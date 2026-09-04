// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)
// Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause

/*
 * GPU donor search for MeshBlock, built on cuBQL (https://github.com/NVIDIA/cuBQL).
 *
 * This is an alternative to MeshBlock::search() with the same inputs and
 * outputs, but a different algorithm:
 *
 *   search()      query OBB -> scan every cell for OBB overlap -> compact the
 *                 survivors into a LIFO list -> recompute their world AABBs ->
 *                 build a 6-D alternating digital tree -> walk the tree once
 *                 per (unique) query point, calling checkContainment().
 *
 *   search_cubql()
 *                 compute one world AABB per cell on the device -> build one
 *                 cuBQL BVH over *all* cells -> one thread per query point,
 *                 traversing the BVH with a degenerate (point) query box and
 *                 running the exact containment test on each candidate.
 *
 * The query-cloud OBB pre-filter is deliberately dropped: it costs a full pass
 * over every cell, which is the same order as building the BVH over every cell
 * in the first place. The BVH is also cached across calls (see gpuMeshDirty),
 * which the ADT is not.
 *
 * Everything geometric is done in double precision, exactly as on the host.
 * Only the BVH bounds are single precision, and those are rounded outward so
 * the broad phase can never produce a false negative.
 */

#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"
#include "cuBQL/builder/cuda/rebinMortonBuilder.h"
#include "cuBQL/traversal/fixedBoxQuery.h"

#include "codetypes.h"
#include "MeshBlock.h"
#include "searchGPU_device.h"

#include <time.h>
#include <stdio.h>
#include <string.h>

#define TIOGA_CUDA_CALL(call)                                           \
  {                                                                     \
    cudaError_t rc__ = (call);                                          \
    if (rc__ != cudaSuccess) {                                          \
      fprintf(stderr,"#tioga: CUDA error %s at %s:%d (%s)\n",           \
              cudaGetErrorString(rc__),__FILE__,__LINE__,#call);        \
      return -1;                                                        \
    }                                                                   \
  }

static inline double gpu_wtime(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC,&ts);
  return (double)ts.tv_sec + 1.0e-9*(double)ts.tv_nsec;
}

/* ================================================================== */
/* persistent device state                                            */
/* ================================================================== */

struct TiogaGpuSearchData
{
  /* --- mesh (rebuilt only when gpuMeshDirty) --- */
  int     nnodes;
  int     ncells;
  int     ntypes;
  int     connSize;
  double *d_x;             /* 3*nnodes                                     */
  int    *d_conn;          /* all types concatenated, already 0-based      */
  int    *d_connOffset;    /* [ntypes]   start of each type inside d_conn  */
  int    *d_cellStart;     /* [ntypes+1] cumulative cell counts            */
  int    *d_nvert;         /* [ntypes]   vertices per cell for each type   */
  double *d_cellRes;       /* [ncells]                                     */
  cuBQL::box3f *d_cellBox; /* [ncells]                                     */
  cuBQL::bvh3f  bvh;
  int     bvhValid;

  /* --- queries (grown on demand) --- */
  int     qCapacity;
  double *d_xsearch;
  int    *d_donorId;

  cudaStream_t stream;
};

/* ================================================================== */
/* kernels                                                            */
/* ================================================================== */

__global__ void k_cellBoxes(cuBQL::box3f *boxes,
                            const double *x, const int *conn,
                            const int *connOffset, const int *cellStart,
                            const int *nvertArr, int ntypes, int ncells)
{
  int cell = blockIdx.x*blockDim.x + threadIdx.x;
  if (cell >= ncells) return;

  int n, i;
  d_cellType(cell, cellStart, ntypes, &n, &i);
  int nvert = nvertArr[n];
  const int *cv = conn + connOffset[n] + (size_t)nvert*i;

  double lo[3], hi[3];
  lo[0] = lo[1] = lo[2] =  BIGVALUE;
  hi[0] = hi[1] = hi[2] = -BIGVALUE;
  for (int m = 0; m < nvert; m++) {
    int i3 = 3*cv[m];
    for (int j = 0; j < 3; j++) {
      double v = x[i3+j];
      lo[j] = fmin(lo[j], v);
      hi[j] = fmax(hi[j], v);
    }
  }
  /* TOL matches the slack the ADT search applies to element bounds */
  cuBQL::box3f b;
  b.lower = cuBQL::vec3f(f_down(lo[0]-TOL), f_down(lo[1]-TOL), f_down(lo[2]-TOL));
  b.upper = cuBQL::vec3f(f_up  (hi[0]+TOL), f_up  (hi[1]+TOL), f_up  (hi[2]+TOL));
  boxes[cell] = b;
}

__global__ __launch_bounds__(128)
void k_search(cuBQL::bvh3f bvh,
              const double *x, const int *conn, const int *connOffset,
              const int *cellStart, const int *nvertArr, int ntypes,
              const double *cellRes,
              const double *xsearch, int *donorId, int nsearch,
              int earlyExit)
{
  int ip = blockIdx.x*blockDim.x + threadIdx.x;
  if (ip >= nsearch) return;

  const double *xp = xsearch + 3*ip;

  cuBQL::box3f qb;
  qb.lower = cuBQL::vec3f(f_down(xp[0]-TOL), f_down(xp[1]-TOL), f_down(xp[2]-TOL));
  qb.upper = cuBQL::vec3f(f_up  (xp[0]+TOL), f_up  (xp[1]+TOL), f_up  (xp[2]+TOL));

  int best = -1;

  auto visit = [&](int cell) -> int {
    int r = d_checkContainment(cell, xp, x, conn, connOffset,
                               cellStart, nvertArr, ntypes, cellRes);
    if (r < 0) return CUBQL_CONTINUE_TRAVERSAL;
    if (r == 0) {
      /* clean hit */
      if (earlyExit) { best = cell; return CUBQL_TERMINATE_TRAVERSAL; }
      best = (best < 0) ? cell : min(best, cell);
      return CUBQL_CONTINUE_TRAVERSAL;
    }
    /* boundary hit on a zero-resolution cell: remember it but keep looking,
       exactly as searchIntersections() does */
    if (best < 0) best = cell;
    return CUBQL_CONTINUE_TRAVERSAL;
  };

  cuBQL::fixedBoxQuery::forEachPrim(visit, bvh, qb);

  donorId[ip] = best;
}

/* ================================================================== */
/* host side                                                          */
/* ================================================================== */

extern "C" {
  void uniquenodes_octree(double *x,int *meshtag,double *rtag,int *itag,int *nn);
}

/* same duplicate-node map as search.C, kept local to avoid a header change */
#ifdef TIOGA_HAS_NODEGID
#include <unordered_map>
#include <algorithm>
namespace {
void gpu_uniquenode_map(uint64_t* node_ids, double* node_res, int* itag, int nnodes)
{
  std::unordered_map<uint64_t,int> lookup;
  for (int i = 0; i < nnodes; i++) {
    auto found = lookup.find(node_ids[i]);
    if (found != lookup.end()) {
      itag[i] = found->second;
      node_res[found->second] = std::max(node_res[found->second], node_res[i]);
    } else {
      lookup[node_ids[i]] = i;
      itag[i] = i;
    }
  }
  for (int i = 0; i < nnodes; i++) node_res[i] = node_res[itag[i]];
}
}
#endif

void MeshBlock::freeCubqlSearchData(void)
{
  if (!gpuData) return;
  TiogaGpuSearchData *g = gpuData;
  if (g->bvhValid) cuBQL::cuda::free(g->bvh, g->stream);
  if (g->d_x)          cudaFree(g->d_x);
  if (g->d_conn)       cudaFree(g->d_conn);
  if (g->d_connOffset) cudaFree(g->d_connOffset);
  if (g->d_cellStart)  cudaFree(g->d_cellStart);
  if (g->d_nvert)      cudaFree(g->d_nvert);
  if (g->d_cellRes)    cudaFree(g->d_cellRes);
  if (g->d_cellBox)    cudaFree(g->d_cellBox);
  if (g->d_xsearch)    cudaFree(g->d_xsearch);
  if (g->d_donorId)    cudaFree(g->d_donorId);
  if (g->stream)       cudaStreamDestroy(g->stream);
  delete g;
  gpuData = NULL;
}

int MeshBlock::search_cubql(void)
{
  double t0, t1, t2;

  searchTimers = SEARCHTIMERS();
  t0 = gpu_wtime();

  if (nsearch == 0) { donorCount = 0; return 0; }

  if (uniform_hex) {
    /* direct logical indexing already beats any tree; leave it on the host */
    search_uniform_hex();
    searchTimers.total = gpu_wtime() - t0;
    return 0;
  }
  if (ihigh != 0) {
    /* high-order containment goes through host callbacks */
    fprintf(stderr,"#tioga: search_cubql() does not support ihigh!=0\n");
    return -2;
  }

  if (!gpuData) {
    gpuData = new TiogaGpuSearchData();
    memset(gpuData, 0, sizeof(TiogaGpuSearchData));
    TIOGA_CUDA_CALL(cudaStreamCreate(&gpuData->stream));
    gpuMeshDirty = 1;
  }
  TiogaGpuSearchData *g = gpuData;
  cudaStream_t s = g->stream;

  /* ---------------- mesh upload (cached) ---------------- */
  t1 = gpu_wtime();
  if (gpuMeshDirty) {
    if (g->bvhValid) { cuBQL::cuda::free(g->bvh, s); g->bvhValid = 0; }
    if (g->d_x)          { cudaFree(g->d_x);          g->d_x = NULL; }
    if (g->d_conn)       { cudaFree(g->d_conn);       g->d_conn = NULL; }
    if (g->d_connOffset) { cudaFree(g->d_connOffset); g->d_connOffset = NULL; }
    if (g->d_cellStart)  { cudaFree(g->d_cellStart);  g->d_cellStart = NULL; }
    if (g->d_nvert)      { cudaFree(g->d_nvert);      g->d_nvert = NULL; }
    if (g->d_cellRes)    { cudaFree(g->d_cellRes);    g->d_cellRes = NULL; }
    if (g->d_cellBox)    { cudaFree(g->d_cellBox);    g->d_cellBox = NULL; }

    /* flatten the per-type connectivity into one 0-based array */
    int *h_connOffset = (int *)malloc(sizeof(int)*ntypes);
    int *h_cellStart  = (int *)malloc(sizeof(int)*(ntypes+1));
    int connSize = 0, cellSum = 0;
    for (int n = 0; n < ntypes; n++) {
      h_connOffset[n] = connSize;
      h_cellStart[n]  = cellSum;
      connSize += nv[n]*nc[n];
      cellSum  += nc[n];
    }
    h_cellStart[ntypes] = cellSum;

    int *h_conn = (int *)malloc(sizeof(int)*(connSize > 0 ? connSize : 1));
    for (int n = 0; n < ntypes; n++) {
      int cnt = nv[n]*nc[n];
      int *dst = h_conn + h_connOffset[n];
      for (int m = 0; m < cnt; m++) dst[m] = vconn[n][m] - BASE;
    }

    g->nnodes   = nnodes;
    g->ncells   = ncells;
    g->ntypes   = ntypes;
    g->connSize = connSize;

    TIOGA_CUDA_CALL(cudaMalloc(&g->d_x,          sizeof(double)*3*nnodes));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_conn,       sizeof(int)*(connSize>0?connSize:1)));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_connOffset, sizeof(int)*ntypes));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellStart,  sizeof(int)*(ntypes+1)));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_nvert,      sizeof(int)*ntypes));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellRes,    sizeof(double)*ncells));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellBox,    sizeof(cuBQL::box3f)*ncells));

    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_x, x, sizeof(double)*3*nnodes,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_conn, h_conn, sizeof(int)*connSize,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_connOffset, h_connOffset, sizeof(int)*ntypes,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_cellStart, h_cellStart, sizeof(int)*(ntypes+1),
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_nvert, nv, sizeof(int)*ntypes,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_cellRes, cellRes, sizeof(double)*ncells,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));

    free(h_conn); free(h_connOffset); free(h_cellStart);
  }
  else if (gpuCoordsDirty) {
    /* Moving mesh: the connectivity is unchanged, so only x[] is re-sent. No
       host staging and no re-walking of the connectivity, which is most of
       what makes a full rebuild expensive. */
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_x, x, sizeof(double)*3*nnodes,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  }

  /* query points */
  if (g->qCapacity < nsearch) {
    if (g->d_xsearch) cudaFree(g->d_xsearch);
    if (g->d_donorId) cudaFree(g->d_donorId);
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_xsearch, sizeof(double)*3*nsearch));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_donorId, sizeof(int)*nsearch));
    g->qCapacity = nsearch;
  }
  TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_xsearch, xsearch, sizeof(double)*3*nsearch,
                                  cudaMemcpyHostToDevice, s));
  TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  searchTimers.transfer = gpu_wtime() - t1;

  /* ---------------- cell AABBs + BVH build (cached) ---------------- */
  t1 = gpu_wtime();
  if (gpuCoordsDirty && g->bvhValid && !gpuMeshDirty && gpuRefit) {
    /* recompute the cell AABBs for the moved coordinates and refit the tree:
       leaf membership is kept, only the bounds move */
    int nb = (ncells + 255)/256;
    k_cellBoxes<<<nb,256,0,s>>>(g->d_cellBox, g->d_x, g->d_conn, g->d_connOffset,
                                g->d_cellStart, g->d_nvert, ntypes, ncells);
    TIOGA_CUDA_CALL(cudaGetLastError());
    cuBQL::cuda::refit(g->bvh, g->d_cellBox, s);
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
    gpuCoordsDirty = 0;
  }
  else if (gpuMeshDirty || gpuCoordsDirty || !g->bvhValid) {
    int nb = (ncells + 255)/256;
    k_cellBoxes<<<nb,256,0,s>>>(g->d_cellBox, g->d_x, g->d_conn, g->d_connOffset,
                                g->d_cellStart, g->d_nvert, ntypes, ncells);
    TIOGA_CUDA_CALL(cudaGetLastError());

    cuBQL::BuildConfig cfg(gpuLeafSize);
    switch (gpuBuilderType) {
    case 1:
      cuBQL::cuda::radixBuilder(g->bvh, g->d_cellBox, (uint32_t)ncells, cfg, s);
      break;
    case 2:
      cuBQL::cuda::rebinRadixBuilder(g->bvh, g->d_cellBox, (uint32_t)ncells, cfg, s);
      break;
    case 3:
      cuBQL::cuda::sahBuilder(g->bvh, g->d_cellBox, (uint32_t)ncells, cfg, s);
      break;
    default:
      cuBQL::gpuBuilder(g->bvh, g->d_cellBox, (uint32_t)ncells, cfg, s);
      break;
    }
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
    g->bvhValid = 1;
    gpuMeshDirty = 0;
    gpuCoordsDirty = 0;
  }
  searchTimers.build = gpu_wtime() - t1;
  searchTimers.candidates = ncells;

  /* ---------------- duplicate query points ---------------- */
  /* Kept on the host so the outputs (xtag, res_search) match search()
     exactly. The GPU searches every point regardless -- duplicates are
     cheaper to just re-search than to gather. */
  t1 = gpu_wtime();
  if (donorId) TIOGA_FREE(donorId);
  donorId = (int *)malloc(sizeof(int)*nsearch);
  if (xtag) TIOGA_FREE(xtag);
  xtag = (int *)malloc(sizeof(int)*nsearch);
#ifdef TIOGA_ENABLE_UNIQUEID
#ifdef TIOGA_HAS_NODEGID
  gpu_uniquenode_map(gid_search.data(), res_search, xtag, nsearch);
#else
  uniquenodes_octree(xsearch, tagsearch, res_search, xtag, &nsearch);
#endif
#else
  for (int i = 0; i < nsearch; i++) xtag[i] = i;
#endif
  searchTimers.dedup = gpu_wtime() - t1;

  /* ---------------- traversal + exact containment ---------------- */
  t1 = gpu_wtime();
  {
    int nb = (nsearch + 127)/128;
    k_search<<<nb,128,0,s>>>(g->bvh, g->d_x, g->d_conn, g->d_connOffset,
                             g->d_cellStart, g->d_nvert, ntypes, g->d_cellRes,
                             g->d_xsearch, g->d_donorId, nsearch, gpuEarlyExit);
    TIOGA_CUDA_CALL(cudaGetLastError());
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  }
  searchTimers.query = gpu_wtime() - t1;

  t2 = gpu_wtime();
  TIOGA_CUDA_CALL(cudaMemcpyAsync(donorId, g->d_donorId, sizeof(int)*nsearch,
                                  cudaMemcpyDeviceToHost, s));
  TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  searchTimers.transfer += gpu_wtime() - t2;

  donorCount = 0;
  ipoint = 0;
  for (int i = 0; i < nsearch; i++) {
    if (donorId[i] > -1) donorCount++;
    ipoint += 3;
  }

  searchTimers.total = gpu_wtime() - t0;
  return 0;
}

/* ================================================================== */
/* One BVH per rank, over the cells of every block                    */
/* ================================================================== */

/*
 * At production block sizes (16^3 active cells) the per-block path above is
 * almost entirely kernel-launch and stream-synchronise latency: 2,000 receptor
 * points is well under a microsecond of real traversal. This batched form
 * builds one BVH over every cell of every block the rank owns and locates
 * every block's points in a single launch.
 *
 * Cells are numbered globally (block 0's cells first, then block 1's, ...).
 * Each cell carries the block it came from and its block-local id, and each
 * query carries the block it belongs to; the traversal visitor rejects any
 * candidate whose block does not match. That keeps donor semantics identical
 * to calling search_cubql() per block even when blocks overlap in space, which
 * is the normal case in overset.
 *
 * Cell connectivity is flattened once into a padded 8-slot vertex list per
 * cell holding indices into the concatenated node array, so the kernel does a
 * direct gather with no per-type lookup.
 */

#define TIOGA_MAXVERT 8

struct TiogaGpuBatchData
{
  int     nblocks;
  int     totalCells;
  int     totalNodes;
  int     totalQueries;
  int     qCapacity;

  double *d_x;          /* 3*totalNodes, all blocks concatenated            */
  int    *d_cellVerts;  /* 8*totalCells, indices into d_x, -1 = unused      */
  int    *d_cellNvert;  /* totalCells                                       */
  int    *d_cellBlock;  /* totalCells                                       */
  int    *d_cellLocal;  /* totalCells, id within that block                 */
  double *d_cellRes;    /* totalCells                                       */
  cuBQL::box3f *d_cellBox;
  cuBQL::bvh3f  bvh;
  int     bvhValid;

  double *d_xsearch;    /* 3*totalQueries                                   */
  int    *d_qBlock;     /* totalQueries                                     */
  int    *d_donorId;    /* totalQueries, block-local cell ids               */

  cudaStream_t stream;
};

static TiogaGpuBatchData *g_batch = NULL;

__global__ void k_cellBoxesBatch(cuBQL::box3f *boxes, const double *x,
                                 const int *cellVerts, const int *cellNvert,
                                 int ncells)
{
  int cell = blockIdx.x*blockDim.x + threadIdx.x;
  if (cell >= ncells) return;

  int nvert = cellNvert[cell];
  const int *cv = cellVerts + TIOGA_MAXVERT*cell;

  double lo[3], hi[3];
  lo[0] = lo[1] = lo[2] =  BIGVALUE;
  hi[0] = hi[1] = hi[2] = -BIGVALUE;
  for (int m = 0; m < nvert; m++) {
    int i3 = 3*cv[m];
    for (int j = 0; j < 3; j++) {
      double v = x[i3+j];
      lo[j] = fmin(lo[j], v);
      hi[j] = fmax(hi[j], v);
    }
  }
  cuBQL::box3f b;
  b.lower = cuBQL::vec3f(f_down(lo[0]-TOL), f_down(lo[1]-TOL), f_down(lo[2]-TOL));
  b.upper = cuBQL::vec3f(f_up  (hi[0]+TOL), f_up  (hi[1]+TOL), f_up  (hi[2]+TOL));
  boxes[cell] = b;
}

__device__ __forceinline__
int d_checkContainmentBatch(int cell, const double *xp, const double *x,
                            const int *cellVerts, const int *cellNvert,
                            const double *cellRes)
{
  double xv[TIOGA_MAXVERT][3], frac[TIOGA_MAXVERT];
  int nvert = cellNvert[cell];
  const int *cv = cellVerts + TIOGA_MAXVERT*cell;
  for (int m = 0; m < nvert; m++) {
    int i3 = 3*cv[m];
    xv[m][0] = x[i3+0]; xv[m][1] = x[i3+1]; xv[m][2] = x[i3+2];
  }
  d_computeNodalWeights(xv, xp, frac, nvert);

  int flagged = 0;
  for (int m = 0; m < nvert; m++) {
    if ((frac[m]+TOL)*(frac[m]-1.0-TOL) > 0.0) return -1;
    if (fabs(frac[m]) < TOL && cellRes[cell] == BIGVALUE) flagged = 1;
  }
  return flagged;
}

__global__ __launch_bounds__(128)
void k_searchBatch(cuBQL::bvh3f bvh, const double *x,
                   const int *cellVerts, const int *cellNvert,
                   const int *cellBlock, const int *cellLocal,
                   const double *cellRes,
                   const double *xsearch, const int *qBlock,
                   int *donorId, int nqueries, int earlyExit)
{
  int ip = blockIdx.x*blockDim.x + threadIdx.x;
  if (ip >= nqueries) return;

  const double *xp = xsearch + 3*ip;
  const int myBlock = qBlock[ip];

  cuBQL::box3f qb;
  qb.lower = cuBQL::vec3f(f_down(xp[0]-TOL), f_down(xp[1]-TOL), f_down(xp[2]-TOL));
  qb.upper = cuBQL::vec3f(f_up  (xp[0]+TOL), f_up  (xp[1]+TOL), f_up  (xp[2]+TOL));

  int best = -1;

  auto visit = [&](int cell) -> int {
    /* a receptor may only be donated to by its own block */
    if (cellBlock[cell] != myBlock) return CUBQL_CONTINUE_TRAVERSAL;
    int r = d_checkContainmentBatch(cell, xp, x, cellVerts, cellNvert, cellRes);
    if (r < 0) return CUBQL_CONTINUE_TRAVERSAL;
    int local = cellLocal[cell];
    if (r == 0) {
      if (earlyExit) { best = local; return CUBQL_TERMINATE_TRAVERSAL; }
      best = (best < 0) ? local : min(best, local);
      return CUBQL_CONTINUE_TRAVERSAL;
    }
    if (best < 0) best = local;
    return CUBQL_CONTINUE_TRAVERSAL;
  };

  cuBQL::fixedBoxQuery::forEachPrim(visit, bvh, qb);
  donorId[ip] = best;
}

void MeshBlock::freeCubqlBatchData(void)
{
  if (!g_batch) return;
  TiogaGpuBatchData *g = g_batch;
  if (g->bvhValid) cuBQL::cuda::free(g->bvh, g->stream);
  if (g->d_x)         cudaFree(g->d_x);
  if (g->d_cellVerts) cudaFree(g->d_cellVerts);
  if (g->d_cellNvert) cudaFree(g->d_cellNvert);
  if (g->d_cellBlock) cudaFree(g->d_cellBlock);
  if (g->d_cellLocal) cudaFree(g->d_cellLocal);
  if (g->d_cellRes)   cudaFree(g->d_cellRes);
  if (g->d_cellBox)   cudaFree(g->d_cellBox);
  if (g->d_xsearch)   cudaFree(g->d_xsearch);
  if (g->d_qBlock)    cudaFree(g->d_qBlock);
  if (g->d_donorId)   cudaFree(g->d_donorId);
  if (g->stream)      cudaStreamDestroy(g->stream);
  delete g;
  g_batch = NULL;
}

int MeshBlock::search_cubql_batch(MeshBlock **blocks, int nblocks,
                                SEARCHTIMERS *timers)
{
  SEARCHTIMERS tm; memset(&tm, 0, sizeof(tm));
  double t0 = gpu_wtime(), t1;

  if (nblocks <= 0) { if (timers) *timers = tm; return 0; }

  for (int b = 0; b < nblocks; b++) {
    if (blocks[b]->ihigh != 0) {
      fprintf(stderr,"#tioga: search_cubql_batch() does not support ihigh!=0\n");
      return -2;
    }
    if (blocks[b]->uniform_hex) {
      fprintf(stderr,"#tioga: search_cubql_batch() does not support uniform_hex "
                     "blocks; call search() on those\n");
      return -2;
    }
  }

  int dirty = 0, coordsDirty = 0;
  if (!g_batch) {
    g_batch = new TiogaGpuBatchData();
    memset(g_batch, 0, sizeof(TiogaGpuBatchData));
    TIOGA_CUDA_CALL(cudaStreamCreate(&g_batch->stream));
    dirty = 1;
  }
  TiogaGpuBatchData *g = g_batch;
  cudaStream_t s = g->stream;

  if (g->nblocks != nblocks) dirty = 1;
  for (int b = 0; b < nblocks; b++) {
    if (blocks[b]->gpuMeshDirty)   dirty = 1;
    if (blocks[b]->gpuCoordsDirty) coordsDirty = 1;
  }
  if (dirty) coordsDirty = 0;   /* a full rebuild subsumes a coordinate update */

  /* ---------------- build the concatenated mesh (cached) ---------------- */
  t1 = gpu_wtime();
  if (dirty || !g->bvhValid) {
    if (g->bvhValid) { cuBQL::cuda::free(g->bvh, s); g->bvhValid = 0; }
    if (g->d_x)         { cudaFree(g->d_x);         g->d_x = NULL; }
    if (g->d_cellVerts) { cudaFree(g->d_cellVerts); g->d_cellVerts = NULL; }
    if (g->d_cellNvert) { cudaFree(g->d_cellNvert); g->d_cellNvert = NULL; }
    if (g->d_cellBlock) { cudaFree(g->d_cellBlock); g->d_cellBlock = NULL; }
    if (g->d_cellLocal) { cudaFree(g->d_cellLocal); g->d_cellLocal = NULL; }
    if (g->d_cellRes)   { cudaFree(g->d_cellRes);   g->d_cellRes = NULL; }
    if (g->d_cellBox)   { cudaFree(g->d_cellBox);   g->d_cellBox = NULL; }

    long totalNodes = 0, totalCells = 0;
    for (int b = 0; b < nblocks; b++) {
      totalNodes += blocks[b]->nnodes;
      totalCells += blocks[b]->ncells;
    }
    g->nblocks    = nblocks;
    g->totalNodes = (int)totalNodes;
    g->totalCells = (int)totalCells;

    double *h_x      = (double *)malloc(sizeof(double)*3*totalNodes);
    int    *h_verts  = (int *)malloc(sizeof(int)*TIOGA_MAXVERT*totalCells);
    int    *h_nvert  = (int *)malloc(sizeof(int)*totalCells);
    int    *h_block  = (int *)malloc(sizeof(int)*totalCells);
    int    *h_local  = (int *)malloc(sizeof(int)*totalCells);
    double *h_res    = (double *)malloc(sizeof(double)*totalCells);

    long nodeOfs = 0, cellOfs = 0;
    for (int b = 0; b < nblocks; b++) {
      MeshBlock *mb = blocks[b];
      memcpy(h_x + 3*nodeOfs, mb->x, sizeof(double)*3*mb->nnodes);
      int lc = 0;
      for (int n = 0; n < mb->ntypes; n++) {
        int nvert = mb->nv[n];
        for (int i = 0; i < mb->nc[n]; i++, lc++) {
          long c = cellOfs + lc;
          int *dst = h_verts + TIOGA_MAXVERT*c;
          for (int m = 0; m < nvert; m++)
            dst[m] = (int)(nodeOfs + mb->vconn[n][nvert*i+m] - BASE);
          for (int m = nvert; m < TIOGA_MAXVERT; m++) dst[m] = -1;
          h_nvert[c] = nvert;
          h_block[c] = b;
          h_local[c] = lc;
          h_res[c]   = mb->cellRes ? mb->cellRes[lc] : 0.0;
        }
      }
      nodeOfs += mb->nnodes;
      cellOfs += mb->ncells;
      mb->gpuMeshDirty = 0;
    }

    TIOGA_CUDA_CALL(cudaMalloc(&g->d_x,         sizeof(double)*3*totalNodes));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellVerts, sizeof(int)*TIOGA_MAXVERT*totalCells));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellNvert, sizeof(int)*totalCells));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellBlock, sizeof(int)*totalCells));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellLocal, sizeof(int)*totalCells));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellRes,   sizeof(double)*totalCells));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_cellBox,   sizeof(cuBQL::box3f)*totalCells));

    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_x, h_x, sizeof(double)*3*totalNodes,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_cellVerts, h_verts,
                                    sizeof(int)*TIOGA_MAXVERT*totalCells,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_cellNvert, h_nvert, sizeof(int)*totalCells,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_cellBlock, h_block, sizeof(int)*totalCells,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_cellLocal, h_local, sizeof(int)*totalCells,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_cellRes, h_res, sizeof(double)*totalCells,
                                    cudaMemcpyHostToDevice, s));
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));

    free(h_x); free(h_verts); free(h_nvert); free(h_block); free(h_local); free(h_res);
  } else if (coordsDirty) {
    /* Moving mesh: the connectivity is unchanged, so only x[] is re-sent.
       Each block's coordinates go straight from its own array into its slice
       of the concatenated buffer -- no host staging, no re-walking the
       connectivity, which is what makes a full rebuild expensive. */
    long nodeOfs = 0;
    for (int b = 0; b < nblocks; b++) {
      MeshBlock *mb = blocks[b];
      TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_x + 3*nodeOfs, mb->x,
                                      sizeof(double)*3*mb->nnodes,
                                      cudaMemcpyHostToDevice, s));
      nodeOfs += mb->nnodes;
      mb->gpuCoordsDirty = 0;
    }
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  }
  tm.transfer = gpu_wtime() - t1;

  /* ---------------- one BVH over every cell (cached) ---------------- */
  t1 = gpu_wtime();
  if (coordsDirty && g->bvhValid && blocks[0]->gpuRefit) {
    /* recompute the cell AABBs for the new coordinates and refit the existing
       tree; leaf membership is kept, only the bounds move */
    int nb = (g->totalCells + 255)/256;
    k_cellBoxesBatch<<<nb,256,0,s>>>(g->d_cellBox, g->d_x, g->d_cellVerts,
                                     g->d_cellNvert, g->totalCells);
    TIOGA_CUDA_CALL(cudaGetLastError());
    cuBQL::cuda::refit(g->bvh, g->d_cellBox, s);
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  } else if (dirty || coordsDirty || !g->bvhValid) {
    if (g->bvhValid) { cuBQL::cuda::free(g->bvh, s); g->bvhValid = 0; }
    int nb = (g->totalCells + 255)/256;
    k_cellBoxesBatch<<<nb,256,0,s>>>(g->d_cellBox, g->d_x, g->d_cellVerts,
                                     g->d_cellNvert, g->totalCells);
    TIOGA_CUDA_CALL(cudaGetLastError());

    cuBQL::BuildConfig cfg(blocks[0]->gpuLeafSize);
    switch (blocks[0]->gpuBuilderType) {
    case 1: cuBQL::cuda::radixBuilder(g->bvh, g->d_cellBox,
                                      (uint32_t)g->totalCells, cfg, s); break;
    case 2: cuBQL::cuda::rebinRadixBuilder(g->bvh, g->d_cellBox,
                                      (uint32_t)g->totalCells, cfg, s); break;
    case 3: cuBQL::cuda::sahBuilder(g->bvh, g->d_cellBox,
                                      (uint32_t)g->totalCells, cfg, s); break;
    default: cuBQL::gpuBuilder(g->bvh, g->d_cellBox,
                                      (uint32_t)g->totalCells, cfg, s); break;
    }
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
    g->bvhValid = 1;
  }
  tm.build = gpu_wtime() - t1;
  tm.candidates = g->totalCells;

  /* ---------------- per-block host bookkeeping + dedup ---------------- */
  t1 = gpu_wtime();
  long totalQ = 0;
  for (int b = 0; b < nblocks; b++) totalQ += blocks[b]->nsearch;
  g->totalQueries = (int)totalQ;

  for (int b = 0; b < nblocks; b++) {
    MeshBlock *mb = blocks[b];
    if (mb->nsearch == 0) continue;
    if (mb->donorId) TIOGA_FREE(mb->donorId);
    mb->donorId = (int *)malloc(sizeof(int)*mb->nsearch);
    if (mb->xtag) TIOGA_FREE(mb->xtag);
    mb->xtag = (int *)malloc(sizeof(int)*mb->nsearch);
#ifdef TIOGA_ENABLE_UNIQUEID
#ifdef TIOGA_HAS_NODEGID
    gpu_uniquenode_map(mb->gid_search.data(), mb->res_search, mb->xtag, mb->nsearch);
#else
    uniquenodes_octree(mb->xsearch, mb->tagsearch, mb->res_search,
                       mb->xtag, &mb->nsearch);
#endif
#else
    for (int i = 0; i < mb->nsearch; i++) mb->xtag[i] = i;
#endif
  }
  tm.dedup = gpu_wtime() - t1;

  if (totalQ == 0) {
    for (int b = 0; b < nblocks; b++) blocks[b]->donorCount = 0;
    tm.total = gpu_wtime() - t0;
    if (timers) *timers = tm;
    return 0;
  }

  /* ---------------- one upload, one launch, one download ------------- */
  t1 = gpu_wtime();
  if (g->qCapacity < (int)totalQ) {
    if (g->d_xsearch) cudaFree(g->d_xsearch);
    if (g->d_qBlock)  cudaFree(g->d_qBlock);
    if (g->d_donorId) cudaFree(g->d_donorId);
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_xsearch, sizeof(double)*3*totalQ));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_qBlock,  sizeof(int)*totalQ));
    TIOGA_CUDA_CALL(cudaMalloc(&g->d_donorId, sizeof(int)*totalQ));
    g->qCapacity = (int)totalQ;
  }

  double *h_q  = (double *)malloc(sizeof(double)*3*totalQ);
  int    *h_qb = (int *)malloc(sizeof(int)*totalQ);
  long qofs = 0;
  for (int b = 0; b < nblocks; b++) {
    MeshBlock *mb = blocks[b];
    if (mb->nsearch == 0) continue;
    memcpy(h_q + 3*qofs, mb->xsearch, sizeof(double)*3*mb->nsearch);
    for (int i = 0; i < mb->nsearch; i++) h_qb[qofs+i] = b;
    qofs += mb->nsearch;
  }
  tm.hostwork += gpu_wtime() - t1;
  t1 = gpu_wtime();
  TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_xsearch, h_q, sizeof(double)*3*totalQ,
                                  cudaMemcpyHostToDevice, s));
  TIOGA_CUDA_CALL(cudaMemcpyAsync(g->d_qBlock, h_qb, sizeof(int)*totalQ,
                                  cudaMemcpyHostToDevice, s));
  TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  tm.transfer += gpu_wtime() - t1;

  t1 = gpu_wtime();
  {
    int nb = ((int)totalQ + 127)/128;
    k_searchBatch<<<nb,128,0,s>>>(g->bvh, g->d_x, g->d_cellVerts, g->d_cellNvert,
                                  g->d_cellBlock, g->d_cellLocal, g->d_cellRes,
                                  g->d_xsearch, g->d_qBlock, g->d_donorId,
                                  (int)totalQ, blocks[0]->gpuEarlyExit);
    TIOGA_CUDA_CALL(cudaGetLastError());
    TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  }
  tm.query = gpu_wtime() - t1;

  t1 = gpu_wtime();
  int *h_donor = (int *)malloc(sizeof(int)*totalQ);
  TIOGA_CUDA_CALL(cudaMemcpyAsync(h_donor, g->d_donorId, sizeof(int)*totalQ,
                                  cudaMemcpyDeviceToHost, s));
  TIOGA_CUDA_CALL(cudaStreamSynchronize(s));
  tm.transfer += gpu_wtime() - t1;
  t1 = gpu_wtime();
  qofs = 0;
  for (int b = 0; b < nblocks; b++) {
    MeshBlock *mb = blocks[b];
    if (mb->nsearch == 0) { mb->donorCount = 0; continue; }
    memcpy(mb->donorId, h_donor + qofs, sizeof(int)*mb->nsearch);
    int dc = 0;
    for (int i = 0; i < mb->nsearch; i++) if (mb->donorId[i] > -1) dc++;
    mb->donorCount = dc;
    mb->ipoint = 3*mb->nsearch;
    qofs += mb->nsearch;
  }
  free(h_donor); free(h_q); free(h_qb);
  tm.hostwork += gpu_wtime() - t1;

  tm.total = gpu_wtime() - t0;
  if (timers) *timers = tm;
  return 0;
}
