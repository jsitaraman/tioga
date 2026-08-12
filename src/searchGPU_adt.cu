// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)
// Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause

/*
 * GPU donor search over tioga's own alternating digital tree.
 *
 * Ported from the exawind-gpu branch of tioga, where the same kernel is
 * inlined into MeshBlock::search(). This is the "native GPU" point of
 * comparison for the cuBQL backend in searchGPU.cu: it changes only where the
 * query loop runs, not what is searched.
 *
 *   search()          host OBB pre-filter -> host buildADT -> host recursive
 *                     tree walk, once per unique query point.
 *
 *   search_adt_gpu()  the same host OBB pre-filter and the same host
 *                     buildADT, then one device thread per query point
 *                     walking the same tree.
 *
 * So the tree, the candidate cell set and the containment math are shared with
 * the host path; only the walk is parallel. That also means the ADT build and
 * the pre-filter stay on the host and stay serial, which is the dominant cost
 * at production block sizes -- exactly the thing this backend exists to
 * measure.
 *
 * Device memory is allocated, filled and released on every call, as in the
 * branch this came from. That is not how one would write it for production,
 * but it keeps the transfer cost visible in the timings instead of hiding it
 * in a cache, and it matches the code being benchmarked.
 *
 * Two differences from the exawind-gpu original, both deliberate:
 *
 *   - The containment test is the one in searchGPU_device.h, which keeps the
 *     "point sits on the boundary of a zero-resolution cell" check that the
 *     original has commented out. Without it the GPU can accept a fringe
 *     donor that the host rejects.
 *
 *   - The stackless fallback tests the point against the subtree bounds in
 *     adtReals before skipping. The original tests it against the bounds of
 *     the single element stored at the node, which is not a bound on the
 *     subtree, so it could skip past cells that do contain the point.
 *
 * Traversal order is breadth first, while the host walk is depth first. Both
 * stop at the first cell that contains the point and is not a boundary hit, so
 * when a point falls inside more than one cell the two can name different
 * (individually valid) donors.
 */

#include "codetypes.h"
#include "MeshBlock.h"
#include "ADT.h"
#include "searchGPU_device.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define TIOGA_CUDA_CALL(call)                                           \
  {                                                                     \
    cudaError_t rc__ = (call);                                          \
    if (rc__ != cudaSuccess) {                                          \
      fprintf(stderr,"#tioga: CUDA error %s at %s:%d (%s)\n",           \
              cudaGetErrorString(rc__),__FILE__,__LINE__,#call);        \
      return -1;                                                        \
    }                                                                   \
  }

/* per-thread node stack, as in the branch this is ported from */
#define ADT_GPU_MAX_STACK 128

static inline double adt_gpu_wtime(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC,&ts);
  return (double)ts.tv_sec + 1.0e-9*(double)ts.tv_nsec;
}

/* Does xp fall inside the box in element[0..ndim)? Element and node bounds are
   stored blocked: the first ndim/2 entries are the lower corner and the rest
   the upper one. */
__device__ __forceinline__
int d_boxContains(const double *element, const double *xp, int ndim)
{
  int flag = 1;
  for (int i = 0; i < ndim/2; i++)
    flag = (flag && (xp[i] >= element[i]-TOL));
  for (int i = ndim/2; i < ndim; i++)
    flag = (flag && (xp[i-ndim/2] <= element[i]+TOL));
  return flag;
}

/* The global extents are stored interleaved instead, as (min,max) per
   dimension, so they need their own test. Mirrors ADT::searchADT(). */
__device__ __forceinline__
int d_extentsContain(const double *adtExtents, const double *xp, int ndim)
{
  int flag = 1;
  for (int i = 0; i < ndim/2; i++)
    flag = (flag && (xp[i] >= adtExtents[2*i]-TOL));
  for (int i = 0; i < ndim/2; i++)
    flag = (flag && (xp[i] <= adtExtents[2*i+1]+TOL));
  return flag;
}

/* ================================================================== */
/* kernel                                                             */
/* ================================================================== */

__global__ __launch_bounds__(128)
void k_adtSearch(const int *adtIntegers, const double *adtReals,
                 const double *adtExtents, const double *coord,
                 const int *elementList,
                 const double *x, const int *conn, const int *connOffset,
                 const int *cellStart, const int *nvertArr, int ntypes,
                 const double *cellRes,
                 const double *xsearch, int *donorId,
                 int nelem, int ndim, int nsearch)
{
  int ip = blockIdx.x*blockDim.x + threadIdx.x;
  if (ip >= nsearch) return;

  const double *xp = xsearch + 3*ip;
  donorId[ip] = -1;

  /* ADT::searchADT() rejects the point against the global extents first */
  if (!d_extentsContain(adtExtents, xp, ndim)) return;

  /* cellIndex[0] is the donor so far, cellIndex[1] flags a boundary hit on a
     zero-resolution cell, which the host keeps searching past. Both follow
     MeshBlock::checkContainment(), including that a later rejected cell
     clears the donor. */
  int cellIndex0 = -1, cellIndex1 = -1;

  int nodeStack[ADT_GPU_MAX_STACK];
  int nstack = 1;
  int bruteSearch = 0;
  nodeStack[0] = 0;

  while (nstack > 0 && !bruteSearch) {
    int mm = 0;
    for (int is = 0; is < nstack && !bruteSearch; is++) {
      int node = nodeStack[is];
      int elem = adtIntegers[4*node];

      if (d_boxContains(coord + ndim*elem, xp, ndim)) {
        int cell = elementList[elem];
        int r = d_checkContainment(cell, xp, x, conn, connOffset,
                                   cellStart, nvertArr, ntypes, cellRes);
        cellIndex0 = (r < 0) ? -1 : cell;
        cellIndex1 = (r < 0) ?  0 : r;
        if (cellIndex0 > -1 && cellIndex1 == 0) { donorId[ip] = cellIndex0; return; }
      }
      /* queue whichever children still bound the point */
      for (int d = 1; d < 3; d++) {
        int nodeChild = adtIntegers[4*node+d];
        if (nodeChild > -1 &&
            d_boxContains(adtReals + ndim*nodeChild, xp, ndim)) {
          if (mm+nstack > ADT_GPU_MAX_STACK-1) {
            /* frontier does not fit: restart as a stackless walk */
            cellIndex0 = cellIndex1 = -1;
            bruteSearch = 1;
            break;
          }
          nodeStack[mm+nstack] = nodeChild;
          mm++;
        }
      }
    }
    if (!bruteSearch) {
      for (int j = 0; j < mm; j++) nodeStack[j] = nodeStack[j+nstack];
      nstack = mm;
    }
  }

  if (bruteSearch) {
    /* Nodes are numbered in pre-order and adtIntegers[4*node+3] holds the
       size of the subtree rooted at node, so advancing by it skips the whole
       subtree. Visiting in this order reproduces the host walk exactly. */
    int node = 0;
    while (node < nelem) {
      if (d_boxContains(adtReals + ndim*node, xp, ndim)) {
        int elem = adtIntegers[4*node];
        if (d_boxContains(coord + ndim*elem, xp, ndim)) {
          int cell = elementList[elem];
          int r = d_checkContainment(cell, xp, x, conn, connOffset,
                                     cellStart, nvertArr, ntypes, cellRes);
          cellIndex0 = (r < 0) ? -1 : cell;
          cellIndex1 = (r < 0) ?  0 : r;
          if (cellIndex0 > -1 && cellIndex1 == 0) { donorId[ip] = cellIndex0; return; }
        }
        node = node + 1;
      }
      else {
        node = node + adtIntegers[4*node+3];
      }
    }
  }

  /* nothing cleanly contained the point: fall back to a boundary hit if the
     last cell tested produced one, matching ADT::searchADT() */
  donorId[ip] = cellIndex0;
}

/* ================================================================== */
/* host side                                                          */
/* ================================================================== */

int MeshBlock::search_adt_gpu(void)
{
  double t1, tq0;

  if (nsearch == 0) return 0;
  if (ihigh != 0) {
    fprintf(stderr,"#tioga: search_adt_gpu() does not support ihigh!=0\n");
    return -2;
  }
  if (!adt) return -2;

  int nelem = adt->get_nelem();
  int ndim  = adt->get_ndim();
  if (nelem == 0) {
    for (int i = 0; i < nsearch; i++) donorId[i] = -1;
    return 0;
  }

  /* ---- flatten the per-type connectivity into one 0-based array ---- */
  t1 = adt_gpu_wtime();
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
  searchTimers.hostwork += adt_gpu_wtime() - t1;

  /* ---- allocate and upload, every call ---- */
  t1 = adt_gpu_wtime();
  int    *d_adtIntegers = NULL, *d_elementList = NULL;
  double *d_adtReals = NULL, *d_adtExtents = NULL, *d_coord = NULL;
  double *d_x = NULL, *d_cellRes = NULL, *d_xsearch = NULL;
  int    *d_conn = NULL, *d_connOffset = NULL, *d_cellStart = NULL;
  int    *d_nvert = NULL, *d_donorId = NULL;

  TIOGA_CUDA_CALL(cudaMalloc(&d_adtIntegers, sizeof(int)*4*nelem));
  TIOGA_CUDA_CALL(cudaMalloc(&d_adtReals,    sizeof(double)*nelem*ndim));
  TIOGA_CUDA_CALL(cudaMalloc(&d_adtExtents,  sizeof(double)*ndim));
  TIOGA_CUDA_CALL(cudaMalloc(&d_coord,       sizeof(double)*nelem*ndim));
  TIOGA_CUDA_CALL(cudaMalloc(&d_elementList, sizeof(int)*nelem));
  TIOGA_CUDA_CALL(cudaMalloc(&d_x,           sizeof(double)*3*nnodes));
  TIOGA_CUDA_CALL(cudaMalloc(&d_conn,        sizeof(int)*(connSize>0?connSize:1)));
  TIOGA_CUDA_CALL(cudaMalloc(&d_connOffset,  sizeof(int)*ntypes));
  TIOGA_CUDA_CALL(cudaMalloc(&d_cellStart,   sizeof(int)*(ntypes+1)));
  TIOGA_CUDA_CALL(cudaMalloc(&d_nvert,       sizeof(int)*ntypes));
  TIOGA_CUDA_CALL(cudaMalloc(&d_cellRes,     sizeof(double)*ncells));
  TIOGA_CUDA_CALL(cudaMalloc(&d_xsearch,     sizeof(double)*3*nsearch));
  TIOGA_CUDA_CALL(cudaMalloc(&d_donorId,     sizeof(int)*nsearch));

  TIOGA_CUDA_CALL(cudaMemcpy(d_adtIntegers, adt->get_Integers(),
                             sizeof(int)*4*nelem, cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_adtReals, adt->get_Reals(),
                             sizeof(double)*nelem*ndim, cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_adtExtents, adt->get_Extents(),
                             sizeof(double)*ndim, cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_coord, elementBbox,
                             sizeof(double)*nelem*ndim, cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_elementList, elementList,
                             sizeof(int)*nelem, cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_x, x, sizeof(double)*3*nnodes,
                             cudaMemcpyHostToDevice));
  if (connSize > 0)
    TIOGA_CUDA_CALL(cudaMemcpy(d_conn, h_conn, sizeof(int)*connSize,
                               cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_connOffset, h_connOffset, sizeof(int)*ntypes,
                             cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_cellStart, h_cellStart, sizeof(int)*(ntypes+1),
                             cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_nvert, nv, sizeof(int)*ntypes,
                             cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_cellRes, cellRes, sizeof(double)*ncells,
                             cudaMemcpyHostToDevice));
  TIOGA_CUDA_CALL(cudaMemcpy(d_xsearch, xsearch, sizeof(double)*3*nsearch,
                             cudaMemcpyHostToDevice));
  searchTimers.transfer += adt_gpu_wtime() - t1;

  /* ---- traverse ---- */
  tq0 = adt_gpu_wtime();
  int block_size = 128;
  int n_blocks = (nsearch + block_size - 1)/block_size;
  k_adtSearch<<<n_blocks,block_size>>>(d_adtIntegers, d_adtReals, d_adtExtents,
                                       d_coord, d_elementList,
                                       d_x, d_conn, d_connOffset, d_cellStart,
                                       d_nvert, ntypes, d_cellRes,
                                       d_xsearch, d_donorId,
                                       nelem, ndim, nsearch);
  TIOGA_CUDA_CALL(cudaGetLastError());
  TIOGA_CUDA_CALL(cudaDeviceSynchronize());
  searchTimers.query += adt_gpu_wtime() - tq0;

  t1 = adt_gpu_wtime();
  TIOGA_CUDA_CALL(cudaMemcpy(donorId, d_donorId, sizeof(int)*nsearch,
                             cudaMemcpyDeviceToHost));
  searchTimers.transfer += adt_gpu_wtime() - t1;

  cudaFree(d_adtIntegers); cudaFree(d_adtReals);  cudaFree(d_adtExtents);
  cudaFree(d_coord);       cudaFree(d_elementList);
  cudaFree(d_x);           cudaFree(d_conn);      cudaFree(d_connOffset);
  cudaFree(d_cellStart);   cudaFree(d_nvert);     cudaFree(d_cellRes);
  cudaFree(d_xsearch);     cudaFree(d_donorId);
  free(h_conn); free(h_connOffset); free(h_cellStart);

  /* ---- duplicate query points inherit their representative's donor, and the
     donor accounting search() would otherwise have done in its own loop ---- */
  t1 = adt_gpu_wtime();
  donorCount = 0;
  ipoint = 0;
  for (int i = 0; i < nsearch; i++) {
    donorId[i] = donorId[xtag[i]];
    if (donorId[i] > -1) donorCount++;
    ipoint += 3;
  }
  searchTimers.hostwork += adt_gpu_wtime() - t1;

  return 0;
}

/*
 * Create the CUDA context up front.
 *
 * The first device call in a process pays for context creation, which is
 * order a second per rank and has nothing to do with searching. Applications
 * call performConnectivity() once per mesh motion, so that cost is amortised
 * away in production but would otherwise land entirely inside the first (and
 * in the test drivers, only) measured search. Doing it here keeps the search
 * timings comparable between backends.
 */
extern "C" void tioga_gpu_context_init(void)
{
  static int done = 0;
  if (done) return;
  done = 1;
  cudaFree(0);
}
