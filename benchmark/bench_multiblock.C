// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

/*
 * Many mesh blocks per rank, one CUDA context per GPU.
 *
 * This is the configuration a production overset run actually has: a rank owns
 * a number of small blocks (16^3 active cells is typical), and calls search()
 * on each of them once per connectivity update. bench_search.C measures one
 * large block, which flatters the GPU -- with small blocks the per-block
 * overhead (H2D upload, BVH build, kernel launch, stream sync) is paid N times
 * and the traversal itself is a rounding error.
 *
 * All blocks live in one process, so there is exactly one CUDA context, and
 * each block keeps its own device state and stream (MeshBlock::gpuData).
 *
 * Reported per run:
 *   cpu_total      sum over blocks of search()          -- one core
 *   gpu_total      sum over blocks of search_gpu()      -- mesh re-uploaded,
 *                                                          BVH rebuilt
 *   gpu_cached     same, with mesh and BVH resident from the previous call
 *
 * Output is one CSV row per configuration; see --header.
 */
#include "meshgen.h"
#include <dirent.h>
#include <unistd.h>

/* Filesystem barrier across independently launched ranks.
   Without this, ranks reach their GPU section at different times, each finds
   an idle GPU, and "total queries / max wall time" measures N x the
   uncontended single-rank rate rather than a sustained throughput. */
static const char *g_barrierDir = NULL;
static int g_nranks = 0, g_rank = 0;

static void barrier(const char *tag)
{
  if (!g_barrierDir || g_nranks <= 1) return;
  char f[1024];
  snprintf(f, sizeof(f), "%s/%s.%d", g_barrierDir, tag, g_rank);
  FILE *fp = fopen(f, "w");
  if (fp) { fputc('x', fp); fclose(fp); }
  const size_t tl = strlen(tag);
  for (;;) {
    int c = 0;
    DIR *d = opendir(g_barrierDir);
    if (d) {
      struct dirent *e;
      while ((e = readdir(d)))
        if (!strncmp(e->d_name, tag, tl) && e->d_name[tl] == '.') c++;
      closedir(d);
    }
    if (c >= g_nranks) break;
    usleep(500);
  }
}

struct Block
{
  Mesh                mesh;
  MeshBlock          *mb;
  std::vector<double> xs;
  std::vector<int>    donorCpu;
};

static void attachQueries(MeshBlock &mb, const std::vector<double> &xs, int nq)
{
  /* MeshBlock's destructor owns these, so they must come from malloc */
  if (mb.xsearch)    TIOGA_FREE(mb.xsearch);
  if (mb.tagsearch)  TIOGA_FREE(mb.tagsearch);
  if (mb.res_search) TIOGA_FREE(mb.res_search);
  mb.nsearch    = nq;
  mb.xsearch    = (double *)malloc(sizeof(double)*3*nq);
  mb.tagsearch  = (int *)malloc(sizeof(int)*nq);
  mb.res_search = (double *)malloc(sizeof(double)*nq);
  memcpy(mb.xsearch, xs.data(), sizeof(double)*3*nq);
  for (int i = 0; i < nq; i++) { mb.tagsearch[i] = 2; mb.res_search[i] = 1.0; }
}

static void usage(const char *p)
{
  fprintf(stderr,
    "usage: %s [options]\n"
    "  --nblocks N     mesh blocks owned by this rank              [64]\n"
    "  --nx N          cells per direction per block (N^3 cells)   [16]\n"
    "  --nq N          query points per block                      [2000]\n"
    "  --eltype T      hex | prism | tet | mixed                   [hex]\n"
    "  --qfrac F       query cloud edge as fraction of a block     [1.0]\n"
    "  --overlap F     block lattice pitch = 1-F, so neighbouring blocks\n"
    "                  overlap in space (exercises the per-block donor\n"
    "                  restriction in --batch)                        [0]\n"
    "  --reps N        timed repetitions (best is reported)        [3]\n"
    "  --seed N        RNG seed                                    [12345]\n"
    "  --builder N     0=spatial median 1=radix 2=rebin 3=SAH      [1]\n"
    "  --skip-dedup    GPU skips the host duplicate-point pass\n"
    "  --batch         one BVH per rank over all blocks, one launch\n"
    "  --moving        emulate a moving mesh: coordinates re-sent and the\n"
    "                  BVH refit each call, connectivity left resident\n"
    "  --no-refit      on --moving, rebuild the BVH instead of refitting\n"
    "  --move F        displace every node by F cell-widths once, before the\n"
    "                  timed phases, so refit runs on coordinates the tree\n"
    "                  was not built for (correctness check)          [0]\n"
    "  --cpu-only      skip the GPU run\n"
    "  --gpu-only      skip the CPU run (no donor cross-check)\n"
    "  --barrier DIR   shared dir for the cross-rank barrier\n"
    "  --nranks N      number of ranks participating in the barrier\n"
    "  --rank N        this process's rank id\n"
    "  --header        print the CSV header and exit\n"
    "  --verbose       print a human readable breakdown too\n", p);
}

static const char *CSV_HEADER =
  "nblocks,nx,eltype,cells_per_block,cells_total,nq_per_block,nq_total,builder,batch,"
  "cpu_total,cpu_filter,cpu_build,cpu_dedup,cpu_query,"
  "gpu_total,gpu_xfer,gpu_build,gpu_dedup,gpu_query,gpu_total_cached,"
  "donors_cpu,donors_gpu,mismatch,mismatch_real,"
  "speedup_total,speedup_cached,"
  "cpu_Mqps,gpu_Mqps,gpu_cached_Mqps,us_per_block_cached,"
  "reps,cpu_window,gpu_window,gpu_cached_window,"
  "cpu_compute,gpu_compute,gpu_compute_cached,"
  "cpu_compute_Mqps,gpu_compute_Mqps,gpu_compute_cached_Mqps,"
  "speedup_compute,speedup_compute_cached,"
  "gpu_hostwork,gpu_compute_plus_hostwork,gpu_compute_plus_hostwork_cached,"
  "speedup_compute_plus_hostwork";

int main(int argc, char **argv)
{
  int nblocks = 64, nx = 16, nq = 2000, reps = 3, builder = 1;
  int cpuOnly = 0, gpuOnly = 0, verbose = 0, moving = 0, refit = 1;
  int eltype = EL_HEX, skipDedup = 0, batch = 0;
  double qfrac = 1.0, perturb = 0.15, overlap = 0.0, moveBy = 0.0;
  uint64_t seed = 12345;

  for (int a = 1; a < argc; a++) {
    if      (!strcmp(argv[a],"--nblocks") && a+1 < argc) nblocks = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--nx")      && a+1 < argc) nx = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--nq")      && a+1 < argc) nq = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--qfrac")   && a+1 < argc) qfrac = atof(argv[++a]);
    else if (!strcmp(argv[a],"--overlap") && a+1 < argc) overlap = atof(argv[++a]);
    else if (!strcmp(argv[a],"--perturb") && a+1 < argc) perturb = atof(argv[++a]);
    else if (!strcmp(argv[a],"--reps")    && a+1 < argc) reps = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--seed")    && a+1 < argc) seed = strtoull(argv[++a],NULL,10);
    else if (!strcmp(argv[a],"--builder") && a+1 < argc) builder = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--eltype")  && a+1 < argc) {
      const char *e = argv[++a];
      if      (!strcmp(e,"hex"))   eltype = EL_HEX;
      else if (!strcmp(e,"prism")) eltype = EL_PRISM;
      else if (!strcmp(e,"tet"))   eltype = EL_TET;
      else if (!strcmp(e,"mixed")) eltype = EL_MIXED;
      else { fprintf(stderr,"unknown --eltype %s\n",e); return 1; }
    }
    else if (!strcmp(argv[a],"--skip-dedup")) skipDedup = 1;
    else if (!strcmp(argv[a],"--batch"))      batch = 1;
    else if (!strcmp(argv[a],"--moving"))    moving = 1;
    else if (!strcmp(argv[a],"--no-refit"))  refit = 0;
    else if (!strcmp(argv[a],"--move") && a+1 < argc) moveBy = atof(argv[++a]);
    else if (!strcmp(argv[a],"--cpu-only"))   cpuOnly = 1;
    else if (!strcmp(argv[a],"--gpu-only"))   gpuOnly = 1;
    else if (!strcmp(argv[a],"--barrier") && a+1 < argc) g_barrierDir = argv[++a];
    else if (!strcmp(argv[a],"--nranks")  && a+1 < argc) g_nranks = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--rank")    && a+1 < argc) g_rank = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--verbose"))    verbose = 1;
    else if (!strcmp(argv[a],"--header")) { printf("%s\n",CSV_HEADER); return 0; }
    else { usage(argv[0]); return 1; }
  }

  /* ---------------- build the blocks ---------------- */
  double tb = wtime();
  std::vector<Block> blk(nblocks);
  int side = 1;
  while (side*side*side < nblocks) side++;      /* lay blocks out in a lattice */
  for (int b = 0; b < nblocks; b++) {
    const double pitch = 1.0 - overlap;
    double org[3] = { pitch*(b % side),
                      pitch*((b / side) % side),
                      pitch*(b / (side*side)) };
    buildMesh(blk[b].mesh, nx, nx, nx, perturb, eltype, org);
    buildQueries(blk[b].mesh, nq, qfrac, 0.0, seed + 7919ull*b, blk[b].xs);
    blk[b].mb = new MeshBlock();
    blk[b].mb->setData(1, blk[b].mesh.nnodes, blk[b].mesh.x, blk[b].mesh.iblank,
                       0, blk[b].mesh.nobc, NULL, blk[b].mesh.obcnode,
                       blk[b].mesh.ntypes, blk[b].mesh.nv, blk[b].mesh.nc,
                       blk[b].mesh.vconn);
    blk[b].mb->preprocess();
    blk[b].mb->gpuBuilderType = builder;
    blk[b].mb->gpuSkipDedup   = skipDedup;
    blk[b].mb->gpuRefit       = refit;
  }
  double tsetup = wtime() - tb;

  /* Displace the mesh once, before anything is searched, so that the CPU
     reference and the GPU both see the moved configuration -- but the GPU's
     BVH gets built on it and then refit against it, which is what a real
     moving-mesh step does. */
  if (moveBy != 0.0) {
    const double d = moveBy / nx;
    for (int b = 0; b < nblocks; b++)
      for (int i = 0; i < blk[b].mesh.nnodes; i++) {
        blk[b].mesh.x[3*i+0] += d;
        blk[b].mesh.x[3*i+1] += 0.5*d;
        blk[b].mesh.x[3*i+2] -= 0.25*d;
      }
  }

  std::vector<MeshBlock *> mbs(nblocks);
  for (int b = 0; b < nblocks; b++) mbs[b] = blk[b].mb;

  const long cellsPer  = blk[0].mesh.ncells;
  const long cellsTot  = cellsPer * nblocks;
  const long nqTot     = (long)nq * nblocks;

  if (verbose)
    printf("# %d blocks x %d^3 (%ld cells, %d queries) each = %ld cells, "
           "%ld queries | setup %.2f s\n",
           nblocks, nx, cellsPer, nq, cellsTot, nqTot, tsetup);

  /* ---------------- CPU: loop over all blocks ---------------- */
  SEARCHTIMERS cpuBest; memset(&cpuBest, 0, sizeof(cpuBest));
  double cpuBestTotal = 1e300, cpuWindow = 0.0, cpuW0 = 0.0;
  int donorsCpu = 0;
  for (int r = 0; r < (gpuOnly ? 0 : reps + 1); r++) {   /* r==0 warm-up */
    if (r == 1) { barrier("cpu"); cpuW0 = wtime(); }
    SEARCHTIMERS acc; memset(&acc, 0, sizeof(acc));
    int donors = 0;
    double t0 = wtime();
    for (int b = 0; b < nblocks; b++) {
      attachQueries(*blk[b].mb, blk[b].xs, nq);
      blk[b].mb->search();
      const SEARCHTIMERS &t = blk[b].mb->searchTimers;
      acc.filter += t.filter; acc.build += t.build;
      acc.dedup  += t.dedup;  acc.query += t.query;
      donors += blk[b].mb->donorCount;
      if (r == 0)
        blk[b].donorCpu.assign(blk[b].mb->donorId, blk[b].mb->donorId + nq);
    }
    acc.total = wtime() - t0;
    if (r == 0) { donorsCpu = donors; continue; }
    if (acc.total < cpuBestTotal) { cpuBestTotal = acc.total; cpuBest = acc; }
  }
  if (gpuOnly) { cpuBestTotal = 0.0; cpuWindow = 0.0; }
  else cpuWindow = wtime() - cpuW0;

  /* ---------------- GPU: loop over all blocks, one context ---------------- */
  SEARCHTIMERS gpuBest; memset(&gpuBest, 0, sizeof(gpuBest));
  double gpuBestTotal = 1e300, gpuBestCached = 1e300;
  int gpuOk = 0, donorsGpu = -1, mismatch = 0, bad = 0;
  double gpuWindow = 0.0, gpuCachedWindow = 0.0;
  SEARCHTIMERS gpuCached; memset(&gpuCached, 0, sizeof(gpuCached));

  if (!cpuOnly) {
    /* warm-up: creates device state, uploads, builds */
    int rc = 0;
    for (int b = 0; b < nblocks; b++) attachQueries(*blk[b].mb, blk[b].xs, nq);
    if (batch) rc = MeshBlock::search_gpu_batch(mbs.data(), nblocks, NULL);
    else for (int b = 0; b < nblocks; b++) rc |= blk[b].mb->search_gpu();
    gpuOk = (rc == 0);

    (void)0;
  }

  /* The barriers below are entered unconditionally -- including by ranks whose
     GPU path failed or which are running --cpu-only -- so that one rank
     bailing out cannot strand the others. */
  {
    barrier("gpu");
    double gpuW0 = wtime();
    if (gpuOk) {
      /* (a) moving mesh: re-upload + rebuild every block, every call */
      for (int r = 0; r < reps; r++) {
        SEARCHTIMERS acc; memset(&acc, 0, sizeof(acc));
        double t0 = wtime();
        for (int b = 0; b < nblocks; b++) {
          attachQueries(*blk[b].mb, blk[b].xs, nq);
          if (moving) blk[b].mb->gpuCoordsDirty = 1;
          else        blk[b].mb->gpuMeshDirty   = 1;
        }
        if (batch) {
          MeshBlock::search_gpu_batch(mbs.data(), nblocks, &acc);
        } else {
          for (int b = 0; b < nblocks; b++) {
            blk[b].mb->search_gpu();
            const SEARCHTIMERS &t = blk[b].mb->searchTimers;
            acc.transfer += t.transfer; acc.build += t.build;
            acc.dedup    += t.dedup;    acc.query += t.query;
            acc.hostwork += t.hostwork;
          }
        }
        acc.total = wtime() - t0;
        if (acc.total < gpuBestTotal) { gpuBestTotal = acc.total; gpuBest = acc; }
      }
    }
    gpuWindow = wtime() - gpuW0;

    barrier("gpuc");
    double gpucW0 = wtime();
    if (gpuOk) {
      /* (b) static mesh: everything resident from the previous call */
      for (int r = 0; r < reps; r++) {
        SEARCHTIMERS cacc; memset(&cacc, 0, sizeof(cacc));
        double t0 = wtime();
        for (int b = 0; b < nblocks; b++) attachQueries(*blk[b].mb, blk[b].xs, nq);
        if (batch) {
          MeshBlock::search_gpu_batch(mbs.data(), nblocks, &cacc);
        } else {
          for (int b = 0; b < nblocks; b++) {
            blk[b].mb->search_gpu();
            const SEARCHTIMERS &t = blk[b].mb->searchTimers;
            cacc.transfer += t.transfer; cacc.build += t.build;
            cacc.dedup    += t.dedup;    cacc.query += t.query;
            cacc.hostwork += t.hostwork;
          }
        }
        double el = wtime() - t0;
        if (el < gpuBestCached) { gpuBestCached = el; gpuCached = cacc; }
      }
    }
    gpuCachedWindow = wtime() - gpucW0;

    if (gpuOk) {
      donorsGpu = 0;
      for (int b = 0; b < nblocks; b++) {
        donorsGpu += blk[b].mb->donorCount;
        if (gpuOnly) continue;
        for (int i = 0; i < nq; i++) {
          int c = blk[b].donorCpu[i], g = blk[b].mb->donorId[i];
          if (c == g) continue;
          mismatch++;
          if (!(donorIsValid(blk[b].mesh, c, &blk[b].xs[3*i]) &&
                donorIsValid(blk[b].mesh, g, &blk[b].xs[3*i]))) bad++;
        }
      }
    }
  }

  const double cpuMq    = cpuBestTotal > 0 ? nqTot / cpuBestTotal / 1e6 : 0.0;
  const double gpuMq    = gpuOk ? nqTot / gpuBestTotal  / 1e6 : 0.0;
  const double gpuCacMq = gpuOk ? nqTot / gpuBestCached / 1e6 : 0.0;
  const double usBlk    = gpuOk ? gpuBestCached*1e6/nblocks : 0.0;

  /* Compute-only comparison: the arithmetic each backend performs, with host
     <-> device transfers and the shared host duplicate-point pass excluded on
     both sides.
       CPU : query-OBB filter + ADT build + tree walk + containment
       GPU : cell-AABB kernel + BVH build/refit + traversal kernel
     The GPU figures are wall time around the kernels including launch and
     stream synchronise, which at these sizes is under 1 % of the total. */
  const double cpuCompute    = cpuBest.filter + cpuBest.build + cpuBest.query;
  const double gpuCompute    = gpuOk ? gpuBest.build   + gpuBest.query   : 0.0;
  const double gpuComputeCac = gpuOk ? gpuCached.build + gpuCached.query : 0.0;
  /* compute plus the host packing/unpacking that sits between the kernels */
  const double gpuComputeHW  = gpuCompute    + (gpuOk ? gpuBest.hostwork   : 0.0);
  const double gpuCompCacHW  = gpuComputeCac + (gpuOk ? gpuCached.hostwork : 0.0);
  const double cpuComputeMq  = cpuCompute    > 0 ? nqTot/cpuCompute/1e6    : 0.0;
  const double gpuComputeMq  = gpuCompute    > 0 ? nqTot/gpuCompute/1e6    : 0.0;
  const double gpuCompCacMq  = gpuComputeCac > 0 ? nqTot/gpuComputeCac/1e6 : 0.0;

  if (verbose) {
    if (!gpuOnly) printf("# CPU  %8.4f s  (filter %6.4f  ADT %6.4f  dedup %6.4f  walk %6.4f)"
           "  donors %d  %.3f M q/s\n",
           cpuBestTotal, cpuBest.filter, cpuBest.build, cpuBest.dedup,
           cpuBest.query, donorsCpu, cpuMq);
    if (gpuOk) {
      printf("# GPU  %8.4f s  (xfer   %6.4f  BVH %6.4f  dedup %6.4f  walk %6.4f)"
             "  donors %d  %.3f M q/s\n",
             gpuBestTotal, gpuBest.transfer, gpuBest.build, gpuBest.dedup,
             gpuBest.query, donorsGpu, gpuMq);
      printf("# + host pack/unpack between kernels: %8.4f s  =>  compute+host "
             "%8.4f s (%.1f M q/s), cached %8.4f s (%.1f M q/s)\n",
             gpuBest.hostwork, gpuComputeHW,
             gpuComputeHW > 0 ? nqTot/gpuComputeHW/1e6 : 0.0,
             gpuCompCacHW, gpuCompCacHW > 0 ? nqTot/gpuCompCacHW/1e6 : 0.0);
      printf("# COMPUTE ONLY (no H2D/D2H, no dedup):  CPU %8.4f s (%.1f M q/s)  "
             "GPU %8.4f s (%.1f M q/s) = %.1fx   |   GPU cached %8.4f s "
             "(%.1f M q/s) = %.1fx\n",
             cpuCompute, cpuComputeMq, gpuCompute, gpuComputeMq,
             gpuCompute > 0 ? cpuCompute/gpuCompute : 0.0,
             gpuComputeCac, gpuCompCacMq,
             gpuComputeCac > 0 ? cpuCompute/gpuComputeCac : 0.0,
         gpuOk ? gpuBest.hostwork : 0.0, gpuComputeHW, gpuCompCacHW,
         gpuComputeHW > 0 ? cpuCompute/gpuComputeHW : 0.0);
      printf("# GPU cached %8.4f s  %.3f M q/s  %.1f us/block  |  "
             "speedup %.2fx / %.2fx cached  |  mismatches %d (%d real)\n",
             gpuBestCached, gpuCacMq, usBlk,
             cpuBestTotal > 0 ? cpuBestTotal/gpuBestTotal  : 0.0,
             cpuBestTotal > 0 ? cpuBestTotal/gpuBestCached : 0.0,
             mismatch, bad);
    }
  }

  printf("%d,%d,%s,%ld,%ld,%d,%ld,%d,%d,"
         "%.6f,%.6f,%.6f,%.6f,%.6f,"
         "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
         "%d,%d,%d,%d,"
         "%.3f,%.3f,"
         "%.4f,%.4f,%.4f,%.2f,"
         "%d,%.6f,%.6f,%.6f,"
         "%.6f,%.6f,%.6f,"
         "%.4f,%.4f,%.4f,"
         "%.3f,%.3f,"
         "%.6f,%.6f,%.6f,%.3f\n",
         nblocks, nx, elTypeName(eltype), cellsPer, cellsTot, nq, nqTot, builder, batch,
         cpuBestTotal, cpuBest.filter, cpuBest.build, cpuBest.dedup, cpuBest.query,
         gpuOk ? gpuBestTotal : 0.0, gpuOk ? gpuBest.transfer : 0.0,
         gpuOk ? gpuBest.build : 0.0, gpuOk ? gpuBest.dedup : 0.0,
         gpuOk ? gpuBest.query : 0.0, gpuOk ? gpuBestCached : 0.0,
         donorsCpu, donorsGpu, gpuOk ? mismatch : -1, gpuOk ? bad : -1,
         (gpuOk && cpuBestTotal > 0) ? cpuBestTotal/gpuBestTotal  : 0.0,
         (gpuOk && cpuBestTotal > 0) ? cpuBestTotal/gpuBestCached : 0.0,
         cpuMq, gpuMq, gpuCacMq, usBlk,
         reps, cpuWindow, gpuWindow, gpuCachedWindow,
         cpuCompute, gpuCompute, gpuComputeCac,
         cpuComputeMq, gpuComputeMq, gpuCompCacMq,
         gpuCompute    > 0 ? cpuCompute/gpuCompute    : 0.0,
         gpuComputeCac > 0 ? cpuCompute/gpuComputeCac : 0.0,
         gpuOk ? gpuBest.hostwork : 0.0, gpuComputeHW, gpuCompCacHW,
         gpuComputeHW > 0 ? cpuCompute/gpuComputeHW : 0.0);
  fflush(stdout);

  MeshBlock::freeGpuBatchData();
  for (int b = 0; b < nblocks; b++) { delete blk[b].mb; freeMesh(blk[b].mesh); }
  return 0;
}
