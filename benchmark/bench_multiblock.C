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
    "  --reps N        timed repetitions (best is reported)        [3]\n"
    "  --seed N        RNG seed                                    [12345]\n"
    "  --builder N     0=spatial median 1=radix 2=rebin 3=SAH      [1]\n"
    "  --skip-dedup    GPU skips the host duplicate-point pass\n"
    "  --cpu-only      skip the GPU run\n"
    "  --header        print the CSV header and exit\n"
    "  --verbose       print a human readable breakdown too\n", p);
}

static const char *CSV_HEADER =
  "nblocks,nx,eltype,cells_per_block,cells_total,nq_per_block,nq_total,builder,"
  "cpu_total,cpu_filter,cpu_build,cpu_dedup,cpu_query,"
  "gpu_total,gpu_xfer,gpu_build,gpu_dedup,gpu_query,gpu_total_cached,"
  "donors_cpu,donors_gpu,mismatch,mismatch_real,"
  "speedup_total,speedup_cached,"
  "cpu_Mqps,gpu_Mqps,gpu_cached_Mqps,us_per_block_cached";

int main(int argc, char **argv)
{
  int nblocks = 64, nx = 16, nq = 2000, reps = 3, builder = 1;
  int cpuOnly = 0, verbose = 0, eltype = EL_HEX, skipDedup = 0;
  double qfrac = 1.0, perturb = 0.15;
  uint64_t seed = 12345;

  for (int a = 1; a < argc; a++) {
    if      (!strcmp(argv[a],"--nblocks") && a+1 < argc) nblocks = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--nx")      && a+1 < argc) nx = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--nq")      && a+1 < argc) nq = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--qfrac")   && a+1 < argc) qfrac = atof(argv[++a]);
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
    else if (!strcmp(argv[a],"--cpu-only"))   cpuOnly = 1;
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
    double org[3] = { (double)(b % side),
                      (double)((b / side) % side),
                      (double)(b / (side*side)) };
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
  }
  double tsetup = wtime() - tb;

  const long cellsPer  = blk[0].mesh.ncells;
  const long cellsTot  = cellsPer * nblocks;
  const long nqTot     = (long)nq * nblocks;

  if (verbose)
    printf("# %d blocks x %d^3 (%ld cells, %d queries) each = %ld cells, "
           "%ld queries | setup %.2f s\n",
           nblocks, nx, cellsPer, nq, cellsTot, nqTot, tsetup);

  /* ---------------- CPU: loop over all blocks ---------------- */
  SEARCHTIMERS cpuBest; memset(&cpuBest, 0, sizeof(cpuBest));
  double cpuBestTotal = 1e300;
  int donorsCpu = 0;
  for (int r = 0; r < reps + 1; r++) {                 /* r==0 warm-up */
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

  /* ---------------- GPU: loop over all blocks, one context ---------------- */
  SEARCHTIMERS gpuBest; memset(&gpuBest, 0, sizeof(gpuBest));
  double gpuBestTotal = 1e300, gpuBestCached = 1e300;
  int gpuOk = 0, donorsGpu = -1, mismatch = 0, bad = 0;

  if (!cpuOnly) {
    /* warm-up: creates every block's context state, uploads, builds */
    int rc = 0;
    for (int b = 0; b < nblocks; b++) {
      attachQueries(*blk[b].mb, blk[b].xs, nq);
      rc |= blk[b].mb->search_gpu();
    }
    gpuOk = (rc == 0);

    if (gpuOk) {
      /* (a) moving mesh: re-upload + rebuild every block, every call */
      for (int r = 0; r < reps; r++) {
        SEARCHTIMERS acc; memset(&acc, 0, sizeof(acc));
        double t0 = wtime();
        for (int b = 0; b < nblocks; b++) {
          attachQueries(*blk[b].mb, blk[b].xs, nq);
          blk[b].mb->gpuMeshDirty = 1;
          blk[b].mb->search_gpu();
          const SEARCHTIMERS &t = blk[b].mb->searchTimers;
          acc.transfer += t.transfer; acc.build += t.build;
          acc.dedup    += t.dedup;    acc.query += t.query;
        }
        acc.total = wtime() - t0;
        if (acc.total < gpuBestTotal) { gpuBestTotal = acc.total; gpuBest = acc; }
      }
      /* (b) static mesh: everything resident from the previous call */
      for (int r = 0; r < reps; r++) {
        double t0 = wtime();
        for (int b = 0; b < nblocks; b++) {
          attachQueries(*blk[b].mb, blk[b].xs, nq);
          blk[b].mb->search_gpu();
        }
        gpuBestCached = std::min(gpuBestCached, wtime() - t0);
      }

      donorsGpu = 0;
      for (int b = 0; b < nblocks; b++) {
        donorsGpu += blk[b].mb->donorCount;
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

  const double cpuMq    = nqTot / cpuBestTotal / 1e6;
  const double gpuMq    = gpuOk ? nqTot / gpuBestTotal  / 1e6 : 0.0;
  const double gpuCacMq = gpuOk ? nqTot / gpuBestCached / 1e6 : 0.0;
  const double usBlk    = gpuOk ? gpuBestCached*1e6/nblocks : 0.0;

  if (verbose) {
    printf("# CPU  %8.4f s  (filter %6.4f  ADT %6.4f  dedup %6.4f  walk %6.4f)"
           "  donors %d  %.3f M q/s\n",
           cpuBestTotal, cpuBest.filter, cpuBest.build, cpuBest.dedup,
           cpuBest.query, donorsCpu, cpuMq);
    if (gpuOk) {
      printf("# GPU  %8.4f s  (xfer   %6.4f  BVH %6.4f  dedup %6.4f  walk %6.4f)"
             "  donors %d  %.3f M q/s\n",
             gpuBestTotal, gpuBest.transfer, gpuBest.build, gpuBest.dedup,
             gpuBest.query, donorsGpu, gpuMq);
      printf("# GPU cached %8.4f s  %.3f M q/s  %.1f us/block  |  "
             "speedup %.2fx / %.2fx cached  |  mismatches %d (%d real)\n",
             gpuBestCached, gpuCacMq, usBlk,
             cpuBestTotal/gpuBestTotal, cpuBestTotal/gpuBestCached,
             mismatch, bad);
    }
  }

  printf("%d,%d,%s,%ld,%ld,%d,%ld,%d,"
         "%.6f,%.6f,%.6f,%.6f,%.6f,"
         "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
         "%d,%d,%d,%d,"
         "%.3f,%.3f,"
         "%.4f,%.4f,%.4f,%.2f\n",
         nblocks, nx, elTypeName(eltype), cellsPer, cellsTot, nq, nqTot, builder,
         cpuBestTotal, cpuBest.filter, cpuBest.build, cpuBest.dedup, cpuBest.query,
         gpuOk ? gpuBestTotal : 0.0, gpuOk ? gpuBest.transfer : 0.0,
         gpuOk ? gpuBest.build : 0.0, gpuOk ? gpuBest.dedup : 0.0,
         gpuOk ? gpuBest.query : 0.0, gpuOk ? gpuBestCached : 0.0,
         donorsCpu, donorsGpu, gpuOk ? mismatch : -1, gpuOk ? bad : -1,
         gpuOk ? cpuBestTotal/gpuBestTotal  : 0.0,
         gpuOk ? cpuBestTotal/gpuBestCached : 0.0,
         cpuMq, gpuMq, gpuCacMq, usBlk);
  fflush(stdout);

  for (int b = 0; b < nblocks; b++) { delete blk[b].mb; freeMesh(blk[b].mesh); }
  return 0;
}
