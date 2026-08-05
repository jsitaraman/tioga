// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

/*
 * Micro-benchmark for MeshBlock::search() (ADT, host) versus
 * MeshBlock::search_gpu() (cuBQL BVH, device).
 *
 * A single unstructured mesh block of perturbed hexahedra stands in for one
 * overset partition, and a cloud of random query points stands in for the
 * receptor points another partition would have sent over. Both searches see
 * exactly the same mesh and the same points, and both fill donorId[].
 *
 * The mesh is deliberately non-affine (interior nodes are displaced by a
 * smooth field) so that the containment test has to run the full Newton
 * iteration for the trilinear map, as it would on a real body-fitted grid.
 * Outer-boundary nodes are declared overset boundary nodes, so tagBoundary()
 * marks the outer few cell layers with cellRes=BIGVALUE, exercising the same
 * donor-rejection logic as a real run.
 *
 * Output is one CSV row per configuration; see --header.
 */
#include "meshgen.h"

struct Result
{
  SEARCHTIMERS t;
  int          donors;
  std::vector<int> donorId;
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

static void capture(MeshBlock &mb, Result &r)
{
  r.t      = mb.searchTimers;
  r.donors = mb.donorCount;
  r.donorId.assign(mb.donorId, mb.donorId + mb.nsearch);
}

static void usage(const char *p)
{
  fprintf(stderr,
    "usage: %s [options]\n"
    "  --nx N          cells per direction (mesh is N^3 hexes)   [32]\n"
    "  --eltype T      hex | prism | tet | mixed                 [hex]\n"
    "  --nq N          number of query points                    [100000]\n"
    "  --qfrac F       query cloud edge as fraction of the block [1.0]\n"
    "                  (>1 puts some query points outside the block)\n"
    "  --dup F         fraction of query points that are exact duplicates [0]\n"
    "  --perturb F     interior node displacement, in cell sizes [0.15]\n"
    "  --reps N        timed repetitions (best is reported)      [3]\n"
    "  --seed N        RNG seed                                  [12345]\n"
    "  --builder N     0=spatial median 1=radix 2=rebin 3=SAH    [0]\n"
    "  --leaf N        cuBQL makeLeafThreshold (0=default)       [0]\n"
    "  --no-early-exit GPU scans all candidates, keeps lowest id\n"
    "  --skip-dedup    GPU skips the host duplicate-point pass\n"
    "  --cpu-only      skip the GPU run\n"
    "  --header        print the CSV header and exit\n"
    "  --verbose       print a human readable breakdown too\n", p);
}

static const char *CSV_HEADER =
  "nx,eltype,ncells,nnodes,nq,qfrac,dup,builder,leaf,"
  "cpu_total,cpu_filter,cpu_build,cpu_dedup,cpu_query,cpu_cand,"
  "gpu_total,gpu_xfer,gpu_build,gpu_dedup,gpu_query,gpu_total_cached,gpu_cand,"
  "donors_cpu,donors_gpu,mismatch,mismatch_real,"
  "speedup_total,speedup_total_cached,speedup_query";

/* No MPI_Init here: MeshBlock::search()/search_gpu() are purely local, and
   keeping the benchmark MPI-free lets it run without an mpirun launcher. */
int main(int argc, char **argv)
{
  int nx = 32, nq = 100000, reps = 3, builder = 0, leaf = 0;
  int cpuOnly = 0, verbose = 0, earlyExit = 1, eltype = EL_HEX, skipDedup = 0;
  double qfrac = 1.0, dupFrac = 0.0, perturb = 0.15;
  uint64_t seed = 12345;

  for (int a = 1; a < argc; a++) {
    if      (!strcmp(argv[a],"--nx")      && a+1 < argc) nx = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--eltype")  && a+1 < argc) {
      const char *e = argv[++a];
      if      (!strcmp(e,"hex"))   eltype = EL_HEX;
      else if (!strcmp(e,"prism")) eltype = EL_PRISM;
      else if (!strcmp(e,"tet"))   eltype = EL_TET;
      else if (!strcmp(e,"mixed")) eltype = EL_MIXED;
      else { fprintf(stderr,"unknown --eltype %s\n",e); return 1; }
    }
    else if (!strcmp(argv[a],"--nq")      && a+1 < argc) nq = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--qfrac")   && a+1 < argc) qfrac = atof(argv[++a]);
    else if (!strcmp(argv[a],"--dup")     && a+1 < argc) dupFrac = atof(argv[++a]);
    else if (!strcmp(argv[a],"--perturb") && a+1 < argc) perturb = atof(argv[++a]);
    else if (!strcmp(argv[a],"--reps")    && a+1 < argc) reps = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--seed")    && a+1 < argc) seed = strtoull(argv[++a],NULL,10);
    else if (!strcmp(argv[a],"--builder") && a+1 < argc) builder = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--leaf")    && a+1 < argc) leaf = atoi(argv[++a]);
    else if (!strcmp(argv[a],"--no-early-exit")) earlyExit = 0;
    else if (!strcmp(argv[a],"--skip-dedup")) skipDedup = 1;
    else if (!strcmp(argv[a],"--cpu-only")) cpuOnly = 1;
    else if (!strcmp(argv[a],"--verbose"))  verbose = 1;
    else if (!strcmp(argv[a],"--header")) { printf("%s\n",CSV_HEADER); return 0; }
    else { usage(argv[0]); return 1; }
  }

  Mesh m;
  double tb = wtime();
  buildMesh(m, nx, nx, nx, perturb, eltype);
  std::vector<double> xs;
  buildQueries(m, nq, qfrac, dupFrac, seed, xs);
  double tmesh = wtime() - tb;

  MeshBlock mb;
  mb.setData(1, m.nnodes, m.x, m.iblank, 0, m.nobc, NULL, m.obcnode,
             m.ntypes, m.nv, m.nc, m.vconn);
  tb = wtime();
  mb.preprocess();
  double tpre = wtime() - tb;

  if (verbose)
    printf("# mesh %d^3 %s = %d cells, %d nodes | build %.3fs preprocess %.3fs\n",
           nx, elTypeName(eltype), m.ncells, m.nnodes, tmesh, tpre);

  /* ---------------- CPU ---------------- */
  Result cpu; SEARCHTIMERS cpuBest; double cpuBestTotal = 1e300;
  memset(&cpuBest, 0, sizeof(cpuBest));
  for (int r = 0; r < reps + 1; r++) {          /* r==0 is a warm-up */
    attachQueries(mb, xs, nq);
    mb.search();
    if (r == 0) continue;
    if (mb.searchTimers.total < cpuBestTotal) {
      cpuBestTotal = mb.searchTimers.total;
      cpuBest = mb.searchTimers;
      capture(mb, cpu);
    }
  }

  /* ---------------- GPU ---------------- */
  Result gpu; SEARCHTIMERS gpuBest; double gpuBestTotal = 1e300;
  double gpuBestCached = 1e300;
  int gpuOk = 0, mismatch = 0, benign = 0, bad = 0;
  memset(&gpuBest, 0, sizeof(gpuBest));
  gpu.donors = -1;

  if (!cpuOnly) {
    mb.gpuBuilderType = builder;
    mb.gpuLeafSize    = leaf;
    mb.gpuEarlyExit   = earlyExit;
    mb.gpuSkipDedup   = skipDedup;

    /* warm-up: creates the context, uploads the mesh, builds the BVH */
    attachQueries(mb, xs, nq);
    int rc = mb.search_gpu();
    gpuOk = (rc == 0);

    if (gpuOk) {
      /* (a) full cost, mesh treated as moving: re-upload + rebuild each call */
      for (int r = 0; r < reps; r++) {
        attachQueries(mb, xs, nq);
        mb.gpuMeshDirty = 1;
        mb.search_gpu();
        if (mb.searchTimers.total < gpuBestTotal) {
          gpuBestTotal = mb.searchTimers.total;
          gpuBest = mb.searchTimers;
          capture(mb, gpu);
        }
      }
      /* (b) static mesh: BVH and mesh copy reused across calls */
      for (int r = 0; r < reps; r++) {
        attachQueries(mb, xs, nq);
        mb.search_gpu();
        gpuBestCached = std::min(gpuBestCached, mb.searchTimers.total);
      }

      /* A mismatch is benign when both answers are genuine donors: the
         point sits on a face shared by two cells and the two traversal
         orders picked different ones. Anything else is a real defect. */
      for (int i = 0; i < nq; i++) {
        if (cpu.donorId[i] == gpu.donorId[i]) continue;
        mismatch++;
        if (donorIsValid(m, cpu.donorId[i], &xs[3*i]) &&
            donorIsValid(m, gpu.donorId[i], &xs[3*i]))
          benign++;
        else if (bad < 5) {
          fprintf(stderr,
                  "# BAD point %d (%.17g %.17g %.17g): cpu donor %d (valid=%d) "
                  "gpu donor %d (valid=%d)\n",
                  i, xs[3*i], xs[3*i+1], xs[3*i+2],
                  cpu.donorId[i], (int)donorIsValid(m, cpu.donorId[i], &xs[3*i]),
                  gpu.donorId[i], (int)donorIsValid(m, gpu.donorId[i], &xs[3*i]));
          bad++;
        } else bad++;
      }
    }
  }

  if (verbose) {
    printf("# CPU  total %8.4f s  (filter %6.4f  buildADT %6.4f  dedup %6.4f  "
           "traverse %6.4f)  cand cells %d  donors %d\n",
           cpuBest.total, cpuBest.filter, cpuBest.build, cpuBest.dedup,
           cpuBest.query, cpuBest.candidates, cpu.donors);
    if (gpuOk)
      printf("# GPU  total %8.4f s  (xfer   %6.4f  buildBVH %6.4f  dedup %6.4f  "
             "traverse %6.4f)  cells %d  donors %d  cached total %8.4f s\n",
             gpuBest.total, gpuBest.transfer, gpuBest.build, gpuBest.dedup,
             gpuBest.query, gpuBest.candidates, gpu.donors, gpuBestCached);
    if (gpuOk)
      printf("# speedup total %6.2fx  cached %6.2fx  traversal-only %6.2fx  "
             "mismatches %d/%d (%d benign, %d real)\n",
             cpuBest.total/gpuBest.total, cpuBest.total/gpuBestCached,
             cpuBest.query/gpuBest.query, mismatch, nq, benign, bad);
  }

  printf("%d,%s,%d,%d,%d,%g,%g,%d,%d,"
         "%.6f,%.6f,%.6f,%.6f,%.6f,%d,"
         "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,"
         "%d,%d,%d,%d,"
         "%.3f,%.3f,%.3f\n",
         nx, elTypeName(eltype), m.ncells, m.nnodes, nq, qfrac, dupFrac,
         builder, leaf,
         cpuBest.total, cpuBest.filter, cpuBest.build, cpuBest.dedup,
         cpuBest.query, cpuBest.candidates,
         gpuOk ? gpuBest.total : 0.0, gpuOk ? gpuBest.transfer : 0.0,
         gpuOk ? gpuBest.build : 0.0, gpuOk ? gpuBest.dedup : 0.0,
         gpuOk ? gpuBest.query : 0.0, gpuOk ? gpuBestCached : 0.0,
         gpuOk ? gpuBest.candidates : 0,
         cpu.donors, gpu.donors, gpuOk ? mismatch : -1, gpuOk ? bad : -1,
         gpuOk ? cpuBest.total/gpuBest.total : 0.0,
         gpuOk ? cpuBest.total/gpuBestCached : 0.0,
         gpuOk ? cpuBest.query/gpuBest.query : 0.0);
  fflush(stdout);

  freeMesh(m);
  return 0;
}
