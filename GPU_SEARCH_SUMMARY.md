# Offloading `MeshBlock::search` to GPU with cuBQL

Branch `beryan/cubql_search`. All work is additive: without `-DTIOGA_ENABLE_CUDA=ON`
the build and behaviour of TIOGA are unchanged.

| file | what |
|---|---|
| `src/searchGPU.cu` | `MeshBlock::search_gpu()` — the cuBQL implementation |
| `src/searchGPU_stub.C` | non-CUDA fallback, returns `-1` |
| `src/search.C`, `src/codetypes.h` | `SEARCHTIMERS` phase breakdown, filled by both backends |
| `src/MeshBlock.{h,C}` | declaration, tuning knobs, device-state lifetime |
| `CMakeLists.txt`, `src/CMakeLists.txt` | `TIOGA_ENABLE_CUDA`, `TIOGA_CUBQL_DIR` |
| `benchmark/` | standalone harness, sweep scripts, CSV results, detailed README |

```
cmake -DTIOGA_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 \
      -DTIOGA_CUBQL_DIR=/scratch/software/cuBQL -DTIOGA_HAS_NODEGID=off ..
```

---

## The benchmark problem

`benchmark/bench_search.C` builds one mesh block and one cloud of query points —
i.e. one partition's worth of the donor search that `tioga::performConnectivity`
drives, isolated from MPI and from everything downstream.

- **Mesh.** `nx^3` cells over the unit cube. Interior nodes are displaced by a
  smooth field, `x += 0.15·dx·sin(πX)sin(πY)sin(πZ)·sin(3πY)` and cyclic
  permutations, so cells are **non-affine**: the containment test must run the
  full Newton iteration for the trilinear map, as on a real body-fitted grid.
  Sizes swept: 4,096 → 1,404,928 cells.
- **Element types.** `--eltype hex | prism | tet | mixed`. Each is carved out of
  the *same* logical hexes (2 prisms, or 6 tets by Kuhn subdivision, per hex; the
  mixed case alternates hex and prism on a checkerboard), so all four describe
  identical geometry and stay conforming. This exercises the `nvert = 8/6/4`
  containment paths and `ntypes > 1`.
- **Boundary conditions.** Outer-boundary nodes are declared overset boundary
  nodes, so `tagBoundary()` marks the outer cell layers with `cellRes=BIGVALUE`
  and the real donor-rejection logic is exercised.
- **Query points.** Uniform random inside a box of edge `qfrac` centred in the
  block. `qfrac=1` covers the whole block; small `qfrac` is the clustered cloud
  the host OBB pre-filter was written for; `qfrac>1` puts points outside the
  block so the "no donor found" path (`donorId = -1`) is exercised too.
  Sizes swept: 1,000 → 4,000,000 points. `--dup` injects exact duplicate points.
- **Not covered:** `uniform_hex` blocks (host direct-indexing path is already
  optimal), `ihigh != 0` (host callbacks), multi-block/multi-rank.

Machine: 1× GH200 120GB (sm_90) + 72-core Grace (Neoverse V2), CUDA 12.9,
nvhpc 25.7, `-O3`. **All CPU numbers are one core** unless stated otherwise —
`search()` has no OpenMP or threading anywhere.

---

## Algorithmic differences

| | `search()` — host ADT | `search_gpu()` — cuBQL |
|---|---|---|
| broad-phase input | cells overlapping the query-cloud OBB | **all** cells in the block |
| pre-filter | full pass over every cell projecting all vertices into the OBB frame; survivors pushed onto a LIFO list embedded in `icell[]`; second pass over their vertices to build world AABBs | none |
| acceleration structure | 6-D ADT over `(xmin,ymin,zmin,xmax,ymax,zmax)` | 3-D BVH over ordinary cell AABBs |
| structure lifetime | **rebuilt on every call** | built once, cached until `gpuMeshDirty` |
| build | serial median-split recursion (`buildADTrecursion`) | parallel GPU builder |
| traversal | recursive, one unique query point at a time | one thread per query point, iterative, 64-entry stack |
| duplicate points | `uniquenodes_octree`, then duplicates copy the representative's donor | same host call; the GPU searches every point anyway, so no gather is needed |
| exact containment | `checkContainment` → `computeNodalWeights` | device port of the same code |
| precision | double throughout | double throughout, **except** BVH bounds |

Three changes matter, in order of impact:

1. **The 6-D ADT is replaced by a cached 3-D BVH.** The ADT rebuild is what makes
   the host cost superlinear — at 1.4M cells it is 76 % of `search()`. Encoding a
   box as a point in 6-D was a way to make a 1-D binary tree answer a box query;
   a BVH does that natively in 3-D. And because the structure is cached, a static
   mesh pays for it once instead of every timestep.

2. **The query-cloud OBB pre-filter is dropped.** It costs a full pass over every
   cell — the same order as just building the BVH over every cell. Measured below:
   even in its best case it never pays back on the GPU path.

3. **The LIFO candidate list is gone.** `icell[p]=iptr; iptr=p` is an inherently
   serial linked list. Nothing replaces it: cell IDs are implicit in the BVH.

**Float bounds are safe.** Only the BVH bounds are single precision, and they are
rounded strictly outward with `nextafterf`, so a float box always contains the
double box it came from and the broad phase cannot produce a false negative. All
geometry — Newton iteration, 3×3 eliminations, the `[-TOL, 1+TOL]` weight test —
is double, from a line-by-line port of `math.c` and `checkContainment.C`.

**Donor selection semantics are preserved.** The device visitor returns
`CUBQL_TERMINATE_TRAVERSAL` on the first cleanly-accepted donor, matching the
ADT's early return, and keeps a boundary hit on a `cellRes==BIGVALUE` cell as a
tentative answer while continuing — exactly what `searchIntersections` does.
`gpuEarlyExit=0` instead scans all candidates and keeps the lowest cell ID, which
is order-independent and therefore deterministic.

### cuBQL features used

- `cuBQL::bvh3f` (`BinaryBVH<float,3>`) as the acceleration structure.
- `cuBQL::gpuBuilder` (adaptive spatial median, default) plus
  `cuBQL::cuda::radixBuilder`, `rebinRadixBuilder` and `sahBuilder`, all selectable
  at runtime via `gpuBuilderType`; `BuildConfig::makeLeafThreshold` via `gpuLeafSize`.
- `cuBQL::fixedBoxQuery::forEachPrim` — the traversal template. The exact
  containment test is the device lambda passed to it, so no candidate list is ever
  materialised in memory; a candidate is tested the moment it is found.
- `cuBQL::box3f` / `cuBQL::vec3f`, `cuBQL::cuda::free`, stream-aware build.
- **Header-only integration.** `searchGPU.cu` defines
  `CUBQL_GPU_BUILDER_IMPLEMENTATION` and instantiates the `float3` builders it
  needs, so no cuBQL library is built, installed or linked — one include path.
- Not used but relevant: `cuBQL::cuda::refit` (refit bounds instead of rebuilding
  for a deforming mesh) and `WideBVH`.

---

## Correctness

Every benchmark run cross-checks the GPU against the CPU **donor by donor** over
the whole query cloud, and classifies any difference:

- The GPU's donor is re-tested independently against `computeNodalWeights` in the
  harness (not via either search path). A difference is **benign** only if *both*
  donors genuinely contain the point — i.e. the point lies on a face shared by two
  cells and the two traversal orders picked different ones. Anything else is
  reported as a **real** mismatch with coordinates.
- Coverage: all four element mixes, `qfrac` from 0.05 to 1.4 (so both the
  "found" and the `donorId = -1` paths), 4k → 1.4M cells, 1k → 4M query points,
  all four builders, leaf sizes 4/8/16.

**Result: no case where the GPU returned a wrong donor**, across ~100
configurations. Two differences turned up:

- 1 in 50,000 (mixed hex/prism): both donors valid, point on a shared face.
- 1 in 4,000,000 (13,824-cell hex): **the CPU returned `-1` where the GPU
  returned a cell the independent check confirms contains the point.** A miss
  in `search()`, not in the port. Flagged, not investigated.

Separately, the full TIOGA regression (`case/run.sh`, sphere + background grid,
16 partitions on 8 ranks) was run before and after the change: interpolation
error statistics are byte-identical, ~1e-16 on every rank.

---

## Timing data

### Shmoo: throughput vs problem size (one CPU core, one GPU)

Millions of query points located per second, one CPU core vs one GH200. The
shared host duplicate-point pass is excluded on both sides (it is identical in
both backends and is not part of the search). Columns:

- `cpu_search` — OBB filter + ADT build + traversal, i.e. all of `search()`
- `gpu_search` — H2D upload + BVH build + traversal (mesh treated as moving)
- `gpu_cached` — mesh and BVH resident on the device (static mesh)
- `cpu_walk` / `gpu_walk` — traversal + exact containment only
- `bad` — GPU donors that do not contain the point (see Correctness)

Raw data: `benchmark/results/shmoo_search.csv` (as measured) and
`benchmark/results/shmoo_throughput.csv` (tidy long form, plots directly).
`benchmark/results/shmoo.csv` is the same grid *including* the host dedup.

| cells | queries | cpu_search | gpu_search | gpu_cached | cpu_walk | gpu_walk | bad |
|---|---|---|---|---|---|---|---|
| 4,096 | 10,000 | 1.045 | 11.9 | 122.0 | 1.391 | 303.0 | 0 |
| 4,096 | 100,000 | 1.322 | 88.8 | 300.3 | 1.443 | 1333.3 | 0 |
| 4,096 | 1,000,000 | 1.350 | 256.1 | 351.7 | 1.439 | 2493.8 | 0 |
| 4,096 | 4,000,000 | 1.359 | 298.6 | 328.7 | 1.446 | 2728.5 | 0 |
| 13,824 | 10,000 | 0.580 | 9.6 | 108.7 | 1.058 | 222.2 | 0 |
| 13,824 | 100,000 | 0.955 | 72.3 | 279.3 | 1.074 | 1075.3 | 0 |
| 13,824 | 1,000,000 | 1.015 | 235.2 | 331.2 | 1.071 | 1785.7 | 0 |
| 13,824 | 4,000,000 | 1.023 | 279.6 | 313.1 | 1.072 | 1874.4 | 1 |
| 32,768 | 10,000 | 0.341 | 9.1 | 122.0 | 0.921 | 285.7 | 0 |
| 32,768 | 100,000 | 0.772 | 59.2 | 297.6 | 0.932 | 1298.7 | 0 |
| 32,768 | 1,000,000 | 0.878 | 230.4 | 342.9 | 0.929 | 2320.2 | 0 |
| 32,768 | 4,000,000 | 0.899 | 282.1 | 325.5 | 0.940 | 2472.2 | 0 |
| 110,592 | 10,000 | 0.122 | 4.8 | 106.4 | 0.728 | 208.3 | 0 |
| 110,592 | 100,000 | 0.493 | 42.6 | 280.1 | 0.765 | 961.5 | 0 |
| 110,592 | 1,000,000 | 0.706 | 187.9 | 324.9 | 0.767 | 1644.7 | 0 |
| 110,592 | 4,000,000 | 0.731 | 256.0 | 307.1 | 0.765 | 1721.2 | 0 |
| 262,144 | 10,000 | 0.055 | 2.8 | 119.0 | 0.687 | 250.0 | 0 |
| 262,144 | 100,000 | 0.319 | 28.5 | 283.3 | 0.710 | 1123.6 | 0 |
| 262,144 | 1,000,000 | 0.619 | 150.2 | 335.6 | 0.713 | 2028.4 | 0 |
| 262,144 | 4,000,000 | 0.674 | 242.1 | 320.0 | 0.715 | 2154.0 | 0 |
| 512,000 | 10,000 | 0.028 | 1.9 | 106.4 | 0.555 | 175.4 | 0 |
| 512,000 | 100,000 | 0.179 | 18.3 | 261.8 | 0.477 | 763.4 | 0 |
| 512,000 | 1,000,000 | 0.475 | 117.6 | 304.0 | 0.583 | 1201.9 | 0 |
| 512,000 | 4,000,000 | 0.541 | 204.9 | 290.3 | 0.582 | 1277.1 | 0 |
| 884,736 | 10,000 | 0.016 | 1.3 | 108.7 | 0.472 | 181.8 | 0 |
| 884,736 | 100,000 | 0.120 | 12.8 | 253.2 | 0.510 | 746.3 | 0 |
| 884,736 | 1,000,000 | 0.385 | 90.0 | 302.1 | 0.520 | 1179.2 | 0 |
| 884,736 | 4,000,000 | 0.354 | 181.2 | 288.2 | 0.381 | 1241.1 | 0 |
| 1,404,928 | 10,000 | 0.010 | 0.9 | 108.7 | 0.385 | 169.5 | 0 |
| 1,404,928 | 100,000 | 0.073 | 8.9 | 251.9 | 0.299 | 757.6 | 0 |
| 1,404,928 | 1,000,000 | 0.292 | 66.9 | 298.2 | 0.426 | 1126.1 | 0 |
| 1,404,928 | 4,000,000 | 0.294 | 153.0 | 285.2 | 0.322 | 1182.7 | 0 |

Reading the table: **`cpu_search` collapses by 100× across the grid** (1.36 →
0.010 M/s) because the ADT is rebuilt every call, so host cost is driven by cell
count while throughput is measured per query. **`gpu_cached` is essentially flat
at ~250-350 M/s** over a 340× range of mesh sizes — with the tree resident, cost
depends only on the query count. `gpu_search` sits between the two: it carries
the H2D upload and BVH rebuild, so it degrades with cell count but far more
slowly than the host. Both GPU columns rise steeply with query count at small
`nq`, where fixed per-call overhead dominates and the GPU is starved.


#### Scaling the CPU column to a socket

The table is one core. Parallel efficiency across 72 cores is **not** constant --
measured with N pinned concurrent processes, `--skip-dedup`:

| block | 1 core | 72 cores | speedup | efficiency |
|---|---|---|---|---|
| 4,096 cells / 1M queries | 1.35 M q/s | 90.8 M q/s | 67.2x | 93 % |
| 262,144 cells / 500k queries | 0.36 | 21.0 | 58.2x | 81 % |
| 1,404,928 cells / 1M queries | 0.28 | 12.3 | 44.0x | 61 % |

Small blocks fit in cache and scale nearly perfectly; a 1.4M-cell block is
memory-bound and loses a third of its scaling. This measures 72 *copies* of the
same mesh, so footprint scales with rank count -- the pessimistic end. A real
MPI run partitions the mesh, so each of 72 ranks holds ~1/72 the cells, which is
the high-efficiency regime, and the ADT build shrinks superlinearly with
per-rank cell count. Use ~90 % for small per-rank blocks, ~60 % for million-cell
blocks.

### Where the time goes

Two representative cases, single core / single GPU:

**Large mesh, modest query cloud** — 1,404,928 cells, 100,000 queries. The
GPU's worst case in the sweep: the tree is expensive to build and there is little
query work to amortise it against.

| phase | CPU | GPU |
|---|---|---|
| OBB filter over all cells | 44.5 ms | — |
| build ADT / build BVH | 985.7 ms | 3.00 ms |
| H2D mesh + queries | — | 5.51 ms |
| dedup (*shared host code*) | 24.3 ms | 24.3 ms |
| traversal + containment | 241.4 ms | 0.14 ms |
| **total** | **1295.9 ms** | **33.1 ms** (39×) |

**Query-dominated** — 262,144 cells, 4,000,000 queries. Closer to a real
overset exchange where many receptors hit one block.

| phase | CPU | GPU |
|---|---|---|
| OBB filter | 176.8 ms | — |
| build ADT / BVH | 159.0 ms | 1.43 ms |
| dedup (*shared*) | 1338.0 ms | 1391.9 ms |
| traversal + containment | 5770.0 ms | 1.89 ms (**3055×**) |
| **total** | **7443.8 ms** | **1406.6 ms** (5.3×) |

The second table is the important one: **the search itself is 3000× faster, but
the shared host dedup (`uniquenodes_octree`, untouched by this work) is now 99 %
of the GPU total.** With it removed (`--skip-dedup`, which is legitimate because
the GPU searches every point regardless), the same case is 466×.

### Does dropping the OBB pre-filter hurt?

884,736 cells, 500,000 queries, `--skip-dedup`. Small `qfrac` is exactly the
regime the filter was written for.

| qfrac | cells surviving filter | CPU filter | CPU ADT build | CPU total | GPU total | GPU cached |
|---|---|---|---|---|---|---|
| 0.05 | 742 | 38.1 | 0.3 | 424.4 | 7.41 | 1.44 |
| 0.25 | 59,338 | 39.1 | 31.1 | 744.1 | 7.49 | 1.50 |
| 1.00 | 884,736 | 50.8 | 598.9 | 1765.2 | 7.57 | 1.65 |

The filter works — at `qfrac=0.05` it cuts the ADT build from 599 ms to 0.3 ms —
but the *scan* that produces it costs 38 ms, five times the GPU's entire
end-to-end search over all 884,736 cells and twenty-six times its cached time.
**No regime in this sweep justifies keeping it on the GPU path.**

### cuBQL builder choice (512,000 cells)

| builder | leaf | BVH build (ms) | traversal (ms) |
|---|---|---|---|
| spatial median (default) | default | 1.71 | 0.132 |
| **radix** | default | **0.78** | **0.101** |
| rebin radix | default | 0.76 | 0.115 |
| SAH | default | 17.22 | 0.101 |
| spatial median | 8 | 1.23 | 0.324 |
| spatial median | 16 | 1.12 | 0.534 |

`radixBuilder` is the right default here: 2.2× faster build, marginally faster
traversal. SAH costs 10× more to build and buys nothing. Larger leaves trade
build for traversal badly — the containment test is too expensive to brute-force.

### Element type (200,000 queries, GPU traversal in ms)

hex 110,592 cells / 0.159 · prism 221,184 / 0.194 · tet 663,552 / 0.226 ·
mixed hex+prism 165,888 / 0.207. Warp divergence from mixed element types is
mild — the mixed mesh costs 30 % more GPU traversal than pure hex at 1.5× the
cells — so bucketing candidates by element type is not worth it here.

### Socket-level: 72 Grace cores vs 1 GH200

The single-core ratios above overstate the practical win. TIOGA parallelises by
MPI rank over partitions, so the honest comparison is a whole 72-core Grace
socket against the one GPU it is paired with. Measured with N pinned concurrent
processes, 1,404,928 cells / 100,000 queries each, in block-searches per second:

| ranks | CPU | GPU, mesh rebuilt each call | GPU, BVH cached |
|---|---|---|---|
| 1 | 0.7 | 29.4 | 40.1 |
| 36 | 22.5 | 149.8 | 1210.9 |
| 72 | 41.1 | 104.3 | 1938.5 |
| 72 + CUDA MPS | 41.4 | **269.9** | **2759.7** |
| **vs 72-core socket** | 1.0× | **6.5×** | **67×** |

So on this configuration the win is **6.5× if the mesh moves every step and is
re-uploaded, 67× if the mesh stays resident on the device.** Two caveats:

- The rebuild column is bound by data movement, not compute: at 72 ranks with MPS
  the phase breakdown is xfer 218 ms, BVH build 33 ms, traversal 17 ms. That is
  72 × ~80 MB moving at ~26 GB/s, far below NVLink-C2C, so it is
  context-scheduling overhead from 72 separate CUDA contexts — an artifact of how
  I measured it, not a property of the algorithm. Without MPS throughput actually
  *drops* from 36→72 ranks as the contexts thrash.
- One CUDA context per rank is the wrong architecture anyway. One rank owning the
  GPU with many blocks in a single context is right, and I have **not** measured it.

---

## Conclusions

1. **The search itself is comprehensively GPU-friendly.** Traversal plus exact
   containment reaches ~2.1 G point-locations/s, against ~0.36 M/s per core and
   ~21 M/s for a full 72-core socket. Correctness is exact.
2. **Caching the acceleration structure is worth more than the traversal speedup**
   for a static or slowly-deforming mesh: 6.5× → 67× at socket scale. Rebuilding
   an ADT every call is the single largest structural cost in the host version.
3. **`uniquenodes_octree` is now the bottleneck**, at up to 99 % of GPU total. It
   is host code shared by both backends and untouched here. The GPU search does
   not need it — it searches every point regardless; it exists only to produce
   `xtag`/`res_search` for the downstream donor exchange (`bookKeeping.C`).
   Cheapest fix: build with `TIOGA_HAS_NODEGID=on` (TIOGA's default, and what
   nalu-wind uses), where dedup keys on 64-bit global node IDs instead of doing
   TOL-proximity matching in an octree — that maps directly onto a sort-by-key
   plus segmented scan.
4. **Drop the OBB pre-filter** on the GPU path. It never pays back.
5. **Use `radixBuilder`** with the default leaf size.

### Next steps, in order of expected value

1. Port the unique-node pass (or switch to `TIOGA_HAS_NODEGID`). Everything else
   is Amdahl-limited until this is done.
2. Keep the mesh device-resident across timesteps and use `cuBQL::cuda::refit`
   for deformation instead of re-uploading and rebuilding.
3. One BVH per rank aggregating many mesh blocks, one CUDA context per GPU. This
   is the configuration a production run would actually use and it is unmeasured.
4. Temporal coherence: test the previous donor and its neighbours before falling
   back to the BVH. For moving overset meshes this can make most searches O(1).

Full methodology, all raw CSVs and the complete tables are in
`benchmark/README.md` and `benchmark/results/`.
