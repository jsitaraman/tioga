# Offloading `MeshBlock::search` to GPU with cuBQL

Branch `beryan/cubql_search`. Additive: without `-DTIOGA_ENABLE_CUDA=ON` the
build and behaviour of TIOGA are unchanged.

```
cmake -DTIOGA_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 \
      -DTIOGA_CUBQL_DIR=/scratch/software/cuBQL -DTIOGA_HAS_NODEGID=off ..
```

| file | what |
|---|---|
| `src/searchGPU.cu` | `search_gpu()` (per block) and `search_gpu_batch()` (one BVH per rank) |
| `src/searchGPU_stub.C` | non-CUDA fallback, returns `-1` |
| `src/search.C`, `src/codetypes.h` | `SEARCHTIMERS` phase breakdown, filled by both backends |
| `benchmark/` | harnesses, sweep scripts, CSV results, methodology notes |

cuBQL is header-only here — `searchGPU.cu` defines
`CUBQL_GPU_BUILDER_IMPLEMENTATION` and instantiates the `float3` builders it
needs. No cuBQL library is built or linked.

Machine: 1× GH200 120GB (sm_90) + 72-core Grace, CUDA 12.9, nvhpc 25.7, `-O3`.

---

## 1. Algorithm

| | `search()` — host ADT | `search_gpu_batch()` — cuBQL |
|---|---|---|
| broad-phase input | cells overlapping the query-cloud OBB | all cells of **every block the rank owns** |
| pre-filter | full pass over every cell, LIFO compaction, second pass for world AABBs | none |
| structure | 6-D ADT over `(xmin..zmax)`, **rebuilt every call** | one 3-D BVH per rank, cached; refit when the mesh moves |
| build | serial median-split recursion | parallel GPU builder |
| traversal | recursive, one unique point at a time | one thread per point, `fixedBoxQuery::forEachPrim` |
| containment | `checkContainment` → `computeNodalWeights` | line-by-line device port, same `TOL`, double precision |

Four changes carry the result:

1. **One BVH per rank, not per block.** At 16³ blocks a per-block BVH is pure
   launch latency. Batching is 5.5–7.4× faster (§4).
2. **The structure is cached.** The ADT is rebuilt on every call; the BVH is
   rebuilt only when connectivity changes, refit when coordinates move.
3. **The OBB pre-filter is dropped.** It costs a full pass over every cell — the
   same order as building the BVH over every cell. Measured: even at
   `qfrac=0.05`, where it cuts the ADT build from 599 ms to 0.3 ms, the *scan*
   costs 38 ms against 7.6 ms for the GPU's entire search. It never pays back.
4. **The serial LIFO candidate list is gone.** Cell IDs are implicit in the BVH.

**Donor semantics are preserved.** The visitor terminates on the first cleanly
accepted donor, matching the ADT's early return, and holds a boundary hit on a
`cellRes==BIGVALUE` cell as tentative while continuing — as
`searchIntersections` does. In the batched path each cell carries its source
block and each query its owning block, and cross-block candidates are rejected,
so overlapping blocks cannot steal each other's receptors.

**Precision.** Only BVH bounds are single precision, rounded strictly outward
with `nextafterf`, so the broad phase cannot produce a false negative. All
geometry — Newton iteration, 3×3 eliminations, the `[-TOL, 1+TOL]` test — is
double.

### cuBQL features used

- `cuBQL::bvh3f` (`BinaryBVH<float,3>`); `box3f`/`vec3f`; `cuda::free`.
- `gpuBuilder` (spatial median) plus `cuda::radixBuilder`, `rebinRadixBuilder`,
  `sahBuilder`, selectable at runtime; `BuildConfig::makeLeafThreshold`.
- `cuBQL::cuda::refit` for the moving-mesh path.
- `cuBQL::fixedBoxQuery::forEachPrim` — the containment test *is* the device
  lambda passed to it, so no candidate list is ever materialised.

---

## 2. The benchmark problem

`benchmark/bench_multiblock.C` (production shape) and `bench_search.C` (single
block). A rank owns N mesh blocks; each is `nx^3` cells with a cloud of receptor
points, isolated from MPI and from everything downstream.

- **Non-affine cells.** Interior nodes are displaced by a smooth field, so the
  containment test runs the full Newton iteration for the trilinear map.
- **Element types.** `hex | prism | tet | mixed`, each carved from the *same*
  logical hexes (2 prisms, or 6 tets by Kuhn subdivision), so all four are the
  same geometry and stay conforming. Exercises `nvert = 8/6/4` and `ntypes > 1`.
- **Boundary conditions.** Outer nodes are declared overset boundary nodes, so
  `tagBoundary()` marks outer cell layers `cellRes=BIGVALUE` and the real
  donor-rejection logic runs.
- **Query cloud.** `qfrac` sets its extent; `qfrac>1` puts points outside the
  block so the `donorId = -1` path runs. `--overlap` makes neighbouring blocks
  overlap in space. `--move` displaces the mesh so refit runs on coordinates the
  tree was not built for.
- **Production point:** 16³ active cells per block, ~2,000 receptor points each,
  many blocks per rank, one CUDA context per GPU.

---

## 3. Correctness

Every run cross-checks GPU against CPU **donor by donor**. The GPU's donor is
re-verified independently against `computeNodalWeights` in the harness. A
difference counts as *benign* only if both donors genuinely contain the point
(shared face, both valid); anything else is reported as real, with coordinates.

Coverage: 4k → 1.4M cells, 1k → 4M points, 1 → 4,608 blocks, all four element
mixes, `qfrac` 0.05–1.4, block overlap 0–0.6, displaced coordinates with refit
and rebuild, all four builders, leaf sizes 4/8/16.

**No case was found where the GPU returned a wrong donor.** Two differences
arose:

- 1 in 50,000 (mixed hex/prism): both donors valid, point on a shared face.
- 1 in 4,000,000 (13,824-cell hex, reproducible): **the CPU returned `-1` where
  the GPU returned a cell that independently checks out.** A miss in `search()`,
  not in the port. Flagged, not investigated.

The full TIOGA regression (`case/run.sh`, sphere + background, 16 partitions on
8 ranks) gives byte-identical interpolation error before and after, ~1e-16.

---

## 4. Results — production configuration

4,608 blocks of 16³ per rank = 18,874,368 cells, 9,216,000 receptor points. One
rank, one CUDA context. CPU baseline is the same total work as 72 ranks × 64
blocks on a full Grace socket, barrier-synchronised.

### Compute versus compute

Host/device transfer and the shared host dedup excluded on both sides.
CPU = OBB filter + ADT build + walk. GPU = cell-AABB kernel + BVH refit +
traversal kernel.

| | seconds | M pts/s | vs 1 core | vs 72-core socket |
|---|---|---|---|---|
| CPU compute, 1 core | 15.787 | 0.58 | 1.0× | — |
| CPU compute, 72-core aggregate | 0.224 | 41.1 | 70× | 1.0× |
| **GPU compute, moving mesh** | 0.0092 | **1007** | **1725×** | **24.5×** |
| GPU compute, static mesh | 0.0052 | 1757 | 3011× | 42.8× |

Neither column contains dedup: `cpu_compute` is `filter + build + query`, and
`dedup` is timed separately on both sides. One asymmetry does remain, though —
the CPU's *query* loop uses `xtag` to skip duplicate points, while the GPU
searches every point. The benchmark generates no duplicates, so as measured both
backends locate the same 9.216M points. At a real duplicate rate `d` the CPU's
walk shrinks as `(1-d)` while its ADT build (52 % of CPU compute) does not, so
the advantage moves only slowly: 24.5× at `d=0`, 22.4× at `d=0.2`, 19.3× at
`d=0.5`.

**No host round-trip inside TIOGA's compute chain** — cell-AABB kernel → refit →
traversal is one dependent chain on one stream. cuBQL's *builders* do round-trip
(`sm_builder` and `radixBuilder` D2H the build state and event-sync before every
level); `refit` does not, which is part of why it wins. Those costs are inside
the figures above.

### Moving mesh

Connectivity does not move, so `gpuCoordsDirty` re-sends only `x[]`, straight
into each block's device slice — no host staging, no re-walking connectivity.
`gpuRefit` (default) refits the tree instead of rebuilding it; always correct,
only query performance degrades as the mesh drifts.

| | xfer | BVH | walk | compute | total | M pts/s |
|---|---|---|---|---|---|---|
| full rebuild, topology re-sent | 0.349 | 0.0084 | 0.0052 | 0.0136 | 0.386 | 23.9 |
| coords only, rebuild | 0.069 | 0.0084 | 0.0052 | 0.0136 | 0.106 | 87.2 |
| **coords only, refit** | 0.069 | 0.0039 | 0.0052 | 0.0091 | **0.101** | **91.4** |
| static, fully cached | — | — | 0.0052 | 0.0052 | 0.061 | 150.9 |

### One BVH per rank vs one per block

16³ blocks, 2,000 points each, cached, dedup excluded:

| blocks | per-block | batched | gain | per-block µs/blk | batched µs/blk |
|---|---|---|---|---|---|
| 64 | 32.5 M pts/s | 239.2 | 7.4× | 61.5 | 8.4 |
| 256 | 30.8 | 192.0 | 6.2× | 64.9 | 10.4 |
| 1024 | 28.6 | 157.5 | 5.5× | 70.1 | 12.7 |

At 64 blocks the BVH build drops 25.1 ms → 0.6 ms and traversal 2.1 → 0.1 ms.

### End-to-end, as currently implemented

| | M pts/s | vs 72-core socket |
|---|---|---|
| GPU compute, moving | 1007 | 24.5× |
| + host pack/unpack between kernels | 230 | 5.6× |
| + H2D/D2H transfer | 91 | 2.2× |
| + host dedup | 5.3 | **0.14×** |

`hostwork` is gathering per-block query arrays into one buffer and scattering
`donorId` back — 3.4 ns/point against 1712 ns/point for the CPU search. It
exists only because TIOGA passes receptor points as per-block host arrays, and
disappears once they are device-resident.

---

## 5. The dedup decides the outcome

`uniquenodes_octree` — duplicate receptor-point detection, unchanged host code
shared by both backends — is 1.37 s of the GPU rank's time against 9 ms of
compute. It is the difference between a 24.5× win and a 6× loss.

It does three jobs. Only the first concerns the search:

1. **Skip redundant searches** (`donorId[i]=donorId[xtag[i]]`). Pure CPU
   optimisation; the GPU searches every point and duplicates cost nothing.
   Droppable.
2. **Max-merge `res_search`** across copies. `nodeRes` can be `BIGVALUE`
   ("mandatory receptor"), and the merge is how that propagates to every copy.
3. **Group-min in `resetCoincident()`** — if any copy of a coincident point
   survives, all copies are un-cancelled.

2 and 3 need the *grouping*, not the search. **Empirically, dropping the
grouping entirely changed nothing** on the standard case: forcing `xtag[i]=i`
gave identical interpolation error to 16 digits and byte-identical donor/fringe
counts on all 16 blocks. That is one topology and one metric, so verify per
configuration — but it suggests the grouping may be droppable outright.

If it is not, it is a sort + segmented reduce (max on resolution, min on cancel)
and belongs in `exchangeSearchData` where the points arrive, not in `search()`.
With `TIOGA_HAS_NODEGID` the key is an exact 64-bit integer, so it is a plain
sort-by-key with no tolerance matching.

---

## 6. Secondary findings

- **Builder:** `cuda::radixBuilder` is the right default — 0.78 ms build vs
  1.71 ms for spatial median at 512k cells, and marginally faster traversal.
  SAH costs 10× more to build and buys nothing. Larger leaves lose: the
  containment test is too expensive to brute-force a leaf.
- **Element type barely matters.** A tet mesh with 6× the cells costs 42 % more
  GPU traversal; mixed hex/prism costs 30 % more than pure hex at 1.5× the
  cells. Not enough to justify bucketing candidates by element type.
- **Single-block scaling:** with a cached BVH, GPU cost tracks query count,
  not mesh size — the host ADT is rebuilt every call, so CPU cost grows with
  both. See `benchmark/results/` from `sweep.sh` / `bench_search`.

---

## 7. Caveats

- `uniform_hex` blocks stay on the host direct-indexing path (already optimal);
  `ihigh != 0` returns `-2` (containment goes through host callbacks).
- The moving-mesh timings re-send coordinates and refit, but those coordinates
  are unchanged, so BVH quality decay under sustained motion is **not**
  captured. Refit correctness under displacement *is* tested (`--move`).
- The compute comparison is at 0 % duplicate receptor points. A real duplicate
  rate discounts the CPU walk but not the ADT build; see §4.
- CPU parallel efficiency is not constant: 93 % for a 4k-cell block, 61 % for a
  1.4M-cell block. The 72-core baseline here replicates the mesh per rank, the
  pessimistic end; a partitioned run puts each rank in the high-efficiency
  regime.
- **A withdrawn measurement:** earlier socket numbers from 72 ranks sharing the
  GPU via MPS (0.62×/5.7×, 6.5×/67×) were invalid and too favourable to the GPU.
  Throughput was total queries over the slowest rank's wall time, but ranks
  reached their GPU phase at different times and each found an idle GPU. The
  tell was an aggregate of 13,572 M pts/s against a ~2,100 M pts/s kernel
  ceiling. Replaced by the barrier-synchronised comparison in §4. (MPS also caps
  at 48 clients on this GPU.)

---

## 8. Next steps, in order

1. **Deal with the dedup.** Either drop the grouping (verify per configuration,
   as in §5) or move it off the critical path. Nothing else changes the outcome.
2. **Keep receptor points device-resident** so the host pack/unpack goes away —
   worth 4× on top (230 → 1007 M pts/s).
3. **Periodic rebuild.** Refit is cheap but tree quality decays under sustained
   motion; measure how often a rebuild is needed on real trajectories.
4. Temporal coherence — test the previous donor and its neighbours before
   falling back to the BVH. For moving overset meshes this can make most
   searches O(1).

Raw data in `benchmark/results/`, methodology in `benchmark/README.md`.
