# GPU donor search (cuBQL) — what matters

Optional GPU path for `MeshBlock` donor location using
[cuBQL](https://github.com/NVIDIA/cuBQL). Host `search()` is unchanged
unless you call the new APIs.

The function that matters for production-shaped TIOGA is
**`MeshBlock::search_gpu_batch()`** — one BVH per rank over every block
that rank owns. `search_gpu()` (one BVH per block) is useful for
micro-benchmarks and as the control that shows why batching is required.

Longer write-up and full tables: `../GPU_SEARCH_SUMMARY.md`.
Raw CSVs: `results/`.

---

## Key takeaways

1. **Batch, don't per-block.** At 16³ cells/block, a BVH + launch per
   block is launch-bound. One BVH per rank is **5.5–7.4×** faster
   (`results/mb_batch_vs_perblock.csv`).
2. **Search compute is fast.** On a GH200, batched GPU compute
   (AABB + refit + walk) is ~**25×** a full Grace socket's host ADT
   compute for 4608×16³ blocks and ~2k points/block — when host dedup
   and transfers are excluded (`results/mb_compute.csv`).
3. **Host dedup decides end-to-end.** With today's
   `uniquenodes_octree` still on the host, batched GPU is
   **~0.14–0.17×** the same socket (`results/mb_production.csv`).
   The GPU searches every point; dedup is for `xtag` / `res_search`
   used later in donor exchange, not for finding donors.
4. **Not wired into TIOGA yet.** `tioga::performConnectivity()` still
   calls `search()` only. This directory exercises the library APIs.
5. **Correct on synthetic donors.** Donor-by-donor vs CPU across
   sweeps; no real GPU false donors. One case: CPU `-1`, GPU a cell
   that independently contains the point.

---

## What was done

| Piece | Role |
|---|---|
| `src/searchGPU.cu` | `search_gpu()`, `search_gpu_batch()`, device containment |
| `src/searchGPU_stub.C` | Non-CUDA build; APIs return `-1` |
| `SEARCHTIMERS` | Phase timing in both backends |
| `bench_multiblock.C` | Many blocks/rank; `--batch`, `--moving` |
| `bench_search.C` | Single-block microbench (flatter GPU) |
| `meshgen.h` | Shared synthetic mesh + query cloud |
| `sweep*.sh` | Reproduce CSV sweeps |

cuBQL is header-only here (`CUBQL_GPU_BUILDER_IMPLEMENTATION` in
`searchGPU.cu`). No separate cuBQL library link.

---

## Why it works (algorithm)

| | Host `search()` | GPU `search_gpu_batch()` |
|---|---|---|
| Broad phase | 6-D ADT on OBB-filtered cells | 3-D BVH on **all** cells of all blocks |
| Structure lifetime | Rebuilt every call | Cached; refit when coords move |
| Traversal | Serial, unique points | One thread per point |
| Containment | `checkContainment` | Same math on device, double |
| Multi-block | One `search()` per block | One launch; reject cross-block hits |

Four design choices carry the result:

- **One BVH per rank** — amortizes build and launch at small blocks.
- **Cache / refit** — ADT rebuild every call is the host scaling wall.
- **Drop query-cloud OBB pre-filter** — the scan costs as much as
  building a BVH over all cells; measured, it never helps the GPU.
- **Float BVH bounds, outward-rounded** — broad phase cannot false-
  negative; exact test stays double.

Donor semantics match the ADT early-exit and `BIGVALUE` fringe
rejection. Batch tags cells and queries by owning block so overlaps
cannot steal receptors.

---

## How to reproduce

```
mkdir build-gpu && cd build-gpu
cmake -DTIOGA_ENABLE_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES=90 \
      -DTIOGA_CUBQL_DIR=/path/to/cuBQL \
      -DTIOGA_HAS_NODEGID=off ..
make -j
```

Without `TIOGA_ENABLE_CUDA`, the library build is unchanged and the
GPU APIs return `-1`.

**Production-shaped run** (start here):

```
./benchmark/bench_multiblock.exe \
  --nblocks 256 --nx 16 --nq 2000 --batch --skip-dedup --verbose

# many-block sweeps → results/mb_*.csv
./benchmark/sweep_multiblock.sh
```

Useful flags: `--batch`, `--moving` / `--move`, `--no-refit`,
`--overlap`, `--eltype hex|prism|tet|mixed`, `--skip-dedup`.

**Single-block microbench** (useful, not the headline):

```
./benchmark/bench_search.exe --nx 64 --nq 500000 --verbose
./benchmark/sweep.sh   # → results/mesh_size.csv, etc.
```

Every run compares GPU `donorId[]` to CPU. A mismatch is *real* only
if the GPU donor fails an independent containment check; shared-face
ambiguities are benign.

Numbers in this repo were taken on **1× GH200 120GB (sm_90) + 72-core
Grace**, CUDA 12.9, nvhpc 25.7, `-O3`.

---

## Results that matter

**Batch vs per-block** (16³, 2k q/block, cached, dedup off):

| blocks | per-block Mqps | batched Mqps | gain |
|---|---|---|---|
| 64 | 32.5 | 239 | 7.4× |
| 256 | 30.8 | 192 | 6.2× |
| 1024 | 28.6 | 158 | 5.5× |

**Same total work, two decompositions** (4608×16³ = 18.9M cells,
9.2M points). GPU: 1 rank, `--batch`. CPU: 72 ranks × 64 blocks,
barrier-synced.

| Path | Mqps | vs 72-core socket |
|---|---|---|
| CPU socket, dedup included | 37.4 | 1.0× |
| GPU batched cached, dedup excluded | 153 | 3.8× |
| GPU batched cached, dedup included | 6.5 | **0.17×** |
| GPU compute only (moving/refit) | 1007 | **24.5×** |

So: kernels are fast; current TIOGA host packing + dedup are not.
Next leverage is dedup / keeping receptors device-resident — not
more traversal tuning. Detail in `GPU_SEARCH_SUMMARY.md` §4–§5.

---

## Caveats

- **Not in the connectivity driver.** Callers must opt into
  `search_gpu_batch()` (or `search_gpu()`).
- **End-to-end with dedup loses today** on the production-shaped case.
- **`ihigh != 0`** returns `-2`. **`uniform_hex`**: per-block falls
  back to host; batch returns `-2`.
- **No MPI / exchangeSearchData** in these harnesses — local synthetic
  query clouds only.
- **Dirty flags:** `setData()` sets `gpuMeshDirty`. Moving coords
  without re-`setData` need `gpuCoordsDirty` (batch path) or a full
  mesh dirty; production does not set these yet.
- **Benches default `TIOGA_HAS_NODEGID=off`.** nalu-wind often builds
  with NODEGID on (different dedup path).
- **Moving-mesh timings** re-send coords and refit, but often on
  undisplaced coordinates for the rate; correctness under `--move`
  displacement is checked separately.
- **Single-block speedups** (hundreds–thousands× walk) are real but
  not the production story — use `bench_multiblock --batch`.

## Knobs (`MeshBlock`)

| member | meaning |
|---|---|
| `gpuMeshDirty` | Re-upload mesh + rebuild BVH (`setData` sets this) |
| `gpuCoordsDirty` | Batch: coords-only upload + refit/rebuild |
| `gpuRefit` | Prefer `cuBQL::cuda::refit` when coords move |
| `gpuBuilderType` | 0 median, 1 radix (best here), 2 rebin, 3 SAH |
| `gpuLeafSize` | Leaf threshold; 0 = builder default |
| `gpuEarlyExit` | 1 = first clean donor (ADT-like); 0 = min cell id |
| `gpuSkipDedup` | Measurement only; breaks exchange-ready `xtag` |
| `searchTimers` | Phase breakdown (both backends) |
