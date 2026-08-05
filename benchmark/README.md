# GPU donor search for TIOGA (cuBQL) — benchmark

`MeshBlock::search_gpu()` is an alternative to `MeshBlock::search()` that
replaces the host-side alternating digital tree with a
[cuBQL](https://github.com/NVIDIA/cuBQL) BVH built and traversed on the GPU.
Both produce the same `donorId[]`, `donorCount` and `xtag[]`.

## What changed

| | `search()` (host) | `search_gpu()` (device) |
|---|---|---|
| broad phase input | cells overlapping the query-cloud OBB | **all** cells of the block |
| pre-filter | full pass over every cell, LIFO compaction, second pass to build world AABBs | none |
| acceleration structure | 6-D ADT, rebuilt every call | 3-D cuBQL BVH, cached until `gpuMeshDirty` |
| traversal | recursive, one unique point at a time | one thread per point, `cuBQL::fixedBoxQuery::forEachPrim` |
| exact containment | `checkContainment` → `computeNodalWeights` | device port of the same code, same `TOL`, double precision |
| duplicate points | `uniquenodes_octree` on the host | same call (see "the remaining bottleneck") |

The query-cloud OBB pre-filter is dropped on purpose: it costs a full pass
over every cell, which is the same order as just building the BVH over every
cell. Measurements below confirm that trade.

Only the BVH bounds are single precision, and they are rounded strictly
outward (`nextafterf`) so the broad phase can never produce a false negative.
All geometry — the trilinear Newton solve, the 3×3 eliminations, the
`[-TOL, 1+TOL]` weight test — runs in double, from a line-by-line port of
`math.c` and `checkContainment.C`.

## Building

```
mkdir build-gpu && cd build-gpu
cmake -DTIOGA_ENABLE_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES=90 \
      -DTIOGA_CUBQL_DIR=/scratch/software/cuBQL \
      -DTIOGA_HAS_NODEGID=off ..
make -j
```

cuBQL is used header-only — `searchGPU.cu` defines
`CUBQL_GPU_BUILDER_IMPLEMENTATION` and instantiates the `float3` builders it
needs, so no cuBQL library is built or linked. Without `TIOGA_ENABLE_CUDA`
the build is unchanged and `search_gpu()` returns `-1`.

## Running

```
./benchmark/bench_search.exe --nx 64 --nq 500000 --verbose
./benchmark/sweep.sh                     # writes benchmark/results/*.csv
./benchmark/summarize.py                 # renders those CSVs as tables
```

The benchmark builds one block of `nx^3` cells with interior nodes displaced
by a smooth field, so the cells are genuinely non-affine and the containment
test has to run the full Newton iteration. Outer-boundary nodes are declared
overset boundary nodes, so `tagBoundary()` marks the outer cell layers with
`cellRes=BIGVALUE` and the donor-rejection logic is exercised. `--eltype`
selects hex / prism / tet / mixed (hex+prism), all carved from the same
logical hexes so the geometry is identical and conforming.

Every run cross-checks GPU against CPU donor-by-donor. A mismatch is only
counted as *real* if the GPU's donor does not actually contain the point
(checked independently against `computeNodalWeights`); a point lying on a
face shared by two cells is legitimately ambiguous and is counted as benign.

## Results

GH200 120GB (sm_90) / 72-core Grace, CUDA 12.9, nvhpc 25.7, `-O3`, single
rank, best of 3 reps. All times in milliseconds.

**Real mismatches across the entire sweep: 0.** (One benign mismatch out of
50,000 appeared in the mixed hex/prism mesh — a point on a shared face where
both answers are valid donors.)

### Growing block, fixed 100k query points

| cells | CPU total | filter | ADT build | dedup | walk | GPU total | BVH | walk | ×total | ×walk |
|---|---|---|---|---|---|---|---|---|---|---|
| 4,096 | 98.3 | 4.4 | 1.8 | 24.2 | 68.0 | 25.3 | 0.66 | 0.079 | 3.9 | 864 |
| 32,768 | 153.4 | 5.1 | 17.0 | 24.3 | 107.0 | 25.8 | 0.87 | 0.084 | 5.9 | 1279 |
| 262,144 | 339.1 | 12.6 | 158.6 | 24.2 | 143.6 | 27.5 | 1.43 | 0.089 | 12.3 | 1617 |
| 884,736 | 852.0 | 29.4 | 599.8 | 24.3 | 198.6 | 30.5 | 2.29 | 0.134 | 27.9 | 1486 |
| 1,404,928 | 1295.9 | 44.5 | 985.7 | 24.3 | 241.4 | 33.1 | 3.00 | 0.137 | 39.2 | 1769 |

The ADT rebuild is what makes the host cost superlinear: at 1.4M cells it is
76 % of `search()`. The cuBQL build over the *same* 1.4M cells takes 3.0 ms —
330× less — and is skipped entirely on the next call if the mesh has not moved.

### Growing query cloud, fixed 262k-cell block

| queries | CPU walk | GPU walk | ×walk | GPU query throughput |
|---|---|---|---|---|
| 16,000 | 23.1 | 0.044 | 523 | 0.36 G/s |
| 256,000 | 363.0 | 0.164 | 2217 | 1.56 G/s |
| 1,000,000 | 1425.4 | 0.506 | 2817 | 1.98 G/s |
| 4,000,000 | 5770.0 | 1.889 | 3055 | 2.12 G/s |

Traversal + exact containment saturates at about **2.1 billion point
locations per second**, against 0.7 M/s on one core.

### Search only (host duplicate-point pass removed, `--skip-dedup`)

| cells | queries | CPU total | GPU total | GPU cached | ×total | ×cached |
|---|---|---|---|---|---|---|
| 4,096 | 500,000 | 502.5 | 2.35 | 1.41 | 213 | 357 |
| 262,144 | 500,000 | 1021.4 | 4.35 | 1.47 | 235 | 697 |
| 1,404,928 | 500,000 | 2562.4 | 10.42 | 1.67 | 246 | 1539 |
| 262,144 | 4,000,000 | 7357.4 | 15.79 | 12.48 | 466 | 590 |

`GPU total` re-uploads the mesh and rebuilds the BVH every call (moving mesh);
`GPU cached` reuses both (static mesh).

### Does dropping the OBB pre-filter hurt?

`qfrac` is the edge of the query cloud as a fraction of the block, so small
`qfrac` is exactly the regime the host pre-filter was written for. 884,736
cells, 500,000 queries, `--skip-dedup`:

| qfrac | cells surviving the host filter | CPU filter | CPU ADT | CPU total | GPU total | GPU cached | ×total |
|---|---|---|---|---|---|---|---|
| 0.05 | 742 | 38.1 | 0.3 | 424.4 | 7.41 | 1.44 | 57 |
| 0.10 | 4,521 | 37.7 | 1.8 | 485.2 | 7.47 | 1.45 | 65 |
| 0.25 | 59,338 | 39.1 | 31.1 | 744.1 | 7.49 | 1.50 | 99 |
| 0.50 | 426,732 | 43.5 | 264.1 | 1183.3 | 7.52 | 1.53 | 157 |
| 1.00 | 884,736 | 50.8 | 598.9 | 1765.2 | 7.57 | 1.65 | 233 |

The filter does its job — at `qfrac=0.05` it cuts the ADT build from 599 ms to
0.3 ms — but the *scan* that produces it costs 38 ms, which is five times the
GPU's entire end-to-end search over all 884,736 cells, and twenty-six times
the GPU's cached-BVH time. There is no regime in this sweep where keeping the
pre-filter would help the GPU path.

### Element type (200k queries)

| type | cells | CPU total | CPU walk | GPU walk | ×walk |
|---|---|---|---|---|---|
| hex | 110,592 | 393.0 | 267.5 | 0.159 | 1686 |
| prism | 221,184 | 520.4 | 322.6 | 0.194 | 1661 |
| tet | 663,552 | 923.7 | 421.9 | 0.226 | 1871 |
| mixed | 165,888 | 475.5 | 311.3 | 0.207 | 1502 |

Element type barely moves GPU traversal: the tet mesh has 6× the cells of the
hex mesh and costs 42 % more traversal time, and the mixed hex/prism mesh
(which does branch on `nvert` within a warp) costs 30 % more than pure hex at
1.5× the cells. Divergence from mixed element types is not large enough here
to justify bucketing candidates by element type.

### cuBQL builder choice (512k cells)

| builder | leaf | BVH build | traversal |
|---|---|---|---|
| spatial median (default) | default | 1.71 | 0.132 |
| radix | default | **0.78** | **0.101** |
| rebin radix | default | 0.76 | 0.115 |
| SAH | default | 17.22 | 0.101 |
| spatial median | 4 | 1.36 | 0.229 |
| spatial median | 8 | 1.23 | 0.324 |
| spatial median | 16 | 1.12 | 0.534 |
| radix | 8 | 0.63 | 0.191 |

`cuBQL::cuda::radixBuilder` (`gpuBuilderType=1`) is the best choice for this
workload: 2.2× faster to build than the default and marginally faster to
traverse. SAH costs 10× more to build and buys nothing. Larger leaves trade
build time for traversal time badly — the containment test is expensive
enough that brute-forcing a leaf does not pay.

## The remaining bottleneck

With the search itself 200–1400× faster, the host duplicate-query-point pass
(`uniquenodes_octree`, unchanged in both backends) becomes the whole cost:

```
262,144 cells, 1,000,000 queries
  CPU  total 1939.5 ms  (filter 51.0  ADT 158.6  dedup 304.4  walk 1425.4)
  GPU  total  308.8 ms  (xfer 2.3  BVH 1.5  dedup 302.8  walk 0.5)
                                             ^^^^^^^^^^^ 98 % of the total
```

End-to-end this is 6.3×; with the dedup removed the same case is 326×. So the
next thing to port is the unique-node pass, not anything in the search. Two
options, in increasing order of effort:

1. Build with `TIOGA_HAS_NODEGID=on` (TIOGA's default, and what nalu-wind
   uses). Dedup then keys on 64-bit global node IDs instead of doing
   TOL-proximity matching in an octree, which maps directly onto a
   sort-by-key plus segmented scan on the device.
2. Port `uniqNodesTree` itself. Harder, because TOL-proximity merging is not
   an equivalence relation and the host result depends on octree leaf
   ordering.

Note that the GPU search does not *need* the dedup at all — it searches every
point regardless, and duplicates cost roughly nothing. The pass only exists
to produce `xtag`/`res_search` for the downstream donor exchange
(`bookKeeping.C`).

## Knobs on MeshBlock

| member | meaning |
|---|---|
| `gpuMeshDirty` | set to 1 whenever `x[]` or connectivity changes; forces re-upload and BVH rebuild. Set automatically by `setData()`. |
| `gpuBuilderType` | 0 = spatial median (cuBQL default), 1 = radix, 2 = rebin radix, 3 = SAH |
| `gpuLeafSize` | cuBQL `makeLeafThreshold`, 0 = builder default |
| `gpuEarlyExit` | 1 = stop at the first accepted donor (matches the ADT), 0 = scan all candidates and keep the lowest cell id (deterministic) |
| `gpuSkipDedup` | 1 = skip the host duplicate-point pass. `donorId` is unaffected; `xtag`/`res_search` are then not set up for the donor exchange. For measurement. |
| `searchTimers` | per-call phase breakdown, filled by both backends |

## Not covered

- `uniform_hex` blocks fall through to `search_uniform_hex()` on the host.
  Direct logical indexing already beats any tree; a BVH would be a regression.
- `ihigh != 0` (high-order) returns `-2`: containment goes through host
  callbacks (`donor_inclusion_test`) that cannot run on the device.
- Multi-block / multi-rank. One BVH per `MeshBlock` per rank; aggregating
  several small blocks into one BVH is the obvious next step for strong
  scaling, and matters most where blocks are small.
- Temporal coherence (test the previous donor, then its neighbours, before
  falling back to the BVH). For moving overset meshes this can turn most
  searches into O(1) and would sit on top of this.
