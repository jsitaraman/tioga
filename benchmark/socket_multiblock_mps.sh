#!/usr/bin/env bash
# The deployment configuration: a GH200 superchip has 72 Grace cores AND one
# Hopper GPU, so the question is not "1 core + GPU vs 72 cores" -- it is
# whether adding the GPU to a fully loaded socket helps.
#
# Runs 72 ranks, each owning NBLOCKS blocks, all sharing the single GPU
# through CUDA MPS (one context per rank, MPS multiplexes them). Compares
# against the same 72 ranks running CPU-only.
#
#   usage: ./socket_multiblock_mps.sh [path/to/bench_multiblock.exe] [outdir]
set -e

BENCH=${1:-$(dirname "$0")/../build-gpu/benchmark/bench_multiblock.exe}
OUT=${2:-$(dirname "$0")/results}
mkdir -p "$OUT"
F="$OUT/mb_socket_mps.csv"
NRANK=${NRANK:-72}
NBLOCKS=${NBLOCKS:-64}
NX=${NX:-16}
NQ=${NQ:-2000}

export CUDA_MPS_PIPE_DIRECTORY=${CUDA_MPS_PIPE_DIRECTORY:-/tmp/tioga-mps}
export CUDA_MPS_LOG_DIRECTORY=${CUDA_MPS_LOG_DIRECTORY:-/tmp/tioga-mps-log}
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
nvidia-cuda-mps-control -d
trap 'echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true' EXIT

echo "mode,ranks,nblocks,nx,nq,cpu_Mqps,gpu_rebuild_Mqps,gpu_cached_Mqps,gpu_rebuild_vs_cpu,gpu_cached_vs_cpu" > "$F"

for mode in "" "--skip-dedup"; do
  rm -f /tmp/mbmps_*.csv
  for ((i=0;i<NRANK;i++)); do
    taskset -c $i "$BENCH" --reps 3 --nblocks $NBLOCKS --nx $NX --nq $NQ $mode \
      > /tmp/mbmps_$i.csv 2>/dev/null &
  done
  wait
  python3 - "$mode" "$NBLOCKS" "$NX" "$NQ" "$F" <<'PY'
import glob, sys
mode = sys.argv[1] or "with-dedup"
nb, nx, nq, out = int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
c = g = gc = []
c, g, gc = [], [], []
for f in glob.glob('/tmp/mbmps_*.csv'):
    r = open(f).read().strip().split(',')
    c.append(float(r[8])); g.append(float(r[13])); gc.append(float(r[18]))
n = len(c); tot = nq * nb * n
mc, mg, mgc = max(c), max(g), max(gc)
row = (f"{mode},{n},{nb},{nx},{nq},{tot/mc/1e6:.2f},{tot/mg/1e6:.2f},"
       f"{tot/mgc/1e6:.2f},{mc/mg:.2f},{mc/mgc:.2f}")
open(out, 'a').write(row + "\n")
print(row)
PY
done

echo "written to $F"
