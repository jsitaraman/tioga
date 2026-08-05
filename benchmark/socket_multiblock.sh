#!/usr/bin/env bash
# Socket-level comparison for the production shape:
#
#   1 rank owning N blocks + 1 GH200      (one CUDA context)
#   vs
#   72 ranks each owning N blocks, CPU only   (one full Grace socket)
#
# Both do the same work per rank, so the CPU column is the aggregate the
# socket can sustain and the GPU column is what one GPU sustains.
#
#   usage: ./socket_multiblock.sh [path/to/bench_multiblock.exe] [outdir]
set -e

BENCH=${1:-$(dirname "$0")/../build-gpu/benchmark/bench_multiblock.exe}
OUT=${2:-$(dirname "$0")/results}
mkdir -p "$OUT"
F="$OUT/mb_socket.csv"
NRANK=${NRANK:-72}

echo "nblocks,nx,nq,ranks,cpu_1core_Mqps,cpu_socket_Mqps,cpu_efficiency,gpu_Mqps,gpu_cached_Mqps,gpu_vs_socket,gpu_cached_vs_socket" > "$F"

for cfg in "64 16 2000" "256 16 2000" "1024 16 2000"; do
  set -- $cfg; nb=$1; nx=$2; nq=$3

  # 1 rank + GPU, one CUDA context
  one=$($BENCH --reps 3 --nblocks $nb --nx $nx --nq $nq)
  cpu1=$(echo "$one" | cut -d, -f26)
  gpu=$(echo  "$one" | cut -d, -f27)
  gpuc=$(echo "$one" | cut -d, -f28)

  # NRANK concurrent CPU-only ranks pinned to distinct cores
  rm -f /tmp/mbsock_*.csv
  for ((i=0;i<NRANK;i++)); do
    taskset -c $i $BENCH --reps 3 --cpu-only --nblocks $nb --nx $nx --nq $nq \
      > /tmp/mbsock_$i.csv 2>/dev/null &
  done
  wait

  python3 - "$nb" "$nx" "$nq" "$NRANK" "$cpu1" "$gpu" "$gpuc" "$F" <<'PY'
import glob, sys
nb, nx, nq, nrank = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
cpu1, gpu, gpuc, out = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
tot = []
for f in glob.glob('/tmp/mbsock_*.csv'):
    tot.append(float(open(f).read().strip().split(',')[8]))   # cpu_total
n = len(tot)
sock = n * int(nq) * int(nb) / max(tot) / 1e6
eff  = sock / cpu1 / n
row = (f"{nb},{nx},{nq},{n},{cpu1:.4f},{sock:.4f},{eff:.3f},"
       f"{gpu:.4f},{gpuc:.4f},{gpu/sock:.2f},{gpuc/sock:.2f}")
open(out, 'a').write(row + "\n")
print(row)
PY
done

echo "written to $F"
