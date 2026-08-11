#!/usr/bin/env bash
# Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause
# 2-D shmoo over (mesh size x query count), for plotting throughput curves.
#
# Runs the grid twice:
#   results/shmoo.csv         as-shipped, including the shared host dedup
#   results/shmoo_search.csv  --skip-dedup, so the GPU phase times are direct
#                             measurements rather than differences of two
#                             larger numbers. Use this one for throughput.
#
#   usage: ./shmoo.sh [path/to/bench_search.exe] [outdir]
set -e

BENCH=${1:-$(dirname "$0")/../build-gpu/benchmark/bench_search.exe}
OUT=${2:-$(dirname "$0")/results}
mkdir -p "$OUT"

NX_LIST="16 24 32 48 64 80 96 112"
NQ_LIST="10000 100000 1000000 4000000"

grid () {                      # grid <outfile> [extra args...]
  local f=$1; shift
  "$BENCH" --header > "$f"
  for nx in $NX_LIST; do
    for nq in $NQ_LIST; do
      "$BENCH" --reps 2 --nx $nx --nq $nq "$@" >> "$f"
      echo "$(basename "$f") nx=$nx nq=$nq done"
    done
  done
  echo "wrote $f"
}

grid "$OUT/shmoo.csv"
grid "$OUT/shmoo_search.csv" --skip-dedup
echo "shmoo complete"
