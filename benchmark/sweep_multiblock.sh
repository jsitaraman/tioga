#!/usr/bin/env bash
# Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause
# Many blocks per rank, one CUDA context per GPU -- the production shape.
# 16^3 active cells per block is typical for the target application.
#
#   usage: ./sweep_multiblock.sh [path/to/bench_multiblock.exe] [outdir]
set -e

BENCH=${1:-$(dirname "$0")/../build-gpu/benchmark/bench_multiblock.exe}
OUT=${2:-$(dirname "$0")/results}
mkdir -p "$OUT"
REPS=3

run_study () {
  local name=$1; shift
  local f="$OUT/$name.csv"
  "$BENCH" --header > "$f"
  echo "### $name -> $f"
  while read -r line; do
    [[ -z "$line" ]] && continue
    # shellcheck disable=SC2086
    "$BENCH" --reps $REPS $line >> "$f"
    tail -1 "$f" | cut -d, -f1,2,5,7,9,14,20,23,24,25,26,27,28,29
  done
}

# 1. how many 16^3 blocks does a rank own? 2000 receptor points each.
run_study mb_nblocks <<'EOF'
--nx 16 --nq 2000 --nblocks 1
--nx 16 --nq 2000 --nblocks 4
--nx 16 --nq 2000 --nblocks 16
--nx 16 --nq 2000 --nblocks 64
--nx 16 --nq 2000 --nblocks 256
--nx 16 --nq 2000 --nblocks 1024
EOF

# 2. same, with the shared host dedup removed, so the per-block GPU
#    overhead is visible on its own
run_study mb_nblocks_nodedup <<'EOF'
--skip-dedup --nx 16 --nq 2000 --nblocks 1
--skip-dedup --nx 16 --nq 2000 --nblocks 4
--skip-dedup --nx 16 --nq 2000 --nblocks 16
--skip-dedup --nx 16 --nq 2000 --nblocks 64
--skip-dedup --nx 16 --nq 2000 --nblocks 256
--skip-dedup --nx 16 --nq 2000 --nblocks 1024
EOF

# 3. receptor points per block, at a fixed 256 blocks
run_study mb_nq <<'EOF'
--nx 16 --nblocks 256 --nq 250
--nx 16 --nblocks 256 --nq 1000
--nx 16 --nblocks 256 --nq 4000
--nx 16 --nblocks 256 --nq 16000
EOF

# 4. block size at fixed total work (~1M cells, ~512k queries)
run_study mb_blocksize <<'EOF'
--nblocks 256 --nx 16 --nq 2000
--nblocks 32  --nx 32 --nq 16000
--nblocks 4   --nx 64 --nq 128000
EOF

echo
echo "written to $OUT"
