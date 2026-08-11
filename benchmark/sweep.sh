#!/usr/bin/env bash
# Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause
# Sweep MeshBlock::search (host ADT) against MeshBlock::search_gpu (cuBQL BVH)
# across problem sizes. Writes one CSV per study into the output directory.
#
#   usage: ./sweep.sh [path/to/bench_search.exe] [outdir]

set -e

BENCH=${1:-$(dirname "$0")/../build-gpu/benchmark/bench_search.exe}
OUT=${2:-$(dirname "$0")/results}

if [[ ! -x "$BENCH" ]]; then
  echo "benchmark executable not found: $BENCH" >&2
  exit 1
fi
mkdir -p "$OUT"

REPS=3

run_study () {           # run_study <name> <fixed args...> -- <varying arg sets>
  local name=$1; shift
  local f="$OUT/$name.csv"
  "$BENCH" --header > "$f"
  echo "### $name -> $f"
  while read -r line; do
    [[ -z "$line" ]] && continue
    # shellcheck disable=SC2086
    "$BENCH" --reps $REPS $line >> "$f"
    tail -1 "$f"
  done
}

# ------------------------------------------------------------------
# 1. mesh size sweep: fixed query cloud, growing block
# ------------------------------------------------------------------
run_study mesh_size <<'EOF'
--nx 16 --nq 100000
--nx 24 --nq 100000
--nx 32 --nq 100000
--nx 40 --nq 100000
--nx 48 --nq 100000
--nx 64 --nq 100000
--nx 80 --nq 100000
--nx 96 --nq 100000
--nx 112 --nq 100000
EOF

# ------------------------------------------------------------------
# 2. query count sweep: fixed block, growing query cloud
# ------------------------------------------------------------------
run_study query_count <<'EOF'
--nx 64 --nq 1000
--nx 64 --nq 4000
--nx 64 --nq 16000
--nx 64 --nq 64000
--nx 64 --nq 256000
--nx 64 --nq 1000000
--nx 64 --nq 4000000
EOF

# ------------------------------------------------------------------
# 3. query cloud extent: the regime the host OBB pre-filter targets.
#    Small qfrac means most of the block is rejected before the ADT is
#    built; qfrac > 1 puts query points outside the block entirely.
# ------------------------------------------------------------------
run_study query_extent <<'EOF'
--nx 64 --nq 100000 --qfrac 0.05
--nx 64 --nq 100000 --qfrac 0.1
--nx 64 --nq 100000 --qfrac 0.25
--nx 64 --nq 100000 --qfrac 0.5
--nx 64 --nq 100000 --qfrac 1.0
--nx 64 --nq 100000 --qfrac 1.4
EOF

# ------------------------------------------------------------------
# 4. element type: exercises the tet / prism / hex containment paths
# ------------------------------------------------------------------
run_study element_type <<'EOF'
--nx 48 --nq 200000 --eltype hex
--nx 48 --nq 200000 --eltype prism
--nx 48 --nq 200000 --eltype tet
--nx 48 --nq 200000 --eltype mixed
EOF

# ------------------------------------------------------------------
# 5. cuBQL builder / leaf size, at a size where the build matters
# ------------------------------------------------------------------
run_study builder <<'EOF'
--nx 80 --nq 100000 --builder 0
--nx 80 --nq 100000 --builder 1
--nx 80 --nq 100000 --builder 2
--nx 80 --nq 100000 --builder 3
--nx 80 --nq 100000 --builder 0 --leaf 4
--nx 80 --nq 100000 --builder 0 --leaf 8
--nx 80 --nq 100000 --builder 0 --leaf 16
--nx 80 --nq 100000 --builder 1 --leaf 8
EOF

# ------------------------------------------------------------------
# 6. search only. The host duplicate-query-point pass (uniquenodes_octree)
#    is identical in both backends and ends up dominating the GPU total,
#    so this study drops it to show what the search itself costs.
#    The GPU searches every point anyway, so donorId is unchanged.
# ------------------------------------------------------------------
run_study search_only <<'EOF'
--skip-dedup --nx 16 --nq 500000
--skip-dedup --nx 32 --nq 500000
--skip-dedup --nx 48 --nq 500000
--skip-dedup --nx 64 --nq 500000
--skip-dedup --nx 80 --nq 500000
--skip-dedup --nx 96 --nq 500000
--skip-dedup --nx 112 --nq 500000
--skip-dedup --nx 64 --nq 16000
--skip-dedup --nx 64 --nq 100000
--skip-dedup --nx 64 --nq 1000000
--skip-dedup --nx 64 --nq 4000000
EOF

# ------------------------------------------------------------------
# 7. does dropping the host OBB pre-filter hurt? Same query-extent sweep
#    as study 3 but on a big block and without the shared dedup, so the
#    filter's cost and benefit are both visible.
# ------------------------------------------------------------------
run_study filter_value <<'EOF'
--skip-dedup --nx 96 --nq 500000 --qfrac 0.05
--skip-dedup --nx 96 --nq 500000 --qfrac 0.1
--skip-dedup --nx 96 --nq 500000 --qfrac 0.25
--skip-dedup --nx 96 --nq 500000 --qfrac 0.5
--skip-dedup --nx 96 --nq 500000 --qfrac 1.0
EOF

echo
echo "all studies written to $OUT"
