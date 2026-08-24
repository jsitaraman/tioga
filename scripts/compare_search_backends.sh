#!/usr/bin/env bash
#
# Copyright TIOGA Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD 3-Clause)
# Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause
#
# Build every donor search backend, run the same case through each of them,
# check that they all pick the same donors as the host search, and print the
# per-phase cost of each.
#
# The backend is a configure time choice (see TIOGA_SEARCH_BACKEND), so the
# backends cannot be compared inside one run. Each gets its own build
# directory and its own run, and the donors each one found are diffed
# afterwards. The query points do not depend on the backend, so at a fixed
# rank count the dumps line up entry for entry.
#
# Usage: scripts/compare_search_backends.sh [-n ranks] [-o outdir]
#                                           [-c cubql_dir] [-a cuda_arch]
#                                           [-u] [-b backends] [-r repeats]
#
#   -n  MPI ranks for the solve (default 8). The grid generator runs on 2n.
#   -o  where builds, dumps and logs go (default build-compare)
#   -c  cuBQL source tree (default /scratch/software/cuBQL)
#   -a  CUDA architecture (default 90)
#   -u  leave host de-duplication on. The GPU backends stand down in that
#       configuration, so this checks that the host path is untouched.
#   -b  space separated backend list (default "cpu adt_gpu cubql cubql_batch")
#   -r  repeat the search this many times and time the last pass (default 1).
#       Passes after the first mark the coordinates dirty, which is what a
#       moving mesh does, so the cuBQL backends refit rather than rebuild.
#       Use this for a steady state number; repeat=1 is the first call only.

set -u

ranks=8
outdir=build-compare
cubql=/scratch/software/cuBQL
arch=90
uniqueid=off
backends="cpu adt_gpu cubql cubql_batch"
repeat=1

while getopts "n:o:c:a:b:r:uh" opt; do
  case $opt in
    n) ranks=$OPTARG ;;
    o) outdir=$OPTARG ;;
    c) cubql=$OPTARG ;;
    a) arch=$OPTARG ;;
    b) backends=$OPTARG ;;
    r) repeat=$OPTARG ;;
    u) uniqueid=on ;;
    h) sed -n '9,33p' "$0"; exit 0 ;;
    *) echo "try -h" >&2; exit 2 ;;
  esac
done

root=$(cd "$(dirname "$0")/.." && pwd)
case_dir=$root/case
outdir=$(mkdir -p "$outdir" && cd "$outdir" && pwd)
dumps=$outdir/dumps
logs=$outdir/logs
mkdir -p "$dumps" "$logs"

if [ ! -d "$case_dir/grid" ]; then
  echo "no $case_dir/grid; nothing to run" >&2
  exit 1
fi

echo "== configuration"
echo "   ranks            $ranks"
echo "   uniqueid (dedup) $uniqueid"
echo "   backends         $backends"
echo "   search repeats   $repeat"
echo "   cuBQL            $cubql"
echo "   output           $outdir"

# ---------------------------------------------------------------- build ----
for be in $backends; do
  bdir=$outdir/build-$be
  echo "== building $be"
  mkdir -p "$bdir"
  cuda=ON
  [ "$be" = cpu ] && cuda=OFF
  if ! cmake -S "$root" -B "$bdir" \
        -DTIOGA_ENABLE_CUDA=$cuda \
        -DCMAKE_CUDA_ARCHITECTURES="$arch" \
        -DTIOGA_CUBQL_DIR="$cubql" \
        -DTIOGA_HAS_NODEGID=off \
        -DTIOGA_ENABLE_UNIQUEID=$uniqueid \
        -DTIOGA_SEARCH_BACKEND=$be > "$logs/cmake-$be.log" 2>&1; then
    echo "   configure failed, see $logs/cmake-$be.log" >&2
    exit 1
  fi
  if ! cmake --build "$bdir" -j "$(nproc)" > "$logs/build-$be.log" 2>&1; then
    echo "   build failed, see $logs/build-$be.log" >&2
    exit 1
  fi
done

# ------------------------------------------------------------------ run ----
# The grid generator is backend independent; run it once, with the first
# build, so every backend searches identical grids.
first=$(echo "$backends" | awk '{print $1}')
echo "== generating grids ($((ranks*2)) parts)"
( cd "$case_dir/grid" && mpirun -np $((ranks*2)) --oversubscribe \
    "$outdir/build-$first/gridGen/buildGrid" ) > "$logs/gridgen.log" 2>&1 \
  || { echo "   gridGen failed, see $logs/gridgen.log" >&2; exit 1; }

for be in $backends; do
  echo "== running $be on $ranks ranks"
  rm -f "$dumps"/donors.$be.*.txt
  ( cd "$case_dir" && TIOGA_DONOR_DUMP="$dumps" TIOGA_SEARCH_TIMERS=1 \
      TIOGA_SEARCH_REPEAT="$repeat" \
      mpirun -np "$ranks" --oversubscribe \
      "$outdir/build-$be/driver/tioga.exe" ) > "$logs/run-$be.log" 2>&1 \
    || { echo "   run failed, see $logs/run-$be.log" >&2; exit 1; }
done

# -------------------------------------------------------------- compare ----
echo
echo "== donor agreement (reference: cpu)"
python3 - "$dumps" "$backends" <<'PY'
import sys, glob, os, re

dumps, backends = sys.argv[1], sys.argv[2].split()
ref = 'cpu' if 'cpu' in backends else backends[0]

def load(be):
    out = {}
    for f in sorted(glob.glob(os.path.join(dumps, f'donors.{be}.*.txt'))):
        key = re.sub(r'^donors\.[^.]+\.', '', os.path.basename(f))
        rows = []
        with open(f) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                p = line.split()
                rows.append((int(p[1]), p[2], p[3], p[4]))
        out[key] = rows
    return out

base = load(ref)
npts = sum(len(v) for v in base.values())
if npts == 0:
    print('  no donor dumps found; did the runs write them?')
    sys.exit(1)
print(f'  {ref}: {len(base)} dumps, {npts} query points, '
      f'{sum(1 for v in base.values() for r in v if r[0] > -1)} with a donor')

bad = 0
for be in backends:
    if be == ref:
        continue
    cur = load(be)
    if set(cur) != set(base):
        print(f'  {be:12s} DUMP SET MISMATCH')
        bad += 1
        continue
    diff = miss = extra = 0
    example = None
    for key in base:
        a, b = base[key], cur[key]
        if len(a) != len(b):
            print(f'  {be:12s} LENGTH MISMATCH in {key}')
            bad += 1
            break
        for i, (x, y) in enumerate(zip(a, b)):
            if x[0] == y[0]:
                continue
            if y[0] == -1:
                miss += 1
            elif x[0] == -1:
                extra += 1
            else:
                diff += 1
            if example is None:
                example = (key, i, x[0], y[0], x[1], x[2], x[3])
    total = diff + miss + extra
    status = 'OK' if total == 0 else 'DIFFERS'
    print(f'  {be:12s} {status:8s} '
          f'different-donor={diff} missed={miss} found-extra={extra} '
          f'({100.0*total/npts:.4f}% of points)')
    if example and total:
        k, i, xa, xb, px, py, pz = example
        print(f'               first at {k}[{i}]: {ref}={xa} {be}={xb} '
              f'point=({px}, {py}, {pz})')
    # a backend that finds no donor where the host does is a real failure;
    # naming a different cell for a point inside several cells is not
    if miss or extra:
        bad += 1
sys.exit(1 if bad else 0)
PY
agree=$?

echo
echo "== search cost per phase (seconds, summed over ranks / slowest rank)"
printf '   %-12s %10s %10s %10s %10s %10s %10s %10s\n' \
  backend total filter build dedup query transfer hostwork
for be in $backends; do
  awk -v be="$be" '
    /^#tioga search: [a-z]+ / {
      split($0,f," "); v[f[3]"_sum"]=f[4]; v[f[3]"_max"]=f[5]
    }
    END {
      printf "   %-12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
        be, v["total_sum"], v["filter_sum"], v["build_sum"], v["dedup_sum"],
        v["query_sum"], v["transfer_sum"], v["hostwork_sum"]
      printf "   %-12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
        "  (slowest)", v["total_max"], v["filter_max"], v["build_max"],
        v["dedup_max"], v["query_max"], v["transfer_max"], v["hostwork_max"]
    }' "$logs/run-$be.log"
done

echo
echo "== interpolation error (should be at round-off for every backend)"
for be in $backends; do
  err=$(awk '/Interpolation error/{f=1;next} f&&NF==2&&$1~/^[0-9]+$/{print $2}' \
        "$logs/run-$be.log" | sort -g | tail -1)
  printf '   %-12s worst rank %s\n' "$be" "${err:-<none found>}"
done

echo
if [ $agree -eq 0 ]; then
  echo "== donors agree with $( [ "$uniqueid" = on ] && echo 'the host search' \
       || echo 'cpu' ) everywhere that matters"
else
  echo "== DONOR MISMATCH, see above" >&2
fi
exit $agree
