#!/usr/bin/env python3
# Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause
"""Turn results/shmoo_search.csv into plot-ready throughput tables.

Writes results/shmoo_throughput.csv (tidy/long form, one row per
cells x queries x backend) and prints a markdown table.

Throughputs are millions of query points located per second:

  cpu_search   nq / (cpu_total - cpu_dedup)   OBB filter + ADT build + walk
  gpu_search   nq / gpu_total                  H2D + BVH build + walk
  gpu_cached   nq / gpu_total_cached           BVH and mesh reused (static mesh)
  cpu_walk     nq / cpu_query                  traversal + containment only
  gpu_walk     nq / gpu_query                  traversal + containment only

Reads the --skip-dedup grid, so the GPU columns are direct measurements. The
shared host duplicate-point pass is excluded everywhere: it is identical in
both backends and is not part of the search. (It cannot be skipped on the CPU
side, so cpu_search still subtracts its measured cost.)

    usage: ./throughput.py [results_dir]
"""
import csv
import os
import sys


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results")
    src = os.path.join(d, "shmoo_search.csv")
    with open(src) as fh:
        rows = list(csv.DictReader(fh))

    out = []
    for r in rows:
        nq = int(r["nq"])
        cells = int(r["ncells"])

        def rate(total, dedup=0.0):
            t = total - dedup
            return nq / t / 1e6 if t > 0 else float("nan")

        out.append(dict(
            cells=cells,
            queries=nq,
            cpu_search=rate(float(r["cpu_total"]), float(r["cpu_dedup"])),
            gpu_search=rate(float(r["gpu_total"])),
            gpu_cached=rate(float(r["gpu_total_cached"])),
            cpu_walk=rate(float(r["cpu_query"])),
            gpu_walk=rate(float(r["gpu_query"])),
            mismatch_real=int(r["mismatch_real"]),
        ))

    dst = os.path.join(d, "shmoo_throughput.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        for o in out:
            w.writerows([{k: (f"{v:.4f}" if isinstance(v, float) else v)
                          for k, v in o.items()}])
    print(f"# wrote {dst}\n")

    hdr = ["cells", "queries", "cpu_search", "gpu_search", "gpu_cached",
           "cpu_walk", "gpu_walk", "bad"]
    body = [[f"{o['cells']:,}", f"{o['queries']:,}",
             f"{o['cpu_search']:.3f}", f"{o['gpu_search']:.1f}",
             f"{o['gpu_cached']:.1f}", f"{o['cpu_walk']:.3f}",
             f"{o['gpu_walk']:.1f}", str(o["mismatch_real"])] for o in out]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "|".join("---" for _ in hdr) + "|")
    for b in body:
        print("| " + " | ".join(b) + " |")


if __name__ == "__main__":
    main()
