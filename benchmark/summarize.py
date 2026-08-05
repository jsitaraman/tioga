#!/usr/bin/env python3
"""Render the sweep CSVs written by sweep.sh as readable tables.

    usage: ./summarize.py [results_dir]
"""
import csv
import os
import sys

COLS = [
    ("cells",      lambda r: f"{int(r['ncells']):,}"),
    ("el",         lambda r: r["eltype"]),
    ("queries",    lambda r: f"{int(r['nq']):,}"),
    ("qfrac",      lambda r: r["qfrac"]),
    ("bld",        lambda r: r["builder"]),
    ("leaf",       lambda r: r["leaf"]),
    ("CPU total",  lambda r: f"{f(r,'cpu_total')*1e3:.1f}"),
    ("  filter",   lambda r: f"{f(r,'cpu_filter')*1e3:.1f}"),
    ("  ADT",      lambda r: f"{f(r,'cpu_build')*1e3:.1f}"),
    ("  dedup",    lambda r: f"{f(r,'cpu_dedup')*1e3:.1f}"),
    ("  walk",     lambda r: f"{f(r,'cpu_query')*1e3:.1f}"),
    ("cand",       lambda r: f"{int(r['cpu_cand']):,}"),
    ("GPU total",  lambda r: f"{f(r,'gpu_total')*1e3:.2f}"),
    ("  xfer",     lambda r: f"{f(r,'gpu_xfer')*1e3:.2f}"),
    ("  BVH",      lambda r: f"{f(r,'gpu_build')*1e3:.2f}"),
    ("  dedup",    lambda r: f"{f(r,'gpu_dedup')*1e3:.2f}"),
    ("  walk",     lambda r: f"{f(r,'gpu_query')*1e3:.3f}"),
    ("GPU cached", lambda r: f"{f(r,'gpu_total_cached')*1e3:.2f}"),
    ("x total",    lambda r: f"{f(r,'speedup_total'):.1f}"),
    ("x cached",   lambda r: f"{f(r,'speedup_total_cached'):.1f}"),
    ("x walk",     lambda r: f"{f(r,'speedup_query'):.0f}"),
    ("bad",        lambda r: r["mismatch_real"]),
]


def f(row, key):
    return float(row[key])


def table(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return
    head = [c[0] for c in COLS]
    body = [[fn(r) for _, fn in COLS] for r in rows]
    w = [max(len(head[i]), *(len(b[i]) for b in body)) for i in range(len(head))]
    print(os.path.basename(path)[:-4] + "   (all times in ms)")
    print("  " + "  ".join(h.rjust(w[i]) for i, h in enumerate(head)))
    print("  " + "  ".join("-" * w[i] for i in range(len(head))))
    for b in body:
        print("  " + "  ".join(b[i].rjust(w[i]) for i in range(len(head))))
    print()


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results")
    order = ["mesh_size", "query_count", "query_extent", "element_type",
             "builder", "search_only"]
    names = sorted(n[:-4] for n in os.listdir(d) if n.endswith(".csv"))
    for n in sorted(names, key=lambda n: (order.index(n) if n in order else 99, n)):
        table(os.path.join(d, n + ".csv"))


if __name__ == "__main__":
    main()
