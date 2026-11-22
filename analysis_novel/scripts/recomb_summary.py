#!/usr/bin/env python3
"""Merge per-method recombination tables into one auditable summary.

A candidate-virus event is 'supported' when >= min-methods independent
scanners agree on the same recombinant/breakpoint. Not-assessable families
carry the reason from recomb_prepare.tsv. No support value is invented.
"""
from __future__ import annotations
import argparse, csv, sys
from collections import defaultdict


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--event-tables", nargs="+", required=True)
    ap.add_argument("--min-methods", type=int, default=2)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    events = defaultdict(set)
    files = []
    for path in a.event_tables:
        rows = read_tsv(path)
        files.append((path, len(rows)))
        for r in rows:
            q = r.get("query", "NA")
            bp = str(r.get("breakpoint", "NA"))
            method = r.get("method") or r.get("scanners", "unknown")
            events[(q, bp)].add(method)

    fields = ["family", "query", "breakpoint", "n_methods", "methods", "status",
              "source_tables"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for (q, bp), methods in sorted(events.items(), key=lambda x: -len(x[1])):
            w.writerow({"family": a.family, "query": q, "breakpoint": bp,
                        "n_methods": len(methods),
                        "methods": ";".join(sorted(methods)),
                        "status": "supported" if len(methods) >= a.min_methods else "weak",
                        "source_tables": ";".join(p for p, n in files)})
    print(f"{a.family}: {len(events)} unique events from "
          f"{sum(n for _, n in files)} method rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
