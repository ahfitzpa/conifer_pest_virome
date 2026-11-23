#!/usr/bin/env python3
"""Run OpenRDP scanners on an alignment and parse events.

OpenRDP 0.1.0 exposes a Python API. We call each requested scanner, collect
the events, and keep those detected by >= min_methods distinct scanners
(the same consensus rule used by Zhao et al. and CaliciTyper). If OpenRDP
cannot be imported we fail loudly -- never emit 'not_assessed' silently.
"""
from __future__ import annotations
import argparse, csv, sys


def parse_detectnum(path, family):
    """Read the eligible file; return (eligible, reference_names, candidate_names)."""
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else None


def iter_fasta(path):
    names = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].strip())
    return names


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--eligible", required=True)
    ap.add_argument("--scanners", required=True)
    ap.add_argument("--min-methods", type=int, default=2)
    ap.add_argument("--bootstrap", type=int, default=100)
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    try:
        from openrdp import Scanner, aliases
    except Exception as exc:
        print(f"ERROR: OpenRDP import failed ({exc}). Check the env "
              "(numpy<2, OpenRDP from GitHub). Failing loudly, not writing not_assessed.",
              file=sys.stderr)
        return 3

    fields = ["query", "donor", "breakpoint", "n_scanners", "scanners", "status"]
    el = parse_detectnum(a.eligible, "")
    if not el or el.get("eligible") != "1":
        print("not assessable: " + (el.get("reason") if el else "no eligibility file") + " -- skipped",
              file=sys.stderr)
        with open(a.out_tsv, "w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=fields, delimiter="\t").writeheader()
        return 0
    names = iter_fasta(a.alignment)
    candidate_names = [n for n in names if n.startswith("cand|")]
    if not candidate_names:
        print("not assessable: no candidate sequences -- skipped", file=sys.stderr)
        with open(a.out_tsv, "w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=fields, delimiter="\t").writeheader()
        return 0

    scanner_list = [s.strip() for s in a.scanners.split() if s.strip()]
    valid_methods = []
    for s in scanner_list:
        key = "threeseq" if s.lower() == "3seq" else s.lower()
        if key not in aliases:
            print(f"WARNING: OpenRDP scanner {s} not found -- skipped", file=sys.stderr)
            continue
        if key not in valid_methods:
            valid_methods.append(key)

    grouped = {}
    if valid_methods:
        scanner = Scanner(methods=tuple(valid_methods))
        results = scanner.run_scans(a.alignment)
        for method in results.keys():
            for d in results[method]:
                q = d.get("recombinant", "NA")
                p1 = d.get("parent1", "-")
                p2 = d.get("parent2", "-")
                dn = f"{p1}/{p2}" if p2 and p2 != "-" else str(p1)
                bp = f"{d.get('start', 'NA')}-{d.get('end', 'NA')}"
                key = (str(q), str(dn), str(bp))
                grouped.setdefault(key, set()).add(method)

    with open(a.out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for (q, dn, bp), methods in sorted(grouped.items(), key=lambda x: -len(x[1])):
            w.writerow({"query": q, "donor": dn, "breakpoint": bp,
                        "n_scanners": len(methods),
                        "scanners": ";".join(sorted(methods)),
                        "status": "supported" if len(methods) >= a.min_methods else "weak"})
    print(f"openrdp: {len(grouped)} candidate events, "
          f"{sum(len(m) >= a.min_methods for m in grouped.values())} supported", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
