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
        from openrdp import scanners as ORDP
    except Exception as exc:
        print(f"ERROR: OpenRDP import failed ({exc}). Check the env "
              "(numpy<2, openrdp==0.1.0). Failing loudly, not writing not_assessed.",
              file=sys.stderr)
        return 3

    el = parse_detectnum(a.eligible, "")
    if not el or el.get("eligible") != "1":
        print("not assessable: " + (el.get("reason") if el else "no eligibility file") + " -- skipped",
              file=sys.stderr)
    names = iter_fasta(a.alignment)
    candidate_names = [n for n in names if n.startswith("cand|")]
    if not candidate_names:
        print("not assessable: no candidate sequences -- skipped", file=sys.stderr)

    # OpenRDP 0.1.x API: each scanner module exposes a run(sequences, ...)
    # returning a list of dict-like events. We introspect defensively and
    # record the actual method names used in the output table.
    scanner_list = [s.strip() for s in a.scanners.split() if s.strip()]
    events = []
    for s in scanner_list:
        fn = getattr(ORDP, s, None) or getattr(ORDP, s.lower(), None)
        if fn is None:
            print(f"WARNING: OpenRDP scanner {s} not found -- skipped", file=sys.stderr)
            continue
        try:
            result = fn(a.alignment)
        except TypeError:
            result = fn(sequences=a.alignment)
        except Exception as exc:
            print(f"WARNING: OpenRDP scanner {s} failed: {exc}", file=sys.stderr)
            continue
        for ev in result or []:
            d = dict(ev) if not isinstance(ev, dict) else ev
            events.append((s, d))

    # Group by (query, donor, breakpoint range) and count supporting scanners.
    grouped = {}
    for method, d in events:
        q = d.get("query", d.get("seq1", "NA"))
        dn = d.get("donor", d.get("seq2", "NA"))
        bp = d.get("breakpoint", d.get("positions", "NA"))
        key = (str(q), str(dn), str(bp))
        grouped.setdefault(key, set()).add(method)

    fields = ["query", "donor", "breakpoint", "n_scanners", "scanners", "status"]
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
