#!/usr/bin/env python3
"""Extract named model/support statistics from an IQ-TREE 2 .iqtree report.

The previous run's distance was unnamed and uninterpretable; every value here
is written with a metric name and a unit-free numeric value so downstream
tables and figures can cite it unambiguously.
"""
from __future__ import annotations
import argparse, csv, re, sys

PATTERNS = {
    "model": r"^Model(?: of evolution)?:\s*(\S+)",
    "n_sites": r"^Number of sites:\s*(\d+)",
    "n_taxa": r"^Number of taxa:\s*(\d+)",
    "log_likelihood": r"^Log likelihood:\s*([-0-9.]+)",
    "ufb_replicates": r"^Ultrafast bootstrap replicates:\s*(\d+)",
    "alrt_replicates": r"^SH-like approximate likelihood ratio test replicates:\s*(\d+)",
    "elapsed_seconds": r"^Total time:\s*([0-9.]+)",
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--iqtree", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    text = open(a.iqtree, encoding="utf-8", errors="ignore").read()
    rows = []
    for metric, pat in PATTERNS.items():
        m = re.search(pat, text, re.M | re.I)
        rows.append({"metric": metric, "value": m.group(1) if m else "NA"})
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["metric", "value"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"{a.out}: {len(rows)} named metrics", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
