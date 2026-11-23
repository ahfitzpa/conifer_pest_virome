#!/usr/bin/env python3
"""Distil the ICTV VMR into a small committed taxon -> lineage table.

scripts/family_taxonomy.py reconciles the three taxonomy callers against this
table. It is committed, unlike references/vmr_msl41.tsv (4.5 MB, gitignored),
so a clone can reconcile family calls without first fetching the VMR.

Regenerate after an MSL bump:
  python3 -I scripts/build_family_lineage.py
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RANKS = ["Phylum", "Class", "Order", "Family"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vmr", default=str(ROOT / "references/vmr_msl41.tsv"))
    ap.add_argument("--out", default=str(ROOT / "common/taxonomy/family_lineage.tsv"))
    a = ap.parse_args()

    # 1: collect each taxon's rank and its ancestors, at every rank we compare on
    lineage: dict[str, dict[str, str]] = {}
    with open(a.vmr, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            path = {r: (row.get(r) or "").strip() for r in RANKS}
            for i, rank in enumerate(RANKS):
                name = path[rank]
                if not name:
                    continue
                rec = lineage.setdefault(name, {"name": name, "rank": rank.lower()})
                for anc in RANKS[:i]:
                    if path[anc]:
                        rec[anc.lower()] = path[anc]

    # 2: write
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["name", "rank", "phylum", "class", "order"]
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for name in sorted(lineage):
            w.writerow({c: lineage[name].get(c, "") for c in cols})
    print(f"[lineage] {len(lineage)} taxa -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
