#!/usr/bin/env python3
"""Summarise one family's ML tree: taxa, model, support and new-lineage count."""
from __future__ import annotations
import argparse, csv, sys
from collections import Counter

from Bio import Phylo


def read_tsv(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--support", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--min-ufb", type=float, default=70)
    ap.add_argument("--min-alrt", type=float, default=90)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tree = Phylo.read(a.tree, "newick")
    tips = [t.name for t in tree.find_clades(terminal=True) if t.name]
    meta = {r["taxon_id"]: r for r in read_tsv(a.metadata)}
    groups = Counter(meta.get(t, {}).get("group", "NA") for t in tips)
    support = {r["metric"]: r["value"] for r in read_tsv(a.support)}

    # IQ-TREE writes SH-aLRT/UFB as an internal-node label, e.g. "85.6/92".
    n_label = n_supported = 0
    for clade in tree.find_clades():
        if clade.is_terminal() or not clade.name:
            continue
        parts = clade.name.split("/")
        if len(parts) != 2:
            continue
        n_label += 1
        try:
            alrt, ufb = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        if alrt >= a.min_alrt and ufb >= a.min_ufb:
            n_supported += 1

    row = {
        "family": a.family,
        "n_taxa": len(tips),
        "n_candidates": groups.get("candidate", 0),
        "n_ictv": groups.get("ICTV", 0),
        "model": support.get("model", "NA"),
        "n_sites": support.get("n_sites", "NA"),
        "log_likelihood": support.get("log_likelihood", "NA"),
        "n_internal_labels": n_label,
        "n_supported_internal": n_supported,
        "support_rule": f"SH-aLRT>={a.min_alrt:g} and UFB>={a.min_ufb:g}",
        "rooting": "midpoint" if "midpoint" in a.tree.lower() else "as_inferred",
        "bootstrap": "not calculated" if "UFB" not in open(a.support).read() else "UFB+SH-aLRT",
    }
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        w.writeheader(); w.writerow(row)
    print(f"{a.out}: {row['n_taxa']} taxa, {n_supported}/{n_label} supported nodes",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
