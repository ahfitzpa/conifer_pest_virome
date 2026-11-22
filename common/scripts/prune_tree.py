#!/usr/bin/env python3
"""Prune a tree to the taxa in an alignment, and strip support labels.

HyPhy and PAML both need a tree over exactly the alignment's taxa. The shared
phylogeny's tree carries reference genomes and any contig the codon QC dropped,
so it has to be pruned rather than passed through.

Support labels are stripped because IQ-TREE writes them as `98.5/95`, which
neither tool will parse. Topology and branch lengths are preserved.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

from Bio import Phylo


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def alignment_ids(path: str) -> list[str]:
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tree", required=True)
    p.add_argument("--alignment", required=True)
    p.add_argument("--out-tree", required=True)
    p.add_argument("--out-report")
    p.add_argument("--format", default="newick")
    p.add_argument("--strip-support", action="store_true")
    p.add_argument("--min-taxa", type=int, default=4)
    args = p.parse_args()

    tree = Phylo.read(args.tree, args.format)
    tips = [t.name for t in tree.get_terminals() if t.name]
    want = alignment_ids(args.alignment)
    keep = [t for t in tips if t in set(want)]
    missing = [w for w in want if w not in set(tips)]
    drop = [t for t in tips if t not in set(want)]

    log(f"[prune] tree {len(tips)} tips, alignment {len(want)} taxa, "
        f"keeping {len(keep)}, pruning {len(drop)}")
    if missing:
        # 1: an alignment taxon absent from the tree is a wiring error, not a
        # filter: say so loudly rather than quietly analysing fewer taxa
        log(f"[prune] WARNING: {len(missing)} alignment taxa are not in the "
            f"tree: {', '.join(missing[:5])}")

    if len(keep) < args.min_taxa:
        log(f"[prune] only {len(keep)} shared taxa; need {args.min_taxa}")
        if args.out_report:
            _write_report(args.out_report, keep, drop, missing)
        return 1

    # 2: prune, then collapse the single-child nodes pruning leaves behind
    for name in drop:
        tree.prune(target=name)

    if args.strip_support:
        # 3: clear internal labels and confidences; terminals keep their names
        for clade in tree.find_clades():
            clade.confidence = None
            if not clade.is_terminal():
                clade.name = None

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tree)) or ".", exist_ok=True)
    Phylo.write(tree, args.out_tree, args.format)

    final = [t.name for t in Phylo.read(args.out_tree, args.format).get_terminals()]
    log(f"[prune] wrote {len(final)} tips -> {args.out_tree}")
    if args.out_report:
        _write_report(args.out_report, final, drop, missing)
    return 0


def _write_report(path: str, kept, dropped, missing):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["taxon", "status"])
        for t in sorted(kept):
            w.writerow([t, "kept"])
        for t in sorted(dropped):
            w.writerow([t, "pruned_not_in_alignment"])
        for t in sorted(missing):
            w.writerow([t, "in_alignment_but_not_in_tree"])


if __name__ == "__main__":
    raise SystemExit(main())
