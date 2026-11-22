#!/usr/bin/env python3
"""Pairwise amino-acid identity between each candidate and each ICTV exemplar.

Uses a local BLOSUM62 alignment (gap open -11, extend -1), the same parameters
used throughout Phase 3. Identity is reported only when >= 50% of the query is
covered, and the alignment length and coverage are written alongside it so the
number is traceable and comparable to the ICTV criterion.
"""
from __future__ import annotations
import argparse, csv, sys

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    seqs = {r.id: str(r.seq).replace("J", "X") for r in SeqIO.parse(a.alignment, "fasta")}
    with open(a.metadata, encoding="utf-8") as fh:
        meta = list(csv.DictReader(fh, delimiter="\t"))
    groups = {r["taxon_id"]: r.get("group", "") for r in meta}
    candidates = [t for t, g in groups.items() if g == "candidate"]
    refs = [t for t, g in groups.items() if g == "ICTV"]

    fields = ["query", "reference", "identity_pct", "query_coverage_pct",
              "aligned_cols", "score", "passes_50pct_coverage"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for q in candidates:
            if q not in seqs:
                continue
            for r in refs:
                if r not in seqs:
                    continue
                aln = aligner.align(seqs[q], seqs[r])[0]
                qa, sa = aln.aligned
                cov = 100.0 * sum(e - s for s, e in qa) / len(seqs[q])
                m = c = 0
                for (qs, qe), (ss, se) in zip(qa, sa):
                    for k in range(qe - qs):
                        c += 1
                        m += seqs[q][qs + k] == seqs[r][ss + k]
                w.writerow({"query": q, "reference": r,
                            "identity_pct": round(100.0 * m / c, 2) if c else 0.0,
                            "query_coverage_pct": round(cov, 2),
                            "aligned_cols": c, "score": float(aln.score),
                            "passes_50pct_coverage": int(cov >= 50)})
    print(f"{len(candidates)} candidates x {len(refs)} references -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
