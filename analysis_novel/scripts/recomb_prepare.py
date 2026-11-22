#!/usr/bin/env python3
"""Gate recombination screening on reference density and alignment length.

Writes the eligibility decision (with the numbers behind it) to a TSV -- never
silently. Families below the threshold produce an empty FASTA and a TSV row
with eligible=0 and a reason, so the summary can report 'not assessable' with
an explicit cause instead of a misleading negative result.
"""
from __future__ import annotations
import argparse, csv, sys


def read_tsv(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--sequences", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--min-refs", type=int, required=True)
    ap.add_argument("--min-aln", type=int, required=True)
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-eligible", required=True)
    a = ap.parse_args()

    meta = read_tsv(a.metadata)
    n_refs = sum(1 for r in meta if r.get("group") == "ICTV")
    n_cand = sum(1 for r in meta if r.get("group") == "candidate")

    seqs, name, chunks = {}, None, []
    with open(a.sequences, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].strip(), []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)

    lengths = {n: len(s.replace("-", "").replace(".", "")) for n, s in seqs.items()}
    aln_len = max(lengths.values()) if lengths else 0
    eligible = int(n_refs >= a.min_refs and aln_len >= a.min_aln and n_cand >= 2)
    reasons = []
    if n_refs < a.min_refs:
        reasons.append(f"only {n_refs} ICTV references (need {a.min_refs})")
    if aln_len < a.min_aln:
        reasons.append(f"alignment {aln_len} nt (need {a.min_aln})")
    if n_cand < 2:
        reasons.append(f"only {n_cand} candidate sequences (need 2)")
    reason = "; ".join(reasons) if reasons else "reference density and length sufficient"

    with open(a.out_faa, "w", encoding="utf-8") as fh:
        if eligible:
            for n, s in seqs.items():
                fh.write(f">{n}\n{s}\n")
    with open(a.out_eligible, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family", "n_ictv_references", "n_candidate_sequences",
                    "alignment_length_nt", "min_references", "min_alignment_nt",
                    "eligible", "reason"])
        w.writerow([a.family, n_refs, n_cand, aln_len, a.min_refs, a.min_aln,
                    eligible, reason])
    print(f"{a.family}: eligible={eligible} ({reason})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
