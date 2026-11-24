#!/usr/bin/env python3
"""Gate nucleotide recombination screening, and emit the alignment to screen.

Input is the scenario phylogeny's trimmed NUCLEOTIDE alignment
(results/phylogeny/<fam>/aln.trimmed.fasta, FAMSA -> ClipKIT) plus its
selected.tsv, which labels every row candidate or reference. Earlier versions
read phylo_prepare's combined.faa, an unaligned amino-acid concatenation:
HyPhy GARD and OpenRDP both parse their input as an aligned nucleotide matrix
and died on the first ragged sequence ("Sequence N is longer than the previous
sequences"), so no family was ever actually screened.

Fragmentary contigs are the other hazard: a 750 nt candidate inside an 8900
column alignment is ~92% gap, and gap-rich rows generate spurious breakpoints
in every scanner. Sequences covering less than --min-coverage of the columns
are dropped before screening, columns left all-gap are then removed, and both
counts are reported in the eligibility table.

Writes the eligibility decision with the numbers behind it -- never silently.
Families below a threshold get an empty FASTA and eligible=0 with a reason, so
the summary reports 'not assessable' with an explicit cause rather than a
misleading negative result.
"""
from __future__ import annotations
import argparse, csv, sys

# Characters carrying no state for coverage purposes.
MISSING = set("-.?nN~")
# OpenRDP accepts only A, T, G, C, N and the gap: "Alignment contains invalid
# characters Y. Sequences can only contain A,T,G,C,-,N." IUPAC ambiguity codes
# are common in these assemblies (Iflaviridae alone carries 52 Y and 19 R), so
# anything outside the accepted set becomes N, which every scanner treats as
# missing data rather than as a state.
ACCEPTED = set("ACGTN-")


def normalise(seq):
    """Upper-case, gaps unified to '-', anything else not ACGTN turned into N."""
    out, replaced = [], 0
    for ch in seq.upper():
        if ch in ".~?":
            out.append("-")
        elif ch in ACCEPTED:
            out.append(ch)
        else:
            out.append("N")
            replaced += 1
    return "".join(out), replaced


def read_tsv(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path):
    """Ordered {name: sequence}; name is the header up to the first space."""
    seqs, name, chunks = {}, None, []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                head = line[1:].strip()
                name, chunks = (head.split()[0] if head else ""), []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def drop_empty_columns(seqs):
    """Remove columns missing in every retained sequence; return (seqs, width)."""
    if not seqs:
        return seqs, 0
    width = max(len(s) for s in seqs.values())
    seqs = {n: s.ljust(width, "-") for n, s in seqs.items()}
    keep = [i for i in range(width)
            if any(s[i] not in MISSING for s in seqs.values())]
    if len(keep) == width:
        return seqs, width
    return {n: "".join(s[i] for i in keep) for n, s in seqs.items()}, len(keep)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--alignment", required=True,
                    help="trimmed nucleotide alignment (aln.trimmed.fasta)")
    ap.add_argument("--selected", required=True,
                    help="selected.tsv, carrying label and origin columns")
    ap.add_argument("--min-refs", type=int, required=True)
    ap.add_argument("--min-aln", type=int, required=True,
                    help="minimum alignment columns (nucleotides)")
    ap.add_argument("--min-coverage", type=float, default=0.3,
                    help="minimum non-gap fraction of columns per sequence")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-eligible", required=True)
    ap.add_argument("--out-candidates", required=True)
    a = ap.parse_args()

    # 1. Origin per sequence, from the phylogeny's own selection table.
    origin = {r["label"]: (r.get("origin") or "").strip()
              for r in read_tsv(a.selected) if r.get("label")}
    seqs = read_fasta(a.alignment)
    widths = {len(s) for s in seqs.values()}
    if len(widths) > 1:
        print(f"ERROR: {a.alignment} is not aligned ({len(widths)} distinct lengths)",
              file=sys.stderr)
        return 2
    n_cols_in = widths.pop() if widths else 0

    # 2. Drop gap-dominated rows, then the columns they alone occupied.
    coverage = {n: (sum(1 for c in s if c not in MISSING) / n_cols_in) if n_cols_in else 0.0
                for n, s in seqs.items()}
    kept = {n: s for n, s in seqs.items() if coverage[n] >= a.min_coverage}
    dropped = [n for n in seqs if n not in kept]
    kept, n_cols = drop_empty_columns(kept)

    # 3. Count what is left, by origin. An unlabelled row is a reference only
    #    if its name carries the ICTV_ prefix the phylogeny writes.
    def is_ref(name):
        o = origin.get(name, "")
        return o == "reference" or (not o and name.startswith("ICTV_"))

    refs = [n for n in kept if is_ref(n)]
    cands = [n for n in kept if not is_ref(n)]
    n_refs, n_cand = len(refs), len(cands)

    eligible = int(n_refs >= a.min_refs and n_cand >= 2 and n_cols >= a.min_aln)
    reasons = []
    if n_refs < a.min_refs:
        reasons.append(f"only {n_refs} ICTV references (need {a.min_refs})")
    if n_cand < 2:
        reasons.append(f"only {n_cand} candidate sequences (need 2)")
    if n_cols < a.min_aln:
        reasons.append(f"alignment {n_cols} nt columns (need {a.min_aln})")
    if dropped:
        reasons.append(f"{len(dropped)} sequence(s) under {a.min_coverage:.0%} coverage removed")
    reason = "; ".join(reasons) if reasons else "reference density, taxa and length sufficient"

    # 4. The screened alignment keeps its gaps: it IS the matrix the scanners
    #    read. Ambiguity codes are folded to N on the way out.
    n_ambiguous = 0
    with open(a.out_fasta, "w", encoding="utf-8") as fh:
        if eligible:
            for n, s in kept.items():
                seq, replaced = normalise(s)
                n_ambiguous += replaced
                fh.write(f">{n}\n{seq}\n")
    if n_ambiguous:
        print(f"{a.family}: {n_ambiguous} ambiguous base(s) written as N",
              file=sys.stderr)
    with open(a.out_candidates, "w", encoding="utf-8") as fh:
        if eligible:
            fh.write("".join(f"{n}\n" for n in cands))
    with open(a.out_eligible, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family", "n_ictv_references", "n_candidate_sequences",
                    "alignment_length_nt", "alignment_columns_input",
                    "n_dropped_low_coverage", "dropped_sequences",
                    "min_references", "min_alignment_nt", "min_coverage",
                    "eligible", "reason"])
        w.writerow([a.family, n_refs, n_cand, n_cols, n_cols_in,
                    len(dropped), ";".join(dropped), a.min_refs, a.min_aln,
                    a.min_coverage, eligible, reason])
    print(f"{a.family}: eligible={eligible} ({reason})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
