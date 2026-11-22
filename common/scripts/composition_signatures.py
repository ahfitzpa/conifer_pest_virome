#!/usr/bin/env python3
"""RSCU and dinucleotide odds ratios per genome, from in-frame CDS.

Both are host-adaptation signals: codon usage and CpG/UpA depletion track the
host's own composition, so they separate genomes that merely look similar from
genomes under similar pressure.

RSCU = observed codon count / mean count across its synonymous family, so 1.0 is
no preference. Odds ratio = f(XY) / (f(X) * f(Y)), so 1.0 is what independence
would give, and CpG depletion shows up below 1.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter

from Bio.Data import CodonTable


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def synonymous_families(code: int) -> dict[str, list[str]]:
    # 1: codon -> its synonymous family, stops excluded (they are not chosen)
    table = CodonTable.unambiguous_dna_by_id[code]
    by_aa: dict[str, list[str]] = {}
    for codon, aa in table.forward_table.items():
        by_aa.setdefault(aa, []).append(codon)
    return {c: sorted(by_aa[aa]) for aa, codons in by_aa.items() for c in codons}


def rscu(seq: str, families: dict[str, list[str]]) -> dict[str, float]:
    # 2: count only clean in-frame codons that belong to a family
    counts = Counter()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if codon in families:
            counts[codon] += 1
    out: dict[str, float] = {}
    for codon, family in families.items():
        total = sum(counts[c] for c in family)
        if total:
            out[codon] = counts[codon] / (total / len(family))
    return out


def dinucleotide_odds(seq: str) -> dict[str, float]:
    # 3: odds ratio over the whole CDS, counting overlapping dinucleotides
    bases = "ACGT"
    seq = "".join(c for c in seq if c in bases)
    if len(seq) < 2:
        return {}
    mono = Counter(seq)
    di = Counter(seq[i:i + 2] for i in range(len(seq) - 1))
    n_mono, n_di = len(seq), len(seq) - 1
    out = {}
    for x in bases:
        for y in bases:
            fx, fy = mono[x] / n_mono, mono[y] / n_mono
            if fx and fy:
                out[x + y] = (di[x + y] / n_di) / (fx * fy)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cds", required=True, help="in-frame CDS FASTA, one per genome")
    p.add_argument("--out-rscu", required=True)
    p.add_argument("--out-dinuc", required=True)
    p.add_argument("--genetic-code", type=int, default=1)
    p.add_argument("--min-codons", type=int, default=100,
                   help="skip a genome with fewer clean codons than this")
    args = p.parse_args()

    families = synonymous_families(args.genetic_code)
    codons = sorted(families)
    dinucs = [x + y for x in "ACGT" for y in "ACGT"]

    rscu_rows, dinuc_rows, skipped = [], [], []
    for name, seq in read_fasta(args.cds):
        seq = seq.upper().replace("U", "T")
        usable = sum(1 for i in range(0, len(seq) - 2, 3)
                     if seq[i:i + 3] in families)
        if usable < args.min_codons:
            skipped.append((name, usable))
            continue
        r = rscu(seq, families)
        d = dinucleotide_odds(seq)
        rscu_rows.append({"genome": name, "n_codons": usable,
                          **{c: f"{r[c]:.4f}" for c in codons if c in r}})
        dinuc_rows.append({"genome": name, "length_nt": len(seq),
                           **{k: f"{d[k]:.4f}" for k in dinucs if k in d}})

    if not rscu_rows:
        log("[composition] no genome had enough clean codons")
        return 1

    os.makedirs(os.path.dirname(os.path.abspath(args.out_rscu)) or ".",
                exist_ok=True)
    with open(args.out_rscu, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "n_codons"] + codons,
                           delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rscu_rows)
    with open(args.out_dinuc, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "length_nt"] + dinucs,
                           delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(dinuc_rows)

    cpg = [float(r["CG"]) for r in dinuc_rows if r.get("CG")]
    log(f"[composition] {len(rscu_rows)} genome(s); "
        f"median CpG odds {sorted(cpg)[len(cpg)//2]:.3f}" if cpg else
        f"[composition] {len(rscu_rows)} genome(s)")
    for name, n in skipped:
        log(f"[composition] skipped {name}: {n} clean codons")
    log(f"[composition] -> {args.out_rscu}, {args.out_dinuc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
