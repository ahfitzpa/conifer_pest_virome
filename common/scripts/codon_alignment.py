#!/usr/bin/env python3
"""Build an in-frame codon alignment from a protein alignment plus its CDS.

The analysis pathways align amino acids, but dN/dS (codeml, HyPhy FUBAR/MEME)
needs nucleotides in frame. This threads each CDS onto its aligned protein, so
every protein column becomes one codon column (pal2nal equivalent).

Writes the codon alignment and a per-sequence QC report.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

GAP = "-"
STOPS_OK_AT_END = True


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def norm_id(name: str, strip_version: bool, strip_suffix: str) -> str:
    # 1: normalise an id so protein and CDS names can meet
    key = name.split()[0]
    if strip_suffix:
        key = re.sub(strip_suffix, "", key)
    if strip_version:
        key = re.sub(r"\.\d+$", "", key)
    return key


def translate(cds: str, table: int) -> str:
    # 2: translate the in-frame part only (ambiguous codons become X)
    usable = len(cds) - (len(cds) % 3)
    return str(Seq(cds[:usable]).translate(table=table, to_stop=False, cds=False))


def compare(aa_ref: str, aa_obs: str) -> tuple[int, int]:
    # 3: count comparable positions and mismatches; X or * on either side is a free pass
    n = mism = 0
    for a, b in zip(aa_ref, aa_obs):
        if a in ("X", "*") or b in ("X", "*"):
            continue
        n += 1
        if a.upper() != b.upper():
            mism += 1
    return n, mism


def thread(aa_aligned: str, cds: str) -> str:
    # 4: one aligned residue consumes one codon; a gap emits a gap codon
    out, i = [], 0
    for residue in aa_aligned:
        if residue == GAP or residue == ".":
            out.append(GAP * 3)
            continue
        codon = cds[i:i + 3]
        if len(codon) < 3:
            codon = codon + "N" * (3 - len(codon))
        out.append(codon)
        i += 3
    return "".join(out)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--aa-aln", required=True, help="aligned protein FASTA")
    p.add_argument("--cds", required=True, help="unaligned nucleotide CDS FASTA")
    p.add_argument("--out", required=True, help="output codon alignment FASTA")
    p.add_argument("--report", help="per-sequence QC report TSV")
    p.add_argument("--genetic-code", type=int, default=1)
    p.add_argument("--max-mismatch-frac", type=float, default=0.02,
                   help="allowed protein/CDS disagreement per sequence")
    p.add_argument("--on-problem", choices=("fail", "drop", "warn"), default="drop",
                   help="what to do with a sequence that fails QC")
    p.add_argument("--min-taxa", type=int, default=4,
                   help="fail if fewer sequences survive (dN/dS needs a few)")
    p.add_argument("--strip-version", action="store_true",
                   help="ignore a trailing .1 when matching ids")
    p.add_argument("--strip-suffix", default="",
                   help="regex removed from ids before matching, e.g. '_ORF[12]$'")
    args = p.parse_args()

    if args.genetic_code not in CodonTable.generic_by_id:
        sys.exit(f"unknown genetic code {args.genetic_code}")

    aa_records = list(SeqIO.parse(args.aa_aln, "fasta"))
    if not aa_records:
        sys.exit(f"no sequences in {args.aa_aln}")
    widths = {len(r.seq) for r in aa_records}
    if len(widths) != 1:
        sys.exit(f"{args.aa_aln} is not an alignment: widths {sorted(widths)}")

    # 5: index the CDS by normalised id, refusing ambiguous duplicates
    cds_by_id: dict[str, str] = {}
    for rec in SeqIO.parse(args.cds, "fasta"):
        key = norm_id(rec.id, args.strip_version, args.strip_suffix)
        if key in cds_by_id and str(rec.seq).upper() != cds_by_id[key]:
            sys.exit(f"duplicate CDS id with different sequence: {key}")
        cds_by_id[key] = str(rec.seq).upper().replace("U", "T")

    rows, kept, problems = [], [], 0
    for rec in aa_records:
        key = norm_id(rec.id, args.strip_version, args.strip_suffix)
        aa_aligned = str(rec.seq).upper()
        aa_plain = aa_aligned.replace(GAP, "").replace(".", "").rstrip("*")
        row = {"sequence": rec.id, "matched_cds": key in cds_by_id,
               "n_residues": len(aa_plain), "n_codons": 0, "cds_len": 0,
               "frame_ok": "", "internal_stop": "", "n_compared": 0,
               "n_mismatch": 0, "mismatch_frac": "", "status": "", "reason": ""}

        if key not in cds_by_id:
            row.update(status="dropped", reason="no CDS with this id")
            rows.append(row)
            problems += 1
            continue

        cds = cds_by_id[key]
        # 6: drop a trailing stop codon only when the CDS is already in frame,
        # so a frameshift is never silently repaired into a valid length
        aa_full = translate(cds, args.genetic_code)
        if STOPS_OK_AT_END and len(cds) % 3 == 0 and aa_full.endswith("*"):
            cds = cds[:-3]
            aa_full = translate(cds, args.genetic_code)
        row["cds_len"] = len(cds)
        row["frame_ok"] = len(cds) % 3 == 0
        row["internal_stop"] = "*" in aa_full[:-1] if aa_full else False

        needed = len(aa_plain) * 3
        if len(cds) < needed:
            row.update(status="dropped",
                       reason=f"CDS too short: {len(cds)} nt for {len(aa_plain)} residues")
            rows.append(row)
            problems += 1
            continue

        n_cmp, mism = compare(aa_plain, aa_full[: len(aa_plain)])
        frac = (mism / n_cmp) if n_cmp else 1.0
        row.update(n_compared=n_cmp, n_mismatch=mism, mismatch_frac=f"{frac:.4f}")

        if not row["frame_ok"] or frac > args.max_mismatch_frac or row["internal_stop"]:
            reason = ("frame not a multiple of 3" if not row["frame_ok"] else
                      "internal stop codon" if row["internal_stop"] else
                      f"protein/CDS disagreement {frac:.1%} > {args.max_mismatch_frac:.1%}")
            if args.on_problem == "fail":
                sys.exit(f"{rec.id}: {reason}")
            row.update(status="dropped" if args.on_problem == "drop" else "kept_with_warning",
                       reason=reason)
            problems += 1
            if args.on_problem == "drop":
                rows.append(row)
                continue
        else:
            row["status"] = "kept"

        codon_aln = thread(aa_aligned, cds[:needed])
        row["n_codons"] = len(codon_aln) // 3
        rows.append(row)
        kept.append((rec.id, rec.description, codon_aln))

    # 7: every kept row must be the same width, or downstream tools read it wrong
    kept_widths = {len(s) for _, _, s in kept}
    if len(kept_widths) > 1:
        sys.exit(f"internal error: ragged codon alignment widths {sorted(kept_widths)}")
    if len(kept) < args.min_taxa:
        log(f"[codon] only {len(kept)} sequence(s) survived QC (need {args.min_taxa})")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        for name, desc, seq in kept:
            header = desc if desc and desc != name else name
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    if args.report:
        cols = ["sequence", "matched_cds", "n_residues", "cds_len", "n_codons",
                "frame_ok", "internal_stop", "n_compared", "n_mismatch",
                "mismatch_frac", "status", "reason"]
        with open(args.report, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    width = next(iter(kept_widths), 0)
    log(f"[codon] {len(kept)}/{len(aa_records)} sequences, {width // 3} codon columns, "
        f"{problems} QC problem(s) -> {args.out}")
    return 0 if len(kept) >= args.min_taxa else 1


if __name__ == "__main__":
    raise SystemExit(main())
