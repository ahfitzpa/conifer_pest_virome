#!/usr/bin/env python3
"""Pick one ORF per genome from Prodigal output: the longest, inside a window.

For Hepeviridae, ORF1 (which carries the RdRp) spans roughly 5 kb of a 7.2 kb
genome, so in any genome assembled far enough to be worth analysing the longest
ORF is ORF1. That is a heuristic, so it is enforced rather than assumed: a
genome whose longest ORF falls outside the window is dropped and reported, not
analysed as if it were ORF1.

Prodigal names ORFs `<contig>_<n>` and records coordinates in the header, so the
parent genome is the name with the final `_<n>` removed.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys


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
                name, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def parent_genome(orf_id: str) -> str:
    # 1: strip Prodigal's trailing _<n> to get the contig it came from
    return re.sub(r"_\d+$", "", orf_id)


def partial_flag(header: str) -> str:
    # 2: Prodigal records partial=00 for a complete ORF, 10/01/11 when an end
    # runs off the contig. A partial ORF1 is still usable but worth recording.
    m = re.search(r"partial=(\d\d)", header)
    return m.group(1) if m else ""


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--faa", required=True, help="Prodigal protein output")
    p.add_argument("--ffn", required=True, help="Prodigal nucleotide output")
    p.add_argument("--out-faa", required=True)
    p.add_argument("--out-ffn", required=True)
    p.add_argument("--out-table", required=True)
    p.add_argument("--min-nt", type=int, default=1500)
    p.add_argument("--max-nt", type=int, default=8000)
    p.add_argument("--min-taxa", type=int, default=4)
    args = p.parse_args()

    cds = {}
    for header, seq in read_fasta(args.ffn):
        cds[header.split()[0]] = seq
    proteins = {}
    headers = {}
    for header, seq in read_fasta(args.faa):
        orf_id = header.split()[0]
        proteins[orf_id] = seq
        headers[orf_id] = header
    log(f"[orf1] {len(proteins)} ORF(s) predicted, {len(cds)} with CDS")

    # 3: group ORFs by genome, then take the longest CDS in each
    by_genome: dict[str, list[tuple[int, str]]] = {}
    for orf_id in proteins:
        if orf_id not in cds:
            continue
        by_genome.setdefault(parent_genome(orf_id), []).append(
            (len(cds[orf_id]), orf_id))

    rows, kept = [], []
    for genome in sorted(by_genome):
        entries = sorted(by_genome[genome], key=lambda e: (-e[0], e[1]))
        length, orf_id = entries[0]
        in_window = args.min_nt <= length <= args.max_nt
        rows.append({
            "genome": genome,
            "orf_id": orf_id,
            "cds_nt": length,
            "protein_aa": len(proteins[orf_id]),
            "n_orfs_in_genome": len(entries),
            "partial": partial_flag(headers[orf_id]),
            "status": "selected" if in_window else "dropped",
            "reason": "" if in_window else
                      f"longest ORF {length} nt outside "
                      f"{args.min_nt}-{args.max_nt} nt",
        })
        if in_window:
            kept.append((genome, orf_id))

    # 4: write one record per genome, named by the genome so the ids match the
    # tree tips rather than Prodigal's ORF numbering
    for out_path, source in ((args.out_faa, proteins), (args.out_ffn, cds)):
        os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
        with open(out_path, "w") as fh:
            for genome, orf_id in kept:
                seq = source[orf_id]
                fh.write(f">{genome}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")

    cols = ["genome", "orf_id", "cds_nt", "protein_aa", "n_orfs_in_genome",
            "partial", "status", "reason"]
    with open(args.out_table, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    dropped = len(rows) - len(kept)
    lengths = sorted(r["cds_nt"] for r in rows if r["status"] == "selected")
    if lengths:
        log(f"[orf1] selected lengths: min {lengths[0]}, "
            f"median {lengths[len(lengths)//2]}, max {lengths[-1]} nt")
    log(f"[orf1] {len(kept)} genome(s) selected, {dropped} dropped "
        f"-> {args.out_faa}")
    if len(kept) < args.min_taxa:
        log(f"[orf1] only {len(kept)} genome(s); need {args.min_taxa}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
