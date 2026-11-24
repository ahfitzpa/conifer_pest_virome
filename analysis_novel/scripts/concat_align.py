#!/usr/bin/env python3
"""Align each conserved-domain partition, then concatenate with gap padding.

MAFFT strips the gaps in its input and realigns across partition boundaries, so
a bare concatenation cannot hold the partitions in register: a candidate
carrying polymerase, nucleoprotein and glycoprotein and a reference carrying
only the polymerase end up with non-homologous columns. Each partition is
therefore aligned on its own and a taxon missing one is padded with that
partition's aligned width in gaps.

  concat.aln.faa          padded concatenation, every row the same length
  concat.aln.trimmed.faa  trimAl of the above, for IQ-TREE 2
  concat.partitions.tsv   column ranges, describing concat.aln.faa (pre-trimAl)
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile


def iter_fasta(path):
    name, chunks = None, []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def tool(name):
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        path = os.path.join(prefix, "bin", name)
        if os.path.exists(path):
            return path
    return shutil.which(name) or name


def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def align(records, strategy, threads, log_fh):
    """One partition through MAFFT; a lone sequence is already its alignment."""
    if len(records) < 2:
        return dict(records)
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "part.faa")
        write_fasta(src, records)
        cmd = [tool("mafft")] + strategy.split() + ["--thread", str(threads), src]
        run = subprocess.run(cmd, capture_output=True, text=True, check=True)
        log_fh.write(run.stderr)
        dst = os.path.join(td, "part.aln.faa")
        with open(dst, "w", encoding="utf-8") as fh:
            fh.write(run.stdout)
        return dict(iter_fasta(dst))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--partitions", required=True)
    ap.add_argument("--out-aln", required=True)
    ap.add_argument("--out-trimmed", required=True)
    ap.add_argument("--out-ranges", required=True)
    ap.add_argument("--mafft-strategy", default="--auto")
    ap.add_argument("--trimal-strategy", default="-automated1")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()

    # 1: group fragments by partition tag written by concat_domains.py
    parts = {}
    for name, seq in iter_fasta(a.fragments):
        tag, _, taxon = name.partition("|")
        if taxon:
            parts.setdefault(tag, []).append((taxon, seq))

    # 2: keep the configured partition order
    order = []
    with open(a.partitions, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            order.append((f"part{r['part']}", r["domain_keyword"]))

    taxa = sorted({t for recs in parts.values() for t, _ in recs})
    if not taxa:
        raise SystemExit(f"no fragments in {a.fragments}")

    # 3: align each partition separately, recording its aligned width
    blocks, ranges, start = {}, [], 1
    with open(a.out_aln + ".log", "w", encoding="utf-8") as log_fh:
        for tag, keyword in order:
            recs = parts.get(tag, [])
            if not recs:
                continue
            aligned = align(recs, a.mafft_strategy, a.threads, log_fh)
            width = len(next(iter(aligned.values())))
            if any(len(s) != width for s in aligned.values()):
                raise SystemExit(f"{tag}: mafft returned ragged lengths")
            blocks[tag] = (aligned, width)
            ranges.append([tag, keyword, start, start + width - 1, len(recs)])
            start += width
        log_fh.write(f"partitions {len(blocks)}, taxa {len(taxa)}\n")
    if not blocks:
        raise SystemExit("every configured partition was empty")

    # 4: concatenate in order, gaps for a taxon with no fragment there
    records = []
    for taxon in taxa:
        records.append((taxon, "".join(
            blocks[tag][0].get(taxon, "-" * blocks[tag][1])
            for tag, _ in order if tag in blocks)))
    write_fasta(a.out_aln, records)

    with open(a.out_ranges, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["part", "domain_keyword", "start", "end", "n_taxa"])
        w.writerows(ranges)

    # 5: trimAl, as the single-marker path does, on the padded concatenation
    with open(a.out_trimmed + ".log", "w", encoding="utf-8") as log_fh:
        cmd = [tool("trimal")] + a.trimal_strategy.split() + [
            "-in", a.out_aln, "-out", a.out_trimmed]
        run = subprocess.run(cmd, capture_output=True, text=True, check=True)
        log_fh.write(run.stdout + run.stderr)

    cols = sum(w for _, w in blocks.values())
    print(f"{len(taxa)} taxa, {len(blocks)} partitions, {cols} columns "
          f"-> {a.out_aln}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
