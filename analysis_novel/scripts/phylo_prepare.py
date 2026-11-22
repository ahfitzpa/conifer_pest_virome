#!/usr/bin/env python3
"""Prepare per-family query/reference FASTA sets for ML phylogenetics.

Reads the per-family bundle written by scripts/ictv_prepare.py (candidate RdRp
domains), an ICTV reference protein FASTA (VMR MSL41), the ICTV metadata table,
and the candidate contig table.  Writes:

  queries.faa       one best RdRp/L domain per candidate contig
  refs.faa          ICTV exemplar RdRp/L proteins for this family
  combined.faa      queries + refs, unique sanitised IDs for MAFFT
  taxa_metadata.tsv taxon_id, group, family, host, accession, source

IDs are rewritten to be safe and consistent across FASTA, metadata and trees:
  candidate -> cand|<sample>|<contig>
  ICTV      -> ICTV|<accession>|<species underscored>
"""
from __future__ import annotations
import argparse, csv, re, sys
from collections import defaultdict

SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def safe(s: str) -> str:
    return SAFE.sub("_", (s or "").strip()).strip("_") or "NA"


def iter_fasta(path):
    if not path:
        return
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


def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def read_tsv(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--rdrp-faa", required=True)
    ap.add_argument("--ictv-faa", required=True)
    ap.add_argument("--ictv-meta", default="")
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out-queries", required=True)
    ap.add_argument("--out-refs", required=True)
    ap.add_argument("--out-combined", required=True)
    ap.add_argument("--out-metadata", required=True)
    a = ap.parse_args()

    fam = a.family
    # best (longest) translated domain per contig
    best = {}
    for name, seq in iter_fasta(a.rdrp_faa):
        parts = name.split("|")
        if len(parts) < 2:
            continue
        sample, contig = parts[0], parts[1]
        key = (sample, contig)
        if key not in best or len(seq) > len(best[key][1]):
            best[key] = (name, seq)

    host_by_contig = {}
    for r in read_tsv(a.candidates):
        cid = (r.get("contig_id") or "").strip()
        if cid:
            host_by_contig[cid] = (r.get("host_organism") or r.get("scientificname") or "Unknown",
                                   r.get("run", ""))

    queries, meta_rows = [], []
    for (sample, contig), (orig, seq) in sorted(best.items()):
        tid = f"cand|{safe(sample)}|{safe(contig)}"
        queries.append((tid, seq))
        host, run = host_by_contig.get(contig, ("Unknown", sample))
        meta_rows.append({"taxon_id": tid, "group": "candidate", "family": fam,
                          "host": host, "accession": run, "source": orig})

    # ICTV references: metadata family column, else family token in header
    wanted_acc = set()
    if a.ictv_meta:
        for r in read_tsv(a.ictv_meta):
            if (r.get("family") or r.get("Family") or "").strip().lower() == fam.lower():
                wanted_acc.add((r.get("accession") or r.get("Accession") or "").strip())
    refs = []
    for name, seq in iter_fasta(a.ictv_faa):
        hit = False
        if wanted_acc:
            hit = any(f"|{acc}|" in name for acc in wanted_acc if acc)
        else:
            hit = f"|{fam}|" in name
        if not hit:
            continue
        parts = name.split("|")
        acc = next((p for p in parts if re.fullmatch(r"[A-Z0-9]{5,12}", p)), "NA")
        species = parts[3] if len(parts) > 3 else "NA"
        tid = f"ICTV|{safe(acc)}|{safe(species)}"
        refs.append((tid, seq))
        meta_rows.append({"taxon_id": tid, "group": "ICTV", "family": fam,
                          "host": "NA", "accession": acc, "source": name})

    write_fasta(a.out_queries, queries)
    write_fasta(a.out_refs, refs)
    write_fasta(a.out_combined, queries + refs)
    with open(a.out_metadata, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["taxon_id", "group", "family", "host",
                                           "accession", "source"], delimiter="\t")
        w.writeheader(); w.writerows(meta_rows)

    print(f"{fam}: {len(queries)} candidate domains, {len(refs)} ICTV references, "
          f"{len(meta_rows)} taxa", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
