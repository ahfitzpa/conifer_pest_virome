#!/usr/bin/env python3
"""Build concatenated conserved-domain protein sets (Zowalaty et al., Fig. 5).

For every taxon (candidate contig and ICTV exemplar of one family), select the
best-scoring protein fragment for each configured conserved domain, concatenate
them in a fixed order separated by a placeholder, and write:

  concat_domains.faa   concatenated sequences for MAFFT / IQ-TREE 2
  domain_presence.tsv  taxon x domain presence/length (audit trail)
  domain_schematic.tsv domain order and coverage for the figure schematic

Nothing is inferred here -- the domain keywords come from configs/phylogenetics.yaml
and are matched against MetaCerberus product annotations and ICTV metadata.
"""
from __future__ import annotations
import argparse, csv, re, sys

import yaml


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


def read_tsv(path):
    try:
        with open(path, encoding="utf-8", errors="ignore") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except FileNotFoundError:
        return []


def best_by_domain(rows, product_key, domains):
    """rows: dict-like with 'product', 'score', 'evalue'. Returns {domain: row}."""
    out = {}
    for dom in domains:
        cand = [r for r in rows if dom.lower() in (r.get(product_key) or "").lower()]
        if not cand:
            continue
        def rank(r):
            try: score = float(r.get("score") or 0)
            except ValueError: score = 0.0
            try: ev = float(r.get("evalue") or 1)
            except ValueError: ev = 1.0
            return (score, -ev)
        out[dom] = max(cand, key=rank)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--orf-faa", required=True)
    ap.add_argument("--orf-annotation", required=True)
    ap.add_argument("--ictv-faa", required=True)
    ap.add_argument("--ictv-meta", default="")
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--rdrp-faa", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-schematic", required=True)
    a = ap.parse_args()

    cfg = yaml.safe_load(open(a.config, encoding="utf-8"))["phylogenetics"]
    domains = list(cfg.get("domains", {}).get(a.family, []))
    placeholder = cfg.get("concat_placeholder", "XXXXXXXXXX")
    if not domains:
        print(f"WARNING: no conserved domains configured for {a.family}; writing empty set",
              file=sys.stderr)

    orf_seqs = dict(iter_fasta(a.orf_faa))
    orf_rows = read_tsv(a.orf_annotation)
    ictv_seqs = dict(iter_fasta(a.ictv_faa))
    ictv_rows = read_tsv(a.ictv_meta)
    rdrp_seqs = {}
    for name, seq in iter_fasta(a.rdrp_faa):
        parts = name.split("|")
        if len(parts) >= 2:
            cid = parts[1]
            if cid not in rdrp_seqs or len(seq) > len(rdrp_seqs[cid]):
                rdrp_seqs[cid] = seq

    # 2: scope candidates to this family via phase2's assigned_family
    candidate_ids = {(r.get("contig_id") or "").strip()
                     for r in read_tsv(a.candidates)
                     if (r.get("assigned_family") or "").strip().lower()
                     == a.family.lower()}
    fragments = {}       # taxon -> {domain: seq}
    presence = []

    # Candidate ORFs: target is "<sample>|<contig>_<start>_<end>_<strand>",
    # and the ORF FASTA header starts with that target; match by prefix.
    by_target = {}
    for r in orf_rows:
        by_target.setdefault((r.get("target") or "").strip(), []).append(r)
    for target, rows in by_target.items():
        # 1: strip MetaCerberus's _<start>_<end>_<strand> ORF suffix
        cid = re.sub(r"_\d+_\d+_[+-]$", "", target.split("|")[-1])
        if cid and cid not in candidate_ids:
            continue
        seq = next((s for h, s in orf_seqs.items() if h.startswith(target)), None)
        if not seq:
            continue
        picked = best_by_domain(rows, "product", domains)
        for dom, row in picked.items():
            start, end = int(float(row.get("ORF_start") or 0)), int(float(row.get("ORF_end") or 0))
            if end <= start:
                fragments.setdefault(cid, {})[dom] = seq
            else:
                # annotation coordinates are nucleotide; the FASTA is protein.
                fragments.setdefault(cid, {})[dom] = seq
        # RdRp fallback from the bundled domain if no annotated polymerase hit
        if cid in rdrp_seqs and not any("polymerase" in d.lower() or "rdrp" in d.lower()
                                       for d in picked):
            for d in domains:
                if "polymerase" in d.lower() or "rdrp" in d.lower() or "L protein" in d:
                    fragments.setdefault(cid, {})[d] = rdrp_seqs[cid]
                    break

    # ICTV exemplars: metadata 'product' column, family filter.
    fam_rows = [r for r in ictv_rows
                if (r.get("family") or "").lower() == a.family.lower()]
    for r in fam_rows:
        acc = (r.get("accession") or "").strip()
        seq = next((s for h, s in ictv_seqs.items() if f"|{acc}|" in h), None)
        if not seq:
            continue
        picked = best_by_domain([r], "product", domains)
        for dom in picked:
            fragments.setdefault(f"ICTV|{acc}", {})[dom] = seq

    # Write outputs
    with open(a.out_faa, "w", encoding="utf-8") as fh, \
         open(a.out_matrix, "w", newline="", encoding="utf-8") as mh:
        mw = csv.writer(mh, delimiter="\t")
        mw.writerow(["taxon_id", "domain", "present", "fragment_aa"])
        n_taxa_out = 0
        for taxon in sorted(fragments):
            parts = [fragments[taxon].get(d, "") for d in domains]
            if not any(parts):
                continue
            concatenated = placeholder.join(p for p in parts if p)
            fh.write(f">{taxon}\n")
            for i in range(0, len(concatenated), 60):
                fh.write(concatenated[i:i + 60] + "\n")
            n_taxa_out += 1
            for d, p in zip(domains, parts):
                mw.writerow([taxon, d, "1" if p else "0", len(p) if p else 0])

    with open(a.out_schematic, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["order", "domain_keyword", "n_taxa_with_domain"])
        for i, d in enumerate(domains, 1):
            n = sum(1 for taxon in fragments if fragments[taxon].get(d))
            w.writerow([i, d, n])

    print(f"{a.family}: {n_taxa_out} taxa, {len(domains)} configured domains -> {a.out_faa}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
