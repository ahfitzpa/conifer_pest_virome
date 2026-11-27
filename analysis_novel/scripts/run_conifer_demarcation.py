#!/usr/bin/env python3
"""Batch ICTV pairwise-identity demarcation across all families under results/phylogeny/.

Stdlib-only (no BioPython/Conda required). Computes exact ungapped-overlap
pairwise identity between every candidate and every ICTV reference directly from
each family's trimmed multiple-sequence alignment (aln.trimmed.fasta or .faa) and
metadata (selected.tsv or taxa_metadata.tsv).

Outputs:
  - <out-dir>/<fam>/pairwise_identity.tsv
  - <out-dir>/<fam>/demarcation_summary.tsv
  - <out-dir>/all_families_demarcation_summary.tsv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    buf: list[str] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf).upper()
                cur_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(buf).upper()
    return seqs


def is_nucleotide(seqs: dict[str, str]) -> bool:
    sample = "".join(list(seqs.values())[:10]).replace("-", "").replace("N", "").replace("X", "")
    if not sample:
        return True
    acgt = sum(1 for c in sample if c in "ACGTU")
    return (acgt / len(sample)) > 0.85


def pairwise_aln_identity(s1: str, s2: str) -> tuple[float, float, int]:
    """Return (identity_pct, query_coverage_pct, overlap_cols) over aligned columns."""
    n = min(len(s1), len(s2))
    q_nongap = 0
    overlap = 0
    matches = 0
    for i in range(n):
        c1, c2 = s1[i], s2[i]
        if c1 not in "-.?":
            q_nongap += 1
            if c2 not in "-.?":
                overlap += 1
                if c1 == c2 and c1 not in "NX":
                    matches += 1
    ident = (100.0 * matches / overlap) if overlap else 0.0
    cov = (100.0 * overlap / q_nongap) if q_nongap else 0.0
    return round(ident, 2), round(cov, 2), overlap


def novelty_tier(ident: float, passes_cov: bool) -> str:
    if not passes_cov:
        return "low_overlap (<50% cov)"
    if ident < 50.0:
        return "putative_novel_genus_or_higher (<50%)"
    if ident < 70.0:
        return "putative_novel_species_or_genus (50-70%)"
    if ident < 90.0:
        return "putative_novel_species (70-90%)"
    return "known_species_like (>=90%)"


def load_meta(path: Path) -> tuple[list[dict[str, str]], list[str], list[str]]:
    rows: list[dict[str, str]] = []
    cands: list[str] = []
    refs: list[str] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tid = (r.get("label") or r.get("taxon_id") or r.get("id") or "").strip()
            if not tid:
                continue
            r["_id"] = tid
            rows.append(r)
            grp = (r.get("origin") or r.get("group") or "").strip().lower()
            if grp == "candidate":
                cands.append(tid)
            elif grp in {"reference", "ictv", "ref"}:
                refs.append(tid)
    return rows, cands, refs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phylogeny-dir", default="results/phylogeny", type=Path)
    ap.add_argument("--out-dir", default="results/demarcation", type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    master_rows: list[dict[str, object]] = []
    n_fams = 0

    for fam_dir in sorted(args.phylogeny_dir.iterdir()):
        if not fam_dir.is_dir():
            continue
        fam = fam_dir.name
        aln_path = next(
            (p for p in (fam_dir / "aln.trimmed.fasta", fam_dir / "aln.trimmed.faa") if p.exists()),
            None,
        )
        meta_path = next(
            (p for p in (fam_dir / "selected.tsv", fam_dir / "taxa_metadata.tsv") if p.exists()),
            None,
        )
        if not aln_path or not meta_path:
            continue

        seqs = read_fasta(aln_path)
        meta_rows, cands, refs = load_meta(meta_path)
        meta_by_id = {r["_id"]: r for r in meta_rows}
        if not cands or not refs:
            continue

        seq_type = "nt" if is_nucleotide(seqs) else "aa"
        fam_out = args.out_dir / fam
        fam_out.mkdir(parents=True, exist_ok=True)

        pw_rows: list[dict[str, object]] = []
        best_rows: list[dict[str, object]] = []

        for q in cands:
            if q not in seqs:
                continue
            best: dict[str, object] | None = None
            for ref in refs:
                if ref not in seqs:
                    continue
                ident, cov, ov = pairwise_aln_identity(seqs[q], seqs[ref])
                rec = {
                    "query": q,
                    "reference": ref,
                    "identity_pct": ident,
                    "query_coverage_pct": cov,
                    "aligned_cols": ov,
                    "seq_type": seq_type,
                    "passes_50pct_coverage": int(cov >= 50.0),
                }
                pw_rows.append(rec)
                if best is None or (rec["passes_50pct_coverage"], rec["identity_pct"], rec["aligned_cols"]) > (
                    best["passes_50pct_coverage"], best["identity_pct"], best["aligned_cols"]
                ):
                    best = rec

            if best is not None:
                m = meta_by_id.get(q, {})
                tier = novelty_tier(float(best["identity_pct"]), bool(best["passes_50pct_coverage"]))
                srow = {
                    "family": fam,
                    "candidate": q,
                    "host": m.get("host") or m.get("scientificname") or "",
                    "bioproject": m.get("bioproject") or "",
                    "checkv_quality": m.get("checkv_quality") or "",
                    "checkv_completeness": m.get("checkv_completeness") or "",
                    "nearest_ictv_ref": best["reference"],
                    "nearest_identity_pct": best["identity_pct"],
                    "query_coverage_pct": best["query_coverage_pct"],
                    "aligned_cols": best["aligned_cols"],
                    "seq_type": seq_type,
                    "novelty_tier": tier,
                }
                best_rows.append(srow)
                master_rows.append(srow)

        if pw_rows:
            with (fam_out / "pairwise_identity.tsv").open("w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(pw_rows[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(pw_rows)
        if best_rows:
            with (fam_out / "demarcation_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(best_rows[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(best_rows)
            n_fams += 1
            print(f"[demarcation] {fam}: {len(best_rows)} candidates x {len(refs)} ICTV refs ({seq_type})")

    if master_rows:
        master_path = args.out_dir / "all_families_demarcation_summary.tsv"
        with master_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(master_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(master_rows)
        print(f"Wrote {master_path} ({len(master_rows)} candidates across {n_fams} families)")
    else:
        print("No families with both candidates and ICTV references found under", args.phylogeny_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
