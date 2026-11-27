#!/usr/bin/env python3
"""Summarize family-unresolved RdRp dark-matter contigs and generate Manuscript Table 1
(top complete/near-complete viral candidates with CheckV %, coverage, length,
multi-run prevalence, and nearest ICTV identity).
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

FOCAL_FAMILIES = {
    "Iflaviridae", "Polycipiviridae", "Amarillovirales",
    "Dicistroviridae", "Partitiviridae", "Phenuiviridae",
}


def as_float(v: str | None) -> float | None:
    try:
        x = float(str(v).strip())
        return x
    except (TypeError, ValueError):
        return None


def parse_len_cov(contig_id: str) -> tuple[float | None, float | None]:
    ml = re.search(r"length_(\d+)", contig_id)
    mc = re.search(r"cov_([0-9.]+)", contig_id)
    length_kb = int(ml.group(1)) / 1000.0 if ml else None
    cov = float(mc.group(1)) if mc else None
    return length_kb, cov


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master-metadata", required=True, type=Path)
    ap.add_argument("--demarcation-summary", type=Path, default=Path("results/demarcation/all_families_demarcation_summary.tsv"))
    ap.add_argument("--out-dir", type=Path, default=Path("results/manuscript_checks"))
    ap.add_argument("--min-checkv", type=float, default=50.0, help="Min CheckV completeness for Table 1")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Load demarcation lookup if available
    demarc: dict[str, dict[str, str]] = {}
    if args.demarcation_summary.exists():
        with args.demarcation_summary.open(encoding="utf-8", errors="ignore") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                demarc[r["candidate"]] = r

    unresolved_by_host: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    fam_runs: dict[tuple[str, str], set[str]] = defaultdict(set)
    all_rows: list[dict[str, str]] = []

    with args.master_metadata.open(newline="", encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh):
            all_rows.append(r)
            fam = (r.get("consensus_family") or "").strip()
            host = (r.get("scientificname") or r.get("family_context") or "Unknown").strip()
            bp = (r.get("bioproject") or "Unknown").strip()
            run = (r.get("run_id") or "").strip()
            if fam and fam.lower() not in {"unclassified", "unknown", "na"}:
                fam_runs[(fam, host)].add(run)
            else:
                ev = as_float(r.get("rdrpcatch_evalue"))
                lkb, cov = parse_len_cov(r.get("contig_id") or "")
                unresolved_by_host[(host, bp)].append({
                    "has_rdrp": ev is not None,
                    "length_kb": lkb or 0.0,
                    "evalue": ev,
                })

    # 1. Unresolved RdRp dark-matter summary table
    unres_path = args.out_dir / "unresolved_rdrp_summary.tsv"
    with unres_path.open("w", newline="", encoding="utf-8") as fh:
        fields = [
            "host", "bioproject", "total_unresolved_contigs", "rdrp_positive_contigs",
            "median_rdrp_length_kb", "max_rdrp_length_kb", "rdrp_contigs_ge_2kb",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        tot_u = tot_r = 0
        for (host, bp), items in sorted(unresolved_by_host.items()):
            rdrps = [x for x in items if x["has_rdrp"]]
            lens = [float(x["length_kb"]) for x in rdrps]
            tot_u += len(items)
            tot_r += len(rdrps)
            w.writerow({
                "host": host,
                "bioproject": bp,
                "total_unresolved_contigs": len(items),
                "rdrp_positive_contigs": len(rdrps),
                "median_rdrp_length_kb": round(statistics.median(lens), 2) if lens else 0.0,
                "max_rdrp_length_kb": round(max(lens), 2) if lens else 0.0,
                "rdrp_contigs_ge_2kb": sum(1 for x in lens if x >= 2.0),
            })
    print(f"Wrote {unres_path} ({tot_u} unresolved contigs, {tot_r} RdRp-positive)")

    # 2. Manuscript Table 1: Top viral candidates (>= min_checkv or focal families with >= 2 kb)
    t1_rows: list[dict[str, object]] = []
    for r in all_rows:
        fam = (r.get("consensus_family") or "").strip()
        if not fam or fam.lower() in {"unclassified", "unknown", "na"}:
            continue
        comp = as_float(r.get("checkv_completeness")) or 0.0
        cid = (r.get("contig_id") or "").strip()
        run = (r.get("run_id") or "").strip()
        host = (r.get("scientificname") or r.get("family_context") or "").strip()
        bp = (r.get("bioproject") or "").strip()
        lkb, cov = parse_len_cov(cid)
        if comp < args.min_checkv and not (fam in FOCAL_FAMILIES and (lkb or 0.0) >= 3.0):
            continue
        # Match demarcation entry by run_nXX or exact contig_id
        dm = demarc.get(cid) or next(
            (v for k, v in demarc.items() if run and k.startswith(run) and (v.get("family") == fam)),
            {},
        )
        t1_rows.append({
            "family": fam,
            "host": host,
            "bioproject": bp,
            "run_id": run,
            "contig_id": cid,
            "length_kb": round(lkb, 2) if lkb is not None else "",
            "spades_cov": round(cov, 1) if cov is not None else "",
            "checkv_quality": (r.get("checkv_quality") or "").strip(),
            "checkv_completeness_pct": round(comp, 1),
            "host_run_prevalence": len(fam_runs.get((fam, host), set())),
            "nearest_ictv_ref": dm.get("nearest_ictv_ref", ""),
            "nearest_identity_pct": dm.get("nearest_identity_pct", ""),
            "novelty_tier": dm.get("novelty_tier", ""),
        })

    t1_rows.sort(key=lambda x: (x["family"] not in FOCAL_FAMILIES, x["family"], -float(x["checkv_completeness_pct"])))
    t1_path = args.out_dir / "table1_top_candidates.tsv"
    if t1_rows:
        with t1_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(t1_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(t1_rows)
        print(f"Wrote {t1_path} ({len(t1_rows)} top candidate rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
