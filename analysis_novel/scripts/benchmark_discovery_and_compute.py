#!/usr/bin/env python3
"""Stage 17.4: Benchmark divergence-bin discovery sensitivity and Snakemake compute cost.

Pure-Python, zero-OpenBLAS script that produces:
  - divergence_bin_sensitivity.tsv  (recovery stratified by ICTV aa ID bins:
                                     <30%, 30-40%, 40-60%, >60%, family_unresolved)
  - compute_cost_summary.tsv        (per-stage wall time, CPU-hours, and peak RSS
                                     aggregated directly from Snakemake benchmarks/)
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import os
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Sequence

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
LOG = logging.getLogger("benchmark_discovery_and_compute")

BINS = ["<30%", "30-40%", "40-60%", ">60%", "family_unresolved"]


def _flt(val: str, default: float = float("nan")) -> float:
    try:
        v = float((val or "").strip())
        return v if math.isfinite(v) else default
    except ValueError:
        return default


def assign_divergence_bin(row: Dict[str, str]) -> str:
    fam = (row.get("consensus_family") or "").strip().lower()
    reason = (row.get("selection_reason") or "").strip().lower()
    if fam in {"", "unresolved", "unassigned", "unknown", "na"} or "unresolved" in reason:
        return "family_unresolved"
    ident = _flt(row.get("nearest_identity_pct", ""))
    if math.isnan(ident):
        return "family_unresolved"
    if ident < 30.0:
        return "<30%"
    if ident < 40.0:
        return "30-40%"
    if ident <= 60.0:
        return "40-60%"
    return ">60%"


def write_tsv_atomic(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def summarize_divergence_bins(
    esm_tsvs: Sequence[Path],
    pebblescout_tsv: Optional[Path],
) -> List[Dict[str, object]]:
    counts = {
        b: {"hmm_gated": 0, "validated_rdrp": 0, "viral_ns3": 0, "rt_purged": 0, "host_purged": 0, "uncertain": 0}
        for b in BINS
    }
    merged_rows: Dict[tuple, Dict[str, str]] = {}
    for esm_tsv in esm_tsvs:
        if not esm_tsv.is_file():
            LOG.warning("Skipping missing ESM validation TSV: %s", esm_tsv)
            continue
        with esm_tsv.open("r", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                run_id = (row.get("run_id") or "").strip()
                contig = (row.get("contig") or "").strip()
                cid = (row.get("candidate_id") or "").strip()
                key = (run_id, contig) if (run_id and contig) else (cid, "")
                merged_rows[key] = row

    for row in merged_rows.values():
        b = assign_divergence_bin(row)
        counts[b]["hmm_gated"] += 1
        verdict = (row.get("structural_verdict") or "").strip().lower()
        target = (row.get("rdrp_target") or "").lower()
        if verdict in {"validated_rdrp_high", "validated_rdrp_moderate"}:
            if ("8gzq" in target or "8gzr" in target) and target.endswith("_b"):
                counts[b]["viral_ns3"] += 1
            else:
                counts[b]["validated_rdrp"] += 1
        elif verdict == "rt_or_non_rdrp_flag":
            counts[b]["rt_purged"] += 1
        elif verdict in {"low_confidence_fold", "no_structural_hit", "fold_failed"}:
            counts[b]["uncertain"] += 1
        else:
            counts[b]["host_purged"] += 1

    # Optional PebbleScout summary check
    ps_max_cov = "NA"
    ps_hits = "NA"
    if pebblescout_tsv and pebblescout_tsv.is_file():
        covs = []
        n_sig = 0
        with pebblescout_tsv.open("r", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cov = _flt(r.get("PctCoverage") or r.get("max_coverage_pct") or r.get("coverage_pct") or "", 0.0)
                covs.append(cov)
                if cov >= 10.0:
                    n_sig += 1
        if covs:
            ps_max_cov = f"{max(covs):.2f}%"
            ps_hits = str(n_sig)

    rows: List[Dict[str, object]] = []
    for b in BINS:
        c = counts[b]
        n = c["hmm_gated"]
        prec_after = (c["validated_rdrp"] / n) if n else 0.0
        rej_rate = ((c["rt_purged"] + c["host_purged"]) / n) if n else 0.0
        rows.append(
            {
                "aa_identity_bin": b,
                "hmm_gated_domains": n,
                "validated_viral_rdrp": c["validated_rdrp"],
                "validated_viral_ns3": c["viral_ns3"],
                "rt_purged": c["rt_purged"],
                "host_other_purged": c["host_purged"],
                "uncertain": c["uncertain"],
                "hmm_raw_rdrp_yield": f"{prec_after:.4f}",
                "structural_rejection_rate": f"{rej_rate:.4f}",
                "pebblescout_sig_hits_ge10pct": ps_hits,
                "pebblescout_max_cov": ps_max_cov,
            }
        )
    return rows


def summarize_snakemake_benchmarks(bench_dirs: Sequence[Path]) -> List[Dict[str, object]]:
    by_rule: Dict[str, List[Dict[str, float]]] = defaultdict(list)
    for bdir in bench_dirs:
        if not bdir.is_dir():
            continue
        for tsv in sorted(bdir.rglob("*.tsv")):
            rule_name = tsv.parent.name
            with tsv.open("r", encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    wall_s = _flt(row.get("s", ""), 0.0)
                    cpu_s = _flt(row.get("cpu_time", ""), wall_s)
                    rss_mb = _flt(row.get("max_rss", ""), 0.0)
                    by_rule[rule_name].append(
                        {"wall_s": wall_s, "cpu_s": cpu_s, "rss_mb": rss_mb}
                    )

    rows: List[Dict[str, object]] = []
    for rule, recs in sorted(by_rule.items()):
        walls = [r["wall_s"] for r in recs]
        cpus = [r["cpu_s"] for r in recs]
        rss = [r["rss_mb"] for r in recs]
        rows.append(
            {
                "rule_name": rule,
                "n_jobs": len(recs),
                "median_wall_min": f"{median(walls) / 60.0:.2f}",
                "max_wall_min": f"{max(walls) / 60.0:.2f}",
                "total_wall_hours": f"{sum(walls) / 3600.0:.3f}",
                "total_cpu_hours": f"{sum(cpus) / 3600.0:.3f}",
                "peak_rss_gb": f"{max(rss) / 1024.0:.2f}",
            }
        )
    return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--esm-validation",
        required=True,
        action="append",
        type=Path,
        help="Stage 15 esm_rdrp_validation.tsv (repeatable).",
    )
    ap.add_argument("--pebblescout-tsv", type=Path, default=None)
    ap.add_argument(
        "--benchmark-dir",
        action="append",
        type=Path,
        default=[],
        help="Snakemake benchmarks/ directory (repeatable).",
    )
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    div_rows = summarize_divergence_bins(args.esm_validation, args.pebblescout_tsv)
    write_tsv_atomic(
        args.out_dir / "divergence_bin_sensitivity.tsv",
        [
            "aa_identity_bin", "hmm_gated_domains", "validated_viral_rdrp",
            "validated_viral_ns3", "rt_purged", "host_other_purged", "uncertain",
            "hmm_raw_rdrp_yield", "structural_rejection_rate",
            "pebblescout_sig_hits_ge10pct", "pebblescout_max_cov",
        ],
        div_rows,
    )

    comp_rows = summarize_snakemake_benchmarks(args.benchmark_dir)
    write_tsv_atomic(
        args.out_dir / "compute_cost_summary.tsv",
        [
            "rule_name", "n_jobs", "median_wall_min", "max_wall_min",
            "total_wall_hours", "total_cpu_hours", "peak_rss_gb",
        ],
        comp_rows,
    )
    LOG.info(
        "Discovery & compute benchmark complete: %d divergence bins, %d Snakemake rules -> %s",
        len(div_rows), len(comp_rows), args.out_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
