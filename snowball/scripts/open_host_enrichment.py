#!/usr/bin/env python3
"""Compute sampling-effort-normalized host enrichment (Fisher's exact test & odds ratios)
for open-host Serratus sOTU detections across 3.8M public SRA runs.

Addresses sampling-bias concerns by normalizing raw SRA run counts per sOTU against
total archived runs in Coleoptera/Curculionidae (and Ips/Hylobius) vs. all other
public SRA runs, and stratifying hits by library source (transcriptomic vs. metagenomic).

Pure-Python statistics using log-gamma hypergeometric exact probabilities (zero SciPy
or OpenBLAS dependency).
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
LOG = logging.getLogger("open_host_enrichment")

BEETLE_RE = re.compile(
    r"(ips |hylobius|dendroctonus|scolyt|curculionid|tomicus|pityogenes|orthotomicus|coleoptera|beetle|weevil)",
    re.IGNORECASE,
)
IPS_RE = re.compile(r"\bips\s+[a-z]+", re.IGNORECASE)
HYLOBIUS_RE = re.compile(r"\bhylobius\b", re.IGNORECASE)
PREDATOR_DIET_RE = re.compile(
    r"(chiroptera|myotis|pipistrellus|eptesicus|vespertilio|rhinolophus|bat |bird |aves|gut metagenome|feces|fecal|stool)",
    re.IGNORECASE,
)

# Conservative Serratus palm_sra2 / srarun archive denominators (3.8M runs indexed in Serratus)
# Used as fallback if --taxon-totals-tsv is not supplied and live SQL is disabled.
DEFAULT_ARCHIVE_TOTALS = {
    "total_sra_runs": 3_837_755,
    "coleoptera_runs": 14_500,
    "curculionidae_runs": 2_400,
    "ips_runs": 320,
    "hylobius_runs": 65,
}


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_log_pmf(a: int, r1: int, c1: int, n: int) -> float:
    return _log_comb(c1, a) + _log_comb(n - c1, r1 - a) - _log_comb(n, r1)


def fisher_exact_enrichment(a: int, b: int, c: int, d: int) -> Tuple[float, float, float, float]:
    """Compute (odds_ratio, ci95_low, ci95_high, p_value_greater) for a 2x2 table:
        [[a, b],
         [c, d]]
    Uses Haldane-Anscombe (+0.5) continuity correction when any cell is 0 for OR/CI,
    and exact log-space hypergeometric summation for one-sided enrichment p-value.
    """
    a, b, c, d = max(0, int(a)), max(0, int(b)), max(0, int(c)), max(0, int(d))
    aa, bb, cc, dd = float(a), float(b), float(c), float(d)
    if min(a, b, c, d) == 0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
    odds_ratio = (aa * dd) / (bb * cc)
    se_log_or = math.sqrt(1.0 / aa + 1.0 / bb + 1.0 / cc + 1.0 / dd)
    log_or = math.log(odds_ratio)
    ci_low = math.exp(log_or - 1.96 * se_log_or)
    ci_high = math.exp(log_or + 1.96 * se_log_or)

    r1 = a + b
    c1 = a + c
    n = a + b + c + d
    if n == 0 or r1 == 0 or c1 == 0:
        return odds_ratio, ci_low, ci_high, 1.0

    upper = min(r1, c1)
    log_probs = [_hypergeom_log_pmf(k, r1, c1, n) for k in range(a, upper + 1)]
    if not log_probs:
        return odds_ratio, ci_low, ci_high, 1.0
    m = max(log_probs)
    p_val = math.exp(m) * sum(math.exp(lp - m) for lp in log_probs)
    p_val = min(1.0, max(0.0, p_val))
    return odds_ratio, ci_low, ci_high, p_val


def load_taxon_totals(path: Optional[Path]) -> Dict[str, int]:
    totals = dict(DEFAULT_ARCHIVE_TOTALS)
    if path is None or not path.is_file():
        return totals
    with path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row.get("taxon_key") or "").strip().lower()
            val = (row.get("n_runs") or "").strip()
            if key and val.isdigit():
                totals[key] = int(val)
    return totals


def query_serratus_denominators(host: str, port: int, dbname: str, user: str, password: str) -> Dict[str, int]:
    """Query public Serratus PostgreSQL RDS srarun table for exact archive denominators."""
    import psycopg2  # optional dependency only when --query-serratus is passed

    conn = psycopg2.connect(host=host, port=port, dbname=dbname, user=user, password=password, connect_timeout=15)
    totals = dict(DEFAULT_ARCHIVE_TOTALS)
    queries = {
        "total_sra_runs": "SELECT COUNT(*) FROM srarun;",
        "coleoptera_runs": (
            "SELECT COUNT(*) FROM srarun WHERE scientific_name ~* "
            "'(ips |hylobius|dendroctonus|scolyt|curculionid|tomicus|pityogenes|orthotomicus|coleoptera|tribolium|anoplophora|agrilus|leptinotarsa|sitophilus)';"
        ),
        "curculionidae_runs": (
            "SELECT COUNT(*) FROM srarun WHERE scientific_name ~* "
            "'(ips |hylobius|dendroctonus|scolyt|curculionid|tomicus|pityogenes|orthotomicus|sitophilus|rhynchophorus)';"
        ),
        "ips_runs": "SELECT COUNT(*) FROM srarun WHERE scientific_name ~* '^ips ';",
        "hylobius_runs": "SELECT COUNT(*) FROM srarun WHERE scientific_name ~* '^hylobius ';",
    }
    try:
        with conn.cursor() as cur:
            for key, sql in queries.items():
                cur.execute(sql)
                res = cur.fetchone()
                if res and res[0] is not None:
                    totals[key] = int(res[0])
    finally:
        conn.close()
    return totals


def write_tsv_atomic(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--open-host-tsv", required=True, type=Path)
    ap.add_argument("--taxon-totals-tsv", type=Path, default=None)
    ap.add_argument("--query-serratus", action="store_true", help="Query live Serratus RDS for exact srarun denominators.")
    ap.add_argument("--db-host", default="serratus-public.c8602223.us-east-1.rds.amazonaws.com")
    ap.add_argument("--db-port", type=int, default=5432)
    ap.add_argument("--db-name", default="serratus")
    ap.add_argument("--db-user", default="serratus_public")
    ap.add_argument("--db-password", default="serratus_public")
    ap.add_argument("--out-enrichment-tsv", required=True, type=Path)
    ap.add_argument("--out-denominators-tsv", type=Path, default=None)
    args = ap.parse_args(argv)

    if args.query_serratus:
        try:
            totals = query_serratus_denominators(args.db_host, args.db_port, args.db_name, args.db_user, args.db_password)
            denom_source = "serratus_rds_live"
        except Exception as exc:
            LOG.warning("Live Serratus denominator query failed (%s); falling back to cached/default totals.", exc)
            totals = load_taxon_totals(args.taxon_totals_tsv)
            denom_source = "default_archive_totals"
    else:
        totals = load_taxon_totals(args.taxon_totals_tsv)
        denom_source = "taxon_totals_tsv" if (args.taxon_totals_tsv and args.taxon_totals_tsv.is_file()) else "default_archive_totals"

    if args.out_denominators_tsv:
        write_tsv_atomic(
            args.out_denominators_tsv,
            ["taxon_key", "n_runs", "source"],
            [{"taxon_key": k, "n_runs": v, "source": denom_source} for k, v in totals.items()],
        )

    by_sotu_hosts: Dict[str, Counter] = defaultdict(Counter)
    by_sotu_sources: Dict[str, Counter] = defaultdict(Counter)

    with args.open_host_tsv.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            host = (
                row.get("scientific_name")
                or row.get("organism")
                or row.get("host_organism")
                or row.get("host")
                or "Unknown"
            ).strip()
            lib_src = (row.get("library_source") or row.get("source") or "UNKNOWN").strip().upper()
            sotus_raw = row.get("matched_sotus") or row.get("sotus") or ""
            for s in sotus_raw.split(";"):
                s = s.strip()
                if s:
                    by_sotu_hosts[s][host] += 1
                    by_sotu_sources[s][lib_src] += 1

    n_total_sra = totals["total_sra_runs"]
    n_coleoptera = totals["coleoptera_runs"]
    n_non_coleoptera = max(1, n_total_sra - n_coleoptera)

    out_rows: List[Dict[str, object]] = []
    for sotu, hosts in sorted(by_sotu_hosts.items(), key=lambda kv: -sum(kv[1].values())):
        total_hits = sum(hosts.values())
        beetle_hits = sum(n for h, n in hosts.items() if BEETLE_RE.search(h))
        non_beetle_hits = total_hits - beetle_hits
        ips_hits = sum(n for h, n in hosts.items() if IPS_RE.search(h))
        hyl_hits = sum(n for h, n in hosts.items() if HYLOBIUS_RE.search(h))
        pred_hits = sum(n for h, n in hosts.items() if PREDATOR_DIET_RE.search(h))

        c_neg_beetle = max(0, n_coleoptera - beetle_hits)
        d_neg_other = max(0, n_non_coleoptera - non_beetle_hits)
        or_val, ci_lo, ci_hi, p_val = fisher_exact_enrichment(beetle_hits, non_beetle_hits, c_neg_beetle, d_neg_other)

        rate_beetle_per_1k = 1000.0 * beetle_hits / max(1, n_coleoptera)
        rate_other_per_1k = 1000.0 * non_beetle_hits / max(1, n_non_coleoptera)

        srcs = by_sotu_sources[sotu]
        n_tx = srcs.get("TRANSCRIPTOMIC", 0)
        n_mg = srcs.get("METAGENOMIC", 0) + srcs.get("METATRANSCRIPTOMIC", 0)

        beetle_frac = beetle_hits / total_hits if total_hits else 0.0
        if beetle_frac >= 0.90 and p_val < 1e-5:
            tier = "beetle_restricted_in_archives"
        elif beetle_frac >= 0.25 and p_val < 1e-3:
            tier = "beetle_enriched"
        elif non_beetle_hits > 0 and (pred_hits / total_hits >= 0.01 or total_hits >= 100):
            tier = "pan_arthropod_generalist"
        else:
            tier = "broad_or_uncertain"

        out_rows.append(
            {
                "sotu": sotu,
                "n_total_runs": total_hits,
                "n_beetle_runs": beetle_hits,
                "n_non_beetle_runs": non_beetle_hits,
                "n_ips_runs": ips_hits,
                "n_hylobius_runs": hyl_hits,
                "pct_beetle": f"{100.0 * beetle_frac:.1f}",
                "rate_beetle_per_1k_sra": f"{rate_beetle_per_1k:.3f}",
                "rate_non_beetle_per_1k_sra": f"{rate_other_per_1k:.5f}",
                "beetle_odds_ratio": f"{or_val:.1f}",
                "or_ci95_low": f"{ci_lo:.1f}",
                "or_ci95_high": f"{ci_hi:.1f}",
                "fisher_exact_p": f"{p_val:.3e}",
                "pct_transcriptomic": f"{100.0 * n_tx / total_hits:.1f}" if total_hits else "0.0",
                "pct_metagenomic": f"{100.0 * n_mg / total_hits:.1f}" if total_hits else "0.0",
                "pct_predator_diet": f"{100.0 * pred_hits / total_hits:.1f}" if total_hits else "0.0",
                "host_association_tier": tier,
                "denominator_source": denom_source,
            }
        )

    write_tsv_atomic(
        args.out_enrichment_tsv,
        [
            "sotu", "n_total_runs", "n_beetle_runs", "n_non_beetle_runs",
            "n_ips_runs", "n_hylobius_runs", "pct_beetle",
            "rate_beetle_per_1k_sra", "rate_non_beetle_per_1k_sra",
            "beetle_odds_ratio", "or_ci95_low", "or_ci95_high", "fisher_exact_p",
            "pct_transcriptomic", "pct_metagenomic", "pct_predator_diet",
            "host_association_tier", "denominator_source",
        ],
        out_rows,
    )
    LOG.info("Open-host enrichment complete: %d sOTU(s) analyzed -> %s", len(out_rows), args.out_enrichment_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
