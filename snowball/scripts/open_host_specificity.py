#!/usr/bin/env python3
"""Open-host Serratus palmprint specificity profiler and Option-B spot-check selector.

Takes an open-host Serratus candidate table (produced by `serratus_snowball.py`
with `host_regex: ""`) and summarises per-sOTU host breadth across the public
SRA archive, classifying each sOTU into:
  - strict_beetle_specialist   (>= 90% Curculionidae/Coleoptera runs)
  - pan_arthropod_generalist   (multi-order arthropods / insectivore diet)
  - broad_low_specificity      (mixed cross-kingdom / environmental hits)

Also emits a compact Phase 2 `manifest_spotcheck.csv` selecting 1 small,
high-coverage representative SRA run per target category (e.g., non-beetle
generalist hosts + newly surfaced sister *Ips* species) so assembly-level
verification can be run cleanly through `Snakefile_phase2`.
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
LOG = logging.getLogger("open_host_specificity")

BEETLE_RE = re.compile(
    r"(ips |hylobius|dendroctonus|pissodes|tomicus|scolytus|curculionidae|coleoptera)",
    re.IGNORECASE,
)
PLANT_FUNGUS_RE = re.compile(
    r"(pinus|picea|abies|larix|pseudotsuga|ophiostoma|grosmannia|ceratocystis|endoconidiophora|fungi|plant )",
    re.IGNORECASE,
)
PREDATOR_DIET_RE = re.compile(
    r"(bat metagenome|bird metagenome|feces metagenome|gut metagenome)",
    re.IGNORECASE,
)
GENERIC_META_RE = re.compile(
    r"(metagenome|unclassified|environmental|synthetic)",
    re.IGNORECASE,
)


def _flt(val: str, default: float = 0.0) -> float:
    try:
        return float((val or "").strip())
    except ValueError:
        return default


def write_tsv_atomic(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def write_csv_atomic(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def classify_sotu(hosts: Counter) -> str:
    total = sum(hosts.values())
    if total == 0:
        return "unclassified"
    beetle_n = sum(n for h, n in hosts.items() if BEETLE_RE.search(h))
    if beetle_n / total >= 0.90:
        return "strict_beetle_specialist"
    non_insect_vert = sum(
        n for h, n in hosts.items()
        if re.search(r"(homo sapiens|mus musculus|human |soil |marine |aquatic )", h, re.I)
    )
    if non_insect_vert / total >= 0.15 and total <= 50:
        return "broad_low_specificity"
    return "pan_arthropod_generalist"


def select_spotcheck_runs(
    by_sotu_rows: Dict[str, List[Dict[str, str]]],
    sotu_classes: Dict[str, str],
    max_size_mb: float = 2500.0,
    max_total_spotchecks: int = 4,
) -> List[Dict[str, object]]:
    """Pick 1 compact, high-coverage spot-check run per generalist/outlier/new-Ips sOTU."""
    chosen: List[Dict[str, object]] = []
    seen_runs = set()

    # Priority order:
    # 1. Generalist sOTUs -> pick 1 named non-beetle insect (preferring small size & high coverage)
    # 2. Specialist sOTUs with a newly surfaced Ips species (e.g. Ips duplicatus) -> pick 1
    # 3. Broad/low-specificity outlier -> pick 1 small run
    def pick_best(rows: List[Dict[str, str]], prefer_non_beetle: bool, prefer_named: bool) -> Optional[Dict[str, str]]:
        cands = []
        for r in rows:
            rid = (r.get("run_id") or "").strip()
            if not rid or rid in seen_runs:
                continue
            size = _flt(r.get("size_mb", ""), 99999.0)
            if size <= 0 or size > max_size_mb:
                continue
            host = (r.get("scientific_name") or r.get("organism") or r.get("host_organism") or "").strip()
            if prefer_non_beetle and BEETLE_RE.search(host):
                continue
            is_named = bool(host) and not GENERIC_META_RE.search(host)
            if prefer_named and not is_named:
                continue
            cov = _flt(r.get("max_coverage", ""), 0.0)
            cands.append((not is_named, -cov, size, r))
        if not cands:
            return None
        cands.sort(key=lambda t: (t[0], t[1], t[2]))
        return cands[0][3]

    # Sort sOTUs by descending run count for deterministic selection
    ordered_sotus = sorted(by_sotu_rows.keys(), key=lambda s: -len(by_sotu_rows[s]))

    for sotu in ordered_sotus:
        if len(chosen) >= max_total_spotchecks:
            break
        cls = sotu_classes.get(sotu, "")
        rows = by_sotu_rows[sotu]
        picked: Optional[Dict[str, str]] = None
        if cls == "pan_arthropod_generalist":
            picked = pick_best(rows, prefer_non_beetle=True, prefer_named=True) or pick_best(
                rows, prefer_non_beetle=True, prefer_named=False
            )
        elif cls == "broad_low_specificity":
            picked = pick_best(rows, prefer_non_beetle=True, prefer_named=False)
        elif cls == "strict_beetle_specialist":
            # Check for a novel sister Ips species (e.g. Ips duplicatus) not yet in Round 1/2
            dup_rows = [
                r for r in rows
                if "duplicatus" in (r.get("scientific_name") or r.get("organism") or r.get("host_organism") or "").lower()
            ]
            if dup_rows:
                picked = pick_best(dup_rows, prefer_non_beetle=False, prefer_named=True)
        if picked is not None:
            rid = (picked.get("run_id") or "").strip()
            seen_runs.add(rid)
            host = (
                picked.get("scientific_name")
                or picked.get("organism")
                or picked.get("host_organism")
                or "Unknown"
            ).strip()
            bp = (picked.get("bioproject") or picked.get("bio_project") or "").strip()
            layout = (picked.get("library_layout") or "paired").strip().lower()
            platform = (picked.get("platform") or "illumina").strip().lower()
            notes = (
                f"spotcheck_sotu={sotu};"
                f"max_cov={picked.get('max_coverage', '')};"
                f"size_mb={picked.get('size_mb', '')}"
            )
            rec: Dict[str, object] = {col: "" for col in PHASE2_MANIFEST_COLS}
            rec.update(
                {
                    "sample_id": rid,
                    "family": "OpenHost_SpotCheck",
                    "source": "sra",
                    "run_id": rid,
                    "platform": platform,
                    "layout": layout,
                    "host_species": host,
                    "notes": notes,
                    "bioproject": bp,
                    "scientificname": host,
                    "scenario": f"spotcheck_{cls}",
                    "host_status": "open_host_spotcheck",
                }
            )
            chosen.append(rec)
    return chosen


PHASE2_MANIFEST_COLS = [
    "sample_id", "family", "source", "run_id", "fq1", "fq2", "platform", "layout",
    "barcode", "library_id", "batch", "condition", "replicate", "host_species",
    "collection_date", "notes", "md5_fq1", "md5_fq2", "adapter", "genome_size",
    "expected_coverage", "bioproject", "scientificname", "ref_genome_path",
    "scenario", "host_status",
]


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--open-host-tsv", required=True, type=Path)
    ap.add_argument("--out-specificity-tsv", required=True, type=Path)
    ap.add_argument("--out-spotcheck-manifest", required=True, type=Path)
    ap.add_argument("--max-spotcheck-size-mb", type=float, default=2500.0)
    ap.add_argument("--max-spotcheck-runs", type=int, default=4)
    args = ap.parse_args(argv)

    by_sotu_hosts: Dict[str, Counter] = defaultdict(Counter)
    by_sotu_rows: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    with args.open_host_tsv.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            host = (
                row.get("scientific_name")
                or row.get("organism")
                or row.get("host_organism")
                or row.get("host")
                or "Unknown"
            ).strip()
            sotus_raw = row.get("matched_sotus") or row.get("sotus") or ""
            for s in sotus_raw.split(";"):
                s = s.strip()
                if s:
                    by_sotu_hosts[s][host] += 1
                    by_sotu_rows[s].append(row)

    spec_rows: List[Dict[str, object]] = []
    sotu_classes: Dict[str, str] = {}
    for sotu, hosts in sorted(by_sotu_hosts.items(), key=lambda kv: -sum(kv[1].values())):
        total = sum(hosts.values())
        beetle_n = sum(n for h, n in hosts.items() if BEETLE_RE.search(h))
        plant_fung_n = sum(n for h, n in hosts.items() if PLANT_FUNGUS_RE.search(h))
        pred_n = sum(n for h, n in hosts.items() if PREDATOR_DIET_RE.search(h))
        cls = classify_sotu(hosts)
        sotu_classes[sotu] = cls
        top_hosts = "; ".join(f"{h} ({n})" for h, n in hosts.most_common(6))
        spec_rows.append(
            {
                "sotu": sotu,
                "n_candidate_runs": total,
                "ecological_class": cls,
                "pct_curculionidae_coleoptera": f"{100.0 * beetle_n / total:.1f}" if total else "0.0",
                "pct_predator_diet_metagenome": f"{100.0 * pred_n / total:.1f}" if total else "0.0",
                "pct_conifer_or_fungus": f"{100.0 * plant_fung_n / total:.1f}" if total else "0.0",
                "top_hosts": top_hosts,
            }
        )

    write_tsv_atomic(
        args.out_specificity_tsv,
        [
            "sotu", "n_candidate_runs", "ecological_class",
            "pct_curculionidae_coleoptera", "pct_predator_diet_metagenome",
            "pct_conifer_or_fungus", "top_hosts",
        ],
        spec_rows,
    )

    spotchecks = select_spotcheck_runs(
        by_sotu_rows,
        sotu_classes,
        max_size_mb=args.max_spotcheck_size_mb,
        max_total_spotchecks=args.max_spotcheck_runs,
    )
    write_csv_atomic(
        args.out_spotcheck_manifest,
        PHASE2_MANIFEST_COLS,
        spotchecks,
    )
    LOG.info(
        "Open-host specificity summary: %d sOTU(s) classified, %d spot-check run(s) selected -> %s",
        len(spec_rows), len(spotchecks), args.out_spotcheck_manifest,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
