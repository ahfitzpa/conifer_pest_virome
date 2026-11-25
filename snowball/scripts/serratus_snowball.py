#!/usr/bin/env python3
"""Query Serratus RDS for SRA runs sharing RdRp sOTUs with validated novel viruses.

Combines:
  1. DIAMOND blastp hits of structurally validated novel RdRps against PalmDB
     (`palmdb_hits.tsv`) -> mapped to Serratus `sotu` IDs (`u...`).
  2. Optional seed-run sOTUs from `palm_sra2` for the runs that yielded the
     validated novel viruses (`include_seed_run_sotus: true`).

Queries the public read-only Serratus PostgreSQL RDS (`palm_sra2` + `srarun`)
using parameterized SQL, excludes already-processed runs from
`--existing-manifest`, enforces Phase 2 compatibility (paired-end Illumina,
size cap, cell-culture exclusion, optional host regex), caps runs per
BioProject for diversity, and writes:
  - `--out-candidates`: ranked TSV of candidate SRA runs & provenance
  - `--out-manifest`:   26-column Phase 2 manifest (`manifest_snowball.csv`)
  - `--out-combined`:   existing manifest + new snowball runs merged cleanly
"""

from __future__ import annotations

import argparse
import csv
import io
import logging
import os
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

__version__ = "1.0.0"
LOG = logging.getLogger("serratus_snowball")

MANIFEST_COLS = [
    "sample_id", "family", "source", "run_id", "fq1", "fq2", "platform", "layout",
    "barcode", "library_id", "batch", "condition", "replicate", "host_species",
    "collection_date", "notes", "md5_fq1", "md5_fq2", "adapter", "genome_size",
    "expected_coverage", "bioproject", "scientificname", "ref_genome_path",
    "scenario", "host_status",
]

CELL_CULTURE_RE = (
    r"(cell[ _-]?culture|cell[ _-]?line|passage|serial[ _-]?passag|"
    r"vero|huh7|hek293|bhk|c6/36|s2 cells|in vitro|synthetic)"
)

SOTU_RE = re.compile(r"\b(u\d+)\b", re.IGNORECASE)


def atomic_write_text(dest: Path, text: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=dest.parent, prefix=f".{dest.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp_path, dest)
    except BaseException:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def atomic_write_csv(dest: Path, cols: list[str], rows: list[dict[str, str]]) -> None:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({c: str(row.get(c, "")) for c in cols})
    atomic_write_text(dest, buf.getvalue())


def extract_sotu_id(target_header: str) -> str:
    """Extract a PalmDB sOTU ID (e.g., u12345) from a DIAMOND subject ID."""
    m = SOTU_RE.search(target_header)
    if m:
        return m.group(1).lower()
    return re.split(r"[|_\s]", target_header.strip())[0]


def parse_palmdb_hits(
    hits_tsv: Path | None,
    min_pident: float = 40.0,
    max_evalue: float = 1e-5,
    max_sotus_per_query: int = 10,
) -> dict[str, set[str]]:
    """Return mapping sotu -> set of matched novel query candidate_ids."""
    sotu_to_queries: dict[str, set[str]] = defaultdict(set)
    if hits_tsv is None or not hits_tsv.is_file() or hits_tsv.stat().st_size == 0:
        return sotu_to_queries

    per_query_count: dict[str, int] = defaultdict(int)
    with hits_tsv.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("qseqid\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            qseqid, sseqid = parts[0].strip(), parts[1].strip()
            try:
                pident = float(parts[2])
                evalue = float(parts[4])
            except ValueError:
                continue
            if pident < min_pident or evalue > max_evalue:
                continue
            if per_query_count[qseqid] >= max_sotus_per_query:
                continue
            sotu = extract_sotu_id(sseqid)
            if sotu:
                sotu_to_queries[sotu].add(qseqid)
                per_query_count[qseqid] += 1

    return sotu_to_queries


def load_existing_manifest(manifest_path: Path | None) -> tuple[set[str], list[dict[str, str]]]:
    if manifest_path is None or not manifest_path.is_file() or manifest_path.stat().st_size == 0:
        return set(), []
    runs: set[str] = set()
    rows: list[dict[str, str]] = []
    with manifest_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            clean = {k: (v or "").strip() for k, v in row.items() if k is not None}
            rows.append(clean)
            for col in ("run_id", "sample_id"):
                if clean.get(col):
                    runs.add(clean[col])
    return runs, rows


def load_seed_runs(validated_tsv: Path | None) -> set[str]:
    runs: set[str] = set()
    if validated_tsv is None or not validated_tsv.is_file() or validated_tsv.stat().st_size == 0:
        return runs
    with validated_tsv.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rid = (row.get("run_id") or "").strip()
            if rid:
                runs.add(rid)
    return runs


def select_diverse_runs(
    candidate_rows: list[dict[str, object]],
    existing_runs: set[str],
    max_runs_per_bioproject: int = 4,
    max_total_runs: int = 12,
) -> list[dict[str, object]]:
    """Filter out existing runs, rank by shared sOTUs + coverage, cap per BioProject."""
    filtered = [r for r in candidate_rows if str(r["run_id"]) not in existing_runs]
    filtered.sort(
        key=lambda r: (
            -int(r.get("n_matched_sotus", 0)),
            -float(r.get("max_coverage", 0.0)),
            -float(r.get("max_pident", 0.0)),
            str(r["run_id"]),
        )
    )

    selected: list[dict[str, object]] = []
    bp_counts: dict[str, int] = defaultdict(int)
    for row in filtered:
        if len(selected) >= max_total_runs:
            break
        bp = str(row.get("bioproject") or "UNKNOWN").strip() or "UNKNOWN"
        if max_runs_per_bioproject > 0 and bp_counts[bp] >= max_runs_per_bioproject:
            continue
        selected.append(row)
        bp_counts[bp] += 1
    return selected


def build_phase2_manifest_rows(
    selected_rows: list[dict[str, object]],
    scenario_label: str = "snowball_rdrp",
    default_family: str = "Snowball_RdRp",
) -> list[dict[str, str]]:
    """Convert selected snowball candidate runs into exact 26-column Phase 2 manifest rows."""
    records: list[dict[str, str]] = []
    for r in selected_rows:
        rid = str(r["run_id"]).strip()
        sci = str(r.get("scientific_name") or "").strip()
        bp = str(r.get("bioproject") or "").strip()
        layout = str(r.get("library_layout") or "PAIRED").strip().lower()
        platform = str(r.get("platform") or "ILLUMINA").strip().lower()
        notes = (
            f"snowball_sotus={r.get('n_matched_sotus', '')};"
            f"max_cov={r.get('max_coverage', '')};"
            f"queries={r.get('matched_queries', '')}"
        )
        rec = {col: "" for col in MANIFEST_COLS}
        rec.update({
            "sample_id": rid,
            "family": default_family,
            "source": "sra",
            "run_id": rid,
            "platform": platform,
            "layout": layout,
            "host_species": sci,
            "notes": notes,
            "bioproject": bp,
            "scientificname": sci,
            "scenario": scenario_label,
            "host_status": "snowball_candidate",
        })
        records.append(rec)
    return records


def query_serratus_rds(
    db_cfg: dict[str, object],
    sotu_to_queries: dict[str, set[str]],
    seed_runs: set[str],
    include_seed_run_sotus: bool,
    host_regex: str,
    min_coverage: int,
    min_run_pident: float,
    paired_only: bool,
    illumina_only: bool,
    exclude_cell_culture: bool,
    max_size_mb: float,
) -> list[dict[str, object]]:
    """Query public read-only Serratus PostgreSQL RDS for matching SRA runs."""
    import psycopg

    conn_cfg = dict(db_cfg)
    if "name" in conn_cfg:
        conn_cfg["dbname"] = conn_cfg.pop("name")
    conn_cfg.setdefault("password", os.environ.get("SERRATUS_DB_PASSWORD", "serratus"))

    all_sotus: dict[str, set[str]] = defaultdict(set)
    for sotu, qset in sotu_to_queries.items():
        all_sotus[sotu].update(qset)

    with psycopg.connect(**conn_cfg) as conn:
        with conn.cursor() as cur:
            if include_seed_run_sotus and seed_runs:
                cur.execute(
                    """
                    SELECT DISTINCT ps.sotu, ps.run_id
                    FROM palm_sra2 ps
                    WHERE ps.run_id = ANY(%s)
                      AND ps.coverage >= %s
                    """,
                    (sorted(seed_runs), max(1, min_coverage)),
                )
                seed_hits = cur.fetchall()
                LOG.info("Seed runs yielded %d (sotu, run_id) links in Serratus", len(seed_hits))
                for sotu, rid in seed_hits:
                    all_sotus[str(sotu)].add(f"seed:{rid}")

            if not all_sotus:
                LOG.warning("No sOTUs found from PalmDB hits or seed runs; returning 0 candidates.")
                return []

            sotu_list = sorted(all_sotus.keys())
            LOG.info("Querying Serratus palm_sra2 + srarun across %d unique sOTU(s)", len(sotu_list))

            clauses = [
                "ps.sotu = ANY(%s)",
                "ps.coverage >= %s",
                "ps.percent_identity >= %s",
            ]
            params: list[object] = [sotu_list, min_coverage, min_run_pident]

            if paired_only:
                clauses.append("UPPER(COALESCE(sr.library_layout, '')) = 'PAIRED'")
            if illumina_only:
                clauses.append("UPPER(COALESCE(sr.platform, '')) = 'ILLUMINA'")
            if max_size_mb > 0:
                clauses.append("(sr.size_mb IS NULL OR sr.size_mb <= %s)")
                params.append(max_size_mb)
            if host_regex:
                clauses.append("sr.scientific_name ~* %s")
                params.append(host_regex)
            if exclude_cell_culture:
                for col in ("sr.scientific_name", "sr.library_name", "sr.sample_name"):
                    clauses.append(f"({col} IS NULL OR {col} !~* %s)")
                    params.append(CELL_CULTURE_RE)

            sql = f"""
                SELECT ps.run_id,
                       COALESCE(sr.bio_project, '')      AS bioproject,
                       COALESCE(sr.scientific_name, '')  AS scientific_name,
                       COALESCE(sr.library_layout, '')   AS library_layout,
                       COALESCE(sr.platform, '')         AS platform,
                       COALESCE(sr.library_strategy, '') AS library_strategy,
                       COALESCE(sr.size_mb, 0)           AS size_mb,
                       ps.sotu,
                       max(ps.coverage)                  AS max_cov,
                       max(ps.percent_identity)          AS max_pid
                FROM palm_sra2 ps
                JOIN srarun sr ON sr.run = ps.run_id
                WHERE {' AND '.join(clauses)}
                GROUP BY ps.run_id, sr.bio_project, sr.scientific_name,
                         sr.library_layout, sr.platform, sr.library_strategy,
                         sr.size_mb, ps.sotu
            """
            cur.execute(sql, tuple(params))
            raw_rows = cur.fetchall()

    by_run: dict[str, dict[str, object]] = {}
    for run_id, bp, sci, layout, plat, strat, size_mb, sotu, cov, pid in raw_rows:
        rec = by_run.get(run_id)
        if rec is None:
            rec = {
                "run_id": str(run_id),
                "bioproject": str(bp),
                "scientific_name": str(sci),
                "library_layout": str(layout),
                "platform": str(plat),
                "library_strategy": str(strat),
                "size_mb": float(size_mb or 0.0),
                "sotus": set(),
                "queries": set(),
                "max_coverage": 0.0,
                "max_pident": 0.0,
            }
            by_run[run_id] = rec
        rec["sotus"].add(str(sotu))
        rec["queries"].update(all_sotus.get(str(sotu), set()))
        rec["max_coverage"] = max(float(rec["max_coverage"]), float(cov or 0.0))
        rec["max_pident"] = max(float(rec["max_pident"]), float(pid or 0.0))

    candidates: list[dict[str, object]] = []
    for rec in by_run.values():
        sotus_sorted = sorted(rec["sotus"])
        queries_sorted = sorted(rec["queries"])
        candidates.append({
            "run_id": rec["run_id"],
            "bioproject": rec["bioproject"],
            "scientific_name": rec["scientific_name"],
            "library_layout": rec["library_layout"],
            "platform": rec["platform"],
            "library_strategy": rec["library_strategy"],
            "size_mb": f"{float(rec['size_mb']):.1f}",
            "n_matched_sotus": len(sotus_sorted),
            "max_coverage": f"{float(rec['max_coverage']):.1f}",
            "max_pident": f"{float(rec['max_pident']):.1f}",
            "matched_sotus": ";".join(sotus_sorted[:20]),
            "matched_queries": ";".join(queries_sorted[:10]),
        })
    return candidates


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, type=Path, help="Path to snowball.yaml")
    ap.add_argument("--palmdb-hits", default=None, type=Path, help="DIAMOND hits TSV against PalmDB")
    ap.add_argument("--validated-tsv", default=None, type=Path, help="Extracted validated RdRps TSV")
    ap.add_argument("--existing-manifest", default=None, type=Path, help="Existing Phase 2 manifest CSV")
    ap.add_argument("--out-candidates", required=True, type=Path, help="Ranked snowball candidate runs TSV")
    ap.add_argument("--out-manifest", required=True, type=Path, help="Phase 2 manifest CSV for new snowball runs")
    ap.add_argument("--out-combined", required=True, type=Path, help="Merged existing + snowball Phase 2 manifest CSV")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv: list[str] | None = None) -> int:
    import yaml

    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8")) or {}
    pdb_cfg = cfg.get("palmdb", {})
    db_cfg = cfg.get("database", {})
    filt = cfg.get("filtering", {})

    sotu_to_queries = parse_palmdb_hits(
        args.palmdb_hits,
        min_pident=float(pdb_cfg.get("min_pident", 40.0)),
        max_evalue=float(pdb_cfg.get("max_evalue", 1e-5)),
        max_sotus_per_query=int(pdb_cfg.get("max_sotus_per_query", 10)),
    )
    LOG.info("Parsed %d unique PalmDB sOTU(s) from DIAMOND hits", len(sotu_to_queries))

    existing_runs, existing_rows = load_existing_manifest(args.existing_manifest)
    seed_runs = load_seed_runs(args.validated_tsv)
    LOG.info("Existing manifest has %d run(s); %d seed run(s) from validated RdRps", len(existing_runs), len(seed_runs))

    candidates = query_serratus_rds(
        db_cfg=db_cfg,
        sotu_to_queries=sotu_to_queries,
        seed_runs=seed_runs,
        include_seed_run_sotus=bool(filt.get("include_seed_run_sotus", True)),
        host_regex=str(filt.get("host_regex", "") or ""),
        min_coverage=int(filt.get("min_coverage", 10)),
        min_run_pident=float(filt.get("min_run_pident", 50.0)),
        paired_only=bool(filt.get("paired_only", True)),
        illumina_only=bool(filt.get("illumina_only", True)),
        exclude_cell_culture=bool(filt.get("exclude_cell_culture", True)),
        max_size_mb=float(filt.get("max_size_mb", 8000)),
    )

    selected = select_diverse_runs(
        candidates,
        existing_runs=existing_runs,
        max_runs_per_bioproject=int(filt.get("max_runs_per_bioproject", 4)),
        max_total_runs=int(filt.get("max_total_runs", 12)),
    )

    cand_cols = [
        "run_id", "bioproject", "scientific_name", "library_layout", "platform",
        "library_strategy", "size_mb", "n_matched_sotus", "max_coverage",
        "max_pident", "matched_sotus", "matched_queries",
    ]
    tsv_lines = ["\t".join(cand_cols)]
    for r in selected:
        tsv_lines.append("\t".join(str(r.get(c, "")) for c in cand_cols))
    atomic_write_text(args.out_candidates, "\n".join(tsv_lines) + "\n")

    snowball_rows = build_phase2_manifest_rows(
        selected,
        scenario_label=str(filt.get("scenario_label", "snowball_rdrp")),
        default_family=str(filt.get("default_family_label", "Snowball_RdRp")),
    )
    atomic_write_csv(args.out_manifest, MANIFEST_COLS, snowball_rows)

    seen: set[str] = set()
    combined_rows: list[dict[str, str]] = []
    for row in existing_rows + snowball_rows:
        rid = row.get("run_id", "")
        if rid and rid in seen:
            continue
        if rid:
            seen.add(rid)
        combined_rows.append(row)
    atomic_write_csv(args.out_combined, MANIFEST_COLS, combined_rows)

    LOG.info(
        "Selected %d new snowball SRA run(s) across %d BioProject(s) -> %s",
        len(selected),
        len({r["bioproject"] for r in selected if r["bioproject"]}),
        args.out_manifest,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
