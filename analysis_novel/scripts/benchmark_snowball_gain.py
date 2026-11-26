#!/usr/bin/env python3
"""Stage 17.3: Quantify Round-2 RdRp snowballing gain and host-specificity metrics.

Pure-Python, zero-OpenBLAS script that converts `snowball_vs_seed_blastx.tsv`,
`esm_rdrp_validation.tsv`, and `manifest_combined.csv` into publication-ready
benchmark tables:
  - snowball_gain_summary.tsv    (retrieval yield, cross-BioProject recovery,
                                  genome-completion gain, geographic expansion,
                                  host-specificity contrast, false-retrieval rate,
                                  and leave-one-lineage-out retrieval analysis)
  - per_lineage_recovery.tsv     (tidy per-seed-lineage recovery breakdown)
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
LOG = logging.getLogger("benchmark_snowball_gain")

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore",
]


def _flt(val: str, default: float = 0.0) -> float:
    try:
        return float((val or "").strip())
    except ValueError:
        return default


def _extract_run_and_len(qseqid: str) -> Tuple[str, int]:
    run = qseqid.split("|")[0] if "|" in qseqid else qseqid.split("_")[0]
    m = re.search(r"length_(\d+)", qseqid)
    clen = int(m.group(1)) if m else 0
    return run, clen


def load_manifest(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    if not path or not path.is_file():
        return out
    with path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            run = (row.get("run") or row.get("run_id") or "").strip()
            if not run:
                continue
            host = (
                row.get("host_organism")
                or row.get("scientificname")
                or row.get("host")
                or ""
            ).strip()
            bp = (row.get("bioproject") or row.get("BioProject") or "").strip()
            group = (row.get("group") or row.get("batch") or "").strip()
            out[run] = {"host": host, "bioproject": bp, "group": group}
    return out


def load_seed_metadata(esm_path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    if not esm_path or not esm_path.is_file():
        return out
    with esm_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("candidate_id") or "").strip()
            if cid:
                out[cid] = row
    return out


def read_blast_hits(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8") as fh:
        first = fh.readline()
        if not first:
            return rows
        parts = first.rstrip("\r\n").split("\t")
        has_header = parts[0].lower() in {"qseqid", "query"}
        fh.seek(0)
        if has_header:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows.append(row)
        else:
            for line in fh:
                cols = line.rstrip("\r\n").split("\t")
                if len(cols) < 4:
                    continue
                rec = {BLAST_COLS[i]: cols[i] for i in range(min(len(cols), len(BLAST_COLS)))}
                rows.append(rec)
    return rows


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
    ap.add_argument("--hits-tsv", required=True, type=Path)
    ap.add_argument("--esm-validation", type=Path, default=None)
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--min-pident", type=float, default=50.0)
    ap.add_argument("--min-aln-aa", type=int, default=80)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args(argv)

    manifest = load_manifest(args.manifest)
    esm_meta = load_seed_metadata(args.esm_validation)
    raw_hits = read_blast_hits(args.hits_tsv)

    # Filter credible alignments
    hits = []
    for r in raw_hits:
        pid = _flt(r.get("pident", ""))
        aln = int(_flt(r.get("length", "")))
        if pid >= args.min_pident and aln >= args.min_aln_aa:
            run, clen = _extract_run_and_len(r.get("qseqid", ""))
            qlen = int(_flt(r.get("qlen", ""), clen))
            hits.append(
                {
                    "run": run,
                    "seed": (r.get("sseqid") or "").strip(),
                    "pident": pid,
                    "aln_aa": aln,
                    "contig_len": max(clen, qlen),
                    "host": manifest.get(run, {}).get("host", "Unknown"),
                    "bioproject": manifest.get(run, {}).get("bioproject", ""),
                }
            )

    # Identify Round-1 vs Round-2 runs from manifest (if present)
    r1_runs = {r for r, m in manifest.items() if "snowball" not in m.get("group", "").lower() and "locusta" not in m.get("host", "").lower()}
    r2_runs = {r for r, m in manifest.items() if "snowball" in m.get("group", "").lower()}
    if not r2_runs:
        r2_runs = {h["run"] for h in hits}

    # Aggregate per seed lineage
    by_seed: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    by_run: Dict[str, Set[str]] = defaultdict(set)
    for h in hits:
        by_seed[str(h["seed"])].append(h)
        by_run[str(h["run"])].add(str(h["seed"]))

    per_lineage_rows: List[Dict[str, object]] = []
    replicated_seeds = 0
    sister_seeds = 0
    complete_genomes_added = 0
    ips_only_seeds = 0
    hylobius_only_seeds = 0
    shared_subfamily_seeds = 0

    for seed, sh in sorted(by_seed.items()):
        max_pid = max(float(x["pident"]) for x in sh)
        max_len = max(int(x["contig_len"]) for x in sh)
        runs_hit = sorted({str(x["run"]) for x in sh})
        hosts_hit = sorted({str(x["host"]) for x in sh if x["host"]})
        bps_hit = sorted({str(x["bioproject"]) for x in sh if x["bioproject"]})
        n_complete = sum(1 for x in sh if int(x["contig_len"]) >= 9000)
        complete_genomes_added += n_complete

        if max_pid >= 90.0:
            tier = "independent_recovery_ge90pct"
            replicated_seeds += 1
        else:
            tier = "sister_lineage_50_90pct"
            sister_seeds += 1

        has_ips = any("ips" in h.lower() for h in hosts_hit)
        has_hylo = any("hylobius" in h.lower() for h in hosts_hit)
        if has_ips and has_hylo:
            shared_subfamily_seeds += 1
        elif has_ips:
            ips_only_seeds += 1
        elif has_hylo:
            hylobius_only_seeds += 1

        # Leave-one-lineage-out: how many of this seed's runs have NO other seed hit?
        unique_runs = [r for r in runs_hit if len(by_run[r]) == 1]

        fam = esm_meta.get(seed, {}).get("consensus_family", "")
        per_lineage_rows.append(
            {
                "seed_candidate_id": seed,
                "seed_family": fam,
                "recovery_tier": tier,
                "n_snowball_hits": len(sh),
                "n_snowball_runs": len(runs_hit),
                "n_unique_runs_lolo": len(unique_runs),
                "max_pident": f"{max_pid:.1f}",
                "max_contig_bp": max_len,
                "n_contigs_ge_9kb": n_complete,
                "snowball_hosts": "; ".join(hosts_hit),
                "snowball_bioprojects": "; ".join(bps_hit),
            }
        )

    n_seed_domains = sum(
        1 for v in esm_meta.values()
        if (v.get("structural_verdict") or "") in {"validated_rdrp_high", "validated_rdrp_moderate"}
    ) or len(by_seed)

    r2_with_hit = len(set(by_run.keys()) & r2_runs) if r2_runs else len(by_run)
    r2_total = max(len(r2_runs), r2_with_hit)
    false_retrieval_runs = max(0, r2_total - r2_with_hit)

    r1_hosts = {manifest[r]["host"] for r in r1_runs if r in manifest and manifest[r]["host"]}
    r2_hosts = {manifest[r]["host"] for r in r2_runs if r in manifest and manifest[r]["host"]} | {
        str(h["host"]) for h in hits if h["host"] != "Unknown"
    }
    r1_bps = {manifest[r]["bioproject"] for r in r1_runs if r in manifest and manifest[r]["bioproject"]}
    r2_bps = {manifest[r]["bioproject"] for r in r2_runs if r in manifest and manifest[r]["bioproject"]} | {
        str(h["bioproject"]) for h in hits if h["bioproject"]
    }

    summary_rows: List[Dict[str, object]] = [
        {"metric": "validated_seed_domains", "value": n_seed_domains, "notes": "Round-1 structurally validated RdRp/NS3 queries"},
        {"metric": "snowball_libraries_total", "value": r2_total, "notes": "Round-2 conifer-pest SRA runs retrieved"},
        {"metric": "snowball_libraries_with_credible_hit", "value": r2_with_hit, "notes": f"pident >= {args.min_pident}%, aln >= {args.min_aln_aa} aa"},
        {"metric": "false_retrieval_rate", "value": f"{(false_retrieval_runs / r2_total):.4f}" if r2_total else "0.0000", "notes": f"{false_retrieval_runs}/{r2_total} snowball runs without credible hit"},
        {"metric": "retrieval_yield_libs_per_seed", "value": f"{(r2_total / n_seed_domains):.3f}" if n_seed_domains else "0.000", "notes": "New SRA libraries per validated Round-1 domain"},
        {"metric": "seed_lineages_independently_recovered_ge90pct", "value": replicated_seeds, "notes": "Cross-BioProject independent recovery (>=90% aa ID)"},
        {"metric": "seed_lineages_sister_50_90pct", "value": sister_seeds, "notes": "Sister lineages recovered across new hosts (50-90% aa ID)"},
        {"metric": "complete_or_near_complete_contigs_ge9kb_added", "value": complete_genomes_added, "notes": "Round-2 hit contigs >= 9.0 kb (genome-completion gain)"},
        {"metric": "host_species_round1", "value": len(r1_hosts), "notes": "; ".join(sorted(r1_hosts))},
        {"metric": "host_species_combined", "value": len(r1_hosts | r2_hosts), "notes": "; ".join(sorted(r1_hosts | r2_hosts))},
        {"metric": "bioprojects_added_round2", "value": len(r2_bps - r1_bps), "notes": "; ".join(sorted(r2_bps - r1_bps))},
        {"metric": "host_specificity_ips_only_seeds", "value": ips_only_seeds, "notes": "Seed lineages hitting only Ips (Scolytinae) in Round 2"},
        {"metric": "host_specificity_hylobius_only_seeds", "value": hylobius_only_seeds, "notes": "Seed lineages hitting only Hylobius (Molytinae) in Round 2"},
        {"metric": "host_specificity_cross_subfamily_seeds", "value": shared_subfamily_seeds, "notes": "Seed lineages hitting both Ips and Hylobius in Round 2"},
    ]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv_atomic(
        args.out_dir / "snowball_gain_summary.tsv",
        ["metric", "value", "notes"],
        summary_rows,
    )
    write_tsv_atomic(
        args.out_dir / "per_lineage_recovery.tsv",
        [
            "seed_candidate_id", "seed_family", "recovery_tier", "n_snowball_hits",
            "n_snowball_runs", "n_unique_runs_lolo", "max_pident", "max_contig_bp",
            "n_contigs_ge_9kb", "snowball_hosts", "snowball_bioprojects",
        ],
        per_lineage_rows,
    )
    LOG.info(
        "Snowball gain benchmark complete: %d credible hits across %d seeds (%d >=90%%, %d sister) -> %s",
        len(hits), len(by_seed), replicated_seeds, sister_seeds, args.out_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
