#!/usr/bin/env python3
"""Week 1 EVE (Retroviridae) and Mycobiome (Partitiviridae) screening helper.

Mode 1 (Extract): Reads phase2 master_metadata.csv and locates assembled contig
FASTAs to extract:
  - results/manuscript_checks/retroviridae_candidates.fasta  (for host-genome EVE check)
  - results/manuscript_checks/partitiviridae_candidates.fasta (for fungal/mycobiome BLASTx/DIAMOND check)
  - results/manuscript_checks/submit_eve_mycobiome_screens.sh (ready-to-submit SLURM script, no --mem=)

Mode 2 (Parse): Given tabular BLAST/DIAMOND outputs (outfmt 6 with stitle),
writes clean manuscript-ready summary tables:
  - results/manuscript_checks/retroviridae_eve_screen.tsv
  - results/manuscript_checks/partitiviridae_mycobiome_screen.tsv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

FUNGAL_KEYWORDS = (
    "fung", "aspergill", "penicill", "fusarium", "botrytis", "sclerotinia",
    "rhizoctonia", "rosellinia", "ophiostoma", "grosmannia", "leptographium",
    "trichoderma", "alternaria", "candida", "saccharomyces", "mycovirus",
    "gammapartitivirus", "alphapartitivirus", "betapartitivirus", "deltapartitivirus",
)
INSECT_KEYWORDS = (
    "coleoptera", "ips ", "hylobius", "locusta", "orthoptera", "insect",
    "beetle", "weevil", "arthropod", "drosophila", "tribolium", "dendroctonus",
)
EVE_KEYWORDS = (
    "chromosome", "scaffold", "genome assembly", "whole genome shotgun",
    "retrotransposon", "gypsy", "copia", "endogenous", "bel/pao", "transposable",
)


def index_fastas(search_dirs: list[Path], target_ids: set[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    if not target_ids:
        return found
    for d in search_dirs:
        if not d.exists():
            continue
        files = [d] if d.is_file() else list(d.rglob("*.fa*"))
        for fp in files:
            if fp.suffix.lower() not in {".fa", ".fna", ".fasta"}:
                continue
            cur_id: str | None = None
            keep = False
            buf: list[str] = []
            with fp.open(encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur_id and keep:
                            found[cur_id] = "".join(buf)
                        hdr = line[1:].strip().split()[0]
                        # Match exact ID or suffix after run_id|
                        cur_id = hdr
                        keep = (hdr in target_ids) or any(hdr.endswith(t) or t.endswith(hdr) for t in target_ids)
                        buf = []
                    elif keep:
                        buf.append(line.strip())
            if cur_id and keep:
                found[cur_id] = "".join(buf)
            if len(found) >= len(target_ids):
                return found
    return found


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for hdr, seq in records:
            fh.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def parse_blast_top_hits(tsv_path: Path) -> dict[str, dict[str, str]]:
    """Parse outfmt '6 qseqid sseqid pident length evalue bitscore stitle'."""
    best: dict[str, dict[str, str]] = {}
    with tsv_path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            q = parts[0]
            if q not in best:
                best[q] = {
                    "qseqid": q,
                    "sseqid": parts[1],
                    "pident": parts[2],
                    "length": parts[3],
                    "evalue": parts[4],
                    "bitscore": parts[5],
                    "stitle": parts[6] if len(parts) > 6 else "",
                }
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master-metadata", required=True, type=Path)
    ap.add_argument("--contigs-dir", nargs="*", type=Path, default=[Path("../results")])
    ap.add_argument("--out-dir", default=Path("results/manuscript_checks"), type=Path)
    ap.add_argument("--parse-eve-blast", type=Path, default=None, help="Tabular BLASTn/DIAMOND output for Retroviridae")
    ap.add_argument("--parse-myco-blast", type=Path, default=None, help="Tabular DIAMOND/BLASTx output for Partitiviridae")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    retro_rows: list[dict[str, str]] = []
    part_rows: list[dict[str, str]] = []

    with args.master_metadata.open(newline="", encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh):
            fam = (r.get("consensus_family") or "").strip()
            if fam == "Retroviridae":
                retro_rows.append(r)
            elif fam == "Partitiviridae":
                part_rows.append(r)

    print(f"[screen] Found {len(retro_rows)} Retroviridae and {len(part_rows)} Partitiviridae rows in metadata")

    # Mode 1: Extract FASTAs & emit SLURM script
    if not args.parse_eve_blast and not args.parse_myco_blast:
        want_ids = {(r.get("contig_id") or "").strip() for r in retro_rows + part_rows if r.get("contig_id")}
        seqs = index_fastas(args.contigs_dir, want_ids)
        print(f"[screen] Extracted {len(seqs)}/{len(want_ids)} target contig sequences")

        def collect(rows: list[dict[str, str]]) -> list[tuple[str, str]]:
            out: list[tuple[str, str]] = []
            seen: set[str] = set()
            for r in rows:
                cid = (r.get("contig_id") or "").strip()
                run = (r.get("run_id") or "").strip()
                host = (r.get("scientificname") or "").strip().replace(" ", "_")
                if not cid or cid in seen:
                    continue
                seq = seqs.get(cid) or next((v for k, v in seqs.items() if k.endswith(cid) or cid.endswith(k)), "")
                if seq:
                    seen.add(cid)
                    out.append((f"{run}|{cid}|{host}", seq))
            return out

        retro_fa = args.out_dir / "retroviridae_candidates.fasta"
        part_fa = args.out_dir / "partitiviridae_candidates.fasta"
        write_fasta(retro_fa, collect(retro_rows))
        write_fasta(part_fa, collect(part_rows))

        slurm = args.out_dir / "submit_eve_mycobiome_screens.sh"
        slurm.write_text(
            "#!/bin/bash\n"
            "#SBATCH --job-name=eve_myco_screen\n"
            "#SBATCH --partition=shared\n"
            "#SBATCH --cpus-per-task=8\n"
            "#SBATCH --time=04:00:00\n"
            "#SBATCH --output=results/manuscript_checks/eve_myco_%j.log\n"
            "# Note: shared partition rejects --mem=; omit explicit --mem= flag.\n\n"
            "set -euo pipefail\n"
            "# Adjust DB paths for your Sonic HPC blast/diamond databases:\n"
            "NT_DB=${NT_DB:-/scratch/databases/nt/nt}\n"
            "NR_DMND=${NR_DMND:-/scratch/databases/diamond/nr.dmnd}\n\n"
            "if [ -s results/manuscript_checks/retroviridae_candidates.fasta ]; then\n"
            "  blastn -query results/manuscript_checks/retroviridae_candidates.fasta \\\n"
            "    -db \"$NT_DB\" -num_threads 8 -max_target_seqs 5 -evalue 1e-5 \\\n"
            "    -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \\\n"
            "    -out results/manuscript_checks/retroviridae_vs_nt.tsv\n"
            "fi\n\n"
            "if [ -s results/manuscript_checks/partitiviridae_candidates.fasta ]; then\n"
            "  diamond blastx -q results/manuscript_checks/partitiviridae_candidates.fasta \\\n"
            "    -d \"$NR_DMND\" -p 8 --max-target-seqs 5 -e 1e-5 --sensitive \\\n"
            "    -f 6 qseqid sseqid pident length evalue bitscore stitle \\\n"
            "    -o results/manuscript_checks/partitiviridae_vs_nr.tsv\n"
            "fi\n",
            encoding="utf-8",
        )
        print(f"Wrote {retro_fa}, {part_fa}, and {slurm}")

    # Mode 2: Parse BLAST/DIAMOND outputs into clean manuscript tables
    if args.parse_eve_blast:
        hits = parse_blast_top_hits(args.parse_eve_blast)
        out_tsv = args.out_dir / "retroviridae_eve_screen.tsv"
        with out_tsv.open("w", newline="", encoding="utf-8") as fh:
            fields = ["qseqid", "top_hit_id", "pident", "aln_len", "evalue", "top_hit_title", "eve_call"]
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for q, h in sorted(hits.items()):
                title_l = h["stitle"].lower()
                is_host = any(k in title_l for k in EVE_KEYWORDS) or any(k in title_l for k in INSECT_KEYWORDS)
                call = "likely_EVE_or_retroelement" if is_host else "putative_exogenous_or_unconfirmed"
                w.writerow({
                    "qseqid": q, "top_hit_id": h["sseqid"], "pident": h["pident"],
                    "aln_len": h["length"], "evalue": h["evalue"],
                    "top_hit_title": h["stitle"], "eve_call": call,
                })
        print(f"Wrote {out_tsv} ({len(hits)} Retroviridae hits classified)")

    if args.parse_myco_blast:
        hits = parse_blast_top_hits(args.parse_myco_blast)
        out_tsv = args.out_dir / "partitiviridae_mycobiome_screen.tsv"
        with out_tsv.open("w", newline="", encoding="utf-8") as fh:
            fields = ["qseqid", "top_hit_id", "pident", "aln_len", "evalue", "top_hit_title", "host_origin_call"]
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for q, h in sorted(hits.items()):
                title_l = h["stitle"].lower()
                if any(k in title_l for k in FUNGAL_KEYWORDS):
                    call = "likely_fungal_mycobiome_associated"
                elif any(k in title_l for k in INSECT_KEYWORDS):
                    call = "invertebrate_associated"
                else:
                    call = "uncertain_host_origin"
                w.writerow({
                    "qseqid": q, "top_hit_id": h["sseqid"], "pident": h["pident"],
                    "aln_len": h["length"], "evalue": h["evalue"],
                    "top_hit_title": h["stitle"], "host_origin_call": call,
                })
        print(f"Wrote {out_tsv} ({len(hits)} Partitiviridae hits classified)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
