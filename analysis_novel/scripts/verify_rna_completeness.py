#!/usr/bin/env python3
"""Verify RNA virus coding-completeness (ORF architecture, intact start/stop
codons, expected ICTV genome length) for focal Coleoptera candidates and emit
a compact FASTA + SLURM remote BLASTx-vs-nr script for GenBank novelty checks.

Replaces DNA/phage-benchmarked CheckV percentages as the primary completeness
evidence for RNA viruses.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

FOCAL_FAMILIES = {
    "Iflaviridae": {"expected_kb": (9.0, 11.5), "expected_orfs": "1 large polyprotein (>=2,500 aa)"},
    "Dicistroviridae": {"expected_kb": (8.5, 10.5), "expected_orfs": "2 non-overlapping ORFs (ORF1 non-structural + ORF2 structural)"},
    "Polycipiviridae": {"expected_kb": (10.0, 12.5), "expected_orfs": "5-6 ORFs (4-5 5'-proximal structural + large 3' polyprotein)"},
    "Aliusviridae": {"expected_kb": (10.0, 14.0), "expected_orfs": "3-5 ORFs including large L (RdRp) ORF"},
    "Phenuiviridae": {"expected_kb": (6.0, 7.8), "expected_orfs": "L segment (~6.5 kb RdRp ORF; segmented genome)"},
    "Amarillovirales": {"expected_kb": (3.0, 12.0), "expected_orfs": "Flavi-like polyprotein or segmented Jingmenvirus-like segment"},
    "Phasmaviridae": {"expected_kb": (6.0, 7.5), "expected_orfs": "L segment (~6.5 kb RdRp ORF; segmented genome)"},
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def find_orfs(seq: str, min_aa: int = 100) -> list[dict[str, object]]:
    """Find all ORFs (>= min_aa) across all 6 frames, noting intact start/stop."""
    seq = seq.upper().replace("U", "T")
    orfs: list[dict[str, object]] = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        n = len(s)
        for frame in range(3):
            aa_list = [CODON_TABLE.get(s[i:i + 3], "X") for i in range(frame, n - 2, 3)]
            aa_str = "".join(aa_list)
            # Split by stop codons
            pos = 0
            for chunk in aa_str.split("*"):
                chunk_len = len(chunk)
                has_stop = (pos + chunk_len) < len(aa_str)
                m_idx = chunk.find("M")
                if m_idx != -1 and (chunk_len - m_idx) >= min_aa:
                    orf_aa = chunk[m_idx:]
                    orfs.append({
                        "strand": strand,
                        "frame": frame + 1,
                        "aa_len": len(orf_aa),
                        "has_start": True,
                        "has_stop": has_stop,
                        "nt_start": frame + (pos + m_idx) * 3 + 1,
                        "nt_end": frame + (pos + chunk_len) * 3 + (3 if has_stop else 0),
                        "aa_seq": orf_aa,
                    })
                elif chunk_len >= min_aa:
                    # 5'-truncated ORF (no upstream start before contig edge)
                    orfs.append({
                        "strand": strand,
                        "frame": frame + 1,
                        "aa_len": chunk_len,
                        "has_start": False,
                        "has_stop": has_stop,
                        "nt_start": frame + pos * 3 + 1,
                        "nt_end": frame + (pos + chunk_len) * 3 + (3 if has_stop else 0),
                        "aa_seq": chunk,
                    })
                pos += chunk_len + 1
    orfs.sort(key=lambda x: -int(x["aa_len"]))
    return orfs


def index_fastas(search_dirs: list[Path], target_ids: set[str]) -> dict[str, str]:
    found: dict[str, str] = {}
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


def classify_completeness(fam: str, length_kb: float, orfs: list[dict[str, object]]) -> str:
    lo, hi = FOCAL_FAMILIES.get(fam, {}).get("expected_kb", (3.0, 15.0))
    intact = [o for o in orfs if o["has_start"] and o["has_stop"]]
    longest = orfs[0] if orfs else None
    if not longest:
        return "fragment"
    if fam in {"Phenuiviridae", "Phasmaviridae"}:
        if length_kb >= lo and longest["has_start"] and longest["has_stop"]:
            return "coding-complete_L_segment"
        return "partial_segment"
    if fam == "Amarillovirales":
        if length_kb >= 8.5 and longest["has_start"] and longest["has_stop"]:
            return "coding-complete_genome"
        return "coding-complete_or_partial_segment" if (longest["has_start"] and longest["has_stop"]) else "partial"
    if length_kb >= lo and longest["has_start"] and longest["has_stop"]:
        if fam == "Dicistroviridae" and len([o for o in intact if int(o["aa_len"]) >= 600]) >= 2:
            return "coding-complete_bicistronic"
        if fam == "Polycipiviridae" and len([o for o in intact if int(o["aa_len"]) >= 150]) >= 4:
            return "coding-complete_polycistronic"
        return "coding-complete"
    if length_kb >= lo * 0.85 and (longest["has_start"] or longest["has_stop"]):
        return "near-coding-complete"
    return "partial"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table1", default=Path("results/manuscript_checks/table1_top_candidates.tsv"), type=Path)
    ap.add_argument("--contigs-dir", nargs="*", type=Path, default=[Path("../results")])
    ap.add_argument("--out-dir", default=Path("results/manuscript_checks"), type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    with args.table1.open(encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fam = r.get("family", "")
            host = (r.get("host") or "").lower()
            if fam in FOCAL_FAMILIES and ("ips" in host or "hylobius" in host):
                rows.append(r)

    want_ids = {r["contig_id"] for r in rows if r.get("contig_id")}
    seqs = index_fastas(args.contigs_dir, want_ids)

    out_rows: list[dict[str, object]] = []
    fa_records: list[tuple[str, str]] = []

    for r in rows:
        cid = r["contig_id"]
        fam = r["family"]
        seq = seqs.get(cid) or next((v for k, v in seqs.items() if k.endswith(cid) or cid.endswith(k)), "")
        if not seq:
            continue
        length_kb = round(len(seq) / 1000.0, 2)
        orfs = find_orfs(seq, min_aa=100)
        intact_orfs = [o for o in orfs if o["has_start"] and o["has_stop"]]
        top3 = "; ".join(f"{o['aa_len']}aa({'ATG' if o['has_start'] else '5trunc'}-{'STOP' if o['has_stop'] else '3trunc'},{o['strand']})" for o in orfs[:3])
        status = classify_completeness(fam, length_kb, orfs)
        lo, hi = FOCAL_FAMILIES[fam]["expected_kb"]
        out_rows.append({
            "family": fam,
            "host": r["host"],
            "bioproject": r["bioproject"],
            "run_id": r["run_id"],
            "contig_id": cid,
            "length_kb": length_kb,
            "expected_ictv_kb": f"{lo}-{hi}",
            "spades_cov": r.get("spades_cov", ""),
            "n_orfs_ge_100aa": len(orfs),
            "n_intact_orfs": len(intact_orfs),
            "longest_orf_aa": orfs[0]["aa_len"] if orfs else 0,
            "top3_orf_architecture": top3,
            "rna_completeness_call": status,
        })
        hdr = f"{fam}|{r['run_id']}|{cid}|{r['host'].replace(' ', '_')}"
        fa_records.append((hdr, seq))

    tsv_path = args.out_dir / "coleoptera_rna_orf_verification.tsv"
    if out_rows:
        with tsv_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(out_rows)
        print(f"Wrote {tsv_path} ({len(out_rows)} Coleoptera focal contigs verified)")

    fa_path = args.out_dir / "coleoptera_top_candidates.fasta"
    with fa_path.open("w", encoding="utf-8") as fh:
        for hdr, seq in fa_records:
            fh.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    print(f"Wrote {fa_path} ({len(fa_records)} sequences for GenBank nr BLASTx check)")

    slurm = args.out_dir / "submit_nr_novelty_check.sh"
    slurm.write_text(
        "#!/bin/bash\n"
        "#SBATCH --job-name=coleop_nr_blastx\n"
        "#SBATCH --partition=shared\n"
        "#SBATCH --cpus-per-task=1\n"
        "#SBATCH --time=04:00:00\n"
        "#SBATCH --output=results/manuscript_checks/coleop_nr_%j.log\n"
        "#SBATCH --error=results/manuscript_checks/coleop_nr_%j.log\n\n"
        "set -euo pipefail\n"
        "module load ncbi-blast/2.16.0\n"
        "echo \"[$(date)] Running remote blastx vs GenBank nr on top Coleoptera RNA virus candidates...\"\n"
        "blastx -remote -db nr \\\n"
        "  -query results/manuscript_checks/coleoptera_top_candidates.fasta \\\n"
        "  -max_target_seqs 5 -evalue 1e-5 \\\n"
        "  -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore stitle' \\\n"
        "  -out results/manuscript_checks/coleoptera_vs_nr_blastx.tsv\n"
        "echo \"[$(date)] Done! ($(wc -l < results/manuscript_checks/coleoptera_vs_nr_blastx.tsv) hits)\"\n",
        encoding="utf-8",
    )
    print(f"Wrote {slurm}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
