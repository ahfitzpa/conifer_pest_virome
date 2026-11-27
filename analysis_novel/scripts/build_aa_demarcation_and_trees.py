#!/usr/bin/env python3
"""Prepare amino-acid (aa) ORF/polyprotein sequences for focal Coleoptera RNA virus
families to eliminate nucleotide substitutional saturation at deep divergence.

Reads results/phylogeny/{fam}/sequences.fasta (or aln.trimmed.fasta) + selected.tsv,
translates the longest ORF per sequence in all 6 frames, writes:
  - results/aa_phylo/{fam}/proteins.faa
  - results/aa_phylo/submit_aa_phylo.sh  (MAFFT -> trimAl -> IQ-TREE2 -m MFP on protein alignments
    + aa pairwise identity via run_conifer_demarcation.py)
"""
from __future__ import annotations

import argparse
from pathlib import Path

FOCAL_FAMILIES = [
    "Iflaviridae",
    "Polycipiviridae",
    "Amarillovirales",
    "Dicistroviridae",
    "Aliusviridae",
    "Phenuiviridae",
]

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


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    cur: str | None = None
    buf: list[str] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    recs.append((cur, "".join(buf)))
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if cur is not None:
        recs.append((cur, "".join(buf)))
    return recs


def longest_orf_aa(nt: str) -> str:
    s = nt.upper().replace("-", "").replace("U", "T")
    best = ""
    for strand_seq in (s, s.translate(COMP)[::-1]):
        n = len(strand_seq)
        for frame in range(3):
            aa = "".join(CODON_TABLE.get(strand_seq[i:i + 3], "X") for i in range(frame, n - 2, 3))
            for chunk in aa.split("*"):
                if len(chunk) > len(best):
                    best = chunk
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phylogeny-dir", default=Path("results/phylogeny"), type=Path)
    ap.add_argument("--out-dir", default=Path("results/aa_phylo"), type=Path)
    ap.add_argument("--families", nargs="+", default=FOCAL_FAMILIES)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ready: list[str] = []

    for fam in args.families:
        fdir = args.phylogeny_dir / fam
        in_fa = next(
            (p for p in (fdir / "sequences.fasta", fdir / "combined.fasta", fdir / "aln.fasta", fdir / "aln.trimmed.fasta") if p.exists()),
            None,
        )
        meta = fdir / "selected.tsv"
        if not in_fa or not meta.exists():
            print(f"[skip] {fam}: missing input fasta or selected.tsv under {fdir}")
            continue
        recs = read_fasta(in_fa)
        fam_out = args.out_dir / fam
        fam_out.mkdir(parents=True, exist_ok=True)
        faa_path = fam_out / "proteins.faa"
        kept = 0
        with faa_path.open("w", encoding="utf-8") as fh:
            for sid, seq in recs:
                aa = longest_orf_aa(seq)
                if len(aa) >= 80:
                    fh.write(f">{sid}\n")
                    for i in range(0, len(aa), 80):
                        fh.write(aa[i:i + 80] + "\n")
                    kept += 1
        # Copy selected.tsv so downstream plotting & demarcation work directly
        (fam_out / "selected.tsv").write_text(meta.read_text(encoding="utf-8"), encoding="utf-8")
        print(f"[aa_phylo] {fam}: translated {kept}/{len(recs)} sequences -> {faa_path}")
        if kept >= 4:
            ready.append(fam)

    slurm = args.out_dir / "submit_aa_phylo.sh"
    fams_str = " ".join(ready)
    slurm.write_text(
        "#!/bin/bash\n"
        "#SBATCH --job-name=coleop_aa_phylo\n"
        "#SBATCH --partition=shared\n"
        "#SBATCH --cpus-per-task=8\n"
        "#SBATCH --time=06:00:00\n"
        "#SBATCH --output=results/aa_phylo/aa_phylo_%j.log\n"
        "#SBATCH --error=results/aa_phylo/aa_phylo_%j.log\n\n"
        "set -euo pipefail\n"
        "# Activate your phylogeny conda env containing mafft, trimal, iqtree2 (edit PHYLO_ENV if needed):\n"
        "PLOT_PY=/scratch/fitzpatria/calicityper_runs/manuscript/figenv/bin/python\n\n"
        f"for FAM in {fams_str}; do\n"
        "  D=results/aa_phylo/${FAM}\n"
        "  echo \"[$(date)] Aligning ${FAM} (amino acid)...\"\n"
        "  mafft --auto --thread 8 ${D}/proteins.faa > ${D}/aln.faa\n"
        "  trimal -in ${D}/aln.faa -out ${D}/aln.trimmed.faa -automated1\n"
        "  echo \"[$(date)] Inferring IQ-TREE2 protein ML phylogeny for ${FAM}...\"\n"
        "  iqtree2 -s ${D}/aln.trimmed.faa -m MFP -B 1000 --alrt 1000 -T 8 --prefix ${D}/${FAM}_aa -redo\n"
        "  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 $PLOT_PY scripts/plot_phylo_tree.py \\\n"
        "    --tree ${D}/${FAM}_aa.treefile --meta ${D}/selected.tsv \\\n"
        "    --out-svg figures/manuscript/${FAM}_aa_tree.svg \\\n"
        "    --out-png figures/manuscript/${FAM}_aa_tree.png --fetch-ncbi-hosts\n"
        "done\n\n"
        "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 $PLOT_PY scripts/run_conifer_demarcation.py \\\n"
        "  --phylogeny-dir results/aa_phylo --out-dir results/aa_demarcation\n"
        "echo \"[$(date)] All amino-acid phylogenies & aa pairwise identities complete!\"\n",
        encoding="utf-8",
    )
    print(f"Wrote {slurm} for families: {fams_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
