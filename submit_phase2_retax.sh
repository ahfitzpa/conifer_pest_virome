#!/bin/bash
#SBATCH --job-name=phase2_retax
#SBATCH --output=logs/phase2_retax_%j.log
#SBATCH --error=logs/phase2_retax_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --partition=shared
#
# Re-derive the taxonomy consensus after a change to the family reconciliation
# (scripts/family_taxonomy.py, scripts/merge_taxonomy_calls.py), then rebuild
# the master database and the reporting tables that every analysis reads.
#
# Forced, because merge_taxonomy_calls.py runs from a `shell:` block: Snakemake
# tracks the command string, not the script body, so editing the script is not
# a rerun trigger and the chain reports "nothing to be done".
#
# No tool is re-run: collect/annotate only re-parse the Kraken2, geNomad and
# RdRpCATCH outputs that are already on disk.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs tmp

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip
# Required, not decorative: without conda.sh the per-rule `conda activate` does
# not update PATH, so a rule's `python` is the base Anaconda interpreter. That
# has pandas but not Biopython, and all 21 collect_viral_positives jobs died on
# `from Bio import SeqIO` even though the rule's env provides it (job 1180235).
# No `conda deactivate` here: on a node with a degraded NFS mount it blocked for
# minutes before Snakemake even started.
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0"

# 1: force the consensus and everything downstream of it
snakemake -s Snakefile_phase2 $SNAKE_OPTS \
    --forcerun collect_viral_positives annotate_viral_positives \
    -- results/reporting/contig_confidence.tsv
rc=$?
echo "[retax] contig_confidence exit $rc"

# 2: the analysis family list, regenerated from the new confidence table
if [ "$rc" -eq 0 ] && [ -s results/reporting/analysis_families.yaml ]; then
    cp results/reporting/analysis_families.yaml analysis_novel/configs/families.yaml
    echo "[retax] families.yaml refreshed: $(grep -c '^    - ' analysis_novel/configs/families.yaml) families"
fi

# 3: what the reconciliation did, for the record
if [ -s results/taxonomy/master_metadata.csv ]; then
    python3 -I - <<'PY'
import csv, collections
c = collections.Counter()
d = collections.Counter()
with open("results/taxonomy/master_metadata.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        c[r.get("family_agreement", "")] += 1
        note = (r.get("family_disagreement") or "").strip()
        if note:
            d[(r.get("family_agreement", ""), r.get("consensus_family", ""), note)] += 1
print("[retax] family_agreement:", dict(c))
print("[retax] recorded disagreements:", sum(d.values()))
for (st, fam, note), n in d.most_common(15):
    print(f"   {n:4} {st:14} -> {fam:20} [{note}]")
PY
fi
exit $rc
