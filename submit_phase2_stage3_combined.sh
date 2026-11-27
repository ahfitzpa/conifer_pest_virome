#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=conifer_s3_comb
#SBATCH --output=logs/stage3_combined_%j.log
#SBATCH --error=logs/stage3_combined_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#
# Rebuild Phase 2 Stage 3 (master_viral_database, metacerberus, contig_confidence,
# contigs_keep, ictv/families.tsv) over ALL 33 runs in manifest_combined.csv
# (21 original conifer runs + 12 snowball runs).
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$PWD/.cache/pip}"
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

export PHASE2_MANIFEST_OVERRIDE="$PWD/snowball/results/snowball/manifest_combined.csv"

echo "[stage3_combined] $(date): rebuilding Stage 3 on $PHASE2_MANIFEST_OVERRIDE"
snakemake -s Snakefile_phase2 \
    --executor slurm --slurm-no-account --jobs 12 --cores 32 \
    --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
    --conda-frontend conda --rerun-incomplete --rerun-triggers mtime \
    --forcerun build_master_database run_metacerberus export_metacerberus_orfs contig_confidence \
    --default-resources slurm_partition=shared runtime=360 mem_mb=0 \
    -- \
    results/master_viral_database.fasta \
    results/taxonomy/master_metadata.csv \
    results/annotation/metacerberus \
    results/reporting/metacerberus_orfs.faa \
    results/reporting/contig_confidence.tsv \
    results/reporting/contigs_keep.txt \
    results/reporting/analysis_families.yaml \
    results/ictv/families.tsv
rc=$?
echo "[stage3_combined] $(date): exit $rc"
if [ "$rc" -eq 0 ] && [ -s results/reporting/analysis_families.yaml ]; then
    cp results/reporting/analysis_families.yaml analysis_novel/configs/families.yaml
    echo "[stage3_combined] synced analysis_families.yaml -> analysis_novel/configs/families.yaml"
fi
exit $rc
