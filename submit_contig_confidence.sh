#!/bin/bash
#SBATCH --job-name=conifer_confidence
#SBATCH --output=logs/contig_confidence_%j.log
#SBATCH --error=logs/contig_confidence_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#
# 1: score contig confidence from MetaCerberus annotation + taxonomy + RdRp
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs
module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$PWD/.cache/pip}"
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true
snakemake -s Snakefile_phase2 --executor slurm --slurm-no-account --jobs 4 \
    --cores 8 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
    --conda-frontend conda --rerun-incomplete --rerun-triggers mtime \
    --default-resources slurm_partition=shared runtime=180 mem_mb=0 \
    -- results/reporting/contig_confidence.tsv
rc=$?
echo "exit $rc"
[ -f results/reporting/contig_confidence_summary.tsv ] && cat results/reporting/contig_confidence_summary.tsv
exit $rc
