#!/usr/bin/env bash
#SBATCH --job-name=open_host_spotcheck
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=logs/spotcheck_%j.out
#SBATCH --error=logs/spotcheck_%j.err
# =============================================================================
# Option-B Open-Host Spot-Check Launcher (for mSystems host-specificity check)
# =============================================================================
# 1. Runs `open_host_specificity` (open-host Serratus SQL + 4 compact spot-check
#    runs in results/open_host/manifest_spotcheck.csv)
# 2. Assembles those 4 small runs via Phase 2 (`Snakefile_phase2`)
# 3. Cross-matches assembled spot-check contigs against Round-1 validated RdRps
#    (`spotcheck_cross_match` -> results/open_host/spotcheck_vs_seed_blastx.tsv)
#
# Submit from snowball/:
#   sbatch submit_spotcheck.sh
# =============================================================================
set -euo pipefail

module load snakemake/8.30.0 2>/dev/null || true
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

SNOWBALL_DIR="$(pwd)"
REPO_ROOT="$(cd "${SNOWBALL_DIR}/.." && pwd)"
mkdir -p "${SNOWBALL_DIR}/logs"

echo "[$(date -Is)] Step 1: Open-host specificity query & spot-check manifest..."
snakemake -s "${SNOWBALL_DIR}/Snakefile_snowball" \
    --cores 4 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix "${SNOWBALL_DIR}/.snakemake/conda" \
    --rerun-triggers mtime \
    open_host_specificity

SPOTCHECK_MANIFEST="${SNOWBALL_DIR}/results/open_host/manifest_spotcheck.csv"
if [ ! -s "${SPOTCHECK_MANIFEST}" ] || [ "$(wc -l < "${SPOTCHECK_MANIFEST}")" -le 1 ]; then
    echo "[$(date -Is)] No spot-check runs selected; exiting."
    exit 0
fi

echo "[$(date -Is)] Selected spot-check runs:"
cat "${SPOTCHECK_MANIFEST}"

echo "[$(date -Is)] Step 2: Running Phase 2 assembly on OpenHost_SpotCheck runs..."
(
    cd "${REPO_ROOT}"
    PHASE2_MANIFEST_OVERRIDE="${SPOTCHECK_MANIFEST}" \
    snakemake -s Snakefile_phase2 \
        --cores 4 \
        --use-conda \
        --conda-frontend conda \
        --conda-prefix "${REPO_ROOT}/.snakemake/conda" \
        --keep-going \
        --rerun-triggers mtime \
        --set-resources filter_contigs:runtime=720
)

echo "[$(date -Is)] Step 3: Cross-matching assembled spot-check contigs against seed RdRps..."
snakemake -s "${SNOWBALL_DIR}/Snakefile_snowball" \
    --cores 4 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix "${SNOWBALL_DIR}/.snakemake/conda" \
    --rerun-triggers mtime \
    spotcheck_cross_match

echo "[$(date -Is)] Done! Results in ${SNOWBALL_DIR}/results/open_host/"
