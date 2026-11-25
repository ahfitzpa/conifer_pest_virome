#!/bin/bash
#SBATCH --job-name=rdrp_snowball
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-04:00:00
#SBATCH --output=logs/snowball_%j.out
#SBATCH --error=logs/snowball_%j.err
#
# Side workflow: Validated Novel RdRps -> PalmDB/Serratus Snowball -> Phase 2
#
# Mode 1 (DEFAULT — safe preview, ~2-5 min, downloads NO SRA reads):
#   sbatch submit_snowball.sh
#   Inspect results/snowball/snowball_candidate_runs.tsv and manifest_snowball.csv
#
# Mode 2 (FULL LOOP — after reviewing manifest_snowball.csv, launch Phase 2
#         download + assembly + taxonomy + analysis_novel on the new runs):
#   RUN_PIPELINE=1 sbatch submit_snowball.sh

set -euo pipefail

REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR:-$PWD}/.." && pwd)"
SNOWBALL_DIR="${REPO_ROOT}/snowball"
cd "${SNOWBALL_DIR}"

mkdir -p logs "${REPO_ROOT}/.cache"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/.cache}"

if command -v module >/dev/null 2>&1; then
    module purge || true
    module load snakemake/8.30.0 || true
fi
export SNAKEMAKE_CONDA_FRONTEND="${CONDA_FRONTEND:-conda}"
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
    conda deactivate 2>/dev/null || true
fi

echo "[snowball] Step 1: Discovering sister-lineage SRA runs via PalmDB + Serratus RDS"
snakemake \
    -s Snakefile_snowball \
    --use-conda \
    --conda-prefix .snakemake/conda \
    --conda-frontend "${CONDA_FRONTEND:-conda}" \
    --cores "${SLURM_CPUS_PER_TASK:-4}" \
    --rerun-incomplete \
    -- snowball_discover

MANIFEST_NEW="${SNOWBALL_DIR}/results/snowball/manifest_snowball.csv"
N_NEW=$(( $(wc -l < "${MANIFEST_NEW}") - 1 ))
echo "[snowball] Discovered ${N_NEW} new candidate SRA run(s) -> ${MANIFEST_NEW}"
echo "[snowball] Review: ${SNOWBALL_DIR}/results/snowball/snowball_candidate_runs.tsv"

if [ "${RUN_PIPELINE:-0}" != "1" ]; then
    echo "[snowball] RUN_PIPELINE=0 (discovery-only). Re-run with RUN_PIPELINE=1 to download & run Phase 2 + analysis_novel."
    exit 0
fi

if [ "${N_NEW}" -le 0 ]; then
    echo "[snowball] 0 new runs to process; exiting."
    exit 0
fi

echo "[snowball] Step 2: Launching general Phase 2 workflow on ${MANIFEST_NEW}"
cd "${REPO_ROOT}"
PHASE2_SNAKEFILE="${PHASE2_SNAKEFILE:-Snakefile_phase2}"
if [ ! -f "${PHASE2_SNAKEFILE}" ] && [ -f "Snakefile" ]; then
    PHASE2_SNAKEFILE="Snakefile"
fi
export PHASE2_MANIFEST_OVERRIDE="${MANIFEST_NEW}"
snakemake \
    -s "${PHASE2_SNAKEFILE}" \
    --executor slurm \
    --slurm-no-account \
    --jobs 12 \
    --cores 32 \
    --latency-wait 300 \
    --use-conda \
    --conda-prefix .snakemake/conda \
    --conda-frontend "${CONDA_FRONTEND:-conda}" \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --default-resources slurm_partition=shared runtime=240 mem_mb=0 \
    --config "manifest=${MANIFEST_NEW}"

echo "[snowball] Step 3: Launching analysis_novel on combined manifest"
cd "${REPO_ROOT}/analysis_novel"
NOVEL_MANIFEST="${SNOWBALL_DIR}/results/snowball/manifest_combined.csv" \
    sbatch submit_novel.sh
echo "[snowball] Submitted analysis_novel follow-up job."
