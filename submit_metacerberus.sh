#!/bin/bash
#SBATCH --job-name=conifer_metac
#SBATCH --output=logs/metac_%j.log
#SBATCH --error=logs/metac_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#
# 1. Rerun MetaCerberus alone. Job 1180648 died at STEP 8 of 10 on 2026-11-23
#    because its parent driver (1180411) was cancelled at 11:43:33, not because
#    anything failed: PHROG/VOG searches were still empty and nothing reached
#    final/. Its final_annotation_summary.tsv is what rules 09 (concat_prepare)
#    and 10 (genome_map_prepare) read, so both die at DAG build without it.
# 2. Scheduler only: the rule itself is submitted as its own SLURM job.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs tmp

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"

# 3. STEP 10 renders figures through kaleido, which needs this Chrome.
CHROME="/scratch/fitzpatria_sh/chrome/chrome-linux64/chrome"
[ -x "$CHROME" ] || { echo "[ERROR] chrome missing at $CHROME" >&2; exit 1; }

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0"

echo "[metac] $(date): results/annotation/metacerberus"
PHASE2_ALLOW_PARTIAL_MASTER=true snakemake -s Snakefile_phase2 $SNAKE_OPTS \
    -- results/annotation/metacerberus
rc=$?
echo "[metac] $(date): exit $rc"

# 4. The one file rules 09 and 10 actually need.
SUMMARY="results/annotation/metacerberus/final/annotations-FragGeneScan_metacerberus_input/final_annotation_summary.tsv"
if [ -s "$SUMMARY" ]; then
    echo "[metac] summary OK: $(wc -l < "$SUMMARY") lines"
else
    echo "[metac] ERROR: $SUMMARY missing or empty" >&2
    rc=1
fi
exit $rc
