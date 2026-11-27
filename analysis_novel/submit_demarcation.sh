#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=conifer_demarc
#SBATCH --output=logs/demarc_%j.log
#SBATCH --error=logs/demarc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#
# 1. Scheduler only: demarcation_all (rules 09 phylo_align -> 12) via --executor slurm.
# 2. Submit with --dependency=afterany:<novel job> so it never shares the workdir lock.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs tmp

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$PWD/.cache/pip}"
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0 slurm_extra='--exclude=sonic59'"

# 3. Target after `--`, as in submit_novel.sh.
echo "[demarc] $(date): demarcation_all"
snakemake -s Snakefile_novel $SNAKE_OPTS -- demarcation_all
rc=$?
echo "[demarc] $(date): demarcation_all exit $rc"

# 4. Report what landed.
echo "[demarc] tables: $(ls results/demarcation/*/demarcation.tsv 2>/dev/null | wc -l)"
exit $rc
