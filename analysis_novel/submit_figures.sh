#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=conifer_figures
#SBATCH --output=logs/figures_%j.log
#SBATCH --error=logs/figures_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#
# 1. Scheduler only: figures_all (rules/15.figures.smk), each figure its own job.
# 2. Run after demarcation and signatures: the overview figures read those
#    tables, and novelty_vs_completeness silently skipped itself on 2026-11-23
#    when results/demarcation was empty.
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

echo "[figures] $(date): figures_all"
snakemake -s Snakefile_novel $SNAKE_OPTS -- figures_all
rc=$?
echo "[figures] $(date): figures_all exit $rc"

# 3. Report what landed.
echo "[figures] overview: $(ls figures/overview/*.png 2>/dev/null | wc -l) png"
echo "[figures] family trees: $(ls figures/*/*_tree_linear.png 2>/dev/null | wc -l)"
exit $rc
