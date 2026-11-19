#!/bin/bash
#SBATCH --job-name=conifer_mc
#SBATCH --output=logs/rerun_metacerberus_%j.log
#SBATCH --error=logs/rerun_metacerberus_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
# Scheduler only: run_metacerberus itself is submitted as its own SLURM job.
# Targeted rerun to confirm the BROWSER_PATH/kaleido fix (chrome_binary in
# configs/paths.yaml) lets MetaCerberus get past STEP 10.
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=mamba
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true
SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 4 --cores 16 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going --default-resources slurm_partition=shared runtime=180 mem_mb=0"
echo "[INFO] Rerunning MetaCerberus only (kaleido/Chrome fix verification)"
PHASE2_ALLOW_PARTIAL_MASTER=true snakemake -s Snakefile_phase2 $SNAKE_OPTS -- results/annotation/metacerberus
echo "[INFO] exit=$?"
