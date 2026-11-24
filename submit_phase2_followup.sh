#!/bin/bash
#SBATCH --job-name=conifer_p2_fu
#SBATCH --output=logs/snakemake_p2_followup_%j.log
#SBATCH --error=logs/snakemake_p2_followup_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=mamba
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$PWD/.cache/pip}"
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda --conda-frontend conda --rerun-incomplete --rerun-triggers mtime params input --keep-going --default-resources slurm_partition=shared runtime=180 mem_mb=0 --resources sra_fetch=4"

# 1: dry run, logged for review
snakemake -s Snakefile_phase2 $SNAKE_OPTS -n --quiet rules

# 2: full workflow (reruns RdRpCATCH timeouts with runtime=1440)
echo "[follow-up] $(date): stage 2"
snakemake -s Snakefile_phase2 $SNAKE_OPTS || true

# 3: master database + MetaCerberus over completed samples
echo "[follow-up] $(date): stage 3"
PHASE2_ALLOW_PARTIAL_MASTER=true snakemake -s Snakefile_phase2 $SNAKE_OPTS \
    --rerun-triggers mtime \
    -- results/master_viral_database.fasta results/taxonomy/excluded_samples.tsv \
    results/annotation/metacerberus
rc=$?

# 4: DAG + report
snakemake -s Snakefile_phase2 --dag --jobs 1 | dot -Tsvg > logs/phase2_dag.svg
snakemake -s Snakefile_phase2 --report report_phase2.html --jobs 1
exit $rc
