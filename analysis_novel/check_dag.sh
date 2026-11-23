#!/bin/bash
#SBATCH --job-name=conifer_dag
#SBATCH --output=logs/dag_%j.log
#SBATCH --error=logs/dag_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#
# 1. Dry-run each target to confirm the rule files still build a DAG.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
module purge
module load snakemake/8.30.0
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export OPENBLAS_NUM_THREADS=1
for t in phylo_all recomb_all demarcation_all signatures_all codon_all genome_map_all figures_all; do
    echo "### $t"
    snakemake -s Snakefile_novel --dryrun -- "$t" 2>&1 \
        | grep -E "^(Job stats|total|[a-z_]+ +[0-9]+)|Error|Exception|MissingInput" | head -15
done
