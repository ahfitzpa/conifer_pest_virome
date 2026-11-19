#!/bin/bash
#SBATCH --job-name=conifer_p3
#SBATCH --output=logs/snakemake_p3_%j.log
#SBATCH --error=logs/snakemake_p3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL

# Phase 3 scheduler job: phylogenetics, conserved domains, diagnostic residues.
# Every rule is submitted as its own SLURM job by --executor slurm.

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip

source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete \
 --rerun-triggers mtime params input --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0"

echo "[SLURM INFO] Phase 3: ML phylogenetics + conserved domains + diagnostic residues"
snakemake -s Snakefile_phase2 $SNAKE_OPTS phylo_all

echo "[SLURM INFO] Phase 3: genome maps"
snakemake -s Snakefile_phase2 $SNAKE_OPTS genome_map_all

echo "[SLURM INFO] Phase 3: recombination screening"
snakemake -s Snakefile_phase2 $SNAKE_OPTS recomb_all

echo "[SLURM INFO] Phase 3: ICTV demarcation tables"
snakemake -s Snakefile_phase2 $SNAKE_OPTS demarcation_all

snakemake -s Snakefile_phase2 --dag --jobs 1 | dot -Tsvg > logs/phase3_dag.svg
snakemake -s Snakefile_phase2 --report report_phase3.html --jobs 1
