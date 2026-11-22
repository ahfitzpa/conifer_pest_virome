#!/bin/bash
#SBATCH --job-name=novel_full_test
#SBATCH --output=/scratch/fitzpatria/Serratus/workflow/analysis_novel/logs/full_test_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#
# End-to-end test of analysis_novel on the 3-run manifest, through the stages
# that do not need phase2's ICTV/reporting outputs (rules 09-12 wait on the
# serratus_port job). Covers: phylogeny -> ORFs -> ORF1 -> codon alignment ->
# tree pruning -> FUBAR/MEME -> selection summary + composition signatures.
set -uo pipefail
cd /scratch/fitzpatria/Serratus/workflow/analysis_novel
mkdir -p logs tmp
module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true
export NOVEL_MANIFEST=../common/test_data/manifest_novel_test.csv

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 20 --cores 32 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0"

for stage in phylogeny_all codon_all signatures_all; do
    echo "[test] $(date): $stage"
    snakemake -s Snakefile_novel $SNAKE_OPTS -- $stage
    echo "[test] $(date): $stage exit $?"
done

echo "### outputs"
find results -type f \( -name "*.tsv" -o -name "*.fasta" -o -name "*.newick" -o -name "*.json" \) | sort
