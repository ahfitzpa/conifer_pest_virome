#!/bin/bash
#SBATCH --job-name=novel_dryrun
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=25
#SBATCH --output=/scratch/fitzpatria/Serratus/workflow/analysis_novel/logs/dryrun_%j.log
# 1: parse analysis_novel and dry-run the phylogeny on the 3-run test manifest
set -uo pipefail
cd /scratch/fitzpatria/Serratus/workflow/analysis_novel
module purge; module load snakemake/8.30.0
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export NOVEL_MANIFEST=../common/test_data/manifest_novel_test.csv
echo "### rule list"
snakemake -s Snakefile_novel --list; echo "rc=$?"
echo "### dry run: phylogeny_all on the 3-run test manifest"
snakemake -s Snakefile_novel -n --quiet rules -- phylogeny_all; echo "rc=$?"
echo "### dry run: novel_all"
snakemake -s Snakefile_novel -n --quiet rules -- novel_all; echo "rc=$?"
