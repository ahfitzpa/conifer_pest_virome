#!/bin/bash
#SBATCH --job-name=novel_0912
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=30
#SBATCH --output=logs/dryrun_0912_%j.log
# 1: now phase2's ICTV/reporting outputs exist, check rules 09-12 resolve here
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
module purge; module load snakemake/8.30.0
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
for t in phylo_all genome_map_all recomb_all demarcation_all novel_all; do
    echo "### $t"
    out=$(snakemake -s Snakefile_novel -n --quiet rules -- "$t" 2>&1)
    rc=$?
    if [ $rc -eq 0 ]; then
        echo "  OK, total jobs: $(echo "$out" | awk '/^total/{print $2}')"
        echo "$out" | awk '/^job +count/,/^total/' | sed 's/^/    /' | head -16
    else
        echo "  FAILED"; echo "$out" | tail -6 | sed 's/^/    /'
    fi
done
