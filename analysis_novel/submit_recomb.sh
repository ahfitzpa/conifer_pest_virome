#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=conifer_recomb
#SBATCH --output=logs/recomb_%j.log
#SBATCH --error=logs/recomb_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#
# 1. Scheduler only: recomb_all on the NUCLEOTIDE alignment, each rule its own job.
# 2. Submit with --dependency=afterany:<previous> so it never shares the workdir lock.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs tmp

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0 slurm_extra='--exclude=sonic59'"

# 3. recomb_prepare reruns on its own: its outputs are renamed (.fasta, candidates.txt).
echo "[recomb] $(date): recomb_all"
snakemake -s Snakefile_novel $SNAKE_OPTS -- recomb_all
rc=$?
echo "[recomb] $(date): recomb_all exit $rc"

# 4. Report eligibility and any breakpoints found.
echo "[recomb] eligible families:"
awk -F'\t' 'FNR==2 && $11==1 {print "  " $1 "  refs=" $2 " cand=" $3 " cols=" $4}' \
    results/recombination/*/eligible.tsv
echo "[recomb] events: $(cat results/recombination/*/recombination_summary.tsv 2>/dev/null | grep -vc '^family') rows"
exit $rc
