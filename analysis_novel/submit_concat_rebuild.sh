#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=conifer_concat
#SBATCH --output=logs/concat_%j.log
#SBATCH --error=logs/concat_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#
# 1. Rebuild the concatenated-domain sets after the concat_domains.py contig-id
#    and assigned_family fix; --rerun-triggers mtime ignores code changes, so
#    the stale sets need --forcerun concat_prepare to be rewritten.
# 2. Submit with --dependency=afterany:<novel> so it never shares the workdir lock.
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
# 3. Target after `--`: --forcerun would otherwise swallow it.
echo "[concat] $(date): phylo_all with forced concat_prepare"
snakemake -s Snakefile_novel $SNAKE_OPTS --forcerun concat_prepare -- phylo_all
rc=$?
echo "[concat] $(date): phylo_all exit $rc"
# 4. Report the taxa each set ended up with, candidates separately: a
#    reference-only set is the failure this rebuild exists to correct.
printf '%-20s %6s %6s %6s\n' family taxa ICTV cand
for f in results/phylo/families/*/concat_domains.faa; do
    fam=$(basename "$(dirname "$f")")
    tot=$(grep -c '^>' "$f" 2>/dev/null || echo 0)
    ictv=$(grep '^>' "$f" 2>/dev/null | grep -c '^>ICTV|' || echo 0)
    printf '%-20s %6s %6s %6s\n' "$fam" "$tot" "$ictv" "$((tot - ictv))"
done
echo "[concat] concat trees: $(ls results/phylo/families/*/concat.tree.newick 2>/dev/null | wc -l)"
exit $rc
