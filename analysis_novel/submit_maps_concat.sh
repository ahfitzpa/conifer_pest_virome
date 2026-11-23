#!/bin/bash
#SBATCH --job-name=conifer_maps
#SBATCH --output=logs/maps_%j.log
#SBATCH --error=logs/maps_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#
# 1. phylo_all (rule 09) and genome_map_all (rule 10), both of which read
#    MetaCerberus's final_annotation_summary.tsv and so died at DAG build until
#    it existed.
# 2. Submit with --dependency=afterok:<metacerberus>,afterany:<recomb> so the
#    summary is there and no two snakemake runs share this workdir lock.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs tmp

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"

# 3. Refuse to start if the blocking input is still absent.
SUMMARY="../results/annotation/metacerberus/final/annotations-FragGeneScan_metacerberus_input/final_annotation_summary.tsv"
[ -s "$SUMMARY" ] || { echo "[maps] ERROR: MetaCerberus summary still missing: $SUMMARY" >&2; exit 1; }

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0"

status=0
for stage in phylo_all genome_map_all; do
    echo "[maps] $(date): $stage"
    snakemake -s Snakefile_novel $SNAKE_OPTS -- "$stage"
    rc=$?
    echo "[maps] $(date): $stage exit $rc"
    [ $rc -eq 0 ] || status=1
done

# 4. Report what landed.
echo "[maps] genome maps: $(ls results/genome_maps/*/genome_map_data.tsv 2>/dev/null | wc -l)"
echo "[maps] domain trees: $(ls results/phylo/families/*/concat_domains.faa 2>/dev/null | wc -l)"
exit $status
