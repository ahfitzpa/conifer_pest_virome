#!/bin/bash
#SBATCH --job-name=novel_analysis
#SBATCH --output=logs/novel_%j.log
#SBATCH --error=logs/novel_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#
# Scheduler only: every rule is submitted as its own SLURM job via
# --executor slurm. Pass NOVEL_MANIFEST to run a different sample set, e.g.
#   sbatch --export=ALL,NOVEL_MANIFEST=../common/test_data/manifest_novel_test.csv submit_novel.sh
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs tmp

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$PWD/.cache/pip}"
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

if [ -n "${PHASE2_DIR:-}" ]; then
    PHASE2="$PHASE2_DIR"
elif [ -f "../configs/paths.yaml" ]; then
    PHASE2=".."
else
    PHASE2="../phase2"
fi
for f in "$PHASE2/results/master_viral_database.fasta" \
         "$PHASE2/results/taxonomy/master_metadata.csv"; do
    [ -e "$f" ] || { echo "[ERROR] phase2 prerequisite missing: $f" >&2; exit 1; }
done

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0"

# Every target is passed after `--`. Without it, --default-resources (and
# --quiet) swallow the target: that is why phase3's submit_phase3.sh failed all
# four stages in 7 seconds with "not enough values to unpack".
run_stage() {
    echo "[novel] $(date): $1"
    snakemake -s Snakefile_novel $SNAKE_OPTS -- "$1"
    local rc=$?
    echo "[novel] $(date): $1 exit $rc"
    return $rc
}

status=0
# 1: per-scenario phylogeny (select -> FAMSA -> ClipKIT -> IQ-TREE2)
run_stage phylogeny_all || status=1
# 2: codon alignments for dN/dS (Prodigal -> ORF1 -> MAFFT -> threading)
run_stage codon_all || status=1
# 3: genomic signatures (FUBAR/MEME, composition)
run_stage signatures_all || status=1
# 4: ported phase3 stages, which read phase2's ICTV and reporting outputs
for stage in phylo_all genome_map_all recomb_all demarcation_all; do
    run_stage "$stage" || status=1
done

snakemake -s Snakefile_novel --report report_novel.html --jobs 1 || true
echo "[novel] $(date): done, status $status"
exit $status
