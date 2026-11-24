#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=conifer_rebuild
#SBATCH --output=logs/rebuild_%j.log
#SBATCH --error=logs/rebuild_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
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
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"

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
 --default-resources slurm_partition=shared runtime=180 mem_mb=0 slurm_extra='--exclude=sonic59'"

# Every target is passed after `--`. Without it, --default-resources (and
# --quiet) swallow the target: that is why phase3's submit_phase3.sh failed all
# four stages in 7 seconds with "not enough values to unpack".
# 1. run_stage <target> [rule-to-force ...]: --rerun-triggers mtime ignores a
#    code change, so the rules whose code was corrected are forced by name.
run_stage() {
    local target="$1"; shift
    local force=""
    [ $# -gt 0 ] && force="--forcerun $*"
    echo "[rebuild] $(date): $target ${force:+($force)}"
    snakemake -s Snakefile_novel $SNAKE_OPTS $force -- "$target"
    local rc=$?
    echo "[rebuild] $(date): $target exit $rc"
    return $rc
}

status=0
# 2. select_contigs is forced once here: every selected.fasta on disk predates
#    the reference-length cap, so Mimiviridae still carried 1.2-1.5 Mb genomes
#    against a 9 kb candidate and killed FAMSA with std::bad_alloc.
run_stage phylogeny_all select_contigs || status=1
# 3. codon and signatures now cover the RNA families only.
run_stage codon_all || status=1
run_stage signatures_all || status=1
# 4. concat_prepare is forced for the same reason: its sets were built before
#    the contig-id and assigned_family fix and hold no candidates.
run_stage phylo_all concat_prepare || status=1
for stage in genome_map_all recomb_all demarcation_all; do
    run_stage "$stage" || status=1
done

snakemake -s Snakefile_novel --report report_novel.html --jobs 1 || true
echo "[rebuild] $(date): done, status $status"
exit $status
