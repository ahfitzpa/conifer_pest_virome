#!/bin/bash
#SBATCH --job-name=poster_figs
#SBATCH --output=logs/poster_figs_%j.log
#SBATCH --error=logs/poster_figs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --partition=shared
#
# Poster figures from whatever the phylogeny stage has finished: one tree per
# family with candidate tips coloured by host, plus a host x family matrix.
# Safe to re-run while the workflow is still going -- families without a tree
# yet are skipped and named, so the set fills in as the run progresses.
#
# SERRATUS_NOVEL points at the other project's analysis_novel so the Hepeviridae
# tree is drawn into its own figures/poster in the same pass.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs

module purge
module load snakemake/8.30.0
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"

# 1: an interpreter with matplotlib + biopython. The rule envs already built by
# this workflow have both; fall back to the declared env if none is present.
PY=""
for e in .snakemake/conda/*_; do
    [ -x "$e/bin/python" ] || continue
    if "$e/bin/python" -c "import matplotlib, Bio.Phylo" 2>/dev/null; then
        PY="$e/bin/python"; break
    fi
done
[ -n "$PY" ] || { echo "[figs] ERROR: no env with matplotlib + biopython" >&2; exit 1; }
echo "[figs] python: $PY"

# 2: this project
"$PY" -I scripts/poster_figures.py \
    --results results/phylogeny --out-dir figures/poster
rc=$?

# 3: the other project's analysis_novel, if it is reachable
SERRATUS_NOVEL="${SERRATUS_NOVEL:-/scratch/fitzpatria/Serratus/workflow/analysis_novel}"
if [ -d "$SERRATUS_NOVEL/results/phylogeny" ]; then
    echo "[figs] --- $SERRATUS_NOVEL ---"
    "$PY" -I scripts/poster_figures.py \
        --results "$SERRATUS_NOVEL/results/phylogeny" \
        --out-dir "$SERRATUS_NOVEL/figures/poster"
    rc2=$?
    [ "$rc" -eq 0 ] && rc=$rc2
fi

echo "[figs] exit $rc"
echo "[figs] conifer:  $(ls figures/poster/*.png 2>/dev/null | wc -l) png"
echo "[figs] serratus: $(ls "$SERRATUS_NOVEL"/figures/poster/*.png 2>/dev/null | wc -l) png"
exit $rc
