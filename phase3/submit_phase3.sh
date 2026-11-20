#!/bin/bash
#SBATCH --job-name=conifer_p3
#SBATCH --output=logs/snakemake_p3_%j.log
#SBATCH --error=logs/snakemake_p3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie

# Phase 3 scheduler job. Every rule is submitted as its own SLURM job by
# --executor slurm, so the real compute happens there, not here.
#
# Run from conifer_pest_virome/phase3. Phase 2 is read as finished output from
# .. (the project root; override with PHASE2_DIR), and nothing here writes into it.

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
# Keep caches off $HOME: it is quota-limited on this cluster and a full home
# stops Snakemake before it starts.
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip

source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 \
 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda \
 --conda-frontend conda --rerun-incomplete \
 --rerun-triggers mtime --keep-going \
 --default-resources slurm_partition=shared runtime=180 mem_mb=0"

# Fail early and loudly if phase2 has not produced what phase3 consumes,
# rather than surfacing it as a confusing MissingInputException per family.
PHASE2="${PHASE2_DIR:-..}"
missing=0
for f in "$PHASE2/results/ictv/families.tsv" \
         "$PHASE2/results/reporting/candidate_contigs.tsv" \
         "$PHASE2/results/reporting/metacerberus_orfs.faa" \
         "$PHASE2/references/ictv_refs.faa"; do
    [ -e "$f" ] || { echo "[ERROR] phase2 prerequisite missing: $f" >&2; missing=1; }
done
[ "$missing" -eq 0 ] || { echo "[ERROR] Run phase2 to completion first." >&2; exit 1; }

echo "[SLURM INFO] Phase 3: ML phylogenetics + conserved domains + diagnostic residues"
snakemake -s Snakefile_phase3 $SNAKE_OPTS phylo_all

echo "[SLURM INFO] Phase 3: genome maps"
snakemake -s Snakefile_phase3 $SNAKE_OPTS genome_map_all

echo "[SLURM INFO] Phase 3: recombination screening"
snakemake -s Snakefile_phase3 $SNAKE_OPTS recomb_all

echo "[SLURM INFO] Phase 3: ICTV demarcation tables"
snakemake -s Snakefile_phase3 $SNAKE_OPTS demarcation_all

snakemake -s Snakefile_phase3 --report report_phase3.html --jobs 1
