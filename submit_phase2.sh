#!/bin/bash
#SBATCH --job-name=conifer_p2
#SBATCH --output=logs/snakemake_p2_%j.log
#SBATCH --error=logs/snakemake_p2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
# Pass your own notification address at submit time (SBATCH directives
# can't read a config file): sbatch --mail-user=you@example.org submit_phase2.sh

# This job is just the Snakemake scheduler -- it submits every rule as
# its own SLURM job via --executor slurm, so the actual compute happens
# in those per-rule jobs, not here. See ENGINEERING_NOTES.md § submit_phase2.sh.

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=mamba
# Keep caches off $HOME: it is quota-limited on this cluster, and a full home
# stops Snakemake before it starts -- its source cache is created under
# ~/.cache/snakemake and the failure is an OSError, not a workflow error.
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip

source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

# sra_fetch=4: caps concurrent SRA downloads regardless of the --jobs
# ceiling above, which is shared across every rule type, not fetch-specific.
# See ENGINEERING_NOTES.md § fetch_and_dump_sra_paired -- this starting
# value hasn't been through the sacct wall-clock-vs-N benchmark the review
# calls for yet.
SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda --conda-frontend conda --rerun-incomplete --rerun-triggers mtime params input --keep-going --default-resources slurm_partition=shared runtime=180 mem_mb=0 --resources sra_fetch=4"

echo "[SLURM INFO] Starting Phase 2: Assembly & Analysis (fanning out across the shared partition via --executor slurm)"

# Stage 1: finish VirSorter2 before anything consumes its score.
#
# annotate_viral_positives takes the VirSorter2 score as a params path rather
# than an input, so that a VirSorter2 failure can no longer block a sample's
# taxonomy_merged.tsv or, through it, MetaCerberus. The cost of that is the
# ordering guarantee: as a params path it no longer forces annotate to wait, so
# annotate could otherwise run first and record NA for a sample whose VirSorter2
# was still running -- and because the merge would then be complete, a later
# rerun would not revisit it and the score would be lost silently.
#
# Running the VirSorter2 targets to completion first restores the ordering
# without restoring the hard gate. It is allowed to fail (|| true): any sample
# VirSorter2 cannot do is listed in logs/virsorter2_failures.tsv and simply gets
# NA in stage 2.
echo "[SLURM INFO] Stage 1/3: VirSorter2 (best-effort; failures are logged, not fatal)"
# --notemp: contig_coverage_stats is outside this partial run, so without it
# Snakemake reclaims each sorted.bam once filter_contigs is done -- which is
# how all 21 BAMs were lost before coverage in run 1173405 (2026-11).
# Stage 2 runs without it and reclaims temp files once every consumer is
# done. rules/99.retention.smk refuses an --until run without --notemp.
snakemake -s Snakefile_phase2 $SNAKE_OPTS --notemp --until run_virsorter2 || true

echo "[SLURM INFO] Stage 2/3: full workflow"
snakemake -s Snakefile_phase2 $SNAKE_OPTS || true

# Stage 3: aggregate whatever stage 2 produced, so one unassemblable sample
# cannot block MetaCerberus. Separate invocation (a DAG is built once) and an
# env var, not --config, which this Snakefile discards.
echo "[SLURM INFO] Stage 3/3: master database + MetaCerberus over the samples that completed"
PHASE2_ALLOW_PARTIAL_MASTER=true snakemake -s Snakefile_phase2 $SNAKE_OPTS \
    --rerun-triggers mtime \
    -- results/master_viral_database.fasta results/taxonomy/excluded_samples.tsv \
    results/annotation/metacerberus

# Reporting (kept local -- trivial, no need to submit these)
snakemake -s Snakefile_phase2 --dag --jobs 1 | dot -Tsvg > logs/phase2_dag.svg
snakemake -s Snakefile_phase2 --report report_phase2.html --jobs 1
