#!/usr/bin/env bash
#SBATCH --job-name=mbio_analyses
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/mbio_analyses_%j.out
#SBATCH --error=logs/mbio_analyses_%j.err
# =============================================================================
# SLURM wrapper for mBio evolutionary, enrichment, benchmark & figure analyses
# =============================================================================
# Runs cleanly on a compute node (avoiding login-node OpenBLAS/thread limits):
#   1. Open-host Fisher's exact enrichment across 3.8M SRA runs (snowball/)
#   2. Pan-Holarctic Ips Polycipiviridae phylogeography & Nei-Gojobori dN/dS
#   3. Stage 17 benchmarks & refreshed manuscript SVG figures (Figs. 1-6)
#
# Submit from analysis_novel/:
#   sbatch submit_mbio_analyses.sh
# =============================================================================
set -euo pipefail

module load snakemake/8.30.0 2>/dev/null || true
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

NOVEL_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
REPO_ROOT="$(cd "${NOVEL_DIR}/.." && pwd)"
SNOWBALL_DIR="${REPO_ROOT}/snowball"
mkdir -p "${NOVEL_DIR}/logs"

echo "[$(date -Is)] Step 1/3: Open-host Fisher enrichment (snowball/)..."
(
    cd "${SNOWBALL_DIR}"
    if [ -s "results/open_host/open_host_candidates.tsv" ]; then
        snakemake -s Snakefile_snowball \
            --use-conda \
            --conda-frontend conda \
            --conda-prefix "${SNOWBALL_DIR}/.snakemake/conda" \
            --cores "${SLURM_CPUS_PER_TASK:-4}" \
            --rerun-triggers mtime \
            --allowed-rules open_host_enrichment \
            -- open_host_enrichment
        echo "--- serratus_denominators.tsv ---"
        cat results/open_host/serratus_denominators.tsv
        echo "--- open_host_enrichment.tsv ---"
        cat results/open_host/open_host_enrichment.tsv
    else
        echo "[WARN] results/open_host/open_host_candidates.tsv not found; skipping Step 1."
    fi
)

echo "[$(date -Is)] Step 2/3: Ips Polycipiviridae phylogeography & dN/dS..."
cd "${NOVEL_DIR}"
snakemake -s Snakefile_novel \
    --cores "${SLURM_CPUS_PER_TASK:-4}" \
    --rerun-triggers mtime \
    --allowed-rules benchmark_polycipiviridae_evolution \
    --forcerun benchmark_polycipiviridae_evolution \
    -- results/benchmarks/05_polycipiviridae_evolution/polycipiviridae_evolution.svg

echo "--- polycipiviridae_genomes_summary.tsv ---"
cat results/benchmarks/05_polycipiviridae_evolution/polycipiviridae_genomes_summary.tsv
echo "--- polycipiviridae_phylogeography_summary.tsv ---"
cat results/benchmarks/05_polycipiviridae_evolution/polycipiviridae_phylogeography_summary.tsv

echo "[$(date -Is)] Step 3/3: Refreshing manuscript SVG figures (Figs. 1-6)..."
snakemake -s Snakefile_novel \
    --cores "${SLURM_CPUS_PER_TASK:-4}" \
    --rerun-triggers mtime \
    --allowed-rules fig1_workflow_and_map fig2_esm_rdrp fig3_snowball_matrix fig4_polycipiviridae fig5_iflaviridae fig6_amarillo_negssrna manuscript_figures_phylo manuscript_figures_complete \
    --forcerun fig3_snowball_matrix fig4_polycipiviridae fig5_iflaviridae fig6_amarillo_negssrna \
    -- manuscript_figures_complete

echo "[$(date -Is)] All mBio analyses complete!"
