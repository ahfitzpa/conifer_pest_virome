#!/bin/bash
#SBATCH --job-name=conifer_esm_rdrp
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-12:00:00
#SBATCH --output=logs/esm_rdrp_%j.out
#SBATCH --error=logs/esm_rdrp_%j.err
#
# GPU SLURM launcher for Stage 15 (ESMFold + Foldseek structural validation of
# medium/high-confidence, deeply divergent or family-unresolved RdRp domains).
#
# Usage (from analysis_novel/):
#   sbatch submit_esm_rdrp_gpu.sh
#
# Optional: override GPU partition at submit time if your cluster uses another
# name (e.g., gpu, gpu-a100, gpu-v100):
#   sbatch --partition=gpu submit_esm_rdrp_gpu.sh
#
# One-time setup of Foldseek's official PDB database (run interactively or via
# sbatch with SETUP_PDB_DB=1):
#   SETUP_PDB_DB=1 sbatch submit_esm_rdrp_gpu.sh

set -euo pipefail

REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR:-$PWD}/.." && pwd)"
ANALYSIS_DIR="${REPO_ROOT}/analysis_novel"
cd "${ANALYSIS_DIR}"

mkdir -p logs "${REPO_ROOT}/.cache/huggingface" "${REPO_ROOT}/.cache/torch" "${REPO_ROOT}/references/foldseek"

# Keep multi-GB ESMFold weights and PyTorch caches in scratch/repo cache, never $HOME
export HF_HOME="${HF_HOME:-${REPO_ROOT}/.cache/huggingface}"
export TORCH_HOME="${TORCH_HOME:-${REPO_ROOT}/.cache/torch}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/.cache}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

export SETUP_PDB_DB="${SETUP_PDB_DB:-0}"
FOLDSEEK_DB="${FOLDSEEK_DB:-${REPO_ROOT}/references/foldseek/pdb}"

echo "[esm_rdrp] Host: $(hostname)"
echo "[esm_rdrp] CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "[esm_rdrp] FOLDSEEK_DB=${FOLDSEEK_DB} (SETUP_PDB_DB=${SETUP_PDB_DB})"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
fi

snakemake \
    -s Snakefile_novel \
    --use-conda \
    --conda-frontend "${CONDA_FRONTEND:-mamba}" \
    --cores "${SLURM_CPUS_PER_TASK:-8}" \
    --resources gpu=1 \
    --rerun-incomplete \
    --config "esm_rdrp={enabled: true, foldseek_db: '${FOLDSEEK_DB}', hf_cache_dir: '${HF_HOME}'}" \
    esm_rdrp_all
