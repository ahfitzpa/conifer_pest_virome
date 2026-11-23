#!/bin/bash
#SBATCH --job-name=conifer_push
#SBATCH --output=logs/push_%j.log
#SBATCH --error=logs/push_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#
# 1. Push from a node: the login node's process limit breaks git's fork().
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
git push origin HEAD
echo "push exit $?"
git log --oneline -3
