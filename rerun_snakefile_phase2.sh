#!/bin/bash
#SBATCH --job-name=conifer_ips
#SBATCH --output=logs/rerun_ips_%j.log
#SBATCH --error=logs/rerun_ips_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=5-00:00:00
# Scheduler only: each rule is submitted as its own SLURM job.
# Targeted run for the 8 Ips_confusus libraries (post dynamic-kmer fix)
# + the 2 Locusta_migratoria libraries whose RdRpCATCH jobs timed out,
# followed by rebuilding the master database and MetaCerberus.

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=mamba
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PIP_CACHE_DIR=/scratch/fitzpatria/.cache/pip

source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda --conda-frontend conda --rerun-incomplete --rerun-triggers mtime --keep-going --default-resources slurm_partition=shared runtime=180 mem_mb=0 --set-resources run_rdrpcatch:runtime=720"

IPS_SAMPLES=(
    SRR12261679 SRR12261680 SRR12261681 SRR12261682
    SRR12261683 SRR12261684 SRR12261685 SRR12261686
)
LOCUSTA_RETRY=(SRR9898676 SRR9898677)

STAGE1_TARGETS=()
STAGE2_TARGETS=()

for s in "${IPS_SAMPLES[@]}"; do
    STAGE1_TARGETS+=("data/03.rnaviralSPAdes/Ips_confusus/${s}/virsorter2/final-viral-score.tsv")
    STAGE2_TARGETS+=(
        "data/03.assembly_stats/Ips_confusus/${s}/${s}_stats.txt"
        "data/05.final_selection/Ips_confusus/${s}/taxonomy_merged.tsv"
        "data/06.rdrp/Ips_confusus/${s}/${s}_rdrp_domains.faa"
        "data/04.taxonomy/Ips_confusus/${s}/.purged"
    )
done

for s in "${LOCUSTA_RETRY[@]}"; do
    STAGE2_TARGETS+=(
        "data/05.final_selection/Locusta_migratoria/${s}/taxonomy_merged.tsv"
        "data/06.rdrp/Locusta_migratoria/${s}/${s}_rdrp_domains.faa"
        "data/04.taxonomy/Locusta_migratoria/${s}/.purged"
    )
done

echo "[INFO] Stage 1/3: Ips_confusus assembly + VirSorter2 (8 samples)"
snakemake -s Snakefile_phase2 $SNAKE_OPTS -- "${STAGE1_TARGETS[@]}" || true

echo "[INFO] Stage 2/3: Ips_confusus taxonomy/RdRp (8) + Locusta RdRpCATCH retry (2)"
snakemake -s Snakefile_phase2 $SNAKE_OPTS -- "${STAGE2_TARGETS[@]}" || true

echo "[INFO] Stage 3/3: Rebuild master database + MetaCerberus over all completed samples"
PHASE2_ALLOW_PARTIAL_MASTER=true snakemake -s Snakefile_phase2 $SNAKE_OPTS \
    -- results/master_viral_database.fasta \
       results/taxonomy/master_metadata.csv \
       results/taxonomy/excluded_samples.tsv \
       results/annotation/metacerberus
echo "[INFO] exit=$?"