#!/usr/bin/env bash
#SBATCH --job-name=conifer_cleanup
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --time=60
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/conifer_pest_virome/logs/cleanup_%j.out
#
# Reclaim regenerable intermediates under data/03.rnaviralSPAdes.
#
#     sbatch cleanup_assembly_intermediates.sh            # dry run
#     sbatch cleanup_assembly_intermediates.sh --apply    # delete
#
# What it removes, and why each is safe:
#
#   sorted.bam / .bai   Already declared temp() in rules/03.assembly.smk. These
#                       are leftovers from runs that predate that annotation,
#                       which Snakemake will not reclaim retrospectively. Only
#                       removed where coverage_stats.tsv exists, because
#                       contig_coverage_stats is the last consumer of the BAM
#                       and nothing downstream can recover it. (This used to
#                       check zero_coverage_contigs.txt, which filter_contigs
#                       writes BEFORE coverage is taken.)
#   contigs.amb/.ann/   bwa-mem2 index of contigs.fasta. contigs.fasta is itself
#   .bwt.2bit.64/.pac   temp() and has already been reclaimed, so the index is
#                       orphaned: it indexes a file that no longer exists.
#   virsorter2/iter-0   VirSorter2's working directory, ~99% of that folder. The
#                       declared outputs (final-viral-score.tsv,
#                       final-viral-combined.fa) and the boundary table, log and
#                       config are kept.
#
# It never touches filtered_contigs.fasta or zero_coverage_contigs.txt, which
# are the actual assembly products. Every candidate is also passed through
# scripts/retention_guard.py filter-deletions, so nothing listed in
# configs/retained_outputs.tsv (or a directory holding one) can be removed.
#
# SAFETY: refuses to delete anything while any job is running in this working
# directory, because the Snakemake driver spawns children that write here.
set -uo pipefail
ROOT=/home/people/fitzpatria/scratch/conifer_pest_virome
ASM="$ROOT/data/03.rnaviralSPAdes"
APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1
mkdir -p "$ROOT/logs"
echo "### host $(hostname) $(date) | apply=$APPLY"

# --- safety gate -----------------------------------------------------------
mine="${SLURM_JOB_ID:-none}"
active=0
for j in $(squeue -u "$USER" -h -o "%i" 2>/dev/null); do
    [ "$j" = "$mine" ] && continue
    wd=$(scontrol show job "$j" 2>/dev/null | sed -n 's/.*WorkDir=\([^ ]*\).*/\1/p')
    case "$wd" in
        *conifer_pest_virome*) echo "  job $j is still running in $wd"; active=$((active+1)) ;;
    esac
done
if [ "$active" -gt 0 ]; then
    echo "### REFUSING: $active job(s) still active in this working directory."
    echo "### Snakemake children write here; deleting now could break the run."
    exit 1
fi
echo "### no jobs active in $ROOT -- safe to proceed"

# --- survey ----------------------------------------------------------------
total=0
# =(): under `set -u`, bash 5.1 treats a declared-but-never-assigned array as
# unbound, so ${#VICTIMS[@]} crashed the script whenever nothing was cleanable.
declare -a VICTIMS=()
for d in "$ASM"/*/*/; do
    [ -d "$d" ] || continue
    # only clean assemblies whose coverage step finished: coverage_stats.tsv is
    # produced by the last consumer of the BAM
    if [ ! -s "$d/coverage_stats.tsv" ]; then
        echo "  SKIP (coverage stats absent, assembly may be incomplete): ${d#$ASM/}"
        continue
    fi
    for item in sorted.bam sorted.bam.bai contigs.amb contigs.ann \
                contigs.bwt.2bit.64 contigs.pac virsorter2/iter-0; do
        p="$d$item"
        if [ -e "$p" ]; then
            kb=$(du -sk "$p" 2>/dev/null | cut -f1)
            total=$((total + kb))
            VICTIMS+=("$p")
        fi
    done
done

# Retention guard: drop anything later analysis reads.
if [ "${#VICTIMS[@]}" -gt 0 ]; then
    SAFE_LIST=$(mktemp "${TMPDIR:-/tmp}/conifer_cleanup_safe.XXXXXX")
    if ! (cd "$ROOT" && printf '%s\0' "${VICTIMS[@]#$ROOT/}" \
            | python3 -I scripts/retention_guard.py filter-deletions) > "$SAFE_LIST"; then
        echo "### retention guard failed -- refusing to delete anything"
        [ -f "$SAFE_LIST" ] && rm -f -- "$SAFE_LIST"
        exit 1
    fi
    SAFE=()
    mapfile -d '' SAFE < "$SAFE_LIST"
    [ -f "$SAFE_LIST" ] && rm -f -- "$SAFE_LIST"
    VICTIMS=()
    for p in "${SAFE[@]}"; do VICTIMS+=("$ROOT/$p"); done
fi

echo "### ${#VICTIMS[@]} items, $((total / 1024 / 1024)) GB (size before retention filtering)"
if [ "${#VICTIMS[@]}" -eq 0 ]; then echo "### nothing to do"; exit 0; fi

if [ "$APPLY" -eq 0 ]; then
    printf '  would remove %s\n' "${VICTIMS[@]}" | sed "s|$ASM/||"
    echo "### dry run; pass --apply to delete"
    exit 0
fi

for p in "${VICTIMS[@]}"; do
    rm -rf "$p" && echo "  removed ${p#$ASM/}"
done
echo "### reclaimed about $((total / 1024 / 1024)) GB"
(cd "$ROOT" && python3 -I scripts/retention_guard.py audit)
echo "### surviving products per assembly:"
ls "$ASM"/*/*/filtered_contigs.fasta 2>/dev/null | wc -l | sed 's/^/  filtered_contigs.fasta: /'
ls "$ASM"/*/*/virsorter2/final-viral-score.tsv 2>/dev/null | wc -l | sed 's/^/  virsorter2 final-viral-score.tsv: /'
echo "### DONE $(date)"
