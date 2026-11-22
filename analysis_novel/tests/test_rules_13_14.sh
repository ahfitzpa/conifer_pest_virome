#!/bin/bash
#SBATCH --job-name=novel_13_14
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=30
#SBATCH --output=/scratch/fitzpatria/Serratus/workflow/analysis_novel/logs/test_13_14_%j.log
# 1: parse analysis_novel with rules 09-14 included
# 2: exercise the new scripts against the host_range fixture (known answers)
set -uo pipefail
W=/scratch/fitzpatria/Serratus/workflow
cd $W/analysis_novel
module purge; module load snakemake/8.30.0
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export NOVEL_MANIFEST=../common/test_data/manifest_novel_test.csv
status=0

echo "### parse + rule list"
snakemake -s Snakefile_novel --list > /tmp/rules.$$ 2>&1
if [ $? -eq 0 ]; then
    echo "  PASS: parses with rules 09-14 ($(grep -c . /tmp/rules.$$) rules)"
    grep -E "^(codon_all|signatures_all|predict_orfs|select_orf1|codon_align_orf1|hyphy_fubar|hyphy_meme|selection_summary|composition_signatures|prune_tree_to_alignment)$" /tmp/rules.$$ | sed 's/^/    /'
else
    echo "  FAIL: parse error"; cat /tmp/rules.$$; status=1
fi
rm -f /tmp/rules.$$

echo "### dry run: codon_all then signatures_all"
for t in codon_all signatures_all novel_all; do
    out=$(snakemake -s Snakefile_novel -n --quiet rules -- $t 2>&1)
    if [ $? -eq 0 ]; then echo "  PASS: $t plans $(echo "$out" | awk '/^total/{print $2}') jobs"
    else echo "  FAIL: $t"; echo "$out" | tail -5; status=1; fi
done

echo "### scripts against the fixture"
export PATH=/scratch/fitzpatria/calicityper_env/bin:$PATH
F=$W/analysis_host_range/tests/fixtures
T=/scratch/fitzpatria/tmp/novel_13_14_${SLURM_JOB_ID:-local}; mkdir -p $T
trap 'rm -rf $T' EXIT

python3 -I $W/common/scripts/prune_tree.py --tree $F/tree_with_support.newick \
    --alignment $F/orf_cds.fasta --out-tree $T/pruned.newick \
    --out-report $T/pruned.tsv --strip-support --min-taxa 4 > $T/prune.log 2>&1
if [ $? -eq 0 ] && ! grep -qE '[0-9]+/[0-9]+' $T/pruned.newick; then
    echo "  PASS: prune_tree stripped 98.5/95 support labels"
else
    echo "  FAIL: support labels survived or prune failed"; tail -3 $T/prune.log; status=1
fi
n=$(grep -o '[A-CU][0-9]:' $T/pruned.newick | wc -l)
[ "$n" -eq 9 ] && echo "  PASS: all 9 fixture tips retained" || { echo "  FAIL: $n tips"; status=1; }

python3 -I $W/common/scripts/composition_signatures.py --cds $F/orf_cds.fasta \
    --out-rscu $T/rscu.tsv --out-dinuc $T/dinuc.tsv --min-codons 50 > $T/comp.log 2>&1
if [ $? -eq 0 ] && [ "$(awk 'NR>1' $T/rscu.tsv | wc -l)" -eq 9 ]; then
    echo "  PASS: composition_signatures wrote 9 genomes"
else
    echo "  FAIL: composition_signatures"; tail -3 $T/comp.log; status=1
fi
# the fixture is almost all Leu CTT with CTC in four tips: RSCU must differ
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} NR>1{print $1, $h["CTT"], $h["CTC"]}' $T/rscu.tsv | head -3 | sed 's/^/    /'

echo "### status $status"
exit $status
