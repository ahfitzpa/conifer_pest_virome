#!/bin/bash
#SBATCH --job-name=fetch_refs
#SBATCH --output=logs/fetch_refs_%j.log
#SBATCH --error=logs/fetch_refs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --partition=shared
#
# Build the ICTV reference sets for the families this pathway analyses.
# Writes references/ictv_refs_all.faa (+ meta) for the RdRp marker phylogeny
# and references/nt/<Family>_refs.fasta for the shared nucleotide phylogeny.
# Re-runs are cheap: raw GenBank is cached under references/ncbi/.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs

FAMILIES="${FAMILIES:-configs/families.yaml}"
PER_FAMILY="${PER_FAMILY:-12}"
REF=../references

# 1: a python with Biopython -- the base module has none, so fall back to an
# env snakemake has already built for this workflow
PY=""
if python3 -c "import Bio" 2>/dev/null; then
    PY="python3"
else
    for e in .snakemake/conda/*_; do
        for cand in "$e"/bin/python3; do
            [ -x "$cand" ] || continue
            if "$cand" -c "import Bio" 2>/dev/null; then PY="$cand"; break 2; fi
        done
    done
fi
if [ -z "$PY" ]; then
    echo "[fetch_refs] ERROR: no python with Biopython found" >&2
    exit 1
fi
echo "[fetch_refs] python: $PY"

# 2: VMR spreadsheet -> TSV, so the fetch needs no openpyxl (skipped if current)
if [ ! -s "$REF/vmr_msl41.tsv" ]; then
    python3 -I scripts/vmr_xlsx_to_tsv.py || {
        echo "[fetch_refs] ERROR: VMR conversion failed" >&2; exit 1; }
fi

# 3: fetch exemplars and extract proteins + genomes
"$PY" -I scripts/fetch_refs_vmr.py \
    --vmr "$REF/vmr_msl41.tsv" \
    --families "$FAMILIES" \
    --per-family "$PER_FAMILY" \
    --out-faa "$REF/ictv_refs_all.faa" \
    --out-meta "$REF/ictv_refs_all_meta.tsv" \
    --out-nt-dir "$REF/nt"
rc=$?

# 4: report
echo "[fetch_refs] exit $rc"
if [ "$rc" -eq 0 ]; then
    echo "[fetch_refs] proteins: $(grep -c '^>' "$REF/ictv_refs_all.faa" || true)"
    echo "[fetch_refs] genome sets: $(ls "$REF/nt" 2>/dev/null | wc -l)"
fi
exit $rc
