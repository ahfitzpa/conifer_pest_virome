# =============================================================================
# Stage 17: Reproducible Method Benchmarks & Out-of-Order Validation
# =============================================================================
# Organises evaluation outputs into numbered subdirectories under
# results/benchmarks/ so each analysis is immediately identifiable:
#
#   results/benchmarks/02_structural_triage/    Ablation, TM sweep, ROC/PR SVG,
#                                               & out-of-order L. migratoria check
#   results/benchmarks/03_snowball_gain/        Retrieval yield, completion gain,
#                                               host specificity & LOLO check
#   results/benchmarks/04_discovery_and_compute/ Divergence-bin sensitivity &
#                                                Snakemake CPU/GPU/RAM summary
#
# Run directly:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 snakemake -s Snakefile_novel \
#     --cores 2 --rerun-triggers mtime -- benchmarks_all
# =============================================================================

import os

BENCH_OUT_DIR = "results/benchmarks"
BENCH_ROUND1_ESM = config.get("benchmarks", {}).get(
    "round1_validation_tsv",
    "benchmark_data/esm_rdrp_validation_round1_284.tsv",
)


rule benchmarks_all:
    input:
        triage_oos = f"{BENCH_OUT_DIR}/02_structural_triage/out_of_order_validation.tsv",
        triage_svg = f"{BENCH_OUT_DIR}/02_structural_triage/roc_pr_curves.svg",
        snowball   = f"{BENCH_OUT_DIR}/03_snowball_gain/snowball_gain_summary.tsv",
        discovery  = f"{BENCH_OUT_DIR}/04_discovery_and_compute/divergence_bin_sensitivity.tsv",
        compute    = f"{BENCH_OUT_DIR}/04_discovery_and_compute/compute_cost_summary.tsv",


rule benchmark_structural_triage:
    """Stage 17.2 + 17.5: Structural triage ablation, TM sweep & L. migratoria OOS check."""
    input:
        validation = f"{ESM_DIR}/esm_rdrp_validation.tsv",
    output:
        cm_col  = f"{BENCH_OUT_DIR}/02_structural_triage/confusion_matrix_coleoptera.tsv",
        cm_loc  = f"{BENCH_OUT_DIR}/02_structural_triage/confusion_matrix_out_of_order_locusta.tsv",
        ablat   = f"{BENCH_OUT_DIR}/02_structural_triage/ablation_summary.tsv",
        sweep   = f"{BENCH_OUT_DIR}/02_structural_triage/threshold_sweep.tsv",
        oos     = f"{BENCH_OUT_DIR}/02_structural_triage/out_of_order_validation.tsv",
        svg     = f"{BENCH_OUT_DIR}/02_structural_triage/roc_pr_curves.svg",
    params:
        out_dir    = f"{BENCH_OUT_DIR}/02_structural_triage",
        manifest   = f"{SNOWBALL_DIR}/manifest_combined.csv",
        round1_tsv = BENCH_ROUND1_ESM,
        gt_tsv     = config.get("esm_rdrp", {}).get("ground_truth_tsv", ""),
    log:
        "logs/benchmarks/02_structural_triage.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.out_dir}"
        VAL_ARGS=""
        [ -s "{params.round1_tsv}" ] && VAL_ARGS="--validation-tsv {params.round1_tsv}"
        VAL_ARGS="$VAL_ARGS --validation-tsv {input.validation}"
        MAN_ARG=""
        [ -s "{params.manifest}" ] && MAN_ARG="--manifest {params.manifest}"
        GT_ARG=""
        [ -n "{params.gt_tsv}" ] && [ -s "{params.gt_tsv}" ] && GT_ARG="--ground-truth-tsv {params.gt_tsv}"
        python3 scripts/benchmark_structural_triage.py \
            $VAL_ARGS \
            $MAN_ARG \
            $GT_ARG \
            --out-dir "{params.out_dir}" \
            > "{log}" 2>&1
        """


rule benchmark_snowball_gain:
    """Stage 17.3: Quantify Round-2 snowballing gain, host specificity, and LOLO."""
    input:
        hits = f"{SNOWBALL_DIR}/snowball_vs_seed_blastx.tsv",
        esm  = f"{ESM_DIR}/esm_rdrp_validation.tsv",
    output:
        summary = f"{BENCH_OUT_DIR}/03_snowball_gain/snowball_gain_summary.tsv",
        lineage = f"{BENCH_OUT_DIR}/03_snowball_gain/per_lineage_recovery.tsv",
    params:
        out_dir  = f"{BENCH_OUT_DIR}/03_snowball_gain",
        manifest = f"{SNOWBALL_DIR}/manifest_combined.csv",
    log:
        "logs/benchmarks/03_snowball_gain.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.out_dir}"
        MAN_ARG=""
        [ -s "{params.manifest}" ] && MAN_ARG="--manifest {params.manifest}"
        python3 scripts/benchmark_snowball_gain.py \
            --hits-tsv "{input.hits}" \
            --esm-validation "{input.esm}" \
            $MAN_ARG \
            --out-dir "{params.out_dir}" \
            > "{log}" 2>&1
        """


rule benchmark_discovery_and_compute:
    """Stage 17.4: Divergence-bin discovery sensitivity & Snakemake compute cost."""
    input:
        esm = f"{ESM_DIR}/esm_rdrp_validation.tsv",
    output:
        div  = f"{BENCH_OUT_DIR}/04_discovery_and_compute/divergence_bin_sensitivity.tsv",
        comp = f"{BENCH_OUT_DIR}/04_discovery_and_compute/compute_cost_summary.tsv",
    params:
        out_dir     = f"{BENCH_OUT_DIR}/04_discovery_and_compute",
        round1_tsv  = BENCH_ROUND1_ESM,
        ps_tsv      = config.get("benchmarks", {}).get("pebblescout_summary_tsv", ""),
        phase2_bdir = os.path.join(PHASE2, "benchmarks"),
    log:
        "logs/benchmarks/04_discovery_and_compute.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.out_dir}"
        VAL_ARGS=""
        [ -s "{params.round1_tsv}" ] && VAL_ARGS="--esm-validation {params.round1_tsv}"
        VAL_ARGS="$VAL_ARGS --esm-validation {input.esm}"
        PS_ARG=""
        [ -n "{params.ps_tsv}" ] && [ -s "{params.ps_tsv}" ] && PS_ARG="--pebblescout-tsv {params.ps_tsv}"
        BDIR_ARGS="--benchmark-dir benchmarks"
        [ -d "{params.phase2_bdir}" ] && BDIR_ARGS="$BDIR_ARGS --benchmark-dir {params.phase2_bdir}"
        python3 scripts/benchmark_discovery_and_compute.py \
            $VAL_ARGS \
            $PS_ARG \
            $BDIR_ARGS \
            --out-dir "{params.out_dir}" \
            > "{log}" 2>&1
        """


rule benchmark_polycipiviridae_evolution:
    """Stage 17.6: Pan-Holarctic Ips Polycipiviridae phylogeography, ORF intactness & Nei-Gojobori dN/dS."""
    output:
        genomes = f"{BENCH_OUT_DIR}/05_polycipiviridae_evolution/polycipiviridae_genomes_summary.tsv",
        pairs   = f"{BENCH_OUT_DIR}/05_polycipiviridae_evolution/polycipiviridae_pairwise_dnds.tsv",
        summary = f"{BENCH_OUT_DIR}/05_polycipiviridae_evolution/polycipiviridae_phylogeography_summary.tsv",
        svg     = f"{BENCH_OUT_DIR}/05_polycipiviridae_evolution/polycipiviridae_evolution.svg",
    params:
        out_dir  = f"{BENCH_OUT_DIR}/05_polycipiviridae_evolution",
        manifest = f"{SNOWBALL_DIR}/manifest_combined.csv",
        fam_dir  = "results/phylogenetics/Polycipiviridae",
    log:
        "logs/benchmarks/05_polycipiviridae_evolution.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.out_dir}"
        FA_ARGS=""
        for f in "{params.fam_dir}"/*novel*.fasta "{params.fam_dir}"/*contigs*.fasta; do
            [ -s "$f" ] && FA_ARGS="$FA_ARGS --contigs-fasta $f"
        done
        if [ -z "$FA_ARGS" ]; then
            echo "[polycipiviridae_evolution] ERROR: No Polycipiviridae contig FASTA found in {params.fam_dir}" > "{log}"
            exit 1
        fi
        MAN_ARG=""
        [ -s "{params.manifest}" ] && MAN_ARG="--manifest {params.manifest}"
        python3 scripts/analyze_polycipiviridae_evolution.py \
            $FA_ARGS \
            $MAN_ARG \
            --out-dir "{params.out_dir}" \
            > "{log}" 2>&1
        """
