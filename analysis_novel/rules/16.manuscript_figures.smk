# =============================================================================
# Stage 16: Reproducible Manuscript Figures (Figs. 1-3 + Structural Gallery)
# =============================================================================
# Declarative, zero-OpenBLAS pure-Python SVG figure rules that rebuild
# automatically whenever upstream ESM-RdRp or Snowball tables change.
#
# Run directly:
#   snakemake -s Snakefile_novel --cores 2 -- manuscript_figures_all
# =============================================================================

import os

MS_FIG_DIR = "figures/manuscript"
SNOWBALL_DIR = os.path.join(PHASE2, "snowball/results/snowball")


rule manuscript_figures_all:
    input:
        fig1 = f"{MS_FIG_DIR}/fig1_workflow_and_map.svg",
        fig2 = f"{MS_FIG_DIR}/fig2_esm_rdrp_validation.svg",
        fig3 = f"{MS_FIG_DIR}/fig3_snowball_sharing_matrix.svg",


rule fig1_workflow_and_map:
    """Fig. 1: Closed-loop 4-stage discovery workflow + Holarctic 5-species map."""
    input:
        host_geo = "configs/host_geography.tsv",
    output:
        svg = f"{MS_FIG_DIR}/fig1_workflow_and_map.svg",
    params:
        manifest = f"{SNOWBALL_DIR}/manifest_combined.csv",
    log:
        "logs/figures/fig1_workflow_and_map.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.svg})
        MAN_ARG=""
        if [ -s "{params.manifest}" ]; then
            MAN_ARG="--manifest {params.manifest}"
        fi
        python3 scripts/plot_fig1_workflow_and_map.py \
            --host-geo "{input.host_geo}" \
            $MAN_ARG \
            --out-svg "{output.svg}" \
            > "{log}" 2>&1
        """


rule fig2_esm_rdrp_validation:
    """Fig. 2: Stage 15 ESMFold + Foldseek structural triage & RT negative control."""
    input:
        validation = f"{ESM_DIR}/esm_rdrp_validation.tsv",
    output:
        svg = f"{MS_FIG_DIR}/fig2_esm_rdrp_validation.svg",
    params:
        manifest = f"{SNOWBALL_DIR}/manifest_combined.csv",
    log:
        "logs/figures/fig2_esm_rdrp_validation.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.svg})
        MAN_ARG=""
        if [ -s "{params.manifest}" ]; then
            MAN_ARG="--manifest {params.manifest}"
        fi
        python3 scripts/plot_fig2_esm_rdrp.py \
            --validation-tsv "{input.validation}" \
            $MAN_ARG \
            --out-validation-svg "{output.svg}" \
            > "{log}" 2>&1
        """


rule fig3_snowball_sharing_matrix:
    """Fig. 3: Snowball replication & pan-Holarctic cross-host sharing matrix."""
    input:
        hits = f"{SNOWBALL_DIR}/snowball_vs_seed_blastx.tsv",
        esm  = f"{ESM_DIR}/esm_rdrp_validation.tsv",
    output:
        svg = f"{MS_FIG_DIR}/fig3_snowball_sharing_matrix.svg",
    params:
        manifest = f"{SNOWBALL_DIR}/manifest_combined.csv",
    log:
        "logs/figures/fig3_snowball_sharing_matrix.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.svg})
        MAN_ARG=""
        if [ -s "{params.manifest}" ]; then
            MAN_ARG="--manifest {params.manifest}"
        fi
        python3 scripts/plot_fig3_snowball_matrix.py \
            --hits-tsv "{input.hits}" \
            --esm-validation "{input.esm}" \
            $MAN_ARG \
            --out-svg "{output.svg}" \
            > "{log}" 2>&1
        """
