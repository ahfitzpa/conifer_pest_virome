# =============================================================================
# Stage 16: Reproducible Manuscript Figures (Figs. 1-6)
# =============================================================================
# Declarative, zero-OpenBLAS pure-Python SVG figure rules that rebuild
# automatically whenever upstream ESM-RdRp, Snowball, or Stage 09/10 outputs
# change.
#
# Run directly:
#   snakemake -s Snakefile_novel --cores 2 -- manuscript_figures_all       # Figs. 1-3
#   snakemake -s Snakefile_novel --cores 2 -- manuscript_figures_phylo     # Figs. 4-6
#   snakemake -s Snakefile_novel --cores 2 -- manuscript_figures_complete  # Figs. 1-6
# =============================================================================

import os

MS_FIG_DIR = "figures/manuscript"
SNOWBALL_DIR = os.path.join(PHASE2, "snowball/results/snowball")


rule manuscript_figures_all:
    input:
        fig1 = f"{MS_FIG_DIR}/fig1_workflow_and_map.svg",
        fig2 = f"{MS_FIG_DIR}/fig2_esm_rdrp_validation.svg",
        fig3 = f"{MS_FIG_DIR}/fig3_snowball_sharing_matrix.svg",


rule manuscript_figures_phylo:
    input:
        fig4 = f"{MS_FIG_DIR}/fig4_polycipiviridae.svg",
        fig5 = f"{MS_FIG_DIR}/fig5_iflaviridae.svg",
        fig6 = f"{MS_FIG_DIR}/fig6_amarillo_negssrna.svg",


rule manuscript_figures_complete:
    input:
        rules.manuscript_figures_all.input,
        rules.manuscript_figures_phylo.input,


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


rule fig4_polycipiviridae:
    """Fig. 4: Pan-Holarctic Ips Polycipiviridae AA phylogeny + ORF architecture."""
    input:
        tree = f"{PHYLO_DIR}/families/Polycipiviridae/tree.newick",
        meta = f"{PHYLO_DIR}/families/Polycipiviridae/taxa_metadata.tsv",
    output:
        svg = f"{MS_FIG_DIR}/fig4_polycipiviridae.svg",
    params:
        genome_map = f"{GM_DIR}/Polycipiviridae/genome_map_data.tsv",
        manifest   = f"{SNOWBALL_DIR}/manifest_combined.csv",
    log:
        "logs/figures/fig4_polycipiviridae.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.svg})
        MAN_ARG=""
        [ -s "{params.manifest}" ] && MAN_ARG="--manifest {params.manifest}"
        GM_ARG=""
        [ -s "{params.genome_map}" ] && GM_ARG="--genome-map {params.genome_map}"
        python3 scripts/plot_phylo_manuscript.py \
            --mode fig4 \
            --tree "{input.tree}" \
            --meta "{input.meta}" \
            $GM_ARG \
            $MAN_ARG \
            --out-svg "{output.svg}" \
            > "{log}" 2>&1
        """


rule fig5_iflaviridae:
    """Fig. 5: Closed-genome conifer-pest Iflaviridae AA phylogeny."""
    input:
        tree = f"{PHYLO_DIR}/families/Iflaviridae/tree.newick",
        meta = f"{PHYLO_DIR}/families/Iflaviridae/taxa_metadata.tsv",
    output:
        svg = f"{MS_FIG_DIR}/fig5_iflaviridae.svg",
    params:
        genome_map = f"{GM_DIR}/Iflaviridae/genome_map_data.tsv",
        manifest   = f"{SNOWBALL_DIR}/manifest_combined.csv",
    log:
        "logs/figures/fig5_iflaviridae.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.svg})
        MAN_ARG=""
        [ -s "{params.manifest}" ] && MAN_ARG="--manifest {params.manifest}"
        GM_ARG=""
        [ -s "{params.genome_map}" ] && GM_ARG="--genome-map {params.genome_map}"
        python3 scripts/plot_phylo_manuscript.py \
            --mode fig5 \
            --tree "{input.tree}" \
            --meta "{input.meta}" \
            $GM_ARG \
            $MAN_ARG \
            --out-svg "{output.svg}" \
            > "{log}" 2>&1
        """


rule fig6_amarillo_negssrna:
    """Fig. 6: Segmented Jingmen-like Amarillovirales & negative-sense RNA viruses."""
    input:
        amar_tree = f"{PHYLO_DIR}/families/Amarillovirales/tree.newick",
        amar_meta = f"{PHYLO_DIR}/families/Amarillovirales/taxa_metadata.tsv",
        phen_tree = f"{PHYLO_DIR}/families/Phenuiviridae/tree.newick",
        phen_meta = f"{PHYLO_DIR}/families/Phenuiviridae/taxa_metadata.tsv",
        aliu_tree = f"{PHYLO_DIR}/families/Aliusviridae/tree.newick",
        aliu_meta = f"{PHYLO_DIR}/families/Aliusviridae/taxa_metadata.tsv",
    output:
        svg = f"{MS_FIG_DIR}/fig6_amarillo_negssrna.svg",
    params:
        manifest = f"{SNOWBALL_DIR}/manifest_combined.csv",
    log:
        "logs/figures/fig6_amarillo_negssrna.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.svg})
        MAN_ARG=""
        [ -s "{params.manifest}" ] && MAN_ARG="--manifest {params.manifest}"
        python3 scripts/plot_phylo_manuscript.py \
            --mode fig6 \
            --panel "Amarillovirales (Jingmen-like)|{input.amar_tree}|{input.amar_meta}" \
            --panel "Phenuiviridae (L RdRp)|{input.phen_tree}|{input.phen_meta}" \
            --panel "Aliusviridae (L RdRp)|{input.aliu_tree}|{input.aliu_meta}" \
            $MAN_ARG \
            --out-svg "{output.svg}" \
            > "{log}" 2>&1
        """
