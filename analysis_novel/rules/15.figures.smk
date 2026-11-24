# Overview and per-family figures.
#
# These scripts existed before this file and were run by hand (SLURM job
# "overview", 1182061/1182120 on 2026-11-23), so nothing regenerated a figure
# when its inputs changed and one of them silently skipped a panel:
# novelty_vs_completeness was dropped because results/demarcation was empty,
# and the log said so where nobody was reading. Declaring the inputs fixes
# both -- a figure is now rebuilt when the tables under it move.
#
# Generic plotters live in common/scripts (any pathway can use them); the
# conifer_* ones hardcode this project's three hosts and stay here.

FIG_DIR = "figures/overview"
FIG_EXT = ["png", "svg", "pdf"]


def _fig(stem, directory=FIG_DIR):
    return expand(directory + "/" + stem + ".{ext}", ext=FIG_EXT)


rule figure_host_family_matrix:
    # 1: host x family matrix and novelty-vs-completeness, across all families
    input:
        phylogeny = rules.phylogeny_all.input,
        demarcation = lambda wc: expand("results/demarcation/{fam}/pairwise_identity.tsv",
                                        fam=get_marker_family_list(wc))
    output:
        matrix = _fig("host_family_matrix"),
        novelty = _fig("novelty_vs_completeness")
    params:
        script = os.path.join(COMMON, "scripts/plot_multifamily_overview.py"),
        phylo_dir = SCENARIO_PHYLO_DIR,
        out_dir = FIG_DIR
    log:
        "logs/figures/host_family_matrix.log"
    benchmark:
        "benchmarks/figures/host_family_matrix.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("PHYLO_PLOT_ENV")
    shell:
        """
        mkdir -p $(dirname {log}) {params.out_dir}
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --phylogeny-dir {params.phylo_dir} \
            --demarcation-dir results/demarcation \
            --out-dir {params.out_dir} > {log} 2>&1
        """


rule figure_complex_heatmap:
    # 2: run x family contig-count heatmap, from phase2's master metadata
    input:
        metadata = MASTER_META,
        phylogeny = rules.phylogeny_all.input
    output:
        _fig("conifer_complex_heatmap")
    params:
        script = "scripts/plot_conifer_complex_heatmap.py",
        out_dir = FIG_DIR
    log:
        "logs/figures/complex_heatmap.log"
    benchmark:
        "benchmarks/figures/complex_heatmap.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("PHYLO_PLOT_ENV")
    shell:
        """
        mkdir -p $(dirname {log}) {params.out_dir}
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --master-metadata {input.metadata} \
            --out-dir {params.out_dir} > {log} 2>&1
        """


rule figure_virome_breadth:
    # 3: breadth panel, PCoA and the RdRp novelty scatter
    input:
        metadata = MASTER_META,
        phylogeny = rules.phylogeny_all.input,
        demarcation = lambda wc: expand("results/demarcation/{fam}/pairwise_identity.tsv",
                                        fam=get_marker_family_list(wc))
    output:
        breadth = _fig("conifer_virome_breadth"),
        pcoa = _fig("conifer_pcoa"),
        novelty = _fig("conifer_rdrp_novelty")
    params:
        script = "scripts/plot_conifer_virome_breadth.py",
        phylo_dir = SCENARIO_PHYLO_DIR,
        out_dir = FIG_DIR
    log:
        "logs/figures/virome_breadth.log"
    benchmark:
        "benchmarks/figures/virome_breadth.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("PHYLO_PLOT_ENV")
    shell:
        """
        mkdir -p $(dirname {log}) {params.out_dir}
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --phylogeny-dir {params.phylo_dir} \
            --demarcation-dir results/demarcation \
            --master-metadata {input.metadata} \
            --out-dir {params.out_dir} > {log} 2>&1
        """


rule figure_host_convergence:
    # 4: compositional convergence across hosts, from the signature tables
    input:
        signatures = rules.signatures_all.input
    output:
        _fig("conifer_host_convergence")
    params:
        script = "scripts/plot_conifer_convergence.py",
        phylo_dir = SCENARIO_PHYLO_DIR,
        out_dir = FIG_DIR
    log:
        "logs/figures/host_convergence.log"
    benchmark:
        "benchmarks/figures/host_convergence.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("PHYLO_PLOT_ENV")
    shell:
        """
        mkdir -p $(dirname {log}) {params.out_dir}
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --signatures-dir results/signatures \
            --phylogeny-dir {params.phylo_dir} \
            --fig-dir {params.out_dir} > {log} 2>&1
        """


rule figure_family_tree:
    # 5: per-family linear phylogram, annotated with host and CheckV quality
    input:
        tree = SCENARIO_PHYLO_DIR + "/{fam}/tree.newick",
        meta = SCENARIO_PHYLO_DIR + "/{fam}/selected.tsv"
    output:
        svg = "figures/{fam}/{fam}_tree_linear.svg",
        png = "figures/{fam}/{fam}_tree_linear.png"
    params:
        script = os.path.join(COMMON, "scripts/plot_phylo_tree.py"),
        title = lambda wc: f"{wc.fam} Maximum-Likelihood Phylogeny"
    log:
        "logs/figures/{fam}_tree.log"
    benchmark:
        "benchmarks/figures/{fam}_tree.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("PHYLO_PLOT_ENV")
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.svg})
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --tree {input.tree} --meta {input.meta} \
            --out-svg {output.svg} --out-png {output.png} \
            --title {params.title:q} > {log} 2>&1
        """


rule figures_all:
    input:
        rules.figure_host_family_matrix.output,
        rules.figure_complex_heatmap.output,
        rules.figure_virome_breadth.output,
        rules.figure_host_convergence.output,
        # the DNA families keep their tree, so this follows phylogeny_all's
        # family set rather than FAMILIES, which would ask for a skipped one
        expand("figures/{fam}/{fam}_tree_linear.png", fam=tree_families())
