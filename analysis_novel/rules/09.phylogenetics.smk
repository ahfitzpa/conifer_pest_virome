# Phase 3: comparative phylogenetics for the candidate virus families.
#
# Depends on Phase 2 outputs from rules/07.reporting.smk (candidate_contigs.tsv,
# metacerberus_orfs.faa) and rules/08.ictv.smk (results/ictv/families.tsv and
# per-family bundles).  NOT in rule all of Phase 2 -- run as a separate phase
# (see SNAKEFILE_INTEGRATION.md).
#
# Three tasks:
#   1. phylo_*      ML tree per family: MAFFT -> trimAl -> IQ-TREE2 with
#                   ModelFinder, UFB and SH-aLRT.
#   2. concat_*     Concatenated conserved-domain phylogeny (Zowalaty Fig. 5).
#   3. diagnostic_* Amino-acid diagnostic residues + conservation annotation
#                   (Zhao Fig. 4).
#
# Families are discovered at runtime from ictv_prepare's index via a checkpoint,
# exactly as the scaffold in rules/08.ictv.smk anticipated.  Nothing here changes
# what Phase 2 calls viral; it only analyses the bundles it produced.

PHYLO_DIR = "results/phylo"
PCFG = config.get("phylogenetics", {})
PHYLO_CONFIG_FILE = "configs/phylogenetics.yaml"


# ---------------------------------------------------------------------------
# Runtime family discovery (checkpoint fan-out over results/ictv/families.tsv)
# ---------------------------------------------------------------------------
checkpoint phylo_family_index:
    input:
        families = ICTV_DIR + "/families.tsv"
    output:
        index = PHYLO_DIR + "/families.tsv"
    benchmark:
        "benchmarks/phylo_family_index/index.tsv"
    threads: 1
    resources: mem_gb=4, runtime=30
    shell:
        """
        mkdir -p $(dirname {output.index})
        cp {input.families} {output.index}
        """


def get_family_list(wildcards):
    checkpoints.phylo_family_index.get()
    families = []
    with open(PHYLO_DIR + "/families.tsv") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fam = line.split("\t")[0].strip()
            if fam in ("assigned_family", "family", ""):
                continue
            families.append(fam)
    wanted = PCFG.get("families", [])
    if wanted:
        families = [f for f in families if f in wanted]
    return families


def get_marker_family_list(wildcards):
    """Families this pathway analyses that actually have an RdRp/L marker.

    The Snakefile sets MARKER_FAMILIES; without it (a pathway including this
    file on its own) every analysed family is used, as before.
    """
    families = get_family_list(wildcards)
    marker = set(globals().get("MARKER_FAMILIES") or families)
    return [f for f in families if f in marker]


def get_concat_families(wildcards):
    domains = PCFG.get("domains", {})
    return [f for f in get_family_list(wildcards) if f in domains]


# ---------------------------------------------------------------------------
# Task 1: single-marker ML phylogeny (RdRp / L domain)
# ---------------------------------------------------------------------------
rule phylo_prepare:
    input:
        rdrp = ICTV_DIR + "/families/{fam}/rdrp_domains.faa",
        ictv_faa = lambda wc: PCFG["ictv_proteins"],
        ictv_meta = lambda wc: PCFG.get("ictv_metadata", ""),
        candidates = lambda wc: PCFG.get("candidates", REPORT_DIR + "/candidate_contigs.tsv")
    output:
        queries = PHYLO_DIR + "/families/{fam}/queries.faa",
        refs = PHYLO_DIR + "/families/{fam}/refs.faa",
        combined = PHYLO_DIR + "/families/{fam}/combined.faa",
        meta = PHYLO_DIR + "/families/{fam}/taxa_metadata.tsv"
    benchmark:
        "benchmarks/phylo_prepare/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/phylo_prepare.py \
            --family {wildcards.fam} \
            --rdrp-faa {input.rdrp} \
            --ictv-faa {input.ictv_faa} \
            --ictv-meta {input.ictv_meta:q} \
            --candidates {input.candidates} \
            --out-queries {output.queries} \
            --out-refs {output.refs} \
            --out-combined {output.combined} \
            --out-metadata {output.meta}
        """


rule phylo_align:
    input:
        combined = PHYLO_DIR + "/families/{fam}/combined.faa"
    output:
        aln = PHYLO_DIR + "/families/{fam}/aln.faa",
        trimmed = PHYLO_DIR + "/families/{fam}/aln.trimmed.faa"
    params:
        strategy = PCFG.get("mafft_strategy", "--auto"),
        trim = PCFG.get("trimal_strategy", "-autom1")
    benchmark:
        "benchmarks/phylo_align/{fam}.tsv"
    threads: 4
    resources: mem_gb=16, runtime=240
    conda: get_env("PHYLO_ENV")
    shell:
        """
        mafft {params.strategy} {input.combined} > {output.aln} 2> {output.aln}.log
        trimal {params.trim} -in {output.aln} -out {output.trimmed} > {output.trimmed}.log 2>&1
        """


# Task 1 ends at the alignment. The ML tree that used to be built here
# (rule phylo_tree -> tree.newick, and phylo_summary over it) was consumed by
# nothing but its own summary table: the tree every downstream step actually
# reads is the scenario phylogeny's nucleotide tree in SCENARIO_PHYLO_DIR,
# which rules 11, 13 and 14 and the figure scripts all use. Two ML trees per
# family were being built for one consumer; aln.trimmed.faa is what rule 12
# (demarcation) needs from this task, so the alignment stays and the second
# IQ-TREE run is gone.


# ---------------------------------------------------------------------------
# Task 2: concatenated conserved-domain phylogeny (Zowalaty et al., Fig. 5)
# ---------------------------------------------------------------------------
rule concat_prepare:
    input:
        orfs = lambda wc: PCFG["orf_proteins"],
        orf_annotation = lambda wc: PCFG["orf_annotation"],
        ictv_faa = lambda wc: PCFG["ictv_proteins"],
        ictv_meta = lambda wc: PCFG.get("ictv_metadata", ""),
        candidates = lambda wc: PCFG.get("candidates", REPORT_DIR + "/candidate_contigs.tsv"),
        rdrp = ICTV_DIR + "/families/{fam}/rdrp_domains.faa"
    output:
        concat = PHYLO_DIR + "/families/{fam}/concat_domains.faa",
        matrix = PHYLO_DIR + "/families/{fam}/domain_presence.tsv",
        schematic = PHYLO_DIR + "/families/{fam}/domain_schematic.tsv"
    params:
        config = PHYLO_CONFIG_FILE
    benchmark:
        "benchmarks/concat_prepare/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/concat_domains.py \
            --family {wildcards.fam} \
            --orf-faa {input.orfs} \
            --orf-annotation {input.orf_annotation} \
            --ictv-faa {input.ictv_faa} \
            --ictv-meta {input.ictv_meta:q} \
            --candidates {input.candidates} \
            --rdrp-faa {input.rdrp} \
            --config {params.config} \
            --out-faa {output.concat} \
            --out-matrix {output.matrix} \
            --out-schematic {output.schematic}
        """


rule concat_align:
    input:
        concat = PHYLO_DIR + "/families/{fam}/concat_domains.faa"
    output:
        aln = PHYLO_DIR + "/families/{fam}/concat.aln.faa",
        trimmed = PHYLO_DIR + "/families/{fam}/concat.aln.trimmed.faa"
    params:
        strategy = PCFG.get("mafft_strategy", "--auto"),
        trim = PCFG.get("trimal_strategy", "-autom1")
    benchmark:
        "benchmarks/concat_align/{fam}.tsv"
    threads: 4
    resources: mem_gb=16, runtime=240
    conda: get_env("PHYLO_ENV")
    shell:
        """
        mafft {params.strategy} {input.concat} > {output.aln} 2> {output.aln}.log
        trimal {params.trim} -in {output.aln} -out {output.trimmed} > {output.trimmed}.log 2>&1
        """


rule concat_tree:
    input:
        aln = PHYLO_DIR + "/families/{fam}/concat.aln.trimmed.faa"
    output:
        tree = PHYLO_DIR + "/families/{fam}/concat.tree.newick",
        iqtree = PHYLO_DIR + "/families/{fam}/concat.iqtree.txt",
        log = PHYLO_DIR + "/families/{fam}/concat.iqtree.log"
    params:
        prefix = PHYLO_DIR + "/families/{fam}/concat_iqtree",
        model = PCFG.get("iqtree_model", "MFP"),
        ufb = PCFG.get("ufb_replicates", 1000),
        alrt = PCFG.get("alrt_replicates", 1000)
    benchmark:
        "benchmarks/concat_tree/{fam}.tsv"
    threads: 8
    resources: mem_gb=32, runtime=2880
    conda: get_env("PHYLO_ENV")
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        iqtree2 -s {input.aln} -m {params.model} \
            -B {params.ufb} -alrt {params.alrt} -T {threads} \
            --prefix {params.prefix} > {output.log} 2>&1
        cp {params.prefix}.treefile {output.tree}
        cp {params.prefix}.iqtree {output.iqtree}
        """


# ---------------------------------------------------------------------------
# Task 3: amino-acid diagnostic residues + conservation (Zhao et al., Fig. 4)
# ---------------------------------------------------------------------------
rule diagnostic_table:
    input:
        aln = PHYLO_DIR + "/families/{fam}/aln.faa",
        meta = PHYLO_DIR + "/families/{fam}/taxa_metadata.tsv"
    output:
        matrix = PHYLO_DIR + "/families/{fam}/diagnostic_residues.tsv",
        conservation = PHYLO_DIR + "/families/{fam}/conservation.tsv"
    params:
        config = PHYLO_CONFIG_FILE
    benchmark:
        "benchmarks/diagnostic_table/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/diagnostic_residues.py \
            --alignment {input.aln} \
            --metadata {input.meta} \
            --config {params.config} \
            --out-matrix {output.matrix} \
            --out-conservation {output.conservation}
        """


rule diagnostic_plot:
    input:
        matrix = PHYLO_DIR + "/families/{fam}/diagnostic_residues.tsv",
        conservation = PHYLO_DIR + "/families/{fam}/conservation.tsv",
        aln = PHYLO_DIR + "/families/{fam}/aln.faa",
        meta = PHYLO_DIR + "/families/{fam}/taxa_metadata.tsv"
    output:
        pdf = "figures/{fam}/diagnostic_{fam}.pdf",
        svg = "figures/{fam}/diagnostic_{fam}.svg",
        png = "figures/{fam}/diagnostic_{fam}.png"
    params:
        config = PHYLO_CONFIG_FILE
    benchmark:
        "benchmarks/diagnostic_plot/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("PHYLO_PLOT_ENV")
    shell:
        """
        mkdir -p figures/{wildcards.fam}
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/plot_diagnostic_matrix.py \
            --matrix {input.matrix} \
            --conservation {input.conservation} \
            --alignment {input.aln} \
            --metadata {input.meta} \
            --family {wildcards.fam} \
            --config {params.config} \
            --out-pdf {output.pdf} --out-svg {output.svg} --out-png {output.png}
        """


# ---------------------------------------------------------------------------
# Phase 3 target: everything, fanned out over discovered families.
# ---------------------------------------------------------------------------
rule phylo_all:
    input:
        alignments = lambda wc: expand(PHYLO_DIR + "/families/{fam}/aln.trimmed.faa",
                                       fam=get_marker_family_list(wc)),
        concat_trees = lambda wc: expand(PHYLO_DIR + "/families/{fam}/concat.tree.newick",
                                         fam=get_concat_families(wc)),
        matrices = lambda wc: expand(PHYLO_DIR + "/families/{fam}/diagnostic_residues.tsv",
                                     fam=get_marker_family_list(wc)),
        plots = lambda wc: expand("figures/{fam}/diagnostic_{fam}.png",
                                  fam=get_marker_family_list(wc))
