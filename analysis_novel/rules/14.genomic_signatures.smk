# Genomic signatures and adaptive evolution (draft 26, split into steps).
#
# Site-level selection from HyPhy on the ORF1 codon alignment, plus compositional
# signatures (RSCU, dinucleotide odds) per genome. One script per step, each
# writing one table, so a failing step cannot take the rest down and every term
# can be recomputed on its own.
#
# HyPhy needs a tree over exactly the taxa in the alignment. The shared
# phylogeny's tree includes reference genomes and any contig the codon QC
# dropped, so it is pruned first rather than passed through.

SIG_DIR = "results/signatures"
SIG_CFG = config.get("signatures", {})


rule prune_tree_to_alignment:
    # 1: restrict the tree to the taxa the codon alignment actually kept, and
    # strip UFBoot/SH-aLRT labels, which HyPhy and PAML will not parse
    input:
        tree=SCENARIO_PHYLO_DIR + "/{fam}/tree.newick",
        aln=CODON_DIR + "/{fam}/orf1.codon.fasta",
    output:
        tree=SIG_DIR + "/{fam}/orf1.pruned.newick",
        report=SIG_DIR + "/{fam}/orf1.pruned.tsv",
    log:
        SIG_DIR + "/{fam}/logs/prune_tree.log",
    params:
        script=os.path.join(COMMON, "scripts/prune_tree.py"),
        min_taxa=SIG_CFG.get("min_taxa", 4),
    conda:
        get_env("COMMON_SEQTOOLS_ENV")
    threads: 1
    resources:
        mem_gb=4,
        runtime=30,
    shell:
        """
        python3 -I {params.script} \
            --tree {input.tree} --alignment {input.aln} \
            --out-tree {output.tree} --out-report {output.report} \
            --strip-support --min-taxa {params.min_taxa} > {log} 2>&1
        """


rule hyphy_fubar:
    # 2: FUBAR, pervasive selection per site
    input:
        aln=CODON_DIR + "/{fam}/orf1.codon.fasta",
        tree=SIG_DIR + "/{fam}/orf1.pruned.newick",
    output:
        json=SIG_DIR + "/{fam}/fubar.json",
    log:
        SIG_DIR + "/{fam}/logs/fubar.log",
    params:
        code=SIG_CFG.get("hyphy_genetic_code", "Universal"),
    conda:
        get_env("GENOMIC_SIGNATURES_ENV")
    threads: 4
    resources:
        mem_gb=16,
        runtime=720,
    shell:
        """
        hyphy CPU={threads} fubar --alignment {input.aln} --tree {input.tree} \
            --code {params.code} --output {output.json} > {log} 2>&1
        """


rule hyphy_meme:
    # 3: MEME, episodic selection per site. Slower than FUBAR and allowed to
    # fail on its own without losing the FUBAR result.
    input:
        aln=CODON_DIR + "/{fam}/orf1.codon.fasta",
        tree=SIG_DIR + "/{fam}/orf1.pruned.newick",
    output:
        json=SIG_DIR + "/{fam}/meme.json",
    log:
        SIG_DIR + "/{fam}/logs/meme.log",
    params:
        code=SIG_CFG.get("hyphy_genetic_code", "Universal"),
    conda:
        get_env("GENOMIC_SIGNATURES_ENV")
    threads: 8
    resources:
        mem_gb=32,
        runtime=1440,
    shell:
        """
        hyphy CPU={threads} meme --alignment {input.aln} --tree {input.tree} \
            --code {params.code} --output {output.json} > {log} 2>&1
        """


rule selection_summary:
    # 4: apply the documented thresholds and categorise each site
    input:
        fubar=SIG_DIR + "/{fam}/fubar.json",
        meme=SIG_DIR + "/{fam}/meme.json",
    output:
        sites=SIG_DIR + "/{fam}/site_selection.tsv",
        summary=SIG_DIR + "/{fam}/selection_summary.tsv",
    log:
        SIG_DIR + "/{fam}/logs/selection_summary.log",
    params:
        script=os.path.join(COMMON, "scripts/selection_summary.py"),
        pos=SIG_CFG.get("omega_positive", 1.5),
        pur=SIG_CFG.get("omega_purifying", 0.5),
        post=SIG_CFG.get("fubar_posterior", 0.9),
        pval=SIG_CFG.get("meme_pvalue", 0.05),
    conda:
        get_env("COMMON_SEQTOOLS_ENV")
    threads: 1
    resources:
        mem_gb=8,
        runtime=60,
    shell:
        """
        python3 -I {params.script} \
            --fubar {input.fubar} --meme {input.meme} \
            --out-sites {output.sites} --out-summary {output.summary} \
            --omega-positive {params.pos} --omega-purifying {params.pur} \
            --fubar-posterior {params.post} --meme-pvalue {params.pval} \
            > {log} 2>&1
        """


rule composition_signatures:
    # 5: RSCU and dinucleotide odds ratios per genome, from the ORF1 CDS
    input:
        cds=CODON_DIR + "/{fam}/orf1.ffn",
    output:
        rscu=SIG_DIR + "/{fam}/rscu.tsv",
        dinuc=SIG_DIR + "/{fam}/dinucleotide_odds.tsv",
    log:
        SIG_DIR + "/{fam}/logs/composition.log",
    params:
        script=os.path.join(COMMON, "scripts/composition_signatures.py"),
        code=SIG_CFG.get("genetic_code", 1),
    conda:
        get_env("COMMON_SEQTOOLS_ENV")
    threads: 1
    resources:
        mem_gb=8,
        runtime=60,
    shell:
        """
        python3 -I {params.script} --cds {input.cds} \
            --out-rscu {output.rscu} --out-dinuc {output.dinuc} \
            --genetic-code {params.code} > {log} 2>&1
        """


rule signatures_all:
    input:
        expand(SIG_DIR + "/{fam}/site_selection.tsv", fam=FAMILIES),
        expand(SIG_DIR + "/{fam}/selection_summary.tsv", fam=FAMILIES),
        expand(SIG_DIR + "/{fam}/rscu.tsv", fam=FAMILIES),
        expand(SIG_DIR + "/{fam}/dinucleotide_odds.tsv", fam=FAMILIES),
