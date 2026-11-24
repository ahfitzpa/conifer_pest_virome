# Codon alignments for dN/dS.
#
# The shared phylogeny aligns whole nucleotide contigs, which is right for a
# tree but wrong for selection analysis: codeml and HyPhy need codons in frame.
# This builds them from the same selected contigs.
#
# ORF1 (which carries the RdRp) is the target, per the project decision. It is
# identified by length: in Hepeviridae ORF1 spans roughly 5 kb of a 7.2 kb
# genome, so the longest Prodigal ORF is ORF1 in any genome assembled far
# enough to be worth analysing. Genomes whose longest ORF falls outside
# `orf_length_window` are dropped rather than analysed as if they were ORF1,
# and the dropped set is reported.

CODON_DIR = "results/codon"
CODON_CFG = config.get("codon", {})


rule predict_orfs:
    # 1: Prodigal in metagenomic mode; viral contigs are not single genomes
    input:
        fasta=SCENARIO_PHYLO_DIR + "/{fam}/selected.fasta",
    output:
        faa=CODON_DIR + "/{fam}/orfs.faa",
        ffn=CODON_DIR + "/{fam}/orfs.ffn",
        gff=CODON_DIR + "/{fam}/orfs.gff",
    log:
        CODON_DIR + "/{fam}/logs/predict_orfs.log",
    params:
        mode=CODON_CFG.get("prodigal_mode", "meta"),
    conda:
        get_env("ORF_TOOLS_ENV")
    threads: 1
    resources:
        mem_gb=8,
        runtime=120,
    shell:
        """
        prodigal -i {input.fasta} -a {output.faa} -d {output.ffn} \
            -o {output.gff} -f gff -p {params.mode} > {log} 2>&1
        """


rule select_orf1:
    # 2: longest ORF per genome, inside the expected ORF1 window
    input:
        faa=CODON_DIR + "/{fam}/orfs.faa",
        ffn=CODON_DIR + "/{fam}/orfs.ffn",
    output:
        faa=CODON_DIR + "/{fam}/orf1.faa",
        ffn=CODON_DIR + "/{fam}/orf1.ffn",
        table=CODON_DIR + "/{fam}/orf1_selection.tsv",
    log:
        CODON_DIR + "/{fam}/logs/select_orf1.log",
    params:
        script=os.path.join(COMMON, "scripts/select_longest_orf.py"),
        lo=CODON_CFG.get("orf_length_window", [1500, 8000])[0],
        hi=CODON_CFG.get("orf_length_window", [1500, 8000])[1],
        min_taxa=CODON_CFG.get("min_taxa", 4),
    conda:
        get_env("COMMON_SEQTOOLS_ENV")
    threads: 1
    resources:
        mem_gb=8,
        runtime=60,
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --faa {input.faa} --ffn {input.ffn} \
            --out-faa {output.faa} --out-ffn {output.ffn} \
            --out-table {output.table} \
            --min-nt {params.lo} --max-nt {params.hi} \
            --min-taxa {params.min_taxa} > {log} 2>&1
        """


rule align_orf1_protein:
    # 3: MAFFT on the proteins; the codon alignment is threaded onto this, so
    # the protein alignment is what defines the columns
    input:
        CODON_DIR + "/{fam}/orf1.faa",
    output:
        CODON_DIR + "/{fam}/orf1.aln.faa",
    log:
        CODON_DIR + "/{fam}/logs/align_orf1.log",
    params:
        extra=CODON_CFG.get("mafft_extra", "--auto --anysymbol"),
    conda:
        get_env("ORF_TOOLS_ENV")
    threads: 8
    resources:
        mem_gb=16,
        runtime=240,
    shell:
        "mafft {params.extra} --thread {threads} {input} > {output} 2> {log}"


rule codon_align_orf1:
    # 4: thread the CDS onto the protein alignment (pal2nal equivalent).
    # on-problem=drop, so a frameshifted or disagreeing sequence is reported and
    # excluded rather than silently repaired into a plausible frame.
    input:
        aln=CODON_DIR + "/{fam}/orf1.aln.faa",
        cds=CODON_DIR + "/{fam}/orf1.ffn",
    output:
        aln=CODON_DIR + "/{fam}/orf1.codon.fasta",
        report=CODON_DIR + "/{fam}/orf1.codon.qc.tsv",
    log:
        CODON_DIR + "/{fam}/logs/codon_align.log",
    params:
        script=os.path.join(COMMON, "scripts/codon_alignment.py"),
        code=CODON_CFG.get("genetic_code", 1),
        max_mismatch=CODON_CFG.get("max_mismatch_frac", 0.02),
        min_taxa=CODON_CFG.get("min_taxa", 4),
    conda:
        get_env("COMMON_SEQTOOLS_ENV")
    threads: 1
    resources:
        mem_gb=8,
        runtime=60,
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --aa-aln {input.aln} --cds {input.cds} \
            --out {output.aln} --report {output.report} \
            --genetic-code {params.code} \
            --max-mismatch-frac {params.max_mismatch} \
            --on-problem drop --min-taxa {params.min_taxa} > {log} 2>&1
        """


rule codon_all:
    input:
        expand(CODON_DIR + "/{fam}/orf1.codon.fasta", fam=rna_analysis_families()),
        expand(CODON_DIR + "/{fam}/orf1.codon.qc.tsv", fam=rna_analysis_families()),
        expand(CODON_DIR + "/{fam}/orf1_selection.tsv", fam=rna_analysis_families()),
