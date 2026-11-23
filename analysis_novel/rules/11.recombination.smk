# Phase 3: recombination screening of candidates with sufficient reference
# density (Zhao et al., Fig. 3 style).
#
# Two screeners, each in its own conda env: OpenRDP (numpy<2, as in
# CaliciTyper) and HyPhy GARD. OpenRDP is the RDP implementation used here --
# the old recomb_rdp5 rule shelled out to the standalone RDP5 binary, which is
# not on bioconda, was never supplied (rdp5_binary was always empty) and so
# never ran. OpenRDP's RDP scanner covers it.
# An alignment is only screened if >= min_references ICTV exemplars exist for
# that family AND it is >= min_alignment_nt long (checked in recomb_prepare).
#
# The screen runs on NUCLEOTIDES: the scenario phylogeny's trimmed alignment
# (FAMSA -> ClipKIT), not rule 09's amino-acid bundle. Every scanner here
# (GARD's ReadNucleotideDataSet, GENECONV, MaxChi, Chimaera, SiScan, 3Seq)
# assumes an aligned nucleotide matrix, and rule 09's combined.faa is neither
# aligned nor nucleotide -- which is why every eligible family died at parse
# with "Sequence N is longer than the previous sequences" on 2026-11-23.

REC_DIR = "results/recombination"
REC_CFG = config["phylogenetics"]["recombination"]
# Set by the pathway's Snakefile; the default keeps this file usable alone.
REC_ALN_DIR = globals().get("SCENARIO_PHYLO_DIR", "results/phylogeny")


rule recomb_prepare:
    input:
        aln = REC_ALN_DIR + "/{fam}/aln.trimmed.fasta",
        selected = REC_ALN_DIR + "/{fam}/selected.tsv"
    output:
        aln = REC_DIR + "/{fam}/recomb.aln.fasta",
        candidates = REC_DIR + "/{fam}/candidates.txt",
        eligible = REC_DIR + "/{fam}/eligible.tsv"
    params:
        min_refs = REC_CFG["min_references"],
        min_aln = REC_CFG["min_alignment_nt"],
        min_cov = REC_CFG.get("min_seq_coverage", 0.3)
    benchmark:
        "benchmarks/recomb_prepare/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/recomb_prepare.py \
            --family {wildcards.fam} \
            --alignment {input.aln} \
            --selected {input.selected} \
            --min-refs {params.min_refs} --min-aln {params.min_aln} \
            --min-coverage {params.min_cov} \
            --out-fasta {output.aln} --out-candidates {output.candidates} \
            --out-eligible {output.eligible}
        """


rule recomb_openrdp:
    input:
        aln = REC_DIR + "/{fam}/recomb.aln.fasta",
        candidates = REC_DIR + "/{fam}/candidates.txt",
        eligible = REC_DIR + "/{fam}/eligible.tsv"
    output:
        events = REC_DIR + "/{fam}/openrdp_events.tsv",
        log = REC_DIR + "/{fam}/openrdp.log"
    params:
        scanners = " ".join(REC_CFG["openrdp_scanners"]),
        min_methods = REC_CFG["openrdp_min_methods"],
        bootstrap = REC_CFG["bootstrap"]
    benchmark:
        "benchmarks/recomb_openrdp/{fam}.tsv"
    threads: 4
    resources: mem_gb=16, runtime=720
    conda: get_env("OPENRDP_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/run_openrdp.py \
            --alignment {input.aln} --eligible {input.eligible} \
            --candidates {input.candidates} \
            --scanners {params.scanners:q} \
            --min-methods {params.min_methods} \
            --bootstrap {params.bootstrap} \
            --out-tsv {output.events} 2> {output.log}
        """


rule recomb_hyphy:
    input:
        aln = REC_DIR + "/{fam}/recomb.aln.fasta",
        candidates = REC_DIR + "/{fam}/candidates.txt",
        eligible = REC_DIR + "/{fam}/eligible.tsv"
    output:
        events = REC_DIR + "/{fam}/hyphy_gard_events.tsv",
        log = REC_DIR + "/{fam}/hyphy_gard.log"
    params:
        method = REC_CFG["hyphy_method"],
        mode = REC_CFG["hyphy_mode"]
    benchmark:
        "benchmarks/recomb_hyphy/{fam}.tsv"
    threads: 4
    resources: mem_gb=16, runtime=1440
    conda: get_env("RECOMB_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/run_hyphy_gard.py \
            --alignment {input.aln} --eligible {input.eligible} \
            --method {params.method} --mode {params.mode} \
            --out-tsv {output.events} 2> {output.log}
        """


rule recomb_summary:
    input:
        events = expand(REC_DIR + "/{{fam}}/{m}_events.tsv",
                        m=["openrdp", "hyphy_gard"])
    output:
        summary = REC_DIR + "/{fam}/recombination_summary.tsv"
    params:
        min_methods = REC_CFG["openrdp_min_methods"]
    benchmark:
        "benchmarks/recomb_summary/{fam}.tsv"
    threads: 1
    resources: mem_gb=4, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/recomb_summary.py \
            --family {wildcards.fam} \
            --event-tables {input.events} \
            --min-methods {params.min_methods} \
            --out {output.summary}
        """


rule recomb_all:
    input:
        summaries = lambda wc: expand(REC_DIR + "/{fam}/recombination_summary.tsv",
                                      fam=get_family_list(wc))
