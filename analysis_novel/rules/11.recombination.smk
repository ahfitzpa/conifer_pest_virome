# Phase 3: recombination screening of candidates with sufficient reference
# density (Zhao et al., Fig. 3 style).
#
# OpenRDP lives in its own conda env (numpy<2, as in CaliciTyper); HyPhy GARD
# in another. RDP5 is a standalone binary supplied via config -- if it is not
# present the rule FAILS LOUDLY rather than writing "not_assessed".
# An alignment is only screened if >= min_references ICTV exemplars exist for
# that family AND it is >= min_alignment_nt long (checked in recomb_prepare).

REC_DIR = "results/recombination"
REC_CFG = config["phylogenetics"]["recombination"]


rule recomb_prepare:
    input:
        combined = PHYLO_DIR + "/families/{fam}/combined.faa",
        metadata = PHYLO_DIR + "/families/{fam}/taxa_metadata.tsv"
    output:
        aln = REC_DIR + "/{fam}/recomb.aln.faa",
        eligible = REC_DIR + "/{fam}/eligible.tsv"
    params:
        min_refs = REC_CFG["min_references"],
        min_aln = REC_CFG["min_alignment_nt"]
    benchmark:
        "benchmarks/recomb_prepare/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        python scripts/recomb_prepare.py \
            --family {wildcards.fam} \
            --sequences {input.combined} \
            --metadata {input.metadata} \
            --min-refs {params.min_refs} --min-aln {params.min_aln} \
            --out-faa {output.aln} --out-eligible {output.eligible}
        """


rule recomb_openrdp:
    input:
        aln = REC_DIR + "/{fam}/recomb.aln.faa",
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
        python scripts/run_openrdp.py \
            --alignment {input.aln} --eligible {input.eligible} \
            --scanners {params.scanners:q} \
            --min-methods {params.min_methods} \
            --bootstrap {params.bootstrap} \
            --out-tsv {output.events} 2> {output.log}
        """


rule recomb_hyphy:
    input:
        aln = REC_DIR + "/{fam}/recomb.aln.faa",
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
        python scripts/run_hyphy_gard.py \
            --alignment {input.aln} --eligible {input.eligible} \
            --method {params.method} --mode {params.mode} \
            --out-tsv {output.events} 2> {output.log}
        """


rule recomb_rdp5:
    input:
        aln = REC_DIR + "/{fam}/recomb.aln.faa",
        eligible = REC_DIR + "/{fam}/eligible.tsv"
    output:
        events = REC_DIR + "/{fam}/rdp5_events.tsv",
        log = REC_DIR + "/{fam}/rdp5.log"
    params:
        binary = REC_CFG.get("rdp5_binary", "")
    benchmark:
        "benchmarks/recomb_rdp5/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=720
    shell:
        """
        python scripts/run_rdp5.py \
            --alignment {input.aln} --eligible {input.eligible} \
            --rdp5-binary {params.binary:q} \
            --out-tsv {output.events} 2> {output.log}
        """


rule recomb_summary:
    input:
        events = expand(REC_DIR + "/{{fam}}/{m}_events.tsv",
                        m=["openrdp", "hyphy_gard"] +
                          (["rdp5"] if REC_CFG.get("rdp5_binary") else []))
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
        python scripts/recomb_summary.py \
            --family {wildcards.fam} \
            --event-tables {input.events} \
            --min-methods {params.min_methods} \
            --out {output.summary}
        """


rule recomb_all:
    input:
        summaries = lambda wc: expand(REC_DIR + "/{fam}/recombination_summary.tsv",
                                      fam=get_family_list(wc))
