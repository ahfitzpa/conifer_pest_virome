# Phase 3: to-scale genome maps of candidate contigs (Zowalaty et al., Fig. 3).
#
# Reads Phase 2 outputs only; changes nothing about what is called viral.
# ORFs come from MetaCerberus product annotations, coordinates from the ORF
# table, contig lengths from the master FASTA, depth from Phase 2 coverage
# stats, and CheckV quality from the candidate table.

GM_DIR = "results/genome_maps"
GM_CFG = config["phylogenetics"]["genome_maps"]


rule genome_map_prepare:
    input:
        candidates = lambda wc: PCFG["candidates"],
        orf_annotation = lambda wc: PCFG["orf_annotation"],
        master = "results/master_viral_database.fasta",
        coverage = get_all_coverage_stats
    output:
        table = GM_DIR + "/{fam}/genome_map_data.tsv"
    params:
        config = PHYLO_CONFIG_FILE
    benchmark:
        "benchmarks/genome_map_prepare/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        python scripts/build_genome_map.py --prepare-only \
            --family {wildcards.fam} \
            --candidates {input.candidates} \
            --orf-annotation {input.orf_annotation} \
            --master-fasta {input.master} \
            --coverage {input.coverage} \
            --config {params.config} \
            --out-table {output.table}
        """


rule genome_map_plot:
    input:
        table = GM_DIR + "/{fam}/genome_map_data.tsv"
    output:
        pdf = "figures/genome_map_{fam}.pdf",
        svg = "figures/genome_map_{fam}.svg",
        png = "figures/genome_map_{fam}.png"
    params:
        config = PHYLO_CONFIG_FILE
    benchmark:
        "benchmarks/genome_map_plot/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("GENOME_MAPS_ENV")
    shell:
        """
        mkdir -p figures
        python scripts/build_genome_map.py --plot-only \
            --family {wildcards.fam} \
            --config {params.config} \
            --in-table {input.table} \
            --out-pdf {output.pdf} --out-svg {output.svg} --out-png {output.png}
        """


rule genome_map_all:
    input:
        tables = lambda wc: expand(GM_DIR + "/{fam}/genome_map_data.tsv",
                                   fam=get_family_list(wc)),
        figures = lambda wc: expand("figures/genome_map_{fam}.png",
                                    fam=get_family_list(wc))
