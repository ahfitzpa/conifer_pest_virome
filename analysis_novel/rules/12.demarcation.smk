# Phase 3: formal ICTV species/genus demarcation tables.
#
# Every threshold comes from configs/phylogenetics.yaml:demarcation.thresholds.
# Families whose criterion is not configured produce an explicit
# `criterion_not_configured` status, never an invented number. When the
# criterion marker (e.g. capsid) differs from the alignment marker (RdRp/L),
# the call is withheld and `criterion_marker_mismatch` is reported.

DM_DIR = "results/demarcation"
DM_CFG = config["phylogenetics"]["demarcation"]


rule demarcation_distances:
    input:
        aln = PHYLO_DIR + "/families/{fam}/aln.trimmed.faa",
        meta = PHYLO_DIR + "/families/{fam}/taxa_metadata.tsv"
    output:
        pairwise = DM_DIR + "/{fam}/pairwise_identity.tsv"
    benchmark:
        "benchmarks/demarcation_distances/{fam}.tsv"
    threads: 4
    resources: mem_gb=16, runtime=240
    conda: get_env("DEMARCATION_ENV")
    shell:
        """
        python scripts/demarcation_distances.py \
            --alignment {input.aln} --metadata {input.meta} \
            --out {output.pairwise}
        """


rule demarcation_table:
    input:
        pairwise = DM_DIR + "/{fam}/pairwise_identity.tsv",
        meta = PHYLO_DIR + "/families/{fam}/taxa_metadata.tsv"
    output:
        table = DM_DIR + "/{fam}/demarcation.tsv"
    params:
        config = PHYLO_CONFIG_FILE
    benchmark:
        "benchmarks/demarcation_table/{fam}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("DEMARCATION_ENV")
    shell:
        """
        python scripts/ictv_demarcation.py \
            --family {wildcards.fam} \
            --pairwise {input.pairwise} --metadata {input.meta} \
            --config {params.config} --out {output.table}
        """


rule demarcation_all:
    input:
        tables = lambda wc: expand(DM_DIR + "/{fam}/demarcation.tsv",
                                   fam=get_family_list(wc))
