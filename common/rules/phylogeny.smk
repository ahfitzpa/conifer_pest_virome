# Shared per-scenario phylogeny: select -> align -> trim -> tree.
#
# Every analysis pathway trees its own manifest's samples, so this file is
# included by each pathway's Snakefile rather than copied into it. Toolchain is
# FAMSA -> ClipKIT -> IQ-TREE2, from phylogenomics/scripts/run_phylogeny_workflow.py.
#
# The including Snakefile must define, before the include:
#   PHASE2          path to the phase2 directory (read-only)
#   INPUTS          resolved by common/rules/inputs.smk (contigs,
#                   contig_metadata, manifest, ...)
#   FAMILIES        list of families to tree
#   SCENARIO_PHYLO_DIR       where outputs go, e.g. "results/phylogeny"
#   COMMON          path to workflow/common
#   get_env(key)    resolves a tools: key to an env yaml

PHYLO = config.get("phylogeny", {})
# Inputs are declared, not assumed: INPUTS comes from common/rules/inputs.smk,
# which the including Snakefile resolves before this file. That is what lets a
# user point the pathway at their own contigs instead of this project's phase2.
MASTER_DB = INPUTS["contigs"]
MASTER_META = INPUTS["contig_metadata"]


def _ref_fasta(wildcards):
    # 1: optional per-family reference set, resolved against COMMON when the
    # config gives a relative path; absent is fine
    path = PHYLO.get("references", {}).get(wildcards.fam, "")
    if not path:
        return []
    if not os.path.isabs(path):
        path = os.path.join(COMMON, path)
    return [path] if os.path.exists(path) else []


rule select_contigs:
    # 2: pick this scenario's contigs for the family, driven by master_metadata
    # (consensus_family), not by the FASTA header's sample-level family
    input:
        fasta=MASTER_DB,
        metadata=MASTER_META,
        manifest=INPUTS["manifest"],
        keep=INPUTS["contig_keep_list"] or [],
        references=_ref_fasta,
    output:
        fasta=SCENARIO_PHYLO_DIR + "/{fam}/selected.fasta",
        table=SCENARIO_PHYLO_DIR + "/{fam}/selected.tsv",
    log:
        SCENARIO_PHYLO_DIR + "/{fam}/logs/select_contigs.log",
    params:
        script=os.path.join(COMMON, "scripts/select_scenario_contigs.py"),
        min_length=PHYLO.get("min_length", 1000),
        max_length=PHYLO.get("max_length", 0),
        agreement=PHYLO.get("agreement", "agree,single_tool"),
        max_per_sample=PHYLO.get("max_per_sample", 1),
        max_taxa=PHYLO.get("max_taxa", 0),
        min_taxa=PHYLO.get("min_taxa", 4),
        refs=lambda w, input: f"--references {input.references[0]}" if input.references else "",
        keep=lambda w, input: f"--keep-contigs {input.keep}" if input.keep else "",
    conda:
        get_env("COMMON_SEQTOOLS_ENV")
    threads: 1
    resources:
        mem_gb=8,
        runtime=60,
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --master-fasta {input.fasta} \
            --metadata {input.metadata} \
            --manifest {input.manifest} \
            --family {wildcards.fam} \
            --out-fasta {output.fasta} \
            --out-table {output.table} \
            --min-length {params.min_length} \
            --max-length {params.max_length} \
            --agreement '{params.agreement}' \
            --max-per-sample {params.max_per_sample} \
            --max-taxa {params.max_taxa} \
            --min-taxa {params.min_taxa} \
            {params.refs} {params.keep} > {log} 2>&1
        """


rule align_famsa:
    # 3: FAMSA; -medoidtree keeps it tractable as the scenarios grow to 329 runs
    input:
        SCENARIO_PHYLO_DIR + "/{fam}/selected.fasta",
    output:
        SCENARIO_PHYLO_DIR + "/{fam}/aln.fasta",
    log:
        SCENARIO_PHYLO_DIR + "/{fam}/logs/align_famsa.log",
    params:
        extra=PHYLO.get("famsa_extra", "-medoidtree -gt nj"),
    conda:
        get_env("SCENARIO_PHYLO_ENV")
    threads: 8
    resources:
        mem_gb=16,
        runtime=240,
    shell:
        "famsa -t {threads} {params.extra} {input} {output} > {log} 2>&1"


rule trim_clipkit:
    # 4: ClipKIT smart-gap, then refuse an alignment too short to tree
    input:
        SCENARIO_PHYLO_DIR + "/{fam}/aln.fasta",
    output:
        SCENARIO_PHYLO_DIR + "/{fam}/aln.trimmed.fasta",
    log:
        SCENARIO_PHYLO_DIR + "/{fam}/logs/trim_clipkit.log",
    params:
        mode=PHYLO.get("clipkit_mode", "smart-gap"),
        min_sites=PHYLO.get("min_sites", 300),
    conda:
        get_env("SCENARIO_PHYLO_ENV")
    threads: 1
    resources:
        mem_gb=8,
        runtime=120,
    shell:
        """
        clipkit {input} -m {params.mode} -o {output} > {log} 2>&1
        sites=$(awk '/^>/{{if (n++) exit; next}} {{w += length($0)}} END {{print w+0}}' {output})
        echo "trimmed width: $sites" >> {log}
        if [ "$sites" -lt {params.min_sites} ]; then
            echo "ERROR: {params.min_sites} sites required, got $sites" >> {log}
            exit 1
        fi
        """


rule tree_iqtree:
    # 5: IQ-TREE2 with model selection, UFBoot and SH-aLRT
    input:
        SCENARIO_PHYLO_DIR + "/{fam}/aln.trimmed.fasta",
    output:
        tree=SCENARIO_PHYLO_DIR + "/{fam}/tree.newick",
        iqtree=SCENARIO_PHYLO_DIR + "/{fam}/iqtree/{fam}.treefile",
    log:
        SCENARIO_PHYLO_DIR + "/{fam}/logs/tree_iqtree.log",
    params:
        prefix=lambda w: f"{SCENARIO_PHYLO_DIR}/{w.fam}/iqtree/{w.fam}",
        model=PHYLO.get("model", "MFP"),
        boot=PHYLO.get("ufboot", 1000),
        alrt=PHYLO.get("alrt", 1000),
        extra=PHYLO.get("iqtree_extra", ""),
    conda:
        get_env("SCENARIO_PHYLO_ENV")
    threads: 8
    resources:
        mem_gb=32,
        runtime=1440,
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        iqtree2 -s {input} -m {params.model} -B {params.boot} \
            -alrt {params.alrt} -T {threads} --prefix {params.prefix} \
            -redo {params.extra} > {log} 2>&1
        cp {output.iqtree} {output.tree}
        """


rule phylogeny_all:
    input:
        expand(SCENARIO_PHYLO_DIR + "/{fam}/tree.newick", fam=FAMILIES),
        expand(SCENARIO_PHYLO_DIR + "/{fam}/selected.tsv", fam=FAMILIES),
