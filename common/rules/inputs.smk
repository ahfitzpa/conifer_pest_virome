# Declared inputs for an analysis pathway.
#
# Each pathway states what data it needs, rather than reaching into this
# project's phase2 directory, so anyone can run it on their own assemblies:
#
#   snakemake -s Snakefile_novel --config \
#       inputs_manifest=my_samples.csv \
#       inputs_contigs=my_contigs.fasta \
#       inputs_contig_metadata=my_contig_table.tsv \
#       -- novel_all
#
# Resolution order, first hit wins:
#   1 --config inputs_<name>=...        (command line, highest)
#   2 ANALYSIS_<NAME> environment variable
#   3 the `inputs:` block of the pathway's config yaml
#   4 the phase2 default, for users running the full project
#
# Required inputs are validated at parse, so a missing one stops the workflow
# with a message naming the input and how to supply it -- not with a
# MissingInputException several rules deep.

import os

INPUT_SPEC = {
    # name:            (required, phase2 default, what it is)
    "manifest":        (True,  "data/test_run/manifest.csv",
                        "sample table; needs sample_id (or run_id) and family, "
                        "and scenario if the pathway routes by it"),
    "contigs":         (True,  "results/master_viral_database.fasta",
                        "nucleotide FASTA of candidate contigs"),
    "contig_metadata": (True,  "results/taxonomy/master_metadata.csv",
                        "per-contig table with run_id, contig_id, "
                        "consensus_family, family_agreement; drives selection"),
    "biosample_metadata": (False, "", "host/geography per BioSample accession"),
    "run_map":         (False, "", "run_id -> biosample, scientificname"),
    "reference_table": (False, "", "reference label -> host, country, year"),
    "contig_keep_list": (False, "results/reporting/contigs_keep.txt",
                        "phase2 contig_confidence keep list; restricts the "
                        "analysis to contigs the annotation supports"),
    "tree":            (False, "", "pre-built tree; skips the phylogeny stage"),
    "codon_alignment": (False, "", "pre-built codon alignment; skips rule 13"),
}


def _env_name(name):
    return "ANALYSIS_" + name.upper()


def resolve_inputs(config, phase2, pathway="analysis", require=None,
                   optional_required=()):
    """Resolve every declared input and validate the required ones.

    `require` overrides which inputs are mandatory for this pathway;
    `optional_required` adds pathway-specific ones (e.g. biosample_metadata for
    host range, which cannot assign host states without it).
    """
    block = (config.get("inputs") or {})
    resolved, sources, problems = {}, {}, []

    for name, (required, default, description) in INPUT_SPEC.items():
        is_required = (name in require) if require is not None else required
        if name in optional_required:
            is_required = True
        cli = config.get(f"inputs_{name}")
        env = os.environ.get(_env_name(name))
        cfg = block.get(name)
        if cli:
            value, source = str(cli), "--config"
        elif env:
            value, source = env, _env_name(name)
        elif cfg:
            value, source = str(cfg), "config yaml"
        elif default:
            # the project's own config knows where its manifest lives
            # (conifer: data/intake/manifest.csv, Serratus: data/test_run/...)
            if name == "manifest" and config.get("manifest"):
                default = config["manifest"]
            value, source = os.path.join(phase2, default), "phase2 default"
        else:
            value, source = "", ""

        explicit = source in ("--config", _env_name(name), "config yaml")
        if value and not os.path.exists(value):
            if not is_required and not explicit:
                # an optional input whose phase2 default has not been produced
                # yet is simply absent -- phase2 writes contigs_keep.txt at the
                # end of its run, and the analysis must still parse before then
                value, source = "", ""
            else:
                problems.append(
                    f"  {name}: '{value}' does not exist (from {source})\n"
                    f"      {description}")
        elif is_required and not value:
            problems.append(
                f"  {name}: required but not set\n"
                f"      {description}\n"
                f"      supply with --config inputs_{name}=PATH, "
                f"{_env_name(name)}=PATH, or inputs: {name}: in the config")
        resolved[name] = value
        sources[name] = source

    if problems:
        raise WorkflowError(
            f"{pathway}: declared inputs are not usable.\n"
            + "\n".join(problems)
            + "\n\nRun with --config inputs_manifest=... inputs_contigs=... "
              "inputs_contig_metadata=... to use your own data instead of this "
              "project's phase2 outputs.")

    for name in sorted(resolved):
        if resolved[name]:
            print(f"[inputs] {name:19} {resolved[name]}  ({sources[name]})")
    return resolved
