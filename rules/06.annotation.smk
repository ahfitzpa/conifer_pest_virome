# Functional annotation of the viral positives with MetaCerberus
# (https://github.com/raw-lab/MetaCerberus). See ENGINEERING_NOTES.md §
# run_metacerberus.
#
# Selection is deliberately loose -- every contig ANY taxonomy tool called
# viral (rules/05.merge_tax.smk), length-filtered only. CheckV/VirSorter2
# support stays an annotation column in master_metadata.csv rather than a
# gate here -- and, since 2026-11-17, not an upstream dependency either:
# annotation reads the per-sample viral_positives.fasta directly rather than
# the master database. See get_annotation_fastas below.

def get_annotation_fastas(wildcards):
    """The viral-positive FASTAs MetaCerberus actually annotates.

    Deliberately NOT results/master_viral_database.fasta. That path pulled in
    build_master_database, and through it every sample's taxonomy_merged.tsv,
    which requires CheckV (and, by submit_phase2.sh's stage ordering,
    VirSorter2). MetaCerberus calls ORFs on contigs and annotates them against
    HMMs: what it needs is the set of contigs the taxonomy tools called viral
    -- collect_viral_positives' output -- not their CheckV quality tier or
    VirSorter2 score. Depending on the merge meant one sample's CheckV failure
    withheld annotation for all of them.

    build_master_database still runs (rule all requests it) and still carries
    the CheckV/VirSorter2 columns into master_metadata.csv; it is simply no
    longer upstream of annotation.
    """
    want = [
        f"data/05.final_selection/{row.family}/{row.sample_id}/viral_positives.fasta"
        for _, row in MANIFEST.iterrows()
    ]
    # Same partial-aggregation switch build_master_database honours, so stage 3
    # of submit_phase2.sh annotates whatever completed instead of blocking.
    if not _allow_partial_master():
        return want
    return [f for f in want if os.path.exists(f)]


rule select_contigs_for_annotation:
    conda: get_env("TAX_CORE_ENV")
    input:
        fastas = get_annotation_fastas
    output:
        fa = "results/annotation/metacerberus_input.fasta"
    params:
        min_len = config["annotation"]["min_contig_len"]
    benchmark:
        "benchmarks/select_contigs_for_annotation/select.tsv"
    shell:
        # Prefix every contig with its sample before concatenating: SPAdes
        # NODE_* names repeat across samples, so without it IDs collide and
        # metacerberus_orfs.faa cannot be traced back to a library. Same
        # "{sample}|{contig}" convention the master database uses.
        """
        mkdir -p "$(dirname {output.fa})"
        TMP="{output.fa}.tmp"
        : > "$TMP"
        for f in {input.fastas}; do
            s=$(basename "$(dirname "$f")")
            seqkit replace -p '^' -r "${{s}}|" "$f" >> "$TMP"
        done
        seqkit seq -m {params.min_len} "$TMP" -o {output.fa}
        rm -f "$TMP"
        """

rule run_metacerberus:
    input:
        fa = "results/annotation/metacerberus_input.fasta",
        db = os.path.join(config["databases"]["metacerberus"], ".done")
    output:
        directory("results/annotation/metacerberus")
    benchmark:
        "benchmarks/run_metacerberus/run.tsv"
    conda: get_env("META_ENV")
    params:
        hmm = " ".join(config["annotation"]["metacerberus_hmm"])
    threads: 16
    # runtime: was 720 (12 h). Job 1164821 reached STEP 10 in 54 min on a
    # 13.8 MB input before dying on the Chrome/kaleido fault, so the work itself
    # is ~1 h here; 1440 is headroom for a larger manifest, and the job exits as
    # soon as it finishes rather than holding the allocation.
    # mem_gb is local bookkeeping only -- this cluster schedules on cores and
    # rejects --mem (see rules/03.assembly.smk).
    resources: mem_gb=64, runtime=1440
    shell:
        # --fraggenescan, not --prodigal: these are mostly eukaryotic RNA
        # viruses, and MetaCerberus documents Prodigal as working poorly
        # outside prokaryotes/phage.
        """
        # STEP 10 renders report figures with plotly -> kaleido, and kaleido >=1
        # needs an external Chrome (its bundled chromium is gone, and plotly 7
        # will not accept kaleido <1). choreographer reads BROWSER_PATH, so
        # pointing it at the shared download avoids an env rebuild. Without it
        # job 1164821 died at STEP 10 with ChromeNotFoundError after 54 min --
        # ORF calling and the HMMER search had finished -- and because `output`
        # is a directory, Snakemake deleted all of that completed work.
        export BROWSER_PATH="{config[chrome_binary]}"
        if [ ! -x "$BROWSER_PATH" ]; then
            echo "ERROR: chrome_binary not executable at $BROWSER_PATH -- MetaCerberus will fail at STEP 10 (report figures). Re-download with 'kaleido_get_chrome --path \\$(dirname \\$(dirname $BROWSER_PATH))'." >&2
            exit 1
        fi
        metacerberus.py \
            --fraggenescan {input.fa} \
            --hmm {params.hmm} \
            --db-path {config[databases][metacerberus]} \
            --dir-out {output} \
            --cpus {threads}
        """
