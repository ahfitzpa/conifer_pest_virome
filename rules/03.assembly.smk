import os

def get_assembly_input(wildcards):
    return {
        "f1": f"data/02.fastp/{wildcards.family}/{wildcards.sample}/{wildcards.sample}_cleaned_1.fastq.gz",
        "f2": f"data/02.fastp/{wildcards.family}/{wildcards.sample}/{wildcards.sample}_cleaned_2.fastq.gz"
    }

# Highmem routing for assemblies that OOM on the shared partition, and the
# known-OOM exception list. Ported from the Hepeviridae run (2026-11-12), where
# seven deep metagenomic libraries failed three times each with "Cannot
# allocate memory" in spades-core.
HIGHMEM_ATTEMPT = 3


def _highmem_samples():
    p = "configs/highmem_assembly_samples.txt"
    if not os.path.exists(p):
        return set()
    with open(p) as fh:
        return {ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")}


HIGHMEM_SAMPLES = _highmem_samples()


def _assembly_highmem(wildcards, attempt):
    return attempt >= HIGHMEM_ATTEMPT or wildcards.sample in HIGHMEM_SAMPLES


rule rnaviral_assembly:
    input:
        unpack(get_assembly_input)
    output:
        # temp(): consumed by filter_contigs AND contig_coverage_stats. Both
        # are declared consumers, so Snakemake holds the BAM until the later
        # of the two finishes. Coverage stats must be taken before the BAM is
        # reclaimed -- nothing downstream can recover it.
        contigs = temp("data/03.rnaviralSPAdes/{family}/{sample}/contigs.fasta"),
        bam = temp("data/03.rnaviralSPAdes/{family}/{sample}/sorted.bam")
    benchmark:
        "benchmarks/rnaviral_assembly/{family}_{sample}.tsv"
    threads: 16
    # No mem_mb: this cluster's `shared` partition schedules on cores
    # only and rejects any --mem request outright; mem_gb below is
    # local-only bookkeeping for the free-memory probe in the shell block.
    resources:
        mem_gb=lambda wildcards, attempt: 256 if _assembly_highmem(wildcards, attempt) else 48,
        runtime=lambda wildcards, attempt: 960 if _assembly_highmem(wildcards, attempt) else 480,
        slurm_extra=lambda wildcards, attempt: ("--constraint=highmem"
                                                if _assembly_highmem(wildcards, attempt) else ""),
        # attempt_no is a resource, and the exception list is read in the shell,
        # because a params: block here would re-run every assembly.
        attempt_no=lambda wildcards, attempt: attempt
    retries: 2
    # Not ASSEMBLY_ENV -- see ENGINEERING_NOTES.md § rnaviral_assembly:
    # bwa-mem3 and hostile=1.1.0's pinned bowtie2 can't coexist in one
    # conda env (zlib version conflict).
    conda: get_env("RNAVIRAL_ASSEMBLY_ENV")
    shell:
        """
        # A listed sample starts on highmem, so a retry would repeat it.
        HIGHMEM_LIST=configs/highmem_assembly_samples.txt
        if [ "{resources.attempt_no}" -gt 1 ] && [ -f "$HIGHMEM_LIST" ] \
           && grep -qx "{wildcards.sample}" "$HIGHMEM_LIST"; then
            echo "ERROR: {wildcards.sample} is in $HIGHMEM_LIST and its highmem attempt already failed." >&2
            echo "Not repeating an identical highmem run. Reduce its input instead." >&2
            exit 1
        fi

        OUT_DIR=$(dirname "{output.contigs}")
        # Node-local disk, not BeeGFS -- see ENGINEERING_NOTES.md §
        # rnaviral_assembly (speed-up plan item 10). SPAdes' own
        # --tmp-dir only covers error-correction files, not the
        # per-k-mer scaffolding graphs that dominate its small-write
        # volume, so -o itself has to move.
        LOCAL_DIR="/tmp/rnaviral_assembly_{wildcards.family}_{wildcards.sample}_${{SLURM_JOB_ID:-$$}}"
        mkdir -p "$LOCAL_DIR"

        # Detect real free memory on whatever node this attempt lands on,
        # rather than trusting the static resources.mem_gb fallback blind.
        AVAIL_MEM_GB=$(( $(free -m | awk '/^Mem:/{{print $7}}') / 1024 ))
        if [ -z "$AVAIL_MEM_GB" ] || [ "$AVAIL_MEM_GB" -le 0 ]; then
            AVAIL_MEM_GB={resources.mem_gb}
        fi
        SPADES_MEM_GB=$(( AVAIL_MEM_GB / 2 ))
        if [ "$SPADES_MEM_GB" -lt 16 ]; then
            SPADES_MEM_GB=16
        fi
        if [ "$SPADES_MEM_GB" -gt {resources.mem_gb} ]; then
            SPADES_MEM_GB={resources.mem_gb}
        fi

        rnaviralspades.py --pe1-1 {input.f1} --pe1-2 {input.f2} \
            -k 21,33,45,55,67,77 -o "$LOCAL_DIR" -t {threads} -m $SPADES_MEM_GB
        cp "$LOCAL_DIR/contigs.fasta" {output.contigs}
        rm -rf "$LOCAL_DIR"

        # See ENGINEERING_NOTES.md § rnaviral_assembly (speed-up plan item 7,
        # bwa-mem3 swap).
        bwa-mem3 index -p "$OUT_DIR/contigs" {output.contigs}
        bwa-mem3 mem -t {threads} --bam=0 "$OUT_DIR/contigs" \
            {input.f1} {input.f2} | \
            samtools view -u -F 4 - | samtools sort -l 1 -m 1G -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule filter_contigs:
    input:
        contigs = "data/03.rnaviralSPAdes/{family}/{sample}/contigs.fasta",
        bam = "data/03.rnaviralSPAdes/{family}/{sample}/sorted.bam"
    output:
        filtered = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    benchmark:
        "benchmarks/filter_contigs/{family}_{sample}.tsv"
    # See ENGINEERING_NOTES.md § filter_contigs for how this figure was
    # picked (not yet sacct-validated -- see HPC_EFFICIENCY_REVIEW.md item
    # 5). samtools idxstats itself isn't threaded -- only seqkit grep uses
    # {threads} below.
    threads: 4
    resources: mem_gb=4, runtime=30
    conda: get_env("ASSEMBLY_ENV")
    shell:
        """
        OUT_DIR=$(dirname "{input.contigs}")
        # See ENGINEERING_NOTES.md § filter_contigs.
        samtools idxstats {input.bam} | awk '$3==0 {{print $1}}' > "$OUT_DIR/zero_coverage_contigs.txt"
        seqkit grep -v -j {threads} -f "$OUT_DIR/zero_coverage_contigs.txt" {input.contigs} -o {output.filtered}
        """


# Item 4 of the conifer_pest_virome capture spec: per-contig coverage mean,
# median and evenness from the assembly BAM.
#
# Declared as a SECOND consumer of the temp() BAM so that the BAM survives
# until these stats exist. Takes filtered_contigs.fasta as well, which both
# orders this rule after filter_contigs and supplies the `retained` flag.
# Stats are emitted for EVERY assembled contig, not just retained ones, so the
# same file serves the stage-count and escalation-fraction reports.
rule contig_coverage_stats:
    input:
        bam = "data/03.rnaviralSPAdes/{family}/{sample}/sorted.bam",
        filtered = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    output:
        tsv = "data/03.rnaviralSPAdes/{family}/{sample}/coverage_stats.tsv",
        idxstats = "data/03.rnaviralSPAdes/{family}/{sample}/idxstats.tsv"
    benchmark:
        "benchmarks/contig_coverage_stats/{family}_{sample}.tsv"
    threads: 2
    resources: mem_gb=8, runtime=60
    conda: get_env("ASSEMBLY_ENV")
    shell:
        # -a: report zero-depth positions too, so breadth and evenness are
        # computed over the full contig length rather than covered bases only.
        # -Q 0 -q 0: no MAPQ/base-quality gate, matching the unfiltered
        # idxstats counts filter_contigs uses.
        #
        # No --max-depth: samtools rewrote `depth` and dropped the depth cap
        # (--max-depth/-d/-m) entirely -- 1.25, which this env ships, errors
        # with "unrecognized option '--max-depth'" and the rule failed in 10 s
        # for every sample (13 failures, coverage_stats at 0/21). Unlimited
        # depth is now the default, so removing the flag keeps the intended
        # "no depth cap" behaviour rather than changing it.
        """
        samtools idxstats {input.bam} > {output.idxstats}
        samtools depth -a -Q 0 -q 0 {input.bam} \
            | python scripts/contig_coverage_stats.py \
                --idxstats {output.idxstats} \
                --depth - \
                --retained {input.filtered} \
                --family {wildcards.family} \
                --sample {wildcards.sample} \
                --out {output.tsv}
        """


# See ENGINEERING_NOTES.md § run_virsorter2 for the thread/runtime tuning
# history. Annotation only since the three-tool taxonomy switch -- its
# score is joined onto every sample's taxonomy_merged.tsv but no longer
# gates what the taxonomy tools see (see ENGINEERING_NOTES.md § three-tool
# taxonomy). The COBRA stage that used to follow it was removed (see
# ENGINEERING_NOTES.md § run_cobra).
rule run_virsorter2:
    # runtime: was 360 (6h), which 142 samples clipped in the 2026-11-02 run and
    # which then aborted the workflow at 86%. Benchmarks over 84 completed runs:
    # median 1.97h, p90 4.86h, max 5.71h -- nothing exceeded 6h, so the limit was
    # marginal rather than the tool hanging. 900 leaves real headroom.
    conda: get_env("TAX_QC_ENV")
    threads: 6
    resources: runtime=900
    retries: 2
    input:
        fasta = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta",
        db = os.path.join(config["databases"]["virsorter2_db"], ".done")
    output:
        tsv = "data/03.rnaviralSPAdes/{family}/{sample}/virsorter2/final-viral-score.tsv",
        fasta = "data/03.rnaviralSPAdes/{family}/{sample}/virsorter2/final-viral-combined.fa",
        dir = directory("data/03.rnaviralSPAdes/{family}/{sample}/virsorter2")
    benchmark:
        "benchmarks/run_virsorter2/{family}_{sample}.tsv"
    shell:
        """
        # VirSorter2 is best-effort: its score is joined on as a column and never
        # removes a contig, but annotate_viral_positives took this file as a hard
        # input, so a failure here blocked that sample's taxonomy_merged.tsv and,
        # through it, MetaCerberus for the whole run. Both declared outputs are
        # therefore always created: on failure they are left empty, and
        # merge_taxonomy_calls.py already treats an empty score file as "no score"
        # and writes NA. Failures are recorded so they can be retried deliberately.
        mkdir -p {output.dir} logs
        # Run under `timeout`, not bare. The handler below only ever caught
        # VirSorter2 exiting non-zero; SLURM's wall-clock kill terminates this
        # shell outright, so a TIMEOUT failed the rule hard and the best-effort
        # contract above did not hold. On the Serratus Hepeviridae run that cost
        # 17 timeouts and 1 d 21 h of a stalled workflow. Owning the kill
        # ourselves means the command returns 124, the handler runs, empty
        # outputs are written and the sample degrades to "no score" exactly as a
        # crash already did. The 600 s margin leaves room for the guarantees below.
        VS2_LIMIT=$(( {resources.runtime} * 60 ))
        # Margin is 600 s, except for short allocations where a flat 600 s would
        # eat most or all of the budget (and go negative under ~10 min), leaving
        # no room for the handler -- there, take a tenth.
        if [ "$VS2_LIMIT" -gt 6000 ]; then VS2_MARGIN=600; else VS2_MARGIN=$(( VS2_LIMIT / 10 )); fi
        VS2_BUDGET=$(( VS2_LIMIT - VS2_MARGIN ))
        VS2_RC=0
        timeout -k 60 "${{VS2_BUDGET}}s" \
            virsorter run -i {input.fasta} -w {output.dir} \
                --db-dir {config[databases][virsorter2_db]} -j {threads} || VS2_RC=$?
        if [ "$VS2_RC" -ne 0 ]; then
            # 124 = timeout sent SIGTERM; 137 = it escalated to SIGKILL.
            if [ "$VS2_RC" -eq 124 ] || [ "$VS2_RC" -eq 137 ]; then
                VS2_REASON="TIMEOUT after ${{VS2_BUDGET}}s"
            else
                VS2_REASON="FAILED rc=$VS2_RC"
            fi
            echo "[virsorter2] $VS2_REASON for {wildcards.family}/{wildcards.sample}; continuing without its score" >&2
            printf '%s\t%s\t%s\t%s\n' "{wildcards.family}" "{wildcards.sample}" "$(date -Is)" "$VS2_REASON" \
                >> logs/virsorter2_failures.tsv
        fi
        # guarantee the declared outputs exist, including after a partial run
        [ -s {output.tsv} ] || : > {output.tsv}
        [ -f {output.fasta} ] || : > {output.fasta}
        """

rule seqkit_stats:
    # Was on raw reads (data/01.reads); moved to the rnaviralSPAdes
    # assembly output per-sample instead -- see ENGINEERING_NOTES.md
    # § seqkit_stats. filtered_contigs.fasta, not contigs.fasta: the
    # latter is temp()-wrapped and deleted once filter_contigs finishes.
    input:
        "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    conda: get_env("ASSEMBLY_ENV")
    output:
        "data/03.assembly_stats/{family}/{sample}/{sample}_stats.txt"
    benchmark:
        "benchmarks/seqkit_stats/{family}_{sample}.tsv"
    shell:
        """
        mkdir -p $(dirname {output})
        seqkit stats {input} > {output}
        """
