# See ENGINEERING_NOTES.md § three-tool taxonomy for why mmseqs (nt and
# GPU/nr), ICCTax and the mmseqs viral/orphan partition were replaced, and
# § run_checkv, § run_rdrpcatch for the fix history carried over here.
#
# Kraken2, geNomad and RdRpCATCH all run independently, in parallel, on
# the same coverage-filtered assembly -- none of them gates another, and
# VirSorter2 no longer pre-filters their input (it dropped RNA viruses the
# other tools catch). Every contig ANY of the three calls viral is kept;
# rules/05.merge_tax.smk unions them and CheckV/VirSorter2 annotate.
#
# Output-tracking note: rules here declare only the result files, never
# directory() outputs, so taxonomy_cleanup can purge tool scratch inside
# these dirs without bumping a tracked output's mtime (which, under
# --rerun-triggers mtime, would re-run everything downstream).

def get_taxonomy_contigs(wildcards):
    return f"data/03.rnaviralSPAdes/{wildcards.family}/{wildcards.sample}/filtered_contigs.fasta"

rule run_kraken2:
    conda: get_env("TAX_CORE_ENV")
    threads: 8
    resources: mem_gb=16, runtime=60
    input:
        contigs=get_taxonomy_contigs,
        db=os.path.join(config["databases"]["kraken2_db"], "hash.k2d")
    output:
        tsv="data/04.taxonomy/{family}/{sample}/kraken2/kraken2_output.tsv",
        report="data/04.taxonomy/{family}/{sample}/kraken2/kraken2_report.txt"
    benchmark:
        "benchmarks/run_kraken2/{family}_{sample}.tsv"
    shell:
        """
        kraken2 --db {config[databases][kraken2_db]} --threads {threads} --use-names \
            --report {output.report} --output {output.tsv} {input.contigs}
        """

rule run_genomad:
    # Pinned to the geNomad version the method comparison and CAMI2
    # benchmark ran (envs/genomad.yaml).
    conda: get_env("GENOMAD_ENV")
    threads: 16
    resources: mem_gb=64, runtime=1440
    input:
        contigs=get_taxonomy_contigs,
        db=os.path.join(config["databases"]["genomad_db"], "version.txt")
    output:
        # geNomad names everything after the input fasta's basename.
        summary="data/04.taxonomy/{family}/{sample}/genomad/filtered_contigs_summary/filtered_contigs_virus_summary.tsv"
    benchmark:
        "benchmarks/run_genomad/{family}_{sample}.tsv"
    shell:
        """
        OUT_DIR=data/04.taxonomy/{wildcards.family}/{wildcards.sample}/genomad
        # geNomad resumes from whatever it finds in OUT_DIR; start clean
        # so a retry can't silently mix in a half-finished earlier attempt.
        rm -rf "$OUT_DIR"
        mkdir -p "$(dirname {output.summary})"

        # geNomad hard-errors on an empty fasta; zero contigs surviving
        # filter_contigs is legitimate, not a failure.
        if [ ! -s {input.contigs} ]; then
            printf "seq_name\\tlength\\ttopology\\tcoordinates\\tn_genes\\tgenetic_code\\tvirus_score\\tfdr\\tn_hallmarks\\tmarker_enrichment\\ttaxonomy\\n" > {output.summary}
            exit 0
        fi

        genomad end-to-end --cleanup --splits 8 --threads {threads} \
            {input.contigs} "$OUT_DIR" {config[databases][genomad_db]}
        """

rule run_rdrpcatch:
    conda: get_env("TAX_CORE_ENV")
    threads: 8
    # 120 -> 240: now scans every coverage-filtered contig, not just the
    # mmseqs orphans it used to get. Not yet sacct-validated.
    resources: mem_gb=16, runtime=240
    input:
        contigs=get_taxonomy_contigs,
        db=os.path.join(config["databases"]["rdrp_catch_hmm"], ".done")
    output:
        tsv="data/04.taxonomy/{family}/{sample}/rdrpcatch/rdrp_annotated.tsv"
    benchmark:
        "benchmarks/run_rdrpcatch/{family}_{sample}.tsv"
    # MAX_LEN below is deliberately a shell variable, NOT a params: entry.
    # submit_phase2.sh runs with --rerun-triggers "mtime params input", so
    # adding a params key re-runs every job of this rule -- measured: all 207,
    # including the 32 already complete (~127 CPU-hours). The shell body is
    # not a rerun trigger, so defining it here changes nothing already done.
    shell:
        """
        OUT_DIR=data/04.taxonomy/{wildcards.family}/{wildcards.sample}/rdrpcatch/scan
        # scan.log lives outside OUT_DIR -- --overwrite makes rdrpcatch
        # recreate that directory mid-run, which would orphan a log file
        # tee'd to a path nested inside it.
        SCAN_LOG="$OUT_DIR.scan.log"
        mkdir -p "$(dirname "$OUT_DIR")"

        # RdRpCATCH hard-errors on any sequence over 300 kb ("Exceeds maximum
        # allowed length: 300000 ... will cause issues with the pyHMMER
        # search"), which kills the scan for the whole sample -- it cost
        # SRR16916432 (362,785 bp contig) and SRR16916437 (309,252 bp) on
        # 2026-11-03. Nothing is lost biologically: the longest known RNA
        # virus genomes are the ~41 kb nidoviruses, so a >300 kb contig in an
        # rnaviralSPAdes assembly is host or bacterial carry-through, never an
        # RdRp candidate. Only this rule needs the cut; Kraken2, geNomad and
        # CheckV all take the contigs whole.
        MAX_LEN=299999
        FILTERED="$(dirname "$OUT_DIR")/contigs_le_${{MAX_LEN}}.fasta"
        seqkit seq -M "$MAX_LEN" {input.contigs} -o "$FILTERED" 2>/dev/null
        N_IN=$(grep -c '^>' {input.contigs} || true)
        N_KEPT=$(grep -c '^>' "$FILTERED" || true)
        echo "[run_rdrpcatch] {wildcards.sample}: $N_KEPT of $N_IN contigs <= $MAX_LEN bp"
        if [ "$N_IN" -ne "$N_KEPT" ]; then
            echo "[run_rdrpcatch] dropped $((N_IN - N_KEPT)) contig(s) over $MAX_LEN bp"
        fi

        # Zero scannable contigs is legitimate (an empty or all-over-length
        # assembly), not a failure -- emit the same header-only placeholder
        # the no-hits branch below uses.
        if [ ! -s "$FILTERED" ]; then
            echo "Contig_name" > {output.tsv}
            rm -f "$FILTERED"
            exit 0
        fi

        set +e
        rdrpcatch scan --input "$FILTERED" --db-dir {config[databases][rdrp_catch_hmm]} \
            --seq-type nuc --cpus {threads} --overwrite --output "$OUT_DIR" 2>&1 | tee "$SCAN_LOG"
        set -e

        FILE=$(find "$OUT_DIR" -name "*_rdrpcatch_output_annotated.tsv" | head -n 1)
        # The HMM-search hits are written before rdrpcatch's final UpSet
        # plot, which can crash on an upsetplot/matplotlib mismatch
        # unrelated to the search itself (seen 2026-10-13 in the method
        # comparison) -- fall back to the pre-plot combined hits table.
        COMBINED=$(find "$OUT_DIR" -path "*tsv_files*" -name "*_combined.tsv" | head -n 1)
        # Strip ANSI escapes and collapse whitespace so rich's line-wrapped
        # "no hits" message can't dodge a plain grep.
        SCAN_LOG_FLAT=$(sed 's/\\x1b\\[[0-9;]*[a-zA-Z]//g' "$SCAN_LOG" | tr -s '[:space:]' ' ')
        if [ -f "$FILE" ]; then
            cp "$FILE" {output.tsv}
        elif [ -f "$COMBINED" ]; then
            cp "$COMBINED" {output.tsv}
        elif echo "$SCAN_LOG_FLAT" | grep -q "No hits found by RdRpCATCH"; then
            # Legitimate "nothing matched any HMM" exit, not a failure.
            echo "Contig_name" > {output.tsv}
        else
            echo "Error: RdRpCATCH output not found, and scan.log doesn't show the known no-hits exit -- check $SCAN_LOG" >&2
            exit 1
        fi
        rm -f "$SCAN_LOG" "$FILTERED"
        """

rule run_checkv:
    # Quality/completeness annotation on the union of all three tools'
    # positives (rules/05.merge_tax.smk: collect_viral_positives) -- not a
    # filter, every positive is kept regardless of CheckV's verdict.
    conda: get_env("TAX_QC_ENV")
    threads: 8
    # See ENGINEERING_NOTES.md § run_checkv for how this figure was picked
    # (not yet sacct-validated -- see HPC_EFFICIENCY_REVIEW.md item 5).
    resources: mem_gb=16, runtime=120
    input:
        fasta="data/05.final_selection/{family}/{sample}/viral_positives.fasta",
        db=os.path.join(config["databases"]["checkv_db"], ".done")
    output:
        tsv="data/04.taxonomy/{family}/{sample}/checkv/quality_summary.tsv"
    benchmark:
        "benchmarks/run_checkv/{family}_{sample}.tsv"
    shell:
        # checkv download_database creates a versioned subdir rather than
        # downloading straight into the configured path; resolve it here.
        """
        OUT_DIR=$(dirname {output.tsv})
        if [ ! -s {input.fasta} ]; then
            mkdir -p "$OUT_DIR"
            printf "contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n" > {output.tsv}
            exit 0
        fi
        DB_DIR=$(find {config[databases][checkv_db]} -maxdepth 1 -type d -name 'checkv-db-*' | head -n1)
        checkv end_to_end {input.fasta} "$OUT_DIR" -d "$DB_DIR" -t {threads} --restart
        """
