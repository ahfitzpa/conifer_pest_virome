rule fastp_cleaning:
    input:
        # Reads straight from staging -- host scrubbing removed from the
        # pipeline entirely, see ENGINEERING_NOTES.md § host-removal pilot.
        unpack(get_fastp_input)
    conda: get_env("ASSEMBLY_ENV")
    output:
        # temp(): rnaviral_assembly is the only consumer; Snakemake deletes
        # these once it finishes. Reports (html/json) are kept.
        fq1_out = temp("data/02.fastp/{family}/{sample}/{sample}_cleaned_1.fastq.gz"),
        fq2_out = temp("data/02.fastp/{family}/{sample}/{sample}_cleaned_2.fastq.gz"),
        report = "data/02.fastp/{family}/{sample}/{sample}_report.html",
        json = "data/02.fastp/{family}/{sample}/{sample}_report.json"
    benchmark:
        "benchmarks/fastp_cleaning/{family}_{sample}.tsv"
    params:
        # Simplified: removed project_id reference
        min_len = 50,
        overlap = 30
    threads: 8
    shell:
        """
        mkdir -p "$(dirname {output.fq1_out})"

        # -z 1: see ENGINEERING_NOTES.md § fastp_cleaning (speed-up plan item 8).
        fastp --in1 {input.fq1} --in2 {input.fq2} \
            --out1 {output.fq1_out} --out2 {output.fq2_out} \
            --html {output.report} \
            --json {output.json} \
            --qualified_quality_phred 20 \
            --length_required {params.min_len} \
            --overlap_len_require {params.overlap} \
            --correction --overlap_diff_percent_limit 2 \
            -z 1 --verbose --thread {threads}
        """