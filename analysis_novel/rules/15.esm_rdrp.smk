# =============================================================================
# Stage 15: ESMFold + Foldseek Structural Validation of Divergent RdRp Domains
# =============================================================================
# Gated strictly to medium/high-confidence novel contigs that are deeply
# divergent (<40% aa identity to nearest ICTV/nr reference) OR family-unresolved
# RdRp-positive contigs (>=200 aa domain).
#
# Run directly on a Sonic GPU node:
#   snakemake -s Snakefile_novel --use-conda --cores 4 esm_rdrp_all
# Or enable inside novel_all by setting `esm_rdrp: enabled: true` in novel.yaml.
# =============================================================================

import os

ESM_CFG = config.get("esm_rdrp", {})
ESM_DIR = ESM_CFG.get("out_dir", "results/esm_rdrp")
ESM_ENV = get_env("ESM_RDRP_ENV")

ESM_CONFIDENCE_TSV = os.path.join(REPORT_DIR, "contig_confidence.tsv")
ESM_RDRP_GLOB = os.path.join(PHASE2, "data/06.rdrp/*/*/*_rdrp_domains.faa")
ESM_DEMARC_GLOBS = ESM_CFG.get(
    "demarcation_tsvs",
    [
        "results/demarcation/all_families_demarcation_summary.tsv",
        "results/demarcation/*/demarcation_summary.tsv",
    ],
)
ESM_FOLDSEEK_DB = ESM_CFG.get("foldseek_db", os.path.join(PHASE2, "references/foldseek/pdb"))
ESM_FOLDSEEK_RT_DB = ESM_CFG.get("foldseek_rt_db", "")
ESM_BFVD_DB = ESM_CFG.get("foldseek_bfvd_db", os.path.join(PHASE2, "references/foldseek/bfvd"))
ESM_DEDICATED_RDRP_DB = bool(ESM_CFG.get("dedicated_rdrp_db", False))


rule esm_rdrp_all:
    input:
        validation = f"{ESM_DIR}/esm_rdrp_validation.tsv",
        summary    = f"{ESM_DIR}/esm_rdrp_summary.tsv",


rule esm_rdrp_bfvd_all:
    """Complementary CPU-only Foldseek search of folded RdRps against BFVD (Kim et al., NAR 2026)."""
    input:
        bfvd_hits = f"{ESM_DIR}/foldseek_bfvd_hits.tsv",


rule esm_rdrp_select:
    """Gate medium/high-confidence, deeply divergent or family-unresolved RdRp domains."""
    input:
        confidence = ESM_CONFIDENCE_TSV,
    output:
        faa = f"{ESM_DIR}/selected_divergent_rdrp.faa",
        tsv = f"{ESM_DIR}/selected_divergent_rdrp.tsv",
    params:
        rdrp_glob        = ESM_RDRP_GLOB,
        demarc_flags     = " ".join(f'--demarcation-tsv "{p}"' for p in ESM_DEMARC_GLOBS),
        keep_confidence  = ESM_CFG.get("keep_confidence", "high,medium"),
        min_aa_len       = int(ESM_CFG.get("min_aa_len", 200)),
        max_aa_len       = int(ESM_CFG.get("max_aa_len", 800)),
        max_identity_pct = float(ESM_CFG.get("max_identity_pct", 40.0)),
        unresolved_flag  = "" if ESM_CFG.get("include_unresolved", True) else "--no-unresolved",
    log:
        "logs/esm_rdrp/select_divergent_rdrp.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.faa})
        python3 scripts/select_divergent_rdrp.py \
            --confidence-tsv "{input.confidence}" \
            --rdrp-faa "{params.rdrp_glob}" \
            {params.demarc_flags} \
            --keep-confidence "{params.keep_confidence}" \
            --min-aa-len {params.min_aa_len} \
            --max-aa-len {params.max_aa_len} \
            --max-identity-pct {params.max_identity_pct} \
            {params.unresolved_flag} \
            --out-faa "{output.faa}" \
            --out-tsv "{output.tsv}" \
            > "{log}" 2>&1
        """


rule esm_rdrp_fold:
    """Fold selected divergent RdRp domains locally on GPU with ESMFold."""
    input:
        faa = f"{ESM_DIR}/selected_divergent_rdrp.faa",
    output:
        metrics = f"{ESM_DIR}/esmfold_metrics.tsv",
        done    = f"{ESM_DIR}/structures/.folded.done",
    params:
        pdb_dir    = f"{ESM_DIR}/structures",
        model_name = ESM_CFG.get("model_name", "facebook/esmfold_v1"),
        cache_dir  = ESM_CFG.get("hf_cache_dir", ""),
        chunk_size = int(ESM_CFG.get("chunk_size", 64)),
    resources:
        gpu = 1,
    threads: 4
    conda:
        ESM_ENV
    log:
        "logs/esm_rdrp/esmfold_rdrp.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.pdb_dir}"
        CACHE_ARG=""
        if [ -n "{params.cache_dir}" ]; then
            CACHE_ARG="--cache-dir {params.cache_dir}"
        elif [ -n "${{HF_HOME:-}}" ]; then
            CACHE_ARG="--cache-dir ${{HF_HOME}}"
        fi
        python3 scripts/run_esmfold_rdrp.py \
            --in-faa "{input.faa}" \
            --out-pdb-dir "{params.pdb_dir}" \
            --out-metrics "{output.metrics}" \
            --model-name "{params.model_name}" \
            --chunk-size {params.chunk_size} \
            $CACHE_ARG \
            > "{log}" 2>&1
        touch "{output.done}"
        """


rule esm_rdrp_foldseek_rdrp:
    """Search predicted RdRp structures against PDB / RdRp Foldseek structure DB."""
    input:
        done = f"{ESM_DIR}/structures/.folded.done",
    output:
        hits = f"{ESM_DIR}/foldseek_rdrp_hits.tsv",
    params:
        pdb_dir = f"{ESM_DIR}/structures",
        db      = ESM_FOLDSEEK_DB,
        tmpdir  = f"{ESM_DIR}/tmp_foldseek_rdrp",
        evalue  = float(ESM_CFG.get("foldseek_evalue", 1e-2)),
    threads: 8
    conda:
        ESM_ENV
    log:
        "logs/esm_rdrp/foldseek_rdrp.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.tmpdir}"
        if [ -z "$(find "{params.pdb_dir}" -maxdepth 1 -name '*.pdb' -size +0c -print -quit)" ]; then
            : > "{output.hits}"
            echo "[foldseek_rdrp] No PDB structures present; wrote empty hits file." > "{log}"
            exit 0
        fi
        if [ ! -e "{params.db}" ] && [ ! -e "{params.db}.dbtype" ]; then
            if [ "${{SETUP_PDB_DB:-0}}" = "1" ]; then
                mkdir -p "$(dirname "{params.db}")" "{params.tmpdir}_dl"
                foldseek databases PDB "{params.db}" "{params.tmpdir}_dl" --threads {threads} > "{log}" 2>&1
            else
                echo "[foldseek_rdrp] ERROR: Foldseek DB not found at {params.db}" > "{log}"
                echo "Re-run with SETUP_PDB_DB=1 or run: foldseek databases PDB {params.db} {params.tmpdir}" >> "{log}"
                exit 1
            fi
        fi
        foldseek easy-search \
            "{params.pdb_dir}" \
            "{params.db}" \
            "{output.hits}" \
            "{params.tmpdir}" \
            --threads {threads} \
            -e {params.evalue} \
            --format-output "query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore,lddt,prob,theader" \
            > "{log}" 2>&1
        """


rule esm_rdrp_foldseek_rt:
    """Optional negative-control Foldseek search against a dedicated RT structure DB."""
    input:
        done = f"{ESM_DIR}/structures/.folded.done",
    output:
        hits = f"{ESM_DIR}/foldseek_rt_hits.tsv",
    params:
        pdb_dir = f"{ESM_DIR}/structures",
        db      = ESM_FOLDSEEK_RT_DB,
        tmpdir  = f"{ESM_DIR}/tmp_foldseek_rt",
        evalue  = float(ESM_CFG.get("foldseek_evalue", 1e-2)),
    threads: 8
    conda:
        ESM_ENV
    log:
        "logs/esm_rdrp/foldseek_rt.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.tmpdir}"
        if [ -z "{params.db}" ]; then
            : > "{output.hits}"
            echo "[foldseek_rt] No dedicated RT DB configured; relying on PDB header RT screen." > "{log}"
            exit 0
        fi
        if [ -z "$(find "{params.pdb_dir}" -maxdepth 1 -name '*.pdb' -size +0c -print -quit)" ]; then
            : > "{output.hits}"
            echo "[foldseek_rt] No PDB structures present; wrote empty hits file." > "{log}"
            exit 0
        fi
        foldseek easy-search \
            "{params.pdb_dir}" \
            "{params.db}" \
            "{output.hits}" \
            "{params.tmpdir}" \
            --threads {threads} \
            -e {params.evalue} \
            --format-output "query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore,lddt,prob,theader" \
            > "{log}" 2>&1
        """


rule esm_rdrp_summarize:
    """Combine ESMFold pLDDT/pTM, catalytic palm motifs, and Foldseek RdRp/RT hits."""
    input:
        selected = f"{ESM_DIR}/selected_divergent_rdrp.tsv",
        metrics  = f"{ESM_DIR}/esmfold_metrics.tsv",
        rdrp     = f"{ESM_DIR}/foldseek_rdrp_hits.tsv",
        rt       = f"{ESM_DIR}/foldseek_rt_hits.tsv",
    output:
        validation = f"{ESM_DIR}/esm_rdrp_validation.tsv",
        summary    = f"{ESM_DIR}/esm_rdrp_summary.tsv",
    params:
        dedicated_flag = "--dedicated-rdrp-db" if ESM_DEDICATED_RDRP_DB else "",
    log:
        "logs/esm_rdrp/summarize_esm_rdrp.log",
    shell:
        """
        mkdir -p $(dirname {log})
        python3 scripts/summarize_esm_rdrp.py \
            --selected-tsv "{input.selected}" \
            --metrics-tsv "{input.metrics}" \
            --foldseek-rdrp "{input.rdrp}" \
            --foldseek-rt "{input.rt}" \
            {params.dedicated_flag} \
            --out-validation "{output.validation}" \
            --out-summary "{output.summary}" \
            > "{log}" 2>&1
        """


rule esm_rdrp_foldseek_bfvd:
    """Complementary CPU-only Foldseek search against BFVD (Big Fantastic Virus Database).

    Retains PDB + RT negative-control triage as the primary filter (since BFVD
    contains only viral structures and cannot flag host helicases or RTs), while
    providing broad ICTV viral-family structural neighbours across predicted
    viral proteomes (Kim et al., NAR 2026).

    Run directly on Sonic without re-folding:
      SETUP_BFVD_DB=1 snakemake -s Snakefile_novel --cores 8 --use-conda \
        --conda-prefix .snakemake/conda --allowed-rules esm_rdrp_foldseek_bfvd esm_rdrp_bfvd_all \
        -- esm_rdrp_bfvd_all
    """
    input:
        done = f"{ESM_DIR}/structures/.folded.done",
    output:
        hits = f"{ESM_DIR}/foldseek_bfvd_hits.tsv",
    params:
        pdb_dir = f"{ESM_DIR}/structures",
        db      = ESM_BFVD_DB,
        tmpdir  = f"{ESM_DIR}/tmp_foldseek_bfvd",
        evalue  = float(ESM_CFG.get("foldseek_evalue", 1e-2)),
    threads: 8
    conda:
        ESM_ENV
    log:
        "logs/esm_rdrp/foldseek_bfvd.log",
    shell:
        """
        mkdir -p $(dirname {log}) "{params.tmpdir}"
        if [ -z "$(find "{params.pdb_dir}" -maxdepth 1 -name '*.pdb' -size +0c -print -quit)" ]; then
            printf "query\ttarget\tfident\talnlen\tevalue\tbits\tqtmscore\tttmscore\talntmscore\tlddt\tprob\ttheader\n" > "{output.hits}"
            echo "[foldseek_bfvd] No PDB structures present; wrote empty hits file." > "{log}"
            exit 0
        fi
        if [ ! -e "{params.db}" ] && [ ! -e "{params.db}.dbtype" ]; then
            if [ "${{SETUP_BFVD_DB:-0}}" = "1" ]; then
                mkdir -p "$(dirname "{params.db}")" "{params.tmpdir}_dl"
                foldseek databases BFVD "{params.db}" "{params.tmpdir}_dl" --threads {threads} > "{log}" 2>&1
                rm -rf "{params.tmpdir}_dl"
            else
                echo "[foldseek_bfvd] ERROR: BFVD Foldseek DB not found at {params.db}" > "{log}"
                echo "Re-run with SETUP_BFVD_DB=1 or run: foldseek databases BFVD {params.db} {params.tmpdir}_dl" >> "{log}"
                exit 1
            fi
        fi
        printf "query\ttarget\tfident\talnlen\tevalue\tbits\tqtmscore\tttmscore\talntmscore\tlddt\tprob\ttheader\n" > "{output.hits}"
        foldseek easy-search \
            "{params.pdb_dir}" \
            "{params.db}" \
            "{params.tmpdir}/raw_hits.tsv" \
            "{params.tmpdir}" \
            --threads {threads} \
            -e {params.evalue} \
            --max-seqs 10 \
            --format-output "query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore,lddt,prob,theader" \
            >> "{log}" 2>&1
        cat "{params.tmpdir}/raw_hits.tsv" >> "{output.hits}"
        rm -rf "{params.tmpdir}"
        """
