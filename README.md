# Conifer Pest Virome Discovery

Reproducible, multi-stage Snakemake workflow for uncovering, structurally validating, and biogeographically profiling divergent RNA viruses across archived transcriptomes of conifer-feeding Coleoptera across the Holarctic (*Ips confusus*, *Ips typographus*, *Ips nitidus*, *Ips subelongatus*, and *Hylobius abietis*), benchmarked against an out-of-order Orthoptera control (*Locusta migratoria*).

Derived from the modular [virus_discovery](https://github.com/ahfitzpa/virus_discovery) pipeline.

---

## End-to-End Analysis Workflow

The repository is organized into **four sequential, reproducible stages**:

| Stage | Entry Point | Purpose | Key Outputs |
| :--- | :--- | :--- | :--- |
| **1. Assembly & Multi-Tool Taxonomy** | `Snakefile_phase2` (`rules/00..08`) | `fastp` QC, `rnaviralSPAdes` assembly, `geNomad`, `RdRpCATCH`, `Kraken2`, `MetaCerberus`, and ICTV `MSL41.v1` classification | `data/03.rnaviralSPAdes/`, `data/06.rdrp/`, `results/contig_confidence.tsv` |
| **2. 3D Structural RdRp Validation** | `analysis_novel/Snakefile_novel` (`rules/15.esm_rdrp.smk`) | Local GPU ESMFold (`facebook/esmfold_v1`) + Foldseek PDB (with RT negative controls) + complementary BFVD v3 viral atlas search | `results/esm_rdrp/esm_rdrp_validation.tsv`, `foldseek_bfvd_hits.tsv`, `structures/*.pdb` |
| **3. RdRp-Guided SRA Retrieval & Open-Host Profiling** | `snowball/Snakefile_snowball` | Maps structurally validated RdRps to PalmDB/Serratus sOTUs to retrieve Round-2 beetle libraries, profiles host breadth across 3.8M SRA runs (Fisher's exact enrichment), and runs compact Option-B spot-check assemblies | `snowball/results/snowball/manifest_combined.csv`, `snowball_vs_seed_blastx.tsv`, `snowball/results/open_host/open_host_enrichment.tsv` |
| **4. Phylogenetics, Evolution, Benchmarks & Figures** | `analysis_novel/Snakefile_novel` (`rules/09..14, 16..17`) | Amino-acid ML phylogenies (`MAFFT`/`trimAl`/`IQ-TREE2`), ICTV demarcation, *Ips* *Polycipiviridae* phylogeography & `dN/dS`, ablation/out-of-order benchmarks, and publication SVGs (Figs. 1–6) | `results/phylogenetics/`, `results/benchmarks/02..05/`, `figures/manuscript/fig1..6.svg` |

---

## Repository Map

```text
.
├── Snakefile_phase2               # Stage 1: Read QC, rnaviralSPAdes assembly, multi-tool classification
├── rules/                         # Stage 1 modular rules (00.setup ... 08.ictv, 99.retention)
├── scripts/                       # Intake, taxonomy consensus, RdRp extraction, retention helpers
├── configs/                       # Portable workflow configs (see paths.example.yaml)
├── envs/                          # Pinned per-rule Conda environments
├── common/                        # Shared ORF, phylogeny, and genomic-signature utilities
│
├── analysis_novel/                # Stages 2 & 4: Structural validation, phylogenetics, benchmarks, figures
│   ├── Snakefile_novel            # Main downstream Snakemake entry point
│   ├── rules/
│   │   ├── 09.phylogenetics.smk   # MAFFT / trimAl / IQ-TREE2 amino-acid trees
│   │   ├── 10.genome_maps.smk     # ORF architecture & domain maps
│   │   ├── 11.recombination.smk   # Recombination screening
│   │   ├── 12.demarcation.smk     # ICTV species-demarcation identity tables
│   │   ├── 13.codon_alignments.smk# Codon alignments & evolutionary constraint
│   │   ├── 14.genomic_signatures.smk # Dinucleotide / CpG / codon-usage signatures
│   │   ├── 15.esm_rdrp.smk        # ESMFold + Foldseek (PDB+RT & BFVD v3)
│   │   ├── 16.manuscript_figures.smk # Reproducible publication SVGs (Figs. 1–6)
│   │   └── 17.benchmarks.smk      # Ablation, OOS L. migratoria, retrieval gain, Polycipiviridae dN/dS
│   ├── benchmark_data/            # Preserved Round-1 284-domain reference table for reproducible merging
│   ├── scripts/                   # Pure-Python validation, benchmarking, evolution, and SVG scripts
│   ├── envs/                      # Pinned Conda environments (incl. ESMFold/Foldseek)
│   └── tests/                     # 26 automated unit tests + Snakemake dry-run checks
│
├── snowball/                      # Stage 3: RdRp-guided SRA retrieval & 3.8M-run open-host profiling
│   ├── Snakefile_snowball         # PalmDB mapping, Serratus SQL retrieval, open-host enrichment, spot-check
│   ├── scripts/                   # serratus_snowball.py, open_host_specificity.py, open_host_enrichment.py
│   ├── configs/                   # Serratus retrieval caps & thresholds
│   ├── submit_snowball.sh         # SLURM wrapper for Round-2 retrieval
│   └── submit_spotcheck.sh        # SLURM wrapper for Option-B open-host spot-check assemblies
│
├── setup_databases.sh             # Idempotent reference database downloader
├── SETUP.md                       # Detailed database & HPC installation guide
├── CITATION.cff                   # Citation metadata
└── LICENSE                        # GPL-3.0 License
```

---

## Reproducible Quickstart

### 0. Verify installation (local, no HPC or databases required)
```bash
python3 -m unittest discover -s analysis_novel/tests -p "test_unit_*.py" -v
```

### 1. Stage 1 — Assembly & multi-tool classification
Provision external databases via `./setup_databases.sh /path/to/db_root` (see [SETUP.md](SETUP.md)), then run:
```bash
sbatch submit_phase2.sh
```

### 2. Stage 2 — Local GPU ESMFold + Foldseek (PDB + RT controls & BFVD v3)
```bash
cd analysis_novel
sbatch submit_esm_rdrp_gpu.sh
# Complementary CPU-only BFVD v3 Foldseek search:
SETUP_BFVD_DB=1 snakemake -s Snakefile_novel --use-conda --cores 4 -- esm_rdrp_bfvd_all
```

### 3. Stage 3 — RdRp-guided SRA retrieval & open-host enrichment
```bash
cd ../snowball
# Preview Round-2 candidate beetle libraries:
snakemake -s Snakefile_snowball --use-conda --cores 4 -- snowball_discover
# Open-host specificity + Fisher's exact enrichment across 3.8M SRA runs:
snakemake -s Snakefile_snowball --use-conda --cores 4 -- open_host_enrichment
# Optional Option-B compact spot-check assemblies:
sbatch submit_spotcheck.sh
```

### 4. Stage 4 — Phylogenetics, benchmarks, Polycipiviridae evolution & Figs. 1–6
```bash
cd ../analysis_novel
sbatch submit_novel.sh
# Benchmarks (structural ablation, OOS L. migratoria, retrieval gain, Polycipiviridae dN/dS):
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 snakemake -s Snakefile_novel --cores 2 --rerun-triggers mtime -- \
  benchmarks_all results/benchmarks/05_polycipiviridae_evolution/polycipiviridae_evolution.svg
# Reproducible manuscript SVG figures (Figs. 1–6):
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 snakemake -s Snakefile_novel --cores 2 --rerun-triggers mtime -- \
  manuscript_figures_complete
```

---

## Reproducibility & Portability

- **No raw data or bloated outputs in git:** FASTQs, reference DBs, container images (`*.sif`), logs, and generated `results/` are gitignored; only the compact 124 KB Round-1 284-domain structural benchmark table is tracked under `analysis_novel/benchmark_data/`.
- **Pinned environments & CI:** Every rule pins its Conda environment, and GitHub Actions (`.github/workflows/ci.yml`) runs all 26 unit tests on every push and pull request.
- **Zero-OpenBLAS pure-Python plotting & statistics:** All manuscript SVG figures, Fisher's exact enrichment tests, and Nei–Gojobori `dN/dS` calculations run in pure Python without BLAS/R runtime issues on HPC login nodes.

## License & Citation

Licensed under GPL-3.0 ([LICENSE](LICENSE)). Citation metadata is provided in [CITATION.cff](CITATION.cff).
