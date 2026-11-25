# Conifer Pest Virome Discovery

Reproducible Snakemake workflow for assembling, classifying, and phylogenetically characterising putative RNA viruses from archived transcriptomes of conifer-feeding Coleoptera (*Ips confusus*, `PRJNA647034`; *Hylobius abietis*, `PRJNA435680` and `PRJNA783427`).

This repository is the study-specific companion workflow derived from the modular [virus_discovery](https://github.com/ahfitzpa/virus_discovery) pipeline.

## Repository Map

```text
.
├── Snakefile_phase2          # Assembly, multi-tool viral classification, and annotation entry point
├── rules/                    # Modular Phase 2 Snakemake rules (00.setup … 08.ictv, 99.retention)
├── scripts/                  # Manifest intake, taxonomy merging, RdRp extraction, retention guard
├── configs/                  # Workflow configs (see paths.example.yaml for portable DB setup)
├── envs/                     # Pinned per-rule Conda environment definitions
├── common/                   # Shared phylogeny, ORF selection, and genomic-signature utilities
├── analysis_novel/           # Downstream characterisation (phylogenetics, recombination, demarcation, signatures)
│   ├── Snakefile_novel       # Downstream Snakemake entry point
│   ├── rules/                # Stages 09–14 (phylogenetics, genome maps, recombination, demarcation, signatures)
│   ├── scripts/              # Phylogeny, ORF verification, and plotting scripts
│   ├── envs/                 # Pinned per-rule Conda environments
│   └── tests/                # Fast unit tests and Snakemake DAG dry-run checks
├── setup_databases.sh        # Idempotent reference database downloader (~35–50 GB)
├── SETUP.md                  # Detailed reference database setup guide
├── CITATION.cff              # Machine-readable citation metadata
└── LICENSE                   # GPL-3.0 License
```

## Quickstart for Reviewers

1. **Provision external reference databases** (not tracked in git; see [SETUP.md](SETUP.md)):
   ```bash
   ./setup_databases.sh /path/to/your/database_root
   ```
   Or copy `configs/paths.example.yaml` to `configs/paths.yaml` and point `databases:` at your local install.
2. **Run Phase 2 assembly & taxonomic classification** (SLURM HPC or local Snakemake 8.x + Conda/Mamba):
   ```bash
   sbatch submit_phase2.sh
   # Optional email notifications: sbatch --mail-user=you@example.org submit_phase2.sh
   ```
3. **Run downstream novel-virus characterisation (`analysis_novel`):**
   ```bash
   cd analysis_novel
   sbatch submit_novel.sh
   ```
4. **Run fast validation checks (no HPC or databases required):**
   ```bash
   python3 -m unittest discover -s analysis_novel/tests -p "test_unit_*.py" -v
   bash analysis_novel/tests/dryrun_novel.sh
   ```

## Reproducibility & Portability

- **External data only:** Raw FASTQ reads, reference databases, container images (`*.sif`), cluster logs (`logs/`), and generated outputs (`results/`) are gitignored and never bundled in the repository.
- **Pinned environments:** Every Snakemake rule declares a version-pinned Conda environment under `envs/`, `common/envs/`, or `analysis_novel/envs/`.
- **Portable cluster scripts:** SLURM wrappers use `${XDG_CACHE_HOME:-$PWD/.cache}`, `${PIP_CACHE_DIR:-$PWD/.cache/pip}`, and `${SLURM_SUBMIT_DIR:-$PWD}` with no hardcoded personal email addresses or home paths.

## License & Citation

Licensed under the GNU General Public License v3.0 ([LICENSE](LICENSE)). Citation metadata is provided in [CITATION.cff](CITATION.cff).
