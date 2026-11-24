# Integrating Phase 3 into the workflow

The new files follow the conventions of `github.com/ahfitzpa/conifer_pest_virome`
(rules with `benchmark`, `threads`, `resources`, `conda: get_env(...)`, helper
scripts under `scripts/`, pinned conda environments under `envs/`).

## 1. Copy files into the repo

```
envs/env_phylo.yaml          -> envs/
envs/env_phylo_plot.yaml     -> envs/
envs/env_genome_maps.yaml    -> envs/
envs/env_recomb.yaml         -> envs/     (HyPhy 2.5.104)
envs/env_openrdp.yaml        -> envs/     (numpy<2 + openrdp==0.1.0, separate)
envs/env_demarcation.yaml    -> envs/
configs/phylogenetics.yaml   -> configs/   (now includes genome_maps,
                                recombination and demarcation blocks)
rules/09.phylogenetics.smk  -> rules/
rules/10.genome_maps.smk    -> rules/
rules/11.recombination.smk  -> rules/
rules/12.demarcation.smk    -> rules/
scripts/phylo_prepare.py     -> scripts/
scripts/phylo_support.py     -> scripts/
scripts/phylo_summary.py     -> scripts/
scripts/concat_domains.py   -> scripts/
scripts/diagnostic_residues.py   -> scripts/
scripts/plot_diagnostic_matrix.py -> scripts/
scripts/build_genome_map.py  -> scripts/
scripts/recomb_prepare.py    -> scripts/
scripts/run_openrdp.py       -> scripts/
scripts/run_hyphy_gard.py    -> scripts/
scripts/run_rdp5.py          -> scripts/
scripts/recomb_summary.py    -> scripts/
scripts/demarcation_distances.py -> scripts/
scripts/ictv_demarcation.py  -> scripts/
submit_phase3.sh             -> repo root
```

## 2. Register the environments in `configs/paths.yaml`

Add to the `tools:` block:

```yaml
  PHYLO_ENV: env_phylo.yaml
  PHYLO_PLOT_ENV: env_phylo_plot.yaml
  GENOME_MAPS_ENV: env_genome_maps.yaml
  REC_ENV: env_recomb.yaml
  OPENRDP_ENV: env_openrdp.yaml
  DEMARCATION_ENV: env_demarcation.yaml
```

## 3. Include the rules in `Snakefile_phase2`

After `include: "rules/08.ictv.smk"` add:

```python
include: "rules/09.phylogenetics.smk"
include: "rules/10.genome_maps.smk"
include: "rules/11.recombination.smk"
include: "rules/12.demarcation.smk"
```

Phase 3 is deliberately **not** added to Phase 2's `rule all`; it is run as a
separate target so that a slow IQ-TREE job cannot block the assembly phase.

## 4. Prerequisite outputs (already in Phase 2)

- `results/ictv/families.tsv` and `results/ictv/families/<family>/` — from `ictv_prepare`
- `results/reporting/candidate_contigs.tsv` — from `candidate_contig_table`
- `results/reporting/metacerberus_orfs.faa` — from `export_metacerberus_orfs`
- `results/annotation/metacerberus/final/annotations-FragGeneScan_metacerberus_input/final_annotation_summary.tsv`

ICTV reference proteins are **not** produced by the pipeline; set their paths in
`configs/phylogenetics.yaml` (defaults point at `references/ictv_refs.faa` and
`references/ictv_refs_meta.tsv`, generated from VMR MSL41.v1).

## 5. Run on Sonic

```bash
# from the repo root, on a login node (the scheduler job submits the rules)
sbatch submit_phase3.sh
```

or interactively for a quick test on one family:

```bash
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=mamba
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
snakemake -s Snakefile_phase2 --cores 8 --use-conda \
  --conda-prefix .snakemake/conda \
  results/phylo/families/Dicistroviridae/phylo_summary.tsv
```

## 6. Outputs

Per family, under `results/phylo/families/<family>/`:
- `queries.faa`, `refs.faa`, `combined.faa`, `taxa_metadata.tsv`
- `aln.faa`, `aln.trimmed.faa`
- `tree.newick`, `iqtree.txt`, `iqtree.log`, `support.tsv`, `phylo_summary.tsv`
- `concat_domains.faa`, `domain_presence.tsv`, `domain_schematic.tsv`, `concat.tree.newick`
- `diagnostic_residues.tsv`, `conservation.tsv`

Figures: `figures/diagnostic_<family>.pdf/.svg/.png`

Named metrics (model, UFB/SH-aLRT support rule, taxa counts) are in
`phylo_summary.tsv` — cite those, not unnamed distances.

## 7. Method notes to carry into the manuscript

- MAFFT `--auto`, trimAl `-automated1`, IQ-TREE 2 `-m MFP` (ModelFinder),
  `-B 1000` (UFB), `-alrt 1000` (SH-aLRT), 8 threads.
- Nodes considered supported at SH-aLRT ≥ 90 **and** UFB ≥ 70 (configurable).
- Rooting: trees are rooted as inferred by IQ-TREE; midpoint root only if no
  outgroup is available, and state this explicitly.
- Concatenated domains are selected from MetaCerberus product annotations using
  the keyword lists in `configs/phylogenetics.yaml`; presence/absence is audited
  in `domain_presence.tsv`.
- Diagnostic columns are those where ≥2 groups are fixed for different residues;
  the full column-level statistics are in `conservation.tsv`.
