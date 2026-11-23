# Script audit — what is in the Snakemake chain, what is not

Started 2026-11-23. A working note, not a finished cleanup: nothing listed
here has been deleted. Verified by grepping every script name against the
`.smk` files, the `Snakefile_*`s, the `submit_*.sh` drivers and the configs,
then against the other scripts for Python-level imports.

## In the chain

`scripts/` is reached from `Snakefile_phase2` (rules 00–08, 99);
`common/scripts/` and `analysis_novel/scripts/` from `Snakefile_novel`
(rules 09–14 plus `common/rules/phylogeny.smk`). Everything not named below
is called by a rule.

Worth knowing: `diagnostic_residues.py` and `plot_diagnostic_matrix.py` are
in the chain, via `phylo_all`, not ad hoc. `concat_domains.py` runs only for
families with an entry under `phylogenetics.domains`, so most families skip
it by design.

## Not referenced by any rule or driver

| Script | Status |
|---|---|
| `scripts/phase1_intake.py` | Phase 1 lives in the Serratus tree; nothing here calls it. |
| `scripts/extract_sra_ids.py` | Superseded by `phase1_intake` / the manifest. |
| `scripts/build_family_lineage.py` | Wrote the family lineage table that rule 05 now derives itself. |
| `scripts/retax_metadata_oneoff.py` | One-off from the 2026-11-23 retax; named in `analysis_novel/configs/families.yaml` as provenance only. |
| `common/scripts/filter_ictv_families.py` | Family filtering now happens in `Snakefile_novel` (`FAMILIES` / `MARKER_FAMILIES`). |
| `common/scripts/_phylo_lib.py` | Imported by nothing. Either dead or an intended shared library that was never wired up — check before deleting. |
| `analysis_novel/scripts/fetch_vmr.py` | Superseded by `vmr_xlsx_to_tsv.py` + `fetch_refs_vmr.py`, which `submit_fetch_refs.sh` does call. |
| `analysis_novel/scripts/poster_figures.py` | Superseded by the five `plot_*.py` figure scripts; its driver `submit_poster_figures.sh` was deleted. |

`analysis_novel/scripts/fetch_refs_extra.py` is not called by a rule but is
imported by `fetch_refs_vmr.py`, so it is live.

## Live, but outside Snakemake

The overview figures are produced by scripts no rule and no committed driver
invokes: `plot_multifamily_overview.py`, `plot_conifer_complex_heatmap.py`,
`plot_conifer_virome_breadth.py`, `plot_conifer_convergence.py`,
`plot_phylo_tree.py`. They were last run by hand as SLURM job "overview"
(1182061, 1182120 on 2026-11-23) and their outputs under `figures/overview/`
are the headline figures. Nothing regenerates them when upstream results
change, and the script that submitted them is gone. They need either a rule
and a `figures_all` target, or a committed `submit_figures.sh`.

## Duplication to decide on

Two trees are built per family, from two pipelines:

- `results/phylogeny/{fam}/` — nucleotide, FAMSA → ClipKIT → IQ-TREE
  (`common/rules/phylogeny.smk`), the scenario phylogeny;
- `results/phylo/families/{fam}/` — amino acid RdRp/L, MAFFT → trimAl →
  IQ-TREE (rule 09), the ported phase3 tree.

`Snakefile_novel` says this is deliberate, and the two answer different
questions, but it doubles the alignment and tree cost for ~35 families and
the near-identical directory names (`phylogeny` vs `phylo`) are a standing
trap — the 2026-11-23 misdiagnosis of the demarcation rule came from exactly
that. Renaming one of them would pay for itself.

`recomb_rdp5` is defined but never runs: `phylogenetics.recombination.rdp5_binary`
is empty, which keeps it out of the `recomb_summary` fan-in. Either supply
the binary or drop the rule.

## Untracked clutter on the cluster only

`rules/03.assembly.smk.bak_2026-11-12`, `rules/05.merge_tax.smk.bak`,
`rules/05.merge_tax.smk.bak_2026-11-12`, `configs/paths.yaml.bak`, and the
`*.pre_retax.*` snapshots under `results/`. All gitignored, none in the repo.
`analysis_novel/errors.log` is HyPhy's own error log and is now gitignored too.
