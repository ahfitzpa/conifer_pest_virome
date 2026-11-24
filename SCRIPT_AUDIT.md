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

`analysis_novel/scripts/fetch_refs_extra.py` is not called by a rule but is
imported by `fetch_refs_vmr.py`, so it is live.

## Live, but outside Snakemake

The overview figures are produced by scripts no rule invokes. They were last
run by hand as SLURM job "overview" (1182061, 1182120 on 2026-11-23) and
their outputs under `figures/overview/` are the headline figures. Nothing
regenerates them when upstream results change.

They are now filed by scope (2026-11-24):

- `common/scripts/plot_multifamily_overview.py`, `common/scripts/plot_phylo_tree.py`
  — generic, no hardcoded hosts, reusable by the surrogate / host-range /
  geography pathways, beside `plot_signatures.py`;
- `analysis_novel/scripts/plot_conifer_complex_heatmap.py`,
  `plot_conifer_virome_breadth.py`, `plot_conifer_convergence.py` — hardcode
  Locusta / Hylobius / Ips 6–13 times each, so they stay with this pathway.

Resolved 2026-11-24: `figures_all` (rules/15.figures.smk) and
`analysis_novel/submit_figures.sh` both exist and are committed, so the
figures now track their inputs.

## Duplication — resolved 2026-11-24

Two ML trees were built per family, from two pipelines:

- `results/phylogeny/{fam}/` — nucleotide, FAMSA → ClipKIT → IQ-TREE
  (`common/rules/phylogeny.smk`), the scenario phylogeny;
- `results/phylo/families/{fam}/` — amino acid RdRp/L, MAFFT → trimAl →
  IQ-TREE (rule 09), the ported phase3 tree.

Tracing the consumers settled it: rule 09's `tree.newick` fed only its own
`phylo_summary.tsv` and nothing else, while the nucleotide tree feeds the
codon alignments (13), the selection tree-pruning (14), recombination (11)
and every figure. Rules `phylo_tree` and `phylo_summary` are gone, with
`scripts/phylo_support.py` and `scripts/phylo_summary.py`. Rule 09 still
produces `aln.trimmed.faa` for demarcation, the concatenated-domain tree and
the diagnostic residues; `phylo_all` now targets the alignment instead of a
tree. One ML tree per family, ~35 IQ-TREE runs saved, no downstream loss.

The directory names `phylo` vs `phylogeny` remain near-identical and remain a
trap — the 2026-11-23 misdiagnosis of the demarcation rule came from exactly
that. Renaming is still worth doing.

`recomb_rdp5` is gone (2026-11-24). It shelled out to the standalone RDP5
binary, which is not on bioconda; `rdp5_binary` was always empty, so the rule
never ran and was excluded from the `recomb_summary` fan-in. OpenRDP's `RDP`
scanner is the RDP implementation this project uses. `scripts/run_rdp5.py`
deleted with it.

## Untracked clutter on the cluster only

`rules/03.assembly.smk.bak_2026-11-12`, `rules/05.merge_tax.smk.bak`,
`rules/05.merge_tax.smk.bak_2026-11-12`, `configs/paths.yaml.bak`, and the
`*.pre_retax.*` snapshots under `results/`. All gitignored, none in the repo.
`analysis_novel/errors.log` is HyPhy's own error log and is now gitignored too.

## Driver scripts (`.sh`) — 2026-11-24

All 17 are tracked and committed; none is stray. Checked three ways: each
driver's `snakemake` invocation and targets against every other driver's, each
file name grepped across the `.md`, `.sh`, `.smk`, `.yaml` and `Snakefile_*`
files, and each `--job-name` against `sacct` from 2026-11-01 on. "runs" below
is that job-name history, with completions in brackets.

| Script | runs (ok) | Referenced by | Verdict |
|---|---|---|---|
| `push_repo.sh` | 9 (9) | — | Keep. The only way to push: the login node's process limit breaks git's `fork()`. |
| `setup_databases.sh` | — | `SETUP.md`, `rules/00.setup.smk` | Keep. One-time database setup. |
| `submit_phase2.sh` | 4 (4) | rules 03, 04, 05, 06 | Keep. The canonical phase2 driver, and four rule files cite it. |
| `submit_phase2_followup.sh` | 1 (1) | — | **Remove.** Superseded by `submit_phase2.sh`. |
| `submit_metacerberus.sh` | 2 (1) | — | Keep. Isolates an expensive step that a cancelled parent leaves half-done. |
| `submit_contig_confidence.sh` | 1 (1) | — | Keep. The only *proven* driver for `contig_confidence.tsv`. |
| `submit_phase2_retax.sh` | 4 (0) | — | Keep, but unproven — see below. |
| `cleanup_assembly_intermediates.sh` | 1 (1) | — | Keep. Maintenance tool, dry-run by default, filtered by the retention guard. |
| `analysis_novel/submit_novel.sh` | 5 (0) | `submit_demarcation.sh`, `submit_novel_rebuild.sh` | Keep. The master driver, seven stages. |
| `analysis_novel/submit_figures.sh` | 1 (0) | `SCRIPT_AUDIT.md` | Keep. `figures_all` is *not* one of `submit_novel.sh`'s stages, so this is its only driver. |
| `analysis_novel/submit_fetch_refs.sh` | 8 (5) | `Snakefile_novel`, `configs/phylogenetics.yaml` | Keep. Builds the ICTV reference sets. |
| `analysis_novel/check_dag.sh` | 2 (1) | — | Keep. Dry-runs every target; cheap smoke test after a rule-file edit. |
| `analysis_novel/submit_recomb.sh` | 2 (0) | `rules/11.recombination.smk` | Keep. The rule's runtime comment depends on this script's 2-day window. |
| `analysis_novel/submit_demarcation.sh` | 2 (1) | — | Keep unless per-stage reruns are unwanted; `demarcation_all` is a `submit_novel.sh` stage. |
| `analysis_novel/submit_maps_concat.sh` | 3 (0) | — | **Remove.** Both targets are `submit_novel.sh` stages. |
| `analysis_novel/submit_concat_rebuild.sh` | 0 (0) | — | **Remove.** Never ran; superseded before it did. |
| `analysis_novel/submit_novel_rebuild.sh` | 1 (—) | — | Remove once its rebuild has landed; its `--forcerun` list is incident-specific. |

### The three that can go

`submit_phase2_followup.sh` is a strict subset of `submit_phase2.sh`: both end
in the same three invocations (full run, then `PHASE2_ALLOW_PARTIAL_MASTER`
over `master_viral_database.fasta` and `results/annotation/metacerberus`, then
the DAG svg and the report), and `submit_phase2.sh` only adds a `--notemp
--until run_virsorter2` pre-stage in front. Re-running `submit_phase2.sh`
does everything the followup does, because Snakemake skips what is already
built and the pre-stage is then a no-op. Nothing references it, and it has run
once, on 2026-11-21.

`analysis_novel/submit_maps_concat.sh` runs `phylo_all` then
`genome_map_all` — both stages of `submit_novel.sh`. Three runs, none
completed. Its only distinct value is rebuilding the maps without the rest,
which `submit_novel.sh` cannot do selectively.

`analysis_novel/submit_concat_rebuild.sh` was written on 2026-11-24 to force
`concat_prepare` after the contig-id fix, queued as job 1186893, and cancelled
before it ran. `submit_novel_rebuild.sh`, written hours later, forces
`concat_prepare` inside its own `phylo_all` stage, so the earlier script is
redundant and has never executed.

### `submit_phase2_retax.sh` has never completed

Four runs on 2026-11-23, no completions: 1180235 and 1180348 failed (1180348
on a `collect_viral_positives` rule failure), 1180411 was cancelled by hand at
11:43, and 1180569 died on a `LockException` — another Snakemake held the
directory. None of that indicts the script, and its purpose is real: forcing
`collect_viral_positives` and `annotate_viral_positives` after an edit to
`family_taxonomy.py` or `merge_taxonomy_calls.py`, which Snakemake cannot
detect because those run from a `shell:` block and it tracks the command
string, not the script body. But the forced path is unproven, so it needs one
successful run before it can be trusted — or removing, since
`submit_contig_confidence.sh` already builds the same target unforced.

### Overlap that is deliberate, not clutter

`submit_demarcation.sh`, `submit_recomb.sh` and `submit_figures.sh` each drive
one stage that `submit_novel.sh` also runs (except `figures_all`, which it
does not). That overlap is what makes a targeted rerun possible without
re-entering the whole chain, and each carries its own wall time and dependency
note. `submit_recomb.sh` in particular cannot be removed casually:
`rules/11.recombination.smk` documents OpenRDP's 36-hour runtime as fitting
"inside the 2-day parent window in submit_recomb.sh", so deleting the script
orphans the constraint that comment records.
