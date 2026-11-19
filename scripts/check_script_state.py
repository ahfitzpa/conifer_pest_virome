from pathlib import Path

R = Path(r"C:\Users\Amy\OneDrive - University College Dublin\Documents\conifer_pest_virome")
print("=== Phase 2-era script state ===")
for s in ["build_host_trees.py", "expand_analysis.py", "make_overview_figure.py",
          "make_msa_figures.py", "make_fig2_paired.py", "functional_profiling.py",
          "analysis.py", "novelty_compare.py", "make_orf_table.py", "overview_analysis.py"]:
    p = R / "scripts" / s
    print(f"{s:24s} {'EXISTS' if p.exists() else 'DELETED'}")
print("\n=== Phase 3 scripts (written, pending HPC run) ===")
for s in ["build_genome_map.py", "recomb_prepare.py", "run_openrdp.py", "run_hyphy_gard.py",
          "run_rdp5.py", "recomb_summary.py", "demarcation_distances.py",
          "ictv_demarcation.py", "phylo_prepare.py", "phylo_summary.py", "phylo_support.py"]:
    p = R / "scripts" / s
    print(f"{s:24s} {'EXISTS' if p.exists() else 'MISSING'}")
print("\n=== Core inputs still present? ===")
for s in ["results/master_viral_database.fasta", "inputs/taxonomy/master_metadata.csv",
          "results/reporting/candidate_contigs.tsv", "results/reporting/metacerberus_orfs.faa"]:
    print(f"{s:50s} {'EXISTS' if (R/s).exists() else 'MISSING'}")
