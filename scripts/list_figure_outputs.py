from pathlib import Path

R = Path(r"C:\Users\Amy\OneDrive - University College Dublin\Documents\conifer_pest_virome")

print("=== Figure/PDF inventory (production figures, not paper_examples) ===")
for d in ["host_trees", "expanded_output", "figure_outputs", "novelty_output", "analysis_output", "figures"]:
    base = R / d
    if not base.exists():
        print(f"{d:16s} MISSING")
        continue
    files = list(base.rglob("*"))
    png = [p for p in files if p.suffix == ".png"]
    pdf = [p for p in files if p.suffix == ".pdf"]
    print(f"{d:16s} png={len(png):3d} pdf={len(pdf):3d}")
    for p in sorted(png):
        print("    PNG", p.name)
    for p in sorted(pdf):
        print("    PDF", p.name)

print("\n=== Phase 3 run outputs (none yet -- scripts written, not executed) ===")
for d in ["genome_maps", "recombination", "demarcation"]:
    base = R / "results" / d
    print(f"-- {d} --")
    for p in base.rglob("*"):
        if p.is_file():
            print("  ", p.relative_to(R))
    if not any(base.rglob("*")):
        print("   (empty)")

print("\n=== Phase 3 conserved-domain concat trees ===")
for p in list(R.rglob("concat*.nwk")) + list(R.rglob("concat.tree.newick")):
    print("  ", p.name)
print("   (none -- Phase 3 scripts written, not executed)")
