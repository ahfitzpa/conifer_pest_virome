"""Smoke tests for the genome-map, recombination and demarcation helpers."""
from pathlib import Path
import subprocess, sys, tempfile, csv

ROOT = Path(r"C:\Users\Amy\OneDrive - University College Dublin\Documents\conifer_pest_virome")
T = Path(tempfile.mkdtemp(prefix="phase3b_smoke_", dir=r"C:\Users\Amy\AppData\Local\Temp\opencode"))
CFG = ROOT / "configs/phylogenetics.yaml"

# genome-map data
(T / "cand.tsv").write_text("contig_id\trun\tassigned_family\thost_organism\tcheckv_quality\tcheckv_completeness\n"
                            "NODE_1\tSRR1\tDicistroviridae\tHylobius abietis\tHigh-quality\t100\n")
(T / "orf.tsv").write_text("target\tproduct\tbest_hit\tevalue\tscore\tORF_start\tORF_end\n"
                           "SRR1|NODE_1\tRNA-dependent RNA polymerase\tVOG1\t1e-50\t200\t100\t1800\n")
(T / "master.fa").write_text(">SRR1|family|scenario|tax|tools|host|NODE_1\n" + "A"*2000 + "\n")
(T / "cov.tsv").write_text("family\tsample\tcontig_id\tcontig_length\tcoverage_mean\n"
                           "Hylobius_abietis\tSRR1\tNODE_1\t2000\t54.5\n")

def run(script, args):
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / script)] + args,
                       capture_output=True, text=True)
    print("OK " if r.returncode == 0 else "FAIL", script, "|", (r.stderr or "")[:220].strip())
    if r.returncode:
        print((r.stdout or "")[:400])
    return r.returncode

rc = 0
rc |= run("build_genome_map.py", ["--prepare-only", "--family", "Dicistroviridae",
    "--candidates", str(T/"cand.tsv"), "--orf-annotation", str(T/"orf.tsv"),
    "--master-fasta", str(T/"master.fa"), "--coverage", str(T/"cov.tsv"),
    "--config", str(CFG), "--out-table", str(T/"gm.tsv")])
rc |= run("build_genome_map.py", ["--plot-only", "--family", "Dicistroviridae",
    "--config", str(CFG), "--in-table", str(T/"gm.tsv"),
    "--out-pdf", str(T/"gm.pdf"), "--out-svg", str(T/"gm.svg"), "--out-png", str(T/"gm.png")])

# recombination prepare (eligible and ineligible)
(T / "meta.tsv").write_text("taxon_id\tgroup\n" + "".join(
    f"cand{i}\tcandidate\n" for i in range(2)) + "".join(
    f"ref{i}\tICTV\n" for i in range(8)))
(T / "seqs.fa").write_text("".join(f">cand{i}\n" + "A"*600 + "\n" for i in range(2)) +
                           "".join(f">ref{i}\n" + "A"*600 + "\n" for i in range(8)))
rc |= run("recomb_prepare.py", ["--family", "Dicistroviridae", "--sequences", str(T/"seqs.fa"),
    "--metadata", str(T/"meta.tsv"), "--min-refs", "8", "--min-aln", "500",
    "--out-faa", str(T/"rec.fa"), "--out-eligible", str(T/"elig.tsv")])
rc |= run("run_rdp5.py", ["--alignment", str(T/"rec.fa"), "--eligible", str(T/"elig.tsv"),
    "--rdp5-binary", "", "--out-tsv", str(T/"rdp5.tsv")])  # expected FAIL (exit 3)
rc |= run("recomb_summary.py", ["--family", "Dicistroviridae",
    "--event-tables", str(T/"rdp5.tsv"), "--min-methods", "2", "--out", str(T/"rs.tsv")])

# demarcation
(T / "pair.tsv").write_text("query\treference\tidentity_pct\tquery_coverage_pct\tpasses_50pct_coverage\n"
    "cand0\tref0\t52.35\t81.99\t1\ncand0\tref1\t40.0\t90\t1\n")
rc |= run("demarcation_distances.py", ["--alignment", str(T/"seqs.fa"),
    "--metadata", str(T/"meta.tsv"), "--out", str(T/"pair2.tsv")])
rc |= run("ictv_demarcation.py", ["--family", "Dicistroviridae",
    "--pairwise", str(T/"pair.tsv"), "--metadata", str(T/"meta.tsv"),
    "--config", str(CFG), "--out", str(T/"dm.tsv")])
print("temp:", T)
sys.exit(0)
