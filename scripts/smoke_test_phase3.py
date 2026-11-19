"""Tiny synthetic inputs to smoke-test the Phase 3 helper scripts."""
from pathlib import Path
import subprocess, sys, tempfile

ROOT = Path(r"C:\Users\Amy\OneDrive - University College Dublin\Documents\conifer_pest_virome")
T = Path(tempfile.mkdtemp(prefix="phase3_smoke_", dir=r"C:\Users\Amy\AppData\Local\Temp\opencode"))

(T / "rdrp.faa").write_text(">SRR1|NODE_1|frame=1|profile=Dicistroviridae.aln\nMKLMNPQRSTV\n>SRR1|NODE_2|frame=2|profile=Dicistroviridae.aln\nMKLMNPQRSTV\n")
(T / "ictv.faa").write_text(">ICTV|Dicistroviridae|Aparavirus|Aparavirus_apisacutum|AF150629|rdrp|1\nMKLMNPQRSTV\n>ICTV|Iflaviridae|Iflavirus|x|AB000001|rdrp|2\nAAAAAAAAAAA\n")
(T / "ictv_meta.tsv").write_text("family\tgenus\tspecies\taccession\tproduct\nDicistroviridae\tAparavirus\tAparavirus apisacutum\tAF150629\tRNA-dependent RNA polymerase\n")
(T / "candidates.tsv").write_text("contig_id\trun\tassigned_family\thost_organism\nNODE_1\tSRR1\tDicistroviridae\tHylobius abietis\nNODE_2\tSRR1\tDicistroviridae\tHylobius abietis\n")
(T / "orfs.faa").write_text(">SRR1|NODE_1|ORF1\nMKLMNPQRSTVVVV\n")
(T / "orf_annotation.tsv").write_text("target\tproduct\tscore\tevalue\tORF_start\tORF_end\nSRR1|NODE_1\tRNA-dependent RNA polymerase\t100\t1e-50\t100\t900\n")
cfg = ROOT / "configs/phylogenetics.yaml"

cmds = [
    [sys.executable, str(ROOT / "scripts/phylo_prepare.py"), "--family", "Dicistroviridae",
     "--rdrp-faa", str(T / "rdrp.faa"), "--ictv-faa", str(T / "ictv.faa"),
     "--ictv-meta", str(T / "ictv_meta.tsv"), "--candidates", str(T / "candidates.tsv"),
     "--out-queries", str(T / "q.faa"), "--out-refs", str(T / "r.faa"),
     "--out-combined", str(T / "c.faa"), "--out-metadata", str(T / "m.tsv")],
    [sys.executable, str(ROOT / "scripts/phylo_support.py"), "--iqtree", str(T / "c.faa"), "--out", str(T / "s.tsv")],
    [sys.executable, str(ROOT / "scripts/concat_domains.py"), "--family", "Dicistroviridae",
     "--orf-faa", str(T / "orfs.faa"), "--orf-annotation", str(T / "orf_annotation.tsv"),
     "--ictv-faa", str(T / "ictv.faa"), "--ictv-meta", str(T / "ictv_meta.tsv"),
     "--candidates", str(T / "candidates.tsv"), "--rdrp-faa", str(T / "rdrp.faa"),
     "--config", str(cfg), "--out-faa", str(T / "concat.faa"),
     "--out-matrix", str(T / "dm.tsv"), "--out-schematic", str(T / "ds.tsv")],
]
for c in cmds:
    r = subprocess.run(c, capture_output=True, text=True)
    print("OK " if r.returncode == 0 else "FAIL ", Path(c[1]).name, "|", (r.stderr or "")[:200].strip())
    if r.returncode:
        print(r.stdout[:300])

# diagnostic test (own alignment + metadata)
(T / "aln.faa").write_text(">cand1\nMKLMNPQRSTV\n>cand2\nMKLMNPQRSTV\n>ICTV1\nMKLMNPQRSTA\n")
(T / "meta.tsv").write_text("taxon_id\tgroup\ncand1\tcandidate\ncand2\tcandidate\nICTV1\tICTV\n")
r = subprocess.run([sys.executable, str(ROOT / "scripts/diagnostic_residues.py"),
    "--alignment", str(T / "aln.faa"), "--metadata", str(T / "meta.tsv"),
    "--config", str(cfg), "--out-matrix", str(T / "diag.tsv"),
    "--out-conservation", str(T / "cons.tsv")], capture_output=True, text=True)
print("OK " if r.returncode == 0 else "FAIL ", "diagnostic_residues.py |", (r.stderr or "")[:200].strip())
if r.returncode:
    print(r.stdout[:300])
print("temp:", T)
