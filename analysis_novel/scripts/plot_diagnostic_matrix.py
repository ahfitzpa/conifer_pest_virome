#!/usr/bin/env python3
"""Plot a Zhao et al. (Fig. 4)-style diagnostic-residue matrix.

One row per diagnostic alignment column, one column per taxon, residues coloured
by amino-acid class, with a group strip above the matrix and a conservation
track.  Vector (PDF/SVG) and 300-dpi PNG output.
"""
from __future__ import annotations
import argparse, csv, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

AA = {"hydrophobic": "#E69F00", "positive": "#0072B2", "negative": "#D55E00",
      "polar": "#009E73", "special": "#CC79A7", "gap": "#e7e3da"}
CLASS = {}
for c in "AVILMFWYC": CLASS[c] = "hydrophobic"
for c in "KRH": CLASS[c] = "positive"
for c in "DE": CLASS[c] = "negative"
for c in "STNQ": CLASS[c] = "polar"
for c in "GP": CLASS[c] = "special"
GROUP_COLOR = {"candidate": "#e34948", "ICTV": "#2a78d6"}


def read_tsv(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def iter_fasta(path):
    name, chunks = None, []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--conservation", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-pdf", required=True)
    ap.add_argument("--out-svg", required=True)
    ap.add_argument("--out-png", required=True)
    a = ap.parse_args()

    matrix = read_tsv(a.matrix)
    cons = {int(r["position"]): float(r["conservation"]) for r in read_tsv(a.conservation)}
    seqs = dict(iter_fasta(a.alignment))
    meta = {r["taxon_id"]: r for r in read_tsv(a.metadata)}
    taxa = sorted(seqs)
    diag = [int(r["position"]) for r in matrix if r.get("diagnostic") == "1"]
    diag = diag[:60]  # keep the panel readable; full table is in the TSV

    fig_w = max(10.0, 0.42 * len(taxa) + 3)
    fig_h = max(4.5, 0.22 * len(diag) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    ax.set_xlim(-0.5, len(taxa) - 0.5)
    ax.set_ylim(-1.9, len(diag) - 0.5 + 1.6)
    ax.set_axis_off()

    for j, t in enumerate(taxa):
        g = meta.get(t, {}).get("group", "NA")
        ax.add_patch(Rectangle((j - 0.45, len(diag) + 0.85), 0.9, 0.55,
                                   facecolor=GROUP_COLOR.get(g, "#8c877d"), edgecolor="none"))

    for i, pos in enumerate(diag):
        y = len(diag) - 1 - i
        for j, t in enumerate(taxa):
            c = seqs[t][pos - 1]
            col = AA.get(CLASS.get(c.upper(), "special"), AA["special"]) if c not in "-Xx" else AA["gap"]
            ax.add_patch(Rectangle((j - 0.45, y - 0.38), 0.9, 0.76,
                                       facecolor=col, edgecolor="white", linewidth=0.4))
            ax.text(j, y, c, ha="center", va="center", fontsize=6.2, color="white"
                    if c not in "-Xx" and CLASS.get(c.upper(), "special") in
                    ("positive", "negative", "special") else "#2a2723")
        ax.text(-0.9, y, str(pos), ha="right", va="center", fontsize=7, color="#2a2723")

    # conservation track
    for i, pos in enumerate(diag):
        y = len(diag) - 1 - i
        v = cons.get(pos, 0.0)
        ax.add_patch(Rectangle((len(taxa) + 0.2, y - 0.38), 2.0 * v, 0.76,
                                   facecolor="#eda100", edgecolor="none"))
    ax.text(len(taxa) + 0.2, len(diag) + 0.2, "conservation", fontsize=7.5, color="#8c877d")

    ax.set_title(f"{a.family}: diagnostic amino-acid residues "
                 f"({len(diag)} columns, full table in diagnostic_residues.tsv)",
                 loc="left", fontweight="bold", fontsize=11)
    handles = [Patch(facecolor=GROUP_COLOR[k], label=k) for k in GROUP_COLOR]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8, ncol=2)
    for ext, path in [("pdf", a.out_pdf), ("svg", a.out_svg), ("png", a.out_png)]:
        fig.savefig(path, bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)
    print(f"plotted {len(diag)} diagnostic columns for {a.family}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
