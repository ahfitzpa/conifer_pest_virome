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
    ap.add_argument("--max-candidates", type=int, default=None,
                    help="Max representative candidates to display (default from config or 25)")
    ap.add_argument("--max-sites", type=int, default=None,
                    help="Max diagnostic columns to display (default from config or 45)")
    a = ap.parse_args()

    cfg: dict = {}
    try:
        import yaml  # type: ignore
        cfg = yaml.safe_load(open(a.config, encoding="utf-8")).get("phylogenetics", {}) or {}
    except ImportError:
        import re as _re
        for line in open(a.config, encoding="utf-8"):
            m = _re.match(r"^\s*(diagnostic_\w+)\s*:\s*([^\s#]+)", line)
            if m:
                cfg[m.group(1)] = m.group(2).strip("\"'")
    group_col = cfg.get("diagnostic_group_column", "group")
    max_candidates = (a.max_candidates if a.max_candidates is not None
                      else int(cfg.get("diagnostic_max_candidates", 25)))
    max_sites = (a.max_sites if a.max_sites is not None
                 else int(cfg.get("diagnostic_max_sites", 45)))

    matrix = read_tsv(a.matrix)
    cons_rows = read_tsv(a.conservation)
    cons = {int(r["position"]): float(r["conservation"]) for r in cons_rows}
    seqs = dict(iter_fasta(a.alignment))
    meta = {r["taxon_id"]: r for r in read_tsv(a.metadata)}

    diag_rows = [r for r in matrix if r.get("diagnostic") == "1"]
    diag = [int(r["position"]) for r in diag_rows][:max_sites]

    # Partition taxa into ICTV references and candidates; subset candidates by alignment coverage
    ictv_taxa = sorted(t for t in seqs if meta.get(t, {}).get(group_col) == "ICTV")
    cand_taxa = [t for t in seqs if meta.get(t, {}).get(group_col) != "ICTV"]
    n_total_cand = len(cand_taxa)
    if len(cand_taxa) > max_candidates:
        cand_taxa.sort(
            key=lambda t: (sum(1 for c in seqs[t] if c not in "-Xx"), t),
            reverse=True,
        )
        cand_taxa = sorted(cand_taxa[:max_candidates])
    else:
        cand_taxa = sorted(cand_taxa)
    shown_taxa = ictv_taxa + cand_taxa

    if not diag or not shown_taxa:
        fig, ax = plt.subplots(figsize=(10.5, 2.8), constrained_layout=True)
        ax.set_axis_off()
        ax.text(0.5, 0.5,
                f"{a.family}: 0 diagnostic columns passed coverage & consensus filters "
                f"(taxa={len(seqs)}; see diagnostic_residues.tsv)",
                ha="center", va="center", fontsize=10, color="#2a2723")
        for ext, path in [("pdf", a.out_pdf), ("svg", a.out_svg), ("png", a.out_png)]:
            fig.savefig(path, bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
        plt.close(fig)
        print(f"plotted 0 diagnostic columns for {a.family}", file=sys.stderr)
        return 0

    # Transposed layout (project_maureen style): taxa as rows, diagnostic sites as columns
    fig_w = min(18.0, max(10.0, 0.28 * len(diag) + 4.8))
    fig_h = max(4.5, 0.24 * len(shown_taxa) + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    ax.set_xlim(-0.5, len(diag) - 0.5)
    ax.set_ylim(-1.6, len(shown_taxa) + 2.2)
    ax.set_axis_off()

    # Top track: per-site conservation bars
    for j, pos in enumerate(diag):
        v = cons.get(pos, 0.0)
        ax.add_patch(Rectangle((j - 0.40, len(shown_taxa) + 0.35), 0.80, 1.15 * v,
                               facecolor="#eda100", edgecolor="none"))
        ax.text(j, len(shown_taxa) + 0.05, str(pos), ha="center", va="bottom",
                fontsize=6.2, rotation=90, color="#2a2723")
    ax.text(-0.8, len(shown_taxa) + 0.9, "cons.", ha="right", va="center",
            fontsize=7.2, color="#8c877d", fontweight="bold")

    # Rows: ICTV at top, candidates below
    for i, t in enumerate(shown_taxa):
        y = len(shown_taxa) - 1 - i
        g = meta.get(t, {}).get(group_col, "NA")
        gcol = GROUP_COLOR.get(g, "#8c877d")
        ax.add_patch(Rectangle((-0.95, y - 0.38), 0.35, 0.76, facecolor=gcol, edgecolor="none"))
        label = t if len(t) <= 32 else t[:29] + "..."
        ax.text(-1.15, y, label, ha="right", va="center", fontsize=6.8, color="#2a2723")

        for j, pos in enumerate(diag):
            c = seqs[t][pos - 1]
            col = AA.get(CLASS.get(c.upper(), "special"), AA["special"]) if c not in "-Xx" else AA["gap"]
            ax.add_patch(Rectangle((j - 0.45, y - 0.38), 0.90, 0.76,
                                   facecolor=col, edgecolor="white", linewidth=0.35))
            ax.text(j, y, c, ha="center", va="center", fontsize=6.0,
                    color="white" if c not in "-Xx" and CLASS.get(c.upper(), "special") in
                    ("positive", "negative", "special") else "#2a2723")

    if ictv_taxa and cand_taxa:
        sep_y = len(cand_taxa) - 0.5
        ax.plot([-1.1, len(diag) - 0.5], [sep_y, sep_y], color="#2a2723", lw=1.0, ls="--")

    subset_note = (f"showing {len(cand_taxa)}/{n_total_cand} candidates + {len(ictv_taxa)} ICTV"
                   if n_total_cand > len(cand_taxa)
                   else f"{len(cand_taxa)} candidates + {len(ictv_taxa)} ICTV")
    ax.set_title(f"{a.family}: diagnostic residues ({len(diag)} sites; {subset_note})",
                 loc="left", fontweight="bold", fontsize=11)
    handles = [Patch(facecolor=GROUP_COLOR[k], label=k) for k in GROUP_COLOR]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8, ncol=2)
    for ext, path in [("pdf", a.out_pdf), ("svg", a.out_svg), ("png", a.out_png)]:
        fig.savefig(path, bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"plotted {len(diag)} diagnostic columns x {len(shown_taxa)} taxa for {a.family}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
