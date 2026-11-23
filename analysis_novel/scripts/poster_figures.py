#!/usr/bin/env python3
"""Poster figures: per-family trees coloured by host, plus a host x family map.

Reads what the shared phylogeny already wrote -- results/phylogeny/<family>/
tree.newick and selected.tsv -- so it needs no re-analysis and can run while
the workflow is still going. Families with no finished tree are skipped and
listed, rather than failing the run.

  --out-dir/<Family>_tree.{png,svg}   one tree per family, candidate tips
                                      coloured by host, ICTV references grey
  --out-dir/host_family_matrix.{png,svg}  candidate contigs per host x family
  --out-dir/figure_manifest.tsv       what was drawn, with taxon counts

Usage:
  python poster_figures.py --results results/phylogeny --out-dir figures/poster
  python poster_figures.py --results results/phylogeny --out-dir figures/poster \
      --families Dicistroviridae Iflaviridae Rhabdoviridae
"""
from __future__ import annotations

import argparse
import csv
import os
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from Bio import Phylo  # noqa: E402

# Colour-blind safe (Okabe-Ito); references stay neutral grey so candidates pop.
HOST_COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                "#56B4E9", "#F0E442"]
REF_COLOUR = "#9A9A9A"
DPI = 300


def read_selected(path):
    """label -> (origin, host); tolerant of a partially written file."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab = (r.get("label") or "").strip()
            if lab:
                out[lab] = ((r.get("origin") or "").strip(),
                            (r.get("host") or "").strip() or "unknown")
    return out


def tidy_host(h):
    return h.replace("_", " ") if h else "unknown"


def draw_tree(fam, tree_path, selected, out_dir, host_colour):
    tree = Phylo.read(tree_path, "newick")
    tips = tree.get_terminals()
    n = len(tips)
    if n < 3:
        return None

    # 1: label and colour every tip from selected.tsv
    n_cand = 0
    for t in tips:
        name = (t.name or "").strip().strip("'\"")
        origin, host = selected.get(name, ("reference", ""))
        if origin == "candidate":
            n_cand += 1
            t.color = host_colour.get(host, "#333333")
            run = name.split("_n")[0]
            t.name = f"{run}  [{tidy_host(host)}]"
        else:
            t.color = REF_COLOUR
            t.name = name.replace("_", " ")[:52]

    # 2: internal nodes carry IQ-TREE's "aLRT/UFboot"; show only strong support
    support = {}
    for cl in tree.get_nonterminals():
        lab = str(cl.confidence) if cl.confidence is not None else (cl.name or "")
        cl.name = None
        try:
            boot = float(str(lab).split("/")[-1])
            if boot >= 90:
                support[id(cl)] = f"{boot:.0f}"
        except (ValueError, IndexError):
            pass

    height = max(4.0, 0.17 * n)
    fig, ax = plt.subplots(figsize=(11, height))
    Phylo.draw(tree, axes=ax, do_show=False,
               show_confidence=False,
               label_colors=lambda s: "#222222",
               branch_labels=lambda c: support.get(id(c), ""))
    ax.set_title(f"{fam} — {n_cand} candidate contigs vs {n - n_cand} "
                 f"ICTV reference genomes", fontsize=13, loc="left", pad=12)
    ax.set_xlabel("substitutions per site", fontsize=9)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_ylabel("")
    ax.set_yticks([])

    present = [h for h in host_colour if any(
        sel_o == "candidate" and sel_h == h for sel_o, sel_h in selected.values())]
    handles = [Line2D([], [], color=host_colour[h], lw=3, label=tidy_host(h))
               for h in present]
    handles.append(Line2D([], [], color=REF_COLOUR, lw=3, label="ICTV reference"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, fontsize=9, title="candidate host",
              title_fontsize=9)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out_dir, f"{fam}_tree.{ext}"), dpi=DPI,
                    bbox_inches="tight")
    plt.close(fig)
    return n, n_cand


def draw_matrix(counts, hosts, fams, out_dir):
    if not counts:
        return
    fig, ax = plt.subplots(figsize=(max(7, 0.5 * len(fams) + 3),
                                    max(3, 0.5 * len(hosts) + 2)))
    grid = [[counts.get((h, f), 0) for f in fams] for h in hosts]
    im = ax.imshow(grid, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(range(len(fams)))
    ax.set_xticklabels(fams, rotation=55, ha="right", fontsize=9)
    ax.set_yticks(range(len(hosts)))
    ax.set_yticklabels([tidy_host(h) for h in hosts], fontsize=10)
    for i, h in enumerate(hosts):
        for j, f in enumerate(fams):
            v = counts.get((h, f), 0)
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=8,
                        color="#222222" if v < max(1, max(map(max, grid))) * 0.6
                        else "white")
    ax.set_title("Candidate contigs per host and virus family", fontsize=12,
                 loc="left", pad=10)
    fig.colorbar(im, ax=ax, shrink=0.8, label="candidate contigs")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out_dir, f"host_family_matrix.{ext}"), dpi=DPI,
                    bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results/phylogeny")
    ap.add_argument("--out-dir", default="figures/poster")
    ap.add_argument("--families", nargs="*", default=None)
    ap.add_argument("--min-candidates", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    fams = sorted(a.families or [d for d in os.listdir(a.results)
                                 if os.path.isdir(os.path.join(a.results, d))])

    # one colour per host, shared across every tree so the poster reads as a set
    host_counts = Counter()
    per_cell = Counter()
    sel_cache = {}
    for fam in fams:
        sel = read_selected(os.path.join(a.results, fam, "selected.tsv"))
        sel_cache[fam] = sel
        for origin, host in sel.values():
            if origin == "candidate":
                host_counts[host] += 1
                per_cell[(host, fam)] += 1
    hosts = [h for h, _ in host_counts.most_common()]
    host_colour = {h: HOST_COLOURS[i % len(HOST_COLOURS)]
                   for i, h in enumerate(hosts)}

    drawn, skipped = [], []
    for fam in fams:
        tree_path = os.path.join(a.results, fam, "tree.newick")
        if not (os.path.exists(tree_path) and os.path.getsize(tree_path) > 0):
            skipped.append((fam, "no tree yet"))
            continue
        sel = sel_cache[fam]
        n_cand = sum(1 for o, _ in sel.values() if o == "candidate")
        if n_cand < a.min_candidates:
            skipped.append((fam, f"{n_cand} candidates"))
            continue
        try:
            got = draw_tree(fam, tree_path, sel, a.out_dir, host_colour)
        except Exception as exc:  # a tree still being written
            skipped.append((fam, f"{type(exc).__name__}: {exc}"))
            continue
        if got is None:
            skipped.append((fam, "fewer than 3 tips"))
            continue
        drawn.append((fam, got[0], got[1]))
        print(f"[fig] {fam:22} {got[0]:4} taxa ({got[1]} candidate)")

    fams_drawn = [f for f, _, _ in drawn]
    draw_matrix(per_cell, hosts, fams_drawn, a.out_dir)

    with open(os.path.join(a.out_dir, "figure_manifest.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "taxa", "candidates", "figure"])
        for fam, n, nc in drawn:
            w.writerow([fam, n, nc, f"{fam}_tree.png"])
    print(f"\n[fig] drew {len(drawn)} tree(s) -> {a.out_dir}")
    if skipped:
        print(f"[fig] skipped {len(skipped)}: "
              + ", ".join(f"{f} ({why})" for f, why in skipped[:12]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
