#!/usr/bin/env python3
"""Per-family figures for selection and composition signatures.

Figure conventions follow project_maureen/scripts/06_plot_qc_and_diagnostics.py,
which the repo notes were themselves adapted from conifer_pest_virome: the
Okabe-Ito colourblind-safe palette, a save_all helper writing PDF, SVG and PNG
together, Agg backend, white facecolour and 300 dpi for the raster.

Three panels, each answering one question:
  A  where selection sits along the ORF -- omega per site on a log axis, points
     coloured by category, so a handful of positively selected sites is visible
     against the purifying background rather than averaged into it
  B  how many sites fall in each category -- a horizontal bar with the count
     labelled directly, because four categories do not need an axis to be read
  C  whether composition shows host adaptation -- dinucleotide odds ratios per
     genome against a reference line at 1.0, with GpC included as the control:
     CpG and UpA depleted while GpC sits near 1 is the expected signal, and all
     three moving together would instead suggest a compositional artefact

No dual axes anywhere, one measure per panel.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Okabe-Ito, as used across project_maureen and conifer. Validated
# colourblind-safe: worst adjacent pair dE 11.0 deutan, 15.6 normal vision.
CAT_COLOURS = {
    "pervasive_purifying": "#0072B2",
    "episodic_positive": "#E69F00",
    "pervasive_positive": "#D55E00",
    "neutral_or_unresolved": "#B8B2A6",
}
CAT_ORDER = ["pervasive_purifying", "neutral_or_unresolved",
             "episodic_positive", "pervasive_positive"]
CAT_LABEL = {
    "pervasive_purifying": "pervasive purifying",
    "neutral_or_unresolved": "neutral / unresolved",
    "episodic_positive": "episodic positive",
    "pervasive_positive": "pervasive positive",
}
INK, INK_SOFT, GRID = "#1F1F1F", "#4A4A4A", "#E7E3DA"
DINUC_COLOUR = {"CG": "#D55E00", "TA": "#E69F00", "GC": "#009E73"}
DINUC_LABEL = {"CG": "CpG", "TA": "UpA", "GC": "GpC (control)"}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def save_all(fig, fig_dir: Path, stem: str) -> None:
    # 1: one figure, three formats, as the other projects' figures do
    fig_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(fig_dir / f"{stem}.{ext}", bbox_inches="tight",
                    dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    log(f"[plot] saved {stem} (.pdf, .svg, .png)")


def read_tsv(path: str) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def as_float(v, default=float("nan")) -> float:
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return default


def panel_sites(ax, sites: list[dict]) -> int:
    # 2: omega per site; log axis, because purifying sites crowd below 1
    drawn = 0
    for cat in CAT_ORDER:
        xs = [int(r["site"]) for r in sites
              if r.get("selection_category") == cat and as_float(r.get("omega")) > 0]
        ys = [as_float(r["omega"]) for r in sites
              if r.get("selection_category") == cat and as_float(r.get("omega")) > 0]
        if not xs:
            continue
        # positives are the point of the panel, so they sit on top and larger
        emphasis = cat in ("episodic_positive", "pervasive_positive")
        ax.scatter(xs, ys, s=26 if emphasis else 9,
                   c=CAT_COLOURS[cat], label=CAT_LABEL[cat],
                   edgecolors="white" if emphasis else "none",
                   linewidths=0.6 if emphasis else 0, zorder=3 if emphasis else 2,
                   alpha=1.0 if emphasis else 0.65)
        drawn += len(xs)
    ax.axhline(1.0, color=INK_SOFT, lw=1.0, ls="--", zorder=1)
    ax.text(0.995, 1.0, "  neutral (omega = 1)", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=7.4, color=INK_SOFT)
    ax.set_yscale("log")
    ax.set_xlabel("codon position in ORF1", fontsize=8.5, color=INK)
    ax.set_ylabel("omega (dN/dS)", fontsize=8.5, color=INK)
    ax.set_title("A  Site-level selection along the ORF", loc="left",
                 fontsize=9.5, fontweight="bold", color=INK)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.0, colors=INK_SOFT)
    if drawn:
        ax.legend(frameon=False, fontsize=7.6, labelcolor=INK_SOFT,
                  loc="upper left", handletextpad=0.4, borderaxespad=0.2)
    return drawn


def panel_counts(ax, sites: list[dict]) -> None:
    # 3: counts per category, labelled directly rather than read off an axis
    counts = {cat: sum(1 for r in sites if r.get("selection_category") == cat)
              for cat in CAT_ORDER}
    ys = np.arange(len(CAT_ORDER))
    widths = [counts[c] for c in CAT_ORDER]
    ax.barh(ys, widths, height=0.6,
            color=[CAT_COLOURS[c] for c in CAT_ORDER], zorder=2)
    total = sum(widths) or 1
    for y, w in zip(ys, widths):
        ax.text(w + total * 0.015, y, f"{w:,}  ({w / total:.0%})",
                va="center", fontsize=7.8, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels([CAT_LABEL[c] for c in CAT_ORDER], fontsize=8.2,
                       color=INK)
    ax.set_xlabel("sites", fontsize=8.5, color=INK)
    ax.set_xlim(0, total * 1.28)
    ax.set_title("B  Sites per selection category", loc="left",
                 fontsize=9.5, fontweight="bold", color=INK)
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.0, colors=INK_SOFT, left=False)


def panel_composition(ax, dinuc: list[dict]) -> None:
    # 4: odds ratios per genome, with GpC as the control for an artefact
    keys = [k for k in ("CG", "TA", "GC") if dinuc and k in dinuc[0]]
    if not keys:
        ax.text(0.5, 0.5, "no dinucleotide table", ha="center", va="center",
                transform=ax.transAxes, fontsize=8.5, color=INK_SOFT)
        ax.set_axis_off()
        return
    rng = np.random.default_rng(0)  # fixed, so the jitter is reproducible
    for i, key in enumerate(keys):
        vals = [as_float(r[key]) for r in dinuc if as_float(r.get(key)) > 0]
        if not vals:
            continue
        jitter = rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=14,
                   c=DINUC_COLOUR[key], alpha=0.55, edgecolors="none", zorder=2)
        med = float(np.median(vals))
        ax.plot([i - 0.28, i + 0.28], [med, med], color=INK, lw=1.8, zorder=4)
        ax.text(i + 0.33, med, f"median {med:.2f}", va="center", fontsize=7.6,
                color=INK)
    ax.axhline(1.0, color=INK_SOFT, lw=1.0, ls="--", zorder=1)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([DINUC_LABEL[k] for k in keys], fontsize=8.2, color=INK)
    ax.set_ylabel("odds ratio (1.0 = independence)", fontsize=8.5, color=INK)
    ax.set_title("C  Dinucleotide composition per genome", loc="left",
                 fontsize=9.5, fontweight="bold", color=INK)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.0, colors=INK_SOFT)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", required=True, help="site_selection.tsv")
    p.add_argument("--summary", help="selection_summary.tsv")
    p.add_argument("--dinucleotide", help="dinucleotide_odds.tsv")
    p.add_argument("--rscu", help="rscu.tsv (row count only, for the caption)")
    p.add_argument("--family", required=True)
    p.add_argument("--fig-dir", required=True)
    p.add_argument("--stem", help="default signatures_<family>")
    args = p.parse_args()

    sites = read_tsv(args.sites)
    dinuc = read_tsv(args.dinucleotide)
    rscu = read_tsv(args.rscu) if args.rscu else []
    if not sites:
        log(f"[plot] {args.sites} has no rows; nothing to draw")
        return 1
    log(f"[plot] {len(sites)} site(s), {len(dinuc)} genome(s) with composition")

    fig = plt.figure(figsize=(13.6, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], width_ratios=[1.0, 1.0])
    panel_sites(fig.add_subplot(gs[0, :]), sites)
    panel_counts(fig.add_subplot(gs[1, 0]), sites)
    panel_composition(fig.add_subplot(gs[1, 1]), dinuc)

    median_omega = ""
    for r in read_tsv(args.summary) if args.summary else []:
        if r.get("metric") == "median_omega":
            median_omega = f", median omega {as_float(r['value']):.3f}"
    fig.suptitle(f"{args.family}: selection and composition signatures",
                 fontsize=11.0, fontweight="bold", color=INK, ha="left", x=0.01)
    fig.text(0.01, -0.012,
             f"{len(sites):,} codon sites{median_omega}; "
             f"{len(dinuc) or len(rscu)} genome(s). "
             f"Counts and thresholds in selection_summary.tsv.",
             fontsize=7.6, color=INK_SOFT, ha="left")

    save_all(fig, Path(args.fig_dir), args.stem or f"signatures_{args.family}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
