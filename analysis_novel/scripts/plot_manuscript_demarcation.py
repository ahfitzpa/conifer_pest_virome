#!/usr/bin/env python3
"""Render Manuscript Fig. 2B: Formal ICTV Demarcation — Nearest ICTV Reference
Identity (%) vs. CheckV Completeness (%) coloured by insect host (title-free).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INK = "#101214"
MUTED = "#505258"
GRID = "#E2E8F0"

HOST_COL = {
    "Locusta migratoria": "#C75300",
    "Hylobius abietis": "#16794F",
    "Ips confusus": "#1868DB",
}
HOST_SHORT = {
    "Locusta migratoria": "L. migratoria",
    "Hylobius abietis": "H. abietis",
    "Ips confusus": "I. confusus",
}
FOCAL = {"Iflaviridae", "Polycipiviridae", "Amarillovirales", "Dicistroviridae", "Partitiviridae", "Phenuiviridae"}


def norm_host(h: str) -> str:
    s = h.strip().lower().replace("_", " ")
    if "locusta" in s:
        return "Locusta migratoria"
    if "hylobius" in s:
        return "Hylobius abietis"
    if "ips" in s:
        return "Ips confusus"
    return "Other"


def as_float(v: str | None) -> float | None:
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--demarcation-summary", default=Path("results/demarcation/all_families_demarcation_summary.tsv"), type=Path)
    ap.add_argument("--out-dir", default=Path("figures/manuscript"), type=Path)
    args = ap.parse_args()

    pts: list[dict[str, object]] = []
    with args.demarcation_summary.open(encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ident = as_float(r.get("nearest_identity_pct"))
            comp = as_float(r.get("checkv_completeness"))
            if ident is None:
                continue
            pts.append({
                "family": r.get("family", ""),
                "candidate": r.get("candidate", ""),
                "host": norm_host(r.get("host", "")),
                "ident": ident,
                "comp": comp if comp is not None else 0.0,
            })

    if not pts:
        raise SystemExit(f"No rows found in {args.demarcation_summary}")

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    for h, col in HOST_COL.items():
        sub = [p for p in pts if p["host"] == h]
        if not sub:
            continue
        ax.scatter(
            [p["comp"] for p in sub],
            [p["ident"] for p in sub],
            c=col, s=56, alpha=0.86, edgecolors="white", linewidths=0.7, zorder=3,
            label=f"{HOST_SHORT[h]} (n = {len(sub)})",
        )

    # Demarcation threshold guides (placed cleanly inside top-left so right margin stays uncluttered)
    for thresh, lbl in ((90, "90% (species-level guide)"), (70, "70% (genus-level guide)"), (50, "50% (deep divergence)")):
        ax.axhline(thresh, color=MUTED, ls="--", lw=0.9, alpha=0.65, zorder=1)
        ax.text(1.5, thresh + 1.4, lbl, va="bottom", ha="left", fontsize=7.2, color=MUTED)

    # Staggered, non-overlapping callouts for the top complete/near-complete focal families (one per family)
    offsets = {
        "Partitiviridae": (-62, 18),
        "Dicistroviridae": (-76, 4),
        "Polycipiviridae": (-96, -10),
        "Iflaviridae": (-50, -32),
        "Phenuiviridae": (-22, 18),
        "Aliusviridae": (-52, 16),
    }
    best_by_fam: dict[str, dict[str, object]] = {}
    for p in pts:
        fam = str(p["family"])
        if fam not in offsets or float(p["comp"]) < 55.0:
            continue
        if fam not in best_by_fam or float(p["comp"]) > float(best_by_fam[fam]["comp"]):
            best_by_fam[fam] = p

    for fam, p in best_by_fam.items():
        dx, dy = offsets[fam]
        ax.annotate(
            fam,
            xy=(float(p["comp"]), float(p["ident"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.6, fontstyle="italic", fontweight="bold", color=INK,
            arrowprops={"arrowstyle": "-", "color": "#64748B", "lw": 0.8, "shrinkA": 2, "shrinkB": 4},
            zorder=5,
        )

    ax.set_xlim(-3, 106)
    ax.set_ylim(0, 103)
    ax.set_xlabel("CheckV genome completeness (%)", fontsize=9.2, color=INK)
    ax.set_ylabel("Nearest ICTV reference identity (%)", fontsize=9.2, color=INK)
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    ax.legend(loc="lower right", frameon=True, facecolor="#F8FAFC", edgecolor="#DDDEE1", fontsize=8.0)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(args.out_dir / f"fig2b_demarcation_novelty.{ext}",
                    bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.out_dir}/fig2b_demarcation_novelty.(pdf|svg|png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
