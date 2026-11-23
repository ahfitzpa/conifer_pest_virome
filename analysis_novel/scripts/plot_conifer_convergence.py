#!/usr/bin/env python3
"""Cross-host & cross-family codon/dinucleotide convergence plot (Conifer Pest Virome).

Reads all `results/signatures/*/all_orf_composition.tsv` (or `dinucleotide_odds.tsv`)
and `results/phylogeny/*/selected.tsv` to plot:
  A  Wright's ENc vs. GC3s with theoretical neutral curve, coloured by insect host
  B  Cross-family CpG vs. UpA odds ratios, coloured by insect host
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK, INK_SOFT, GRID = "#101214", "#505258", "#E2E8F0"

HOST_COLOURS = {
    "Ips confusus": "#1868DB",
    "Hylobius abietis": "#16794F",
    "Locusta migratoria": "#C75300",
    "ICTV Reference": "#64748B",
    "Other": "#8F4FDB",
}


def norm_host(h: str, origin: str = "") -> str:
    if origin == "reference":
        return "ICTV Reference"
    s = h.strip().lower().replace("_", " ")
    if "ips" in s and "confusus" in s:
        return "Ips confusus"
    if "hylobius" in s and "abietis" in s:
        return "Hylobius abietis"
    if "locusta" in s and "migratoria" in s:
        return "Locusta migratoria"
    return "Other"


def load_host_lookup(phylo_dir: Path) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for tsv in sorted(phylo_dir.glob("*/selected.tsv")):
        with tsv.open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                lbl = (r.get("label") or "").strip()
                if not lbl:
                    continue
                grp = norm_host(r.get("host", ""), r.get("origin", ""))
                lookup[lbl] = grp
                lookup[lbl.split("_")[0]] = grp
    return lookup


def as_float(v: str | None) -> float:
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signatures-dir", default="results/signatures", type=Path)
    ap.add_argument("--phylogeny-dir", default="results/phylogeny", type=Path)
    ap.add_argument("--fig-dir", default="figures/overview", type=Path)
    args = ap.parse_args()

    host_map = load_host_lookup(args.phylogeny_dir)
    rows: list[dict] = []

    for tsv in sorted(args.signatures_dir.glob("*/all_orf_composition.tsv")):
        fam = tsv.parent.name
        with tsv.open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = (r.get("genome") or "").strip()
                grp = host_map.get(g) or host_map.get(g.split("_")[0]) or norm_host(g)
                rows.append({
                    "family": fam,
                    "host": grp,
                    "gc3s": as_float(r.get("GC3s")),
                    "enc": as_float(r.get("ENc")),
                    "cpg": as_float(r.get("CpG")),
                    "upa": as_float(r.get("UpA")),
                })

    if not rows:
        for tsv in sorted(args.signatures_dir.glob("*/dinucleotide_odds.tsv")):
            fam = tsv.parent.name
            with tsv.open(newline="", encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    g = (r.get("genome") or "").strip()
                    grp = host_map.get(g) or host_map.get(g.split("_")[0]) or norm_host(g)
                    rows.append({
                        "family": fam,
                        "host": grp,
                        "gc3s": float("nan"),
                        "enc": float("nan"),
                        "cpg": as_float(r.get("CG")),
                        "upa": as_float(r.get("TA")),
                    })

    if not rows:
        raise SystemExit(f"No composition tables found under {args.signatures_dir}")

    fig, (ax_enc, ax_di) = plt.subplots(1, 2, figsize=(12.4, 5.2), constrained_layout=True)

    # Panel A: Wright's ENc vs GC3s
    xs_null = np.linspace(0.05, 0.95, 200)
    enc_null = 2.0 + xs_null + 29.0 / (xs_null ** 2 + (1.0 - xs_null) ** 2)
    ax_enc.plot(xs_null, np.clip(enc_null, 20, 61), color=INK, lw=1.4, ls="--", label="Wright neutral curve", zorder=1)
    for host, col in HOST_COLOURS.items():
        pts = [r for r in rows if r["host"] == host and np.isfinite(r["gc3s"]) and np.isfinite(r["enc"])]
        if not pts:
            continue
        ax_enc.scatter([p["gc3s"] for p in pts], [p["enc"] for p in pts],
                       s=20, c=col, alpha=0.75, edgecolors="white", linewidths=0.3,
                       label=f"{host} (n={len(pts)})", zorder=3)
    ax_enc.set_xlabel("GC3s", fontsize=9, color=INK)
    ax_enc.set_ylabel("Effective Number of Codons (ENc)", fontsize=9, color=INK)
    ax_enc.set_title("A  Codon usage vs. GC3s across all mapped ORFs", loc="left", fontsize=10, fontweight="bold", color=INK)
    ax_enc.set_ylim(20, 63)
    ax_enc.grid(color=GRID, lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax_enc.spines[sp].set_visible(False)
    ax_enc.legend(frameon=False, fontsize=7.5, loc="lower left")

    # Panel B: CpG vs UpA across hosts
    for host, col in HOST_COLOURS.items():
        pts = [r for r in rows if r["host"] == host and np.isfinite(r["cpg"]) and np.isfinite(r["upa"])]
        if not pts:
            continue
        ax_di.scatter([p["cpg"] for p in pts], [p["upa"] for p in pts],
                      s=20, c=col, alpha=0.75, edgecolors="white", linewidths=0.3,
                      label=f"{host} (n={len(pts)})", zorder=3)
    ax_di.axvline(1.0, color=INK_SOFT, lw=0.9, ls="--", zorder=1)
    ax_di.axhline(1.0, color=INK_SOFT, lw=0.9, ls="--", zorder=1)
    ax_di.set_xlabel("CpG odds ratio (1.0 = independence)", fontsize=9, color=INK)
    ax_di.set_ylabel("UpA odds ratio (1.0 = independence)", fontsize=9, color=INK)
    ax_di.set_title("B  Cross-family CpG vs. UpA host convergence", loc="left", fontsize=10, fontweight="bold", color=INK)
    ax_di.grid(color=GRID, lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax_di.spines[sp].set_visible(False)
    ax_di.legend(frameon=False, fontsize=7.5, loc="best")

    fig.suptitle("Cross-host codon & dinucleotide convergence across viral families",
                 fontsize=11.5, fontweight="bold", color=INK, ha="left", x=0.01)

    args.fig_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(args.fig_dir / f"conifer_host_convergence.{ext}",
                    bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.fig_dir}/conifer_host_convergence.(pdf|svg|png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
