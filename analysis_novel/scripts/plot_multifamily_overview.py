#!/usr/bin/env python3
"""Multi-family overview plots for poster & manuscript reporting.

Reads all `results/phylogeny/*/selected.tsv` tables (available as soon as
`select_contigs` / `phylogeny_all` runs) and, when present,
`results/demarcation/*/pairwise_identity.tsv` (from `demarcation_all`),
producing two Atlassian-themed figures (PDF, SVG, 300-dpi PNG):

  1. host_family_matrix      Host/BioProject x Viral Family bubble matrix
                             (size = contig count, colour = median CheckV
                             completeness) + right-hand contig length range
                             strip (min-median-max in kb).
  2. novelty_vs_completeness Nearest-ICTV AA % identity vs CheckV completeness
                             across all families, coloured by host and sized
                             by contig length (written when pairwise_identity.tsv
                             files exist).
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

# ggsci::pal_atlassian("categorical8") + Atlassian Design System tokens
ATLASSIAN_HOSTS = [
    "#1868DB",  # Atlassian Blue
    "#C75300",  # chart-orange
    "#22A06B",  # chart-green
    "#803FA5",  # chart-purple-bold
    "#2898BD",  # chart-teal
    "#AE2E24",  # chart-red
    "#1C2B42",  # brand-boldest navy
    "#505258",  # neutral ink
]
INK = "#101214"
INK_SOFT = "#505258"
GRID = "#DDDEE1"
COMP_CMAP = LinearSegmentedColormap.from_list(
    "atlassian_comp", ["#DDDEE1", "#2898BD", "#1868DB", "#1C2B42"]
)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", errors="ignore", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def as_float(v, default=None):
    try:
        s = str(v).strip()
        if not s or s.lower() in ("na", "nan", "none"):
            return default
        return float(s)
    except (TypeError, ValueError):
        return default


def save_all(fig, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(
            out_dir / f"{stem}.{ext}",
            bbox_inches="tight",
            dpi=300 if ext == "png" else None,
            facecolor="white",
        )
    plt.close(fig)
    print(f"[overview] saved {out_dir / stem}.{{pdf,svg,png}}", file=sys.stderr)


def load_candidates(phylogeny_dir: Path, group_col: str) -> list[dict]:
    rows = []
    for tsv in sorted(phylogeny_dir.glob("*/selected.tsv")):
        fam_dir = tsv.parent.name
        for r in read_tsv(tsv):
            origin = (r.get("origin") or "candidate").strip().lower()
            if origin != "candidate":
                continue
            fam = (r.get("consensus_family") or fam_dir).strip()
            host = (r.get(group_col) or "").strip()
            if not host:
                host = (r.get("bioproject") or r.get("run_id") or "Unknown").strip()
            length = as_float(r.get("length"), 0.0) or 0.0
            comp = as_float(r.get("checkv_completeness"), None)
            rows.append({
                "label": (r.get("label") or "").strip(),
                "family": fam,
                "host": host,
                "bioproject": (r.get("bioproject") or "").strip(),
                "run_id": (r.get("run_id") or "").strip(),
                "length": length,
                "completeness": comp,
                "agreement": (r.get("family_agreement") or "").strip(),
            })
    return rows


def load_nearest_ictv(demarcation_dir: Path) -> dict[str, dict]:
    """Return max AA % identity to an ICTV exemplar per candidate query."""
    best: dict[str, dict] = {}
    if not demarcation_dir.exists():
        return best
    for tsv in sorted(demarcation_dir.glob("*/pairwise_identity.tsv")):
        for r in read_tsv(tsv):
            if str(r.get("passes_50pct_coverage") or "0").strip() != "1":
                continue
            q = (r.get("query") or "").strip()
            ident = as_float(r.get("identity_pct"), None)
            if not q or ident is None:
                continue
            prev = best.get(q)
            if prev is None or ident > prev["identity_pct"]:
                best[q] = {
                    "identity_pct": ident,
                    "reference": (r.get("reference") or "").strip(),
                    "query_coverage_pct": as_float(r.get("query_coverage_pct"), 0.0) or 0.0,
                }
    return best


def plot_host_family_matrix(rows: list[dict], out_dir: Path, max_families: int) -> None:
    fam_counts = defaultdict(int)
    host_counts = defaultdict(int)
    cell_rows: dict[tuple[str, str], list[dict]] = defaultdict(list)
    fam_lengths: dict[str, list[float]] = defaultdict(list)

    for r in rows:
        fam_counts[r["family"]] += 1
        host_counts[r["host"]] += 1
        cell_rows[(r["family"], r["host"])].append(r)
        if r["length"] > 0:
            fam_lengths[r["family"]].append(r["length"] / 1000.0)

    fam_n_hosts = {
        f: sum(1 for h in host_counts if cell_rows.get((f, h)))
        for f in fam_counts
    }
    families = sorted(fam_counts, key=lambda f: (-fam_n_hosts[f], -fam_counts[f], f))[:max_families]
    families = list(reversed(families))  # multi-host & largest at top
    hosts = sorted(host_counts, key=lambda h: (-host_counts[h], h))

    fig_h = max(5.2, 0.28 * len(families) + 2.2)
    fig_w = max(8.8, 0.95 * len(hosts) + 6.0)
    fig, (ax_mat, ax_len) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [max(2.6, len(hosts) * 0.95), 2.1], "wspace": 0.14},
        sharey=True,
    )

    max_cnt = max((len(cell_rows[(f, h)]) for f in families for h in hosts), default=1)
    xs, ys, sizes, cols = [], [], [], []
    for yi, fam in enumerate(families):
        for xi, h in enumerate(hosts):
            grp = cell_rows.get((fam, h), [])
            if not grp:
                continue
            xs.append(xi)
            ys.append(yi)
            sizes.append(45 + 320 * (len(grp) / max_cnt))
            comps = [g["completeness"] for g in grp if g["completeness"] is not None]
            cols.append(statistics.median(comps) if comps else 0.0)

    sc = ax_mat.scatter(
        xs, ys, s=sizes, c=cols, cmap=COMP_CMAP, vmin=0, vmax=100,
        edgecolors=INK, linewidths=0.7, zorder=3,
    )
    for x, y, fam, h in [
        (xi, yi, fam, h)
        for yi, fam in enumerate(families)
        for xi, h in enumerate(hosts)
        if cell_rows.get((fam, h))
    ]:
        cnt = len(cell_rows[(fam, h)])
        comps = [g["completeness"] for g in cell_rows[(fam, h)] if g["completeness"] is not None]
        med_c = statistics.median(comps) if comps else 0.0
        txt_col = "white" if med_c >= 55 else INK
        ax_mat.text(x, y, str(cnt), ha="center", va="center",
                    fontsize=7.2, fontweight="bold", color=txt_col, zorder=4)

    ax_mat.set_xticks(range(len(hosts)))
    ax_mat.set_xticklabels([h.replace("_", " ") for h in hosts], rotation=28, ha="right", fontsize=8.5, fontstyle="italic", color=INK)
    ax_mat.set_xlabel("Host species", fontsize=8.8, color=INK)
    ax_mat.set_yticks(range(len(families)))
    ax_mat.set_yticklabels(families, fontsize=8.5, fontstyle="italic", color=INK)
    ax_mat.set_ylabel("Viral family", fontsize=8.8, color=INK)
    ax_mat.set_xlim(-0.6, len(hosts) - 0.4)
    ax_mat.set_ylim(-0.6, len(families) - 0.4)
    ax_mat.grid(True, color=GRID, linewidth=0.6, zorder=0)
    for sp in ("top", "right"):
        ax_mat.spines[sp].set_visible(False)
    ax_mat.set_title("A", loc="left", fontsize=10.5, fontweight="bold", color=INK)

    cbar = fig.colorbar(sc, ax=[ax_mat, ax_len], fraction=0.025, pad=0.02, aspect=32)
    cbar.set_label("Median CheckV completeness (%)", fontsize=8.2, color=INK_SOFT)
    cbar.ax.tick_params(labelsize=7.5, colors=INK_SOFT)

    # Right panel: contig length min–median–max (kb)
    for yi, fam in enumerate(families):
        ls = fam_lengths.get(fam, [])
        if not ls:
            continue
        lo, med, hi = min(ls), statistics.median(ls), max(ls)
        ax_len.plot([lo, hi], [yi, yi], color="#2898BD", lw=2.2, zorder=2,
                    label="Min–max length" if yi == 0 else None)
        ax_len.scatter([med], [yi], color="#1C2B42", s=26, zorder=3,
                       label="Median length" if yi == 0 else None)
    ax_len.set_xlabel("Contig length (kb)", fontsize=8.8, color=INK)
    ax_len.set_title("B", loc="left", fontsize=10.5, fontweight="bold", color=INK)
    ax_len.grid(axis="x", color=GRID, linewidth=0.6, zorder=0)
    for sp in ("top", "right"):
        ax_len.spines[sp].set_visible(False)
    ax_len.tick_params(labelsize=8, colors=INK_SOFT)
    ax_len.legend(frameon=False, fontsize=7.6, loc="lower right")

    save_all(fig, out_dir, "host_family_matrix")


def plot_novelty_vs_completeness(rows: list[dict], nearest: dict[str, dict], out_dir: Path) -> None:
    pts = []
    for r in rows:
        hit = nearest.get(r["label"])
        if not hit or r["completeness"] is None:
            continue
        pts.append({
            **r,
            "identity_pct": hit["identity_pct"],
        })
    if not pts:
        print("[overview] no overlapping candidates with both CheckV & demarcation identity; "
              "skipping novelty_vs_completeness", file=sys.stderr)
        return

    hosts = sorted({p["host"] for p in pts})
    host_col = {h: ATLASSIAN_HOSTS[i % len(ATLASSIAN_HOSTS)] for i, h in enumerate(hosts)}

    fig, ax = plt.subplots(figsize=(8.6, 5.4), constrained_layout=True)
    for h in hosts:
        sub = [p for p in pts if p["host"] == h]
        ax.scatter(
            [p["identity_pct"] for p in sub],
            [p["completeness"] for p in sub],
            s=[max(24, min(220, 18 * (p["length"] / 1000.0))) for p in sub],
            c=host_col[h], label=f"{h} (n={len(sub)})",
            alpha=0.82, edgecolors="white", linewidths=0.6, zorder=3,
        )

    ax.axvline(50.0, color=INK_SOFT, ls="--", lw=0.9, zorder=1)
    ax.axhline(80.0, color=INK_SOFT, ls=":", lw=0.9, zorder=1)
    ax.set_xlim(0, 102)
    ax.set_ylim(-3, 105)
    ax.set_xlabel("Max AA identity to nearest ICTV reference (%)", fontsize=9.5, color=INK)
    ax.set_ylabel("CheckV completeness (%)", fontsize=9.5, color=INK)
    ax.set_title("Novelty vs. Genome Completeness Across Families",
                 loc="left", fontsize=11, fontweight="bold", color=INK)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right", title="Host / Source", title_fontsize=8.5)
    save_all(fig, out_dir, "novelty_vs_completeness")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phylogeny-dir", default="results/phylogeny")
    ap.add_argument("--demarcation-dir", default="results/demarcation")
    ap.add_argument("--out-dir", default="figures/overview")
    ap.add_argument("--group-by", default="host",
                    help="Column in selected.tsv for x-axis grouping (host or bioproject)")
    ap.add_argument("--max-families", type=int, default=30,
                    help="Max families to show in host_family_matrix")
    a = ap.parse_args()

    rows = load_candidates(Path(a.phylogeny_dir), a.group_by)
    if not rows:
        print(f"[overview] no candidate rows found under {a.phylogeny_dir}/*/selected.tsv",
              file=sys.stderr)
        return 1

    out_dir = Path(a.out_dir)
    plot_host_family_matrix(rows, out_dir, a.max_families)

    nearest = load_nearest_ictv(Path(a.demarcation_dir))
    if nearest:
        plot_novelty_vs_completeness(rows, nearest, out_dir)
    else:
        print(f"[overview] no demarcation pairwise_identity.tsv found yet under {a.demarcation_dir}; "
              "skipping novelty_vs_completeness", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
