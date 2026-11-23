#!/usr/bin/env python3
"""ComplexHeatmap-style multi-track virome overview for the Conifer Pest Virome.

Reads phase2 `master_metadata.csv` (and optionally `results/phylogeny/*/selected.tsv`)
and renders a cohesive multi-track figure (PDF, SVG, 300-dpi PNG):
  - Top column tracks: Viral richness per run, BioProject strip, Host species strip
  - Left row track:    Viral genome class (+ssRNA, -ssRNA, dsRNA, RT, DNA)
  - Main matrix:       Top viral families x 21 runs (fill = max CheckV completeness %,
                       cell text = contig count per run)
  - Right row tracks:  (1) Best RdRp-CATCH support (-log10 E-value),
                       (2) Contig length range (kb, min-median-max),
                       (3) Run prevalence (n / 21 runs)
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

INK, INK_SOFT, GRID = "#101214", "#505258", "#E2E8F0"

COMP_CMAP = LinearSegmentedColormap.from_list(
    "atlassian_comp", ["#F1F5F9", "#60A5FA", "#1868DB", "#1C2B42"]
)

HOST_ORDER = ["Locusta migratoria", "Hylobius abietis", "Ips confusus"]
HOST_SHORT = {
    "Locusta migratoria": "L. migratoria",
    "Hylobius abietis": "H. abietis",
    "Ips confusus": "I. confusus",
}
HOST_COL = {
    "Locusta migratoria": "#C75300",
    "Hylobius abietis": "#16794F",
    "Ips confusus": "#1868DB",
}

BP_PALETTE = ["#2898BD", "#8F4FDB", "#D97706", "#0D9488", "#E11D48", "#4F46E5", "#65A30D"]

BALTIMORE_MAP = {
    "Iflaviridae": "ssRNA", "Dicistroviridae": "ssRNA", "Flaviviridae": "ssRNA",
    "Amarillovirales": "ssRNA", "Polycipiviridae": "ssRNA", "Tombusviridae": "ssRNA",
    "Solemoviridae": "ssRNA", "Steitzviridae": "ssRNA", "Narnaviridae": "ssRNA",
    "Picornaviridae": "ssRNA", "Hepeviridae": "ssRNA", "Virgaviridae": "ssRNA",
    "Tymoviridae": "ssRNA", "Nodaviridae": "ssRNA", "Botourmiaviridae": "ssRNA",
    "Mitoviridae": "ssRNA", "Mesoniviridae": "ssRNA",
    "Phenuiviridae": "ssRNA", "Xinmoviridae": "ssRNA", "Nyamiviridae": "ssRNA",
    "Orthomyxoviridae": "ssRNA", "Aliusviridae": "ssRNA", "Phasmaviridae": "ssRNA",
    "Rhabdoviridae": "ssRNA", "Chuviridae": "ssRNA", "Qinviridae": "ssRNA",
    "Peribunyaviridae": "ssRNA", "Arenaviridae": "ssRNA",
    "Partitiviridae": "dsRNA", "Spinareoviridae": "dsRNA", "Duplornaviricota": "dsRNA",
    "Totiviridae": "dsRNA", "Sedoreoviridae": "dsRNA", "Chrysoviridae": "dsRNA",
    "Retroviridae": "RT", "Metaviridae": "RT", "Pseudoviridae": "RT",
    "Peduoviridae": "DNA", "Zimmerviridae": "DNA", "Adintoviridae": "DNA",
    "Straboviridae": "DNA", "Chimalliviridae": "DNA", "Vilmaviridae": "DNA",
    "Demerecviridae": "DNA", "Iridoviridae": "DNA", "Phycodnaviridae": "DNA",
    "Baculoviridae": "DNA", "Duplopiviricetes": "DNA", "Mimiviridae": "DNA",
}

BALT_COL = {
    "ssRNA": "#1868DB",
    "dsRNA": "#16794F",
    "RT": "#8F4FDB",
    "DNA": "#2898BD",
    "Other": "#94A3B8",
}


def norm_host(h: str) -> str:
    s = h.strip().lower().replace("_", " ")
    if "locusta" in s:
        return "Locusta migratoria"
    if "hylobius" in s:
        return "Hylobius abietis"
    if "ips" in s:
        return "Ips confusus"
    return h.strip() or "Other"


def as_float(v: str | None) -> float | None:
    try:
        x = float(str(v).strip())
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def parse_len_kb(contig_id: str) -> float | None:
    m = re.search(r"length_(\d+)", contig_id)
    return (int(m.group(1)) / 1000.0) if m else None


def load_master(master_csv: Path) -> list[dict]:
    rows: list[dict] = []
    with master_csv.open(newline="", encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh):
            fam = (r.get("consensus_family") or "").strip()
            run = (r.get("run_id") or "").strip()
            if not fam or not run or fam.lower() in {"unclassified", "unknown", "na"}:
                continue
            host = norm_host(r.get("scientificname") or r.get("family_context") or "")
            if host not in HOST_ORDER:
                continue
            bp = (r.get("bioproject") or "").strip() or "Unknown"
            comp = as_float(r.get("checkv_completeness"))
            ev = as_float(r.get("rdrpcatch_evalue"))
            neg_log_ev = float(-np.log10(max(ev, 1e-300))) if ev is not None else None
            length_kb = parse_len_kb(r.get("contig_id") or "")
            rows.append({
                "family": fam,
                "run": run,
                "host": host,
                "bioproject": bp,
                "completeness": comp if comp is not None else 0.0,
                "neg_log_ev": neg_log_ev,
                "length_kb": length_kb,
            })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master-metadata", required=True, type=Path, help="phase2 master_metadata.csv")
    ap.add_argument("--out-dir", default="figures/overview", type=Path)
    ap.add_argument("--top-families", type=int, default=22)
    args = ap.parse_args()

    rows = load_master(args.master_metadata)
    if not rows:
        raise SystemExit(f"No valid viral rows parsed from {args.master_metadata}")

    run_host: dict[str, str] = {}
    run_bp: dict[str, str] = {}
    run_fams: dict[str, set[str]] = defaultdict(set)
    fam_hosts: dict[str, set[str]] = defaultdict(set)
    fam_runs: dict[str, set[str]] = defaultdict(set)
    fam_cnt: dict[str, int] = defaultdict(int)
    fam_evs: dict[str, list[float]] = defaultdict(list)
    fam_lens: dict[str, list[float]] = defaultdict(list)

    for r in rows:
        run_host[r["run"]] = r["host"]
        run_bp[r["run"]] = r["bioproject"]
        run_fams[r["run"]].add(r["family"])
        fam_hosts[r["family"]].add(r["host"])
        fam_runs[r["family"]].add(r["run"])
        fam_cnt[r["family"]] += 1
        if r["neg_log_ev"] is not None:
            fam_evs[r["family"]].append(r["neg_log_ev"])
        if r["length_kb"] is not None:
            fam_lens[r["family"]].append(r["length_kb"])

    ordered_runs = sorted(
        run_host.keys(),
        key=lambda u: (HOST_ORDER.index(run_host[u]), run_bp[u], u),
    )
    top_fams = sorted(
        fam_runs.keys(),
        key=lambda f: (-len(fam_hosts[f]), -len(fam_runs[f]), -fam_cnt[f], f),
    )[: args.top_families]

    n_f, n_r = len(top_fams), len(ordered_runs)
    f_idx = {f: i for i, f in enumerate(top_fams)}
    r_idx = {u: j for j, u in enumerate(ordered_runs)}

    cell_comp = np.full((n_f, n_r), np.nan)
    cell_n = np.zeros((n_f, n_r), dtype=int)
    for r in rows:
        fi = f_idx.get(r["family"])
        rj = r_idx.get(r["run"])
        if fi is None or rj is None:
            continue
        cell_n[fi, rj] += 1
        cur = cell_comp[fi, rj]
        cell_comp[fi, rj] = r["completeness"] if np.isnan(cur) else max(cur, r["completeness"])

    bps = sorted(set(run_bp.values()))
    bp_col = {bp: BP_PALETTE[i % len(BP_PALETTE)] for i, bp in enumerate(bps)}

    # Layout: 4 rows (richness, bioproject strip, host strip, main heatmap)
    #         x 4 cols (genome-class strip, main 21-run matrix, Length kb, Prevalence)
    fig = plt.figure(figsize=(14.6, max(7.4, 0.32 * n_f + 3.0)))
    gs = fig.add_gridspec(
        4, 4,
        height_ratios=[0.75, 0.22, 0.24, n_f * 0.30],
        width_ratios=[0.22, 5.4, 1.35, 1.15],
        hspace=0.06, wspace=0.09,
    )

    # Top track 1: Viral family richness per run
    ax_rich = fig.add_subplot(gs[0, 1])
    rich_vals = [len(run_fams[u]) for u in ordered_runs]
    ax_rich.bar(range(n_r), rich_vals, width=0.72, color=[HOST_COL[run_host[u]] for u in ordered_runs], alpha=0.85)
    ax_rich.set_xlim(-0.5, n_r - 0.5)
    ax_rich.set_xticks([])
    ax_rich.set_ylabel("Families", fontsize=7.5, color=INK)
    ax_rich.tick_params(axis="y", labelsize=7.0, colors=INK_SOFT)
    ax_rich.grid(axis="y", color=GRID, lw=0.5)
    for sp in ("top", "right"):
        ax_rich.spines[sp].set_visible(False)

    # Top track 2: BioProject strip
    ax_bp = fig.add_subplot(gs[1, 1])
    for j, u in enumerate(ordered_runs):
        ax_bp.add_patch(plt.Rectangle((j - 0.48, 0), 0.96, 1, facecolor=bp_col[run_bp[u]], edgecolor="white", lw=0.6))
    ax_bp.set_xlim(-0.5, n_r - 0.5)
    ax_bp.set_ylim(0, 1)
    ax_bp.set_xticks([])
    ax_bp.set_yticks([0.5])
    ax_bp.set_yticklabels(["BioProject"], fontsize=7.5, fontweight="bold", color=INK)
    ax_bp.tick_params(length=0)
    for sp in ax_bp.spines.values():
        sp.set_visible(False)

    # Top track 3: Host species strip
    ax_host = fig.add_subplot(gs[2, 1])
    host_spans: list[tuple[str, int, int]] = []
    cur_h, start_j = None, 0
    for j, u in enumerate(ordered_runs):
        h = run_host[u]
        ax_host.add_patch(plt.Rectangle((j - 0.48, 0), 0.96, 1, facecolor=HOST_COL[h], edgecolor="white", lw=0.6))
        if h != cur_h:
            if cur_h is not None:
                host_spans.append((cur_h, start_j, j - 1))
            cur_h, start_j = h, j
    if cur_h is not None:
        host_spans.append((cur_h, start_j, n_r - 1))
    for h, j0, j1 in host_spans:
        ax_host.text((j0 + j1) / 2.0, 0.5, HOST_SHORT[h], ha="center", va="center",
                     fontsize=7.6, fontweight="bold", fontstyle="italic", color="white")
    ax_host.set_xlim(-0.5, n_r - 0.5)
    ax_host.set_ylim(0, 1)
    ax_host.set_xticks([])
    ax_host.set_yticks([0.5])
    ax_host.set_yticklabels(["Host"], fontsize=7.5, fontweight="bold", color=INK)
    ax_host.tick_params(length=0)
    for sp in ax_host.spines.values():
        sp.set_visible(False)

    # Left row track: Genome class strip
    ax_cls = fig.add_subplot(gs[3, 0])
    for i, f in enumerate(top_fams):
        b = BALTIMORE_MAP.get(f, "Other")
        ax_cls.add_patch(plt.Rectangle((0, i - 0.46), 1, 0.92, facecolor=BALT_COL[b], edgecolor="white", lw=0.6))
    ax_cls.set_xlim(0, 1)
    ax_cls.set_ylim(n_f - 0.5, -0.5)
    ax_cls.set_xticks([0.5])
    ax_cls.set_xticklabels(["Class"], fontsize=7.2, fontweight="bold", color=INK, rotation=45, ha="right")
    ax_cls.set_yticks(range(n_f))
    ax_cls.set_yticklabels(top_fams, fontsize=8.3, fontstyle="italic", color=INK)
    ax_cls.tick_params(length=0)
    for sp in ax_cls.spines.values():
        sp.set_visible(False)

    # Main heatmap matrix (21 runs x top families)
    ax_hm = fig.add_subplot(gs[3, 1])
    im = ax_hm.imshow(cell_comp, aspect="auto", cmap=COMP_CMAP, vmin=0, vmax=100, interpolation="nearest")
    ax_hm.set_xticks(range(n_r))
    ax_hm.set_xticklabels(ordered_runs, rotation=60, ha="right", fontsize=7.0, color=INK)
    ax_hm.set_yticks([])
    ax_hm.set_xticks(np.arange(-0.5, n_r, 1), minor=True)
    ax_hm.set_yticks(np.arange(-0.5, n_f, 1), minor=True)
    ax_hm.grid(which="minor", color="#FFFFFF", lw=0.8)
    ax_hm.tick_params(which="minor", length=0)
    for _, j0, _ in host_spans[1:]:
        ax_hm.axvline(j0 - 0.5, color=INK, lw=1.3, ls="--")
    for i in range(n_f):
        for j in range(n_r):
            if cell_n[i, j] > 0:
                val = cell_comp[i, j]
                ax_hm.text(j, i, str(cell_n[i, j]), ha="center", va="center",
                           fontsize=6.6, fontweight="bold", color="white" if val >= 55 else INK)

    # Right track 1: Contig length range (min-median-max kb)
    ax_len = fig.add_subplot(gs[3, 2], sharey=ax_hm)
    for i, f in enumerate(top_fams):
        ls = fam_lens[f]
        if not ls:
            continue
        lo, med, hi = min(ls), statistics.median(ls), max(ls)
        ax_len.plot([lo, hi], [i, i], color="#2898BD", lw=2.0, zorder=2)
        ax_len.scatter([med], [i], color="#1C2B42", s=22, zorder=3)
    ax_len.set_xlabel("Contig length\n(kb)", fontsize=7.8, color=INK)
    ax_len.tick_params(axis="x", labelsize=7.0, colors=INK_SOFT)
    ax_len.tick_params(axis="y", left=False, labelleft=False)
    ax_len.grid(axis="x", color=GRID, lw=0.5)
    for sp in ("top", "right"):
        ax_len.spines[sp].set_visible(False)

    # Right track 2: Run prevalence (n / 21 runs)
    ax_prev = fig.add_subplot(gs[3, 3], sharey=ax_hm)
    prev_vals = [len(fam_runs[f]) for f in top_fams]
    ax_prev.barh(range(n_f), prev_vals, height=0.68, color="#475569", alpha=0.85)
    for i, v in enumerate(prev_vals):
        ax_prev.text(v + 0.35, i, f"{v}/{n_r}", va="center", fontsize=6.8, fontweight="bold", color=INK)
    ax_prev.set_xlim(0, n_r + 5)
    ax_prev.set_xlabel("Prevalence\n(runs)", fontsize=7.8, color=INK)
    ax_prev.tick_params(axis="x", labelsize=7.0, colors=INK_SOFT)
    ax_prev.tick_params(axis="y", left=False, labelleft=False)
    ax_prev.grid(axis="x", color=GRID, lw=0.5)
    for sp in ("top", "right"):
        ax_prev.spines[sp].set_visible(False)

    # Reserve generous bottom margin so legends never overlap x-axis labels or each other
    fig.subplots_adjust(bottom=0.24, top=0.96, left=0.14, right=0.97)
    ax_hm.set_xlabel("Sequencing run (grouped by host & BioProject; cell number = contig count)",
                     fontsize=8.2, color=INK, labelpad=6)

    cax = fig.add_axes([0.14, 0.085, 0.22, 0.016])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("Max CheckV completeness (%) per run", fontsize=7.6, color=INK)
    cbar.ax.tick_params(labelsize=7.0, colors=INK_SOFT)

    host_handles = [Patch(facecolor=HOST_COL[h], label=HOST_SHORT[h]) for h in HOST_ORDER]
    class_handles = [Patch(facecolor=BALT_COL[k], label=k) for k in ("ssRNA", "dsRNA", "RT", "DNA", "Other")]
    bp_handles = [Patch(facecolor=bp_col[bp], label=bp) for bp in bps]
    fig.legend(handles=host_handles, title="Host species", loc="lower left",
               bbox_to_anchor=(0.42, 0.055), ncol=3, frameon=False, fontsize=7.4, title_fontsize=7.8)
    fig.legend(handles=bp_handles, title="BioProject", loc="lower left",
               bbox_to_anchor=(0.72, 0.055), ncol=2, frameon=False, fontsize=7.4, title_fontsize=7.8)
    fig.legend(handles=class_handles, title="Genome class", loc="lower center",
               bbox_to_anchor=(0.55, 0.002), ncol=5, frameon=False, fontsize=7.4, title_fontsize=7.8)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(args.out_dir / f"conifer_complex_heatmap.{ext}",
                    bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.out_dir}/conifer_complex_heatmap.(pdf|svg|png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
