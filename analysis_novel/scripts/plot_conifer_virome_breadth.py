#!/usr/bin/env python3
"""Publication-grade Conifer Pest Virome ordination (t-SNE / PCoA), run-level
prevalence, and Shi et al. (2016)-style RdRp novelty vs. CheckV completeness plots.

Reads `results/phylogeny/*/selected.tsv` and (optionally)
`results/demarcation/*/pairwise_identity.tsv` or `master_metadata.csv`.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

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
BP_MARKERS = ["o", "s", "^", "D", "P", "v", "X", "*"]


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


def load_candidates(phylo_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for tsv in sorted(phylo_dir.glob("*/selected.tsv")):
        fam = tsv.parent.name
        with tsv.open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("origin") == "reference":
                    continue
                lbl = (r.get("label") or "").strip()
                run = lbl.split("_")[0] if lbl else "unknown"
                host = norm_host(r.get("host", ""))
                bp = (r.get("bioproject") or "").strip() or "Unknown BP"
                comp = as_float(r.get("checkv_completeness"))
                rows.append({
                    "family": fam,
                    "label": lbl,
                    "run": run,
                    "host": host,
                    "bioproject": bp,
                    "completeness": comp if comp is not None else 0.0,
                })
    return rows


def load_rdrp_identity(demarc_dir: Path | None, master_meta: Path | None) -> dict[str, tuple[float, float]]:
    """Return label/run -> (identity_pct, coverage_pct) from demarcation or master_metadata."""
    best: dict[str, tuple[float, float]] = {}
    if demarc_dir and demarc_dir.exists():
        for tsv in sorted(demarc_dir.glob("*/pairwise_identity.tsv")):
            with tsv.open(newline="", encoding="utf-8", errors="ignore") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    cand = (r.get("candidate") or "").strip()
                    ident = as_float(r.get("identity_pct"))
                    cov = as_float(r.get("query_coverage_pct")) or 0.0
                    if not cand or ident is None:
                        continue
                    keys = [cand]
                    for part in cand.split("|"):
                        if part and part not in {"cand", "ICTV"}:
                            keys.append(part)
                    for k in keys:
                        if k not in best or ident > best[k][0]:
                            best[k] = (ident, cov)
    if master_meta and master_meta.exists():
        with master_meta.open(newline="", encoding="utf-8", errors="ignore") as fh:
            for r in csv.DictReader(fh):
                ident = None
                for col in ("rdrp_pident", "pident", "best_pident", "identity_pct", "diamond_pident"):
                    if col in r:
                        ident = as_float(r.get(col))
                        if ident is not None:
                            break
                if ident is None:
                    continue
                cov = as_float(r.get("query_coverage_pct") or r.get("qcovhsp") or "100") or 100.0
                for col in ("contig_name", "contig", "label", "seq_id"):
                    val = (r.get(col) or "").strip()
                    if val and (val not in best or ident > best[val][0]):
                        best[val] = (ident, cov)
    return best


def compute_pcoa(X: np.ndarray) -> tuple[np.ndarray, tuple[float, float]]:
    """Bray-Curtis dissimilarity + classical MDS (PCoA) with % variance explained."""
    n = X.shape[0]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(X[i] - X[j]).sum()
            den = (X[i] + X[j]).sum()
            d = (num / den) if den > 0 else 0.0
            D[i, j] = D[j, i] = d
    D2 = D ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D2 @ H
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1]
    pos_sum = float(np.sum(vals[vals > 0])) or 1.0
    v1 = max(float(vals[idx[0]]), 0.0)
    v2 = max(float(vals[idx[1]]), 0.0) if n > 1 else 0.0
    coords = np.column_stack([vecs[:, idx[0]] * np.sqrt(v1), vecs[:, idx[1]] * np.sqrt(v2)])
    return coords, (100.0 * v1 / pos_sum, 100.0 * v2 / pos_sum)


def load_community_from_master(master_meta: Path) -> tuple[dict[str, str], dict[str, str], dict[tuple[str, str], float]]:
    """Aggregate log1p(sum of SPAdes contig coverage) per (run, consensus_family) from master_metadata.csv."""
    import re
    run_host: dict[str, str] = {}
    run_bp: dict[str, str] = {}
    abund: dict[tuple[str, str], float] = defaultdict(float)
    with master_meta.open(newline="", encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh):
            run = (r.get("run_id") or "").strip()
            fam = (r.get("consensus_family") or "").strip()
            if not run or not fam or fam.lower() in {"unclassified", "unknown", "na"}:
                continue
            host = norm_host(r.get("scientificname") or r.get("family_context") or "")
            if host not in HOST_ORDER:
                continue
            run_host[run] = host
            run_bp[run] = (r.get("bioproject") or "").strip() or "Unknown BP"
            cid = r.get("contig_id") or ""
            m = re.search(r"cov_([0-9.]+)", cid)
            cov = float(m.group(1)) if m else 1.0
            abund[(run, fam)] += cov
    return run_host, run_bp, abund


def load_rdrp_evalue(master_meta: Path | None) -> list[tuple[float, float, str, str]]:
    """Return list of (-log10(rdrpcatch_evalue), checkv_completeness, host, family)."""
    pts: list[tuple[float, float, str, str]] = []
    if not master_meta or not master_meta.exists():
        return pts
    with master_meta.open(newline="", encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh):
            ev = as_float(r.get("rdrpcatch_evalue"))
            comp = as_float(r.get("checkv_completeness"))
            if ev is None or comp is None:
                continue
            neg_log_ev = -np.log10(max(ev, 1e-300))
            host = norm_host(r.get("scientificname") or r.get("family_context") or "")
            fam = (r.get("consensus_family") or r.get("rdrpcatch_family") or "Unknown").strip()
            if host in HOST_ORDER:
                pts.append((float(neg_log_ev), float(comp), host, fam))
    return pts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phylogeny-dir", default="results/phylogeny", type=Path)
    ap.add_argument("--demarcation-dir", default="results/demarcation", type=Path)
    ap.add_argument("--master-metadata", default=None, type=Path,
                    help="Phase2 master_metadata.csv for abundance-weighted PCoA & RdRp E-values")
    ap.add_argument("--out-dir", default="figures/overview", type=Path)
    ap.add_argument("--top-families", type=int, default=16)
    args = ap.parse_args()

    rows = load_candidates(args.phylogeny_dir)
    if not rows:
        raise SystemExit(f"No candidate rows found under {args.phylogeny_dir}/*/selected.tsv")

    run_host: dict[str, str] = {}
    run_bp: dict[str, str] = {}
    fam_runs: dict[str, set[str]] = defaultdict(set)
    fam_cnt: dict[str, int] = defaultdict(int)
    for r in rows:
        if r["host"] in HOST_ORDER:
            run_host[r["run"]] = r["host"]
            run_bp[r["run"]] = r["bioproject"]
        fam_runs[r["family"]].add(r["run"])
        fam_cnt[r["family"]] += 1

    ordered_runs = sorted(run_host.keys(), key=lambda u: (HOST_ORDER.index(run_host[u]), u))
    run_idx = {u: i for i, u in enumerate(ordered_runs)}

    # Abundance-weighted community matrix from master_metadata.csv (all viral contigs)
    if args.master_metadata and args.master_metadata.exists():
        m_host, m_bp, abund = load_community_from_master(args.master_metadata)
        for u in m_host:
            run_host.setdefault(u, m_host[u])
            run_bp[u] = m_bp[u]
        pcoa_runs = sorted(m_host.keys(), key=lambda u: (HOST_ORDER.index(m_host[u]), u))
        pcoa_fams = sorted({f for (_, f) in abund})
        X = np.zeros((len(pcoa_runs), len(pcoa_fams)), dtype=float)
        p_ri = {u: i for i, u in enumerate(pcoa_runs)}
        p_fi = {f: j for j, f in enumerate(pcoa_fams)}
        for (u, f), cov in abund.items():
            X[p_ri[u], p_fi[f]] = np.log1p(cov)
        pcoa_subtitle = f"Bray-Curtis PCoA (all contigs, {len(pcoa_fams)} families)"
    else:
        pcoa_runs = ordered_runs
        pcoa_fams = sorted(fam_cnt.keys())
        X = np.zeros((len(pcoa_runs), len(pcoa_fams)), dtype=float)
        p_fi = {f: j for j, f in enumerate(pcoa_fams)}
        for r in rows:
            ri = run_idx.get(r["run"])
            fi = p_fi.get(r["family"])
            if ri is not None and fi is not None:
                X[ri, fi] += np.log1p(r["completeness"] + 1.0)
        pcoa_subtitle = f"Bray-Curtis PCoA ({len(pcoa_fams)} families)"

    coords, (var1, var2) = compute_pcoa(X)

    # Top families for prevalence heatmap
    top_fams = sorted(fam_runs.keys(), key=lambda f: (-len(fam_runs[f]), -fam_cnt[f], f))[: args.top_families]
    fam_idx_top = {f: i for i, f in enumerate(top_fams)}
    cell_comp = np.full((len(top_fams), len(ordered_runs)), np.nan)
    cell_n = np.zeros((len(top_fams), len(ordered_runs)), dtype=int)
    for r in rows:
        fi = fam_idx_top.get(r["family"])
        rj = run_idx.get(r["run"])
        if fi is None or rj is None:
            continue
        cell_n[fi, rj] += 1
        cur = cell_comp[fi, rj]
        cell_comp[fi, rj] = r["completeness"] if np.isnan(cur) else max(cur, r["completeness"])

    # Main 2-panel figure: A = Abundance-Weighted Virome PCoA, B = Run-level Prevalence
    fig = plt.figure(figsize=(14.4, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.45])

    ax_a = fig.add_subplot(gs[0, 0])
    bps = sorted(set(run_bp.get(u, "Unknown BP") for u in pcoa_runs))
    bp_marker = {bp: BP_MARKERS[i % len(BP_MARKERS)] for i, bp in enumerate(bps)}
    seen_labels = set()
    for h in HOST_ORDER:
        idxs = [i for i, u in enumerate(pcoa_runs) if run_host[u] == h]
        if not idxs:
            continue
        cx, cy = float(np.mean(coords[idxs, 0])), float(np.mean(coords[idxs, 1]))
        for i in idxs:
            ax_a.plot([cx, coords[i, 0]], [cy, coords[i, 1]], color=HOST_COL[h], lw=0.8, alpha=0.35, zorder=2)
            bp = run_bp.get(pcoa_runs[i], "Unknown BP")
            lbl = f"{HOST_SHORT[h]} ({bp})" if len(bps) > 1 else HOST_SHORT[h]
            ax_a.scatter(
                coords[i, 0], coords[i, 1],
                s=75, c=HOST_COL[h], marker=bp_marker[bp],
                edgecolors="white", linewidths=0.8, alpha=0.92, zorder=3,
                label=lbl if lbl not in seen_labels else None,
            )
            seen_labels.add(lbl)
        ax_a.scatter([cx], [cy], s=130, c=HOST_COL[h], marker="X", edgecolors="white", linewidths=1.0, zorder=4)
    ax_a.set_xlabel(f"PCoA 1 ({var1:.1f}% variance)", fontsize=9, color=INK)
    ax_a.set_ylabel(f"PCoA 2 ({var2:.1f}% variance)", fontsize=9, color=INK)
    ax_a.set_title("A", loc="left", fontsize=11.0, fontweight="bold", color=INK)
    ax_a.grid(color=GRID, lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax_a.spines[sp].set_visible(False)
    ax_a.legend(frameon=False, fontsize=7.8, loc="best")

    ax_b = fig.add_subplot(gs[0, 1])
    im = ax_b.imshow(cell_comp, aspect="auto", cmap=COMP_CMAP, vmin=0, vmax=100, interpolation="nearest")
    ax_b.set_xticks(range(len(ordered_runs)))
    ax_b.set_xticklabels(ordered_runs, rotation=60, ha="right", fontsize=7.0, color=INK)
    ax_b.set_yticks(range(len(top_fams)))
    ax_b.set_yticklabels([f"{f} ({len(fam_runs[f])}/{len(ordered_runs)})" for f in top_fams], fontsize=8.2, fontstyle="italic", color=INK)
    ax_b.set_ylabel("Viral family (run prevalence)", fontsize=9, color=INK)

    for i in range(len(top_fams)):
        for j in range(len(ordered_runs)):
            if cell_n[i, j] > 0:
                val = cell_comp[i, j]
                ax_b.text(j, i, str(cell_n[i, j]), ha="center", va="center",
                          fontsize=6.8, fontweight="bold",
                          color="white" if val >= 55 else INK)

    host_spans: list[tuple[str, int, int]] = []
    cur_h, start_j = None, 0
    for j, u in enumerate(ordered_runs):
        h = run_host[u]
        if h != cur_h:
            if cur_h is not None:
                host_spans.append((cur_h, start_j, j - 1))
                ax_b.axvline(j - 0.5, color=INK, lw=1.2, ls="--")
            cur_h, start_j = h, j
    if cur_h is not None:
        host_spans.append((cur_h, start_j, len(ordered_runs) - 1))

    for h, j0, j1 in host_spans:
        ax_b.plot([j0 - 0.4, j1 + 0.4], [-0.72, -0.72], color=HOST_COL[h], lw=4.0, clip_on=False)
        ax_b.text((j0 + j1) / 2.0, -1.12, HOST_SHORT[h], ha="center", va="bottom",
                  fontsize=8.2, fontweight="bold", fontstyle="italic", color=HOST_COL[h], clip_on=False)

    ax_b.set_title("B", loc="left", fontsize=11.0, fontweight="bold", color=INK, pad=22)
    ax_b.set_xlabel("Sequencing run (grouped by insect host)", fontsize=9, color=INK)
    ax_b.set_xticks(np.arange(-0.5, len(ordered_runs), 1), minor=True)
    ax_b.set_yticks(np.arange(-0.5, len(top_fams), 1), minor=True)
    ax_b.grid(which="minor", color="#FFFFFF", lw=0.8)
    ax_b.tick_params(which="minor", length=0)

    cbar = fig.colorbar(im, ax=ax_b, fraction=0.03, pad=0.02)
    cbar.set_label("Max CheckV completeness (%)", fontsize=8.2, color=INK_SOFT)
    cbar.ax.tick_params(labelsize=7.5, colors=INK_SOFT)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(args.out_dir / f"conifer_virome_breadth.{ext}",
                    bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.out_dir}/conifer_virome_breadth.(pdf|svg|png)")

    # Standalone title-free Bray-Curtis PCoA panel
    fig_p, ax_p = plt.subplots(figsize=(5.6, 4.6), constrained_layout=True)
    seen_p = set()
    for h in HOST_ORDER:
        idxs = [i for i, u in enumerate(pcoa_runs) if run_host[u] == h]
        if not idxs:
            continue
        cx, cy = float(np.mean(coords[idxs, 0])), float(np.mean(coords[idxs, 1]))
        for i in idxs:
            ax_p.plot([cx, coords[i, 0]], [cy, coords[i, 1]], color=HOST_COL[h], lw=0.9, alpha=0.38, zorder=2)
            bp = run_bp.get(pcoa_runs[i], "Unknown BP")
            lbl = f"{HOST_SHORT[h]} ({bp})" if len(bps) > 1 else HOST_SHORT[h]
            ax_p.scatter(
                coords[i, 0], coords[i, 1],
                s=82, c=HOST_COL[h], marker=bp_marker[bp],
                edgecolors="white", linewidths=0.8, alpha=0.94, zorder=3,
                label=lbl if lbl not in seen_p else None,
            )
            seen_p.add(lbl)
        ax_p.scatter([cx], [cy], s=140, c=HOST_COL[h], marker="X", edgecolors="white", linewidths=1.0, zorder=4)
    ax_p.set_xlabel(f"PCoA 1 ({var1:.1f}% variance)", fontsize=9.5, color=INK)
    ax_p.set_ylabel(f"PCoA 2 ({var2:.1f}% variance)", fontsize=9.5, color=INK)
    ax_p.grid(color=GRID, lw=0.6, zorder=0)
    for sp in ("top", "right"):
        ax_p.spines[sp].set_visible(False)
    ax_p.tick_params(labelsize=8.2, colors=INK_SOFT)
    ax_p.legend(frameon=False, fontsize=8.0, loc="best")
    for ext in ("pdf", "svg", "png"):
        fig_p.savefig(args.out_dir / f"conifer_pcoa.{ext}",
                      bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig_p)
    print(f"Wrote {args.out_dir}/conifer_pcoa.(pdf|svg|png)")

    # RdRp-CATCH -log10(E-value) vs CheckV completeness (%) from master_metadata.csv
    ev_pts = load_rdrp_evalue(args.master_metadata)
    if ev_pts:
        fig2, ax2 = plt.subplots(figsize=(7.8, 5.2), constrained_layout=True)
        for h in HOST_ORDER:
            pts = [p for p in ev_pts if p[2] == h]
            if not pts:
                continue
            ax2.scatter([p[0] for p in pts], [p[1] for p in pts], s=28, c=HOST_COL[h],
                        alpha=0.78, edgecolors="white", linewidths=0.4,
                        label=f"{HOST_SHORT[h]} (n={len(pts)})", zorder=3)
        ax2.axvline(5.0, color="#C75300", lw=1.0, ls="--", label="E = 1e-5 threshold")
        ax2.axhline(50.0, color="#1868DB", lw=1.0, ls=":", label="50% CheckV completeness")
        ax2.set_xlabel("RdRp-CATCH homology support (−log10 E-value)", fontsize=9, color=INK)
        ax2.set_ylabel("CheckV genome completeness (%)", fontsize=9, color=INK)
        ax2.set_title("Conifer Pest Virome: RdRp homology support vs. genome completeness",
                      loc="left", fontsize=10.2, fontweight="bold", color=INK)
        ax2.grid(color=GRID, lw=0.6, zorder=0)
        for sp in ("top", "right"):
            ax2.spines[sp].set_visible(False)
        ax2.legend(frameon=False, fontsize=7.8, loc="best")
        for ext in ("pdf", "svg", "png"):
            fig2.savefig(args.out_dir / f"conifer_rdrp_novelty.{ext}",
                         bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
        plt.close(fig2)
        print(f"Wrote {args.out_dir}/conifer_rdrp_novelty.(pdf|svg|png)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
