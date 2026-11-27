#!/usr/bin/env python3
"""Render Manuscript Fig. 5: Cross-family evolutionary constraint & dinucleotide
host-composition signatures across all analysed viral families (title-free).

Panel A: Median dN/dS (omega) and MEME episodic positive site counts ordered by
         genome class (dsRNA, ssRNA, RT, DNA/phage), showing pervasive purifying
         selection (0 pervasive positive sites across all families).
Panel B: Paired Candidate vs. ICTV Reference CpG O/E dumbbell plot across
         multi-candidate families, highlighting arthropod RNA-virus concordance
         and the Partitiviridae mycobiome CpG gap.
"""
from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

INK = "#101214"
MUTED = "#505258"
GRID = "#E2E8F0"

GENOME_CLASS = {
    "Partitiviridae": "dsRNA", "Spinareoviridae": "dsRNA", "Sedoreoviridae": "dsRNA",
    "Totiviridae": "dsRNA", "Chrysoviridae": "dsRNA", "Duplornaviricota": "dsRNA",
    "Duplopiviricetes": "dsRNA",
    "Iflaviridae": "ssRNA", "Dicistroviridae": "ssRNA", "Flaviviridae": "ssRNA",
    "Amarillovirales": "ssRNA", "Polycipiviridae": "ssRNA", "Solemoviridae": "ssRNA",
    "Phenuiviridae": "ssRNA", "Aliusviridae": "ssRNA", "Phasmaviridae": "ssRNA",
    "Rhabdoviridae": "ssRNA", "Nyamiviridae": "ssRNA", "Orthomyxoviridae": "ssRNA",
    "Chuviridae": "ssRNA", "Xinmoviridae": "ssRNA", "Peribunyaviridae": "ssRNA",
    "Steitzviridae": "ssRNA", "Narnaviridae": "ssRNA", "Tombusviridae": "ssRNA",
    "Roniviridae": "ssRNA", "Tolivirales": "ssRNA",
    "Retroviridae": "RT", "Metaviridae": "RT", "Pseudoviridae": "RT",
    "Zimmerviridae": "DNA/Phage", "Peduoviridae": "DNA/Phage", "Adintoviridae": "DNA/Phage",
    "Baculoviridae": "DNA/Phage", "Phycodnaviridae": "DNA/Phage", "Mimiviridae": "DNA/Phage",
    "Chimalliviridae": "DNA/Phage", "Straboviridae": "DNA/Phage", "Demerecviridae": "DNA/Phage",
    "Vilmaviridae": "DNA/Phage", "Sarkviridae": "DNA/Phage", "Zierdtviridae": "DNA/Phage",
    "Iridoviridae": "DNA/Phage",
}
CLASS_ORDER = ["dsRNA", "ssRNA", "RT", "DNA/Phage", "Other"]
CLASS_COL = {
    "dsRNA": "#16794F",
    "ssRNA": "#1868DB",
    "RT": "#8F4FDB",
    "DNA/Phage": "#C75300",
    "Other": "#94A3B8",
}


def as_float(v: str | None) -> float | None:
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def load_selection_summaries(sig_dir: Path) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for fam_dir in sorted(sig_dir.iterdir()):
        if not fam_dir.is_dir():
            continue
        fam = fam_dir.name
        sel_sum = fam_dir / "selection_summary.tsv"
        omega: float | None = None
        n_episodic = 0
        n_pervasive = 0
        if sel_sum.exists():
            kv: dict[str, str] = {}
            with sel_sum.open(encoding="utf-8", errors="ignore") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                fnames = reader.fieldnames or []
                if "metric" in fnames and "value" in fnames:
                    for r in reader:
                        k = (r.get("metric") or "").strip()
                        v = (r.get("value") or "").strip()
                        if k:
                            kv[k] = v
                else:
                    for r in reader:
                        kv.update({k: v for k, v in r.items() if k and v is not None})
                        break
            omega = as_float(kv.get("median_omega") or kv.get("median_dnds") or kv.get("omega_median"))
            n_episodic = int(float(kv.get("n_episodic_positive") or kv.get("meme_positive") or 0))
            n_pervasive = int(float(kv.get("n_pervasive_positive") or kv.get("fubar_positive") or 0))
        if omega is not None:
            cls = GENOME_CLASS.get(fam, "Other")
            out.append({
                "family": fam, "class": cls, "omega": omega,
                "n_episodic": n_episodic, "n_pervasive": n_pervasive,
            })
    out.sort(key=lambda x: (CLASS_ORDER.index(str(x["class"])), float(x["omega"])))
    return out


def load_cpg_pairs(sig_dir: Path, min_candidates: int = 2) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for fam_dir in sorted(sig_dir.iterdir()):
        if not fam_dir.is_dir():
            continue
        fam = fam_dir.name
        comp_tsv = next(
            (p for p in (
                fam_dir / "dinucleotide_odds.tsv",
                fam_dir / "all_orf_dinuc.tsv",
                fam_dir / "all_orf_composition.tsv",
                fam_dir / "composition_signatures.tsv",
            ) if p.exists()),
            None,
        )
        if not comp_tsv:
            continue
        cand_vals: list[float] = []
        ref_vals: list[float] = []
        with comp_tsv.open(encoding="utf-8", errors="ignore") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cpg = as_float(r.get("CG") or r.get("CpG") or r.get("cpg_oe") or r.get("CpG_OE"))
                if cpg is None:
                    continue
                tid = (r.get("genome") or r.get("taxon_id") or r.get("label") or r.get("seq_id") or "").strip()
                orig = (r.get("origin") or r.get("group") or "").strip().lower()
                is_ref = (orig in {"reference", "ictv", "ref"}) or tid.startswith("ICTV_") or (
                    not any(tid.startswith(p) for p in ("SRR", "ERR", "DRR"))
                )
                if is_ref:
                    ref_vals.append(cpg)
                else:
                    cand_vals.append(cpg)
        if len(cand_vals) >= min_candidates and len(ref_vals) >= 1:
            out.append({
                "family": fam,
                "class": GENOME_CLASS.get(fam, "Other"),
                "cand_cpg": statistics.median(cand_vals),
                "ref_cpg": statistics.median(ref_vals),
                "n_cand": len(cand_vals),
                "n_ref": len(ref_vals),
                "gap": abs(statistics.median(cand_vals) - statistics.median(ref_vals)),
            })
    out.sort(key=lambda x: float(x["cand_cpg"]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signatures-dir", default=Path("results/signatures"), type=Path)
    ap.add_argument("--out-dir", default=Path("figures/manuscript"), type=Path)
    ap.add_argument("--min-cpg-candidates", type=int, default=2)
    args = ap.parse_args()

    sel = load_selection_summaries(args.signatures_dir)
    cpg = load_cpg_pairs(args.signatures_dir, min_candidates=args.min_cpg_candidates)
    if not sel and not cpg:
        raise SystemExit(f"No selection or composition summaries found under {args.signatures_dir}")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.5, max(5.6, 0.25 * max(len(sel), len(cpg), 12) + 1.8)))
    plt.subplots_adjust(wspace=0.42, bottom=0.16, top=0.95, left=0.15, right=0.97)

    # Panel A: Median dN/dS (omega) by genome class
    if sel:
        y = list(range(len(sel)))
        cols = [CLASS_COL[str(r["class"])] for r in sel]
        vals = [float(r["omega"]) for r in sel]
        ax_a.barh(y, vals, color=cols, height=0.72, alpha=0.88)
        for i, r in enumerate(sel):
            ep = int(r["n_episodic"])
            ax_a.text(float(r["omega"]) + 0.012, i, f"ω={float(r['omega']):.3f} ({ep} ep.)",
                      va="center", fontsize=6.6, color=INK)
        ax_a.set_yticks(y)
        ax_a.set_yticklabels([str(r["family"]) for r in sel], fontsize=7.6, fontstyle="italic", color=INK)
        ax_a.invert_yaxis()
        ax_a.set_xlim(0, max(vals) * 1.38)
    ax_a.axvline(1.0, color="#DC2626", ls="--", lw=0.9, alpha=0.7)
    ax_a.set_xlabel("Median dN/dS (ω; FUBAR) · parentheses = MEME episodic sites", fontsize=8.4, color=INK)
    ax_a.grid(axis="x", color=GRID, lw=0.6)
    for sp in ("top", "right"):
        ax_a.spines[sp].set_visible(False)
    ax_a.text(-0.18, 1.01, "A", transform=ax_a.transAxes, fontsize=12, fontweight="bold", color=INK)

    # Panel B: Candidate vs. ICTV Reference CpG O/E dumbbell plot
    if cpg:
        y2 = list(range(len(cpg)))
        for i, r in enumerate(cpg):
            c_val, r_val = float(r["cand_cpg"]), float(r["ref_cpg"])
            ax_b.plot([c_val, r_val], [i, i], color="#94A3B8", lw=1.8, zorder=2)
            ax_b.scatter([c_val], [i], color="#1868DB", s=46, zorder=3, edgecolors="white", lw=0.6)
            ax_b.scatter([r_val], [i], color="#64748B", s=42, marker="D", zorder=3, edgecolors="white", lw=0.6)
        ax_b.set_yticks(y2)
        ax_b.set_yticklabels([f"{r['family']} (n={r['n_cand']})" for r in cpg], fontsize=7.8, fontstyle="italic", color=INK)
        ax_b.invert_yaxis()
    ax_b.axvline(1.0, color=MUTED, ls=":", lw=1.0)
    ax_b.set_xlabel("Median CpG Observed / Expected ratio (O/E)", fontsize=8.4, color=INK)
    ax_b.grid(axis="x", color=GRID, lw=0.6)
    for sp in ("top", "right"):
        ax_b.spines[sp].set_visible(False)
    ax_b.text(-0.18, 1.01, "B", transform=ax_b.transAxes, fontsize=12, fontweight="bold", color=INK)

    # Legends cleanly placed below axes (no overlap, no figure title)
    cls_handles = [Patch(facecolor=CLASS_COL[k], label=k) for k in CLASS_ORDER if any(r["class"] == k for r in sel)]
    if cls_handles:
        ax_a.legend(handles=cls_handles, title="Genome class", loc="upper center",
                    bbox_to_anchor=(0.48, -0.10), ncol=len(cls_handles), frameon=False, fontsize=7.6, title_fontsize=8.0)
    cpg_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#1868DB", markersize=7, label="Candidate median"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="#64748B", markersize=6.5, label="ICTV Reference median"),
    ]
    ax_b.legend(handles=cpg_handles, loc="upper center", bbox_to_anchor=(0.48, -0.10),
                ncol=2, frameon=False, fontsize=7.8)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(args.out_dir / f"fig5_evolutionary_constraint.{ext}",
                    bbox_inches="tight", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.out_dir}/fig5_evolutionary_constraint.(pdf|svg|png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
