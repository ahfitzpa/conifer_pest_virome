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


HOST_PALETTE = {
    "ICTV Vertebrate Ref": "#475569",
    "ICTV Reference": "#475569",
    "Vertebrate / Zoonotic": "#1868DB",
    "Ips confusus": "#1868DB",
    "Hylobius abietis": "#16794F",
    "Locusta migratoria": "#C75300",
    "Crustacean": "#C75300",
    "Invertebrate / Insect": "#8F4FDB",
    "Marine / Coastal": "#2898BD",
    "Freshwater": "#16794F",
    "Soil / Wastewater": "#8F4FDB",
    "Other / Unknown": "#94A3B8",
}

BIOPROJECT_GROUP = {
    "PRJEB6505": "Vertebrate / Zoonotic",
    "PRJEB26687": "Vertebrate / Zoonotic",
    "PRJEB17935": "Vertebrate / Zoonotic",
    "PRJNA751908": "Vertebrate / Zoonotic",
    "PRJDB7591": "Vertebrate / Zoonotic",
    "PRJNA683917": "Vertebrate / Zoonotic",
    "PRJNA655509": "Crustacean",
    "PRJNA702424": "Crustacean",
    "PRJNA669582": "Crustacean",
    "PRJNA754774": "Crustacean",
    "PRJEB34490": "Marine / Coastal",
    "PRJNA576751": "Marine / Coastal",
    "PRJNA555787": "Marine / Coastal",
    "PRJNA665326": "Marine / Coastal",
    "PRJNA654739": "Marine / Coastal",
    "PRJNA640812": "Freshwater",
    "PRJNA691101": "Freshwater",
    "PRJNA657047": "Freshwater",
    "PRJNA663589": "Freshwater",
    "PRJNA716119": "Soil / Wastewater",
    "PRJNA616041": "Soil / Wastewater",
    "PRJNA571457": "Soil / Wastewater",
    "PRJNA673783": "Soil / Wastewater",
}

RUN_TO_BIOPROJECT = {
    "DRR155082": "PRJDB7591",
    "ERR1301444": "PRJEB6505", "ERR1301533": "PRJEB6505", "ERR1303023": "PRJEB6505",
    "ERR2001273": "PRJEB6505", "ERR2001276": "PRJEB6505", "ERR3179191": "PRJEB6505",
    "ERR3179286": "PRJEB6505", "ERR3180651": "PRJEB6505", "ERR3180737": "PRJEB6505",
    "ERR1742975": "PRJEB17935", "ERR1742977": "PRJEB17935",
    "ERR2737741": "PRJEB26687", "ERR2737742": "PRJEB26687",
    "ERR3660995": "PRJEB34490", "ERR3660996": "PRJEB34490", "ERR3661002": "PRJEB34490",
    "ERR3661003": "PRJEB34490", "ERR3661004": "PRJEB34490", "ERR3661007": "PRJEB34490",
    "ERR3661008": "PRJEB34490", "ERR3661010": "PRJEB34490", "ERR3661012": "PRJEB34490",
    "ERR3661014": "PRJEB34490", "ERR3661016": "PRJEB34490", "ERR3661018": "PRJEB34490",
    "ERR3661021": "PRJEB34490",
    "SRR10402494": "PRJNA571457",
    "SRR11450579": "PRJNA616041", "SRR11450586": "PRJNA616041",
    "SRR12063560": "PRJNA640812", "SRR12063601": "PRJNA640812", "SRR12063613": "PRJNA640812",
    "SRR12660293": "PRJNA654739",
    "SRR12671109": "PRJNA655509", "SRR12671112": "PRJNA655509", "SRR12671136": "PRJNA655509",
    "SRR12685647": "PRJNA663589",
    "SRR12720641": "PRJNA665326", "SRR12720672": "PRJNA665326", "SRR12720673": "PRJNA665326",
    "SRR12835363": "PRJNA669582",
    "SRR12978839": "PRJNA673783",
    "SRR13259316": "PRJNA683917",
    "SRR13450175": "PRJNA691101", "SRR13450176": "PRJNA691101",
    "SRR13450956": "PRJNA657047",
    "SRR13729885": "PRJNA702424", "SRR13729895": "PRJNA702424",
    "SRR14039655": "PRJNA716119", "SRR14039765": "PRJNA716119", "SRR14039946": "PRJNA716119",
    "SRR14857291": "PRJNA576751", "SRR14857293": "PRJNA576751", "SRR14857296": "PRJNA576751",
    "SRR14857297": "PRJNA576751", "SRR14857300": "PRJNA576751",
    "SRR15340407": "PRJNA751908", "SRR15340415": "PRJNA751908",
    "SRR15462253": "PRJNA754774",
    "SRR16347784": "PRJNA555787", "SRR16347789": "PRJNA555787",
    "SRR16347793": "PRJNA555787", "SRR16347798": "PRJNA555787",
}

HEPE_ORF1_DOMAINS = [
    ("Met", 0.00, 0.14, "#E2E8F0"),
    ("Y", 0.14, 0.25, "#F1F5F9"),
    ("PCP", 0.25, 0.35, "#E2E8F0"),
    ("HVR", 0.35, 0.46, "#FEF3C7"),
    ("X/Macro", 0.46, 0.56, "#F1F5F9"),
    ("Hel", 0.56, 0.72, "#E2E8F0"),
    ("RdRp", 0.72, 1.00, "#DBEAFE"),
]

HEPE_ORF2_DOMAINS = [
    ("N-term / ARM", 0.00, 0.18, "#E2E8F0"),
    ("S (Shell)", 0.18, 0.48, "#DBEAFE"),
    ("M (Middle)", 0.48, 0.68, "#F1F5F9"),
    ("P (Protruding)", 0.68, 1.00, "#FEF3C7"),
]


def _classify_broad_host(host: str) -> str:
    h = host.strip().lower().replace("_", " ")
    if not h or h in {"na", "unknown", "not applicable", "none"}:
        return ""
    if "ips" in h and "confusus" in h:
        return "Ips confusus"
    if "hylobius" in h and "abietis" in h:
        return "Hylobius abietis"
    if "locusta" in h and "migratoria" in h:
        return "Locusta migratoria"
    if any(k in h for k in ("ips", "hylobius", "locusta", "beetle", "insect", "weevil")):
        return "Invertebrate / Insect"
    if any(k in h for k in ("homarus", "lobster", "crayfish", "pacifastacus", "crustace")):
        return "Crustacean"
    if any(k in h for k in ("marine", "coastal", "ocean", "sea")):
        return "Marine / Coastal"
    if any(k in h for k in ("freshwater", "lake", "river")):
        return "Freshwater"
    if any(k in h for k in ("soil", "sludge", "wastewater", "paddy")):
        return "Soil / Wastewater"
    return "Vertebrate / Zoonotic"


def build_host_map(meta_rows: list[dict], family: str = "") -> dict[str, str]:
    """Map genome/taxon_id/accession -> broad Host/Source group."""
    out: dict[str, str] = {}
    for run, bp in RUN_TO_BIOPROJECT.items():
        if bp in BIOPROJECT_GROUP:
            out[run] = BIOPROJECT_GROUP[bp]

    ref_label = "ICTV Vertebrate Ref" if family == "Hepeviridae" else "ICTV Reference"
    id_cols = ("taxon_id", "label", "genome", "taxon", "id", "seq_id", "contig", "name", "run_accession", "accession", "sample")
    for r in meta_rows:
        origin = (r.get("origin") or r.get("group") or "").strip().lower()
        bp = (r.get("bioproject") or r.get("BioProject") or "").strip()
        acc = (r.get("accession") or "").strip()
        host = (r.get("host") or r.get("Host") or r.get("species") or "").strip()
        if not bp and acc in RUN_TO_BIOPROJECT:
            bp = RUN_TO_BIOPROJECT[acc]

        if origin in {"reference", "ictv"}:
            grp = ref_label
        elif bp in BIOPROJECT_GROUP:
            grp = BIOPROJECT_GROUP[bp]
        else:
            grp = _classify_broad_host(host) or _classify_broad_host(r.get("sample", "")) or "Other / Unknown"

        for col in id_cols:
            val = (r.get(col) or "").strip()
            if not val:
                continue
            out.setdefault(val, grp)
            for part in val.split("|"):
                part = part.strip()
                if not part:
                    continue
                out.setdefault(part, grp)
                head = part.split("_")[0]
                if head and head not in {"cand", "ICTV"}:
                    out.setdefault(head, grp)
    return out


def lookup_host(genome: str, host_map: dict[str, str], family: str = "") -> str:
    if genome in host_map:
        return host_map[genome]
    if genome.startswith("ICTV|"):
        return "ICTV Vertebrate Ref" if family == "Hepeviridae" else "ICTV Reference"
    for part in genome.split("|"):
        if part in host_map:
            return host_map[part]
        head = part.split("_")[0]
        if head in RUN_TO_BIOPROJECT:
            return BIOPROJECT_GROUP[RUN_TO_BIOPROJECT[head]]
        if head in host_map:
            return host_map[head]
    # Check conifer host tokens directly in genome ID if metadata was partial
    broad = _classify_broad_host(genome)
    if broad in {"Ips confusus", "Hylobius abietis", "Locusta migratoria"}:
        return broad
    # GenBank accessions (not ERR/SRR/DRR) default to ICTV Ref
    head = genome.split("_")[0]
    if head and not head.startswith(("ERR", "SRR", "DRR", "cand")):
        return "ICTV Vertebrate Ref" if family == "Hepeviridae" else "ICTV Reference"
    return "Other / Unknown"


def panel_sites(ax, sites: list[dict], family: str = "", locus: str = "ORF1") -> int:
    # 2: omega per site; log axis, because purifying sites crowd below 1
    drawn = 0
    max_site = 1
    for cat in CAT_ORDER:
        xs = [int(r["site"]) for r in sites
              if r.get("selection_category") == cat and as_float(r.get("omega")) > 0]
        ys = [as_float(r["omega"]) for r in sites
              if r.get("selection_category") == cat and as_float(r.get("omega")) > 0]
        if not xs:
            continue
        max_site = max(max_site, max(xs))
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
    if family == "Hepeviridae" and max_site >= 300:
        domains = HEPE_ORF2_DOMAINS if locus.upper() == "ORF2" else HEPE_ORF1_DOMAINS
        for dname, f0, f1, bg in domains:
            x0, x1 = f0 * max_site, f1 * max_site
            ax.axvspan(x0, x1, ymin=0.90, ymax=1.0, color=bg, alpha=0.85, zorder=1)
            ax.text((x0 + x1) / 2.0, 0.95, dname, transform=ax.get_xaxis_transform(),
                    ha="center", va="center", fontsize=6.8, fontweight="bold", color=INK)
    ax.set_xlabel(f"codon position in {locus}", fontsize=8.5, color=INK)
    ax.set_ylabel("omega (dN/dS)", fontsize=8.5, color=INK)
    ax.set_title(f"A  Site-level selection along {locus}", loc="left",
                 fontsize=9.5, fontweight="bold", color=INK)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.0, colors=INK_SOFT)
    if drawn:
        ax.legend(frameon=False, fontsize=7.6, labelcolor=INK_SOFT,
                  loc="lower left", ncol=4, handletextpad=0.3, borderaxespad=0.2)
    return drawn


def panel_cpg_vs_upa(ax, dinuc: list[dict], host_map: dict[str, str]) -> None:
    """Host-stratified CpG vs UpA odds ratios to separate vertebrate vs environmental signatures."""
    groups_seen: dict[str, list[tuple[float, float]]] = {}
    for r in dinuc:
        cg, ta = as_float(r.get("CG")), as_float(r.get("TA"))
        if not (cg > 0 and ta > 0):
            continue
        grp = lookup_host(r.get("genome", ""), host_map)
        groups_seen.setdefault(grp, []).append((cg, ta))
    for grp, pts in groups_seen.items():
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.scatter(xs, ys, s=22, c=HOST_PALETTE.get(grp, "#94A3B8"),
                   label=f"{grp} (n={len(pts)})", alpha=0.82,
                   edgecolors="white", linewidths=0.4, zorder=3)
    ax.axvline(1.0, color=INK_SOFT, lw=0.9, ls="--", zorder=1)
    ax.axhline(1.0, color=INK_SOFT, lw=0.9, ls="--", zorder=1)
    ax.set_xlabel("CpG odds ratio (1.0 = independence)", fontsize=8.5, color=INK)
    ax.set_ylabel("UpA odds ratio (1.0 = independence)", fontsize=8.5, color=INK)
    ax.set_title("B  Host-stratified CpG vs. UpA depletion", loc="left",
                 fontsize=9.5, fontweight="bold", color=INK)
    ax.grid(color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.0, colors=INK_SOFT)
    ax.legend(frameon=False, fontsize=7.0, labelcolor=INK_SOFT, loc="best", handletextpad=0.3)


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


def panel_composition(ax, dinuc: list[dict], host_map: dict[str, str] | None = None) -> None:
    # 4: odds ratios per genome, coloured by host/source when metadata is available
    keys = [k for k in ("CG", "TA", "GC") if dinuc and k in dinuc[0]]
    if not keys:
        ax.text(0.5, 0.5, "no dinucleotide table", ha="center", va="center",
                transform=ax.transAxes, fontsize=8.5, color=INK_SOFT)
        ax.set_axis_off()
        return
    rng = np.random.default_rng(0)  # fixed, so the jitter is reproducible
    for i, key in enumerate(keys):
        rows = [r for r in dinuc if as_float(r.get(key)) > 0]
        if not rows:
            continue
        vals = [as_float(r[key]) for r in rows]
        cols = [
            HOST_PALETTE.get(lookup_host(r.get("genome", ""), host_map), DINUC_COLOUR[key])
            if host_map else DINUC_COLOUR[key]
            for r in rows
        ]
        jitter = rng.uniform(-0.14, 0.14, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=16,
                   c=cols, alpha=0.72, edgecolors="white", linewidths=0.3, zorder=2)
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
    p.add_argument("--metadata", help="optional taxa_metadata.tsv or selected.tsv for host stratification")
    p.add_argument("--locus", default="ORF1", help="Target locus label (e.g. ORF1, ORF2)")
    p.add_argument("--family", required=True)
    p.add_argument("--fig-dir", required=True)
    p.add_argument("--stem", help="default signatures_<family>")
    args = p.parse_args()

    sites = read_tsv(args.sites)
    dinuc = read_tsv(args.dinucleotide)
    rscu = read_tsv(args.rscu) if args.rscu else []
    meta_rows = read_tsv(args.metadata) if args.metadata else []
    host_map = build_host_map(meta_rows, family=args.family) if meta_rows else {}
    if not sites:
        log(f"[plot] {args.sites} has no rows; nothing to draw")
        return 1
    log(f"[plot] {len(sites)} site(s), {len(dinuc)} genome(s) with composition")

    fig = plt.figure(figsize=(13.6, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], width_ratios=[1.0, 1.0])
    panel_sites(fig.add_subplot(gs[0, :]), sites, family=args.family, locus=args.locus)
    if host_map and dinuc and "CG" in dinuc[0] and "TA" in dinuc[0]:
        panel_cpg_vs_upa(fig.add_subplot(gs[1, 0]), dinuc, host_map)
    else:
        panel_counts(fig.add_subplot(gs[1, 0]), sites)
    panel_composition(fig.add_subplot(gs[1, 1]), dinuc, host_map=host_map or None)

    median_omega = ""
    for r in read_tsv(args.summary) if args.summary else []:
        if r.get("metric") == "median_omega":
            median_omega = f", median omega {as_float(r['value']):.3f}"
    fig.suptitle(rf"$\it{{{args.family}}}$ ({args.locus}): selection and host-composition signatures",
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
