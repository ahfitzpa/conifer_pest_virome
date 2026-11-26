#!/usr/bin/env python3
"""Generate Figure 1: Closed-loop petabase discovery workflow & Holarctic sampling map.

Pure-Python zero-dependency SVG generator (no matplotlib/OpenBLAS required).
Reads `host_geography.tsv` and optional `manifest_combined.csv` to plot exact
per-species / per-BioProject library counts across North America, Europe, and Asia.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import tempfile
from collections import Counter, defaultdict
from html import escape
from pathlib import Path

__version__ = "1.0.0"
LOG = logging.getLogger("plot_fig1_workflow_and_map")


def atomic_write_text(dest: Path, text: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=dest.parent, prefix=f".{dest.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp_path, dest)
    except BaseException:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def load_host_geography(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def count_manifest_runs(manifest_csv: Path | None) -> tuple[Counter[str], Counter[str]]:
    """Return (by_host, by_bioproject) counts from manifest_combined.csv."""
    by_host: Counter[str] = Counter()
    by_bp: Counter[str] = Counter()
    if manifest_csv is None or not manifest_csv.is_file() or manifest_csv.stat().st_size == 0:
        return by_host, by_bp
    with manifest_csv.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            host = (row.get("scientificname") or row.get("host_species") or "").strip().replace("_", " ")
            bp = (row.get("bioproject") or "").strip()
            if host.lower() == "locusta migratoria":
                continue
            if host:
                by_host[host] += 1
            if bp:
                by_bp[bp] += 1
    return by_host, by_bp


def build_svg(geo_rows: list[dict[str, str]], manifest_csv: Path | None) -> str:
    by_host, by_bp = count_manifest_runs(manifest_csv)

    # Aggregate unique species metadata
    species_meta: dict[str, dict[str, object]] = {}
    for r in geo_rows:
        sp = r["host_species"].strip()
        if sp not in species_meta:
            species_meta[sp] = {
                "subfamily": r["subfamily"].strip(),
                "continent": r["continent"].strip(),
                "region": r["region"].strip(),
                "color": r["color_hex"].strip(),
                "x": float(r["map_x"]),
                "y": float(r["map_y"]),
                "bps": [],
                "stages": set(),
            }
        species_meta[sp]["bps"].append(r["bioproject"].strip())
        species_meta[sp]["stages"].add(r["stage"].strip())

    W, H = 1400, 680
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<defs>',
        '  <style>',
        '    text { font-family: "Inter", "Helvetica Neue", Arial, sans-serif; fill: #0F172A; }',
        '    .ptitle { font-size: 16px; font-weight: 700; }',
        '    .box-t { font-size: 13px; font-weight: 700; }',
        '    .box-s { font-size: 11px; fill: #334155; }',
        '    .badge { font-size: 10.5px; font-weight: 700; }',
        '    .sp-lbl { font-size: 12px; font-weight: 700; font-style: italic; }',
        '    .sub-lbl { font-size: 10.5px; fill: #475569; }',
        '  </style>',
        '  <marker id="arr-navy" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto">',
        '    <path d="M 0 1 L 10 5 L 0 9 z" fill="#1E293B"/>',
        '  </marker>',
        '  <marker id="arr-teal" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto">',
        '    <path d="M 0 1 L 10 5 L 0 9 z" fill="#0F766E"/>',
        '  </marker>',
        '</defs>',
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="#FFFFFF"/>',
    ]

    # =========================================================================
    # PANEL A: 4-Stage Closed-Loop Petabase Discovery Framework (top row)
    # =========================================================================
    svg.append('<text x="28" y="34" class="ptitle">a   Closed-loop protein- and structure-guided petabase virus discovery</text>')

    stages = [
        {
            "x": 28, "w": 300, "col": "#1D4ED8", "bg": "#EFF6FF", "bd": "#BFDBFE",
            "step": "STAGE 1 · SEED DISCOVERY",
            "title": "Serratus Protein-Space Index",
            "lines": [
                "• Translated RdRp palmprints surface",
                "  divergent beetle RNA viruses (25–52% aa ID)",
                "• Nucleotide k-mers (PebbleScout/MetaGraph)",
                "  return 0 discovery hits from ICTV references",
            ],
        },
        {
            "x": 366, "w": 300, "col": "#0F766E", "bg": "#F0FDFA", "bd": "#99F6E4",
            "step": "STAGE 2 · ASSEMBLY & TRIAGE",
            "title": "Phase 2 + ORF Verification",
            "lines": [
                "• fastp → rnaviralSPAdes → multi-tool",
                "  (geNomad, RdRpCATCH, Kraken2, MetaCerberus)",
                "• Direct 6-frame ORF completeness check",
                "  (replaces CheckV for segmented RNA)",
            ],
        },
        {
            "x": 704, "w": 310, "col": "#7C3AED", "bg": "#F5F3FF", "bd": "#DDD6FE",
            "step": "STAGE 3 · 3D STRUCTURE FILTER",
            "title": "ESMFold + Foldseek (A100 GPU)",
            "lines": [
                "• Folds gated <40% ID & unresolved domains",
                "• Validates 19 RdRps (TM 0.69–0.90) + 4 NS3",
                "• Purges 79.5% HMM false positives (38 RTs",
                "  + 78 host DEAH-box helicases / MTases)",
            ],
        },
        {
            "x": 1052, "w": 320, "col": "#D97706", "bg": "#FFFBEB", "bd": "#FDE68A",
            "step": "STAGE 4 · AUTOMATED SNOWBALL",
            "title": "Validated RdRp → Serratus Loop",
            "lines": [
                "• DIAMOND vs PalmDB → Serratus RDS query",
                "• Recovers 12 new Holarctic bark-beetle runs",
                "• 7 viruses replicated (93.5–100% aa ID)",
                "• 8 complete ~12.3 kb Ips Polycipiviridae",
            ],
        },
    ]

    sy, sh = 52, 166
    for i, st in enumerate(stages):
        x, w = st["x"], st["w"]
        svg.append(f'<rect x="{x}" y="{sy}" width="{w}" height="{sh}" rx="10" fill="{st["bg"]}" stroke="{st["bd"]}" stroke-width="1.5"/>')
        svg.append(f'<rect x="{x}" y="{sy}" width="{w}" height="28" rx="10" fill="{st["col"]}"/>')
        svg.append(f'<rect x="{x}" y="{sy+18}" width="{w}" height="10" fill="{st["col"]}"/>')
        svg.append(f'<text x="{x+14}" y="{sy+19}" class="badge" style="fill:#FFFFFF">{escape(st["step"])}</text>')
        svg.append(f'<text x="{x+14}" y="{sy+50}" class="box-t">{escape(st["title"])}</text>')
        for j, ln in enumerate(st["lines"]):
            svg.append(f'<text x="{x+14}" y="{sy+72 + j*19}" class="box-s">{escape(ln)}</text>')
        if i < len(stages) - 1:
            nx = stages[i + 1]["x"]
            svg.append(f'<line x1="{x+w+4}" y1="{sy+sh/2:.1f}" x2="{nx-6}" y2="{sy+sh/2:.1f}" stroke="#1E293B" stroke-width="2.0" marker-end="url(#arr-navy)"/>')

    # Closed feedback arc from Stage 4 back to Stage 2 (combined 24-library phylogenetics)
    svg.append(
        f'<path d="M 1212 {sy+sh+2} L 1212 {sy+sh+24} L 516 {sy+sh+24} L 516 {sy+sh+6}" '
        'fill="none" stroke="#0F766E" stroke-width="2.0" stroke-dasharray="6,4" marker-end="url(#arr-teal)"/>'
    )
    svg.append(
        f'<rect x="742" y="{sy+sh+12}" width="286" height="22" rx="11" fill="#FFFFFF" stroke="#99F6E4"/>'
    )
    svg.append(
        f'<text x="885" y="{sy+sh+27}" text-anchor="middle" class="badge" style="fill:#0F766E">'
        'Closed loop → Combined 24-library AA phylogenetics</text>'
    )

    # =========================================================================
    # PANEL B: Holarctic 3-Continent Conifer-Pest Sampling Map (bottom row)
    # =========================================================================
    my0 = 278
    svg.append(f'<text x="28" y="{my0}" class="ptitle">b   Pan-Holarctic conifer-feeding Curculionidae sampling (5 species · 3 continents)</text>')

    mx0, mw, mh = 28, 860, 370
    svg.append(f'<rect x="{mx0}" y="{my0+14}" width="{mw}" height="{mh}" rx="12" fill="#F8FAFC" stroke="#CBD5E1" stroke-width="1.3"/>')

    # Clean stylized Holarctic landmass silhouettes (North America, Europe, Asia)
    # Coordinates inside map box (offset by mx0, my0+14)
    ox, oy = mx0, my0 + 14
    landmasses = [
        # North America (boreal / western conifer belt)
        "M 60,95 L 135,60 L 235,68 L 285,115 L 255,175 L 210,225 L 155,245 L 115,205 L 80,150 Z",
        # Europe
        "M 390,85 L 460,55 L 535,70 L 545,145 L 495,180 L 420,175 L 385,135 Z",
        # Asia / Siberia / Qinghai-Tibet / East Asia
        "M 545,70 L 660,52 L 785,78 L 820,145 L 775,215 L 705,240 L 615,220 L 555,175 Z",
    ]
    for path_d in landmasses:
        svg.append(f'<g transform="translate({ox},{oy})"><path d="{path_d}" fill="#E2E8F0" stroke="#94A3B8" stroke-width="1.3"/></g>')

    # Continent headers
    svg.append(f'<text x="{ox+165}" y="{oy+42}" text-anchor="middle" class="badge" style="fill:#475569">NORTH AMERICA</text>')
    svg.append(f'<text x="{ox+465}" y="{oy+42}" text-anchor="middle" class="badge" style="fill:#475569">EUROPE</text>')
    svg.append(f'<text x="{ox+690}" y="{oy+42}" text-anchor="middle" class="badge" style="fill:#475569">ASIA</text>')

    # Snowball connection arcs showing trans-continental virus sharing
    # I. confusus (185,165) -> I. typographus (480,138) -> I. nitidus (675,175) -> I. subelongatus (725,152)
    arcs = [
        (185, 165, 480, 138, "#1D4ED8", "Polycipiviridae (88.9% aa ID)"),
        (480, 138, 675, 175, "#D97706", "Polycipiviridae & Spinareoviridae"),
        (185, 165, 440, 120, "#0F766E", "Iflaviridae & NS3 (61–63%)"),
    ]
    for x1, y1, x2, y2, col, _ in arcs:
        cx, cy = (x1 + x2) / 2.0, min(y1, y2) - 42
        svg.append(
            f'<g transform="translate({ox},{oy})">'
            f'<path d="M {x1},{y1} Q {cx:.1f},{cy:.1f} {x2},{y2}" fill="none" stroke="{col}" '
            f'stroke-width="2.0" stroke-dasharray="5,4" opacity="0.75"/></g>'
        )

    # Species pins
    for sp, meta in species_meta.items():
        px, py = ox + float(meta["x"]), oy + float(meta["y"])
        col = str(meta["color"])
        n_runs = by_host.get(sp, 0)
        is_seed = "seed" in meta["stages"]
        r = 11 if is_seed else 9
        if is_seed:
            svg.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{r+4}" fill="{col}" fill-opacity="0.18"/>')
        svg.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{r}" fill="{col}" stroke="#FFFFFF" stroke-width="2"/>')
        lbl_y = py + 24 if sp != "Ips typographus" else py + 24
        tag = f" ({n_runs} runs)" if n_runs else ""
        svg.append(f'<text x="{px:.1f}" y="{lbl_y:.1f}" text-anchor="middle" class="sp-lbl">{escape(sp)}</text>')
        if tag:
            svg.append(f'<text x="{px:.1f}" y="{lbl_y+14:.1f}" text-anchor="middle" class="sub-lbl">{escape(tag.strip())}</text>')

    # Map legend bottom-left
    svg.append(f'<rect x="{ox+16}" y="{oy+mh-48}" width="405" height="34" rx="8" fill="#FFFFFF" stroke="#CBD5E1"/>')
    svg.append(f'<circle cx="{ox+34}" cy="{oy+mh-31}" r="7" fill="#1D4ED8" stroke="#FFFFFF" stroke-width="1.5"/>')
    svg.append(f'<text x="{ox+48}" y="{oy+mh-27}" class="sub-lbl">Seed hosts (12 libraries)</text>')
    svg.append(f'<circle cx="{ox+205}" cy="{oy+mh-31}" r="6" fill="#D97706" stroke="#FFFFFF" stroke-width="1.5"/>')
    svg.append(f'<text x="{ox+218}" y="{oy+mh-27}" class="sub-lbl">Snowball hosts (+12 libraries)</text>')

    # Right summary table card (Panel b sidebar)
    tx, tw = mx0 + mw + 22, W - (mx0 + mw + 22) - 28
    svg.append(f'<rect x="{tx}" y="{my0+14}" width="{tw}" height="{mh}" rx="12" fill="#F8FAFC" stroke="#CBD5E1" stroke-width="1.3"/>')
    svg.append(f'<text x="{tx+18}" y="{my0+42}" class="box-t">Dataset Expansion Summary</text>')
    svg.append(f'<line x1="{tx+16}" y1="{my0+52}" x2="{tx+tw-16}" y2="{my0+52}" stroke="#CBD5E1"/>')

    yrow = my0 + 78
    for sp, meta in species_meta.items():
        col = str(meta["color"])
        bps = ", ".join(dict.fromkeys(meta["bps"]))
        stage_lbl = "Seed + Snowball" if len(meta["stages"]) > 1 else ("Seed" if "seed" in meta["stages"] else "Snowball")
        n_runs = by_host.get(sp, 0)
        run_txt = f"{n_runs} runs" if n_runs else f"{len(meta['bps'])} BP(s)"
        svg.append(f'<circle cx="{tx+24}" cy="{yrow-4}" r="6" fill="{col}"/>')
        svg.append(f'<text x="{tx+38}" y="{yrow}" class="sp-lbl">{escape(sp)}</text>')
        svg.append(f'<text x="{tx+tw-18}" y="{yrow}" text-anchor="end" class="badge" style="fill:{col}">{escape(stage_lbl)}</text>')
        svg.append(f'<text x="{tx+38}" y="{yrow+16}" class="sub-lbl">{escape(str(meta["continent"]))} · {escape(str(meta["subfamily"]))} · {escape(run_txt)}</text>')
        svg.append(f'<text x="{tx+38}" y="{yrow+30}" class="sub-lbl" style="font-family:monospace;font-size:9.8px;">{escape(bps)}</text>')
        yrow += 58

    svg.append('</svg>')
    return "\n".join(svg)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--host-geo", required=True, type=Path, help="Path to configs/host_geography.tsv")
    ap.add_argument("--manifest", default=None, type=Path, help="Optional combined Phase 2 manifest CSV")
    ap.add_argument("--out-svg", required=True, type=Path, help="Output SVG path")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    geo_rows = load_host_geography(args.host_geo)
    svg_text = build_svg(geo_rows, args.manifest)
    atomic_write_text(args.out_svg, svg_text)
    LOG.info("Wrote Fig. 1 workflow + Holarctic map -> %s", args.out_svg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
