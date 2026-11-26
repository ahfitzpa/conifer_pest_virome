#!/usr/bin/env python3
"""Generate Figure 3: Snowball replication & pan-Holarctic cross-host sharing matrix.

Reads `snowball_vs_seed_blastx.tsv` (from `snowball/Snakefile_snowball`),
`esm_rdrp_validation.tsv` (for original lineage family labels), and
`manifest_combined.csv` (for snowball run host species & BioProject metadata),
producing a clean publication-ready SVG dot-matrix separating:
  - >= 90% aa identity: independent same-virus / strain replications across BioProjects
  - 40-89% aa identity: full-length trans-continental sister lineages across beetle species
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import sys
import tempfile
from collections import defaultdict
from html import escape
from pathlib import Path

__version__ = "1.0.0"
LOG = logging.getLogger("plot_fig3_snowball_matrix")


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


def load_run_metadata(manifest_csv: Path | None) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    if manifest_csv is None or not manifest_csv.is_file():
        return meta
    with manifest_csv.open("r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rid = (r.get("run_id") or r.get("sample_id") or "").strip()
            if not rid:
                continue
            meta[rid] = {
                "host": (r.get("scientificname") or r.get("host_species") or "Unknown").strip().replace("_", " "),
                "bioproject": (r.get("bioproject") or "").strip(),
            }
    return meta


def load_seed_labels(esm_tsv: Path | None) -> dict[str, str]:
    labels: dict[str, str] = {}
    if esm_tsv is None or not esm_tsv.is_file():
        return labels
    with esm_tsv.open("r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cid = (r.get("candidate_id") or "").strip()
            fam = (r.get("consensus_family") or "Unresolved").strip()
            tgt = (r.get("rdrp_target") or "").strip()
            if tgt in ("8gzq-assembly1_B", "8gzr-assembly1_B"):
                fam = "Viral NS3 helicase"
            short = cid.split("__", 1)[-1] if "__" in cid else cid
            labels[cid] = f"{fam} ({short})"
    return labels


def identity_color(pident: float) -> str:
    if pident >= 95.0:
        return "#0F766E"  # exact strain replication: deep teal
    if pident >= 90.0:
        return "#0D9488"  # same virus replication: teal
    if pident >= 80.0:
        return "#2563EB"  # close sister species: royal blue
    if pident >= 60.0:
        return "#7C3AED"  # divergent sister lineage: violet
    return "#D97706"      # deep homolog: amber


def parse_blastx_matrix(
    hits_tsv: Path,
    min_pident: float = 40.0,
    min_alnlen: int = 100,
) -> dict[tuple[str, str], dict[str, float]]:
    """Aggregate best hit per (seed_id, snowball_run_id)."""
    cells: dict[tuple[str, str], dict[str, float]] = {}
    if not hits_tsv.is_file() or hits_tsv.stat().st_size == 0:
        return cells
    with hits_tsv.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("qseqid\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            qseqid, sseqid = parts[0].strip(), parts[1].strip()
            try:
                pident = float(parts[2])
                alnlen = float(parts[3])
                qlen = float(parts[4])
                bitscore = float(parts[7])
            except ValueError:
                continue
            if pident < min_pident or alnlen < min_alnlen:
                continue
            run_id = qseqid.split("|", 1)[0] if "|" in qseqid else qseqid.split("_", 1)[0]
            # Extract contig length from SPAdes header if qlen is missing
            if qlen <= 0:
                m = re.search(r"length_(\d+)", qseqid)
                if m:
                    qlen = float(m.group(1))
            key = (sseqid, run_id)
            prev = cells.get(key)
            if prev is None or (pident, bitscore) > (prev["pident"], prev["bitscore"]):
                cells[key] = {
                    "pident": pident,
                    "alnlen": alnlen,
                    "qlen_kb": qlen / 1000.0,
                    "bitscore": bitscore,
                }
    return cells


def build_svg(
    cells: dict[tuple[str, str], dict[str, float]],
    run_meta: dict[str, dict[str, str]],
    seed_labels: dict[str, str],
) -> str:
    seeds = sorted(
        {k[0] for k in cells},
        key=lambda s: (
            -max(cells[(s, r)]["pident"] for (ss, r) in cells if ss == s),
            -sum(1 for (ss, _) in cells if ss == s),
            s,
        ),
    )
    runs = sorted(
        {k[1] for k in cells},
        key=lambda r: (run_meta.get(r, {}).get("host", "ZZZ"), run_meta.get(r, {}).get("bioproject", ""), r),
    )

    left_pad = 300
    top_pad = 128
    col_w = 56
    row_h = 34
    right_pad = 40
    bot_pad = 80

    W = max(920, left_pad + max(len(runs), 6) * col_w + right_pad)
    H = max(360, top_pad + max(len(seeds), 4) * row_h + bot_pad)

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<defs>',
        '  <style>',
        '    text { font-family: "Inter", "Helvetica Neue", Arial, sans-serif; fill: #0F172A; }',
        '    .title { font-size: 15px; font-weight: 700; }',
        '    .row-lbl { font-size: 11.5px; font-weight: 600; fill: #1E293B; }',
        '    .col-lbl { font-size: 10.5px; font-weight: 600; fill: #334155; }',
        '    .host-hdr { font-size: 11px; font-weight: 700; font-style: italic; fill: #0F172A; }',
        '    .cell-txt { font-size: 9px; font-weight: 700; }',
        '    .leg { font-size: 11px; fill: #334155; }',
        '  </style>',
        '</defs>',
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="#FFFFFF"/>',
        '<text x="24" y="30" class="title">Snowball replication &amp; pan-Holarctic cross-host virus sharing (DIAMOND blastx)</text>',
    ]

    if not seeds or not runs:
        svg.append('<text x="24" y="80" class="row-lbl">No snowball cross-match hits passing thresholds.</text>')
        svg.append('</svg>')
        return "\n".join(svg)

    # Top host grouping brackets
    groups: list[tuple[str, int, int]] = []
    cur_h, start_i = None, 0
    for idx, rid in enumerate(runs):
        h = run_meta.get(rid, {}).get("host", "Unknown")
        if h != cur_h:
            if cur_h is not None:
                groups.append((cur_h, start_i, idx - 1))
            cur_h, start_i = h, idx
    if cur_h is not None:
        groups.append((cur_h, start_i, len(runs) - 1))

    for hname, i0, i1 in groups:
        gx0 = left_pad + i0 * col_w + 6
        gx1 = left_pad + (i1 + 1) * col_w - 6
        svg.append(f'<line x1="{gx0}" y1="{top_pad-58}" x2="{gx1}" y2="{top_pad-58}" stroke="#334155" stroke-width="1.5"/>')
        svg.append(f'<text x="{(gx0+gx1)/2:.1f}" y="{top_pad-66}" text-anchor="middle" class="host-hdr">{escape(hname)}</text>')

    # Column run labels
    for j, rid in enumerate(runs):
        cx = left_pad + j * col_w + col_w / 2.0
        svg.append(f'<text x="{cx:.1f}" y="{top_pad-14}" text-anchor="end" class="col-lbl" transform="rotate(-40 {cx:.1f} {top_pad-14})">{escape(rid)}</text>')

    # Grid & dots
    for i, seed in enumerate(seeds):
        cy = top_pad + i * row_h + row_h / 2.0
        if i % 2 == 0:
            svg.append(f'<rect x="{left_pad}" y="{top_pad + i*row_h}" width="{len(runs)*col_w}" height="{row_h}" fill="#F8FAFC"/>')
        lbl = seed_labels.get(seed, seed)
        svg.append(f'<text x="{left_pad-12}" y="{cy+4:.1f}" text-anchor="end" class="row-lbl">{escape(lbl)}</text>')

        for j, rid in enumerate(runs):
            cx = left_pad + j * col_w + col_w / 2.0
            hit = cells.get((seed, rid))
            if not hit:
                svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="2.2" fill="#E2E8F0"/>')
                continue
            pid = hit["pident"]
            qkb = hit["qlen_kb"]
            r = max(6.5, min(14.5, 5.5 + qkb * 0.72))
            col = identity_color(pid)
            svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{col}" stroke="#FFFFFF" stroke-width="1.2"/>')
            svg.append(f'<text x="{cx:.1f}" y="{cy+3.2:.1f}" text-anchor="middle" class="cell-txt" style="fill:#FFFFFF">{pid:.0f}</text>')

    # Bottom legend
    ly = top_pad + len(seeds) * row_h + 38
    leg_items = [
        ("#0F766E", "≥95% aa ID (same strain)"),
        ("#0D9488", "90–95% (same virus)"),
        ("#2563EB", "80–90% (close sister)"),
        ("#7C3AED", "60–80% (sister lineage)"),
        ("#D97706", "40–60% (distant homolog)"),
    ]
    lx = 24
    for col, txt in leg_items:
        svg.append(f'<circle cx="{lx+6}" cy="{ly-4}" r="6" fill="{col}"/>')
        svg.append(f'<text x="{lx+18}" y="{ly}" class="leg">{escape(txt)}</text>')
        lx += len(txt) * 6.5 + 34
    svg.append(f'<text x="24" y="{ly+22}" class="leg" style="font-style:italic;fill:#64748B;">Dot radius scales with assembled contig length (1.5–12.4 kb); numbers inside dots show % amino-acid identity.</text>')

    svg.append('</svg>')
    return "\n".join(svg)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hits-tsv", required=True, type=Path, help="snowball_vs_seed_blastx.tsv")
    ap.add_argument("--manifest", default=None, type=Path, help="manifest_combined.csv")
    ap.add_argument("--esm-validation", default=None, type=Path, help="esm_rdrp_validation.tsv")
    ap.add_argument("--min-pident", type=float, default=40.0)
    ap.add_argument("--min-alnlen", type=int, default=100)
    ap.add_argument("--out-svg", required=True, type=Path)
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    cells = parse_blastx_matrix(args.hits_tsv, args.min_pident, args.min_alnlen)
    run_meta = load_run_metadata(args.manifest)
    seed_labels = load_seed_labels(args.esm_validation)
    svg = build_svg(cells, run_meta, seed_labels)
    atomic_write_text(args.out_svg, svg)
    LOG.info("Wrote Fig. 3 snowball sharing matrix (%d hit cells) -> %s", len(cells), args.out_svg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
