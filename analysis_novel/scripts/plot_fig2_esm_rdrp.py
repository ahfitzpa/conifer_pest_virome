#!/usr/bin/env python3
"""Generate Figure 2: Stage 15 ESMFold + Foldseek structural validation & 3D RdRp folds.

Pure-Python zero-dependency SVG generator (no matplotlib/PyMOL/OpenBLAS required).
Reads `results/esm_rdrp/esm_rdrp_validation.tsv` and optional showcase `.pdb`
structures from `results/esm_rdrp/structures/`, writing:
  - `--out-validation-svg`: 3-panel structural validation space, RT negative-control
                            diagonal, and host triage breakdown
  - `--out-3d-svg` (optional): 3-panel 3D pLDDT-coloured RdRp backbone gallery with
                               catalytic Motif A & Motif C highlights
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import re
import sys
import tempfile
from collections import Counter
from html import escape
from pathlib import Path

__version__ = "1.0.0"
LOG = logging.getLogger("plot_fig2_esm_rdrp")

COLORS = {
    "Validated RdRp (high)":    "#0F766E",
    "Validated RdRp (mod.)":    "#38BDF8",
    "Viral NS3 helicase":       "#7C3AED",
    "RT / retroelement flag":   "#DC2626",
    "Host helicase / non-RdRp": "#94A3B8",
    "Low pLDDT / no hit":       "#CBD5E1",
}
ORDER = [
    "Validated RdRp (high)",
    "Validated RdRp (mod.)",
    "Viral NS3 helicase",
    "RT / retroelement flag",
    "Host helicase / non-RdRp",
    "Low pLDDT / no hit",
]

AA3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


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


def safe_float(v: str | None, default: float = 0.0) -> float:
    try:
        return float(v) if v else default
    except ValueError:
        return default


def host_of(run_id: str, manifest_hosts: dict[str, str]) -> str | None:
    if run_id in manifest_hosts:
        h = manifest_hosts[run_id]
        return None if h.lower() == "locusta migratoria" else h
    if run_id.startswith("SRR122616"):
        return "Ips confusus"
    if run_id.startswith(("SRR67659", "SRR170354", "SRR63238")):
        return "Hylobius abietis"
    if run_id.startswith("SRR98986"):
        return None
    return "Other Coleoptera"


def clean_category(row: dict[str, str]) -> str:
    v = (row.get("structural_verdict") or "").strip()
    tgt = (row.get("rdrp_target") or "").strip()
    if v in ("validated_rdrp_high", "validated_rdrp_moderate") and tgt in ("8gzq-assembly1_B", "8gzr-assembly1_B"):
        return "Viral NS3 helicase"
    if v == "validated_rdrp_high":
        return "Validated RdRp (high)"
    if v == "validated_rdrp_moderate":
        return "Validated RdRp (mod.)"
    if v == "rt_or_non_rdrp_flag":
        return "RT / retroelement flag"
    if v == "non_rdrp_structural_hit":
        return "Host helicase / non-RdRp"
    return "Low pLDDT / no hit"


def load_manifest_hosts(manifest_csv: Path | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if manifest_csv is None or not manifest_csv.is_file():
        return out
    with manifest_csv.open("r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rid = (r.get("run_id") or r.get("sample_id") or "").strip()
            h = (r.get("scientificname") or r.get("host_species") or "").strip().replace("_", " ")
            if rid and h:
                out[rid] = h
    return out


def build_validation_svg(validation_tsv: Path, manifest_hosts: dict[str, str]) -> str:
    records = []
    with validation_tsv.open("r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            h = host_of(r.get("run_id", ""), manifest_hosts)
            if h is None:
                continue
            records.append({
                "cid": r.get("candidate_id", ""),
                "host": h,
                "cat": clean_category(r),
                "plddt": safe_float(r.get("mean_plddt")),
                "tm": safe_float(r.get("rdrp_tmscore")),
                "rt_tm": safe_float(r.get("rt_tmscore")),
                "aa_len": int(safe_float(r.get("folded_aa_len"), 300)),
                "unres": (r.get("consensus_family") or "").strip().lower() in ("", "unresolved", "unclassified"),
            })

    n_total = len(records)
    counts = Counter(r["cat"] for r in records)
    W, H = 1420, 590

    def map_x(val, vmin, vmax, x0, x1):
        return x0 + (val - vmin) / (vmax - vmin) * (x1 - x0)

    def map_y(val, vmin, vmax, y0, y1):
        return y1 - (val - vmin) / (vmax - vmin) * (y1 - y0)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<defs><style>',
        '  text { font-family: "Inter", "Helvetica Neue", Arial, sans-serif; fill: #0F172A; }',
        '  .title { font-size: 15px; font-weight: 700; }',
        '  .axis-label { font-size: 13px; font-weight: 600; fill: #1E293B; }',
        '  .tick { font-size: 11px; fill: #475569; }',
        '  .annot { font-size: 11px; font-weight: 600; }',
        '  .annot-sub { font-size: 10.5px; font-style: italic; }',
        '  .legend { font-size: 11.5px; font-weight: 600; fill: #1E293B; }',
        '</style></defs>',
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="#FFFFFF"/>',
    ]

    lx, ly = 55, 34
    for cat in ORDER:
        label = f"{cat} ({counts[cat]})"
        svg.append(f'<rect x="{lx}" y="{ly-10}" width="13" height="13" rx="3" fill="{COLORS[cat]}"/>')
        svg.append(f'<text x="{lx+18}" y="{ly+1}" class="legend">{escape(label)}</text>')
        lx += len(label) * 7.2 + 34

    # Panel A
    ax0, ax1, ay0, ay1 = 74, 540, 92, 505
    svg.append(f'<text x="{ax0-20}" y="{ay0-18}" class="title">a   ESMFold + Foldseek validation (Coleoptera, n = {n_total})</text>')
    zx0 = map_x(0.50, 0.0, 1.0, ax0, ax1)
    zy1 = map_y(70.0, 20.0, 100.0, ay0, ay1)
    svg.append(f'<rect x="{zx0:.1f}" y="{ay0}" width="{ax1-zx0:.1f}" height="{zy1-ay0:.1f}" fill="#CCFBF1" opacity="0.38"/>')

    for yv in [20, 40, 50, 60, 70, 80, 100]:
        yy = map_y(yv, 20, 100, ay0, ay1)
        if yv == 70:
            svg.append(f'<line x1="{ax0}" y1="{yy:.1f}" x2="{ax1}" y2="{yy:.1f}" stroke="#0F766E" stroke-width="1.2" stroke-dasharray="6,4"/>')
        elif yv == 50:
            svg.append(f'<line x1="{ax0}" y1="{yy:.1f}" x2="{ax1}" y2="{yy:.1f}" stroke="#64748B" stroke-width="0.9" stroke-dasharray="3,4"/>')
        else:
            svg.append(f'<line x1="{ax0}" y1="{yy:.1f}" x2="{ax1}" y2="{yy:.1f}" stroke="#F1F5F9" stroke-width="1.0"/>')
        svg.append(f'<text x="{ax0-10}" y="{yy+4:.1f}" text-anchor="end" class="tick">{yv}</text>')

    for xv in [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]:
        xx = map_x(xv, 0.0, 1.0, ax0, ax1)
        if abs(xv - 0.5) < 1e-6:
            svg.append(f'<line x1="{xx:.1f}" y1="{ay0}" x2="{xx:.1f}" y2="{ay1}" stroke="#0F766E" stroke-width="1.2" stroke-dasharray="6,4"/>')
        svg.append(f'<text x="{xx:.1f}" y="{ay1+20}" text-anchor="middle" class="tick">{xv:.1f}</text>')

    svg.append(f'<line x1="{ax0}" y1="{ay0}" x2="{ax0}" y2="{ay1}" stroke="#334155" stroke-width="1.3"/>')
    svg.append(f'<line x1="{ax0}" y1="{ay1}" x2="{ax1}" y2="{ay1}" stroke="#334155" stroke-width="1.3"/>')

    for cat in reversed(ORDER):
        for rec in records:
            if rec["cat"] != cat:
                continue
            cx = map_x(min(rec["tm"], 0.99), 0.0, 1.0, ax0, ax1)
            cy = map_y(max(21.0, min(rec["plddt"], 99.0)), 20.0, 100.0, ay0, ay1)
            r = 3.5 + (min(rec["aa_len"], 800) - 200) / 170.0
            col = COLORS[cat]
            op = "0.92" if cat.startswith("Validated") or cat == "Viral NS3 helicase" else "0.68"
            if rec["unres"]:
                svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{col}" fill-opacity="{op}" stroke="#FFFFFF" stroke-width="0.7"/>')
            else:
                pts = f"{cx:.1f},{cy-r*1.25:.1f} {cx-r*1.1:.1f},{cy+r*0.95:.1f} {cx+r*1.1:.1f},{cy+r*0.95:.1f}"
                svg.append(f'<polygon points="{pts}" fill="{col}" fill-opacity="{op}" stroke="#FFFFFF" stroke-width="0.7"/>')

    svg.append(f'<text x="{(ax0+ax1)/2:.1f}" y="{ay1+44}" text-anchor="middle" class="axis-label">Best Foldseek Structural Alignment TM-score</text>')
    svg.append(f'<text x="24" y="{(ay0+ay1)/2:.1f}" text-anchor="middle" class="axis-label" transform="rotate(-90 24 {(ay0+ay1)/2:.1f})">ESMFold Mean pLDDT</text>')

    # Panel B
    bx0, bx1, by0, by1 = 626, 976, 92, 505
    svg.append(f'<text x="{bx0-16}" y="{by0-18}" class="title">b   RdRp vs. RT negative-control fold</text>')
    svg.append(f'<polygon points="{bx0},{by0} {bx0},{map_y(0.05, 0.0, 1.0, by0, by1):.1f} {map_x(0.95, 0.0, 1.0, bx0, bx1):.1f},{by0}" fill="#FEE2E2" opacity="0.42"/>')
    svg.append(f'<line x1="{bx0}" y1="{by1}" x2="{bx1}" y2="{by0}" stroke="#475569" stroke-width="1.2" stroke-dasharray="5,4"/>')
    for tv in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        xx = map_x(tv, 0.0, 1.0, bx0, bx1)
        yy = map_y(tv, 0.0, 1.0, by0, by1)
        svg.append(f'<text x="{xx:.1f}" y="{by1+20}" text-anchor="middle" class="tick">{tv:.1f}</text>')
        svg.append(f'<text x="{bx0-10}" y="{yy+4:.1f}" text-anchor="end" class="tick">{tv:.1f}</text>')
    svg.append(f'<line x1="{bx0}" y1="{by0}" x2="{bx0}" y2="{by1}" stroke="#334155" stroke-width="1.3"/>')
    svg.append(f'<line x1="{bx0}" y1="{by1}" x2="{bx1}" y2="{by1}" stroke="#334155" stroke-width="1.3"/>')

    for cat in reversed(ORDER):
        for rec in records:
            if rec["cat"] != cat:
                continue
            cx = map_x(min(rec["tm"], 0.99), 0.0, 1.0, bx0, bx1)
            cy = map_y(min(rec["rt_tm"], 0.99), 0.0, 1.0, by0, by1)
            svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="3.8" fill="{COLORS[cat]}" fill-opacity="0.85" stroke="#FFFFFF" stroke-width="0.6"/>')

    svg.append(f'<text x="{(bx0+bx1)/2:.1f}" y="{by1+44}" text-anchor="middle" class="axis-label">Best Foldseek RdRp / PDB TM-score</text>')
    svg.append(f'<text x="{bx0-44}" y="{(by0+by1)/2:.1f}" text-anchor="middle" class="axis-label" transform="rotate(-90 {bx0-44} {(by0+by1)/2:.1f})">Foldseek RT Control TM-score</text>')

    # Panel C: Stacked bars by host
    cx0, cx1, cy0, cy1 = 1055, 1385, 92, 505
    svg.append(f'<text x="{cx0-18}" y="{cy0-18}" class="title">c   Structural triage by beetle host</text>')
    host_list = [h for h, _ in Counter(r["host"] for r in records).most_common(5)]
    bar_h = min(68, int(240 / max(1, len(host_list))))
    step = bar_h + 44
    for idx, hk in enumerate(host_list):
        by = cy0 + 42 + idx * step
        sub = [r for r in records if r["host"] == hk]
        tot = len(sub)
        svg.append(f'<text x="{cx0}" y="{by-10}" class="annot" font-style="italic">{escape(hk)} <tspan font-style="normal" font-weight="400" fill="#475569">— n = {tot}</tspan></text>')
        xcur = float(cx0)
        for cat in ORDER:
            cnt = sum(1 for r in sub if r["cat"] == cat)
            if cnt == 0:
                continue
            pct = 100.0 * cnt / tot
            bw = (cx1 - cx0) * (pct / 100.0)
            svg.append(f'<rect x="{xcur:.1f}" y="{by}" width="{bw:.1f}" height="{bar_h}" fill="{COLORS[cat]}" stroke="#FFFFFF" stroke-width="1.2"/>')
            if pct >= 10.0:
                tcol = "#FFFFFF" if cat in ("Validated RdRp (high)", "RT / retroelement flag", "Viral NS3 helicase") else "#0F172A"
                svg.append(f'<text x="{xcur + bw/2:.1f}" y="{by + bar_h/2 + 4:.1f}" text-anchor="middle" font-size="10.5" font-weight="700" style="fill:{tcol}">{cnt} ({pct:.0f}%)</text>')
            xcur += bw

    svg.append('</svg>')
    return "\n".join(svg)


def _load_ca(pdb_path: Path):
    pts = []
    seq = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
            resnum = int(line[22:26])
            aa = AA3TO1.get(line[17:20].strip(), "X")
            pts.append([resnum, float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[60:66])])
            seq.append(aa)
    if pts and max(p[4] for p in pts) <= 1.0:
        for p in pts:
            p[4] *= 100.0
    return pts, "".join(seq)


def _plddt_rgb(val: float) -> tuple[int, int, int]:
    if val >= 90.0:
        return (0, 83, 214)
    if val >= 70.0:
        return (101, 203, 243)
    if val >= 50.0:
        return (255, 219, 19)
    return (255, 125, 69)


def _shade(rgb: tuple[int, int, int], factor: float, fog: float) -> str:
    r = int(max(0, min(255, rgb[0] * factor * (1 - fog) + 255 * fog)))
    g = int(max(0, min(255, rgb[1] * factor * (1 - fog) + 255 * fog)))
    b = int(max(0, min(255, rgb[2] * factor * (1 - fog) + 255 * fog)))
    return f"#{r:02x}{g:02x}{b:02x}"


def _jacobi_pca(coords: list[tuple[float, float, float]]):
    n = len(coords)
    mx = sum(c[0] for c in coords) / n
    my = sum(c[1] for c in coords) / n
    mz = sum(c[2] for c in coords) / n
    centered = [(c[0] - mx, c[1] - my, c[2] - mz) for c in coords]
    cov = [[0.0] * 3 for _ in range(3)]
    for x, y, z in centered:
        v = (x, y, z)
        for i in range(3):
            for j in range(3):
                cov[i][j] += v[i] * v[j]
    V = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    A = [row[:] for row in cov]
    for _ in range(40):
        p, q, max_v = 0, 1, abs(A[0][1])
        for i, j in ((0, 2), (1, 2)):
            if abs(A[i][j]) > max_v:
                p, q, max_v = i, j, abs(A[i][j])
        if max_v < 1e-9:
            break
        theta = 0.5 * math.atan2(2 * A[p][q], A[q][q] - A[p][p])
        c, s = math.cos(theta), math.sin(theta)
        for i in range(3):
            A[i][p], A[i][q] = c * A[i][p] - s * A[i][q], s * A[i][p] + c * A[i][q]
        for i in range(3):
            A[p][i], A[q][i] = c * A[p][i] - s * A[q][i], s * A[p][i] + c * A[q][i]
        for i in range(3):
            V[i][p], V[i][q] = c * V[i][p] - s * V[i][q], s * V[i][p] + c * V[i][q]
    order = sorted(range(3), key=lambda k: -A[k][k])
    basis = [[V[row][k] for row in range(3)] for k in order]
    return [
        (
            sum(c[j] * basis[0][j] for j in range(3)),
            sum(c[j] * basis[1][j] for j in range(3)),
            sum(c[j] * basis[2][j] for j in range(3)),
        )
        for c in centered
    ]


def build_3d_gallery_svg(pdb_paths: list[Path]) -> str:
    valid_pdbs = [p for p in pdb_paths if p.is_file() and p.stat().st_size > 100][:3]
    W, H = 1440, 480
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<defs><style>',
        '  text { font-family: "Inter", "Helvetica Neue", Arial, sans-serif; fill: #0F172A; }',
        '  .ptitle { font-size: 14px; font-weight: 700; }',
        '  .psub { font-size: 11px; fill: #475569; }',
        '  .legend { font-size: 11.5px; font-weight: 600; fill: #1E293B; }',
        '</style></defs>',
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="#FFFFFF"/>',
    ]
    if not valid_pdbs:
        svg.append('<text x="28" y="60" class="ptitle">No showcase PDB files provided.</text></svg>')
        return "\n".join(svg)

    panel_w, gap, start_x, card_y, card_h = 450, 22, 22, 50, 410
    for idx, pdb_path in enumerate(valid_pdbs):
        px = start_x + idx * (panel_w + gap)
        ca, seq = _load_ca(pdb_path)
        mean_b = sum(c[4] for c in ca) / max(1, len(ca))
        svg.append(f'<rect x="{px}" y="{card_y}" width="{panel_w}" height="{card_h}" rx="12" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="1.2"/>')
        svg.append(f'<text x="{px+18}" y="{card_y+28}" class="ptitle">{chr(97+idx)}   {escape(pdb_path.stem)}</text>')
        svg.append(f'<text x="{px+18}" y="{card_y+46}" class="psub">{len(ca)} aa · mean pLDDT {mean_b:.1f}</text>')
        if len(ca) < 4:
            continue

        proj = _jacobi_pca([(c[1], c[2], c[3]) for c in ca])
        mc = re.search(r"(GDD|SDD|GDN)", seq)
        mc_lo, mc_hi = (ca[mc.start()][0], ca[mc.end() - 1][0]) if mc else (-1, -1)

        cx_mid, cy_mid = px + panel_w / 2.0, card_y + 224.0
        max_r = max(1e-6, max(max(abs(x), abs(y)) for x, y, _ in proj))
        scale = 165.0 / max_r
        z_min, z_max = min(z for _, _, z in proj), max(z for _, _, z in proj)
        z_span = max(1e-6, z_max - z_min)

        items = []
        for i in range(len(proj) - 1):
            x1, y1, z1 = proj[i]
            x2, y2, z2 = proj[i + 1]
            zm = 0.5 * (z1 + z2)
            dn = (zm - z_min) / z_span
            pA = (cx_mid + x1 * scale, cy_mid - y1 * scale)
            pB = (cx_mid + x2 * scale, cy_mid - y2 * scale)
            is_mc = mc_lo <= ca[i][0] <= mc_hi
            items.append((zm, pA, pB, _plddt_rgb(0.5 * (ca[i][4] + ca[i + 1][4])), dn, is_mc))

        items.sort(key=lambda it: it[0])
        for _, pA, pB, rgb, dn, is_mc in items:
            col = "#E11D48" if is_mc else _shade(rgb, 0.90 + 0.14 * dn, fog=(1.0 - dn) * 0.40)
            w = (6.0 if is_mc else 3.8) + 2.2 * dn
            svg.append(f'<line x1="{pA[0]:.1f}" y1="{pA[1]:.1f}" x2="{pB[0]:.1f}" y2="{pB[1]:.1f}" stroke="{col}" stroke-width="{w:.1f}" stroke-linecap="round"/>')

    svg.append('</svg>')
    return "\n".join(svg)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--validation-tsv", required=True, type=Path, help="results/esm_rdrp/esm_rdrp_validation.tsv")
    ap.add_argument("--manifest", default=None, type=Path, help="Optional manifest CSV for host names")
    ap.add_argument("--out-validation-svg", required=True, type=Path, help="Output 3-panel structural triage SVG")
    ap.add_argument("--pdb", action="append", default=[], type=Path, help="Optional showcase PDB file(s) for 3D gallery")
    ap.add_argument("--out-3d-svg", default=None, type=Path, help="Optional output 3D RdRp gallery SVG")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    hosts = load_manifest_hosts(args.manifest)
    val_svg = build_validation_svg(args.validation_tsv, hosts)
    atomic_write_text(args.out_validation_svg, val_svg)
    LOG.info("Wrote Fig. 2 structural validation SVG -> %s", args.out_validation_svg)

    if args.out_3d_svg is not None:
        gallery_svg = build_3d_gallery_svg(args.pdb)
        atomic_write_text(args.out_3d_svg, gallery_svg)
        LOG.info("Wrote Fig. 2 3D RdRp gallery SVG -> %s", args.out_3d_svg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
