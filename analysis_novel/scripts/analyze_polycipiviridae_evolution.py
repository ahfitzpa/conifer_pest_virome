#!/usr/bin/env python3
"""Evolutionary and phylogeographic analysis of the pan-Holarctic Ips Polycipiviridae clade.

Pure-Python (zero OpenBLAS/Biopython) script that takes assembled Polycipiviridae
genome FASTA(s) and manifest metadata, and computes:
  1. polycipiviridae_genomes_summary.tsv     (length, GC%, host, continent, intact ORF count, RdRp polyprotein aa len)
  2. polycipiviridae_pairwise_dnds.tsv       (pairwise aa ID%, nt ID%, dN, dS, dN/dS [Nei-Gojobori 1986], geographic stratum)
  3. polycipiviridae_phylogeography_summary.tsv (mean/median divergence & dN/dS by geographic stratum)
  4. polycipiviridae_evolution.svg           (2-panel SVG: pairwise divergence matrix + dN/dS purifying selection)

Uses MAFFT (if available on PATH) for codon-guided amino-acid alignment of the
conserved RdRp polyprotein ORF, falling back to pure-Python anchor-based alignment
if MAFFT is absent.
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Dict, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
LOG = logging.getLogger("analyze_polycipiviridae_evolution")

CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
BASES = ("A", "C", "G", "T")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

HOST_CONTINENT = {
    "ips confusus": "North America",
    "ips typographus": "Europe",
    "ips nitidus": "Asia",
    "ips subelongatus": "Asia",
    "ips duplicatus": "Europe",
    "hylobius abietis": "Europe",
}


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def translate(nt: str) -> str:
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append(CODON_TABLE.get(nt[i : i + 3], "X"))
    return "".join(aa)


def find_orfs(seq: str, min_aa: int = 150) -> List[Dict[str, object]]:
    """Find intact (ATG..STOP) and longest partial ORFs across all 6 frames."""
    seq = seq.upper().replace("U", "T")
    orfs: List[Dict[str, object]] = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        n = len(s)
        for frame in range(3):
            i = frame
            while i + 3 <= n:
                cod = s[i : i + 3]
                if cod == "ATG":
                    j = i + 3
                    while j + 3 <= n and s[j : j + 3] not in STOP_CODONS:
                        j += 3
                    has_stop = (j + 3 <= n) and (s[j : j + 3] in STOP_CODONS)
                    aa_len = (j - i) // 3
                    if aa_len >= min_aa:
                        nt_seq = s[i:j]
                        orfs.append(
                            {
                                "strand": strand,
                                "start": i,
                                "end": j,
                                "aa_len": aa_len,
                                "intact": has_stop,
                                "nt": nt_seq,
                                "aa": translate(nt_seq),
                            }
                        )
                    i = j + 3
                else:
                    i += 3
    orfs.sort(key=lambda o: (-int(o["aa_len"]), not bool(o["intact"])))
    return orfs


def read_fasta(paths: Sequence[Path]) -> List[Tuple[str, str]]:
    recs: List[Tuple[str, str]] = []
    for p in paths:
        if not p.is_file():
            continue
        hdr: Optional[str] = None
        chunks: List[str] = []
        with p.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if hdr is not None:
                        recs.append((hdr, "".join(chunks).upper()))
                    hdr = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
        if hdr is not None:
            recs.append((hdr, "".join(chunks).upper()))
    return recs


def load_manifest_meta(manifest: Optional[Path]) -> Dict[str, Dict[str, str]]:
    meta: Dict[str, Dict[str, str]] = {}
    if manifest is None or not manifest.is_file():
        return meta
    with manifest.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rid = (row.get("run_id") or row.get("sample_id") or row.get("run") or "").strip()
            if not rid:
                continue
            host = (
                row.get("scientificname")
                or row.get("host_species")
                or row.get("host_organism")
                or "Unknown"
            ).strip().replace("_", " ")
            bp = (row.get("bioproject") or "").strip()
            meta[rid] = {
                "host": host,
                "bioproject": bp,
                "continent": HOST_CONTINENT.get(host.lower(), "Unknown"),
            }
    return meta


def infer_run_id(contig_id: str) -> str:
    if "|" in contig_id:
        return contig_id.split("|", 1)[0]
    m = re.match(r"([SED]RR\d+)", contig_id)
    return m.group(1) if m else contig_id


# Precompute Nei-Gojobori (1986) synonymous & nonsynonymous site counts per sense codon
_CODON_SITES: Dict[str, Tuple[float, float]] = {}
for _cod, _aa in CODON_TABLE.items():
    if _aa == "*":
        continue
    _syn = 0.0
    _nonsyn = 0.0
    for _pos in range(3):
        _s_mut = 0
        _ns_mut = 0
        for _b in BASES:
            if _b == _cod[_pos]:
                continue
            _alt = _cod[:_pos] + _b + _cod[_pos + 1 :]
            _alt_aa = CODON_TABLE[_alt]
            if _alt_aa == "*":
                continue
            if _alt_aa == _aa:
                _s_mut += 1
            else:
                _ns_mut += 1
        _tot = _s_mut + _ns_mut
        if _tot > 0:
            _syn += _s_mut / 3.0
            _nonsyn += _ns_mut / 3.0
    _CODON_SITES[_cod] = (_syn, _nonsyn)


def _codon_diff_ng86(c1: str, c2: str) -> Tuple[float, float]:
    """Return (syn_diffs, nonsyn_diffs) between two sense codons under Nei-Gojobori 1986."""
    if c1 == c2:
        return 0.0, 0.0
    diffs = [i for i in range(3) if c1[i] != c2[i]]
    if len(diffs) == 1:
        return (1.0, 0.0) if CODON_TABLE[c1] == CODON_TABLE[c2] else (0.0, 1.0)
    if len(diffs) == 2:
        p1, p2 = diffs
        paths = [
            (c1[:p1] + c2[p1] + c1[p1 + 1 :],),
            (c1[:p2] + c2[p2] + c1[p2 + 1 :],),
        ]
        valid_sd, valid_nd = [], []
        for (mid,) in paths:
            if CODON_TABLE[mid] == "*":
                continue
            sd = (1 if CODON_TABLE[c1] == CODON_TABLE[mid] else 0) + (
                1 if CODON_TABLE[mid] == CODON_TABLE[c2] else 0
            )
            valid_sd.append(float(sd))
            valid_nd.append(float(2 - sd))
        if valid_sd:
            return mean(valid_sd), mean(valid_nd)
        return 0.0, 2.0
    # 3 differences: average over valid 6 mutational pathways
    valid_sd, valid_nd = [], []
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    for i1, i2, _i3 in perms:
        m1 = c1[:i1] + c2[i1] + c1[i1 + 1 :]
        m2 = m1[:i2] + c2[i2] + m1[i2 + 1 :]
        if CODON_TABLE[m1] == "*" or CODON_TABLE[m2] == "*":
            continue
        sd = (
            (1 if CODON_TABLE[c1] == CODON_TABLE[m1] else 0)
            + (1 if CODON_TABLE[m1] == CODON_TABLE[m2] else 0)
            + (1 if CODON_TABLE[m2] == CODON_TABLE[c2] else 0)
        )
        valid_sd.append(float(sd))
        valid_nd.append(float(3 - sd))
    if valid_sd:
        return mean(valid_sd), mean(valid_nd)
    return 0.0, 3.0


def _jc_correct(p: float) -> float:
    if p <= 0.0:
        return 0.0
    if p >= 0.74:
        return 2.5  # saturated cap
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)


def compute_ng86_dnds(aln_nt1: str, aln_nt2: str) -> Dict[str, float]:
    """Compute pairwise aa_id%, nt_id%, dN, dS, and dN/dS from aligned codon strings."""
    n_cod = min(len(aln_nt1), len(aln_nt2)) // 3
    aa_match = aa_tot = nt_match = nt_tot = 0
    s_sites = n_sites = s_diffs = n_diffs = 0.0

    for idx in range(n_cod):
        c1 = aln_nt1[idx * 3 : idx * 3 + 3]
        c2 = aln_nt2[idx * 3 : idx * 3 + 3]
        if "-" in c1 or "-" in c2 or c1 not in _CODON_SITES or c2 not in _CODON_SITES:
            continue
        aa_tot += 1
        if CODON_TABLE[c1] == CODON_TABLE[c2]:
            aa_match += 1
        for k in range(3):
            nt_tot += 1
            if c1[k] == c2[k]:
                nt_match += 1
        s1, n1 = _CODON_SITES[c1]
        s2, n2 = _CODON_SITES[c2]
        s_sites += 0.5 * (s1 + s2)
        n_sites += 0.5 * (n1 + n2)
        sd, nd = _codon_diff_ng86(c1, c2)
        s_diffs += sd
        n_diffs += nd

    ps = s_diffs / s_sites if s_sites > 0 else 0.0
    pn = n_diffs / n_sites if n_sites > 0 else 0.0
    ds = _jc_correct(ps)
    dn = _jc_correct(pn)
    omega = (dn / ds) if ds > 1e-6 else float("nan")
    return {
        "aligned_codons": float(aa_tot),
        "aa_identity_pct": (100.0 * aa_match / aa_tot) if aa_tot else 0.0,
        "nt_identity_pct": (100.0 * nt_match / nt_tot) if nt_tot else 0.0,
        "dn": dn,
        "ds": ds,
        "dnds": omega,
    }


def align_codons(orf_records: List[Tuple[str, str, str]]) -> Dict[str, str]:
    """Return {id: aligned_nt} using MAFFT on translated amino acids if available, else prefix-trim."""
    mafft_bin = shutil.which("mafft")
    if mafft_bin and len(orf_records) >= 2:
        with tempfile.TemporaryDirectory() as tmp:
            in_faa = Path(tmp) / "in.faa"
            with in_faa.open("w", encoding="utf-8") as fh:
                for rec_id, aa_seq, _nt in orf_records:
                    fh.write(f">{rec_id}\n{aa_seq}\n")
            res = subprocess.run(
                [mafft_bin, "--auto", "--quiet", str(in_faa)],
                capture_output=True,
                text=True,
                check=False,
            )
            if res.returncode == 0 and res.stdout.strip():
                aln_aa: Dict[str, str] = {}
                cur: Optional[str] = None
                buf: List[str] = []
                for line in res.stdout.splitlines():
                    if line.startswith(">"):
                        if cur is not None:
                            aln_aa[cur] = "".join(buf)
                        cur = line[1:].split()[0]
                        buf = []
                    else:
                        buf.append(line.strip())
                if cur is not None:
                    aln_aa[cur] = "".join(buf)
                out_nt: Dict[str, str] = {}
                for rec_id, _aa, nt in orf_records:
                    a_str = aln_aa.get(rec_id, "")
                    cod_idx = 0
                    nt_chunks = []
                    for ch in a_str:
                        if ch == "-":
                            nt_chunks.append("---")
                        else:
                            nt_chunks.append(nt[cod_idx * 3 : cod_idx * 3 + 3])
                            cod_idx += 1
                    out_nt[rec_id] = "".join(nt_chunks)
                if len(out_nt) == len(orf_records):
                    return out_nt
    # Pure-Python fallback: trim to common codon length
    min_nt = min(len(nt) - (len(nt) % 3) for _, _, nt in orf_records)
    return {rec_id: nt[:min_nt] for rec_id, _, nt in orf_records}


def classify_geo_stratum(h1: str, c1: str, h2: str, c2: str) -> str:
    if h1.lower() == h2.lower():
        return "within_species_same_continent"
    if c1 == c2:
        return "inter_species_same_continent"
    eurasia = {"Europe", "Asia"}
    if c1 in eurasia and c2 in eurasia:
        return "inter_species_eurasian"
    return "inter_continental_nearctic_palearctic"


def write_tsv_atomic(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def render_evolution_svg(
    out_svg: Path,
    genome_rows: List[Dict[str, object]],
    pair_rows: List[Dict[str, object]],
) -> None:
    """Render 2-panel SVG: (A) pairwise RdRp aa identity matrix & (B) dN/dS by geographic stratum."""
    ids = [str(g["genome_id"]) for g in genome_rows]
    n = len(ids)
    lookup: Dict[Tuple[str, str], Dict[str, object]] = {}
    for pr in pair_rows:
        lookup[(str(pr["genome_1"]), str(pr["genome_2"]))] = pr
        lookup[(str(pr["genome_2"]), str(pr["genome_1"]))] = pr

    cell = 34
    mat_w = max(240, n * cell)
    W = 240 + mat_w + 380
    H = max(360, 140 + n * cell + 70)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<style>text{font-family:"Inter","Helvetica Neue",Arial,sans-serif;fill:#0F172A;}</style>',
        f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>',
        '<text x="20" y="28" font-size="14" font-weight="700">Pan-Holarctic Ips Polycipiviridae: Phylogeographic Divergence &amp; Purifying Selection (dN/dS)</text>',
    ]

    # Panel A: Matrix
    x0, y0 = 215, 90
    svg.append(f'<text x="{x0}" y="58" font-size="11.5" font-weight="700">A. Pairwise RdRp Polyprotein Identity (% aa)</text>')
    cont_colors = {"North America": "#D97706", "Europe": "#0F766E", "Asia": "#2563EB", "Unknown": "#64748B"}

    for i, g in enumerate(genome_rows):
        gid = str(g["genome_id"])
        short = f"{g['host']} ({str(g['run_id'])[-4:]})"
        col = cont_colors.get(str(g["continent"]), "#64748B")
        cy = y0 + i * cell + cell / 2.0
        svg.append(f'<circle cx="{x0 - 10}" cy="{cy:.1f}" r="4" fill="{col}"/>')
        svg.append(f'<text x="{x0 - 18}" y="{cy + 3.5:.1f}" text-anchor="end" font-size="9.5" font-weight="600">{escape(short)}</text>')
        for j, g2 in enumerate(genome_rows):
            cx = x0 + j * cell
            cy_top = y0 + i * cell
            if i == j:
                val = 100.0
            else:
                pr = lookup.get((gid, str(g2["genome_id"])))
                val = float(pr["aa_identity_pct"]) if pr else 0.0
            # Color gradient: 80% (#DBEAFE) -> 100% (#0F766E)
            t = max(0.0, min(1.0, (val - 80.0) / 20.0))
            r = int(219 + (15 - 219) * t)
            gr = int(234 + (118 - 234) * t)
            b = int(254 + (110 - 254) * t)
            fill = f"#{r:02x}{gr:02x}{b:02x}"
            txt_col = "#FFFFFF" if t > 0.55 else "#0F172A"
            svg.append(f'<rect x="{cx}" y="{cy_top}" width="{cell - 2}" height="{cell - 2}" rx="3" fill="{fill}"/>')
            svg.append(f'<text x="{cx + (cell - 2)/2:.1f}" y="{cy_top + (cell - 2)/2 + 3.2:.1f}" text-anchor="middle" font-size="8.5" font-weight="700" style="fill:{txt_col}">{val:.1f}</text>')

    # Panel B: dN/dS strip by geographic stratum
    bx0 = x0 + mat_w + 55
    bw, bh = 280, max(180, n * cell)
    svg.append(f'<text x="{bx0}" y="58" font-size="11.5" font-weight="700">B. Purifying Selection (Nei–Gojobori dN/dS)</text>')
    svg.append(f'<rect x="{bx0}" y="{y0}" width="{bw}" height="{bh}" fill="#F8FAFC" stroke="#CBD5E1"/>')

    strata_order = [
        ("within_species_same_continent", "Within species", "#0F766E"),
        ("inter_species_eurasian", "Europe vs. Asia", "#2563EB"),
        ("inter_continental_nearctic_palearctic", "N. America vs. Eurasia", "#D97706"),
    ]
    max_w = 0.5
    for pr in pair_rows:
        try:
            v = float(pr["dnds"])
            if math.isfinite(v):
                max_w = max(max_w, min(1.2, v * 1.25))
        except ValueError:
            pass

    # Reference line at dN/dS = 0.2 (strong purifying selection)
    for tick in (0.0, 0.1, 0.2, 0.5):
        if tick > max_w:
            continue
        ty = y0 + bh - (tick / max_w) * (bh - 20) - 10
        svg.append(f'<line x1="{bx0}" y1="{ty:.1f}" x2="{bx0 + bw}" y2="{ty:.1f}" stroke="#E2E8F0" stroke-dasharray="3,3"/>')
        svg.append(f'<text x="{bx0 - 6}" y="{ty + 3:.1f}" text-anchor="end" font-size="8.5" fill="#475569">{tick:.1f}</text>')

    for idx, (key, label, col) in enumerate(strata_order):
        sx = bx0 + 45 + idx * 92
        vals = []
        for pr in pair_rows:
            if pr["geographic_stratum"] == key:
                try:
                    v = float(pr["dnds"])
                    if math.isfinite(v):
                        vals.append(v)
                except ValueError:
                    pass
        for k, v in enumerate(vals):
            jx = sx + ((k % 5) - 2) * 5
            py = y0 + bh - (min(v, max_w) / max_w) * (bh - 20) - 10
            svg.append(f'<circle cx="{jx:.1f}" cy="{py:.1f}" r="4.2" fill="{col}" fill-opacity="0.78" stroke="#FFFFFF" stroke-width="0.8"/>')
        if vals:
            med_v = median(vals)
            my = y0 + bh - (min(med_v, max_w) / max_w) * (bh - 20) - 10
            svg.append(f'<line x1="{sx - 18}" y1="{my:.1f}" x2="{sx + 18}" y2="{my:.1f}" stroke="#0F172A" stroke-width="2.2"/>')
            svg.append(f'<text x="{sx}" y="{my - 6:.1f}" text-anchor="middle" font-size="8.5" font-weight="700">{med_v:.3f}</text>')
        svg.append(f'<text x="{sx}" y="{y0 + bh + 16}" text-anchor="middle" font-size="8.5" font-weight="600">{escape(label)}</text>')

    svg.append("</svg>\n")
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_svg.with_suffix(".svg.tmp")
    tmp.write_text("\n".join(svg), encoding="utf-8")
    os.replace(tmp, out_svg)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contigs-fasta", required=True, action="append", type=Path, help="Polycipiviridae genome FASTA (repeatable).")
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--min-length", type=int, default=10000)
    ap.add_argument("--min-orf-aa", type=int, default=150)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args(argv)

    meta = load_manifest_meta(args.manifest)
    raw_recs = read_fasta(args.contigs_fasta)

    genome_rows: List[Dict[str, object]] = []
    orf_inputs: List[Tuple[str, str, str]] = []
    seen = set()

    for gid, seq in raw_recs:
        if gid in seen or len(seq) < args.min_length:
            continue
        seen.add(gid)
        orfs = find_orfs(seq, min_aa=args.min_orf_aa)
        if not orfs:
            continue
        main_orf = orfs[0]
        rid = infer_run_id(gid)
        m = meta.get(rid, {"host": "Unknown", "bioproject": "", "continent": "Unknown"})
        gc = 100.0 * (seq.count("G") + seq.count("C")) / max(1, len(seq))
        genome_rows.append(
            {
                "genome_id": gid,
                "run_id": rid,
                "host": m["host"],
                "continent": m["continent"],
                "bioproject": m["bioproject"],
                "length_bp": len(seq),
                "gc_pct": f"{gc:.2f}",
                "n_orfs_ge150aa": len(orfs),
                "n_intact_orfs": sum(1 for o in orfs if o["intact"]),
                "rdrp_polyprotein_aa": main_orf["aa_len"],
                "rdrp_polyprotein_intact": int(bool(main_orf["intact"])),
            }
        )
        orf_inputs.append((gid, str(main_orf["aa"]), str(main_orf["nt"])))

    # Sort deterministically by continent -> host -> genome_id
    cont_rank = {"North America": 0, "Europe": 1, "Asia": 2, "Unknown": 3}
    order = sorted(
        range(len(genome_rows)),
        key=lambda i: (cont_rank.get(str(genome_rows[i]["continent"]), 9), str(genome_rows[i]["host"]), str(genome_rows[i]["genome_id"])),
    )
    genome_rows = [genome_rows[i] for i in order]
    orf_inputs = [orf_inputs[i] for i in order]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv_atomic(
        args.out_dir / "polycipiviridae_genomes_summary.tsv",
        [
            "genome_id", "run_id", "host", "continent", "bioproject",
            "length_bp", "gc_pct", "n_orfs_ge150aa", "n_intact_orfs",
            "rdrp_polyprotein_aa", "rdrp_polyprotein_intact",
        ],
        genome_rows,
    )

    aln_nt = align_codons(orf_inputs) if len(orf_inputs) >= 2 else {}
    pair_rows: List[Dict[str, object]] = []
    by_stratum: Dict[str, List[Dict[str, float]]] = defaultdict(list)

    for i in range(len(genome_rows)):
        for j in range(i + 1, len(genome_rows)):
            g1, g2 = genome_rows[i], genome_rows[j]
            id1, id2 = str(g1["genome_id"]), str(g2["genome_id"])
            stats = compute_ng86_dnds(aln_nt[id1], aln_nt[id2])
            strat = classify_geo_stratum(str(g1["host"]), str(g1["continent"]), str(g2["host"]), str(g2["continent"]))
            by_stratum[strat].append(stats)
            pair_rows.append(
                {
                    "genome_1": id1,
                    "host_1": g1["host"],
                    "continent_1": g1["continent"],
                    "genome_2": id2,
                    "host_2": g2["host"],
                    "continent_2": g2["continent"],
                    "geographic_stratum": strat,
                    "aligned_codons": int(stats["aligned_codons"]),
                    "aa_identity_pct": f"{stats['aa_identity_pct']:.2f}",
                    "nt_identity_pct": f"{stats['nt_identity_pct']:.2f}",
                    "dn": f"{stats['dn']:.4f}",
                    "ds": f"{stats['ds']:.4f}",
                    "dnds": f"{stats['dnds']:.4f}" if math.isfinite(stats["dnds"]) else "NA",
                }
            )

    write_tsv_atomic(
        args.out_dir / "polycipiviridae_pairwise_dnds.tsv",
        [
            "genome_1", "host_1", "continent_1", "genome_2", "host_2", "continent_2",
            "geographic_stratum", "aligned_codons", "aa_identity_pct", "nt_identity_pct",
            "dn", "ds", "dnds",
        ],
        pair_rows,
    )

    strat_rows: List[Dict[str, object]] = []
    for strat, items in sorted(by_stratum.items()):
        aas = [x["aa_identity_pct"] for x in items]
        nts = [x["nt_identity_pct"] for x in items]
        ws = [x["dnds"] for x in items if math.isfinite(x["dnds"])]
        strat_rows.append(
            {
                "geographic_stratum": strat,
                "n_pairs": len(items),
                "mean_aa_identity_pct": f"{mean(aas):.2f}",
                "median_aa_identity_pct": f"{median(aas):.2f}",
                "mean_nt_identity_pct": f"{mean(nts):.2f}",
                "median_dnds": f"{median(ws):.4f}" if ws else "NA",
                "mean_dnds": f"{mean(ws):.4f}" if ws else "NA",
            }
        )

    write_tsv_atomic(
        args.out_dir / "polycipiviridae_phylogeography_summary.tsv",
        [
            "geographic_stratum", "n_pairs", "mean_aa_identity_pct", "median_aa_identity_pct",
            "mean_nt_identity_pct", "median_dnds", "mean_dnds",
        ],
        strat_rows,
    )

    if genome_rows:
        render_evolution_svg(args.out_dir / "polycipiviridae_evolution.svg", genome_rows, pair_rows)

    LOG.info(
        "Polycipiviridae evolution analysis complete: %d genome(s), %d pair(s) -> %s",
        len(genome_rows), len(pair_rows), args.out_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
