#!/usr/bin/env python3
"""Render a linear (rectangular) publication-grade ML phylogram SVG from a Newick tree
and selected.tsv metadata, with colour-coded Host/Source, CheckV Quality, and CheckV
Completeness heatmap columns for both ICTV references and Serratus candidates.
"""
from __future__ import annotations

import argparse
import csv
import html
import re
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

INK = "#101214"
MUTED = "#505258"
REF_GREY = "#94A3B8"
CAND_HIGH = "#1868DB"
CAND_MED = "#C75300"
CAND_LOW = "#2898BD"
SUPP_FILL = "#101214"

# BioProject -> (Category, HexColour, DisplayLabel) for Serratus candidates
BIOPROJECT_SOURCE: dict[str, tuple[str, str, str]] = {
    # Vertebrate / Zoonotic
    "PRJEB6505": ("Vertebrate / Zoonotic", "#1868DB", "Vert / Zoon."),
    "PRJEB26687": ("Vertebrate / Zoonotic", "#1868DB", "Vert / Zoon."),
    "PRJEB17935": ("Vertebrate / Zoonotic", "#1868DB", "Vert / Zoon."),
    "PRJNA751908": ("Vertebrate / Zoonotic", "#1868DB", "Vert / Zoon."),
    "PRJDB7591": ("Vertebrate / Zoonotic", "#1868DB", "Vert / Zoon."),
    "PRJNA683917": ("Vertebrate / Zoonotic", "#1868DB", "Vert / Zoon."),
    # Crustacean (Lobster / Crayfish)
    "PRJNA655509": ("Crustacean (Lobster/Crayfish)", "#C75300", "Crustacean"),
    "PRJNA702424": ("Crustacean (Lobster/Crayfish)", "#C75300", "Crustacean"),
    "PRJNA669582": ("Crustacean (Lobster/Crayfish)", "#C75300", "Crustacean"),
    "PRJNA754774": ("Crustacean (Lobster/Crayfish)", "#C75300", "Crustacean"),
    # Marine / Coastal
    "PRJEB34490": ("Marine / Coastal", "#2898BD", "Marine"),
    "PRJNA576751": ("Marine / Coastal", "#2898BD", "Marine"),
    "PRJNA555787": ("Marine / Coastal", "#2898BD", "Marine"),
    "PRJNA665326": ("Marine / Coastal", "#2898BD", "Marine"),
    "PRJNA654739": ("Marine / Coastal", "#2898BD", "Marine"),
    # Freshwater
    "PRJNA640812": ("Freshwater", "#16794F", "Freshwater"),
    "PRJNA691101": ("Freshwater", "#16794F", "Freshwater"),
    "PRJNA657047": ("Freshwater", "#16794F", "Freshwater"),
    "PRJNA663589": ("Freshwater", "#16794F", "Freshwater"),
    # Soil / Wastewater
    "PRJNA716119": ("Soil / Wastewater", "#8F4FDB", "Soil / WWTP"),
    "PRJNA616041": ("Soil / Wastewater", "#8F4FDB", "Soil / WWTP"),
    "PRJNA571457": ("Soil / Wastewater", "#8F4FDB", "Soil / WWTP"),
    "PRJNA673783": ("Soil / Wastewater", "#8F4FDB", "Soil / WWTP"),
}


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    x: float = 0.0
    y: float = 0.0


def parse_newick(s: str) -> Node:
    s = s.strip().rstrip(";")
    root = Node()
    cur = root
    i, n = 0, len(s)
    while i < n:
        ch = s[i]
        if ch == "(":
            c = Node(parent=cur)
            cur.children.append(c)
            cur = c
            i += 1
        elif ch == ",":
            assert cur.parent is not None
            sib = Node(parent=cur.parent)
            cur.parent.children.append(sib)
            cur = sib
            i += 1
        elif ch == ")":
            assert cur.parent is not None
            cur = cur.parent
            i += 1
            j = i
            while j < n and s[j] not in ",():;":
                j += 1
            tok = s[i:j].strip()
            if tok:
                cur.name = tok
            i = j
        elif ch == ":":
            i += 1
            j = i
            while j < n and s[j] not in ",();":
                j += 1
            try:
                cur.length = max(0.0, float(s[i:j]))
            except ValueError:
                cur.length = 0.0
            i = j
        else:
            j = i
            while j < n and s[j] not in ",():;":
                j += 1
            cur.name = s[i:j].strip()
            i = j
    return root


def subtree_size(node: Node) -> int:
    if not node.children:
        return 1
    return sum(subtree_size(c) for c in node.children)


def ladderize(node: Node, reverse: bool = False) -> None:
    """Sort children by subtree size so clades cascade cleanly."""
    if not node.children:
        return
    for c in node.children:
        ladderize(c, reverse=reverse)
    node.children.sort(key=subtree_size, reverse=reverse)


def get_leaves(node: Node) -> list[Node]:
    if not node.children:
        return [node]
    out: list[Node] = []
    for c in node.children:
        out.extend(get_leaves(c))
    return out


def assign_xy(node: Node, leaf_y: dict[int, float], depth: float = 0.0) -> tuple[float, float]:
    node.x = depth
    if not node.children:
        y = leaf_y[id(node)]
        node.y = y
        return y, y
    ys: list[float] = []
    for c in node.children:
        y0, y1 = assign_xy(c, leaf_y, depth + c.length)
        ys.extend([y0, y1])
    node.y = (min(ys) + max(ys)) / 2.0
    return min(ys), max(ys)


def parse_support(label: str) -> tuple[float | None, float | None]:
    if "/" in label:
        parts = label.split("/", 1)
        try:
            return float(parts[0]), float(parts[1])
        except ValueError:
            return None, None
    return None, None


def extract_accession(label: str) -> str:
    parts = [p for p in re.split(r"[_|]", label) if p]
    if parts and parts[0].upper() == "ICTV":
        parts = parts[1:]
    if len(parts) >= 2 and parts[0] in {"NC", "AC", "NG", "NM", "NP", "NR", "XM", "XP"}:
        return f"{parts[0]}_{parts[1]}"
    return parts[0] if parts else label


def fetch_ncbi_hosts(accessions: list[str]) -> dict[str, tuple[str, str]]:
    """Query NCBI Entrez nuccore GenBank records for /host and ORGANISM."""
    out: dict[str, tuple[str, str]] = {}
    batch_size = 100
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nuccore&id={','.join(batch)}&rettype=gb&retmode=text"
        )
        text = urllib.request.urlopen(url, timeout=30).read().decode("utf-8", errors="replace")
        for rec in text.split("\n//\n"):
            m_acc = re.search(r"^ACCESSION\s+(\S+)", rec, re.M)
            if not m_acc:
                continue
            acc = m_acc.group(1)
            m_host = re.search(r'/host="([^"]+)"', rec)
            m_org = re.search(r"^\s+ORGANISM\s+(.+)", rec, re.M)
            out[acc] = (
                m_host.group(1).strip() if m_host else "",
                m_org.group(1).strip() if m_org else "",
            )
    return out


def classify_ictv_host(host_str: str, org_str: str, label: str) -> tuple[str, str, str]:
    """Map an ICTV reference's host/organism/label to (Category, HexColour, DisplayLabel)."""
    s = f"{host_str} {org_str} {label}".lower()
    # Vertebrate hosts
    if any(k in s for k in ("homo sapiens", "human")):
        return ("Human", "#1C2B42", "Human")
    if any(k in s for k in ("sus scrofa", "swine", "pig", "camelus", "camel", "alces", "moose", "rabbit", "oryctolagus")):
        return ("Swine / Ungulate", "#B34B00", "Swine / Ung.")
    if any(k in s for k in ("mustela", "ferret")):
        return ("Mustelid (Ferret)", "#9E3B6B", "Ferret")
    if any(k in s for k in ("rattus", "apodemus", "eothenomys", "microtus", "cricetulus", "crocidura", "ochotona", "rat", "vole", "shrew", "rodent")):
        return ("Rodent / Shrew", "#6B4FBB", "Rodent")
    if any(k in s for k in ("eptesicus", "rhinolophus", "desmodus", "bat", "chirohepevirus")):
        return ("Bat", "#42526E", "Bat")
    if any(k in s for k in ("chicken", "egretta", "egret", "avian", "avihepevirus", "gallus")):
        return ("Bird", "#0E7490", "Bird")
    if any(k in s for k in ("trout", "oncorhynchus", "piscihepevirus", "fish")):
        return ("Fish", "#15803D", "Fish")
    # Insect / Arthropod hosts (Iflaviridae, Dicistroviridae, Polycipiviridae, etc.)
    if any(k in s for k in ("apis", "bombus", "vespa", "solenopsis", "formica", "bee", "wasp", "ant", "hymenoptera", "aladeformis", "sacbrood", "varroa")):
        return ("Hymenoptera (Bee/Ant)", "#B34B00", "Bee / Hym.")
    if any(k in s for k in ("bombyx", "perina", "ectropis", "spodoptera", "helicoverpa", "lymantria", "moth", "silkworm", "lepidoptera", "flacherie", "betanuda")):
        return ("Lepidoptera (Moth)", "#6B4FBB", "Moth / Lep.")
    if any(k in s for k in ("nilaparvata", "brevicoryne", "halyomorpha", "aphid", "planthopper", "bug", "hemiptera", "lygus")):
        return ("Hemiptera (Bug/Aphid)", "#0E7490", "Bug / Hem.")
    if any(k in s for k in ("acheta", "locusta", "gryllus", "cricket", "locust", "orthoptera")):
        return ("Orthoptera (Cricket)", "#9E3B6B", "Cricket / Ort.")
    if any(k in s for k in ("coleoptera", "beetle", "weevil", "dinocampus", "tribolium", "ips", "hylobius")):
        return ("Coleoptera (Beetle)", "#16794F", "Beetle / Col.")
    if any(k in s for k in ("drosophila", "aedes", "culex", "anopheles", "musca", "fly", "mosquito", "diptera")):
        return ("Diptera (Fly/Mosquito)", "#42526E", "Fly / Dip.")
    if any(k in s for k in ("iflavirus", "dicistrovirus", "insect", "arthropod")):
        return ("Other Insect", "#64748B", "Insect")
    return ("Other Reference", "#64748B", "Other Ref")


def short_label(name: str, is_ref: bool) -> str:
    if is_ref and ("_" in name or "|" in name):
        parts = [p for p in re.split(r"[_|]", name) if p]
        if parts and parts[0].upper() == "ICTV":
            parts = parts[1:]
        if parts and parts[-1].upper() == "ICTV":
            parts = parts[:-1]
        if len(parts) >= 3 and parts[0] in {"NC", "AC"}:
            acc = f"{parts[0]}_{parts[1]}"
            sp = " ".join(parts[2:])
        elif len(parts) >= 2:
            acc = parts[0]
            sp = " ".join(parts[1:])
        else:
            return name
        for g in ("Paslahepevirus ", "Rocahepevirus ", "Avihepevirus ", "Chirohepevirus ", "Piscihepevirus "):
            if sp.startswith(g):
                sp = g[0] + ". " + sp[len(g) :]
        return f"{sp} | {acc}" if len(sp) <= 25 else f"{sp[:23]}… | {acc}"
    return name


def comp_colour(pct: float | None) -> str:
    if pct is None:
        return "#E2E8F0"
    if pct >= 90:
        return "#1C2B42"
    if pct >= 70:
        return "#1868DB"
    if pct >= 50:
        return "#2898BD"
    return "#94A3B8"


def format_title_with_italic_family(title: str) -> str:
    """Italicize viral family names ending in -viridae inside SVG <text>."""
    tokens = re.split(r"(\b[A-Z][a-z]+viridae\b)", title)
    out: list[str] = []
    for tok in tokens:
        if re.fullmatch(r"[A-Z][a-z]+viridae", tok):
            out.append(f'<tspan font-style="italic">{html.escape(tok)}</tspan>')
        else:
            out.append(html.escape(tok))
    return "".join(out)


def render_tree_svg(
    tree_path: Path,
    meta_path: Path,
    out_svg: Path,
    title: str = "Hepeviridae Maximum-Likelihood Phylogeny",
    fetch_ncbi: bool = False,
) -> None:
    tree = parse_newick(tree_path.read_text(encoding="utf-8"))
    ladderize(tree, reverse=False)
    meta: dict[str, dict[str, str]] = {}
    with meta_path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["label"]] = r

    leaves = get_leaves(tree)
    n_leaves = len(leaves)

    ref_hosts: dict[str, tuple[str, str]] = {}
    if fetch_ncbi:
        ref_accs = [
            extract_accession(lf.name)
            for lf in leaves
            if meta.get(lf.name, {}).get("origin") == "reference"
        ]
        if ref_accs:
            ref_hosts = fetch_ncbi_hosts(ref_accs)

    y_top = 26.0
    dy = 10.2
    leaf_y = {id(lf): y_top + idx * dy for idx, lf in enumerate(leaves)}
    assign_xy(tree, leaf_y, 0.0)
    max_x = max(lf.x for lf in leaves) or 1.0

    w = 840
    legend_h = 116
    h = int(y_top + n_leaves * dy + legend_h)
    x0 = 24.0
    tree_w = 390.0
    x_scale = tree_w / max_x
    label_x = x0 + tree_w + 12.0
    tile_host = w - 156.0   # width 68
    tile_qual = w - 82.0    # width 30
    tile_comp = w - 46.0    # width 32

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Inter,sans-serif">',
        f'<rect width="{w}" height="{h}" fill="#FFFFFF"/>',
        f'<text x="{tile_host + 34}" y="14" text-anchor="middle" font-size="8.5" font-weight="700" fill="{INK}">Host / Source</text>',
        f'<text x="{tile_qual + 15}" y="14" text-anchor="middle" font-size="8.5" font-weight="700" fill="{INK}">Origin</text>',
        f'<text x="{tile_comp + 16}" y="14" text-anchor="middle" font-size="8.5" font-weight="700" fill="{INK}">CheckV</text>',
        f'<line x1="16" y1="18" x2="{w - 14}" y2="18" stroke="#E2E8F0" stroke-width="1"/>',
    ]

    # Subtle zebra striping across label + heatmap columns for horizontal readability
    for idx, lf in enumerate(leaves):
        if idx % 2 == 1:
            parts.append(f'<rect x="{label_x - 4:.1f}" y="{lf.y - dy / 2:.2f}" width="{w - label_x - 10:.1f}" height="{dy:.2f}" fill="#F8FAFC"/>')

    def px(val: float) -> float:
        return x0 + val * x_scale

    def draw_node(node: Node) -> None:
        if not node.children:
            return
        nx = px(node.x)
        y_min = min(c.y for c in node.children)
        y_max = max(c.y for c in node.children)
        parts.append(f'<line x1="{nx:.2f}" y1="{y_min:.2f}" x2="{nx:.2f}" y2="{y_max:.2f}" stroke="{INK}" stroke-width="1.05" stroke-linecap="round"/>')
        for c in node.children:
            cx = px(c.x)
            parts.append(f'<line x1="{nx:.2f}" y1="{c.y:.2f}" x2="{cx:.2f}" y2="{c.y:.2f}" stroke="{INK}" stroke-width="1.05" stroke-linecap="round"/>')
            draw_node(c)
        alrt, ufb = parse_support(node.name)
        if alrt is not None and ufb is not None and alrt >= 80 and ufb >= 95:
            parts.append(f'<circle cx="{nx:.2f}" cy="{node.y:.2f}" r="2.2" fill="{SUPP_FILL}"/>')

    draw_node(tree)

    cand_counts: dict[tuple[str, str], int] = {}
    ref_counts: dict[tuple[str, str], int] = {}
    n_refs = 0

    for lf in leaves:
        m = meta.get(lf.name, {})
        is_ref = m.get("origin") == "reference"
        q = m.get("checkv_quality", "")
        bp = m.get("bioproject", "")
        if is_ref:
            n_refs += 1
            col = REF_GREY
            tcol = MUTED
            fw = "400"
            fs = "italic"
            qual_txt = "Ref"
            qual_tcol = "#FFFFFF"
            comp_str = "Ref"
            comp_val = None
            acc = extract_accession(lf.name)
            h_str, o_str = ref_hosts.get(acc, (m.get("host", ""), ""))
            h_cat, h_col, h_txt = classify_ictv_host(h_str, o_str, lf.name)
            ref_counts[(h_cat, h_col)] = ref_counts.get((h_cat, h_col), 0) + 1
            h_tip = f"{h_cat}: {h_str or o_str or acc}"
        else:
            col = CAND_HIGH if "High" in q else (CAND_MED if "Medium" in q else CAND_LOW)
            tcol = INK
            fw = "700"
            fs = "normal"
            qual_txt = "High" if "High" in q else ("Med" if "Medium" in q else "Low")
            qual_tcol = "#FFFFFF"
            try:
                comp_val = float(m.get("checkv_completeness") or "")
                comp_str = f"{comp_val:.0f}%"
            except ValueError:
                comp_val = None
                comp_str = "—"
            if bp in BIOPROJECT_SOURCE:
                h_cat, h_col, h_txt = BIOPROJECT_SOURCE[bp]
            else:
                hl = (m.get("host") or lf.name).lower().replace("_", " ")
                if "ips" in hl and "confusus" in hl:
                    h_cat, h_col, h_txt = ("Ips confusus", "#1868DB", "I. confusus")
                elif "hylobius" in hl and "abietis" in hl:
                    h_cat, h_col, h_txt = ("Hylobius abietis", "#16794F", "H. abietis")
                elif "locusta" in hl and "migratoria" in hl:
                    h_cat, h_col, h_txt = ("Locusta migratoria", "#C75300", "L. migratoria")
                else:
                    h_cat, h_col, h_txt = ("Other", "#64748B", m.get("host", "Other")[:11] or "Other")
            cand_counts[(h_cat, h_col)] = cand_counts.get((h_cat, h_col), 0) + 1
            h_tip = f"{h_cat} ({bp or m.get('host', 'unknown')})"

        tx = px(lf.x)
        parts.append(f'<line x1="{tx + 2.8:.2f}" y1="{lf.y:.2f}" x2="{label_x - 4:.2f}" y2="{lf.y:.2f}" stroke="#CBD5E1" stroke-width="0.65" stroke-dasharray="1.5 2"/>')
        parts.append(f'<circle cx="{tx:.2f}" cy="{lf.y:.2f}" r="{2.3 if is_ref else 2.8}" fill="{col}" stroke="#FFFFFF" stroke-width="0.5"/>')
        lbl = short_label(lf.name, is_ref)
        parts.append(f'<text x="{label_x:.2f}" y="{lf.y + 2.5:.2f}" font-size="7.3" font-weight="{fw}" font-style="{fs}" fill="{tcol}">{html.escape(lbl)}</text>')

        # Host / Source tile (both ICTV References and Candidates)
        parts.append(f'<rect x="{tile_host}" y="{lf.y - 4.3:.2f}" width="68" height="8.6" rx="1.5" fill="{h_col}"><title>{html.escape(h_tip)}</title></rect>')
        parts.append(f'<text x="{tile_host + 34}" y="{lf.y + 2.3:.2f}" text-anchor="middle" font-size="6.4" font-weight="700" fill="#FFFFFF">{html.escape(h_txt)}</text>')

        # Origin / Quality tile
        parts.append(f'<rect x="{tile_qual}" y="{lf.y - 4.3:.2f}" width="30" height="8.6" rx="1.5" fill="{col}"/>')
        parts.append(f'<text x="{tile_qual + 15}" y="{lf.y + 2.3:.2f}" text-anchor="middle" font-size="6.2" font-weight="700" fill="{qual_tcol}">{qual_txt}</text>')

        # CheckV % tile
        c_bg = comp_colour(comp_val)
        c_txt = "#FFFFFF" if (comp_val is not None and comp_val >= 50) else INK
        parts.append(f'<rect x="{tile_comp}" y="{lf.y - 4.3:.2f}" width="32" height="8.6" rx="1.5" fill="{c_bg}"/>')
        parts.append(f'<text x="{tile_comp + 16}" y="{lf.y + 2.3:.2f}" text-anchor="middle" font-size="6.3" font-weight="700" fill="{c_txt}">{comp_str}</text>')

    # Scale bar
    sb_y = h - legend_h + 18
    sb_len = 0.5 * x_scale
    parts.append(f'<line x1="{x0}" y1="{sb_y}" x2="{x0 + sb_len:.1f}" y2="{sb_y}" stroke="{INK}" stroke-width="2.0" stroke-linecap="round"/>')
    parts.append(f'<text x="{x0 + sb_len / 2:.1f}" y="{sb_y + 13}" text-anchor="middle" font-size="8.8" font-weight="700" fill="{INK}">0.5 substitutions / site</text>')

    # Clean, spacious legend box at bottom (dynamically populated from tree taxa)
    box_y = h - legend_h + 38
    parts.append(f'<rect x="16" y="{box_y}" width="{w - 32}" height="{legend_h - 46}" rx="6" fill="#F8FAFC" stroke="#DDDEE1" stroke-width="1"/>')

    # Row 1: Origin / CheckV Quality + Node Support
    r1_y = box_y + 17
    parts.append(f'<text x="26" y="{r1_y}" font-size="8.2" font-weight="800" fill="{INK}">Origin &amp; Quality:</text>')
    leg1 = [
        ("Candidate (High-quality)", CAND_HIGH),
        ("Candidate (Medium-quality)", CAND_MED),
        ("Candidate (Low-quality)", CAND_LOW),
        (f"ICTV Reference (n = {n_refs})", REF_GREY),
        ("Node support: SH-aLRT ≥ 80 & UFB ≥ 95", SUPP_FILL),
    ]
    lx = 122
    for lbl, col in leg1:
        parts.append(f'<circle cx="{lx}" cy="{r1_y - 3}" r="3.5" fill="{col}"/>')
        parts.append(f'<text x="{lx + 7}" y="{r1_y}" font-size="8.2" fill="{INK}">{html.escape(lbl)}</text>')
        lx += len(lbl) * 4.5 + 22

    # Row 2: Candidate Host / Source
    r2_y = box_y + 37
    parts.append(f'<text x="26" y="{r2_y}" font-size="8.2" font-weight="800" fill="{INK}">Candidate Host / Source:</text>')
    lx = 160
    for (cat, col), cnt in sorted(cand_counts.items(), key=lambda kv: (-kv[1], kv[0][0])):
        lbl = f"{cat} (n = {cnt})"
        parts.append(f'<rect x="{lx - 3}" y="{r2_y - 7}" width="8" height="8" rx="1.5" fill="{col}"/>')
        parts.append(f'<text x="{lx + 8}" y="{r2_y}" font-size="8.2" fill="{INK}">{html.escape(lbl)}</text>')
        lx += len(lbl) * 4.35 + 18

    # Row 3: ICTV Reference Hosts
    r3_y = box_y + 57
    parts.append(f'<text x="26" y="{r3_y}" font-size="8.2" font-weight="800" fill="{INK}">ICTV Reference Host (GenBank):</text>')
    lx = 175
    for (cat, col), cnt in sorted(ref_counts.items(), key=lambda kv: (-kv[1], kv[0][0])):
        lbl = f"{cat} (n = {cnt})"
        parts.append(f'<rect x="{lx - 3}" y="{r3_y - 7}" width="8" height="8" rx="1.5" fill="{col}"/>')
        parts.append(f'<text x="{lx + 8}" y="{r3_y}" font-size="8.2" fill="{INK}">{html.escape(lbl)}</text>')
        lx += len(lbl) * 4.35 + 16

    parts.append("</svg>")
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_svg.write_text("\n".join(parts), encoding="utf-8")


def svg_to_png(svg_path: Path, png_path: Path, dpi: int = 300) -> None:
    import shutil, subprocess
    png_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import cairosvg  # type: ignore
        cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), dpi=dpi)
        return
    except Exception:
        pass
    if shutil.which("rsvg-convert"):
        subprocess.run(["rsvg-convert", "-d", str(dpi), "-p", str(dpi), "-o", str(png_path), str(svg_path)], check=True)
        return
    for cmd in ("magick", "convert"):
        if shutil.which(cmd):
            subprocess.run([cmd, "-density", str(dpi), str(svg_path), str(png_path)], check=True)
            return
    raise SystemExit("ERROR: --out-png requires cairosvg, rsvg-convert, or ImageMagick (magick/convert)")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree", required=True, type=Path, help="Input Newick tree file")
    p.add_argument("--meta", required=True, type=Path, help="Input selected.tsv metadata")
    p.add_argument("--out-svg", type=Path, help="Output SVG path")
    p.add_argument("--out-png", type=Path, help="Optional 300-dpi PNG output path")
    p.add_argument("--title", default="Hepeviridae Maximum-Likelihood Phylogeny")
    p.add_argument("--fetch-ncbi-hosts", action="store_true", help="Fetch /host qualifiers for ICTV references from NCBI")
    args = p.parse_args()
    if not args.out_svg and not args.out_png:
        p.error("Provide --out-svg and/or --out-png")
    svg_out = args.out_svg or args.out_png.with_suffix(".svg")
    render_tree_svg(args.tree, args.meta, svg_out, title=args.title, fetch_ncbi=args.fetch_ncbi_hosts)
    print(f"Wrote {svg_out}")
    if args.out_png:
        svg_to_png(svg_out, args.out_png)
        print(f"Wrote {args.out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
