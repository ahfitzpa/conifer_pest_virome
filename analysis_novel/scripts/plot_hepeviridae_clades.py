#!/usr/bin/env python3
"""Generate ICTV genus- and subfamily-specific circular (ggtree fan-style) ML
phylograms for Hepeviridae from a family-level IQ-TREE Newick file and
selected.tsv metadata.

Features:
  - Outgroup-roots on Piscihepevirus (Parahepevirinae) so classic ICTV genera
    (Rocahepevirus, Paslahepevirus, Chirohepevirus, Avihepevirus) and the
    Parahepevirinae / Aquatic-Environmental clade separate cleanly.
  - Assigns vertebrate candidates strictly to their nearest Orthohepevirinae
    ICTV genus by exact patristic distance (no cross-genus reference mixing).
  - Renders clean ggtree-style circular fan phylograms (no heatmap, no title),
    colour-coded by granular host (fetched from NCBI SRA/BioSample for SRR/ERR/DRR
    runs and NCBI GenBank /host for ICTV references when --fetch-ncbi-hosts is set).
  - Writes .treefile, .svg, and 300-dpi .png to workflow/analysis_novel/figures/Hepeviridae/.
"""
from __future__ import annotations

import argparse
import copy
import csv
import html
import json
import math
import re
import sys
import urllib.request
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import plot_phylo_tree as pt  # noqa: E402

INK = "#101214"
MUTED = "#505258"
BRANCH_COL = "#334155"

VERT_BIOPROJECTS = {
    bp for bp, (cat, _, _) in pt.BIOPROJECT_SOURCE.items()
    if "Vertebrate" in cat
}
AQUATIC_ENV_BIOPROJECTS = {
    bp for bp, (cat, _, _) in pt.BIOPROJECT_SOURCE.items()
    if any(k in cat for k in ("Marine", "Freshwater", "Crustacean", "Soil"))
}

ORTHO_GENERA: dict[str, tuple[str, ...]] = {
    "Rocahepevirus": ("rocahepevirus",),
    "Paslahepevirus": ("paslahepevirus",),
    "Chirohepevirus": ("chirohepevirus",),
    "Avihepevirus": ("avihepevirus",),
}

# Granular fallback host labels per Serratus BioProject (overridden by live SRA/BioSample host when fetched)
BIOPROJECT_GRANULAR_HOST: dict[str, tuple[str, str]] = {
    "PRJEB6505": ("Rodent / Rat (VIZIONS)", "#1868DB"),
    "PRJEB26687": ("Treeshrew / Rodent", "#2898BD"),
    "PRJEB17935": ("Rodent (Rat / Shrew)", "#1868DB"),
    "PRJNA751908": ("Bat", "#6E5DC6"),
    "PRJDB7591": ("Rodent (Wild)", "#1868DB"),
    "PRJNA683917": ("Rodent / Small Mammal", "#1868DB"),
    "PRJNA655509": ("Crustacean (Lobster/Crayfish)", "#C75300"),
    "PRJNA702424": ("Crustacean (Lobster/Crayfish)", "#C75300"),
    "PRJNA669582": ("Crustacean (Lobster/Crayfish)", "#C75300"),
    "PRJNA754774": ("Crustacean (Lobster/Crayfish)", "#C75300"),
    "PRJEB34490": ("Marine / Coastal", "#2898BD"),
    "PRJNA576751": ("Marine / Coastal", "#2898BD"),
    "PRJNA555787": ("Marine / Coastal", "#2898BD"),
    "PRJNA665326": ("Marine / Coastal", "#2898BD"),
    "PRJNA654739": ("Marine / Coastal", "#2898BD"),
    "PRJNA640812": ("Freshwater", "#16794F"),
    "PRJNA691101": ("Freshwater", "#16794F"),
    "PRJNA657047": ("Freshwater", "#16794F"),
    "PRJNA663589": ("Freshwater", "#16794F"),
    "PRJNA716119": ("Soil / Wastewater", "#8F4FDB"),
    "PRJNA616041": ("Soil / Wastewater", "#8F4FDB"),
    "PRJNA571457": ("Soil / Wastewater", "#8F4FDB"),
    "PRJNA673783": ("Soil / Wastewater", "#8F4FDB"),
}

GRANULAR_HOST_RULES: list[tuple[tuple[str, ...], str, str]] = [
    (("rattus norvegicus", "norway rat", "brown rat"), "Rat (R. norvegicus)", "#1868DB"),
    (("rattus tanezumi", "rattus rattus", "rattus", "bandicota", "niviventer", "maxomys", "berylmys"), "Rat (Rattus / Muridae)", "#1868DB"),
    (("suncus", "crocidura", "sorex", "shrew", "eulipotyphla"), "Shrew (Soricidae)", "#0D9488"),
    (("tupaia", "treeshrew", "scandentia"), "Treeshrew (Tupaia)", "#2898BD"),
    (("apodemus", "microtus", "myodes", "eothenomys", "mus ", "mouse", "vole", "hamster", "cricet", "gerbil", "rodent"), "Other Rodent", "#2563EB"),
    (("mustela", "ferret", "mink", "fox", "vulpes", "carniv", "cat", "dog"), "Carnivore / Ferret", "#D97706"),
    (("homo sapiens", "human"), "Human", "#E11D48"),
    (("sus scrofa", "swine", "pig", "boar"), "Swine", "#C75300"),
    (("oryctolagus", "rabbit", "lepus", "camel", "cervus", "deer", "macaque", "macaca"), "Other Mammal", "#B45309"),
    (("bat", "chiroptera", "eptesicus", "rhinolophus", "pipistrellus", "myotis", "desmodus", "rousettus"), "Bat", "#6E5DC6"),
    (("gallus", "chicken", "avian", "bird", "egretta", "sparrow", "pigeon", "duck", "goose", "falco", "kestrel"), "Bird / Avian", "#16794F"),
    (("oncorhynchus", "trout", "salmon", "fish", "pisci", "teleost", "cyprinus"), "Fish (Piscihepevirus)", "#0284C7"),
    (("homarus", "procambarus", "cherax", "lobster", "crayfish", "shrimp", "crab", "crustace", "decapod", "penaeus"), "Crustacean", "#C75300"),
    (("marine", "seawater", "coastal", "ocean", "estuar", "sediment", "bivalve", "oyster", "mussel"), "Marine / Coastal", "#2898BD"),
    (("freshwater", "river", "lake", "pond", "stream", "wetland"), "Freshwater", "#16794F"),
    (("soil", "wastewater", "sewage", "sludge", "compost"), "Soil / Wastewater", "#8F4FDB"),
]


def all_leaves(node: pt.Node) -> list[pt.Node]:
    if not node.children:
        return [node]
    out: list[pt.Node] = []
    for c in node.children:
        out.extend(all_leaves(c))
    return out


def set_parents(node: pt.Node, parent: pt.Node | None = None) -> None:
    node.parent = parent
    for c in node.children:
        set_parents(c, node)


def patristic_distances_from(root: pt.Node, start_leaf: pt.Node) -> dict[str, float]:
    set_parents(root)
    dists: dict[str, float] = {}
    stack: list[tuple[pt.Node, pt.Node | None, float]] = [(start_leaf, None, 0.0)]
    while stack:
        cur, prev, d = stack.pop()
        if not cur.children and cur is not start_leaf:
            dists[cur.name] = d
        if cur.parent is not None and cur.parent is not prev:
            stack.append((cur.parent, cur, d + cur.length))
        for c in cur.children:
            if c is not prev:
                stack.append((c, cur, d + c.length))
    return dists


def reroot_on_leaf(root: pt.Node, outgroup_name: str) -> pt.Node:
    set_parents(root)
    target = next((lf for lf in all_leaves(root) if lf.name == outgroup_name), None)
    if target is None or target.parent is None:
        return root
    if target.parent is root and len(root.children) == 2:
        return root

    path: list[pt.Node] = []
    cur: pt.Node | None = target
    while cur is not None:
        path.append(cur)
        cur = cur.parent

    p = path[1]
    p.children.remove(target)
    half = target.length / 2.0

    prev = p
    prev_len = p.length
    prev_name = p.name
    for idx in range(2, len(path)):
        anc = path[idx]
        anc.children.remove(prev)
        prev.children.append(anc)
        next_len = anc.length
        next_name = anc.name
        anc.length = prev_len
        anc.name = prev_name
        prev = anc
        prev_len = next_len
        prev_name = next_name

    old_root = path[-1]
    if len(old_root.children) == 1 and old_root is not p:
        only = old_root.children[0]
        only.length += old_root.length
        parent_of_old = path[-2]
        parent_of_old.children.remove(old_root)
        parent_of_old.children.append(only)

    new_root = pt.Node(name="", length=0.0, children=[target, p])
    target.length = half
    p.length = half
    p.name = ""
    set_parents(new_root)
    return new_root


def prune_tree(node: pt.Node, keep: set[str]) -> pt.Node | None:
    if not node.children:
        return pt.Node(name=node.name, length=node.length) if node.name in keep else None
    kept: list[pt.Node] = []
    for c in node.children:
        pc = prune_tree(c, keep)
        if pc is not None:
            kept.append(pc)
    if not kept:
        return None
    if len(kept) == 1:
        only = kept[0]
        only.length += node.length
        return only
    out = pt.Node(name=node.name, length=node.length, children=kept)
    for c in out.children:
        c.parent = out
    return out


def to_newick(node: pt.Node, is_root: bool = True) -> str:
    if not node.children:
        s = node.name
    else:
        inner = ",".join(to_newick(c, is_root=False) for c in node.children)
        s = f"({inner}){node.name}"
    if not is_root:
        s += f":{node.length:.6f}"
    return s + (";" if is_root else "")


def load_meta(meta_path: Path) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    with meta_path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = r.get("label") or r.get("taxon_id") or r.get("id") or ""
            if key:
                meta[key] = r
    return meta


def ref_genus(name: str) -> str | None:
    low = name.lower()
    if "piscihepevirus" in low or "cutthroat" in low or "parahepe" in low:
        return "Piscihepevirus"
    for genus, keys in ORTHO_GENERA.items():
        if any(k in low for k in keys):
            return genus
    return None


def is_candidate(name: str, meta: dict[str, dict[str, str]]) -> bool:
    return meta.get(name, {}).get("origin") != "reference"


def extract_sra_run(label: str) -> str:
    m = re.search(r"([SED]RR\d+)", label)
    return m.group(1) if m else ""


def fetch_sra_hosts(runs: list[str]) -> dict[str, str]:
    """Query NCBI SRA eSearch + eSummary for exact BioSample/Organism host strings."""
    uniq = sorted({r for r in runs if r})
    if not uniq:
        return {}
    out: dict[str, str] = {}
    for i in range(0, len(uniq), 25):
        batch = uniq[i:i + 25]
        term = "+OR+".join(f"{r}[Accession]" for r in batch)
        try:
            url1 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode=json&retmax=100&term={term}"
            with urllib.request.urlopen(url1, timeout=20) as resp:
                ids = json.loads(resp.read().decode("utf-8")).get("esearchresult", {}).get("idlist", [])
            if not ids:
                continue
            url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&retmode=json&id={','.join(ids)}"
            with urllib.request.urlopen(url2, timeout=25) as resp:
                data = json.loads(resp.read().decode("utf-8")).get("result", {})
            for uid in ids:
                rec = data.get(uid, {})
                expxml = html.unescape(rec.get("expxml", ""))
                runs_xml = html.unescape(rec.get("runs", ""))
                m_org = re.search(r'ScientificName="([^"]+)"', expxml)
                m_title = re.search(r"<Title>(.*?)</Title>", expxml)
                org = m_org.group(1).strip() if m_org else ""
                title = m_title.group(1).strip() if m_title else ""
                # Prefer specific organism; if metagenome, append title clue
                desc = f"{org} {title}".strip() if "metagenome" in org.lower() else (org or title)
                for r_acc in re.findall(r'acc="([SED]RR\d+)"', runs_xml):
                    if desc and r_acc not in out:
                        out[r_acc] = desc
        except Exception:
            continue
    return out


def classify_granular_host(
    label: str,
    m: dict[str, str],
    is_ref: bool,
    ref_hosts: dict[str, tuple[str, str]],
    sra_hosts: dict[str, str],
) -> tuple[str, str, str]:
    """Return (host_category, hex_colour, short_tag_for_tip)."""
    if is_ref:
        acc = pt.extract_accession(label)
        h_str, o_str = ref_hosts.get(acc, (m.get("host", ""), ""))
        blob = f"{h_str} {o_str} {label}".lower()
        for keys, cat, col in GRANULAR_HOST_RULES:
            if any(k in blob for k in keys):
                return cat, col, h_str or cat
        return "Other Reference", "#64748B", h_str or "Ref"

    run = extract_sra_run(label)
    sra_desc = sra_hosts.get(run, "")
    bp = m.get("bioproject", "")
    meta_host = m.get("host", "") or m.get("scientificname", "")
    blob = f"{sra_desc} {meta_host} {label}".lower()
    for keys, cat, col in GRANULAR_HOST_RULES:
        if any(k in blob for k in keys):
            # Build concise specific tag if SRA gave a binomial species name
            tag = sra_desc.split()[0] + " " + sra_desc.split()[1] if len(sra_desc.split()) >= 2 and "metagenome" not in sra_desc.lower() else cat
            return cat, col, tag
    if bp in BIOPROJECT_GRANULAR_HOST:
        cat, col = BIOPROJECT_GRANULAR_HOST[bp]
        return cat, col, cat
    return "Other / Unclassified", "#64748B", meta_host or "Other"


def classify_ictv_clades(root: pt.Node, meta: dict[str, dict[str, str]]) -> tuple[dict[str, set[str]], pt.Node]:
    leaves = all_leaves(root)
    pisci_refs = [lf.name for lf in leaves if not is_candidate(lf.name, meta) and ref_genus(lf.name) == "Piscihepevirus"]
    if pisci_refs:
        root = reroot_on_leaf(root, sorted(pisci_refs)[0])
        leaves = all_leaves(root)

    genus_refs: dict[str, set[str]] = {g: set() for g in ORTHO_GENERA}
    ref_to_genus: dict[str, str] = {}
    vert_tips: set[str] = set()
    aquatic_tips: set[str] = set()
    vert_cand_leaves: list[pt.Node] = []

    for lf in leaves:
        m = meta.get(lf.name, {})
        if not is_candidate(lf.name, meta):
            g = ref_genus(lf.name)
            if g in ORTHO_GENERA:
                genus_refs[g].add(lf.name)
                ref_to_genus[lf.name] = g
                vert_tips.add(lf.name)
            else:
                aquatic_tips.add(lf.name)
        else:
            bp = m.get("bioproject", "")
            if bp in VERT_BIOPROJECTS:
                vert_tips.add(lf.name)
                vert_cand_leaves.append(lf)
            elif bp in AQUATIC_ENV_BIOPROJECTS:
                aquatic_tips.add(lf.name)
            else:
                hl = (m.get("host") or "").lower()
                if any(k in hl for k in ("rat", "rodent", "mouse", "bat", "bird", "swine", "pig", "human", "treeshrew", "mammal")):
                    vert_tips.add(lf.name)
                    vert_cand_leaves.append(lf)
                else:
                    aquatic_tips.add(lf.name)

    genus_cands: dict[str, set[str]] = {g: set() for g in ORTHO_GENERA}
    for lf in vert_cand_leaves:
        dists = patristic_distances_from(root, lf)
        nearest_ref = min(
            (r for r in ref_to_genus if r in dists),
            key=lambda r: dists[r],
            default=None,
        )
        if nearest_ref is not None:
            genus_cands[ref_to_genus[nearest_ref]].add(lf.name)

    clades: dict[str, set[str]] = {}
    for genus, refs in genus_refs.items():
        if not refs:
            continue
        clades[f"hepeviridae_genus_{genus}"] = set(refs) | genus_cands[genus]

    clades["hepeviridae_subfamily_Orthohepevirinae"] = vert_tips
    clades["hepeviridae_subfamily_Parahepevirinae_aquatic"] = aquatic_tips
    return clades, root


def annular_wedge_path(cx: float, cy: float, r_in: float, r_out: float, a0: float, a1: float) -> str:
    """Return SVG path for an annular sector (ggtreeExtra geom_fruit ring tile)."""
    x0_in, y0_in = cx + r_in * math.cos(a0), cy + r_in * math.sin(a0)
    x1_in, y1_in = cx + r_in * math.cos(a1), cy + r_in * math.sin(a1)
    x0_out, y0_out = cx + r_out * math.cos(a0), cy + r_out * math.sin(a0)
    x1_out, y1_out = cx + r_out * math.cos(a1), cy + r_out * math.sin(a1)
    large = 1 if abs(a1 - a0) > math.pi else 0
    return (
        f"M {x0_out:.2f} {y0_out:.2f} "
        f"A {r_out:.2f} {r_out:.2f} 0 {large} 1 {x1_out:.2f} {y1_out:.2f} "
        f"L {x1_in:.2f} {y1_in:.2f} "
        f"A {r_in:.2f} {r_in:.2f} 0 {large} 0 {x0_in:.2f} {y0_in:.2f} Z"
    )


def render_circular_tree_svg(
    tree_path: Path,
    meta: dict[str, dict[str, str]],
    out_svg: Path,
    ref_hosts: dict[str, tuple[str, str]],
    sra_hosts: dict[str, str],
    show_labels: bool = False,
) -> None:
    """Render a ggtree / ggtreeExtra-style open-fan circular phylogram SVG with
    concentric annular Host & Origin rings, 3-o'clock ring labels, and a clean
    right-side vertical legend (no title, no CheckV heatmap, no tip labels by default).
    """
    tree = pt.parse_newick(tree_path.read_text(encoding="utf-8"))
    pt.ladderize(tree, reverse=False)
    leaves = pt.get_leaves(tree)
    n_leaves = len(leaves)

    # Open fan centered at 0 deg (3 o'clock) so ring labels sit cleanly in the gap
    open_deg = 22.0
    start_deg = open_deg / 2.0
    span_deg = 360.0 - open_deg
    step_rad = math.radians(span_deg / max(1, n_leaves))
    leaf_angle: dict[int, float] = {}
    for idx, lf in enumerate(leaves):
        leaf_angle[id(lf)] = math.radians(start_deg) + (idx + 0.5) * step_rad

    node_r: dict[int, float] = {}
    node_a: dict[int, float] = {}

    def assign_polar(node: pt.Node, depth: float = 0.0) -> tuple[float, float]:
        r = depth + node.length
        node_r[id(node)] = r
        if not node.children:
            a = leaf_angle[id(node)]
            node_a[id(node)] = a
            return a, a
        angles = [assign_polar(c, r) for c in node.children]
        lo = min(x[0] for x in angles)
        hi = max(x[1] for x in angles)
        node_a[id(node)] = 0.5 * (lo + hi)
        return lo, hi

    assign_polar(tree, 0.0)
    max_r = max(node_r[id(lf)] for lf in leaves) or 1.0

    # Concentric radii (clean ggtreeExtra proportions without label clutter)
    r_inner = 12.0
    tree_r = 215.0 if n_leaves <= 35 else 235.0
    r_scale = (tree_r - r_inner) / max_r
    outer_tree_r = r_inner + max_r * r_scale

    ring1_in = outer_tree_r + 10.0
    ring1_out = ring1_in + 22.0   # Ring 1: Granular Host / Source
    ring2_in = ring1_out + 4.0
    ring2_out = ring2_in + 14.0   # Ring 2: Origin (Candidate vs. ICTV Reference)
    label_r = ring2_out + 8.0
    outer_pad = 148.0 if show_labels else 24.0

    cx = ring2_out + outer_pad + 18.0
    cy = ring2_out + outer_pad + 14.0
    legend_w = 245.0
    w = int(cx + ring2_out + outer_pad + legend_w)
    h = int(2 * cy + 20)

    def polar_xy(depth: float, ang: float) -> tuple[float, float]:
        rr = r_inner + depth * r_scale
        return cx + rr * math.cos(ang), cy + rr * math.sin(ang)

    # Pre-classify all tips
    tip_info: dict[int, tuple[bool, str, str, str]] = {}
    host_counts: dict[tuple[str, str], int] = {}
    n_cand = 0
    n_ref = 0
    for lf in leaves:
        m = meta.get(lf.name, {})
        is_ref = m.get("origin") == "reference"
        if is_ref:
            n_ref += 1
        else:
            n_cand += 1
        h_cat, h_col, h_tag = classify_granular_host(lf.name, m, is_ref, ref_hosts, sra_hosts)
        host_counts[(h_cat, h_col)] = host_counts.get((h_cat, h_col), 0) + 1
        base_lbl = pt.short_label(lf.name, is_ref)
        if not is_ref and h_tag and h_tag.lower() not in base_lbl.lower():
            disp = f"{base_lbl} [{h_tag[:18]}]"
        else:
            disp = base_lbl
        tip_info[id(lf)] = (is_ref, h_cat, h_col, disp)

    sorted_hosts = sorted(host_counts.items(), key=lambda kv: (-kv[1], kv[0][0]))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Inter,Arial,Helvetica,sans-serif">',
        f'<rect width="{w}" height="{h}" fill="#FFFFFF"/>',
    ]

    # 1. Draw circular tree branches (clean ggtree black strokes)
    def draw_node(node: pt.Node) -> None:
        r0 = node_r[id(node)]
        rr0 = r_inner + r0 * r_scale
        if node.children:
            child_angs = [node_a[id(c)] for c in node.children]
            a_min, a_max = min(child_angs), max(child_angs)
            x1, y1 = cx + rr0 * math.cos(a_min), cy + rr0 * math.sin(a_min)
            x2, y2 = cx + rr0 * math.cos(a_max), cy + rr0 * math.sin(a_max)
            large_arc = 1 if (a_max - a_min) > math.pi else 0
            if rr0 > 1.0 and abs(a_max - a_min) > 1e-4:
                parts.append(
                    f'<path d="M {x1:.2f} {y1:.2f} A {rr0:.2f} {rr0:.2f} 0 {large_arc} 1 {x2:.2f} {y2:.2f}" '
                    f'fill="none" stroke="{INK}" stroke-width="1.25" stroke-linecap="round"/>'
                )
            for c in node.children:
                ac = node_a[id(c)]
                rc = node_r[id(c)]
                cx1, cy1 = polar_xy(r0, ac)
                cx2, cy2 = polar_xy(rc, ac)
                parts.append(
                    f'<line x1="{cx1:.2f}" y1="{cy1:.2f}" x2="{cx2:.2f}" y2="{cy2:.2f}" '
                    f'stroke="{INK}" stroke-width="1.25" stroke-linecap="round"/>'
                )
                draw_node(c)
            alrt, ufb = pt.parse_support(node.name)
            if alrt is not None and ufb is not None and alrt >= 80 and ufb >= 95:
                nx, ny = polar_xy(r0, node_a[id(node)])
                parts.append(f'<circle cx="{nx:.2f}" cy="{ny:.2f}" r="2.6" fill="{INK}" stroke="#FFFFFF" stroke-width="0.5"/>')

    draw_node(tree)

    # 2. Draw tip dots, dashed guides, concentric Host & Origin rings (optional outer labels)
    font_sz = 8.2 if n_leaves <= 28 else (7.0 if n_leaves <= 55 else 5.8)
    half_step = step_rad * 0.485
    cand_ring_col = "#1C2B42"
    ref_ring_col = "#CBD5E1"

    for lf in leaves:
        is_ref, h_cat, h_col, disp = tip_info[id(lf)]
        ang = node_a[id(lf)]
        deg = (math.degrees(ang)) % 360.0
        tx, ty = polar_xy(node_r[id(lf)], ang)

        gx, gy = cx + (ring1_in - 2.0) * math.cos(ang), cy + (ring1_in - 2.0) * math.sin(ang)
        parts.append(
            f'<line x1="{tx:.2f}" y1="{ty:.2f}" x2="{gx:.2f}" y2="{gy:.2f}" '
            f'stroke="#CBD5E1" stroke-width="0.7" stroke-dasharray="1.6 2.0"/>'
        )
        if is_ref:
            parts.append(f'<circle cx="{tx:.2f}" cy="{ty:.2f}" r="2.7" fill="{h_col}" stroke="#FFFFFF" stroke-width="0.5"/>')
        else:
            parts.append(f'<circle cx="{tx:.2f}" cy="{ty:.2f}" r="4.0" fill="{h_col}" stroke="#FFFFFF" stroke-width="0.8"/>')

        a0, a1 = ang - half_step, ang + half_step
        parts.append(
            f'<path d="{annular_wedge_path(cx, cy, ring1_in, ring1_out, a0, a1)}" '
            f'fill="{h_col}" stroke="#FFFFFF" stroke-width="0.6"/>'
        )
        orig_col = ref_ring_col if is_ref else cand_ring_col
        parts.append(
            f'<path d="{annular_wedge_path(cx, cy, ring2_in, ring2_out, a0, a1)}" '
            f'fill="{orig_col}" stroke="#FFFFFF" stroke-width="0.6"/>'
        )

        if show_labels:
            lx, ly = cx + label_r * math.cos(ang), cy + label_r * math.sin(ang)
            if 90.0 < deg < 270.0:
                rot = deg - 180.0
                anchor = "end"
            else:
                rot = deg
                anchor = "start"
            fw = "400" if is_ref else "700"
            fs = "italic" if is_ref else "normal"
            tcol = MUTED if is_ref else INK
            parts.append(
                f'<text x="{lx:.2f}" y="{ly:.2f}" transform="rotate({rot:.2f} {lx:.2f} {ly:.2f})" '
                f'text-anchor="{anchor}" dominant-baseline="central" font-size="{font_sz}" '
                f'font-weight="{fw}" font-style="{fs}" fill="{tcol}">{html.escape(disp)}</text>'
            )

    # Ring column labels inside the 3-o'clock open wedge
    r1_mid = 0.5 * (ring1_in + ring1_out)
    r2_mid = 0.5 * (ring2_in + ring2_out)
    parts.append(
        f'<text x="{cx + r1_mid:.1f}" y="{cy:.1f}" transform="rotate(-90 {cx + r1_mid:.1f} {cy:.1f})" '
        f'text-anchor="middle" dominant-baseline="central" font-size="9.5" font-weight="800" fill="{INK}">Host</text>'
    )
    parts.append(
        f'<text x="{cx + r2_mid:.1f}" y="{cy:.1f}" transform="rotate(-90 {cx + r2_mid:.1f} {cy:.1f})" '
        f'text-anchor="middle" dominant-baseline="central" font-size="9.5" font-weight="800" fill="{INK}">Origin</text>'
    )

    # Adaptive scale bar at bottom-left of circle (left-aligned so text never clips)
    sb_val = 0.1 if max_r < 0.45 else (0.2 if max_r < 1.1 else (0.5 if max_r < 3.0 else 1.0))
    sb_px = sb_val * r_scale
    sb_x = 24.0
    sb_y = h - 30.0
    parts.append(f'<line x1="{sb_x:.1f}" y1="{sb_y:.1f}" x2="{sb_x + sb_px:.1f}" y2="{sb_y:.1f}" stroke="{INK}" stroke-width="2.2" stroke-linecap="round"/>')
    parts.append(f'<text x="{sb_x:.1f}" y="{sb_y + 15:.1f}" text-anchor="start" font-size="9.4" font-weight="700" fill="{INK}">{sb_val:g} substitutions / site</text>')

    # 3. Clean ggplot2 / ggtree right-side vertical legend
    leg_x = cx + ring2_out + outer_pad + 16.0
    leg_y = max(60.0, cy - (len(sorted_hosts) * 19 + 145) / 2.0)

    parts.append(f'<text x="{leg_x:.1f}" y="{leg_y:.1f}" font-size="11.0" font-weight="800" fill="{INK}">Host / Source</text>')
    cur_y = leg_y + 20.0
    for (cat, col), cnt in sorted_hosts:
        parts.append(f'<rect x="{leg_x:.1f}" y="{cur_y - 9.5:.1f}" width="12" height="12" rx="2" fill="{col}"/>')
        parts.append(f'<text x="{leg_x + 18:.1f}" y="{cur_y:.1f}" font-size="9.6" fill="{INK}">{html.escape(cat)} (n = {cnt})</text>')
        cur_y += 19.0

    cur_y += 15.0
    parts.append(f'<text x="{leg_x:.1f}" y="{cur_y:.1f}" font-size="11.0" font-weight="800" fill="{INK}">Origin (Outer Ring)</text>')
    cur_y += 20.0
    parts.append(f'<rect x="{leg_x:.1f}" y="{cur_y - 9.5:.1f}" width="12" height="12" rx="2" fill="{cand_ring_col}"/>')
    parts.append(f'<text x="{leg_x + 18:.1f}" y="{cur_y:.1f}" font-size="9.6" font-weight="700" fill="{INK}">Candidate (n = {n_cand})</text>')
    cur_y += 19.0
    parts.append(f'<rect x="{leg_x:.1f}" y="{cur_y - 9.5:.1f}" width="12" height="12" rx="2" fill="{ref_ring_col}"/>')
    parts.append(f'<text x="{leg_x + 18:.1f}" y="{cur_y:.1f}" font-size="9.6" font-style="italic" fill="{MUTED}">ICTV Reference (n = {n_ref})</text>')

    cur_y += 23.0
    parts.append(f'<text x="{leg_x:.1f}" y="{cur_y:.1f}" font-size="11.0" font-weight="800" fill="{INK}">Node Support</text>')
    cur_y += 18.0
    parts.append(f'<circle cx="{leg_x + 6.0:.1f}" cy="{cur_y - 3.5:.1f}" r="3.2" fill="{INK}"/>')
    parts.append(f'<text x="{leg_x + 18:.1f}" y="{cur_y:.1f}" font-size="9.4" fill="{INK}">SH-aLRT ≥ 80 &amp; UFB ≥ 95</text>')

    parts.append("</svg>")
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_svg.write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", required=True, type=Path, help="Full Hepeviridae IQ-TREE .treefile")
    ap.add_argument("--meta", required=True, type=Path, help="Hepeviridae selected.tsv metadata")
    ap.add_argument(
        "--out-dir",
        default="workflow/analysis_novel/figures/Hepeviridae",
        type=Path,
        help="Output directory (default: workflow/analysis_novel/figures/Hepeviridae)",
    )
    ap.add_argument("--layout", choices=("circular", "linear"), default="circular", help="Tree layout (default: circular)")
    ap.add_argument("--show-labels", action="store_true", help="Show outer radial tip labels on circular trees")
    ap.add_argument("--fetch-ncbi-hosts", action="store_true", help="Fetch ICTV GenBank /host and SRA BioSample hosts from NCBI")
    ap.add_argument("--candidates-only", action="store_true", help="Only render clades containing >= 1 candidate")
    ap.add_argument("--no-png", action="store_true", help="Skip 300-dpi PNG conversion (write SVG only)")
    args = ap.parse_args()

    raw_root = pt.parse_newick(args.tree.read_text(encoding="utf-8"))
    meta = load_meta(args.meta)
    clades, rooted = classify_ictv_clades(raw_root, meta)

    # Fetch NCBI hosts once for all ICTV refs and SRA candidate runs
    ref_hosts: dict[str, tuple[str, str]] = {}
    sra_hosts: dict[str, str] = {}
    if args.fetch_ncbi_hosts:
        leaves = all_leaves(rooted)
        ref_accs = [pt.extract_accession(lf.name) for lf in leaves if not is_candidate(lf.name, meta)]
        sra_runs = [extract_sra_run(lf.name) for lf in leaves if is_candidate(lf.name, meta)]
        if ref_accs:
            ref_hosts = pt.fetch_ncbi_hosts(ref_accs)
        if sra_runs:
            sra_hosts = fetch_sra_hosts(sra_runs)
            print(f"[ncbi] resolved granular SRA/BioSample host for {len(sra_hosts)}/{len(set(r for r in sra_runs if r))} candidate run(s)")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for stem, keep in clades.items():
        n_cand = sum(1 for t in keep if is_candidate(t, meta))
        n_ref = len(keep) - n_cand
        if len(keep) < 3:
            print(f"[skip] {stem}: only {len(keep)} tip(s)")
            continue
        if args.candidates_only and n_cand == 0:
            print(f"[skip] {stem}: 0 candidates ({n_ref} ICTV refs)")
            continue
        pruned = prune_tree(copy.deepcopy(rooted), keep)
        if pruned is None:
            continue
        pruned.length = 0.0
        nwk_path = args.out_dir / f"{stem}.treefile"
        svg_path = args.out_dir / f"{stem}.svg"
        png_path = args.out_dir / f"{stem}.png"
        nwk_path.write_text(to_newick(pruned) + "\n", encoding="utf-8")
        if args.layout == "circular":
            render_circular_tree_svg(nwk_path, meta, svg_path, ref_hosts, sra_hosts, show_labels=args.show_labels)
        else:
            pt.render_tree_svg(nwk_path, args.meta, svg_path, title="", fetch_ncbi=args.fetch_ncbi_hosts)
        print(f"Wrote {svg_path} ({len(keep)} taxa: {n_cand} candidates, {n_ref} ICTV refs)")
        if not args.no_png:
            pt.svg_to_png(svg_path, png_path, dpi=300)
            print(f"Wrote {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
