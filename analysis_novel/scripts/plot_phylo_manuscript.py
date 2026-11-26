#!/usr/bin/env python3
"""Render publication-grade phylogenetic manuscript figures (Figs. 4, 5, 6) as pure SVG.

Zero-OpenBLAS, pure-Python stdlib renderer designed to run safely on HPC login
nodes and produce clean vector SVGs from Stage 09 (`results/phylo/families/*`)
and Stage 10 (`results/genome_maps/*/genome_map_data.tsv`) outputs.

Modes:
  fig4  Pan-Holarctic Ips Polycipiviridae: midpoint-rooted AA ML phylogeny
        paired with a to-scale ORF genome-architecture track.
  fig5  Closed-Genome Iflaviridae Clade: midpoint-rooted AA ML phylogeny
        coloured by host species and continent.
  fig6  Segmented Jingmen-like Amarillovirales & Negative-Sense RNA Viruses:
        3-panel composite AA ML phylogeny (Amarillovirales, Phenuiviridae,
        Aliusviridae).
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
LOG = logging.getLogger("plot_phylo_manuscript")

# Okabe-Ito-inspired host palette consistent across Figs. 1-6.
HOST_STYLE: Dict[str, Tuple[str, str, str]] = {
    "Ips confusus": ("#0072B2", "North America", "I. confusus (NA)"),
    "Hylobius abietis": ("#009E73", "Europe", "H. abietis (EU)"),
    "Ips typographus": ("#56B4E9", "Europe", "I. typographus (EU)"),
    "Ips nitidus": ("#D55E00", "East Asia", "I. nitidus (EA)"),
    "Ips subelongatus": ("#E69F00", "East Asia", "I. subelongatus (EA)"),
}
ICTV_COLOR = "#64748B"
OTHER_HOST_COLOR = "#7C3AED"

ORF_COLOURS: List[Tuple[str, str, str]] = [
    ("rdrp|polymerase|replicase|rna-dependent", "#DC2626", "RdRp / Polymerase"),
    ("helicase|ntpase|2c", "#2563EB", "Helicase"),
    ("protease|peptidase|3c", "#D97706", "Protease"),
    ("capsid|coat|vp1|vp2|vp3|vp4|structural|glycoprotein|nucleocapsid", "#059669", "Structural / Capsid"),
    ("methyltransferase|mtase|capping", "#7C3AED", "Methyltransferase"),
]
ORF_DEFAULT_COLOUR = "#94A3B8"


# ---------------------------------------------------------------------------
# Newick parser, midpoint rooter, and ladderizer (pure stdlib)
# ---------------------------------------------------------------------------
@dataclass
class TreeNode:
    name: str = ""
    length: float = 0.0
    children: List["TreeNode"] = field(default_factory=list)
    parent: Optional["TreeNode"] = None
    support_alrt: Optional[float] = None
    support_ufb: Optional[float] = None

    def is_leaf(self) -> bool:
        return len(self.children) == 0


def _parse_support(label: str) -> Tuple[Optional[float], Optional[float]]:
    label = (label or "").strip().strip("'\"")
    if not label:
        return None, None
    if "/" in label:
        parts = label.split("/", 1)
        try:
            return float(parts[0]), float(parts[1])
        except ValueError:
            return None, None
    try:
        val = float(label)
        return None, val
    except ValueError:
        return None, None


def parse_newick(text: str) -> TreeNode:
    """Parse a Newick string (including IQ-TREE dual support labels) into a TreeNode."""
    s = text.strip()
    if not s:
        raise ValueError("Empty Newick tree")
    if s.endswith(";"):
        s = s[:-1].strip()

    tokens: List[str] = []
    i, n = 0, len(s)
    while i < n:
        ch = s[i]
        if ch in "(),:":
            tokens.append(ch)
            i += 1
        elif ch in ("'", '"'):
            quote = ch
            j = i + 1
            buf = []
            while j < n and s[j] != quote:
                buf.append(s[j])
                j += 1
            tokens.append("".join(buf))
            i = j + 1
        elif ch.isspace():
            i += 1
        else:
            j = i
            while j < n and s[j] not in "(),:" and not s[j].isspace():
                j += 1
            tokens.append(s[i:j])
            i = j

    pos = 0

    def parse_subtree(parent: Optional[TreeNode]) -> TreeNode:
        nonlocal pos
        node = TreeNode(parent=parent)
        if pos < len(tokens) and tokens[pos] == "(":
            pos += 1
            while True:
                child = parse_subtree(node)
                node.children.append(child)
                if pos < len(tokens) and tokens[pos] == ",":
                    pos += 1
                    continue
                if pos < len(tokens) and tokens[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Malformed Newick near token {pos}")
        if pos < len(tokens) and tokens[pos] not in "(),:":
            node.name = tokens[pos]
            pos += 1
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            if pos < len(tokens):
                try:
                    node.length = max(0.0, float(tokens[pos]))
                except ValueError:
                    node.length = 0.0
                pos += 1
        if not node.is_leaf():
            node.support_alrt, node.support_ufb = _parse_support(node.name)
        return node

    root = parse_subtree(None)
    return root


def _all_leaves(node: TreeNode) -> List[TreeNode]:
    if node.is_leaf():
        return [node]
    out: List[TreeNode] = []
    for c in node.children:
        out.extend(_all_leaves(c))
    return out


def _prune_leaves(root: TreeNode, drop_names: set) -> TreeNode:
    """Remove leaves whose names match `drop_names`, collapsing unary internal nodes."""
    if not drop_names:
        return root

    def recurse(node: TreeNode) -> Optional[TreeNode]:
        if node.is_leaf():
            return None if node.name in drop_names else node
        kept = []
        for c in node.children:
            res = recurse(c)
            if res is not None:
                kept.append(res)
        if not kept:
            return None
        if len(kept) == 1:
            only = kept[0]
            only.length += node.length
            return only
        node.children = kept
        for c in kept:
            c.parent = node
        return node

    pruned = recurse(root)
    if pruned is None:
        return root
    pruned.parent = None
    return pruned


def midpoint_root(root: TreeNode) -> TreeNode:
    """Midpoint-root a tree by finding the longest tip-to-tip path."""
    leaves = _all_leaves(root)
    if len(leaves) < 3:
        return root

    # Build undirected adjacency list
    adj: Dict[int, List[Tuple[int, float]]] = defaultdict(list)
    id_to_node: Dict[int, TreeNode] = {}

    def walk(u: TreeNode) -> None:
        id_to_node[id(u)] = u
        for c in u.children:
            w = max(c.length, 1e-6)
            adj[id(u)].append((id(c), w))
            adj[id(c)].append((id(u), w))
            walk(c)

    walk(root)

    def farthest(start_id: int) -> Tuple[int, float, Dict[int, Tuple[Optional[int], float]]]:
        dist: Dict[int, float] = {start_id: 0.0}
        prev: Dict[int, Tuple[Optional[int], float]] = {start_id: (None, 0.0)}
        stack = [start_id]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + w
                    prev[v] = (u, w)
                    stack.append(v)
        best = max(dist.items(), key=lambda kv: kv[1])
        return best[0], best[1], prev

    end1, _, _ = farthest(id(leaves[0]))
    end2, max_dist, prev = farthest(end1)
    if max_dist <= 0:
        return root

    # Reconstruct path end1 -> end2
    path: List[Tuple[int, float]] = []
    cur: Optional[int] = end2
    while cur is not None:
        p, w = prev[cur]
        path.append((cur, w))
        cur = p
    path.reverse()

    half = max_dist / 2.0
    acc = 0.0
    u_id = path[0][0]
    v_id = path[1][0]
    edge_w = path[1][1]
    split_from_u = half
    for idx in range(1, len(path)):
        nid, w = path[idx]
        if acc + w >= half:
            u_id = path[idx - 1][0]
            v_id = nid
            edge_w = w
            split_from_u = half - acc
            break
        acc += w

    # Rebuild rooted tree from a new root inserted between u_id and v_id
    new_root = TreeNode(name="root", length=0.0)

    def build_directed(curr_id: int, blocked_id: int, parent_node: TreeNode, branch_len: float) -> TreeNode:
        orig = id_to_node[curr_id]
        copy = TreeNode(
            name=orig.name,
            length=branch_len,
            parent=parent_node,
            support_alrt=orig.support_alrt,
            support_ufb=orig.support_ufb,
        )
        for nxt_id, w in adj[curr_id]:
            if nxt_id == blocked_id:
                continue
            child = build_directed(nxt_id, curr_id, copy, w)
            copy.children.append(child)
        # Collapse degree-2 old root if it has no name/leaf status
        if len(copy.children) == 1 and not orig.is_leaf():
            only = copy.children[0]
            only.length += copy.length
            only.parent = parent_node
            return only
        return copy

    left = build_directed(u_id, v_id, new_root, max(0.0, split_from_u))
    right = build_directed(v_id, u_id, new_root, max(0.0, edge_w - split_from_u))
    new_root.children = [left, right]
    return new_root


def ladderize(node: TreeNode) -> int:
    """Sort children by descending leaf count for clean visual presentation."""
    if node.is_leaf():
        return 1
    counts = [(ladderize(c), c) for c in node.children]
    counts.sort(key=lambda t: t[0], reverse=True)
    node.children = [c for _, c in counts]
    return sum(n for n, _ in counts)


# ---------------------------------------------------------------------------
# Metadata & ORF loaders
# ---------------------------------------------------------------------------
def load_manifest_hosts(manifest_path: Optional[Path]) -> Dict[str, str]:
    hosts: Dict[str, str] = {}
    if not manifest_path or not manifest_path.is_file():
        return hosts
    with manifest_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            run = (row.get("run") or "").strip()
            host = (row.get("host_organism") or row.get("host") or "").strip()
            if run and host:
                hosts[run] = host
    return hosts


def load_taxa_metadata(meta_path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    if not meta_path or not meta_path.is_file():
        return out
    with meta_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tid = (row.get("taxon_id") or "").strip()
            if tid:
                out[tid] = row
    return out


def load_genome_map_data(
    table_path: Optional[Path],
) -> Dict[str, Tuple[int, List[Dict[str, str]]]]:
    """Return contig_id -> (contig_length, list of ORF dicts)."""
    by_contig: Dict[str, Tuple[int, List[Dict[str, str]]]] = {}
    if not table_path or not table_path.is_file():
        return by_contig
    with table_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("contig_id") or "").strip()
            if not cid:
                continue
            length = int(float(row.get("contig_length") or 0))
            cur_len, orfs = by_contig.get(cid, (length, []))
            cur_len = max(cur_len, length)
            start = int(float(row.get("ORF_start") or 0))
            end = int(float(row.get("ORF_end") or 0))
            if end > start:
                orfs.append(row)
            by_contig[cid] = (cur_len, orfs)
    return by_contig


def normalise_host(raw: str) -> str:
    low = (raw or "").strip().lower()
    for sp in HOST_STYLE:
        if sp.lower() in low:
            return sp
    return (raw or "").strip()


def classify_tip(
    tip_name: str,
    meta: Dict[str, Dict[str, str]],
    manifest_hosts: Dict[str, str],
) -> Tuple[bool, str, str, str, str]:
    """Return (is_candidate, host_species, display_label, color, contig_key)."""
    row = meta.get(tip_name, {})
    group = (row.get("group") or "").strip().lower()
    is_cand = group in {"candidate", "query", "novel"} or tip_name.startswith("cand|")

    # Extract run & contig from cand|<run>|<contig> or <run>_NODE_...
    run_id = ""
    contig_id = ""
    if tip_name.startswith("cand|"):
        parts = tip_name.split("|")
        run_id = parts[1] if len(parts) > 1 else ""
        contig_id = parts[2] if len(parts) > 2 else ""
    else:
        m = re.match(r"((?:SRR|ERR|DRR)\d+)[_|](.+)", tip_name)
        if m:
            run_id, contig_id = m.group(1), m.group(2)
            is_cand = True

    if is_cand:
        host = normalise_host(row.get("host") or manifest_hosts.get(run_id, ""))
        color, region, short_host = HOST_STYLE.get(
            host, (OTHER_HOST_COLOR, "", host or "Candidate")
        )
        # Compact, human-readable label (no ugly full NODE_..._cov_... strings)
        node_m = re.search(r"NODE_(\d+)_length_(\d+)", contig_id)
        if node_m:
            kb = int(node_m.group(2)) / 1000.0
            label = f"{short_host} · {run_id} ({kb:.1f} kb)"
        elif run_id:
            label = f"{short_host} · {run_id}"
        else:
            label = tip_name
        return True, host, label, color, contig_id or tip_name

    # ICTV reference
    if tip_name.startswith("ICTV|"):
        parts = tip_name.split("|")
        sp = parts[2].replace("_", " ") if len(parts) > 2 else tip_name
        acc = parts[1] if len(parts) > 1 else ""
        label = f"{sp} ({acc})" if acc else sp
    else:
        label = tip_name.replace("_", " ")
    if len(label) > 48:
        label = label[:45] + "..."
    return False, "ICTV reference", label, ICTV_COLOR, ""


def orf_colour(product: str) -> str:
    low = (product or "").lower()
    for pat, col, _ in ORF_COLOURS:
        if re.search(pat, low):
            return col
    return ORF_DEFAULT_COLOUR


# ---------------------------------------------------------------------------
# SVG panel renderer for a phylogenetic tree (+ optional ORF track)
# ---------------------------------------------------------------------------
def render_tree_panel(
    tree: TreeNode,
    meta: Dict[str, Dict[str, str]],
    manifest_hosts: Dict[str, str],
    *,
    x0: float,
    y0: float,
    tree_w: float,
    label_w: float,
    row_h: float,
    panel_letter: str,
    panel_title: str,
    panel_subtitle: str,
    orf_map: Optional[Dict[str, Tuple[int, List[Dict[str, str]]]]] = None,
    orf_x0: float = 0.0,
    orf_w: float = 0.0,
    max_genome_bp: int = 13000,
) -> Tuple[List[str], float, List[str]]:
    """Render one phylogenetic tree panel and return (svg_lines, bottom_y, hosts_seen)."""
    leaves = _all_leaves(tree)
    n_leaves = len(leaves)
    svg: List[str] = []
    hosts_seen: List[str] = []

    # Panel header
    svg.append(
        f'<text x="{x0:.1f}" y="{y0 + 18:.1f}" font-size="16" font-weight="700" '
        f'fill="#0F172A">{escape(panel_letter)}</text>'
    )
    svg.append(
        f'<text x="{x0 + 22:.1f}" y="{y0 + 18:.1f}" font-size="13.5" font-weight="700" '
        f'fill="#0F172A">{escape(panel_title)}</text>'
    )
    if panel_subtitle:
        svg.append(
            f'<text x="{x0 + 22:.1f}" y="{y0 + 34:.1f}" font-size="10.5" '
            f'fill="#475569">{escape(panel_subtitle)}</text>'
        )

    top_y = y0 + 50.0

    # Assign leaf Y positions and cumulative X distances
    leaf_y: Dict[int, float] = {}
    for idx, lf in enumerate(leaves):
        leaf_y[id(lf)] = top_y + idx * row_h

    node_x_raw: Dict[int, float] = {}
    node_y: Dict[int, float] = {}

    def compute_coords(node: TreeNode, depth: float) -> float:
        x_val = depth + node.length
        node_x_raw[id(node)] = x_val
        if node.is_leaf():
            y_val = leaf_y[id(node)]
            node_y[id(node)] = y_val
            return y_val
        ys = [compute_coords(c, x_val) for c in node.children]
        y_val = (min(ys) + max(ys)) / 2.0
        node_y[id(node)] = y_val
        return y_val

    compute_coords(tree, 0.0)
    max_x = max(node_x_raw.values()) if node_x_raw else 1.0
    if max_x <= 1e-9:
        max_x = 1.0

    def sx(raw_x: float) -> float:
        return x0 + 10.0 + (raw_x / max_x) * (tree_w - 20.0)

    # Subtle alternating row bands & ORF panel header if requested
    if orf_map is not None and orf_w > 0:
        svg.append(
            f'<text x="{orf_x0:.1f}" y="{y0 + 18:.1f}" font-size="16" font-weight="700" '
            f'fill="#0F172A">B</text>'
        )
        svg.append(
            f'<text x="{orf_x0 + 22:.1f}" y="{y0 + 18:.1f}" font-size="13.5" font-weight="700" '
            f'fill="#0F172A">Genome &amp; ORF Architecture</text>'
        )
        svg.append(
            f'<text x="{orf_x0 + 22:.1f}" y="{y0 + 34:.1f}" font-size="10.5" '
            f'fill="#475569">To-scale predicted ORFs aligned to tree tips (0–{max_genome_bp/1000:.0f} kb)</text>'
        )
        # Vertical kb gridlines
        step_kb = 2 if max_genome_bp <= 16000 else 5
        for kb in range(0, int(max_genome_bp // 1000) + 1, step_kb):
            gx = orf_x0 + (kb * 1000.0 / max_genome_bp) * orf_w
            svg.append(
                f'<line x1="{gx:.1f}" y1="{top_y - 8:.1f}" '
                f'x2="{gx:.1f}" y2="{top_y + max(0, n_leaves - 1) * row_h + 10:.1f}" '
                f'stroke="#E2E8F0" stroke-width="0.8" stroke-dasharray="2,2"/>'
            )
            svg.append(
                f'<text x="{gx:.1f}" y="{top_y + max(0, n_leaves - 1) * row_h + 24:.1f}" '
                f'text-anchor="middle" font-size="9" fill="#64748B">{kb} kb</text>'
            )

    # Draw branches
    def draw_branches(node: TreeNode) -> None:
        nx = sx(node_x_raw[id(node)])
        ny = node_y[id(node)]
        if not node.is_leaf():
            cys = [node_y[id(c)] for c in node.children]
            svg.append(
                f'<line x1="{nx:.2f}" y1="{min(cys):.2f}" x2="{nx:.2f}" y2="{max(cys):.2f}" '
                f'stroke="#334155" stroke-width="1.35" stroke-linecap="round"/>'
            )
            for c in node.children:
                cx = sx(node_x_raw[id(c)])
                cy = node_y[id(c)]
                svg.append(
                    f'<line x1="{nx:.2f}" y1="{cy:.2f}" x2="{cx:.2f}" y2="{cy:.2f}" '
                    f'stroke="#334155" stroke-width="1.35" stroke-linecap="round"/>'
                )
                draw_branches(c)
            # Support marker (UFB >= 80)
            ufb = node.support_ufb
            if ufb is not None and ufb >= 80.0 and node.parent is not None:
                dot_col = "#0F172A" if ufb >= 95.0 else "#64748B"
                svg.append(
                    f'<circle cx="{nx:.2f}" cy="{ny:.2f}" r="2.6" fill="{dot_col}"/>'
                )
                if ufb >= 90.0:
                    sup_txt = (
                        f"{node.support_alrt:.0f}/{ufb:.0f}"
                        if node.support_alrt is not None
                        else f"{ufb:.0f}"
                    )
                    svg.append(
                        f'<text x="{nx - 4:.1f}" y="{ny - 3.5:.1f}" text-anchor="end" '
                        f'font-size="8" fill="#475569">{escape(sup_txt)}</text>'
                    )

    draw_branches(tree)

    # Draw tip dots, labels, and optional aligned ORF tracks
    tip_align_x = x0 + tree_w + 4.0
    for lf in leaves:
        lx = sx(node_x_raw[id(lf)])
        ly = node_y[id(lf)]
        is_cand, host, label, color, contig_key = classify_tip(lf.name, meta, manifest_hosts)
        if is_cand and host and host not in hosts_seen:
            hosts_seen.append(host)

        # Dotted guide from tip to aligned label column
        if tip_align_x - lx > 4.0:
            svg.append(
                f'<line x1="{lx + 3:.2f}" y1="{ly:.2f}" x2="{tip_align_x - 3:.2f}" y2="{ly:.2f}" '
                f'stroke="#CBD5E1" stroke-width="0.7" stroke-dasharray="2,2"/>'
            )

        if is_cand:
            svg.append(
                f'<circle cx="{lx:.2f}" cy="{ly:.2f}" r="3.8" fill="{color}" '
                f'stroke="#FFFFFF" stroke-width="1.0"/>'
            )
            svg.append(
                f'<text x="{tip_align_x:.1f}" y="{ly + 3.5:.1f}" font-size="10" '
                f'font-weight="700" fill="{color}">{escape(label)}</text>'
            )
        else:
            svg.append(
                f'<circle cx="{lx:.2f}" cy="{ly:.2f}" r="2.2" fill="{ICTV_COLOR}"/>'
            )
            svg.append(
                f'<text x="{tip_align_x:.1f}" y="{ly + 3.5:.1f}" font-size="9.5" '
                f'font-style="italic" fill="#475569">{escape(label)}</text>'
            )

        # Optional aligned ORF track on the right
        if orf_map is not None and orf_w > 0 and is_cand:
            entry = orf_map.get(contig_key)
            if entry is None:
                # Fallback substring match on contig_key
                for k, v in orf_map.items():
                    if contig_key and (contig_key in k or k in contig_key):
                        entry = v
                        break
            if entry is not None:
                clen, orfs = entry
                bar_len = min(orf_w, (clen / max_genome_bp) * orf_w) if clen > 0 else 0.0
                if bar_len > 0:
                    svg.append(
                        f'<line x1="{orf_x0:.1f}" y1="{ly:.2f}" '
                        f'x2="{orf_x0 + bar_len:.1f}" y2="{ly:.2f}" '
                        f'stroke="#475569" stroke-width="1.6"/>'
                    )
                for orf in orfs:
                    st = max(0, int(float(orf.get("ORF_start") or 0)))
                    en = min(max_genome_bp, int(float(orf.get("ORF_end") or 0)))
                    if en <= st:
                        continue
                    ox = orf_x0 + (st / max_genome_bp) * orf_w
                    ow = max(3.0, ((en - st) / max_genome_bp) * orf_w)
                    col = orf_colour(orf.get("product", ""))
                    svg.append(
                        f'<rect x="{ox:.1f}" y="{ly - 5.0:.1f}" width="{ow:.1f}" height="10.0" '
                        f'rx="2" fill="{col}" stroke="#FFFFFF" stroke-width="0.6"/>'
                    )
                if clen > 0:
                    svg.append(
                        f'<text x="{orf_x0 + bar_len + 5:.1f}" y="{ly + 3.2:.1f}" '
                        f'font-size="8.5" fill="#475569">{clen/1000:.1f} kb</text>'
                    )

    # Scale bar
    scale_y = top_y + max(0, n_leaves - 1) * row_h + 22.0
    # Choose nice scale bar length (~15-25% of max_x)
    raw_step = max_x * 0.2
    mag = 10 ** math.floor(math.log10(max(raw_step, 1e-6)))
    nice = round(raw_step / mag) * mag
    if nice <= 0:
        nice = mag
    sb_w = (nice / max_x) * (tree_w - 20.0)
    sb_x = x0 + 12.0
    svg.append(
        f'<line x1="{sb_x:.1f}" y1="{scale_y:.1f}" x2="{sb_x + sb_w:.1f}" y2="{scale_y:.1f}" '
        f'stroke="#0F172A" stroke-width="1.8"/>'
    )
    svg.append(
        f'<text x="{sb_x + sb_w / 2:.1f}" y="{scale_y + 13:.1f}" text-anchor="middle" '
        f'font-size="9" fill="#334155">{nice:g} subs/site</text>'
    )

    return svg, scale_y + 24.0, hosts_seen


def build_legend(x0: float, y0: float, include_orfs: bool) -> Tuple[List[str], float]:
    svg: List[str] = []
    cx = x0
    svg.append(
        f'<text x="{cx:.1f}" y="{y0 + 11:.1f}" font-size="9.5" font-weight="700" '
        f'fill="#334155">Host / Continent:</text>'
    )
    cx += 96.0
    for sp, (col, _reg, short) in HOST_STYLE.items():
        svg.append(f'<circle cx="{cx:.1f}" cy="{y0 + 7.5:.1f}" r="4.0" fill="{col}"/>')
        svg.append(
            f'<text x="{cx + 8:.1f}" y="{y0 + 11:.1f}" font-size="9.5" fill="#1E293B">{escape(short)}</text>'
        )
        cx += 118.0
    svg.append(f'<circle cx="{cx:.1f}" cy="{y0 + 7.5:.1f}" r="3.0" fill="{ICTV_COLOR}"/>')
    svg.append(
        f'<text x="{cx + 8:.1f}" y="{y0 + 11:.1f}" font-size="9.5" font-style="italic" '
        f'fill="#475569">ICTV reference</text>'
    )
    y_next = y0 + 22.0
    if include_orfs:
        ox = x0
        svg.append(
            f'<text x="{ox:.1f}" y="{y_next + 11:.1f}" font-size="9.5" font-weight="700" '
            f'fill="#334155">ORF Annotation:</text>'
        )
        ox += 96.0
        for _pat, col, lbl in ORF_COLOURS:
            svg.append(
                f'<rect x="{ox:.1f}" y="{y_next + 3:.1f}" width="11" height="9" rx="2" fill="{col}"/>'
            )
            svg.append(
                f'<text x="{ox + 15:.1f}" y="{y_next + 11:.1f}" font-size="9.5" fill="#1E293B">{escape(lbl)}</text>'
            )
            ox += 118.0
        svg.append(
            f'<rect x="{ox:.1f}" y="{y_next + 3:.1f}" width="11" height="9" rx="2" fill="{ORF_DEFAULT_COLOUR}"/>'
        )
        svg.append(
            f'<text x="{ox + 15:.1f}" y="{y_next + 11:.1f}" font-size="9.5" fill="#475569">Hypothetical / other</text>'
        )
        y_next += 22.0
    return svg, y_next


def _prepare_tree(
    tree_path: Path,
    meta: Dict[str, Dict[str, str]],
    manifest_hosts: Dict[str, str],
    exclude_hosts: Sequence[str],
) -> TreeNode:
    tree = parse_newick(tree_path.read_text(encoding="utf-8"))
    if exclude_hosts:
        ex_low = [h.lower() for h in exclude_hosts if h]
        drop = set()
        for lf in _all_leaves(tree):
            is_cand, host, _lbl, _col, _ck = classify_tip(lf.name, meta, manifest_hosts)
            if is_cand and any(ex in host.lower() for ex in ex_low):
                drop.add(lf.name)
        if drop:
            tree = _prune_leaves(tree, drop)
    tree = midpoint_root(tree)
    ladderize(tree)
    return tree


def write_svg_atomic(out_path: Path, width: float, height: float, body: List[str]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.0f} {height:.0f}" '
        f'width="{width:.0f}" height="{height:.0f}" '
        f'font-family="Inter, Helvetica Neue, Arial, sans-serif">',
        f'<rect width="{width:.0f}" height="{height:.0f}" fill="#FFFFFF"/>',
    ]
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    tmp.write_text("\n".join(header + body + ["</svg>\n"]), encoding="utf-8")
    os.replace(tmp, out_path)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", required=True, choices=["fig4", "fig5", "fig6"])
    ap.add_argument("--tree", type=Path, help="Primary Newick treefile (fig4 / fig5).")
    ap.add_argument("--meta", type=Path, help="Primary taxa_metadata.tsv (fig4 / fig5).")
    ap.add_argument("--genome-map", type=Path, help="Optional genome_map_data.tsv (fig4 / fig5).")
    ap.add_argument(
        "--panel",
        action="append",
        default=[],
        help="For fig6: repeatable 'FamilyLabel|tree.newick|taxa_metadata.tsv' spec.",
    )
    ap.add_argument("--manifest", type=Path, help="Combined manifest CSV for run->host lookup.")
    ap.add_argument(
        "--exclude-hosts",
        default="Locusta migratoria",
        help="Comma-separated host names to exclude from candidate tips.",
    )
    ap.add_argument("--out-svg", required=True, type=Path)
    args = ap.parse_args(argv)

    manifest_hosts = load_manifest_hosts(args.manifest)
    exclude_hosts = [h.strip() for h in args.exclude_hosts.split(",") if h.strip()]

    body: List[str] = []

    if args.mode == "fig4":
        if not args.tree or not args.tree.is_file():
            raise SystemExit(f"Missing --tree for fig4: {args.tree}")
        meta = load_taxa_metadata(args.meta)
        tree = _prepare_tree(args.tree, meta, manifest_hosts, exclude_hosts)
        orf_map = load_genome_map_data(args.genome_map)
        max_bp = max((v[0] for v in orf_map.values()), default=13000)
        max_bp = max(13000, int(math.ceil(max_bp / 1000.0) * 1000))

        row_h = 18.0
        width = 1180.0
        panel_svg, bottom_y, _ = render_tree_panel(
            tree,
            meta,
            manifest_hosts,
            x0=24.0,
            y0=18.0,
            tree_w=380.0,
            label_w=270.0,
            row_h=row_h,
            panel_letter="A",
            panel_title="Pan-Holarctic Ips Polycipiviridae — RdRp Amino-Acid Phylogeny",
            panel_subtitle="Midpoint-rooted IQ-TREE2 ML tree (SH-aLRT / UFB ≥ 90% shown; dark dots UFB ≥ 95%)",
            orf_map=orf_map,
            orf_x0=700.0,
            orf_w=430.0,
            max_genome_bp=max_bp,
        )
        body.extend(panel_svg)
        leg_svg, final_y = build_legend(24.0, bottom_y + 6.0, include_orfs=True)
        body.extend(leg_svg)
        write_svg_atomic(args.out_svg, width, final_y + 18.0, body)
        LOG.info("Wrote Fig. 4 Polycipiviridae SVG (%d leaves) -> %s", len(_all_leaves(tree)), args.out_svg)
        return 0

    if args.mode == "fig5":
        if not args.tree or not args.tree.is_file():
            raise SystemExit(f"Missing --tree for fig5: {args.tree}")
        meta = load_taxa_metadata(args.meta)
        tree = _prepare_tree(args.tree, meta, manifest_hosts, exclude_hosts)
        orf_map = load_genome_map_data(args.genome_map) if args.genome_map else None
        has_orfs = bool(orf_map)
        width = 1140.0 if has_orfs else 860.0
        panel_svg, bottom_y, _ = render_tree_panel(
            tree,
            meta,
            manifest_hosts,
            x0=24.0,
            y0=18.0,
            tree_w=460.0 if not has_orfs else 380.0,
            label_w=310.0,
            row_h=17.5,
            panel_letter="A",
            panel_title="Closed-Genome Conifer-Pest Iflaviridae — RdRp Amino-Acid Phylogeny",
            panel_subtitle="Midpoint-rooted IQ-TREE2 ML tree across North American, European, and East Asian hosts",
            orf_map=orf_map if has_orfs else None,
            orf_x0=710.0,
            orf_w=380.0,
            max_genome_bp=11000,
        )
        body.extend(panel_svg)
        leg_svg, final_y = build_legend(24.0, bottom_y + 6.0, include_orfs=has_orfs)
        body.extend(leg_svg)
        write_svg_atomic(args.out_svg, width, final_y + 18.0, body)
        LOG.info("Wrote Fig. 5 Iflaviridae SVG (%d leaves) -> %s", len(_all_leaves(tree)), args.out_svg)
        return 0

    # fig6: multi-panel composite (Amarillovirales, Phenuiviridae, Aliusviridae)
    if not args.panel:
        raise SystemExit("fig6 requires at least one --panel 'Title|tree.newick|taxa_metadata.tsv'")
    width = 920.0
    cur_y = 18.0
    letters = ["A", "B", "C", "D", "E"]
    for idx, spec in enumerate(args.panel):
        parts = spec.split("|")
        fam_title = parts[0].strip()
        t_path = Path(parts[1].strip())
        m_path = Path(parts[2].strip()) if len(parts) > 2 and parts[2].strip() else None
        if not t_path.is_file():
            LOG.warning("Skipping missing panel tree: %s", t_path)
            continue
        meta = load_taxa_metadata(m_path)
        tree = _prepare_tree(t_path, meta, manifest_hosts, exclude_hosts)
        p_svg, cur_y, _ = render_tree_panel(
            tree,
            meta,
            manifest_hosts,
            x0=24.0,
            y0=cur_y,
            tree_w=490.0,
            label_w=340.0,
            row_h=17.0,
            panel_letter=letters[idx % len(letters)],
            panel_title=f"{fam_title} — RdRp Amino-Acid Phylogeny",
            panel_subtitle="Midpoint-rooted IQ-TREE2 ML tree (SH-aLRT / UFB ≥ 90% shown)",
        )
        body.extend(p_svg)
        cur_y += 14.0

    leg_svg, final_y = build_legend(24.0, cur_y + 4.0, include_orfs=False)
    body.extend(leg_svg)
    write_svg_atomic(args.out_svg, width, final_y + 18.0, body)
    LOG.info("Wrote Fig. 6 multi-panel SVG (%d panels) -> %s", len(args.panel), args.out_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
