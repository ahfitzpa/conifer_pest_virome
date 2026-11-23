#!/usr/bin/env python3
"""Reconcile the three taxonomy callers' family names against ICTV lineage.

The old rule was a plain string comparison: any two different names were a
`conflict`, and conflicts are excluded from every downstream analysis. That
discarded contigs for two reasons that are not disagreements at all:

  * mixed rank. geNomad reports a family, RdRpCATCH's profile sets are named
    after whatever rank the profile covers, so Flaviviridae vs Amarillovirales
    and Orthomyxoviridae vs Articulavirales are the same answer at two ranks.
  * a superseded name. Reoviridae was split into Sedoreoviridae and
    Spinareoviridae, so Spinareoviridae vs Reoviridae is one taxon, not two.

Where the names really are sister families (Picornaviridae vs Iflaviridae,
Rhabdoviridae vs Xinmoviridae), RdRpCATCH wins for RNA viruses: it assigns
from the RdRp itself, which is the family-level marker, where geNomad votes
over whole-genome protein content. The rejected calls are recorded in
`family_disagreement` so nothing is silently overwritten.

Statuses returned, in `family_agreement`:

  single_tool    one tool resolved a family
  agree          every tool that resolved one gave the same name
  agree_rank     same lineage at different ranks, or a superseded name;
                 consensus_family is the most specific valid name
  rdrp_resolved  sister-family disagreement settled by RdRpCATCH
  conflict       names disagree and RdRpCATCH has no call to arbitrate
  unresolved     no tool resolved a family
"""
from __future__ import annotations

import csv
import os

RANK_DEPTH = {"phylum": 0, "class": 1, "order": 2, "family": 3}

# common/ sits one level up in conifer (scripts/) and two in Serratus
# (phase2/scripts/), so the table is searched for rather than assumed.
# FAMILY_LINEAGE_TSV overrides, for a layout neither pattern covers.
_HERE = os.path.dirname(os.path.abspath(__file__))
_CANDIDATES = [
    os.environ.get("FAMILY_LINEAGE_TSV", ""),
    os.path.join(_HERE, os.pardir, "common", "taxonomy", "family_lineage.tsv"),
    os.path.join(_HERE, os.pardir, os.pardir, "common", "taxonomy", "family_lineage.tsv"),
    os.path.join(_HERE, os.pardir, "taxonomy", "family_lineage.tsv"),
]
LINEAGE_TSV = next((os.path.normpath(p) for p in _CANDIDATES
                    if p and os.path.exists(p)),
                   os.path.normpath(_CANDIDATES[1]))

# Taxa the callers still emit that ICTV has renamed or split. A name maps to
# the successor(s) that replaced it; a call is compatible with any of them.
# Extend this rather than special-casing in the comparison.
SUPERSEDED = {
    "reoviridae": {"sedoreoviridae", "spinareoviridae"},
    "bunyaviridae": {"peribunyaviridae", "phenuiviridae", "hantaviridae",
                     "nairoviridae"},
    "luteoviridae": {"solemoviridae", "tombusviridae"},
}

_LINEAGE: dict[str, dict[str, str]] | None = None


def load_lineage(path: str = LINEAGE_TSV) -> dict[str, dict[str, str]]:
    """taxon (lowercased) -> {rank, phylum, class, order}; empty if absent."""
    global _LINEAGE
    if _LINEAGE is not None:
        return _LINEAGE
    table: dict[str, dict[str, str]] = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                name = (row.get("name") or "").strip()
                if name:
                    table[name.lower()] = {
                        "rank": (row.get("rank") or "").strip().lower(),
                        "phylum": (row.get("phylum") or "").strip(),
                        "class": (row.get("class") or "").strip(),
                        "order": (row.get("order") or "").strip(),
                    }
    _LINEAGE = table
    return table


def _depth(name: str, table) -> int:
    return RANK_DEPTH.get(table.get(name.lower(), {}).get("rank", ""), -1)


def _is_ancestor(ancestor: str, descendant: str, table) -> bool:
    """True when `ancestor` is one of `descendant`'s higher-rank lineage names."""
    rec = table.get(descendant.lower())
    if not rec:
        return False
    anc_rank = table.get(ancestor.lower(), {}).get("rank", "")
    if anc_rank not in ("phylum", "class", "order"):
        return False
    return rec.get(anc_rank, "").lower() == ancestor.lower()


def _compatible(a: str, b: str, table) -> bool:
    """Same taxon at a different rank, or a superseded name and its successor."""
    if a.lower() == b.lower():
        return True
    if _is_ancestor(a, b, table) or _is_ancestor(b, a, table):
        return True
    for old, new in ((a, b), (b, a)):
        successors = SUPERSEDED.get(old.lower())
        if successors and (new.lower() in successors
                           or any(_is_ancestor(new, s, table) for s in successors)):
            return True
    return False


def reconcile(fams: dict[str, str], tool_priority, rdrp_tool: str = "rdrpcatch",
              lineage_path: str = LINEAGE_TSV):
    """Return (consensus_family, family_agreement, family_disagreement).

    `fams` maps tool name -> family call; blank calls are ignored by the caller.
    """
    table = load_lineage(lineage_path)
    calls = {t: f for t, f in fams.items() if isinstance(f, str) and f.strip()}
    if not calls:
        return "", "unresolved", ""

    names = list(dict.fromkeys(calls.values()))
    if len(names) == 1:
        return names[0], ("agree" if len(calls) > 1 else "single_tool"), ""

    note = ";".join(f"{t}={calls[t]}" for t in tool_priority if t in calls)

    # 1: one lineage at mixed ranks? take the most specific known name
    best = max(names, key=lambda n: (_depth(n, table), -names.index(n)))
    if all(_compatible(best, other, table) for other in names if other != best):
        if _depth(best, table) < 0:
            # nothing in the table resolved; fall through to arbitration
            pass
        else:
            return best, "agree_rank", note

    # 2: genuine disagreement -- RdRpCATCH decides for RNA viruses, because its
    # call comes from the RdRp, the family-level marker
    rdrp_call = calls.get(rdrp_tool, "").strip()
    if rdrp_call:
        return rdrp_call, "rdrp_resolved", note

    # 3: no RdRp call to arbitrate: keep the old majority/priority winner
    from collections import Counter
    counts = Counter(calls.values())
    top = max(counts.values())
    leaders = {f for f, c in counts.items() if c == top}
    winner = next(calls[t] for t in tool_priority if calls.get(t) in leaders)
    return winner, "conflict", note
