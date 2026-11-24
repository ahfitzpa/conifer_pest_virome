#!/usr/bin/env python3
"""Apply family-specific ICTV demarcation thresholds to candidate viruses.

Thresholds come only from configs/phylogenetics.yaml:demarcation.thresholds.
If a family has no configured threshold, or the criterion marker differs from
the marker actually evaluated (RdRp/L alignment), the row is labelled
explicitly and NO species/genus call is made.
"""
from __future__ import annotations
import argparse, csv, sys

import yaml


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def marker_matches(criterion, evaluated, aliases=None):
    """How the criterion's marker relates to the one the alignment measures.

    Returns "direct" for a substring match either way, so "rdrp" matches
    "RdRp" and "RdRP amino acid"; "alias" when the config declares the two
    interchangeable under demarcation.marker_aliases; None otherwise. A
    criterion defined on a different protein (capsid, ORF5, concatenated core
    genes) does NOT match an RdRp alignment and must not be applied to it.

    An alias is an approximation the caller has chosen to accept, so every row
    that used one says so in its note.
    """
    if not criterion or not evaluated:
        return None
    c, e = criterion.lower(), evaluated.lower()
    if c in e or e in c:
        return "direct"
    for alias in (aliases or {}).get(e, []):
        if alias.lower() == c:
            return "alias"
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--pairwise", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cfg = yaml.safe_load(open(a.config, encoding="utf-8"))["phylogenetics"]
    rule = cfg["demarcation"]["thresholds"].get(a.family, {})
    criterion_marker = rule.get("species_marker")
    species_max = rule.get("species_max_pct")
    genus_marker = rule.get("genus_marker")
    genus_max = rule.get("genus_max_pct")
    evaluated_marker = cfg.get("marker", "rdrp")  # the phylo alignment marker
    aliases = cfg["demarcation"].get("marker_aliases") or {}

    rows = [r for r in read_tsv(a.pairwise)]
    best = {}
    for r in rows:
        if r.get("passes_50pct_coverage") != "1":
            continue
        q = r["query"]
        if q not in best or float(r["identity_pct"]) > float(best[q]["identity_pct"]):
            best[q] = r

    fields = ["family", "candidate", "closest_ictv", "identity_pct", "query_coverage_pct",
              "alignment_marker", "criterion_marker", "species_threshold_pct",
              "below_species_threshold", "genus_threshold_pct", "below_genus_threshold",
              "call", "status", "note"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for q, r in sorted(best.items()):
            base = {"family": a.family, "candidate": q,
                    "closest_ictv": r["reference"],
                    "identity_pct": r["identity_pct"],
                    "query_coverage_pct": r["query_coverage_pct"],
                    "alignment_marker": evaluated_marker,
                    "criterion_marker": criterion_marker or "NA",
                    "species_threshold_pct": species_max if species_max is not None else "NA",
                    "below_species_threshold": "NA",
                    "genus_threshold_pct": genus_max if genus_max is not None else "NA",
                    "below_genus_threshold": "NA",
                    "call": "withheld", "status": "ok", "note": rule.get("note", "")}
            # Species and genus are judged independently. Partitiviridae and
            # Spinareoviridae publish a numeric GENUS criterion and no species
            # one; evaluating genus only when a species threshold existed threw
            # those away and withheld every row.
            identity = float(r["identity_pct"])
            reasons, below, below_g = [], None, None

            s_match = marker_matches(criterion_marker, evaluated_marker, aliases)
            g_match = marker_matches(genus_marker, evaluated_marker, aliases)

            if species_max is None:
                reasons.append("species criterion not configured")
            elif s_match is None:
                reasons.append(f"species criterion uses {criterion_marker} but "
                               f"alignment is {evaluated_marker}")
            else:
                below = identity < float(species_max)
                base["below_species_threshold"] = int(below)
                if s_match == "alias":
                    reasons.append(f"species threshold defined on {criterion_marker}, "
                                   f"applied to the {evaluated_marker} alignment by "
                                   "configured alias (approximation)")

            if genus_max is None:
                reasons.append("genus criterion not configured")
            elif g_match is None:
                reasons.append(f"genus criterion uses {genus_marker} but "
                               f"alignment is {evaluated_marker}")
            else:
                below_g = identity < float(genus_max)
                base["below_genus_threshold"] = int(below_g)
                if g_match == "alias":
                    reasons.append(f"genus threshold defined on {genus_marker}, "
                                   f"applied to the {evaluated_marker} alignment by "
                                   "configured alias (approximation)")

            if below_g:
                base["call"] = "novel genus"
            elif below is not None:
                base["call"] = "novel species" if below else "existing species"
            elif below_g is False:
                base["call"] = "existing genus"
            mismatch = ((species_max is not None and s_match is None)
                        or (genus_max is not None and g_match is None))
            if base["call"] == "withheld":
                base["status"] = ("criterion_marker_mismatch" if mismatch
                                  else "criterion_not_configured")
            elif "alias" in (s_match, g_match):
                # The call stands, but it rests on a marker equivalence the
                # config asserts rather than one ICTV states. Filterable.
                base["status"] = "ok_marker_alias"
            if reasons:
                base["note"] = (base["note"] + " | " + "; ".join(reasons)).strip(" |")
            w.writerow(base)
    print(f"{a.family}: {len(best)} candidates evaluated -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
