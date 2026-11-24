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


def marker_matches(criterion, evaluated):
    """True when the criterion's marker is the one the alignment measures.

    Substring either way, so "rdrp" matches "RdRp" and "RdRP amino acid". A
    criterion defined on a different protein (capsid, ORF5, the full-length L,
    concatenated core genes) does NOT match an RdRp alignment and must not be
    applied to it.
    """
    if not criterion or not evaluated:
        return False
    c, e = criterion.lower(), evaluated.lower()
    return c in e or e in c


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

            if species_max is None:
                reasons.append("species criterion not configured")
            elif not marker_matches(criterion_marker, evaluated_marker):
                reasons.append(f"species criterion uses {criterion_marker} but "
                               f"alignment is {evaluated_marker}")
            else:
                below = identity < float(species_max)
                base["below_species_threshold"] = int(below)

            if genus_max is None:
                reasons.append("genus criterion not configured")
            elif not marker_matches(genus_marker, evaluated_marker):
                reasons.append(f"genus criterion uses {genus_marker} but "
                               f"alignment is {evaluated_marker}")
            else:
                below_g = identity < float(genus_max)
                base["below_genus_threshold"] = int(below_g)

            if below_g:
                base["call"] = "novel genus"
            elif below is not None:
                base["call"] = "novel species" if below else "existing species"
            elif below_g is False:
                base["call"] = "existing genus"
            if base["call"] == "withheld":
                base["status"] = ("criterion_marker_mismatch"
                                  if any("alignment is" in x for x in reasons)
                                  else "criterion_not_configured")
            if reasons:
                base["note"] = (base["note"] + " | " + "; ".join(reasons)).strip(" |")
            w.writerow(base)
    print(f"{a.family}: {len(best)} candidates evaluated -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
