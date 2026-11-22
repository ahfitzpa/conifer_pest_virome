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
            if species_max is None:
                base["status"] = "criterion_not_configured"
                base["note"] = (base["note"] + " | species criterion not configured; "
                                "no call made").strip(" |")
            elif criterion_marker and evaluated_marker and \
                    criterion_marker.lower() not in evaluated_marker.lower() and \
                    evaluated_marker.lower() not in criterion_marker.lower():
                base["status"] = "criterion_marker_mismatch"
                base["note"] = (base["note"] + f" | criterion uses {criterion_marker} but "
                                f"alignment is {evaluated_marker}; no call made").strip(" |")
            else:
                below = float(r["identity_pct"]) < float(species_max)
                base["below_species_threshold"] = int(below)
                if genus_max is not None and genus_marker:
                    below_g = float(r["identity_pct"]) < float(genus_max)
                    base["below_genus_threshold"] = int(below_g)
                    base["call"] = ("novel genus" if below_g else
                                    "novel species" if below else "existing species")
                else:
                    base["call"] = "novel species" if below else "existing species"
            w.writerow(base)
    print(f"{a.family}: {len(best)} candidates evaluated -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
