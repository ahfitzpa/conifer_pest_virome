#!/usr/bin/env python3
"""Decide which ICTV families are worth analysing, and say why for each one.

`ictv_prepare.py` writes a bundle for every family with at least one candidate
contig -- 122 in conifer, 318 in Serratus -- and `phylogenetics.families: []`
means "all of them". Most of those families are not findings: in conifer, 27
rest on a single contig from a single run, and 70 come from one BioProject.
Running ML phylogenetics, genome maps, recombination and demarcation over all of
them is mostly compute spent on noise.

The criteria below are evidence-based and each is recorded per family, so the
decision is auditable rather than a threshold someone has to trust:

  in_ictv        the family name appears in the ICTV VMR. `_unresolved` and
                 anything misspelled or renamed fails here.
  n_contigs      at least `--min-contigs` contigs assigned to it.
  n_runs         seen in at least `--min-runs` sequencing runs. One run cannot
                 distinguish a virus from a one-off misassignment.
  n_bioprojects  seen in at least `--min-bioprojects` BioProjects. This is the
                 criterion phase1 already applies to sOTUs
                 (`min_bioprojects_per_sotu`), and it is what separates a real
                 signal from one study's systematic artefact.
  rdrp           an RdRp domain is required ONLY where the VMR says the family
                 is an RNA virus and not reverse-transcribing. Demanding RdRp of
                 DNA families would wrongly discard legitimate ones -- in these
                 datasets Baculoviridae, Kyanoviridae and other dsDNA families
                 carry no RdRp by definition.

Target families named with `--target-family` always pass: the scenario's family
of interest must never be filtered out by a generic rule.

Tiers: A passes everything; B fails only replication (one run or one
BioProject) while having the contigs and the expected replicase, so it is a
candidate worth a look but not a conclusion; C is excluded.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

NULL = {"", "nan", "na", "none", "null", "-"}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def clean(v) -> str:
    v = str(v or "").strip()
    return "" if v.lower() in NULL else v


def as_int(v) -> int:
    try:
        return int(float(clean(v)))
    except ValueError:
        return 0


def load_vmr(path: str) -> dict[str, str]:
    """family -> genome composition, from the ICTV VMR (xlsx or csv)."""
    if not path or not os.path.exists(path):
        return {}
    rows = []
    if path.lower().endswith((".xlsx", ".xls")):
        try:
            import pandas as pd
        except ImportError:
            log("[families] pandas unavailable; pass the VMR as csv instead")
            return {}
        sheet = None
        xl = pd.ExcelFile(path)
        for name in xl.sheet_names:
            if "vmr" in name.lower() or "msl" in name.lower():
                sheet = name
                break
        df = pd.read_excel(path, sheet_name=sheet or 0)
        rows = df.to_dict("records")
    else:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))

    counts: dict[str, dict[str, int]] = {}
    for r in rows:
        fam = clean(r.get("Family") or r.get("family"))
        comp = clean(r.get("Genome") or r.get("genome"))
        if fam and comp:
            counts.setdefault(fam, {}).setdefault(comp, 0)
            counts[fam][comp] += 1
    # 1: a family's composition is its commonest entry (16 families list more
    # than one across their isolates)
    return {fam: max(c, key=c.get) for fam, c in counts.items()}


def expects_rdrp(composition: str) -> bool:
    # 2: RNA viruses encode an RdRp; reverse-transcribing ones use RT instead
    c = composition.upper()
    return "RNA" in c and "RT" not in c


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--families-tsv", required=True, help="ictv_prepare output")
    p.add_argument("--vmr", help="ICTV VMR xlsx/csv; without it in_ictv is skipped")
    p.add_argument("--out-tsv", required=True, help="decisions per family")
    p.add_argument("--out-families", help="newline-separated passing families")
    p.add_argument("--out-yaml", help="YAML fragment for phylogenetics.families")
    p.add_argument("--min-contigs", type=int, default=2)
    p.add_argument("--min-runs", type=int, default=2)
    p.add_argument("--min-bioprojects", type=int, default=2)
    p.add_argument("--require-rdrp-for-rna", action="store_true", default=True)
    p.add_argument("--no-require-rdrp-for-rna", dest="require_rdrp_for_rna",
                   action="store_false")
    p.add_argument("--allow-unlisted", action="store_true",
                   help="keep families absent from the VMR")
    p.add_argument("--target-family", action="append", default=[],
                   help="always keep (repeatable)")
    p.add_argument("--tiers", default="A",
                   help="comma-separated tiers to pass, e.g. A or A,B")
    args = p.parse_args()

    vmr = load_vmr(args.vmr) if args.vmr else {}
    if args.vmr:
        log(f"[families] VMR: {len(vmr)} families, "
            f"{sum(1 for c in vmr.values() if expects_rdrp(c))} RNA")
    keep_tiers = {t.strip().upper() for t in args.tiers.split(",") if t.strip()}
    targets = set(args.target_family)

    with open(args.families_tsv, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    log(f"[families] {len(rows)} family bundle(s) in {args.families_tsv}")

    out, counts = [], {"A": 0, "B": 0, "C": 0}
    for r in rows:
        fam = clean(r.get("assigned_family"))
        n_contigs, n_runs = as_int(r.get("n_contigs")), as_int(r.get("n_runs"))
        n_bp = as_int(r.get("n_bioprojects"))
        n_rdrp = as_int(r.get("n_rdrp_domains"))
        composition = vmr.get(fam, "")
        in_ictv = bool(composition) or not vmr

        fails, notes = [], []
        if fam in targets:
            notes.append("target family, always kept")
        else:
            if vmr and not in_ictv and not args.allow_unlisted:
                fails.append("not an ICTV family name")
            if n_contigs < args.min_contigs:
                fails.append(f"{n_contigs} contig(s) < {args.min_contigs}")
            if n_runs < args.min_runs:
                fails.append(f"{n_runs} run(s) < {args.min_runs}")
            if n_bp < args.min_bioprojects:
                fails.append(f"{n_bp} bioproject(s) < {args.min_bioprojects}")
            if (args.require_rdrp_for_rna and expects_rdrp(composition)
                    and n_rdrp == 0):
                fails.append(f"no RdRp domain, but {composition} expects one")

        # 3: tier by WHICH criteria failed, so replication-only failures stay
        # visible instead of being lumped in with name and evidence failures
        replication_only = bool(fails) and all(
            "run(s)" in f or "bioproject(s)" in f for f in fails)
        if not fails:
            tier = "A"
        elif replication_only and n_contigs >= args.min_contigs:
            tier = "B"
        else:
            tier = "C"
        counts[tier] += 1

        out.append({
            "assigned_family": fam, "tier": tier,
            "passes": "1" if tier in keep_tiers else "0",
            "n_contigs": n_contigs, "n_runs": n_runs, "n_bioprojects": n_bp,
            "n_rdrp_domains": n_rdrp,
            "genome_composition": composition,
            "in_ictv": "1" if in_ictv else "0",
            "expects_rdrp": "1" if expects_rdrp(composition) else "0",
            "failed_criteria": "; ".join(fails),
            "notes": "; ".join(notes),
        })

    out.sort(key=lambda r: (r["tier"], -r["n_contigs"], r["assigned_family"]))
    cols = ["assigned_family", "tier", "passes", "n_contigs", "n_runs",
            "n_bioprojects", "n_rdrp_domains", "genome_composition", "in_ictv",
            "expects_rdrp", "failed_criteria", "notes"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)) or ".", exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    passing = [r["assigned_family"] for r in out if r["passes"] == "1"]
    if args.out_families:
        with open(args.out_families, "w") as fh:
            fh.write("".join(f"{f}\n" for f in sorted(passing)))
    if args.out_yaml:
        with open(args.out_yaml, "w") as fh:
            fh.write("# Generated by filter_ictv_families.py -- do not hand-edit.\n")
            fh.write(f"# Criteria: >={args.min_contigs} contigs, "
                     f">={args.min_runs} runs, >={args.min_bioprojects} "
                     f"bioprojects, ICTV-listed, RdRp where the VMR expects it. "
                     f"Tiers kept: {','.join(sorted(keep_tiers))}.\n")
            fh.write("phylogenetics:\n  families:\n")
            for f in sorted(passing):
                fh.write(f"    - {f}\n")

    log(f"[families] tier A {counts['A']}, tier B {counts['B']}, "
        f"tier C {counts['C']}; {len(passing)} pass (tiers "
        f"{','.join(sorted(keep_tiers))})")
    for r in out[:8]:
        if r["tier"] != "A":
            break
        log(f"    A  {r['assigned_family']:22} {r['n_contigs']:5} contigs, "
            f"{r['n_runs']:3} runs, {r['genome_composition']}")
    log(f"[families] -> {args.out_tsv}")
    return 0 if passing else 1


if __name__ == "__main__":
    raise SystemExit(main())
