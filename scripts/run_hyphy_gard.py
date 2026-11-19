#!/usr/bin/env python3
"""Run a HyPhy recombination screen (GARD, or Walrus if configured).

Uses the HyPhy CLI (`hyphy gard --alignment ... --output ...`). Parses the
resulting JSON for breakpoints and supports, and writes a tidy event table.
Fails loudly if HyPhy is absent or the run does not produce a result file.
"""
from __future__ import annotations
import argparse, csv, json, subprocess, sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--eligible", required=True)
    ap.add_argument("--method", default="gard")
    ap.add_argument("--mode", default="Faster")
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    with open(a.eligible, encoding="utf-8") as fh:
        el = next(csv.DictReader(fh, delimiter="\t"), {})
    if el.get("eligible") != "1":
        print("not assessable: " + el.get("reason", "unknown") + " -- skipped", file=sys.stderr)
        fields = ["method", "breakpoint", "support", "status"]
        with open(a.out_tsv, "w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=fields, delimiter="\t").writeheader()
        return 0

    out_json = str(Path(a.out_tsv).with_suffix(".gard.json"))
    cmd = ["hyphy", a.method, "--alignment", a.alignment,
           "--mode", a.mode, "--output", out_json]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        print("ERROR: hyphy not found -- check env_recomb.yaml", file=sys.stderr)
        return 3
    if r.returncode != 0 or not Path(out_json).exists():
        print(f"ERROR: hyphy failed (rc={r.returncode}): {r.stderr[-500:]}", file=sys.stderr)
        return 4

    try:
        data = json.loads(Path(out_json).read_text(encoding="utf-8", errors="ignore"))
    except Exception as exc:
        print(f"ERROR: cannot parse HyPhy JSON: {exc}", file=sys.stderr)
        return 5

    # GARD JSON stores breakpoints under 'best' / 'breaks' / 'evidence'
    breaks = data.get("best", {}).get("breaks", []) or data.get("breaks", [])
    fields = ["method", "breakpoint", "support", "status"]
    with open(a.out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for b in breaks:
            pos = b.get("position", b.get("breakpoint", "NA"))
            sup = b.get("support", b.get("AIC", "NA"))
            w.writerow({"method": a.method, "breakpoint": pos,
                        "support": sup, "status": "GARD_inferred"})
    print(f"hyphy {a.method}: {len(breaks)} breakpoints", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
