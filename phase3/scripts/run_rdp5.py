#!/usr/bin/env python3
"""Drive RDP5 (standalone binary) on an alignment and parse its events.

RDP5 is not on bioconda. The binary path must be supplied in config. If it is
missing we FAIL LOUDLY (exit 3) rather than writing a 'not_assessed' table that
could be mistaken for a real negative result -- the same rule used by the
OpenRDP wrapper.
"""
from __future__ import annotations
import argparse, csv, re, subprocess, sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--eligible", required=True)
    ap.add_argument("--rdp5-binary", default="")
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    with open(a.eligible, encoding="utf-8") as fh:
        el = next(csv.DictReader(fh, delimiter="\t"), {})
    if el.get("eligible") != "1":
        print("not assessable: " + el.get("reason", "unknown") + " -- skipped", file=sys.stderr)
        fields = ["query", "donor", "breakpoint", "method", "status"]
        with open(a.out_tsv, "w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=fields, delimiter="\t").writeheader()
        return 0

    binary = Path(a.rdp5_binary) if a.rdp5_binary else None
    if binary is None or not binary.exists():
        print("ERROR: RDP5 binary not configured or not found "
              f"({a.rdp5_binary or 'empty'}). Set phylogenetics.recombination.rdp5_binary "
              "in configs/phylogenetics.yaml, or remove 'rdp5' from methods. "
              "Failing loudly.", file=sys.stderr)
        return 3

    work = Path(a.out_tsv).parent / "rdp5_work"
    work.mkdir(parents=True, exist_ok=True)
    # RDP5 uses an interactive command file; write the standard batch for a
    # single alignment scan with all seven methods.
    batch = work / "commands.txt"
    batch.write_text(
        f'"{a.alignment}"\n1\n2\n3\n4\n5\n6\n7\nY\n'
        f'"{work / "result"}"\n', encoding="utf-8")
    try:
        r = subprocess.run([str(binary)], stdin=open(batch), capture_output=True,
                           text=True, cwd=work, timeout=3600)
    except Exception as exc:
        print(f"ERROR: RDP5 run failed: {exc}", file=sys.stderr)
        return 4
    result_file = next(work.glob("result*"), None)
    if result_file is None:
        print(f"ERROR: RDP5 produced no result file. stderr: {r.stderr[-500:]}", file=sys.stderr)
        return 5

    fields = ["query", "donor", "breakpoint", "method", "status"]
    rows = []
    for line in result_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        # RDP5 report lines name the recombinant and the two parents, with
        # breakpoint positions in the first columns. Parse defensively.
        parts = re.split(r"\s+", line.strip())
        if len(parts) >= 6 and any(c.isdigit() for c in parts[0]):
            rows.append({"query": parts[1], "donor": f"{parts[2]}/{parts[3]}",
                         "breakpoint": parts[0], "method": "RDP5", "status": "RDP5_reported"})
    with open(a.out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"RDP5: {len(rows)} events", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
