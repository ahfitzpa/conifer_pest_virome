#!/usr/bin/env python3
"""Flatten the ICTV VMR spreadsheet to a TSV.

The TSV is the input to fetch_refs_vmr.py, so that script needs only the
standard library and Biopython -- no openpyxl in the workflow envs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT.parent / "references"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", default=str(REF / "VMR_MSL41.v1.20260729.xlsx"))
    ap.add_argument("--sheet", default="VMR MSL41")
    ap.add_argument("--out", default=str(REF / "vmr_msl41.tsv"))
    a = ap.parse_args()

    # 1: read the sheet
    df = pd.read_excel(a.xlsx, sheet_name=a.sheet, header=0)
    # 2: exemplar isolates only -- additional isolates duplicate the species
    col = "Exemplar or additional isolate"
    df = df[df[col].astype(str).str.strip().str.upper().str.startswith("E")]
    # 3: write
    df.to_csv(a.out, sep="\t", index=False)
    print(f"[vmr] {len(df)} exemplar rows -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
