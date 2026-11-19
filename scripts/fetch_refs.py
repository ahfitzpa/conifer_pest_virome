#!/usr/bin/env python3
"""Fetch ICTV VMR exemplar genomes from NCBI and extract RdRp/L and capsid proteins.

Caches raw GenBank records under reference_data/ncbi/ so reruns are offline.
Outputs:
  reference_data/ictv_refs.faa       all extracted reference proteins
  reference_data/ictv_refs_meta.tsv  one row per extracted protein
"""
from __future__ import annotations

import io
import re
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(r"C:\Users\Amy\OneDrive - University College Dublin\Documents\conifer_pest_virome")
REF = ROOT / "reference_data"
CACHE = REF / "ncbi"
CACHE.mkdir(parents=True, exist_ok=True)
FAMILIES = ["Dicistroviridae", "Chuviridae"]

RDRP_KEYS = [
    "rna-directed rna polymerase", "rna-dependent rna polymerase", "rdrp",
    "replicase", "l protein", "nonstructural polyprotein",
    "non-structural polyprotein", "polymerase",
]
CAPSID_KEYS = ["capsid", "structural polyprotein"]


def ascii_safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-]+", "_", str(s)).strip("_")


def fetch_genbank(acc: str):
    path = CACHE / f"{acc}.gb"
    if path.exists():
        return path.read_text(encoding="utf-8", errors="ignore")
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=gb&retmode=text")
    txt = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "ignore")
    path.write_text(txt, encoding="utf-8")
    time.sleep(0.4)  # stay under NCBI's 3 requests/sec guideline
    return txt


def classify(product: str) -> str | None:
    p = product.lower()
    if any(k in p for k in RDRP_KEYS):
        return "rdrp"
    if any(k in p for k in CAPSID_KEYS):
        return "capsid"
    return None


def main():
    vmr = pd.read_excel(REF / "VMR_MSL41.v1.20260729.xlsx",
                        sheet_name="VMR MSL41", header=0)
    sel = vmr[(vmr["Family"].isin(FAMILIES)) &
              (vmr["Exemplar or additional isolate"] == "E")].copy()
    print("ICTV exemplar species:", len(sel))

    records = []
    seen = set()
    for _, row in sel.iterrows():
        accfield = str(row["Virus GENBANK accession"])
        accessions = re.findall(r"[A-Z]{1,2}\d{5,9}", accfield)
        for acc in accessions:
            if acc in seen:
                continue
            seen.add(acc)
            txt = fetch_genbank(acc)
            rec = SeqIO.read(io.StringIO(txt), "genbank")
            for f in rec.features:
                if f.type != "CDS":
                    continue
                product = f.qualifiers.get("product", [""])[0]
                kind = classify(product)
                if kind is None:
                    continue
                seq = f.extract(rec.seq).translate(to_stop=False)
                seq = str(seq).split("*")[0]  # stop at first stop codon
                if len(seq) < 100:
                    continue
                records.append({
                    "family": row["Family"],
                    "genus": row["Genus"],
                    "species": row["Species"],
                    "virus_name": row["Virus name(s)"],
                    "accession": acc,
                    "product": product,
                    "kind": kind,
                    "length": len(seq),
                    "seq": seq,
                })
        print(f"  {row['Family'][:4]} {row['Species']!r:55s} "
              f"acc={accessions} -> {sum(r['accession'] in accessions for r in records)} proteins")

    meta = pd.DataFrame(records)
    with open(REF / "ictv_refs.faa", "w", encoding="utf-8") as fh:
        for i, r in enumerate(meta.itertuples()):
            hid = (f"ICTV|{ascii_safe(r.family)}|{ascii_safe(r.genus)}|"
                   f"{ascii_safe(r.species)}|{r.accession}|"
                   f"{ascii_safe(r.kind)}|{i}")
            fh.write(f">{hid}\n{r.seq}\n")
    meta.drop(columns=["seq"]).to_csv(REF / "ictv_refs_meta.tsv", sep="\t", index=False)
    print("\nextracted proteins:", len(meta))
    print(meta.groupby(["family", "kind"]).size().to_string())


if __name__ == "__main__":
    main()
