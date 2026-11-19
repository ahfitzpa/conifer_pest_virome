#!/usr/bin/env python3
"""Fetch ICTV exemplars for the additional candidate families (reuses ncbi cache)."""
import io, re, sys, time
from pathlib import Path
import pandas as pd
from Bio import SeqIO

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"C:\Users\Amy\OneDrive - University College Dublin\Documents\conifer_pest_virome")
REF = ROOT / "reference_data"
CACHE = REF / "ncbi"
CACHE.mkdir(exist_ok=True)
FAMS = ["Aliusviridae", "Phasmaviridae", "Orthomyxoviridae",
        "Iflaviridae", "Tombusviridae"]
RDRP_KEYS = ["rna-directed rna polymerase", "rna-dependent rna polymerase", "rdrp",
             "replicase", "l protein", "nonstructural polyprotein",
             "non-structural polyprotein", "polymerase"]
CAPSID_KEYS = ["capsid", "structural polyprotein"]


def ascii_safe(s):
    return re.sub(r"[^A-Za-z0-9_.\-]+", "_", str(s)).strip("_")


def fetch(acc):
    p = CACHE / f"{acc}.gb"
    if p.exists():
        return p.read_text(encoding="utf-8", errors="ignore")
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=gb&retmode=text")
    txt = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "ignore")
    p.write_text(txt, encoding="utf-8")
    time.sleep(0.4)
    return txt


import urllib.request

vmr = pd.read_excel(REF / "VMR_MSL41.v1.20260729.xlsx", sheet_name="VMR MSL41", header=0)
sel = vmr[vmr.Family.isin(FAMS) & (vmr["Exemplar or additional isolate"] == "E")]
out, seen = [], set()
for _, row in sel.iterrows():
    accs = re.findall(r"[A-Z]{1,2}\d{5,9}", str(row["Virus GENBANK accession"]))
    for acc in accs:
        if acc in seen:
            continue
        seen.add(acc)
        rec = SeqIO.read(io.StringIO(fetch(acc)), "genbank")
        for f in rec.features:
            if f.type != "CDS":
                continue
            prod = f.qualifiers.get("product", [""])[0]
            pl = prod.lower()
            kind = ("rdrp" if any(k in pl for k in RDRP_KEYS)
                    else "capsid" if any(k in pl for k in CAPSID_KEYS) else None)
            # Iflaviridae/Tombusviridae annotate the single ORF simply as
            # "polyprotein"; its C-terminal region contains the RdRp.
            if kind is None and pl.strip() == "polyprotein" and \
                    row.Family in ("Iflaviridae", "Tombusviridae"):
                kind = "rdrp"
            if kind is None:
                continue
            seq = str(f.extract(rec.seq).translate(to_stop=False)).split("*")[0]
            if len(seq) < 100:
                continue
            out.append({"family": row.Family, "genus": row.Genus, "species": row.Species,
                        "virus_name": row["Virus name(s)"], "accession": acc,
                        "product": prod, "kind": kind, "length": len(seq), "seq": seq})
    print(f"{row.Family[:9]:10s} {str(row.Species)[:45]:45s} {accs}")

meta = pd.DataFrame(out)
with open(REF / "ictv_refs_extra.faa", "w", encoding="utf-8") as fh:
    for i, r in enumerate(meta.itertuples()):
        hid = (f"ICTV|{ascii_safe(r.family)}|{ascii_safe(r.genus)}|{ascii_safe(r.species)}|"
               f"{r.accession}|{ascii_safe(r.kind)}|{i}")
        fh.write(f">{hid}\n{r.seq}\n")
meta.drop(columns=["seq"]).to_csv(REF / "ictv_refs_extra_meta.tsv", sep="\t", index=False)
print("\nproteins:", len(meta))
print(meta.groupby(["family", "kind"]).size().to_string())
