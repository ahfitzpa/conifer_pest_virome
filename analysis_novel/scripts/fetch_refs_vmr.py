#!/usr/bin/env python3
"""Fetch ICTV exemplar RdRp/capsid proteins for the families being analysed.

Replaces fetch_refs.py and fetch_refs_extra.py, whose family lists were
hardcoded to an earlier project: 24 of this project's 25 analysed families had
no references at all, so every tree was candidates-only and demarcation had
nothing to compare against.

Two things it does that the old scripts did not:

  * matches a target against ANY VMR rank, not just Family. The generated
    family list carries orders and classes (Amarillovirales, Tolivirales,
    Duplopiviricetes, Duplornaviricota) because the annotation tools report
    at mixed ranks, and a Family-only match silently returned nothing.
  * caps the exemplars per target, picked round-robin across genera, so a
    large family does not dominate (Steitzviridae alone has 1289 exemplars).

The `family` column of the metadata holds the TARGET label, not the VMR
family, because scripts/phylo_prepare.py selects references by matching that
column against the analysed family name.

Outputs (default): references/ictv_refs_all.faa, ictv_refs_all_meta.tsv
Caches raw GenBank under references/ncbi/ so reruns are offline.
"""
from __future__ import annotations

import argparse
import csv
import io
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT.parent / "references"
CACHE = REF / "ncbi"

RANKS = ["Family", "Order", "Class", "Phylum", "Genus", "Species"]
RDRP_KEYS = [
    "rna-directed rna polymerase", "rna-dependent rna polymerase", "rdrp",
    "replicase", "l protein", "nonstructural polyprotein",
    "non-structural polyprotein", "polymerase",
    # mature-peptide names for the RdRp of the single-polyprotein families:
    # flaviviruses call it NS5, picorna-like viruses 3D.
    "ns5", "3d", "protein 3d",
]
CAPSID_KEYS = ["capsid", "structural polyprotein", "vp1", "coat protein"]
# Last resort for a record whose only CDS is an unannotated polyprotein. Kept
# distinct from rdrp because it is the whole ORF, not the marker domain.
POLY_KEYS = ["polyprotein", "orf1", "orf 1"]
# "polymerase" alone has to stay in RDRP_KEYS -- several RNA families annotate
# their RdRp with nothing else (Chuviridae, Nyamiviridae). It also matches the
# DNA polymerase of every large DNA virus and phage, which is why an earlier
# run extracted 4019 "rdrp" proteins for Mimiviridae and 2299 for
# Chimalliviridae. These veto a product before it can be called an RdRp.
NOT_RDRP_KEYS = [
    "dna polymerase", "dna-directed", "dna directed", "dna pol", "primase",
    "reverse transcriptase", "rna polymerase ii", "poly(a)", "polyadenylate",
    "polymerase accessory", "polymerase-associated", "processivity",
    "terminase", "topoisomerase", "nucleotidyltransferase",
    "rna polymerase beta", "polymerase binding", "ribosylase",
]

# Matched on word boundaries, not as bare substrings: "l protein" (the
# mononegavirus RdRp) otherwise matches "hypotheticaL PROTEIN", which labelled
# 4646 Mimiviridae and 2009 Chimalliviridae hypothetical proteins as RdRp.
_RDRP_RE = re.compile(r"\b(?:" + "|".join(re.escape(k) for k in RDRP_KEYS) + r")\b")
_CAPSID_RE = re.compile(r"\b(?:" + "|".join(re.escape(k) for k in CAPSID_KEYS) + r")\b")
_POLY_RE = re.compile(r"\b(?:" + "|".join(re.escape(k) for k in POLY_KEYS) + r")\b")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def ascii_safe(s) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-]+", "_", str(s)).strip("_")


def classify(product: str):
    p = product.lower().strip()
    if any(k in p for k in NOT_RDRP_KEYS):
        pass
    elif _RDRP_RE.search(p):
        return "rdrp"
    if _CAPSID_RE.search(p):
        return "capsid"
    if _POLY_RE.search(p):
        return "polyprotein"
    return None


def extract_proteins(rec, min_aa):
    """Classified proteins from one GenBank record, preferring mature peptides.

    The single-polyprotein families (Flaviviridae/Amarillovirales,
    Picornaviridae, Potyviridae, Iflaviridae) annotate one CDS called
    'polyprotein', which matched no keyword at all -- those families came back
    with zero references. Their RdRp is annotated as a mat_peptide (NS5, 3D),
    which is also the right size for a marker alignment, so mat_peptide is
    read first and the whole polyprotein is kept only when nothing better
    exists.
    """
    found = {"mat_peptide": [], "CDS": []}
    for feat in rec.features:
        if feat.type not in ("CDS", "mat_peptide"):
            continue
        product = (feat.qualifiers.get("product")
                   or feat.qualifiers.get("gene") or [""])[0]
        kind = classify(product)
        if kind is None:
            continue
        try:
            aa = str(feat.extract(rec.seq).translate(to_stop=False))
        except Exception:  # malformed or trans-spliced coordinates
            continue
        aa = aa.split("*")[0]
        if len(aa) < min_aa:
            continue
        found[feat.type].append((kind, product, aa))
    # mature peptides give the domain; fall back to CDS when there are none
    out = [x for x in found["mat_peptide"] if x[0] in ("rdrp", "capsid")]
    if not any(k == "rdrp" for k, _, _ in out):
        out += [x for x in found["CDS"] if x[0] == "rdrp"]
    if not any(k == "capsid" for k, _, _ in out):
        out += [x for x in found["CDS"] if x[0] == "capsid"]
    if not out:
        out = [x for x in found["CDS"] + found["mat_peptide"]
               if x[0] == "polyprotein"]
    return out


def fetch_genbank(acc: str, retries: int = 3):
    """Return the GenBank text for an accession, caching it under references/ncbi."""
    path = CACHE / f"{acc}.gb"
    if path.exists() and path.stat().st_size > 0:
        return path.read_text(encoding="utf-8", errors="ignore")
    url = f"{EUTILS}?db=nuccore&id={acc}&rettype=gb&retmode=text"
    for attempt in range(retries):
        try:
            txt = urllib.request.urlopen(url, timeout=120).read().decode("utf-8", "ignore")
        except (urllib.error.URLError, OSError) as exc:
            if attempt == retries - 1:
                print(f"    [warn] {acc}: {exc}", file=sys.stderr)
                return None
            time.sleep(2 * (attempt + 1))
            continue
        if not txt.strip() or "LOCUS" not in txt[:2000]:
            return None
        path.write_text(txt, encoding="utf-8")
        time.sleep(0.4)  # NCBI asks for <=3 requests/sec
        return txt
    return None


def pick_exemplars(rows, cap):
    """Up to `cap` rows, round-robin across genera so one genus cannot fill the set."""
    if cap <= 0 or len(rows) <= cap:
        return rows
    by_genus = defaultdict(list)
    for r in rows:
        by_genus[r.get("Genus") or "NA"].append(r)
    picked, genera = [], sorted(by_genus)
    while len(picked) < cap:
        progressed = False
        for g in genera:
            if by_genus[g]:
                picked.append(by_genus[g].pop(0))
                progressed = True
                if len(picked) == cap:
                    break
        if not progressed:
            break
    return picked


def resolve_target(name, rows_by_rank):
    """Find the rank at which a target name matches, Family first."""
    for rank in RANKS:
        hits = rows_by_rank[rank].get(name.lower())
        if hits:
            return rank, hits
    return None, []


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vmr", default=str(REF / "vmr_msl41.tsv"),
                    help="VMR TSV from vmr_xlsx_to_tsv.py")
    ap.add_argument("--families", required=True,
                    help="comma-separated targets, or a YAML/text file listing them")
    ap.add_argument("--per-family", type=int, default=12,
                    help="max exemplar isolates per target (0 = no cap)")
    ap.add_argument("--out-faa", default=str(REF / "ictv_refs_all.faa"))
    ap.add_argument("--out-meta", default=str(REF / "ictv_refs_all_meta.tsv"))
    ap.add_argument("--out-nt-dir", default="",
                    help="also write <Target>_refs.fasta genome nucleotide sets "
                         "here, for the shared nucleotide phylogeny")
    ap.add_argument("--min-aa", type=int, default=100)
    a = ap.parse_args()

    CACHE.mkdir(parents=True, exist_ok=True)

    # 1: target list, from a file (YAML list or one per line) or a comma string
    p = Path(a.families)
    if p.exists():
        targets = [ln.strip().lstrip("- ").strip().strip("'\"")
                   for ln in p.read_text().splitlines()]
        targets = [t for t in targets
                   if t and not t.startswith("#") and not t.endswith(":")]
    else:
        targets = [t.strip() for t in a.families.split(",") if t.strip()]
    targets = list(dict.fromkeys(targets))
    print(f"[refs] {len(targets)} target(s)")

    # 2: index the VMR by every rank we are willing to match on
    with open(a.vmr, encoding="utf-8", newline="") as fh:
        vmr = list(csv.DictReader(fh, delimiter="\t"))
    rows_by_rank = {r: defaultdict(list) for r in RANKS}
    for row in vmr:
        for rank in RANKS:
            v = (row.get(rank) or "").strip()
            if v:
                rows_by_rank[rank][v.lower()].append(row)

    # 3: per target, resolve the rank, pick exemplars, extract proteins
    records, unresolved, empty = [], [], []
    genomes: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for target in targets:
        rank, hits = resolve_target(target, rows_by_rank)
        if not hits:
            unresolved.append(target)
            print(f"[refs] {target:24} -- not in VMR MSL41")
            continue
        chosen = pick_exemplars(hits, a.per_family)
        n_before = len(records)
        seen = set()
        for row in chosen:
            accs = re.findall(r"[A-Z]{1,2}\d{5,9}",
                              str(row.get("Virus GENBANK accession") or ""))
            for acc in accs:
                if acc in seen:
                    continue
                seen.add(acc)
                txt = fetch_genbank(acc)
                if not txt:
                    continue
                try:
                    rec = SeqIO.read(io.StringIO(txt), "genbank")
                except (ValueError, OSError) as exc:
                    print(f"    [warn] {acc}: unreadable ({exc})", file=sys.stderr)
                    continue
                if a.out_nt_dir and len(rec.seq):
                    genomes[target].append((
                        f"ICTV_{acc}_{ascii_safe(row.get('Species') or 'NA')}",
                        str(rec.seq)))
                for kind, product, aa in extract_proteins(rec, a.min_aa):
                    records.append({
                        "family": target,          # the analysed label
                        "matched_rank": rank,
                        "vmr_family": row.get("Family") or "NA",
                        "genus": row.get("Genus") or "NA",
                        "species": row.get("Species") or "NA",
                        "virus_name": row.get("Virus name(s)") or "NA",
                        "accession": acc,
                        "product": product,
                        "kind": kind,
                        "length": len(aa),
                        "seq": aa,
                    })
        got = len(records) - n_before
        if got == 0:
            empty.append(target)
        print(f"[refs] {target:24} rank={rank:8} isolates={len(chosen):4} proteins={got}")

    if not records:
        print("[refs] ERROR: no proteins extracted", file=sys.stderr)
        return 1

    # 4: write FASTA; headers keep |family| and |accession| so phylo_prepare matches
    with open(a.out_faa, "w", encoding="utf-8") as fh:
        for i, r in enumerate(records):
            fh.write(f">ICTV|{ascii_safe(r['family'])}|{ascii_safe(r['genus'])}|"
                     f"{ascii_safe(r['species'])}|{r['accession']}|"
                     f"{ascii_safe(r['kind'])}|{i}\n{r['seq']}\n")
    cols = ["family", "matched_rank", "vmr_family", "genus", "species",
            "virus_name", "accession", "product", "kind", "length"]
    with open(a.out_meta, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in records:
            w.writerow({k: r[k] for k in cols})

    # 5: per-target genome nucleotide sets for the shared phylogeny
    if a.out_nt_dir:
        nt_dir = Path(a.out_nt_dir)
        nt_dir.mkdir(parents=True, exist_ok=True)
        for target, recs in sorted(genomes.items()):
            with open(nt_dir / f"{ascii_safe(target)}_refs.fasta", "w",
                      encoding="utf-8") as fh:
                for name, seq in recs:
                    fh.write(f">{name}\n{seq}\n")
        print(f"[refs] wrote {len(genomes)} genome set(s) to {nt_dir}")

    # 6: summary
    per_fam = defaultdict(lambda: defaultdict(int))
    for r in records:
        per_fam[r["family"]][r["kind"]] += 1
    print(f"\n[refs] {len(records)} proteins across {len(per_fam)} target(s)")
    for fam in sorted(per_fam):
        k = per_fam[fam]
        print(f"  {fam:24} rdrp={k['rdrp']:4} capsid={k['capsid']:4} "
              f"polyprotein={k['polyprotein']:4}")
    if unresolved:
        print(f"\n[refs] not in VMR MSL41 ({len(unresolved)}): {', '.join(unresolved)}")
    if empty:
        print(f"[refs] matched but no RdRp/capsid CDS ({len(empty)}): {', '.join(empty)}")
    print(f"\n[refs] wrote {a.out_faa}\n[refs] wrote {a.out_meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
