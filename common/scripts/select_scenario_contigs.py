#!/usr/bin/env python3
"""Select one scenario's contigs for a family from phase2's master database.

Shared by every analysis pathway: each one trees its own manifest's samples, so
the selection rule is the same and only the manifest differs.

Selection is driven by phase2's master_metadata.csv, NOT by the FASTA header.
The header's `family` field is the SAMPLE's family (every row of a Hepeviridae
manifest carries it), so filtering on it keeps essentially every contig in the
run. The per-contig assignment is `consensus_family`, with `family_agreement`
recording whether the tools agreed.

Emits a FASTA with short tree-safe labels plus a table mapping each label back to
its run, host, evidence and original contig name. Nothing downstream should parse
the header again: join on the table instead.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

NULL = {"", "nan", "na", "none", "null", "n/a", "-"}
AGREEMENT_RANK = {"agree": 0, "single_tool": 1, "conflict": 2}
HEADER_FIELDS = ["run_id", "family", "bioproject", "family_context",
                 "detected_by", "host", "contig_name"]
TABLE_COLUMNS = ["label", "run_id", "contig_name", "length", "consensus_family",
                 "family_agreement", "detected_by", "n_tools", "checkv_quality",
                 "checkv_completeness", "host", "bioproject", "origin"]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def clean(v: str) -> str:
    v = (v or "").strip()
    return "" if v.lower() in NULL else v


def as_float(v: str) -> float:
    try:
        return float(clean(v))
    except ValueError:
        return -1.0


def read_fasta(path):
    # 1: stream the fasta; the master database is large, so never hold it all
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def parse_header(header: str) -> dict[str, str]:
    # 2: split on pipes; the contig name may itself contain one
    parts = header.split("|")
    row = {f: clean(parts[i]) if i < len(parts) else ""
           for i, f in enumerate(HEADER_FIELDS)}
    if len(parts) > len(HEADER_FIELDS):
        row["contig_name"] = "|".join(parts[len(HEADER_FIELDS) - 1:])
    return row


def manifest_runs(path: str, columns: tuple[str, ...]) -> set[str]:
    # 3: accept either a manifest csv/tsv or a plain list of run ids
    with open(path, newline="") as fh:
        first = fh.readline()
        fh.seek(0)
        if "," not in first and "\t" not in first:
            return {ln.split(":")[0].strip() for ln in fh if ln.strip()}
        delim = "\t" if "\t" in first else ","
        reader = csv.DictReader(fh, delimiter=delim)
        present = [c for c in columns if c in (reader.fieldnames or [])]
        if not present:
            sys.exit(f"none of {columns} in {path}; found: {reader.fieldnames}")
        runs = set()
        for r in reader:
            for c in present:
                if clean(r.get(c, "")):
                    runs.add(clean(r[c]))
                    break
        return runs


def load_contig_evidence(path: str, runs: set[str], family: str) -> dict:
    # 4: per-contig evidence for this family only, keyed by (run_id, contig_id)
    keep = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("run_id") not in runs:
                continue
            if clean(r.get("consensus_family", "")) != family:
                continue
            keep[(r["run_id"], clean(r.get("contig_id", "")))] = {
                "family_agreement": clean(r.get("family_agreement", "")),
                "detected_by": clean(r.get("detected_by", "")),
                "n_tools": clean(r.get("n_tools", "")),
                "checkv_quality": clean(r.get("checkv_quality", "")),
                "checkv_completeness": clean(r.get("checkv_completeness", "")),
            }
    return keep


def safe_label(run_id: str, contig_name: str, taken: set[str]) -> str:
    # 5: newick, PAML and HyPhy dislike punctuation and long names, so build a
    # short label from the run id plus the contig's NODE number
    m = re.search(r"NODE_(\d+)", contig_name) or re.search(r"(\d+)", contig_name)
    base = re.sub(r"[^A-Za-z0-9_]", "_", f"{run_id}_n{m.group(1) if m else '1'}")[:40]
    label, i = base, 1
    while label in taken:
        i += 1
        label = f"{base}_{i}"
    taken.add(label)
    return label


def write_fasta_record(fh, label: str, seq: str):
    fh.write(f">{label}\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i + 60] + "\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--master-fasta", required=True)
    p.add_argument("--metadata", required=True,
                   help="phase2 results/taxonomy/master_metadata.csv")
    p.add_argument("--manifest", required=True)
    p.add_argument("--manifest-columns", default="sample_id,run_id,run")
    p.add_argument("--family", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-table", required=True)
    p.add_argument("--references", help="extra FASTA appended as origin=reference")
    p.add_argument("--min-length", type=int, default=1000)
    p.add_argument("--max-length", type=int, default=0,
                   help="0 = no cap; set it to exclude host/misassembly contigs "
                        "far above the expected genome size")
    p.add_argument("--agreement", default="agree,single_tool",
                   help="family_agreement values accepted; 'conflict' is excluded "
                        "by default")
    p.add_argument("--max-per-sample", type=int, default=1,
                   help="keep the N best contigs per run (0 = no cap)")
    p.add_argument("--max-taxa", type=int, default=0)
    p.add_argument("--min-taxa", type=int, default=4)
    p.add_argument("--keep-contigs",
                   help="phase2 contigs_keep.txt: only these contigs are "
                        "eligible, so host genome the annotation rejected can "
                        "never reach an analysis")
    args = p.parse_args()

    accepted_agreement = {a.strip() for a in args.agreement.split(",") if a.strip()}
    runs = manifest_runs(args.manifest, tuple(c.strip() for c in
                                              args.manifest_columns.split(",")))
    if not runs:
        sys.exit(f"no run ids in {args.manifest}")

    # 4b: the phase2 confidence filter is applied before anything else, so a
    # contig the annotation called host genome is not eligible at all
    allowed: set[str] | None = None
    if args.keep_contigs and os.path.exists(args.keep_contigs):
        with open(args.keep_contigs) as fh:
            allowed = set()
            for line in fh:
                name = line.strip()
                if name:
                    allowed.add(name)
                    allowed.add(name.split("|")[-1])
        log(f"[select] confidence filter: {len(allowed)} eligible contig key(s)")

    evidence = load_contig_evidence(args.metadata, runs, args.family)
    if allowed is not None:
        before = len(evidence)
        evidence = {k: v for k, v in evidence.items()
                    if k[1] in allowed or f"{k[0]}|{k[1]}" in allowed}
        log(f"[select] {before} -> {len(evidence)} contig(s) after the "
            f"confidence filter")
    log(f"[select] {len(runs)} run(s) in the manifest; {len(evidence)} contig(s) "
        f"assigned to {args.family} in the metadata")
    if not evidence:
        log(f"[select] no contig has consensus_family == {args.family}")
        return 1

    # 6: one pass over the master database, keeping only assigned contigs
    stats = {"scanned": 0, "assigned": 0, "agreement": 0, "length": 0}
    per_run: dict[str, list] = {}
    for header, seq in read_fasta(args.master_fasta):
        stats["scanned"] += 1
        meta = parse_header(header)
        ev = evidence.get((meta["run_id"], meta["contig_name"]))
        if ev is None:
            continue
        stats["assigned"] += 1
        if accepted_agreement and ev["family_agreement"] not in accepted_agreement:
            continue
        stats["agreement"] += 1
        if len(seq) < args.min_length:
            continue
        if args.max_length and len(seq) > args.max_length:
            continue
        stats["length"] += 1
        per_run.setdefault(meta["run_id"], []).append((meta, ev, seq))
    log(f"[select] scanned {stats['scanned']}, matched {stats['assigned']}, "
        f"passed agreement {stats['agreement']}, passed length {stats['length']}")

    # 7: rank within a run by evidence, not by length alone -- the longest contig
    # in a run is often host sequence or a misassembly
    def rank_key(item):
        meta, ev, seq = item
        return (AGREEMENT_RANK.get(ev["family_agreement"], 9),
                -as_float(ev["checkv_completeness"]),
                -as_float(ev["n_tools"]),
                -len(seq),
                meta["contig_name"])

    chosen = []
    for run_id in sorted(per_run):
        entries = sorted(per_run[run_id], key=rank_key)
        chosen.extend(entries[: args.max_per_sample] if args.max_per_sample > 0
                      else entries)

    if args.max_taxa and len(chosen) > args.max_taxa:
        chosen = sorted(chosen, key=rank_key)[: args.max_taxa]
        log(f"[select] capped to the {args.max_taxa} best-supported contigs")
    chosen.sort(key=lambda e: (e[0]["run_id"], e[0]["contig_name"]))

    taken: set[str] = set()
    rows: list[dict[str, str]] = []
    os.makedirs(os.path.dirname(os.path.abspath(args.out_fasta)) or ".", exist_ok=True)
    with open(args.out_fasta, "w") as out:
        for meta, ev, seq in chosen:
            label = safe_label(meta["run_id"], meta["contig_name"], taken)
            rows.append({"label": label, "run_id": meta["run_id"],
                         "contig_name": meta["contig_name"], "length": str(len(seq)),
                         "consensus_family": args.family, "host": meta["host"],
                         "bioproject": meta["bioproject"], "origin": "candidate",
                         **ev})
            write_fasta_record(out, label, seq)

        n_ref = 0
        if args.references and os.path.exists(args.references):
            for header, seq in read_fasta(args.references):
                label = re.sub(r"[^A-Za-z0-9_]", "_", header.split()[0])[:40]
                while label in taken:
                    label += "_r"
                taken.add(label)
                rows.append({"label": label, "run_id": "", "contig_name": header,
                             "length": str(len(seq)), "consensus_family": args.family,
                             "origin": "reference"})
                write_fasta_record(out, label, seq)
                n_ref += 1
            log(f"[select] added {n_ref} reference sequence(s)")

    with open(args.out_table, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TABLE_COLUMNS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n_cand = sum(1 for r in rows if r["origin"] == "candidate")
    lengths = sorted(int(r["length"]) for r in rows if r["origin"] == "candidate")
    if lengths:
        log(f"[select] candidate lengths: min {lengths[0]}, "
            f"median {lengths[len(lengths)//2]}, max {lengths[-1]}")
    log(f"[select] {n_cand} candidate + {len(rows) - n_cand} reference "
        f"= {len(rows)} taxa -> {args.out_fasta}")
    if n_cand < args.min_taxa:
        log(f"[select] only {n_cand} candidate(s); need {args.min_taxa} for a tree")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
