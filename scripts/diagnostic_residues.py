#!/usr/bin/env python3
"""Amino-acid diagnostic residues and per-column conservation (Zhao et al., Fig. 4).

For a protein alignment annotated with group labels (e.g. candidate vs ICTV),
identify columns where groups are fixed for different residues, and report
convolution-free, named statistics:

  diagnostic_residues.tsv  position, per-group residue/fixed flag, diagnostic flag
  conservation.tsv         position, consensus, conservation fraction, Shannon entropy
"""
from __future__ import annotations
import argparse, csv, math, sys
from collections import Counter


def iter_fasta(path):
    name, chunks = None, []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_tsv(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-conservation", required=True)
    a = ap.parse_args()

    import yaml
    cfg = yaml.safe_load(open(a.config, encoding="utf-8"))["phylogenetics"]
    group_col = cfg.get("diagnostic_group_column", "group")
    groups = list(cfg.get("diagnostic_groups", ["candidate", "ICTV"]))

    seqs = dict(iter_fasta(a.alignment))
    if not seqs:
        raise SystemExit("empty alignment")
    L = len(next(iter(seqs.values())))
    if any(len(s) != L for s in seqs.values()):
        raise SystemExit("sequences are not aligned to equal length")

    meta = {r["taxon_id"]: r for r in read_tsv(a.metadata)}
    taxon_group = {t: (meta.get(t, {}).get(group_col) or "NA") for t in seqs}

    matrix, cons_rows = [], []
    for i in range(L):
        col = {t: s[i] for t, s in seqs.items()}
        residues = [c for c in col.values() if c not in "-Xx"]
        n = len(residues)
        counts = Counter(residues)
        consensus, cons_frac = ("NA", 0.0)
        if counts:
            consensus, top = counts.most_common(1)[0]
            cons_frac = round(top / n, 4)
        entropy = -sum((c / n) * math.log2(c / n) for c in counts.values()) if n else 0.0
        cons_rows.append({
            "position": i + 1, "n_non_gap": n, "n_gaps": len(col) - n,
            "consensus_aa": consensus, "conservation": cons_frac,
            "shannon_entropy": round(entropy, 4),
        })

        row = {"position": i + 1, "consensus_aa": consensus,
               "conservation": cons_frac}
        fixed = {}
        for g in groups:
            gr = [c for t, c in col.items() if taxon_group[t] == g and c not in "-Xx"]
            gcounts = Counter(gr)
            if gr and len(gcounts) == 1:
                fixed[g] = True
                row[f"{g}_aa"] = gr[0]
            else:
                fixed[g] = False
                row[f"{g}_aa"] = gcounts.most_common(1)[0][0] if gcounts else "NA"
            row[f"{g}_fixed"] = "1" if fixed[g] else "0"
        present = {row[f"{g}_aa"] for g in groups if fixed[g] and row[f"{g}_aa"] != "NA"}
        row["diagnostic"] = "1" if len(present) >= 2 else "0"
        row["n_fixed_groups"] = len(present)
        matrix.append(row)

    fields = ["position", "consensus_aa", "conservation"]
    for g in groups:
        fields += [f"{g}_aa", f"{g}_fixed"]
    fields += ["diagnostic", "n_fixed_groups"]

    with open(a.out_matrix, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(matrix)
    with open(a.out_conservation, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cons_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(cons_rows)

    n_diag = sum(r["diagnostic"] == "1" for r in matrix)
    print(f"alignment {L} aa; {n_diag} diagnostic columns", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
