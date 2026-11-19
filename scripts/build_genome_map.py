#!/usr/bin/env python3
"""To-scale genome maps of candidate contigs (Zowalaty et al., Fig. 3 style).

Two modes keep data preparation (tax_core env) separate from plotting
(genome_maps env):

  --prepare-only  ORF annotation + contig lengths + coverage + CheckV
                  -> one tidy genome_map_data.tsv (audit trail)
  --plot-only     TSV -> PDF/SVG/300-dpi PNG, ORFs offset by reading frame
                  and coloured by product keyword.
"""
from __future__ import annotations
import argparse, csv, sys
from collections import defaultdict

import yaml


def read_tsv(path):
    try:
        with open(path, encoding="utf-8", errors="ignore") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except FileNotFoundError:
        return []


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


def colour_for(product, cfg):
    for key, colour in cfg["product_colours"].items():
        if key.lower() in (product or "").lower():
            return colour
    return cfg.get("hypothetical_colour", "#e7e3da")


def abbreviate(product):
    p = (product or "").lower()
    for key, short in [("rna-directed rna polymerase", "L"), ("rna-dependent rna polymerase", "RdRp"),
                       ("replicase", "Rep"), ("nucleoprotein", "N"), ("glycoprotein", "GP"),
                       ("structural", "Cap"), ("capsid", "Cap"), ("polyprotein", "PP"),
                       ("helicase", "Hel"), ("protease", "Pro")]:
        if key in p:
            return short
    return "ORF"


def prepare(a):
    cfg = yaml.safe_load(open(a.config, encoding="utf-8"))["phylogenetics"]["genome_maps"]
    candidates = [r for r in read_tsv(a.candidates)
                  if (r.get("assigned_family") or "").strip() == a.family]
    lengths = {name.split("|")[-1]: len(seq) for name, seq in iter_fasta(a.master_fasta)}
    coverage = {}
    for path in a.coverage:
        for r in read_tsv(path):
            cid = r.get("contig_id", "")
            if cid:
                coverage[cid] = r

    orf_by_contig = defaultdict(list)
    for r in read_tsv(a.orf_annotation):
        target = (r.get("target") or "").strip()
        cid = target.split("|")[-1]
        if cid:
            orf_by_contig[cid].append(r)

    fields = ["family", "sample", "contig_id", "contig_length", "coverage_mean",
              "checkv_quality", "checkv_completeness", "host", "ORF_start", "ORF_end",
              "frame", "product", "best_hit", "evalue", "product_colour"]
    rows = []
    for c in candidates:
        cid = (c.get("contig_id") or "").strip()
        sample = c.get("run", "")
        L = int(float(coverage.get(cid, {}).get("contig_length") or lengths.get(cid) or 0))
        cov = coverage.get(cid, {}).get("coverage_mean", "")
        orfs = orf_by_contig.get(cid, [])
        if not orfs:
            rows.append(dict(family=a.family, sample=sample, contig_id=cid,
                             contig_length=L, coverage_mean=cov,
                             checkv_quality=c.get("checkv_quality", ""),
                             checkv_completeness=c.get("checkv_completeness", ""),
                             host=c.get("host_organism", ""), ORF_start=0, ORF_end=0,
                             frame=0, product="(no ORF call)", best_hit="", evalue="",
                             product_colour=cfg["hypothetical_colour"]))
            continue
        for r in orfs:
            start = int(float(r.get("ORF_start") or 0))
            end = int(float(r.get("ORF_end") or 0))
            prod = r.get("product") or "Hypothetical"
            rows.append(dict(family=a.family, sample=sample, contig_id=cid,
                             contig_length=L, coverage_mean=cov,
                             checkv_quality=c.get("checkv_quality", ""),
                             checkv_completeness=c.get("checkv_completeness", ""),
                             host=c.get("host_organism", ""),
                             ORF_start=start, ORF_end=end, frame=start % 3,
                             product=prod, best_hit=r.get("best_hit", ""),
                             evalue=r.get("evalue", ""),
                             product_colour=colour_for(prod, cfg)))
    with open(a.out_table, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"{a.family}: {len(rows)} map rows for "
          f"{len({r['contig_id'] for r in rows})} contigs", file=sys.stderr)


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    cfg = yaml.safe_load(open(a.config, encoding="utf-8"))["phylogenetics"]["genome_maps"]
    rows = read_tsv(a.in_table)
    if not rows:
        raise SystemExit("empty genome map table")
    by_sample = defaultdict(lambda: defaultdict(list))
    meta = {}
    for r in rows:
        by_sample[r["sample"]][r["contig_id"]].append(r)
        meta[(r["sample"], r["contig_id"])] = r

    samples = sorted(by_sample)
    fig, axes = plt.subplots(len(samples), 1, figsize=(13, max(3.4, 2.3 * len(samples))),
                             squeeze=False, constrained_layout=True)
    max_len = max(float(r["contig_length"]) for r in rows)
    for si, sample in enumerate(samples):
        ax = axes[si][0]
        contigs = sorted(by_sample[sample])
        for ci, cid in enumerate(contigs):
            y = len(contigs) - 1 - ci
            L = float(meta[(sample, cid)]["contig_length"])
            ax.add_patch(Rectangle((0, y - 0.28), L, 0.56, facecolor="#e7e3da",
                                   edgecolor="#cbc6bb", linewidth=0.5))
            labelled = 0
            for r in sorted(by_sample[sample][cid], key=lambda x: int(x["ORF_start"] or 0)):
                s, e = float(r["ORF_start"]), float(r["ORF_end"])
                if e <= s:
                    continue
                # reading-frame offset, as in Zowalaty et al. Fig. 3
                fy = y - 0.22 + (int(r["frame"]) % 3) * 0.14
                ax.add_patch(Rectangle((s, fy), e - s, 0.14, facecolor=r["product_colour"],
                                       edgecolor="none"))
                if labelled < cfg.get("max_orf_label", 8) and (e - s) > max_len * 0.012:
                    ax.text(s + (e - s) / 2, y + 0.30, abbreviate(r["product"]),
                            ha="center", va="bottom", fontsize=6, color="#2a2723")
                    labelled += 1
            info = meta[(sample, cid)]
            ax.text(-max_len * 0.012, y, f"{cid}  ({info['checkv_quality'] or 'NA'}; "
                    f"{info['coverage_mean'] or 'NA'}x)", ha="right", va="center", fontsize=7)
        ax.set_xlim(-max_len * 0.34, max_len * 1.02)
        ax.set_ylim(-0.7, len(contigs) + 0.5)
        ax.set_yticks([])
        ax.set_title(sample, loc="left", fontsize=9, fontweight="bold")
        ax.grid(axis="x", visible=False)
        if si == len(samples) - 1:
            ax.set_xlabel("Genome position (nt)")
        # scale bar
        sb = max_len * 0.1
        ax.plot([max_len * 0.55, max_len * 0.55 + sb], [-0.45, -0.45], color="#2a2723", lw=1.4)
        ax.text(max_len * 0.55 + sb / 2, -0.62, f"{sb/1000:.1f} kb", ha="center", fontsize=7)
    fig.suptitle(f"{a.family}: genome organisation of candidate contigs "
                 "(ORFs offset by reading frame; colours = product)",
                 fontsize=11, fontweight="bold", x=0.01, ha="left")
    for ext, path in [("pdf", a.out_pdf), ("svg", a.out_svg), ("png", a.out_png)]:
        fig.savefig(path, bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)
    print(f"plotted {len(samples)} libraries", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--candidates")
    ap.add_argument("--orf-annotation")
    ap.add_argument("--master-fasta")
    ap.add_argument("--coverage", nargs="*", default=[])
    ap.add_argument("--out-table")
    ap.add_argument("--in-table")
    ap.add_argument("--out-pdf")
    ap.add_argument("--out-svg")
    ap.add_argument("--out-png")
    a = ap.parse_args()
    if a.prepare_only:
        prepare(a)
    elif a.plot_only:
        plot(a)
    else:
        raise SystemExit("choose --prepare-only or --plot-only")
    return 0


if __name__ == "__main__":
    sys.exit(main())
