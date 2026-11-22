#!/usr/bin/env python3
"""Turn HyPhy FUBAR and MEME output into one per-site table and a summary.

Thresholds are explicit arguments, and each site records which test supported
it, so a category can always be traced back to the number that produced it.

HyPhy's JSON puts the numbers in `MLE.content` and their names in `MLE.headers`,
so columns are matched by name rather than by position: the order has changed
between HyPhy versions, and reading beta from the wrong column would invert the
result silently.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load_mle(path: str) -> tuple[list[str], list[list]]:
    # 1: headers + rows, tolerating a missing or unfinished file
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        log(f"[selection] no usable file at {path}")
        return [], []
    try:
        with open(path) as fh:
            doc = json.load(fh)
    except (json.JSONDecodeError, OSError) as e:
        log(f"[selection] could not parse {path}: {e}")
        return [], []
    mle = doc.get("MLE") or {}
    raw_headers = mle.get("headers") or []
    headers = [str(h[0]) if isinstance(h, (list, tuple)) and h else str(h)
               for h in raw_headers]
    content = mle.get("content") or {}
    rows: list[list] = []
    for key in sorted(content, key=lambda k: (len(str(k)), str(k))):
        block = content[key]
        if isinstance(block, list):
            rows.extend(block)
    return headers, rows


def column(headers: list[str], *wanted: str) -> int:
    # 2: first header whose name contains one of the wanted substrings
    low = [h.lower() for h in headers]
    for w in wanted:
        for i, h in enumerate(low):
            if w in h:
                return i
    return -1


def value(row: list, idx: int):
    if idx < 0 or idx >= len(row):
        return None
    try:
        return float(row[idx])
    except (TypeError, ValueError):
        return None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fubar")
    p.add_argument("--meme")
    p.add_argument("--out-sites", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--omega-positive", type=float, default=1.5)
    p.add_argument("--omega-purifying", type=float, default=0.5)
    p.add_argument("--fubar-posterior", type=float, default=0.9)
    p.add_argument("--meme-pvalue", type=float, default=0.05)
    p.add_argument("--min-sites", type=int, default=1)
    args = p.parse_args()

    f_head, f_rows = load_mle(args.fubar)
    m_head, m_rows = load_mle(args.meme)
    log(f"[selection] FUBAR {len(f_rows)} site(s), MEME {len(m_rows)} site(s)")
    if not f_rows and not m_rows:
        log("[selection] neither test produced sites; nothing to summarise")
        return 1

    i_alpha = column(f_head, "alpha", "syn")
    i_beta = column(f_head, "beta", "non-syn", "nonsyn")
    i_post = column(f_head, "prob[beta>alpha]", "posterior", "prob")
    j_pval = column(m_head, "p-value", "pvalue", "p value")
    j_beta_plus = column(m_head, "beta+", "beta +")
    for name, head, idx in (("FUBAR beta", f_head, i_beta),
                            ("MEME p-value", m_head, j_pval)):
        if head and idx < 0:
            log(f"[selection] WARNING: no {name} column in {head}")

    n = max(len(f_rows), len(m_rows))
    rows = []
    counts = {"pervasive_positive": 0, "episodic_positive": 0,
              "pervasive_purifying": 0, "neutral_or_unresolved": 0}
    for site in range(1, n + 1):
        f = f_rows[site - 1] if site <= len(f_rows) else []
        m = m_rows[site - 1] if site <= len(m_rows) else []
        alpha, beta = value(f, i_alpha), value(f, i_beta)
        post = value(f, i_post)
        pval = value(m, j_pval)
        omega = (beta / alpha) if (alpha and beta is not None and alpha > 0) else None

        # 3: pervasive first (FUBAR), then episodic (MEME), then purifying
        if (omega is not None and omega >= args.omega_positive
                and post is not None and post >= args.fubar_posterior):
            category, support = "pervasive_positive", "fubar"
        elif pval is not None and pval <= args.meme_pvalue:
            category, support = "episodic_positive", "meme"
        elif omega is not None and omega <= args.omega_purifying:
            category, support = "pervasive_purifying", "fubar"
        else:
            category, support = "neutral_or_unresolved", ""
        counts[category] += 1
        rows.append({
            "site": site,
            "alpha_syn": "" if alpha is None else f"{alpha:.5f}",
            "beta_nonsyn": "" if beta is None else f"{beta:.5f}",
            "omega": "" if omega is None else f"{omega:.5f}",
            "fubar_posterior": "" if post is None else f"{post:.5f}",
            "meme_pvalue": "" if pval is None else f"{pval:.5g}",
            "meme_beta_plus": ("" if value(m, j_beta_plus) is None
                               else f"{value(m, j_beta_plus):.5f}"),
            "selection_category": category,
            "supported_by": support,
        })

    cols = ["site", "alpha_syn", "beta_nonsyn", "omega", "fubar_posterior",
            "meme_pvalue", "meme_beta_plus", "selection_category", "supported_by"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out_sites)) or ".",
                exist_ok=True)
    with open(args.out_sites, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    omegas = [float(r["omega"]) for r in rows if r["omega"]]
    with open(args.out_summary, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["n_sites", len(rows)])
        for k, v in counts.items():
            w.writerow([f"n_{k}", v])
        if omegas:
            w.writerow(["median_omega", f"{sorted(omegas)[len(omegas)//2]:.5f}"])
            w.writerow(["mean_omega", f"{sum(omegas)/len(omegas):.5f}"])
        w.writerow(["omega_positive_threshold", args.omega_positive])
        w.writerow(["omega_purifying_threshold", args.omega_purifying])
        w.writerow(["fubar_posterior_threshold", args.fubar_posterior])
        w.writerow(["meme_pvalue_threshold", args.meme_pvalue])
        w.writerow(["fubar_used", bool(f_rows)])
        w.writerow(["meme_used", bool(m_rows)])

    log(f"[selection] {len(rows)} sites: " +
        ", ".join(f"{k} {v}" for k, v in counts.items()))
    log(f"[selection] -> {args.out_sites}")
    return 0 if len(rows) >= args.min_sites else 1


if __name__ == "__main__":
    raise SystemExit(main())
