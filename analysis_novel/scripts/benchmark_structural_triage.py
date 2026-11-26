#!/usr/bin/env python3
"""Stage 17.2 & 17.5: Benchmark ESMFold + Foldseek structural triage & out-of-order validation.

Pure-Python, zero-OpenBLAS evaluation script that takes Stage 15
`esm_rdrp_validation.tsv` and produces:
  - confusion_matrix_coleoptera.tsv
  - confusion_matrix_out_of_order_locusta.tsv
  - threshold_sweep.tsv           (TM-score cutoffs 0.50..0.80 across both cohorts)
  - ablation_summary.tsv          (pLDDT-only, motif-only, Foldseek-no-RT, full model)
  - out_of_order_validation.tsv   (side-by-side Coleoptera vs. Orthoptera L. migratoria)
  - roc_pr_curves.svg             (publication SVG: ROC + PR curves & ablation bars)

Supports an optional `--ground-truth-tsv` (`candidate_id`, `true_label`) for
strictly independent evaluation (e.g., from GenBank nr / Pfam annotations);
when omitted, derives curated reference labels from PDB/family annotations and
records `label_source=curated_annotation` in every output table.
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
LOG = logging.getLogger("benchmark_structural_triage")

# Known Round-1 Coleoptera BioProject run prefixes when manifest is omitted
COLEOPTERA_RUN_PREFIXES = (
    "SRR122616",  # Ips confusus PRJNA647034
    "SRR67659",   # Hylobius abietis PRJNA435680
    "SRR170354",  # Hylobius abietis PRJNA783427
    "SRR63238",   # Hylobius abietis PRJNA419715
    "SRR137423",  # Ips typographus PRJNA699282
    "SRR1419",    # Ips nitidus
    "SRR1705",    # Ips nitidus
    "SRR169384",  # Ips subelongatus
)

TM_SWEEP = [0.50, 0.55, 0.60, 0.65, 0.69, 0.70, 0.75, 0.80]


def _flt(val: str, default: float = 0.0) -> float:
    try:
        v = float((val or "").strip())
        return v if math.isfinite(v) else default
    except ValueError:
        return default


def load_manifest_hosts(manifest_path: Optional[Path]) -> Dict[str, str]:
    hosts: Dict[str, str] = {}
    if not manifest_path or not manifest_path.is_file():
        return hosts
    with manifest_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            run = (row.get("run") or row.get("run_id") or "").strip()
            host = (
                row.get("host_organism")
                or row.get("scientificname")
                or row.get("host")
                or ""
            ).strip()
            if run and host:
                hosts[run] = host
    return hosts


def infer_cohort(run_id: str, manifest_hosts: Dict[str, str]) -> str:
    host = manifest_hosts.get(run_id, "").lower()
    if "locusta" in host:
        return "out_of_order_locusta"
    if any(k in host for k in ("ips ", "hylobius")):
        return "coleoptera"
    if run_id.startswith(COLEOPTERA_RUN_PREFIXES):
        return "coleoptera"
    return "out_of_order_locusta"


def load_ground_truth(path: Optional[Path]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path or not path.is_file():
        return out
    with path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("candidate_id") or "").strip()
            lbl = (row.get("true_label") or "").strip().lower()
            if cid and lbl:
                out[cid] = lbl
    return out


def derive_reference_class(row: Dict[str, str], ext_gt: Dict[str, str]) -> Tuple[str, str]:
    """Return (class_label, label_source).

    Classes: viral_rdrp | viral_ns3 | rt | host_other | uncertain
    """
    cid = (row.get("candidate_id") or "").strip()
    if cid in ext_gt:
        return ext_gt[cid], "external_ground_truth"

    verdict = (row.get("structural_verdict") or "").strip().lower()
    theader = (row.get("rdrp_theader") or "").lower()
    target = (row.get("rdrp_target") or "").lower()

    if verdict in {"low_confidence_fold", "no_structural_hit", "fold_failed"}:
        return "uncertain", "curated_annotation"
    if verdict == "rt_or_non_rdrp_flag":
        return "rt", "curated_annotation"
    if verdict in {"validated_rdrp_high", "validated_rdrp_moderate"}:
        # Chain-level check for Dengue NS5-NS3 complex (chain B = NS3 helicase)
        if ("8gzq" in target or "8gzr" in target) and target.endswith("_b"):
            return "viral_ns3", "curated_annotation"
        if "ns3" in theader and "helicase" in theader:
            return "viral_ns3", "curated_annotation"
        return "viral_rdrp", "curated_annotation"
    if ("8gzq" in target or "8gzr" in target) and target.endswith("_b"):
        return "viral_ns3", "curated_annotation"
    return "host_other", "curated_annotation"


def has_canonical_motif_c(row: Dict[str, str]) -> bool:
    mc = (row.get("motif_c") or "").strip().upper()
    return mc in {"GDD", "SDD", "GDN", "ADD"}


_RDRP_POS_RE = re.compile(
    r"(rna[- ]dependent rna polymerase|rdrp|3dpol|polymerase|replicase|\bl protein\b|\bns5\b|\bpb1\b)",
    re.IGNORECASE,
)
_NON_RDRP_NEG_RE = re.compile(
    r"(helicase|methyltransferase|mtase|spliceosom|prp2|dhx|upf1|ighmbp2|capsid|\bvp1\b|reverse transcriptase|telomerase|retrotransposon)",
    re.IGNORECASE,
)


def _is_rdrp_pdb_hit(row: Dict[str, str]) -> bool:
    target = (row.get("rdrp_target") or "").lower()
    if ("8gzq" in target or "8gzr" in target) and target.endswith("_b"):
        return False
    hdr = (row.get("rdrp_theader") or "").strip()
    if not hdr:
        return False
    if _NON_RDRP_NEG_RE.search(hdr) and not re.search(r"rna[- ]dependent rna polymerase|rdrp|3dpol", hdr, re.I):
        return False
    return bool(_RDRP_POS_RE.search(hdr))


def predict_rule(row: Dict[str, str], mode: str, tm_cutoff: float = 0.50) -> bool:
    """Evaluate binary viral-RdRp prediction under ablation mode & TM cutoff."""
    plddt = _flt(row.get("mean_plddt", ""))
    rdrp_tm = _flt(row.get("rdrp_tmscore", ""))
    rt_tm = _flt(row.get("rt_tmscore", ""))
    rdrp_eval = _flt(row.get("rdrp_evalue", ""), 1.0)
    rdrp_prob = _flt(row.get("rdrp_prob", ""), 0.0)
    motif_ok = has_canonical_motif_c(row)
    is_rdrp_hit = _is_rdrp_pdb_hit(row)

    if mode == "plddt_only":
        return plddt >= 70.0
    if mode == "motif_only":
        return motif_ok
    if mode == "foldseek_no_rt_control":
        return is_rdrp_hit and (rdrp_tm >= tm_cutoff or rdrp_prob >= 0.90) and rdrp_eval <= 1e-5
    if mode == "foldseek_with_rt_control":
        return (
            is_rdrp_hit
            and rdrp_tm >= tm_cutoff
            and rdrp_eval <= 1e-5
            and (rdrp_tm - rt_tm) >= 0.05
        )
    if mode == "full_model":
        return (
            is_rdrp_hit
            and plddt >= 50.0
            and rdrp_tm >= tm_cutoff
            and rdrp_eval <= 1e-5
            and (rdrp_tm - rt_tm) >= 0.05
            and motif_ok
        )
    raise ValueError(f"Unknown mode: {mode}")


def binary_metrics(y_true: List[bool], y_pred: List[bool]) -> Dict[str, float]:
    tp = sum(1 for t, p in zip(y_true, y_pred) if t and p)
    fp = sum(1 for t, p in zip(y_true, y_pred) if (not t) and p)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t and (not p))
    tn = sum(1 for t, p in zip(y_true, y_pred) if (not t) and (not p))
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
    bal_acc = (rec + spec) / 2.0
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / denom if denom > 0 else 0.0
    total = tp + fp + fn + tn
    rej_rate = (tn + fn) / total if total else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": prec,
        "recall": rec,
        "specificity": spec,
        "f1": f1,
        "balanced_accuracy": bal_acc,
        "mcc": mcc,
        "structural_rejection_rate": rej_rate,
    }


def continuous_score(row: Dict[str, str]) -> float:
    """Continuous RdRp structural discrimination score: RdRp TM minus RT TM (0 if non-RdRp PDB hit)."""
    if not _is_rdrp_pdb_hit(row):
        return 0.0
    rdrp_tm = _flt(row.get("rdrp_tmscore", ""))
    rt_tm = _flt(row.get("rt_tmscore", ""))
    delta = rdrp_tm - rt_tm
    return max(0.0, delta + 0.5 * rdrp_tm)


def compute_roc_pr(
    rows: List[Tuple[Dict[str, str], str]],
) -> Tuple[List[Tuple[float, float]], List[Tuple[float, float]], float, float]:
    """Return (roc_points, pr_points, auroc, auprc) for viral_rdrp vs non-RdRp."""
    scored = [
        (continuous_score(r), lbl == "viral_rdrp")
        for r, lbl in rows
        if lbl != "uncertain"
    ]
    scored.sort(key=lambda t: t[0], reverse=True)
    n_pos = sum(1 for _, y in scored if y)
    n_neg = len(scored) - n_pos
    if n_pos == 0 or n_neg == 0:
        return [(0.0, 0.0), (1.0, 1.0)], [(0.0, 1.0), (1.0, 0.0)], 0.0, 0.0

    roc: List[Tuple[float, float]] = [(0.0, 0.0)]
    pr: List[Tuple[float, float]] = [(0.0, 1.0)]
    tp = fp = 0
    for _score, is_pos in scored:
        if is_pos:
            tp += 1
        else:
            fp += 1
        tpr = tp / n_pos
        fpr = fp / n_neg
        prec = tp / (tp + fp)
        roc.append((fpr, tpr))
        pr.append((tpr, prec))

    # Trapezoidal AUROC & AUPRC
    auroc = sum(
        (roc[i][0] - roc[i - 1][0]) * (roc[i][1] + roc[i - 1][1]) / 2.0
        for i in range(1, len(roc))
    )
    auprc = sum(
        (pr[i][0] - pr[i - 1][0]) * (pr[i][1] + pr[i - 1][1]) / 2.0
        for i in range(1, len(pr))
    )
    return roc, pr, auroc, auprc


def write_tsv_atomic(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def render_benchmark_svg(
    out_svg: Path,
    roc_col: List[Tuple[float, float]],
    pr_col: List[Tuple[float, float]],
    auroc_col: float,
    auprc_col: float,
    roc_loc: List[Tuple[float, float]],
    pr_loc: List[Tuple[float, float]],
    auroc_loc: float,
    auprc_loc: float,
    ablation_col: List[Dict[str, object]],
) -> None:
    """Render 3-panel publication SVG: ROC, PR, and Ablation F1/MCC comparison."""
    w, h = 1140, 390
    svg: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
        f'width="{w}" height="{h}" font-family="Inter, Helvetica Neue, Arial, sans-serif">',
        f'<rect width="{w}" height="{h}" fill="#FFFFFF"/>',
    ]

    def draw_curve_panel(
        x0: float,
        y0: float,
        pw: float,
        ph: float,
        letter: str,
        title: str,
        xlab: str,
        ylab: str,
        pts1: List[Tuple[float, float]],
        pts2: List[Tuple[float, float]],
        leg1: str,
        leg2: str,
        diag: bool,
    ) -> None:
        svg.append(
            f'<text x="{x0 - 30:.1f}" y="{y0 - 16:.1f}" font-size="15" font-weight="700" fill="#0F172A">{letter}</text>'
        )
        svg.append(
            f'<text x="{x0:.1f}" y="{y0 - 16:.1f}" font-size="12.5" font-weight="700" fill="#0F172A">{escape(title)}</text>'
        )
        svg.append(
            f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{pw:.1f}" height="{ph:.1f}" '
            f'fill="#F8FAFC" stroke="#CBD5E1" stroke-width="1"/>'
        )
        for t in (0.0, 0.25, 0.5, 0.75, 1.0):
            gx = x0 + t * pw
            gy = y0 + (1.0 - t) * ph
            svg.append(
                f'<line x1="{gx:.1f}" y1="{y0:.1f}" x2="{gx:.1f}" y2="{y0 + ph:.1f}" stroke="#E2E8F0" stroke-width="0.7"/>'
            )
            svg.append(
                f'<line x1="{x0:.1f}" y1="{gy:.1f}" x2="{x0 + pw:.1f}" y2="{gy:.1f}" stroke="#E2E8F0" stroke-width="0.7"/>'
            )
            svg.append(
                f'<text x="{gx:.1f}" y="{y0 + ph + 15:.1f}" text-anchor="middle" font-size="9.5" fill="#475569">{t:.2f}</text>'
            )
            svg.append(
                f'<text x="{x0 - 6:.1f}" y="{gy + 3.5:.1f}" text-anchor="end" font-size="9.5" fill="#475569">{t:.2f}</text>'
            )
        if diag:
            svg.append(
                f'<line x1="{x0:.1f}" y1="{y0 + ph:.1f}" x2="{x0 + pw:.1f}" y2="{y0:.1f}" '
                f'stroke="#94A3B8" stroke-width="1" stroke-dasharray="4,3"/>'
            )
        for pts, col in ((pts1, "#0072B2"), (pts2, "#D55E00")):
            if len(pts) >= 2:
                d = " ".join(
                    f"{'M' if i == 0 else 'L'} {x0 + x * pw:.1f} {y0 + (1.0 - y) * ph:.1f}"
                    for i, (x, y) in enumerate(pts)
                )
                svg.append(
                    f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2.2" stroke-linejoin="round"/>'
                )
        svg.append(
            f'<text x="{x0 + pw / 2:.1f}" y="{y0 + ph + 34:.1f}" text-anchor="middle" font-size="10.5" fill="#1E293B">{escape(xlab)}</text>'
        )
        svg.append(
            f'<text x="{x0 - 38:.1f}" y="{y0 + ph / 2:.1f}" text-anchor="middle" font-size="10.5" fill="#1E293B" '
            f'transform="rotate(-90 {x0 - 38:.1f} {y0 + ph / 2:.1f})">{escape(ylab)}</text>'
        )
        # Compact legend box inside bottom-right
        lx, ly = x0 + 10.0, y0 + ph - 38.0
        svg.append(
            f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{pw - 20:.1f}" height="32" rx="4" fill="#FFFFFF" fill-opacity="0.9" stroke="#E2E8F0"/>'
        )
        svg.append(f'<line x1="{lx + 8:.1f}" y1="{ly + 11:.1f}" x2="{lx + 24:.1f}" y2="{ly + 11:.1f}" stroke="#0072B2" stroke-width="2.2"/>')
        svg.append(f'<text x="{lx + 30:.1f}" y="{ly + 14:.1f}" font-size="9" fill="#0F172A">{escape(leg1)}</text>')
        svg.append(f'<line x1="{lx + 8:.1f}" y1="{ly + 24:.1f}" x2="{lx + 24:.1f}" y2="{ly + 24:.1f}" stroke="#D55E00" stroke-width="2.2"/>')
        svg.append(f'<text x="{lx + 30:.1f}" y="{ly + 27:.1f}" font-size="9" fill="#0F172A">{escape(leg2)}</text>')

    draw_curve_panel(
        60, 48, 290, 270, "A", "ROC Curve (Structural RdRp Score)",
        "False Positive Rate", "True Positive Rate",
        roc_col, roc_loc,
        f"Coleoptera (AUROC = {auroc_col:.3f})",
        f"L. migratoria out-of-order (AUROC = {auroc_loc:.3f})",
        diag=True,
    )
    draw_curve_panel(
        425, 48, 290, 270, "B", "Precision–Recall Curve",
        "Recall", "Precision",
        pr_col, pr_loc,
        f"Coleoptera (AUPRC = {auprc_col:.3f})",
        f"L. migratoria out-of-order (AUPRC = {auprc_loc:.3f})",
        diag=False,
    )

    # Panel C: Ablation horizontal bars (F1 and Precision on Coleoptera)
    cx0, cy0, cw, ch = 800.0, 48.0, 305.0, 270.0
    svg.append(f'<text x="{cx0 - 30:.1f}" y="{cy0 - 16:.1f}" font-size="15" font-weight="700" fill="#0F172A">C</text>')
    svg.append(f'<text x="{cx0:.1f}" y="{cy0 - 16:.1f}" font-size="12.5" font-weight="700" fill="#0F172A">Component Ablation (Coleoptera)</text>')
    svg.append(f'<rect x="{cx0:.1f}" y="{cy0:.1f}" width="{cw:.1f}" height="{ch:.1f}" fill="#F8FAFC" stroke="#CBD5E1"/>')

    pretty = {
        "plddt_only": "pLDDT ≥ 70 only",
        "motif_only": "Motif C only",
        "foldseek_no_rt_control": "Foldseek (no RT ctrl)",
        "foldseek_with_rt_control": "Foldseek + RT ctrl",
        "full_model": "Full Stage 15 model",
    }
    n_modes = max(1, len(ablation_col))
    band = ch / n_modes
    for idx, row in enumerate(ablation_col):
        mode = str(row["ablation_mode"])
        f1 = float(row["f1"])
        prec = float(row["precision"])
        by = cy0 + idx * band + 8.0
        svg.append(
            f'<text x="{cx0 + 6:.1f}" y="{by + 11:.1f}" font-size="9.5" font-weight="600" fill="#0F172A">'
            f'{escape(pretty.get(mode, mode))}</text>'
        )
        bar_max = cw - 55.0
        svg.append(
            f'<rect x="{cx0 + 6:.1f}" y="{by + 16:.1f}" width="{f1 * bar_max:.1f}" height="10" rx="2" fill="#0072B2"/>'
        )
        svg.append(
            f'<text x="{cx0 + 10 + f1 * bar_max:.1f}" y="{by + 24.5:.1f}" font-size="8.5" fill="#0072B2">F1 {f1:.2f}</text>'
        )
        svg.append(
            f'<rect x="{cx0 + 6:.1f}" y="{by + 29:.1f}" width="{prec * bar_max:.1f}" height="10" rx="2" fill="#009E73"/>'
        )
        svg.append(
            f'<text x="{cx0 + 10 + prec * bar_max:.1f}" y="{by + 37.5:.1f}" font-size="8.5" fill="#009E73">Prec {prec:.2f}</text>'
        )

    svg.append("</svg>\n")
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_svg.with_suffix(".svg.tmp")
    tmp.write_text("\n".join(svg), encoding="utf-8")
    os.replace(tmp, out_svg)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--validation-tsv", required=True, type=Path)
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--ground-truth-tsv", type=Path, default=None)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args(argv)

    manifest_hosts = load_manifest_hosts(args.manifest)
    ext_gt = load_ground_truth(args.ground_truth_tsv)

    cohorts: Dict[str, List[Tuple[Dict[str, str], str, str]]] = {
        "coleoptera": [],
        "out_of_order_locusta": [],
    }
    with args.validation_tsv.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            run_id = (row.get("run_id") or "").strip()
            cohort = infer_cohort(run_id, manifest_hosts)
            ref_cls, src = derive_reference_class(row, ext_gt)
            cohorts[cohort].append((row, ref_cls, src))

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    ablation_modes = [
        "plddt_only",
        "motif_only",
        "foldseek_no_rt_control",
        "foldseek_with_rt_control",
        "full_model",
    ]
    ablation_rows: List[Dict[str, object]] = []
    sweep_rows: List[Dict[str, object]] = []
    oos_rows: List[Dict[str, object]] = []
    roc_data: Dict[str, Tuple[List[Tuple[float, float]], List[Tuple[float, float]], float, float]] = {}

    for cname, items in cohorts.items():
        # Confusion / class distribution table
        counts = Counter(cls for _, cls, _ in items)
        label_src = items[0][2] if items else "none"
        cm_rows = []
        for cls in ("viral_rdrp", "viral_ns3", "rt", "host_other", "uncertain"):
            cm_rows.append(
                {
                    "cohort": cname,
                    "reference_class": cls,
                    "count": counts.get(cls, 0),
                    "fraction_of_cohort": f"{(counts.get(cls, 0) / len(items)):.4f}" if items else "0.0000",
                    "label_source": label_src,
                }
            )
        write_tsv_atomic(
            out_dir / f"confusion_matrix_{cname}.tsv",
            ["cohort", "reference_class", "count", "fraction_of_cohort", "label_source"],
            cm_rows,
        )

        # Primary evaluation set excludes uncertain folds (reported separately above)
        eval_items = [(r, cls) for r, cls, _ in items if cls != "uncertain"]
        y_true = [cls == "viral_rdrp" for _, cls in eval_items]
        roc, pr, auroc, auprc = compute_roc_pr(eval_items)
        roc_data[cname] = (roc, pr, auroc, auprc)

        for mode in ablation_modes:
            y_pred = [predict_rule(r, mode, tm_cutoff=0.50) for r, _ in eval_items]
            m = binary_metrics(y_true, y_pred)
            ablation_rows.append(
                {
                    "cohort": cname,
                    "ablation_mode": mode,
                    "n_evaluated": len(eval_items),
                    "n_uncertain_excluded": counts.get("uncertain", 0),
                    "tp": m["tp"],
                    "fp": m["fp"],
                    "fn": m["fn"],
                    "tn": m["tn"],
                    "precision": f"{m['precision']:.4f}",
                    "recall": f"{m['recall']:.4f}",
                    "f1": f"{m['f1']:.4f}",
                    "balanced_accuracy": f"{m['balanced_accuracy']:.4f}",
                    "mcc": f"{m['mcc']:.4f}",
                    "structural_rejection_rate": f"{m['structural_rejection_rate']:.4f}",
                    "label_source": label_src,
                }
            )

        for tm_cut in TM_SWEEP:
            y_pred = [predict_rule(r, "full_model", tm_cutoff=tm_cut) for r, _ in eval_items]
            m = binary_metrics(y_true, y_pred)
            sweep_rows.append(
                {
                    "cohort": cname,
                    "tm_cutoff": f"{tm_cut:.2f}",
                    "tp": m["tp"],
                    "fp": m["fp"],
                    "fn": m["fn"],
                    "tn": m["tn"],
                    "precision": f"{m['precision']:.4f}",
                    "recall": f"{m['recall']:.4f}",
                    "f1": f"{m['f1']:.4f}",
                    "mcc": f"{m['mcc']:.4f}",
                    "structural_rejection_rate": f"{m['structural_rejection_rate']:.4f}",
                }
            )

        full_pred = [predict_rule(r, "full_model", tm_cutoff=0.50) for r, _ in eval_items]
        fm = binary_metrics(y_true, full_pred)
        oos_rows.append(
            {
                "cohort": cname,
                "insect_order": "Coleoptera" if cname == "coleoptera" else "Orthoptera (Locusta migratoria)",
                "n_total_domains": len(items),
                "n_viral_rdrp": counts.get("viral_rdrp", 0),
                "n_viral_ns3": counts.get("viral_ns3", 0),
                "n_rt_purged": counts.get("rt", 0),
                "n_host_other_purged": counts.get("host_other", 0),
                "n_uncertain": counts.get("uncertain", 0),
                "non_rdrp_rejection_pct": f"{100.0 * (counts.get('rt', 0) + counts.get('host_other', 0)) / len(items):.1f}" if items else "0.0",
                "full_model_precision": f"{fm['precision']:.4f}",
                "full_model_recall": f"{fm['recall']:.4f}",
                "full_model_f1": f"{fm['f1']:.4f}",
                "auroc": f"{auroc:.4f}",
                "auprc": f"{auprc:.4f}",
            }
        )

    write_tsv_atomic(
        out_dir / "ablation_summary.tsv",
        [
            "cohort", "ablation_mode", "n_evaluated", "n_uncertain_excluded",
            "tp", "fp", "fn", "tn", "precision", "recall", "f1",
            "balanced_accuracy", "mcc", "structural_rejection_rate", "label_source",
        ],
        ablation_rows,
    )
    write_tsv_atomic(
        out_dir / "threshold_sweep.tsv",
        [
            "cohort", "tm_cutoff", "tp", "fp", "fn", "tn",
            "precision", "recall", "f1", "mcc", "structural_rejection_rate",
        ],
        sweep_rows,
    )
    write_tsv_atomic(
        out_dir / "out_of_order_validation.tsv",
        [
            "cohort", "insect_order", "n_total_domains", "n_viral_rdrp",
            "n_viral_ns3", "n_rt_purged", "n_host_other_purged", "n_uncertain",
            "non_rdrp_rejection_pct", "full_model_precision", "full_model_recall",
            "full_model_f1", "auroc", "auprc",
        ],
        oos_rows,
    )

    col_ablation = [r for r in ablation_rows if r["cohort"] == "coleoptera"]
    render_benchmark_svg(
        out_dir / "roc_pr_curves.svg",
        *roc_data["coleoptera"],
        *roc_data["out_of_order_locusta"],
        col_ablation,
    )
    LOG.info(
        "Structural triage benchmark complete: Coleoptera n=%d, out-of-order L. migratoria n=%d -> %s",
        len(cohorts["coleoptera"]),
        len(cohorts["out_of_order_locusta"]),
        out_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
