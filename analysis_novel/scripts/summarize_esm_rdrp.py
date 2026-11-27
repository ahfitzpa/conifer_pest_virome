#!/usr/bin/env python3
"""Integrate ESMFold metrics with Foldseek RdRp & RT negative-control searches.

Produces:
  1. `--out-validation`: per-candidate structural validation table combining
     candidate metadata, ESMFold pLDDT/pTM, canonical palm motifs, best RdRp
     Foldseek structural hit, best Reverse Transcriptase (RT) negative-control
     hit, and a conservative structural verdict.
  2. `--out-summary`: aggregate counts by verdict tier and consensus family.

Foldseek searches are expected to be run with:
  --format-output query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore,lddt,prob,theader
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path

__version__ = "1.0.0"

LOG = logging.getLogger("summarize_esm_rdrp")

FOLDSEEK_COLS = [
    "query",
    "target",
    "fident",
    "alnlen",
    "evalue",
    "bits",
    "qtmscore",
    "ttmscore",
    "alntmscore",
    "lddt",
    "prob",
    "theader",
]

RT_PATTERNS = re.compile(
    r"(reverse\s+transcriptase|telomerase|retrotransposon|retrovirus|"
    r"group\s+ii\s+intron|maturase|hiv-?1\s+rt|m-?mlv|rna-directed\s+dna\s+polymerase)",
    re.IGNORECASE,
)

NON_RDRP_VIRAL_PATTERNS = re.compile(
    r"(capsid|coat\s+protein|envelope|glycoprotein|nucleocapsid|matrix\s+protein|"
    r"protease|helicase|methyltransferase|spliceosom|deah-box|dead-box|ribosom)",
    re.IGNORECASE,
)

RDRP_PATTERNS = re.compile(
    r"(rna-dependent\s+rna\s+polymerase|rna-directed\s+rna\s+polymerase|"
    r"\brdrp\b|replicase|polyprotein|\bnsp12\b|\bns5\b|\bns5b\b|\bpb1\b|\b3dpol\b|"
    r"\bl\s+protein\b|large\s+structural\s+protein|viral\s+polymerase)",
    re.IGNORECASE,
)

EXPLICIT_RDRP_PATTERNS = re.compile(
    r"(rna-dependent\s+rna\s+polymerase|rna-directed\s+rna\s+polymerase|\brdrp\b|\b3dpol\b)",
    re.IGNORECASE,
)


def atomic_write_text(dest: Path, text: str) -> None:
    """Write text atomically via a sibling temporary file and os.replace()."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=dest.parent, prefix=f".{dest.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp_path, dest)
    except BaseException:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def safe_float(val: str | None, default: float = 0.0) -> float:
    try:
        return float(val) if val not in (None, "") else default
    except ValueError:
        return default


def normalize_query_id(raw_query: str) -> str:
    """Strip Foldseek's .pdb / chain suffix (e.g., cand_0001__foo.pdb -> cand_0001__foo)."""
    q = raw_query.strip()
    if q.endswith(".pdb"):
        q = q[:-4]
    if q.endswith(".cif"):
        q = q[:-4]
    return q


def classify_target_header(header: str) -> str:
    """Classify a Foldseek target description as 'rt', 'rdrp', or 'other'."""
    if RT_PATTERNS.search(header):
        return "rt"
    if NON_RDRP_VIRAL_PATTERNS.search(header) and not EXPLICIT_RDRP_PATTERNS.search(header):
        return "other"
    if RDRP_PATTERNS.search(header):
        return "rdrp"
    return "other"


def load_foldseek_hits(
    path: Path | None,
    assume_all_rdrp: bool = False,
    assume_all_rt: bool = False,
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Parse Foldseek tabular output into (best_rdrp_hit, best_rt_hit) per query.

    When searching a general PDB database (`assume_all_rdrp=False`), target
    headers are classified into RdRp vs. RT hits so a single PDB Foldseek search
    still provides an internal RT negative control. When searching a dedicated
    RdRp DB (`assume_all_rdrp=True`) or dedicated RT DB (`assume_all_rt=True`),
    hits are assigned directly.
    """
    best_rdrp: dict[str, dict[str, str]] = {}
    best_rt: dict[str, dict[str, str]] = {}
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return best_rdrp, best_rt

    def _better(new: dict[str, str], old: dict[str, str] | None) -> bool:
        if old is None:
            return True
        new_key = (safe_float(new["alntmscore"]), safe_float(new["prob"]), safe_float(new["bits"]))
        old_key = (safe_float(old["alntmscore"]), safe_float(old["prob"]), safe_float(old["bits"]))
        return new_key > old_key

    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            # Skip optional header line if present
            if parts[0] == "query" and len(parts) >= 2 and parts[1] == "target":
                continue
            while len(parts) < len(FOLDSEEK_COLS):
                parts.append("")
            rec = dict(zip(FOLDSEEK_COLS, parts[: len(FOLDSEEK_COLS)]))
            qid = normalize_query_id(rec["query"])
            rec["query_id"] = qid

            if assume_all_rt:
                if _better(rec, best_rt.get(qid)):
                    best_rt[qid] = rec
                continue

            role = "rdrp" if assume_all_rdrp else classify_target_header(rec["theader"])
            if role == "rt":
                if _better(rec, best_rt.get(qid)):
                    best_rt[qid] = rec
            elif role == "rdrp":
                if _better(rec, best_rdrp.get(qid)):
                    best_rdrp[qid] = rec
            else:
                # For unclassified PDB headers, keep only if no explicit RdRp hit exists
                # and mark header_class='other' so downstream logic can inspect it.
                rec["header_class"] = "other"
                if qid not in best_rdrp:
                    best_rdrp[qid] = rec

    return best_rdrp, best_rt


def assign_verdict(
    mean_plddt: float,
    fold_status: str,
    rdrp_hit: dict[str, str] | None,
    rt_hit: dict[str, str] | None,
    high_plddt: float = 70.0,
    min_plddt: float = 50.0,
    high_tm: float = 0.50,
    mod_tm: float = 0.40,
    rt_margin: float = 0.05,
) -> tuple[str, str]:
    """Return (structural_verdict, rationale) for a single candidate."""
    if fold_status not in ("ok", "cached"):
        return "folding_failed", f"ESMFold status={fold_status}"
    if mean_plddt < min_plddt:
        return "low_confidence_fold", f"mean pLDDT {mean_plddt:.1f} < {min_plddt:g}"

    rdrp_tm = safe_float(rdrp_hit.get("alntmscore")) if rdrp_hit else 0.0
    rdrp_prob = safe_float(rdrp_hit.get("prob")) if rdrp_hit else 0.0
    rdrp_eval = safe_float(rdrp_hit.get("evalue"), 1.0) if rdrp_hit else 1.0
    rdrp_is_explicit = (rdrp_hit is not None) and (rdrp_hit.get("header_class", "rdrp") == "rdrp")

    rt_tm = safe_float(rt_hit.get("alntmscore")) if rt_hit else 0.0
    rt_eval = safe_float(rt_hit.get("evalue"), 1.0) if rt_hit else 1.0

    # Negative-control RT flag: RT structural match beats RdRp match
    if rt_hit and rt_eval <= 1e-3 and (rt_tm >= mod_tm) and (rt_tm > rdrp_tm + rt_margin or not rdrp_is_explicit):
        return (
            "rt_or_non_rdrp_flag",
            f"RT/retroelement structural hit (TM={rt_tm:.3f}, E={rt_eval:.2e}) exceeds RdRp (TM={rdrp_tm:.3f})",
        )

    if rdrp_hit and rdrp_is_explicit:
        if (
            mean_plddt >= high_plddt
            and (rdrp_tm >= high_tm or rdrp_prob >= 0.90)
            and rdrp_eval <= 1e-5
            and (not rt_hit or rdrp_tm >= rt_tm + rt_margin)
        ):
            return (
                "validated_rdrp_high",
                f"pLDDT={mean_plddt:.1f}, Foldseek RdRp TM={rdrp_tm:.3f}, prob={rdrp_prob:.2f}, E={rdrp_eval:.2e}",
            )
        if (rdrp_tm >= mod_tm or rdrp_prob >= 0.50) and rdrp_eval <= 1e-3 and rdrp_tm >= rt_tm:
            return (
                "validated_rdrp_moderate",
                f"pLDDT={mean_plddt:.1f}, Foldseek RdRp TM={rdrp_tm:.3f}, prob={rdrp_prob:.2f}, E={rdrp_eval:.2e}",
            )

    if rdrp_hit and not rdrp_is_explicit and rdrp_tm >= mod_tm and rdrp_eval <= 1e-3:
        return (
            "non_rdrp_structural_hit",
            f"Top PDB hit is non-RdRp ({rdrp_hit.get('target', '')}: TM={rdrp_tm:.3f}, E={rdrp_eval:.2e})",
        )

    return "no_structural_hit", f"pLDDT={mean_plddt:.1f}, no passing Foldseek RdRp hit (best TM={rdrp_tm:.3f})"


def load_tsv_by_key(path: Path, key_col: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if not path.is_file() or path.stat().st_size == 0:
        return out
    with path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            k = (row.get(key_col) or "").strip()
            if k:
                out[k] = row
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selected-tsv", required=True, type=Path,
                    help="Candidate manifest from select_divergent_rdrp.py")
    ap.add_argument("--metrics-tsv", required=True, type=Path,
                    help="ESMFold metrics TSV from run_esmfold_rdrp.py")
    ap.add_argument("--foldseek-rdrp", required=True, type=Path,
                    help="Foldseek easy-search TSV against RdRp/PDB structure DB")
    ap.add_argument("--foldseek-rt", default=None, type=Path,
                    help="Optional Foldseek easy-search TSV against dedicated RT negative-control DB")
    ap.add_argument("--dedicated-rdrp-db", action="store_true",
                    help="Treat all hits in --foldseek-rdrp as RdRp (e.g., Rider/Viro3D RdRp-only DB)")
    ap.add_argument("--out-validation", required=True, type=Path,
                    help="Output per-candidate structural validation TSV")
    ap.add_argument("--out-summary", required=True, type=Path,
                    help="Output summary counts TSV")
    ap.add_argument("--log-level", default="INFO",
                    choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    selected = load_tsv_by_key(args.selected_tsv, "candidate_id")
    metrics = load_tsv_by_key(args.metrics_tsv, "candidate_id")
    rdrp_hits, pdb_rt_hits = load_foldseek_hits(
        args.foldseek_rdrp,
        assume_all_rdrp=args.dedicated_rdrp_db,
    )
    _, dedicated_rt_hits = load_foldseek_hits(args.foldseek_rt, assume_all_rt=True)

    # Merge RT negative-control hits (dedicated RT DB takes precedence if better TM-score)
    rt_hits = dict(pdb_rt_hits)
    for qid, hit in dedicated_rt_hits.items():
        prev = rt_hits.get(qid)
        if prev is None or safe_float(hit["alntmscore"]) > safe_float(prev["alntmscore"]):
            rt_hits[qid] = hit

    out_cols = [
        "candidate_id", "contig", "run_id", "confidence", "consensus_family",
        "nearest_ictv_ref", "nearest_identity_pct", "selection_reason",
        "domain_aa_len", "folded_aa_len", "was_trimmed",
        "fold_status", "mean_plddt", "median_plddt", "frac_plddt_ge_70", "ptm",
        "motif_a", "motif_b", "motif_c", "palm_motifs_ordered",
        "rdrp_target", "rdrp_tmscore", "rdrp_prob", "rdrp_evalue", "rdrp_fident", "rdrp_theader",
        "rt_target", "rt_tmscore", "rt_evalue", "rt_theader",
        "structural_verdict", "verdict_rationale",
    ]

    rows: list[dict[str, str]] = []
    verdict_counts: Counter[str] = Counter()

    for cand_id, sel in selected.items():
        met = metrics.get(cand_id, {})
        rh = rdrp_hits.get(cand_id)
        th = rt_hits.get(cand_id)
        mean_plddt = safe_float(met.get("mean_plddt"), 0.0)
        fold_status = met.get("status", "missing")

        verdict, rationale = assign_verdict(
            mean_plddt=mean_plddt,
            fold_status=fold_status,
            rdrp_hit=rh,
            rt_hit=th,
        )
        verdict_counts[verdict] += 1

        rows.append({
            "candidate_id": cand_id,
            "contig": sel.get("contig", ""),
            "run_id": sel.get("run_id", ""),
            "confidence": sel.get("confidence", ""),
            "consensus_family": sel.get("consensus_family", ""),
            "nearest_ictv_ref": sel.get("nearest_ictv_ref", ""),
            "nearest_identity_pct": sel.get("nearest_identity_pct", ""),
            "selection_reason": sel.get("selection_reason", ""),
            "domain_aa_len": sel.get("domain_aa_len", ""),
            "folded_aa_len": sel.get("folded_aa_len", ""),
            "was_trimmed": sel.get("was_trimmed", ""),
            "fold_status": fold_status,
            "mean_plddt": met.get("mean_plddt", ""),
            "median_plddt": met.get("median_plddt", ""),
            "frac_plddt_ge_70": met.get("frac_plddt_ge_70", ""),
            "ptm": met.get("ptm", ""),
            "motif_a": met.get("motif_a", ""),
            "motif_b": met.get("motif_b", ""),
            "motif_c": met.get("motif_c", ""),
            "palm_motifs_ordered": met.get("palm_motifs_ordered", ""),
            "rdrp_target": rh.get("target", "") if rh else "",
            "rdrp_tmscore": rh.get("alntmscore", "") if rh else "",
            "rdrp_prob": rh.get("prob", "") if rh else "",
            "rdrp_evalue": rh.get("evalue", "") if rh else "",
            "rdrp_fident": rh.get("fident", "") if rh else "",
            "rdrp_theader": rh.get("theader", "") if rh else "",
            "rt_target": th.get("target", "") if th else "",
            "rt_tmscore": th.get("alntmscore", "") if th else "",
            "rt_evalue": th.get("evalue", "") if th else "",
            "rt_theader": th.get("theader", "") if th else "",
            "structural_verdict": verdict,
            "verdict_rationale": rationale,
        })

    val_lines = ["\t".join(out_cols)]
    for r in rows:
        val_lines.append("\t".join(str(r.get(c, "")) for c in out_cols))
    atomic_write_text(args.out_validation, "\n".join(val_lines) + "\n")

    sum_lines = ["metric\tcount"]
    sum_lines.append(f"total_selected\t{len(rows)}")
    for verdict in (
        "validated_rdrp_high",
        "validated_rdrp_moderate",
        "rt_or_non_rdrp_flag",
        "non_rdrp_structural_hit",
        "low_confidence_fold",
        "no_structural_hit",
        "folding_failed",
    ):
        sum_lines.append(f"{verdict}\t{verdict_counts.get(verdict, 0)}")
    atomic_write_text(args.out_summary, "\n".join(sum_lines) + "\n")

    LOG.info(
        "Summarized %d candidate(s): %d high, %d moderate, %d RT-flagged -> %s",
        len(rows),
        verdict_counts.get("validated_rdrp_high", 0),
        verdict_counts.get("validated_rdrp_moderate", 0),
        verdict_counts.get("rt_or_non_rdrp_flag", 0),
        args.out_validation,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
