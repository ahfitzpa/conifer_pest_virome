#!/usr/bin/env python3
"""Select medium-to-high confidence, deeply divergent RdRp domains for ESMFold.

Gating criteria (all configurable via CLI flags):
  1. Contig confidence in `contig_confidence.tsv` is in `--keep-confidence`
     (default: `high,medium`) AND `has_rdrp == 1` (or an extracted RdRp domain
     exists for that contig).
  2. Extracted RdRp amino-acid domain length >= `--min-aa-len` (default: 200 aa).
  3. At least one divergence condition holds:
       (a) Nearest ICTV/nr amino-acid identity < `--max-identity-pct`
           (default: 40.0%) in demarcation / candidate summary table(s); OR
       (b) Family-unresolved (`consensus_family` is blank, `unclassified`,
           `unresolved`, `unknown`, or `NA`) when `--include-unresolved` is set.

When a contig has multiple RdRpCATCH domain hits, the highest-scoring (or
longest) domain is retained by default so ESMFold does not redundantly fold
overlapping hits from the same locus. Domains longer than `--max-aa-len`
(default: 800 aa) are center-trimmed to fit GPU memory without OOM.
"""

from __future__ import annotations

import argparse
import csv
import glob
import logging
import os
import re
import sys
import tempfile
from pathlib import Path

__version__ = "1.0.0"

LOG = logging.getLogger("select_divergent_rdrp")

UNRESOLVED_TOKENS = {
    "",
    "unclassified",
    "unresolved",
    "unknown",
    "none",
    "na",
    "n/a",
    "other",
    "unassigned",
}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


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


def sanitize_sequence(seq: str) -> str:
    """Strip stop codons/whitespace and replace non-canonical residues with X."""
    cleaned = re.sub(r"[\s*]+", "", seq.upper())
    return "".join(ch if ch in VALID_AA else "X" for ch in cleaned)


def trim_centered(seq: str, max_len: int) -> tuple[str, bool]:
    """Center-trim a sequence exceeding max_len to preserve the RdRp palm core."""
    if max_len <= 0 or len(seq) <= max_len:
        return seq, False
    start = (len(seq) - max_len) // 2
    return seq[start : start + max_len], True


def canonical_contig_key(raw: str) -> str:
    """Extract a normalized (run_id, node_short) or contig token for joining."""
    s = raw.strip()
    if "|" in s:
        parts = s.split("|")
        run_id, contig = parts[0], parts[1]
        m = re.match(r"NODE_(\d+)_", contig)
        if m:
            return f"{run_id}_n{m.group(1)}"
        return f"{run_id}|{contig}"
    return s


def parse_rdrp_header(header: str) -> dict[str, str]:
    """Parse headers written by scripts/extract_rdrp_domains.py.

    Format:
      {sample}|{contig_id}|frame={frame}|aa={a}-{b}|profile={p}|db={d}|evalue={e}|norm_bitscore_profile={nbp}
    """
    token = header.strip().split()[0]
    parts = token.split("|")
    info: dict[str, str] = {
        "raw_header": token,
        "sample": parts[0] if len(parts) >= 1 else "",
        "contig_id": parts[1] if len(parts) >= 2 else parts[0],
        "full_contig": f"{parts[0]}|{parts[1]}" if len(parts) >= 2 else parts[0],
        "frame": "",
        "aa_range": "",
        "profile": "",
        "db": "",
        "evalue": "",
        "norm_bitscore_profile": "",
    }
    for field in parts[2:]:
        if "=" in field:
            k, v = field.split("=", 1)
            if k == "aa":
                info["aa_range"] = v
            elif k in info:
                info[k] = v
    return info


def iter_fasta(path: Path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def load_demarcation_identities(paths: list[str]) -> dict[str, tuple[float, str, str]]:
    """Map candidate/contig keys -> (nearest_identity_pct, nearest_ref, family)."""
    out: dict[str, tuple[float, str, str]] = {}
    expanded: list[Path] = []
    for pat in paths:
        matches = sorted(glob.glob(pat))
        if matches:
            expanded.extend(Path(m) for m in matches)
        elif Path(pat).is_file():
            expanded.append(Path(pat))

    for p in expanded:
        if not p.is_file() or p.stat().st_size == 0:
            continue
        with p.open("r", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                id_raw = (
                    row.get("nearest_identity_pct")
                    or row.get("best_identity_pct")
                    or row.get("pident")
                    or ""
                ).strip()
                if not id_raw:
                    continue
                try:
                    ident = float(id_raw)
                except ValueError:
                    continue
                ref = (row.get("nearest_ictv_ref") or row.get("nearest_ref") or "").strip()
                fam = (row.get("family") or "").strip()

                keys: set[str] = set()
                for col in ("candidate", "contig", "contig_id", "query"):
                    val = (row.get(col) or "").strip()
                    if val:
                        keys.add(val)
                        keys.add(canonical_contig_key(val))
                run_id = (row.get("run_id") or "").strip()
                contig_id = (row.get("contig_id") or "").strip()
                if run_id and contig_id:
                    full = f"{run_id}|{contig_id}"
                    keys.add(full)
                    keys.add(canonical_contig_key(full))

                for k in keys:
                    # Keep the highest known identity if duplicate entries exist
                    # (conservative: do not mislabel a close match as divergent).
                    if k not in out or ident > out[k][0]:
                        out[k] = (ident, ref, fam)
    return out


def load_confidence_table(path: Path) -> dict[str, dict[str, str]]:
    """Load contig_confidence.tsv indexed by both full contig and canonical key."""
    index: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            if not contig:
                continue
            index[contig] = row
            index[canonical_contig_key(contig)] = row
            if "|" in contig:
                index[contig.split("|", 1)[1]] = row
    return index


def is_unresolved_family(family: str) -> bool:
    return family.strip().lower() in UNRESOLVED_TOKENS


def safe_float(val: str, default: float = -1.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def select_candidates(
    confidence_tsv: Path,
    rdrp_faa_patterns: list[str],
    demarcation_patterns: list[str],
    keep_confidence: set[str],
    min_aa_len: int,
    max_aa_len: int,
    max_identity_pct: float,
    include_unresolved: bool,
    dedup_per_contig: bool,
) -> tuple[list[dict[str, str]], list[tuple[str, str]]]:
    """Filter RdRp domains and return (metadata_rows, fasta_records)."""
    conf_idx = load_confidence_table(confidence_tsv)
    demarc_idx = load_demarcation_identities(demarcation_patterns)
    LOG.info(
        "Loaded %d confidence index keys and %d demarcation identity keys",
        len(conf_idx),
        len(demarc_idx),
    )

    faa_files: list[Path] = []
    seen_paths: set[Path] = set()
    for pat in rdrp_faa_patterns:
        matches = sorted(glob.glob(pat))
        candidates = [Path(m) for m in matches] if matches else ([Path(pat)] if Path(pat).is_file() else [])
        for p in candidates:
            rp = p.resolve()
            if rp not in seen_paths and p.is_file() and p.stat().st_size > 0:
                seen_paths.add(rp)
                faa_files.append(p)

    LOG.info("Scanning %d RdRp domain FASTA file(s)", len(faa_files))

    # Collect passing hits; optionally deduplicate per contig keeping best score/length.
    per_contig: dict[str, tuple[tuple[float, int], dict[str, str], str]] = {}
    multi_records: list[tuple[dict[str, str], str]] = []

    for faa in faa_files:
        for header, raw_seq in iter_fasta(faa):
            info = parse_rdrp_header(header)
            seq = sanitize_sequence(raw_seq)
            aa_len = len(seq)
            if aa_len < min_aa_len:
                continue

            full_contig = info["full_contig"]
            canon = canonical_contig_key(full_contig)
            conf_row = (
                conf_idx.get(full_contig)
                or conf_idx.get(canon)
                or conf_idx.get(info["contig_id"])
            )
            if conf_row is None:
                continue

            conf_level = (conf_row.get("confidence") or "").strip().lower()
            if conf_level not in keep_confidence:
                continue

            cons_family = (conf_row.get("consensus_family") or "").strip()
            unresolved = is_unresolved_family(cons_family)

            demarc_hit = (
                demarc_idx.get(full_contig)
                or demarc_idx.get(canon)
                or demarc_idx.get(info["contig_id"])
            )
            ident_pct: float | None = demarc_hit[0] if demarc_hit else None
            nearest_ref = demarc_hit[1] if demarc_hit else ""

            reasons: list[str] = []
            if ident_pct is not None and ident_pct < max_identity_pct:
                reasons.append(f"divergent_aa_id<{max_identity_pct:g}pct")
            if include_unresolved and unresolved:
                reasons.append("family_unresolved_rdrp")

            if not reasons:
                continue

            folded_seq, was_trimmed = trim_centered(seq, max_aa_len)
            score_tuple = (safe_float(info["norm_bitscore_profile"], -1.0), aa_len)

            meta = {
                "contig": full_contig,
                "canonical_id": canon,
                "run_id": (conf_row.get("run_id") or info["sample"]).strip(),
                "confidence": conf_level,
                "consensus_family": cons_family if cons_family else "Unresolved",
                "nearest_ictv_ref": nearest_ref,
                "nearest_identity_pct": f"{ident_pct:.2f}" if ident_pct is not None else "",
                "selection_reason": ";".join(reasons),
                "rdrp_profile": info["profile"],
                "rdrp_db": info["db"],
                "rdrp_evalue": info["evalue"],
                "norm_bitscore_profile": info["norm_bitscore_profile"],
                "aa_range": info["aa_range"],
                "domain_aa_len": str(aa_len),
                "folded_aa_len": str(len(folded_seq)),
                "was_trimmed": "1" if was_trimmed else "0",
                "source_faa": str(faa),
            }

            if dedup_per_contig:
                prev = per_contig.get(full_contig)
                if prev is None or score_tuple > prev[0]:
                    per_contig[full_contig] = (score_tuple, meta, folded_seq)
            else:
                multi_records.append((meta, folded_seq))

    if dedup_per_contig:
        ordered = sorted(
            [(m, s) for _, m, s in per_contig.values()],
            key=lambda x: (
                0 if x[0]["confidence"] == "high" else 1,
                x[0]["consensus_family"],
                x[0]["contig"],
            ),
        )
    else:
        ordered = multi_records

    rows: list[dict[str, str]] = []
    fasta_records: list[tuple[str, str]] = []
    for idx, (meta, folded_seq) in enumerate(ordered, start=1):
        safe_canon = re.sub(r"[^A-Za-z0-9_.-]+", "_", meta["canonical_id"])[:60]
        cand_id = f"cand_{idx:04d}__{safe_canon}"
        row = {"candidate_id": cand_id, **meta}
        rows.append(row)
        fasta_records.append((cand_id, folded_seq))

    return rows, fasta_records


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--confidence-tsv", required=True, type=Path,
                    help="Path to results/reporting/contig_confidence.tsv")
    ap.add_argument("--rdrp-faa", action="append", required=True,
                    help="Glob or path to *_rdrp_domains.faa (repeatable)")
    ap.add_argument("--demarcation-tsv", action="append", default=[],
                    help="Optional glob or path to demarcation summary TSV(s)")
    ap.add_argument("--keep-confidence", default="high,medium",
                    help="Comma-separated confidence tiers to keep (default: high,medium)")
    ap.add_argument("--min-aa-len", type=int, default=200,
                    help="Minimum RdRp domain amino-acid length (default: 200)")
    ap.add_argument("--max-aa-len", type=int, default=800,
                    help="Maximum aa length passed to ESMFold; longer domains are center-trimmed (default: 800)")
    ap.add_argument("--max-identity-pct", type=float, default=40.0,
                    help="Select candidates with nearest ICTV/nr aa identity strictly below this %% (default: 40.0)")
    ap.add_argument("--no-unresolved", action="store_true",
                    help="Exclude family-unresolved RdRp-positive contigs")
    ap.add_argument("--keep-all-hits", action="store_true",
                    help="Keep all RdRpCATCH hits per contig instead of deduplicating to the best hit")
    ap.add_argument("--out-faa", required=True, type=Path,
                    help="Output FASTA of selected RdRp domains for ESMFold")
    ap.add_argument("--out-tsv", required=True, type=Path,
                    help="Output TSV manifest of selected candidates")
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

    keep_conf = {x.strip().lower() for x in args.keep_confidence.split(",") if x.strip()}
    rows, fasta_records = select_candidates(
        confidence_tsv=args.confidence_tsv,
        rdrp_faa_patterns=args.rdrp_faa,
        demarcation_patterns=args.demarcation_tsv,
        keep_confidence=keep_conf,
        min_aa_len=args.min_aa_len,
        max_aa_len=args.max_aa_len,
        max_identity_pct=args.max_identity_pct,
        include_unresolved=not args.no_unresolved,
        dedup_per_contig=not args.keep_all_hits,
    )

    cols = [
        "candidate_id", "contig", "canonical_id", "run_id", "confidence",
        "consensus_family", "nearest_ictv_ref", "nearest_identity_pct",
        "selection_reason", "rdrp_profile", "rdrp_db", "rdrp_evalue",
        "norm_bitscore_profile", "aa_range", "domain_aa_len",
        "folded_aa_len", "was_trimmed", "source_faa",
    ]
    tsv_lines = ["\t".join(cols)]
    for r in rows:
        tsv_lines.append("\t".join(str(r.get(c, "")) for c in cols))
    atomic_write_text(args.out_tsv, "\n".join(tsv_lines) + "\n")

    faa_lines: list[str] = []
    for cand_id, seq in fasta_records:
        faa_lines.append(f">{cand_id}")
        for i in range(0, len(seq), 60):
            faa_lines.append(seq[i : i + 60])
    atomic_write_text(args.out_faa, ("\n".join(faa_lines) + "\n") if faa_lines else "")

    LOG.info(
        "Selected %d divergent RdRp candidate(s) -> %s, %s",
        len(rows),
        args.out_faa,
        args.out_tsv,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
