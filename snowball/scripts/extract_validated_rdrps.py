#!/usr/bin/env python3
"""Extract structurally validated novel RdRp domains from Stage 15 ESM-RdRp outputs.

Reads `esm_rdrp_validation.tsv` and `selected_divergent_rdrp.faa`, filters to
requested structural verdict tiers (default: validated_rdrp_high,
validated_rdrp_moderate), and writes:
  - `--out-faa`: clean protein FASTA ready for PalmDB/DIAMOND or palmID
  - `--out-tsv`: compact provenance table (candidate_id, contig, run_id,
                 consensus_family, mean_plddt, rdrp_tmscore, structural_verdict)
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import tempfile
from pathlib import Path

__version__ = "1.0.0"
LOG = logging.getLogger("extract_validated_rdrps")


def atomic_write_text(dest: Path, text: str) -> None:
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


def iter_fasta(path: Path):
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header.split()[0], "".join(chunks)


def extract_validated(
    validation_tsv: Path,
    selected_faa: Path,
    allowed_verdicts: set[str],
) -> tuple[list[dict[str, str]], list[tuple[str, str]]]:
    kept_meta: dict[str, dict[str, str]] = {}
    with validation_tsv.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            verdict = (row.get("structural_verdict") or "").strip()
            if verdict in allowed_verdicts:
                cid = (row.get("candidate_id") or "").strip()
                if cid:
                    kept_meta[cid] = row

    seqs = dict(iter_fasta(selected_faa)) if selected_faa.is_file() else {}

    out_rows: list[dict[str, str]] = []
    out_fasta: list[tuple[str, str]] = []
    for cid, row in kept_meta.items():
        seq = seqs.get(cid, "")
        if not seq:
            LOG.warning("Validated candidate %s not found in %s; skipping", cid, selected_faa)
            continue
        out_rows.append(row)
        out_fasta.append((cid, seq))

    return out_rows, out_fasta


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--validation-tsv", required=True, type=Path)
    ap.add_argument("--selected-faa", required=True, type=Path)
    ap.add_argument("--verdicts", default="validated_rdrp_high,validated_rdrp_moderate",
                    help="Comma-separated structural_verdict values to keep")
    ap.add_argument("--out-faa", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
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

    allowed = {v.strip() for v in args.verdicts.split(",") if v.strip()}
    rows, records = extract_validated(args.validation_tsv, args.selected_faa, allowed)

    cols = [
        "candidate_id", "contig", "run_id", "confidence", "consensus_family",
        "nearest_ictv_ref", "nearest_identity_pct", "mean_plddt",
        "rdrp_target", "rdrp_tmscore", "structural_verdict",
    ]
    tsv_lines = ["\t".join(cols)]
    for r in rows:
        tsv_lines.append("\t".join(str(r.get(c, "")) for c in cols))
    atomic_write_text(args.out_tsv, "\n".join(tsv_lines) + "\n")

    faa_lines: list[str] = []
    for cid, seq in records:
        faa_lines.append(f">{cid}")
        for i in range(0, len(seq), 60):
            faa_lines.append(seq[i : i + 60])
    atomic_write_text(args.out_faa, ("\n".join(faa_lines) + "\n") if faa_lines else "")

    LOG.info("Extracted %d structurally validated RdRp(s) -> %s", len(records), args.out_faa)
    return 0


if __name__ == "__main__":
    sys.exit(main())
