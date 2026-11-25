#!/usr/bin/env python3
"""Run local GPU ESMFold on selected divergent RdRp protein domains.

Uses HuggingFace `transformers` (`facebook/esmfold_v1`, pure PyTorch) so it runs
cleanly on Linux/SLURM GPU nodes without requiring legacy `openfold` CUDA C++
extensions. Existing non-empty PDB files in `--out-pdb-dir` are reused so long
HPC jobs are safely resumable after preemption or walltime limits.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import statistics
import sys
import tempfile
import time
from pathlib import Path

__version__ = "1.0.0"

LOG = logging.getLogger("run_esmfold_rdrp")

# Canonical RdRp catalytic palm motifs (Wolf et al. 2018; Olendraite et al. 2023)
MOTIF_A_RE = re.compile(r"D[A-Z]{3,5}D")
MOTIF_B_RE = re.compile(r"[ST]G[A-Z]{2,4}[TN][A-Z]{2,4}[NS]")
MOTIF_C_RE = re.compile(r"(GDD|SDD|GDN|IDD|ADD|VDD)")


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


def iter_fasta(path: Path):
    """Yield (record_id, sequence) tuples from a FASTA file."""
    rec_id: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rec_id is not None:
                    yield rec_id, "".join(chunks)
                rec_id = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line)
    if rec_id is not None:
        yield rec_id, "".join(chunks)


def scan_palm_motifs(seq: str) -> dict[str, str]:
    """Scan an RdRp amino-acid sequence for canonical palm motifs A, B, and C."""
    s = seq.upper()
    ma = MOTIF_A_RE.search(s)
    mb = MOTIF_B_RE.search(s)
    mc = MOTIF_C_RE.search(s)
    # Check canonical N->C ordering A < B < C when multiple motifs are present
    ordered = False
    if ma and mb and mc:
        ordered = ma.start() < mb.start() < mc.start()
    elif ma and mc:
        ordered = ma.start() < mc.start()
    return {
        "motif_a": ma.group(0) if ma else "",
        "motif_b": mb.group(0) if mb else "",
        "motif_c": mc.group(1) if mc else "",
        "palm_motifs_ordered": "1" if ordered else "0",
    }


def parse_pdb_plddt(pdb_text: str) -> dict[str, float]:
    """Extract per-residue CA pLDDT values from a PDB string (0-100 scale)."""
    scores: list[float] = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  "):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            bfac = float(line[60:66].strip())
            scores.append(bfac)
        except ValueError:
            continue

    if not scores:
        return {
            "n_residues": 0.0,
            "mean_plddt": 0.0,
            "median_plddt": 0.0,
            "frac_plddt_ge_70": 0.0,
            "frac_plddt_ge_50": 0.0,
        }

    # Normalize 0-1 scale (some transformers versions) to 0-100 scale
    if max(scores) <= 1.0:
        scores = [x * 100.0 for x in scores]

    n = len(scores)
    return {
        "n_residues": float(n),
        "mean_plddt": statistics.fmean(scores),
        "median_plddt": statistics.median(scores),
        "frac_plddt_ge_70": sum(1 for x in scores if x >= 70.0) / n,
        "frac_plddt_ge_50": sum(1 for x in scores if x >= 50.0) / n,
    }


def fold_sequences(
    records: list[tuple[str, str]],
    out_pdb_dir: Path,
    model_name: str,
    cache_dir: str | None,
    chunk_size: int,
    use_fp16: bool,
    allow_cpu: bool,
) -> list[dict[str, str]]:
    """Fold protein sequences with local ESMFold and return per-sequence metrics."""
    out_pdb_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    if not records:
        LOG.warning("No candidate sequences to fold; writing empty metrics table.")
        return rows

    # Check which records need folding before loading the multi-GB model weights
    to_fold: list[tuple[str, str]] = []
    cached_rows: dict[str, dict[str, str]] = {}
    for cand_id, seq in records:
        pdb_path = out_pdb_dir / f"{cand_id}.pdb"
        motifs = scan_palm_motifs(seq)
        if pdb_path.is_file() and pdb_path.stat().st_size > 100:
            stats = parse_pdb_plddt(pdb_path.read_text(encoding="utf-8"))
            cached_rows[cand_id] = {
                "candidate_id": cand_id,
                "aa_len": str(len(seq)),
                "status": "cached",
                "pdb_path": str(pdb_path),
                "mean_plddt": f"{stats['mean_plddt']:.2f}",
                "median_plddt": f"{stats['median_plddt']:.2f}",
                "frac_plddt_ge_70": f"{stats['frac_plddt_ge_70']:.3f}",
                "frac_plddt_ge_50": f"{stats['frac_plddt_ge_50']:.3f}",
                "ptm": "",
                "runtime_sec": "0.0",
                **motifs,
            }
        else:
            to_fold.append((cand_id, seq))

    LOG.info(
        "Found %d total candidate(s): %d cached PDB(s), %d to fold",
        len(records),
        len(cached_rows),
        len(to_fold),
    )

    folded_rows: dict[str, dict[str, str]] = {}
    if to_fold:
        import torch
        from transformers import AutoTokenizer, EsmForProteinFolding

        has_cuda = torch.cuda.is_available()
        if not has_cuda and not allow_cpu:
            raise RuntimeError(
                f"CUDA GPU is required for local ESMFold on Sonic GPU nodes, "
                f"but torch.cuda.is_available() is False "
                f"(torch={torch.__version__}, torch.version.cuda={torch.version.cuda}). "
                f"Submit to a GPU partition (e.g., --gres=gpu:1) or pass --allow-cpu for testing."
            )

        device = torch.device("cuda" if has_cuda else "cpu")
        if has_cuda:
            gpu_name = torch.cuda.get_device_name(0)
            vram_gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
            LOG.info("Using GPU: %s (%.1f GiB VRAM)", gpu_name, vram_gb)
        else:
            LOG.warning("Running ESMFold on CPU (--allow-cpu set); this will be slow.")

        LOG.info("Loading ESMFold model %s (cache_dir=%s)", model_name, cache_dir)
        tokenizer = AutoTokenizer.from_pretrained(model_name, cache_dir=cache_dir)
        model = EsmForProteinFolding.from_pretrained(model_name, cache_dir=cache_dir)
        model = model.to(device)
        model.eval()
        if chunk_size > 0 and hasattr(model, "trunk") and hasattr(model.trunk, "set_chunk_size"):
            model.trunk.set_chunk_size(chunk_size)
        if has_cuda and use_fp16 and hasattr(model, "esm"):
            model.esm = model.esm.half()

        for idx, (cand_id, seq) in enumerate(to_fold, start=1):
            t0 = time.monotonic()
            motifs = scan_palm_motifs(seq)
            pdb_path = out_pdb_dir / f"{cand_id}.pdb"
            LOG.info("[%d/%d] Folding %s (%d aa)", idx, len(to_fold), cand_id, len(seq))
            try:
                tokenized = tokenizer(
                    [seq],
                    return_tensors="pt",
                    add_special_tokens=False,
                )
                tokenized = {k: v.to(device) for k, v in tokenized.items()}
                with torch.no_grad():
                    output = model(**tokenized)
                pdb_list = model.output_to_pdb(output)
                pdb_text = pdb_list[0]
                atomic_write_text(pdb_path, pdb_text)

                stats = parse_pdb_plddt(pdb_text)
                ptm_val = ""
                if hasattr(output, "ptm") and output.ptm is not None:
                    ptm_val = f"{float(output.ptm.mean().item()):.4f}"
                elapsed = time.monotonic() - t0
                folded_rows[cand_id] = {
                    "candidate_id": cand_id,
                    "aa_len": str(len(seq)),
                    "status": "ok",
                    "pdb_path": str(pdb_path),
                    "mean_plddt": f"{stats['mean_plddt']:.2f}",
                    "median_plddt": f"{stats['median_plddt']:.2f}",
                    "frac_plddt_ge_70": f"{stats['frac_plddt_ge_70']:.3f}",
                    "frac_plddt_ge_50": f"{stats['frac_plddt_ge_50']:.3f}",
                    "ptm": ptm_val,
                    "runtime_sec": f"{elapsed:.1f}",
                    **motifs,
                }
            except RuntimeError as exc:
                elapsed = time.monotonic() - t0
                is_oom = "out of memory" in str(exc).lower()
                LOG.error(
                    "Failed folding %s (%d aa) after %.1fs: %s",
                    cand_id,
                    len(seq),
                    elapsed,
                    exc,
                )
                folded_rows[cand_id] = {
                    "candidate_id": cand_id,
                    "aa_len": str(len(seq)),
                    "status": "oom_error" if is_oom else "runtime_error",
                    "pdb_path": "",
                    "mean_plddt": "0.00",
                    "median_plddt": "0.00",
                    "frac_plddt_ge_70": "0.000",
                    "frac_plddt_ge_50": "0.000",
                    "ptm": "",
                    "runtime_sec": f"{elapsed:.1f}",
                    **motifs,
                }
            finally:
                if has_cuda:
                    torch.cuda.empty_cache()

    for cand_id, _ in records:
        if cand_id in cached_rows:
            rows.append(cached_rows[cand_id])
        elif cand_id in folded_rows:
            rows.append(folded_rows[cand_id])

    return rows


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-faa", required=True, type=Path,
                    help="Input FASTA of selected divergent RdRp domains")
    ap.add_argument("--out-pdb-dir", required=True, type=Path,
                    help="Directory to write predicted .pdb structure files")
    ap.add_argument("--out-metrics", required=True, type=Path,
                    help="Output TSV of per-structure ESMFold pLDDT/pTM and palm motif metrics")
    ap.add_argument("--model-name", default="facebook/esmfold_v1",
                    help="HuggingFace model ID or local path (default: facebook/esmfold_v1)")
    ap.add_argument("--cache-dir", default=None,
                    help="HuggingFace cache directory (recommended: scratch-backed dir on HPC)")
    ap.add_argument("--chunk-size", type=int, default=64,
                    help="ESMFold trunk axial attention chunk size to reduce VRAM (default: 64)")
    ap.add_argument("--no-fp16", action="store_true",
                    help="Keep ESM-2 language model backbone in fp32 instead of fp16")
    ap.add_argument("--allow-cpu", action="store_true",
                    help="Allow CPU fallback if no CUDA GPU is available (testing only)")
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

    records = list(iter_fasta(args.in_faa)) if args.in_faa.is_file() else []
    rows = fold_sequences(
        records=records,
        out_pdb_dir=args.out_pdb_dir,
        model_name=args.model_name,
        cache_dir=args.cache_dir,
        chunk_size=args.chunk_size,
        use_fp16=not args.no_fp16,
        allow_cpu=args.allow_cpu,
    )

    cols = [
        "candidate_id", "aa_len", "status", "pdb_path",
        "mean_plddt", "median_plddt", "frac_plddt_ge_70", "frac_plddt_ge_50",
        "ptm", "runtime_sec", "motif_a", "motif_b", "motif_c", "palm_motifs_ordered",
    ]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in cols))
    atomic_write_text(args.out_metrics, "\n".join(lines) + "\n")
    LOG.info("Wrote ESMFold metrics for %d candidate(s) -> %s", len(rows), args.out_metrics)
    return 0


if __name__ == "__main__":
    sys.exit(main())
