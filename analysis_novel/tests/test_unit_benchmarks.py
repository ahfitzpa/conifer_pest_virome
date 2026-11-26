"""Unit tests for Stage 17 benchmark & out-of-order validation scripts.

Run with:
    python3 -m unittest discover -s analysis_novel/tests -p "test_unit_*.py" -v
"""
from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ANALYSIS_DIR / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import benchmark_discovery_and_compute  # noqa: E402
import benchmark_snowball_gain  # noqa: E402
import benchmark_structural_triage  # noqa: E402


class TestStage17Benchmarks(unittest.TestCase):
    def test_structural_triage_and_out_of_order_split(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            val = tdir / "esm_rdrp_validation.tsv"
            val.write_text(
                "candidate_id\tcontig\trun_id\tconfidence\tconsensus_family\tnearest_identity_pct\tselection_reason\tmean_plddt\tmotif_c\trdrp_target\trdrp_tmscore\trdrp_prob\trdrp_evalue\trdrp_theader\trt_tmscore\tstructural_verdict\n"
                "c1\tSRR12261680|N1\tSRR12261680\thigh\tPolycipiviridae\t24.9\tdivergent_aa_id<40pct\t84.5\tGDD\t3nky_A\t0.88\t1.0\t1e-20\tFMDV 3Dpol RdRp\t0.62\tvalidated_rdrp_high\n"
                "c2\tSRR6765939|N2\tSRR6765939\tmedium\tUnresolved\t\tunresolved\t82.0\t\t6id1_A\t0.91\t1.0\t1e-25\tDEAH-box helicase Prp2\t0.0\thost_non_rdrp\n"
                "c3\tSRR6765939|N3\tSRR6765939\tmedium\tUnresolved\t\tunresolved\t78.0\tGDD\t1rtd_A\t0.85\t1.0\t1e-18\tHIV-1 reverse transcriptase\t0.92\trt_or_non_rdrp_flag\n"
                "c4\tSRR9990001|N4\tSRR9990001\tmedium\tUnresolved\t\tunresolved\t80.0\t\t6id1_A\t0.90\t1.0\t1e-22\tDEAH-box helicase DHX15\t0.0\thost_non_rdrp\n",
                encoding="utf-8",
            )
            out_dir = tdir / "02_structural_triage"
            rc = benchmark_structural_triage.main([
                "--validation-tsv", str(val),
                "--out-dir", str(out_dir),
            ])
            self.assertEqual(rc, 0)
            self.assertTrue((out_dir / "ablation_summary.tsv").is_file())
            self.assertTrue((out_dir / "threshold_sweep.tsv").is_file())
            self.assertTrue((out_dir / "out_of_order_validation.tsv").is_file())
            svg = (out_dir / "roc_pr_curves.svg").read_text(encoding="utf-8")
            self.assertTrue(svg.startswith("<svg "))

            with (out_dir / "out_of_order_validation.tsv").open("r", encoding="utf-8") as fh:
                rows = {r["cohort"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(int(rows["coleoptera"]["n_viral_rdrp"]), 1)
            self.assertEqual(int(rows["out_of_order_locusta"]["n_total_domains"]), 1)
            self.assertEqual(int(rows["out_of_order_locusta"]["n_viral_rdrp"]), 0)

    def test_snowball_gain_benchmark_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            hits = tdir / "hits.tsv"
            hits.write_text(
                "qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore\n"
                "SRR13742307|NODE_1_length_12390\tcand_0018\t88.9\t512\t12390\t512\t1e-200\t950\n"
                "SRR6323862|NODE_5_length_6550\tcand_0017\t100.0\t604\t6550\t604\t0.0\t1200\n",
                encoding="utf-8",
            )
            man = tdir / "manifest.csv"
            man.write_text(
                "run_id,scientificname,bioproject,group\n"
                "SRR12261680,Ips confusus,PRJNA647034,seed\n"
                "SRR13742307,Ips typographus,PRJNA699282,Snowball_RdRp\n"
                "SRR6323862,Hylobius abietis,PRJNA419715,Snowball_RdRp\n",
                encoding="utf-8",
            )
            out_dir = tdir / "03_snowball_gain"
            rc = benchmark_snowball_gain.main([
                "--hits-tsv", str(hits),
                "--manifest", str(man),
                "--out-dir", str(out_dir),
            ])
            self.assertEqual(rc, 0)
            with (out_dir / "snowball_gain_summary.tsv").open("r", encoding="utf-8") as fh:
                metrics = {r["metric"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(int(metrics["seed_lineages_independently_recovered_ge90pct"]), 1)
            self.assertEqual(int(metrics["seed_lineages_sister_50_90pct"]), 1)
            self.assertEqual(int(metrics["complete_or_near_complete_contigs_ge9kb_added"]), 1)

    def test_discovery_and_compute_benchmark_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            val = tdir / "val.tsv"
            val.write_text(
                "candidate_id\tconsensus_family\tnearest_identity_pct\tselection_reason\trdrp_target\tstructural_verdict\n"
                "c1\tPolycipiviridae\t24.9\tdivergent_aa_id<40pct\t3nky_A\tvalidated_rdrp_high\n"
                "c2\tUnresolved\t\tunresolved\t6id1_A\thost_non_rdrp\n",
                encoding="utf-8",
            )
            bdir = tdir / "benchmarks" / "esm_rdrp_fold"
            bdir.mkdir(parents=True)
            (bdir / "run.tsv").write_text(
                "s\th:m:s\tmax_rss\tcpu_time\n"
                "120.0\t0:02:00\t4096.0\t240.0\n",
                encoding="utf-8",
            )
            out_dir = tdir / "04_discovery_and_compute"
            rc = benchmark_discovery_and_compute.main([
                "--esm-validation", str(val),
                "--benchmark-dir", str(tdir / "benchmarks"),
                "--out-dir", str(out_dir),
            ])
            self.assertEqual(rc, 0)
            with (out_dir / "compute_cost_summary.tsv").open("r", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["rule_name"], "esm_rdrp_fold")
            self.assertEqual(rows[0]["peak_rss_gb"], "4.00")


if __name__ == "__main__":
    unittest.main()
