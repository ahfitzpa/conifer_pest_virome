"""Unit tests for the RdRp -> PalmDB/Serratus snowball side workflow.

Run with:
    python3 -m unittest discover -s analysis_novel/tests -p "test_unit_*.py" -v
"""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

SNOWBALL_SCRIPTS = Path(__file__).resolve().parents[2] / "snowball" / "scripts"
if str(SNOWBALL_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SNOWBALL_SCRIPTS))

import extract_validated_rdrps  # noqa: E402
import serratus_snowball  # noqa: E402


class TestSnowballWorkflow(unittest.TestCase):
    def test_extract_validated_rdrps_filters_verdicts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            val = tdir / "val.tsv"
            val.write_text(
                "candidate_id\tcontig\trun_id\tconfidence\tconsensus_family\tnearest_ictv_ref\tnearest_identity_pct\tmean_plddt\trdrp_target\trdrp_tmscore\tstructural_verdict\n"
                "cand_0001\tSRR1|N1\tSRR1\thigh\tPolycipiviridae\tREF1\t24.9\t84.5\t3nky\t0.881\tvalidated_rdrp_high\n"
                "cand_0002\tSRR1|N2\tSRR1\tmedium\tUnresolved\t\t\t72.0\t1rtd\t0.610\trt_or_non_rdrp_flag\n",
                encoding="utf-8",
            )
            faa = tdir / "sel.faa"
            faa.write_text(">cand_0001\nMKTAYIAK\n>cand_0002\nGGGGGGGG\n", encoding="utf-8")

            rows, recs = extract_validated_rdrps.extract_validated(
                val, faa, {"validated_rdrp_high", "validated_rdrp_moderate"}
            )
            self.assertEqual(len(rows), 1)
            self.assertEqual(recs, [("cand_0001", "MKTAYIAK")])

    def test_parse_palmdb_hits_and_sotu_extraction(self) -> None:
        self.assertEqual(serratus_snowball.extract_sotu_id("u12345_RdRp_palm"), "u12345")
        self.assertEqual(serratus_snowball.extract_sotu_id("u98765"), "u98765")

        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(
                "cand_0001\tu12345\t62.5\t110\t1e-25\t180\n"
                "cand_0001\tu99999\t31.0\t90\t1e-10\t95\n"
                "cand_0002\tu12345\t55.0\t105\t1e-18\t150\n",
                encoding="utf-8",
            )
            mapping = serratus_snowball.parse_palmdb_hits(hits, min_pident=40.0, max_evalue=1e-5)
            self.assertEqual(set(mapping.keys()), {"u12345"})
            self.assertEqual(mapping["u12345"], {"cand_0001", "cand_0002"})

    def test_select_diverse_runs_excludes_existing_and_caps_bioproject(self) -> None:
        cands = [
            {"run_id": "SRR_OLD", "bioproject": "PRJ1", "n_matched_sotus": 5, "max_coverage": 200.0, "max_pident": 90.0},
            {"run_id": "SRR_A1",  "bioproject": "PRJ2", "n_matched_sotus": 3, "max_coverage": 150.0, "max_pident": 85.0},
            {"run_id": "SRR_A2",  "bioproject": "PRJ2", "n_matched_sotus": 3, "max_coverage": 120.0, "max_pident": 82.0},
            {"run_id": "SRR_A3",  "bioproject": "PRJ2", "n_matched_sotus": 2, "max_coverage": 110.0, "max_pident": 80.0},
            {"run_id": "SRR_B1",  "bioproject": "PRJ3", "n_matched_sotus": 2, "max_coverage": 95.0,  "max_pident": 78.0},
        ]
        sel = serratus_snowball.select_diverse_runs(
            cands,
            existing_runs={"SRR_OLD"},
            max_runs_per_bioproject=2,
            max_total_runs=3,
        )
        self.assertEqual([r["run_id"] for r in sel], ["SRR_A1", "SRR_A2", "SRR_B1"])

        manifest_rows = serratus_snowball.build_phase2_manifest_rows(sel)
        self.assertEqual(len(manifest_rows), 3)
        self.assertEqual(list(manifest_rows[0].keys()), serratus_snowball.MANIFEST_COLS)
        self.assertEqual(manifest_rows[0]["source"], "sra")
        self.assertEqual(manifest_rows[0]["run_id"], "SRR_A1")


    def test_open_host_specificity_and_spotcheck_selection(self) -> None:
        import open_host_specificity

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            inp = tdir / "open_host.tsv"
            inp.write_text(
                "run_id\tbioproject\tscientific_name\tmatched_sotus\tmax_coverage\tsize_mb\n"
                "SRR1\tPRJNA1\tIps nitidus\tu993423\t80\t900\n"
                "SRR2\tPRJNA2\tIps typographus\tu993423\t60\t800\n"
                "SRR3\tPRJNA3\tIps duplicatus\tu1029776\t55\t700\n"
                "SRR4\tPRJNA4\tVespula vulgaris\tu653762\t120\t600\n"
                "SRR5\tPRJNA5\tGraminella nigrifrons\tu653762\t90\t1200\n"
                "SRR6\tPRJNA6\tMus musculus\tu761589\t30\t500\n"
                "SRR7\tPRJNA7\tHomo sapiens\tu761589\t25\t450\n",
                encoding="utf-8",
            )
            spec = tdir / "spec.tsv"
            man = tdir / "spotcheck.csv"
            rc = open_host_specificity.main([
                "--open-host-tsv", str(inp),
                "--out-specificity-tsv", str(spec),
                "--out-spotcheck-manifest", str(man),
            ])
            self.assertEqual(rc, 0)
            spec_txt = spec.read_text(encoding="utf-8")
            self.assertIn("strict_beetle_specialist", spec_txt)
            self.assertIn("pan_arthropod_generalist", spec_txt)
            self.assertIn("broad_low_specificity", spec_txt)
            man_txt = man.read_text(encoding="utf-8")
            self.assertIn("SRR4", man_txt)  # Vespula vulgaris generalist
            self.assertIn("SRR3", man_txt)  # Ips duplicatus sister specialist
            self.assertIn("OpenHost_SpotCheck", man_txt)
            with man.open("r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh)
                self.assertEqual(reader.fieldnames, serratus_snowball.MANIFEST_COLS)
                rows = list(reader)
            self.assertGreaterEqual(len(rows), 2)
            self.assertEqual(rows[0]["source"], "sra")
            self.assertEqual(rows[0]["sample_id"], rows[0]["run_id"])


if __name__ == "__main__":
    unittest.main()
