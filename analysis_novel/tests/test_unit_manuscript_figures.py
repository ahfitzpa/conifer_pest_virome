"""Unit tests for Stage 16 reproducible manuscript figure generators (Figs. 1-3).

Run with:
    python3 -m unittest discover -s analysis_novel/tests -p "test_unit_*.py" -v
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ANALYSIS_DIR / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import plot_fig1_workflow_and_map  # noqa: E402
import plot_fig2_esm_rdrp  # noqa: E402
import plot_fig3_snowball_matrix  # noqa: E402


class TestManuscriptFigureGenerators(unittest.TestCase):
    def test_fig1_workflow_and_map_cli(self) -> None:
        host_geo = ANALYSIS_DIR / "configs" / "host_geography.tsv"
        self.assertTrue(host_geo.is_file())
        with tempfile.TemporaryDirectory() as tmp:
            out_svg = Path(tmp) / "fig1.svg"
            rc = plot_fig1_workflow_and_map.main([
                "--host-geo", str(host_geo),
                "--out-svg", str(out_svg),
            ])
            self.assertEqual(rc, 0)
            text = out_svg.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("<svg "))
            self.assertIn("Ips confusus", text)
            self.assertIn("Ips nitidus", text)
            self.assertIn("Hylobius abietis", text)

    def test_fig2_esm_rdrp_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            val = tdir / "val.tsv"
            val.write_text(
                "candidate_id\tcontig\trun_id\tconfidence\tconsensus_family\tfolded_aa_len\tmean_plddt\trdrp_target\trdrp_tmscore\trt_tmscore\tstructural_verdict\n"
                "cand_1\tSRR12261680|N2\tSRR12261680\thigh\tPolycipiviridae\t512\t84.5\t3nky_A\t0.881\t0.650\tvalidated_rdrp_high\n"
                "cand_2\tSRR6765939|N9\tSRR6765939\tmedium\tUnresolved\t300\t75.0\t1rtd_A\t0.890\t0.910\trt_or_non_rdrp_flag\n",
                encoding="utf-8",
            )
            out_svg = tdir / "fig2.svg"
            rc = plot_fig2_esm_rdrp.main([
                "--validation-tsv", str(val),
                "--out-validation-svg", str(out_svg),
            ])
            self.assertEqual(rc, 0)
            text = out_svg.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("<svg "))
            self.assertIn("n = 2", text)

    def test_fig3_snowball_matrix_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            hits = tdir / "hits.tsv"
            hits.write_text(
                "qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore\n"
                "SRR13742307|NODE_1_length_12390\tcand_0018__SRR12261680_n2\t88.9\t512\t12390\t512\t1e-200\t950\n"
                "SRR6323862|NODE_5_length_6550\tcand_0017__SRR6765939_n252\t100.0\t604\t6550\t604\t0.0\t1200\n",
                encoding="utf-8",
            )
            man = tdir / "manifest.csv"
            man.write_text(
                "run_id,scientificname,bioproject\n"
                "SRR13742307,Ips typographus,PRJNA699282\n"
                "SRR6323862,Hylobius abietis,PRJNA419715\n",
                encoding="utf-8",
            )
            out_svg = tdir / "fig3.svg"
            rc = plot_fig3_snowball_matrix.main([
                "--hits-tsv", str(hits),
                "--manifest", str(man),
                "--out-svg", str(out_svg),
            ])
            self.assertEqual(rc, 0)
            text = out_svg.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("<svg "))
            self.assertIn("Ips typographus", text)
            self.assertIn("Hylobius abietis", text)


if __name__ == "__main__":
    unittest.main()
