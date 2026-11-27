"""Unit tests for Stage 16 reproducible manuscript figure generators (Figs. 1-6).

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
import plot_phylo_manuscript  # noqa: E402


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

    def test_fig4_polycipiviridae_with_orfs_and_host_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            tree = tdir / "tree.newick"
            tree.write_text(
                "(((cand|SRR12261680|NODE_2_length_12300:0.12,"
                "cand|SRR13742307|NODE_1_length_12390:0.11)98.5/99:0.25,"
                "cand|SRR9999999|NODE_9_length_5000:0.40)92/95:0.30,"
                "(ICTV|NC_035455|Sinv3_virus:0.55,ICTV|NC_040651|Sinv5_virus:0.48)100/100:0.35);\n",
                encoding="utf-8",
            )
            meta = tdir / "meta.tsv"
            meta.write_text(
                "taxon_id\tgroup\tfamily\thost\taccession\tsource\n"
                "cand|SRR12261680|NODE_2_length_12300\tcandidate\tPolycipiviridae\tIps confusus\t\tassembly\n"
                "cand|SRR13742307|NODE_1_length_12390\tcandidate\tPolycipiviridae\tIps typographus\t\tassembly\n"
                "cand|SRR9999999|NODE_9_length_5000\tcandidate\tPolycipiviridae\tLocusta migratoria\t\tassembly\n"
                "ICTV|NC_035455|Sinv3_virus\tictv\tPolycipiviridae\t\tNC_035455\tICTV\n"
                "ICTV|NC_040651|Sinv5_virus\tictv\tPolycipiviridae\t\tNC_040651\tICTV\n",
                encoding="utf-8",
            )
            gm = tdir / "gm.tsv"
            gm.write_text(
                "family\tsample\tcontig_id\tcontig_length\tORF_start\tORF_end\tproduct\n"
                "Polycipiviridae\tSRR12261680\tNODE_2_length_12300\t12300\t500\t2200\tCapsid VP1\n"
                "Polycipiviridae\tSRR12261680\tNODE_2_length_12300\t12300\t6500\t11800\tRNA-dependent RNA polymerase\n",
                encoding="utf-8",
            )
            esm_val = tdir / "esm.tsv"
            esm_val.write_text(
                "candidate_id\tcontig\tstructural_verdict\trdrp_tmscore\n"
                "c1\tSRR12261680|NODE_2_length_12300\tvalidated_rdrp_high\t0.881\n",
                encoding="utf-8",
            )
            sb_hits = tdir / "sb.tsv"
            sb_hits.write_text(
                "qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore\n"
                "SRR13742307|NODE_1_length_12390\tc1\t88.9\t512\t12390\t512\t1e-200\t950\n",
                encoding="utf-8",
            )
            out_svg = tdir / "fig4.svg"
            rc = plot_phylo_manuscript.main([
                "--mode", "fig4",
                "--tree", str(tree),
                "--meta", str(meta),
                "--genome-map", str(gm),
                "--esm-validation", str(esm_val),
                "--snowball-hits", str(sb_hits),
                "--out-svg", str(out_svg),
            ])
            self.assertEqual(rc, 0)
            text = out_svg.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("<svg "))
            self.assertIn("I. confusus (NA)", text)
            self.assertIn("I. typographus (EU)", text)
            self.assertIn("Genome &amp; ORF Architecture", text)
            self.assertIn("RdRp 3D", text)
            self.assertIn("3D 0.88", text)
            self.assertIn("89%", text)
            self.assertNotIn("Locusta migratoria", text)
            self.assertNotIn("SRR9999999", text)

    def test_fig5_iflaviridae_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            tree = tdir / "tree.newick"
            tree.write_text(
                "((cand|SRR16938486|NODE_3_length_9800:0.15,cand|SRR6765939|NODE_4_length_9940:0.18)95/97:0.25,"
                "ICTV|NC_002066|Infectious_flacherie_virus:0.50);\n",
                encoding="utf-8",
            )
            meta = tdir / "meta.tsv"
            meta.write_text(
                "taxon_id\tgroup\tfamily\thost\taccession\tsource\n"
                "cand|SRR16938486|NODE_3_length_9800\tcandidate\tIflaviridae\tIps subelongatus\t\tassembly\n"
                "cand|SRR6765939|NODE_4_length_9940\tcandidate\tIflaviridae\tHylobius abietis\t\tassembly\n"
                "ICTV|NC_002066|Infectious_flacherie_virus\tictv\tIflaviridae\t\tNC_002066\tICTV\n",
                encoding="utf-8",
            )
            out_svg = tdir / "fig5.svg"
            rc = plot_phylo_manuscript.main([
                "--mode", "fig5",
                "--tree", str(tree),
                "--meta", str(meta),
                "--out-svg", str(out_svg),
            ])
            self.assertEqual(rc, 0)
            text = out_svg.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("<svg "))
            self.assertIn("I. subelongatus (EA)", text)
            self.assertIn("H. abietis (EU)", text)

    def test_fig6_multipanel_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            t1 = tdir / "amar.nwk"
            t2 = tdir / "phen.nwk"
            for tp in (t1, t2):
                tp.write_text(
                    "((cand|SRR6765940|NODE_1_length_3100:0.12,ICTV|NC_012345|Jingmen_tick_virus:0.30)99/100:0.2,"
                    "ICTV|NC_054321|Alongshan_virus:0.35);\n",
                    encoding="utf-8",
                )
            m1 = tdir / "meta.tsv"
            m1.write_text(
                "taxon_id\tgroup\tfamily\thost\taccession\tsource\n"
                "cand|SRR6765940|NODE_1_length_3100\tcandidate\tAmarillovirales\tHylobius abietis\t\tassembly\n",
                encoding="utf-8",
            )
            out_svg = tdir / "fig6.svg"
            rc = plot_phylo_manuscript.main([
                "--mode", "fig6",
                "--panel", f"Amarillovirales|{t1}|{m1}",
                "--panel", f"Phenuiviridae|{t2}|{m1}",
                "--out-svg", str(out_svg),
            ])
            self.assertEqual(rc, 0)
            text = out_svg.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("<svg "))
            self.assertIn("Amarillovirales", text)
            self.assertIn("Phenuiviridae", text)


if __name__ == "__main__":
    unittest.main()
