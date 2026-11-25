"""Unit tests for Stage 15 ESMFold + Foldseek RdRp structural validation.

Run with:
    python3 -m unittest discover -s analysis_novel/tests -p "test_unit_*.py" -v
"""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import run_esmfold_rdrp  # noqa: E402
import select_divergent_rdrp  # noqa: E402
import summarize_esm_rdrp  # noqa: E402


class TestSelectDivergentRdRp(unittest.TestCase):
    def test_gating_confidence_identity_unresolved_and_dedup(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            conf_tsv = tdir / "contig_confidence.tsv"
            conf_tsv.write_text(
                "contig\trun_id\tconfidence\tkeep\thas_rdrp\tconsensus_family\n"
                "SRR1|NODE_1_length_5000_cov_50\tSRR1\thigh\t1\t1\tPolycipiviridae\n"
                "SRR1|NODE_2_length_5000_cov_50\tSRR1\thigh\t1\t1\tIflaviridae\n"
                "SRR1|NODE_3_length_2500_cov_20\tSRR1\tmedium\t1\t1\t\n"
                "SRR1|NODE_4_length_1200_cov_10\tSRR1\tlow\t0\t1\t\n"
                "SRR1|NODE_5_length_900_cov_10\tSRR1\thigh\t1\t1\tUnclassified\n",
                encoding="utf-8",
            )

            demarc_tsv = tdir / "all_families_demarcation_summary.tsv"
            demarc_tsv.write_text(
                "family\tcandidate\tnearest_ictv_ref\tnearest_identity_pct\n"
                "Polycipiviridae\tSRR1_n1\tICTV_REF1\t24.9\n"
                "Iflaviridae\tSRR1_n2\tICTV_REF2\t52.0\n",
                encoding="utf-8",
            )

            faa = tdir / "SRR1_rdrp_domains.faa"
            seq250 = "ACDEFGHIKLMNPQRSTVWY" * 13  # 260 aa
            seq150 = "ACDEFGHIKLMNPQRSTVWY" * 7   # 140 aa (<200 aa)
            seq900 = "A" * 900                    # >800 aa -> should center-trim
            faa.write_text(
                f">SRR1|NODE_1_length_5000_cov_50|frame=1|aa=10-270|profile=P1|db=RVMT|evalue=1e-20|norm_bitscore_profile=0.8\n{seq250}\n"
                f">SRR1|NODE_1_length_5000_cov_50|frame=1|aa=15-260|profile=P2|db=RVMT|evalue=1e-10|norm_bitscore_profile=0.4\n{seq250}\n"
                f">SRR1|NODE_2_length_5000_cov_50|frame=1|aa=10-270|profile=P1|db=RVMT|evalue=1e-30|norm_bitscore_profile=0.9\n{seq250}\n"
                f">SRR1|NODE_3_length_2500_cov_20|frame=2|aa=1-900|profile=P3|db=NeoRdRp|evalue=1e-12|norm_bitscore_profile=0.6\n{seq900}\n"
                f">SRR1|NODE_4_length_1200_cov_10|frame=1|aa=1-260|profile=P3|db=NeoRdRp|evalue=1e-12|norm_bitscore_profile=0.6\n{seq250}\n"
                f">SRR1|NODE_5_length_900_cov_10|frame=1|aa=1-140|profile=P3|db=NeoRdRp|evalue=1e-12|norm_bitscore_profile=0.6\n{seq150}\n",
                encoding="utf-8",
            )

            rows, recs = select_divergent_rdrp.select_candidates(
                confidence_tsv=conf_tsv,
                rdrp_faa_patterns=[str(faa)],
                demarcation_patterns=[str(demarc_tsv)],
                keep_confidence={"high", "medium"},
                min_aa_len=200,
                max_aa_len=800,
                max_identity_pct=40.0,
                include_unresolved=True,
                dedup_per_contig=True,
            )

            # Expect NODE_1 (high, 24.9% < 40%, deduped to profile P1)
            # and NODE_3 (medium, unresolved, 900aa trimmed to 800aa).
            # Excluded: NODE_2 (52.0% >= 40%), NODE_4 (low), NODE_5 (140aa < 200aa).
            self.assertEqual(len(rows), 2)
            self.assertEqual(len(recs), 2)

            by_contig = {r["contig"]: r for r in rows}
            self.assertIn("SRR1|NODE_1_length_5000_cov_50", by_contig)
            self.assertEqual(by_contig["SRR1|NODE_1_length_5000_cov_50"]["rdrp_profile"], "P1")
            self.assertIn("divergent_aa_id<40pct", by_contig["SRR1|NODE_1_length_5000_cov_50"]["selection_reason"])

            self.assertIn("SRR1|NODE_3_length_2500_cov_20", by_contig)
            self.assertEqual(by_contig["SRR1|NODE_3_length_2500_cov_20"]["folded_aa_len"], "800")
            self.assertEqual(by_contig["SRR1|NODE_3_length_2500_cov_20"]["was_trimmed"], "1")
            self.assertIn("family_unresolved_rdrp", by_contig["SRR1|NODE_3_length_2500_cov_20"]["selection_reason"])


class TestEsmFoldHelpers(unittest.TestCase):
    def test_palm_motif_scan(self) -> None:
        seq = "AAAAADAAADAAAAAASGAATAAAANAAAAAGDDAAAAA"
        res = run_esmfold_rdrp.scan_palm_motifs(seq)
        self.assertEqual(res["motif_a"], "DAAAD")
        self.assertEqual(res["motif_c"], "GDD")
        self.assertEqual(res["palm_motifs_ordered"], "1")

    def test_pdb_plddt_scales(self) -> None:
        # 0-1 scale B-factors (older transformers)
        pdb_unit = (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.80           C\n"
            "ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00  0.60           C\n"
        )
        stats_unit = run_esmfold_rdrp.parse_pdb_plddt(pdb_unit)
        self.assertAlmostEqual(stats_unit["mean_plddt"], 70.0, places=2)
        self.assertAlmostEqual(stats_unit["frac_plddt_ge_70"], 0.5, places=3)

        # 0-100 scale B-factors
        pdb_pct = (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 85.00           C\n"
            "ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 75.00           C\n"
        )
        stats_pct = run_esmfold_rdrp.parse_pdb_plddt(pdb_pct)
        self.assertAlmostEqual(stats_pct["mean_plddt"], 80.0, places=2)
        self.assertAlmostEqual(stats_pct["frac_plddt_ge_70"], 1.0, places=3)


class TestSummarizeEsmRdRp(unittest.TestCase):
    def test_verdict_classification_and_rt_negative_control(self) -> None:
        # High-confidence RdRp validation
        v_hi, _ = summarize_esm_rdrp.assign_verdict(
            mean_plddt=82.5,
            fold_status="ok",
            rdrp_hit={"alntmscore": "0.68", "prob": "0.99", "evalue": "1e-14", "target": "7BV2_A"},
            rt_hit={"alntmscore": "0.31", "evalue": "1e-2", "target": "1RTD_A"},
        )
        self.assertEqual(v_hi, "validated_rdrp_high")

        # RT negative-control flag beats RdRp
        v_rt, _ = summarize_esm_rdrp.assign_verdict(
            mean_plddt=78.0,
            fold_status="ok",
            rdrp_hit={"alntmscore": "0.32", "prob": "0.40", "evalue": "1e-2", "target": "7BV2_A"},
            rt_hit={"alntmscore": "0.64", "evalue": "1e-11", "target": "1RTD_A"},
        )
        self.assertEqual(v_rt, "rt_or_non_rdrp_flag")

        # Low-confidence fold
        v_lo, _ = summarize_esm_rdrp.assign_verdict(
            mean_plddt=41.2,
            fold_status="ok",
            rdrp_hit={"alntmscore": "0.55", "prob": "0.80", "evalue": "1e-6", "target": "7BV2_A"},
            rt_hit=None,
        )
        self.assertEqual(v_lo, "low_confidence_fold")

    def test_end_to_end_summary_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            sel = tdir / "selected.tsv"
            sel.write_text(
                "candidate_id\tcontig\trun_id\tconfidence\tconsensus_family\tnearest_ictv_ref\tnearest_identity_pct\tselection_reason\tdomain_aa_len\tfolded_aa_len\twas_trimmed\n"
                "cand_0001__SRR1_n1\tSRR1|NODE_1\tSRR1\thigh\tPolycipiviridae\tICTV_1\t24.9\tdivergent_aa_id<40pct\t320\t320\t0\n"
                "cand_0002__SRR1_n3\tSRR1|NODE_3\tSRR1\tmedium\tUnresolved\t\t\tfamily_unresolved_rdrp\t280\t280\t0\n",
                encoding="utf-8",
            )
            met = tdir / "metrics.tsv"
            met.write_text(
                "candidate_id\taa_len\tstatus\tpdb_path\tmean_plddt\tmedian_plddt\tfrac_plddt_ge_70\tfrac_plddt_ge_50\tptm\truntime_sec\tmotif_a\tmotif_b\tmotif_c\tpalm_motifs_ordered\n"
                "cand_0001__SRR1_n1\t320\tok\tp1.pdb\t81.2\t83.0\t0.85\t0.95\t0.77\t12.1\tDAAAD\t\tGDD\t1\n"
                "cand_0002__SRR1_n3\t280\tok\tp2.pdb\t74.5\t76.0\t0.72\t0.91\t0.65\t10.4\t\t\t\t0\n",
                encoding="utf-8",
            )
            fs = tdir / "foldseek_rdrp.tsv"
            fs.write_text(
                "cand_0001__SRR1_n1.pdb\t7BV2_A\t0.21\t290\t1e-15\t420\t0.67\t0.65\t0.66\t0.71\t0.99\tRNA-dependent RNA polymerase\n"
                "cand_0002__SRR1_n3.pdb\t1RTD_A\t0.24\t250\t1e-12\t380\t0.62\t0.60\t0.61\t0.68\t0.97\tHIV-1 Reverse Transcriptase\n",
                encoding="utf-8",
            )
            out_val = tdir / "val.tsv"
            out_sum = tdir / "sum.tsv"

            rc = summarize_esm_rdrp.main([
                "--selected-tsv", str(sel),
                "--metrics-tsv", str(met),
                "--foldseek-rdrp", str(fs),
                "--out-validation", str(out_val),
                "--out-summary", str(out_sum),
            ])
            self.assertEqual(rc, 0)

            with out_val.open("r", encoding="utf-8") as fh:
                rows = {r["candidate_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(rows["cand_0001__SRR1_n1"]["structural_verdict"], "validated_rdrp_high")
            self.assertEqual(rows["cand_0002__SRR1_n3"]["structural_verdict"], "rt_or_non_rdrp_flag")


if __name__ == "__main__":
    unittest.main()
