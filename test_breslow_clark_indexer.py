"""
Unit Test Suite for Breslow Depth & Clark Level Melanoma Indexer
================================================================
Comprehensive verification of AJCC 8th Edition pT criteria, boundary conditions,
Clark invasion microstaging, SLNB recommendations, surgical margins, and CLI.
"""

import csv
import json
import os
import tempfile
import unittest

from breslow_clark_indexer import (
    ClarkLevel,
    AnatomicSite,
    TilCategory,
    SlnbRecommendation,
    MelanomaSpecimenInput,
    TStageAssignment,
    SlnbRiskPrediction,
    SurgicalMarginGuideline,
    MelanomaStagingReport,
    BreslowClarkMelanomaIndexer,
    format_melanoma_report,
)
import cli


class TestAJCC8TStaging(unittest.TestCase):
    """Test Pathologic T category assignments per AJCC 8th Edition rules."""

    def test_tis_in_situ(self):
        stage = BreslowClarkMelanomaIndexer.assign_t_stage(0.0, ulcerated=False, is_in_situ=True)
        self.assertEqual(stage.category, "Tis")
        self.assertTrue(stage.is_in_situ)

    def test_pt1a_boundary_thin_non_ulcerated(self):
        # 0.79 mm, non-ulcerated -> pT1a
        stage = BreslowClarkMelanomaIndexer.assign_t_stage(0.79, ulcerated=False)
        self.assertEqual(stage.category, "pT1a")
        self.assertEqual(stage.base_category, "pT1")
        self.assertEqual(stage.substage_modifier, "a")

    def test_pt1b_under_08_with_ulceration(self):
        # 0.50 mm, ulcerated -> pT1b
        stage = BreslowClarkMelanomaIndexer.assign_t_stage(0.50, ulcerated=True)
        self.assertEqual(stage.category, "pT1b")
        self.assertEqual(stage.substage_modifier, "b")

    def test_pt1b_exact_08_without_ulceration(self):
        # 0.80 mm, non-ulcerated -> pT1b in AJCC 8
        stage = BreslowClarkMelanomaIndexer.assign_t_stage(0.80, ulcerated=False)
        self.assertEqual(stage.category, "pT1b")

    def test_pt1b_exact_10_without_ulceration(self):
        # 1.00 mm, non-ulcerated -> pT1b
        stage = BreslowClarkMelanomaIndexer.assign_t_stage(1.00, ulcerated=False)
        self.assertEqual(stage.category, "pT1b")

    def test_pt2a_and_pt2b(self):
        # 1.01 mm non-ulcerated -> pT2a
        s2a = BreslowClarkMelanomaIndexer.assign_t_stage(1.01, ulcerated=False)
        self.assertEqual(s2a.category, "pT2a")

        # 2.00 mm ulcerated -> pT2b
        s2b = BreslowClarkMelanomaIndexer.assign_t_stage(2.00, ulcerated=True)
        self.assertEqual(s2b.category, "pT2b")

    def test_pt3a_and_pt3b(self):
        # 2.50 mm non-ulcerated -> pT3a
        s3a = BreslowClarkMelanomaIndexer.assign_t_stage(2.50, ulcerated=False)
        self.assertEqual(s3a.category, "pT3a")

        # 4.00 mm ulcerated -> pT3b
        s3b = BreslowClarkMelanomaIndexer.assign_t_stage(4.00, ulcerated=True)
        self.assertEqual(s3b.category, "pT3b")

    def test_pt4a_and_pt4b(self):
        # 4.01 mm non-ulcerated -> pT4a
        s4a = BreslowClarkMelanomaIndexer.assign_t_stage(4.01, ulcerated=False)
        self.assertEqual(s4a.category, "pT4a")

        # 6.50 mm ulcerated -> pT4b
        s4b = BreslowClarkMelanomaIndexer.assign_t_stage(6.50, ulcerated=True)
        self.assertEqual(s4b.category, "pT4b")


class TestSlnbGuidelinesAndNomogram(unittest.TestCase):
    """Test Sentinel Lymph Node Biopsy decision rules and probability calculations."""

    def test_slnb_pt1a_low_risk_not_recommended(self):
        # 0.40 mm, age 60, no ulceration, no LVI, mitoses 0 -> Not recommended
        slnb = BreslowClarkMelanomaIndexer.calculate_slnb_risk(
            breslow_mm=0.40,
            patient_age=60,
            ulcerated=False,
            lvi=False,
            site=AnatomicSite.EXTREMITY,
            mitotic_rate=0.0
        )
        self.assertEqual(slnb.recommendation, SlnbRecommendation.NOT_RECOMMENDED)
        self.assertLess(slnb.probability_pct, 10.0)

    def test_slnb_pt1a_with_high_mitotic_rate_discuss(self):
        # 0.50 mm, mitoses 3.0/mm2 -> Discuss/consider
        slnb = BreslowClarkMelanomaIndexer.calculate_slnb_risk(
            breslow_mm=0.50,
            patient_age=55,
            ulcerated=False,
            lvi=False,
            site=AnatomicSite.TRUNK,
            mitotic_rate=3.0
        )
        self.assertEqual(slnb.recommendation, SlnbRecommendation.DISCUSS_CONSIDER)

    def test_slnb_pt1a_with_lvi_discuss(self):
        slnb = BreslowClarkMelanomaIndexer.calculate_slnb_risk(
            breslow_mm=0.45,
            patient_age=50,
            ulcerated=False,
            lvi=True,
            site=AnatomicSite.TRUNK
        )
        self.assertEqual(slnb.recommendation, SlnbRecommendation.DISCUSS_CONSIDER)

    def test_slnb_pt1b_discuss_and_consider(self):
        # 0.90 mm -> pT1b -> Discuss/consider
        slnb = BreslowClarkMelanomaIndexer.calculate_slnb_risk(
            breslow_mm=0.90,
            patient_age=50,
            ulcerated=False,
            lvi=False,
            site=AnatomicSite.TRUNK
        )
        self.assertEqual(slnb.recommendation, SlnbRecommendation.DISCUSS_CONSIDER)

    def test_slnb_pt2_plus_recommended(self):
        # 1.50 mm -> Recommended
        slnb = BreslowClarkMelanomaIndexer.calculate_slnb_risk(
            breslow_mm=1.50,
            patient_age=50,
            ulcerated=False,
            lvi=False,
            site=AnatomicSite.TRUNK
        )
        self.assertEqual(slnb.recommendation, SlnbRecommendation.RECOMMENDED)

    def test_slnb_microsatellites_forces_recommended(self):
        slnb = BreslowClarkMelanomaIndexer.calculate_slnb_risk(
            breslow_mm=0.50,
            patient_age=50,
            ulcerated=False,
            lvi=False,
            site=AnatomicSite.TRUNK,
            microsatellites=True
        )
        self.assertEqual(slnb.recommendation, SlnbRecommendation.RECOMMENDED)

    def test_nomogram_monotonicity(self):
        # Increasing thickness increases predicted probability
        p1 = BreslowClarkMelanomaIndexer.calculate_slnb_risk(1.0, 50, False, False, AnatomicSite.TRUNK).probability_pct
        p2 = BreslowClarkMelanomaIndexer.calculate_slnb_risk(3.0, 50, False, False, AnatomicSite.TRUNK).probability_pct
        p3 = BreslowClarkMelanomaIndexer.calculate_slnb_risk(3.0, 50, True, True, AnatomicSite.TRUNK).probability_pct
        self.assertLess(p1, p2)
        self.assertLess(p2, p3)


class TestSurgicalMargins(unittest.TestCase):
    """Test Wide Local Excision margin guidance."""

    def test_margin_in_situ(self):
        m = BreslowClarkMelanomaIndexer.get_surgical_margins(0.0, is_in_situ=True)
        self.assertIn("0.5 cm to 1.0 cm", m.recommended_clinical_margin_cm)

    def test_margin_thin_melanoma(self):
        m = BreslowClarkMelanomaIndexer.get_surgical_margins(0.75)
        self.assertEqual(m.recommended_clinical_margin_cm, "1.0 cm")

    def test_margin_intermediate_melanoma(self):
        m = BreslowClarkMelanomaIndexer.get_surgical_margins(1.50)
        self.assertEqual(m.recommended_clinical_margin_cm, "1.0 cm to 2.0 cm")

    def test_margin_thick_melanoma(self):
        m = BreslowClarkMelanomaIndexer.get_surgical_margins(3.20)
        self.assertEqual(m.recommended_clinical_margin_cm, "2.0 cm")


class TestValidationAndFullPipeline(unittest.TestCase):
    """Test specimen input validation, full pipeline execution, and report structure."""

    def test_negative_depth_raises_value_error(self):
        specimen = MelanomaSpecimenInput(specimen_id="ERR-1", breslow_depth_mm=-1.5)
        with self.assertRaises(ValueError):
            BreslowClarkMelanomaIndexer.stage_melanoma(specimen)

    def test_invalid_age_raises_value_error(self):
        specimen = MelanomaSpecimenInput(specimen_id="ERR-2", breslow_depth_mm=1.0, patient_age=150)
        with self.assertRaises(ValueError):
            BreslowClarkMelanomaIndexer.stage_melanoma(specimen)

    def test_full_pipeline_thick_ulcerated(self):
        specimen = MelanomaSpecimenInput(
            specimen_id="SPEC-FULL-01",
            patient_age=65,
            breslow_depth_mm=4.5,
            ulcerated=True,
            clark_level=ClarkLevel.LEVEL_IV,
            lymphovascular_invasion=True,
            perineural_invasion=True,
            microsatellitosis=True,
            anatomic_site=AnatomicSite.HEAD_NECK
        )
        rep = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
        self.assertEqual(rep.t_stage.category, "pT4b")
        self.assertEqual(rep.slnb_evaluation.recommendation, SlnbRecommendation.RECOMMENDED)
        self.assertGreaterEqual(rep.adverse_features_count, 3)
        self.assertTrue(any("Microsatellitosis" in a for a in rep.critical_alerts))
        self.assertIn("pT4b", rep.clinical_summary)

    def test_clark_level_string_mapping(self):
        specimen = MelanomaSpecimenInput(
            specimen_id="SPEC-CLARK",
            breslow_depth_mm=1.2,
            clark_level=ClarkLevel.LEVEL_III
        )
        rep = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
        self.assertIn("Clark Level III", rep.clark_level)


class TestFormattingAndCLI(unittest.TestCase):
    """Test JSON/text formatting and CLI integration."""

    def test_json_and_dict_serialization(self):
        specimen = MelanomaSpecimenInput(
            specimen_id="SPEC-JSON",
            breslow_depth_mm=0.85,
            ulcerated=False
        )
        report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
        d = report.to_dict()
        self.assertEqual(d["specimen_id"], "SPEC-JSON")
        self.assertEqual(d["t_stage"]["category"], "pT1b")

        js = report.to_json()
        parsed = json.loads(js)
        self.assertEqual(parsed["t_stage"]["category"], "pT1b")

    def test_text_report_rendering(self):
        specimen = MelanomaSpecimenInput(
            specimen_id="SPEC-RENDER",
            breslow_depth_mm=1.5,
            ulcerated=True
        )
        report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
        txt = format_melanoma_report(report)
        self.assertIn("CUTANEOUS MELANOMA HISTOPATHOLOGIC STAGING DOSSIER", txt)
        self.assertIn("SPEC-RENDER", txt)
        self.assertIn("pT2b", txt)

    def test_cli_demo_modes(self):
        self.assertEqual(cli.main(["--demo", "in_situ"]), 0)
        self.assertEqual(cli.main(["--demo", "t1a"]), 0)
        self.assertEqual(cli.main(["--demo", "t1b"]), 0)
        self.assertEqual(cli.main(["--demo", "t2b"]), 0)
        self.assertEqual(cli.main(["--demo", "t4b"]), 0)

    def test_cli_batch_csv_processing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_in = os.path.join(tmpdir, "input.csv")
            csv_out = os.path.join(tmpdir, "output.csv")
            with open(csv_in, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["specimen_id", "patient_age", "breslow_depth_mm", "ulcerated"])
                writer.writeheader()
                writer.writerow({"specimen_id": "MEL-01", "patient_age": "50", "breslow_depth_mm": "0.65", "ulcerated": "false"})
                writer.writerow({"specimen_id": "MEL-02", "patient_age": "60", "breslow_depth_mm": "2.20", "ulcerated": "true"})

            ret = cli.main(["--batch-csv", csv_in, "--output", csv_out])
            self.assertEqual(ret, 0)
            self.assertTrue(os.path.exists(csv_out))


if __name__ == "__main__":
    unittest.main()
