#!/usr/bin/env python3
"""
Command-Line Interface for Breslow Depth & Clark Level Melanoma Indexer
======================================================================
Supports interactive entry, direct argument parsing, batch CSV evaluation,
benchmark clinical scenarios, and JSON output formatting.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import List, Optional

from breslow_clark_indexer import (
    ClarkLevel,
    AnatomicSite,
    TilCategory,
    MelanomaSpecimenInput,
    BreslowClarkMelanomaIndexer,
    format_melanoma_report,
)


def run_demo(scenario: str = "all") -> int:
    """Runs validated clinical benchmark scenarios."""
    scenarios = {
        "in_situ": MelanomaSpecimenInput(
            specimen_id="DEMO-TIS-01",
            patient_age=48,
            breslow_depth_mm=0.0,
            is_in_situ=True,
            clark_level=ClarkLevel.LEVEL_I,
            anatomic_site=AnatomicSite.TRUNK
        ),
        "t1a": MelanomaSpecimenInput(
            specimen_id="DEMO-T1A-01",
            patient_age=52,
            breslow_depth_mm=0.45,
            ulcerated=False,
            clark_level=ClarkLevel.LEVEL_II,
            mitotic_rate_per_mm2=0.0,
            anatomic_site=AnatomicSite.EXTREMITY
        ),
        "t1b": MelanomaSpecimenInput(
            specimen_id="DEMO-T1B-01",
            patient_age=42,
            breslow_depth_mm=0.85,
            ulcerated=False,
            clark_level=ClarkLevel.LEVEL_III,
            mitotic_rate_per_mm2=2.5,
            anatomic_site=AnatomicSite.TRUNK
        ),
        "t2b": MelanomaSpecimenInput(
            specimen_id="DEMO-T2B-01",
            patient_age=63,
            breslow_depth_mm=1.65,
            ulcerated=True,
            clark_level=ClarkLevel.LEVEL_IV,
            lymphovascular_invasion=True,
            anatomic_site=AnatomicSite.HEAD_NECK
        ),
        "t4b": MelanomaSpecimenInput(
            specimen_id="DEMO-T4B-01",
            patient_age=70,
            breslow_depth_mm=5.20,
            ulcerated=True,
            clark_level=ClarkLevel.LEVEL_V,
            microsatellitosis=True,
            perineural_invasion=True,
            anatomic_site=AnatomicSite.ACRAL
        )
    }

    selected = scenarios.items() if scenario == "all" else [(scenario, scenarios[scenario])] if scenario in scenarios else []
    if not selected:
        print(f"Unknown scenario: {scenario}. Choose from: {list(scenarios.keys())} or 'all'")
        return 1

    for name, specimen in selected:
        report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
        print(format_melanoma_report(report))
        print("\n")
    return 0


def interactive_mode() -> int:
    """Guides the user through entering melanoma pathology parameters."""
    print("=" * 60)
    print(" Cutaneous Melanoma Microstaging - Interactive Entry")
    print("=" * 60)
    try:
        specimen_id = input("Enter Specimen ID [MEL-2026-001]: ").strip() or "MEL-2026-001"
        age_str = input("Patient Age in years [55]: ").strip() or "55"
        age = int(age_str)

        in_situ_str = input("Is this Melanoma in situ (Tis)? (y/n) [n]: ").strip().lower()
        is_in_situ = in_situ_str in ("y", "yes", "true", "1")

        if is_in_situ:
            depth = 0.0
        else:
            depth_str = input("Breslow Depth in mm [0.75]: ").strip() or "0.75"
            depth = float(depth_str)

        ulc_str = input("Primary tumor ulceration present? (y/n) [n]: ").strip().lower()
        ulcerated = ulc_str in ("y", "yes", "true", "1")

        clark_str = input("Clark Invasion Level (1-5, or enter to skip) [3]: ").strip() or "3"
        clark = ClarkLevel(int(clark_str)) if clark_str.isdigit() and int(clark_str) in [1, 2, 3, 4, 5] else None

        mit_str = input("Mitotic Rate (mitoses/mm2) [0.0]: ").strip() or "0.0"
        mitoses = float(mit_str)

        lvi_str = input("Lymphovascular invasion (LVI)? (y/n) [n]: ").strip().lower()
        lvi = lvi_str in ("y", "yes", "true", "1")

        site_str = input("Anatomic Site (trunk/extremity/head_neck/acral/mucosal) [trunk]: ").strip().lower() or "trunk"
        site = AnatomicSite(site_str) if site_str in [s.value for s in AnatomicSite] else AnatomicSite.TRUNK

        specimen = MelanomaSpecimenInput(
            specimen_id=specimen_id,
            patient_age=age,
            breslow_depth_mm=depth,
            is_in_situ=is_in_situ,
            ulcerated=ulcerated,
            clark_level=clark,
            mitotic_rate_per_mm2=mitoses,
            lymphovascular_invasion=lvi,
            anatomic_site=site
        )

        report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
        print("\n" + format_melanoma_report(report))
        return 0

    except Exception as e:
        print(f"Error during interactive entry: {e}", file=sys.stderr)
        return 1


def process_batch_csv(input_csv: str, output_csv: Optional[str] = None) -> int:
    """Processes a CSV file containing melanoma pathology records."""
    try:
        with open(input_csv, mode="r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        if not rows:
            print("Error in batch processing: input CSV contains no data rows", file=sys.stderr)
            return 1

        results = []
        for r in rows:
            specimen_id = (
                r.get("specimen_id") or r.get("id") or r.get("patient_id")
                or r.get("Patient_ID") or r.get("specimen") or "SPEC-001"
            )
            age_raw = r.get("patient_age") or r.get("age") or "55"
            try:
                age = int(float(age_raw))
            except (ValueError, TypeError):
                age = 55

            depth_raw = (
                r.get("breslow_depth_mm") or r.get("breslow_thickness_mm")
                or r.get("breslow_mm") or r.get("depth") or r.get("breslow")
                or r.get("v1") or "0.0"
            )
            try:
                depth = float(depth_raw)
            except (ValueError, TypeError):
                depth = 0.0

            in_situ_raw = str(r.get("is_in_situ", "")).strip().lower()
            in_situ = in_situ_raw in ("true", "1", "yes", "in_situ", "in situ")

            ulc_raw = str(
                r.get("ulcerated") or r.get("ulceration") or r.get("ulceration_status") or "false"
            ).strip().lower()
            ulcerated = ulc_raw in ("true", "1", "yes", "present")

            clark_raw = str(r.get("clark_level") or r.get("clark") or "").strip().upper()
            clark_val = None
            clark_map = {
                "1": ClarkLevel.LEVEL_I, "I": ClarkLevel.LEVEL_I, "LEVEL_I": ClarkLevel.LEVEL_I,
                "2": ClarkLevel.LEVEL_II, "II": ClarkLevel.LEVEL_II, "LEVEL_II": ClarkLevel.LEVEL_II,
                "3": ClarkLevel.LEVEL_III, "III": ClarkLevel.LEVEL_III, "LEVEL_III": ClarkLevel.LEVEL_III,
                "4": ClarkLevel.LEVEL_IV, "IV": ClarkLevel.LEVEL_IV, "LEVEL_IV": ClarkLevel.LEVEL_IV,
                "5": ClarkLevel.LEVEL_V, "V": ClarkLevel.LEVEL_V, "LEVEL_V": ClarkLevel.LEVEL_V,
            }
            if clark_raw in clark_map:
                clark_val = clark_map[clark_raw]

            mitotic_raw = (
                r.get("mitotic_rate_per_mm2") or r.get("mitotic_rate")
                or r.get("mitoses") or "0.0"
            )
            try:
                mitoses = float(mitotic_raw)
            except (ValueError, TypeError):
                mitoses = 0.0

            lvi_raw = str(
                r.get("lymphovascular_invasion") or r.get("lvi") or "false"
            ).strip().lower()
            lvi = lvi_raw in ("true", "1", "yes", "present")

            perineural_raw = str(
                r.get("perineural_invasion") or r.get("neurotropism")
                or r.get("perineural") or "false"
            ).strip().lower()
            perineural = perineural_raw in ("true", "1", "yes", "present")

            micro_raw = str(
                r.get("microsatellitosis") or r.get("microscopic_satellitosis")
                or r.get("satellitosis") or "false"
            ).strip().lower()
            microsatellites = micro_raw in ("true", "1", "yes", "present")

            site_val = str(r.get("anatomic_site") or r.get("site") or "trunk").lower()
            site = AnatomicSite(site_val) if site_val in [s.value for s in AnatomicSite] else AnatomicSite.TRUNK

            specimen = MelanomaSpecimenInput(
                specimen_id=specimen_id,
                patient_age=age,
                breslow_depth_mm=depth,
                is_in_situ=in_situ,
                ulcerated=ulcerated,
                clark_level=clark_val,
                mitotic_rate_per_mm2=mitoses,
                lymphovascular_invasion=lvi,
                perineural_invasion=perineural,
                microsatellitosis=microsatellites,
                anatomic_site=site
            )
            rep = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
            row_res = dict(r)
            row_res["t_stage"] = rep.t_stage.category
            row_res["clark_level_staged"] = rep.clark_level
            row_res["slnb_recommendation"] = rep.slnb_evaluation.recommendation.value
            row_res["slnb_probability_pct"] = rep.slnb_evaluation.probability_pct
            row_res["recommended_margins"] = rep.surgical_margins.recommended_clinical_margin_cm
            row_res["margin_rationale"] = rep.surgical_margins.guideline_rationale
            row_res["adverse_features_count"] = rep.adverse_features_count
            row_res["adverse_features"] = "; ".join(rep.adverse_features_list) if rep.adverse_features_list else "None"
            results.append(row_res)

        if output_csv:
            with open(output_csv, mode="w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
                writer.writeheader()
                writer.writerows(results)
            print(f"Successfully processed {len(results)} records -> {output_csv}")
        else:
            print(json.dumps(results, indent=2))
        return 0

    except Exception as e:
        print(f"Error in batch processing: {e}", file=sys.stderr)
        return 1


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    # Check if 'batch' subcommand is used
    if argv and argv[0] == "batch":
        batch_parser = argparse.ArgumentParser(
            prog="cli.py batch",
            description="Batch evaluate melanoma pathology CSV records"
        )
        batch_parser.add_argument("-i", "--input", "--batch-csv", dest="batch_csv", required=True,
                                  help="Input CSV file for batch processing")
        batch_parser.add_argument("-o", "--output", dest="output", default=None,
                                  help="Output CSV or JSON file path")
        args = batch_parser.parse_args(argv[1:])
        return process_batch_csv(args.batch_csv, args.output)

    parser = argparse.ArgumentParser(
        description="Cutaneous Melanoma Histopathologic Staging & SLNB Indexer (AJCC 8th Ed / NCCN)"
    )
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive microstaging mode")
    parser.add_argument("--demo", choices=["in_situ", "t1a", "t1b", "t2b", "t4b", "all"], help="Run benchmark demo scenario")
    parser.add_argument("--specimen-id", default="MEL-001", help="Pathology specimen accession ID")
    parser.add_argument("--age", type=int, default=55, help="Patient age in years")
    parser.add_argument("--depth", type=float, default=0.75, help="Breslow depth in millimeters (mm)")
    parser.add_argument("--in-situ", action="store_true", help="Melanoma in situ (Tis)")
    parser.add_argument("--ulcerated", action="store_true", help="Ulceration present")
    parser.add_argument("--clark-level", type=int, choices=[1, 2, 3, 4, 5], help="Clark Level of invasion (I-V)")
    parser.add_argument("--mitoses", type=float, default=0.0, help="Mitotic rate (mitoses/mm2)")
    parser.add_argument("--lvi", action="store_true", help="Lymphovascular invasion present")
    parser.add_argument("--perineural", action="store_true", help="Perineural invasion present")
    parser.add_argument("--microsatellites", action="store_true", help="Microsatellitosis present")
    parser.add_argument("--tils", choices=["absent", "non-brisk", "brisk"], default="non-brisk", help="Tumor-infiltrating lymphocytes")
    parser.add_argument("--regression", action="store_true", help="Histologic regression present")
    parser.add_argument("--site", choices=["trunk", "extremity", "head_neck", "acral"], default="trunk", help="Cutaneous primary site")

    parser.add_argument("--batch-csv", help="Input CSV file for batch processing")
    parser.add_argument("--output", "-o", help="Output file path (CSV or JSON)")
    parser.add_argument("--file", "-f", help="Load specimen input JSON file")
    parser.add_argument("--json", "-j", action="store_true", help="Output result as JSON")

    args = parser.parse_args(argv)

    if args.interactive:
        return interactive_mode()

    if args.demo:
        return run_demo(args.demo)

    if args.batch_csv:
        return process_batch_csv(args.batch_csv, args.output)

    if args.file:
        with open(args.file, "r") as fp:
            data = json.load(fp)
        clark_val = ClarkLevel(data["clark_level"]) if "clark_level" in data and data["clark_level"] else None
        specimen = MelanomaSpecimenInput(
            specimen_id=data.get("specimen_id", "FILE-SPEC"),
            patient_age=data.get("patient_age", 55),
            breslow_depth_mm=data.get("breslow_depth_mm", 0.75),
            is_in_situ=data.get("is_in_situ", False),
            ulcerated=data.get("ulcerated", False),
            clark_level=clark_val,
            mitotic_rate_per_mm2=data.get("mitotic_rate_per_mm2", 0.0),
            lymphovascular_invasion=data.get("lymphovascular_invasion", False),
            perineural_invasion=data.get("perineural_invasion", False),
            microsatellitosis=data.get("microsatellitosis", False),
            anatomic_site=AnatomicSite(data.get("anatomic_site", "trunk"))
        )
    else:
        clark = ClarkLevel(args.clark_level) if args.clark_level else None
        til_map = {"absent": TilCategory.ABSENT, "non-brisk": TilCategory.NON_BRISK, "brisk": TilCategory.BRISK}
        specimen = MelanomaSpecimenInput(
            specimen_id=args.specimen_id,
            patient_age=args.age,
            breslow_depth_mm=0.0 if args.in_situ else args.depth,
            is_in_situ=args.in_situ,
            ulcerated=args.ulcerated,
            clark_level=clark,
            mitotic_rate_per_mm2=args.mitoses,
            lymphovascular_invasion=args.lvi,
            perineural_invasion=args.perineural,
            microsatellitosis=args.microsatellites,
            tumor_infiltrating_lymphocytes=til_map.get(args.tils, TilCategory.NON_BRISK),
            regression_present=args.regression,
            anatomic_site=AnatomicSite(args.site)
        )

    report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)

    if args.json:
        if args.output:
            with open(args.output, "w") as out_fp:
                out_fp.write(report.to_json())
        else:
            print(report.to_json())
    else:
        formatted = format_melanoma_report(report)
        if args.output:
            with open(args.output, "w") as out_fp:
                out_fp.write(formatted)
        else:
            print(formatted)

    return 0


if __name__ == "__main__":
    sys.exit(main())

