#!/usr/bin/env python3
"""
Microscopic Staging Features for Breslow-Clark Melanoma Indexer.
Extracts and quantifies microscopic features beyond Breslow/Clark
for enhanced melanoma staging and prognostication.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass


@dataclass
class MicroscopicFeatures:
    """Microscopic staging features."""
    tumor_infiltrating_lymphocytes: str = "absent"  # absent, brisk, non-brisk
    regression: bool = False
    solar_elastosis: bool = False
    vertical_growth_phase: bool = True
    horizontal_growth_phase: bool = False
    microsatellites: bool = False
    in_transit_metastases: bool = False
    neurotropism: bool = False
    pre_existing_nevus: bool = False
    tumor_boarder: str = "well_circumscribed"  # well_circumscribed, poorly_defined, infiltrative
    maturation_pattern: str = "present"  # present, absent, partial
    mitotic_figure_position: str = "superficial"  # superficial, deep, diffuse


def extract_microscopic_features(features: MicroscopicFeatures) -> Dict[str, Any]:
    """Extract and quantify microscopic staging features."""
    extracted = {
        "til_status": features.tumor_infiltrating_lymphocytes,
        "til_score": 0,
        "regression_present": features.regression,
        "growth_phase": "VGP" if features.vertical_growth_phase else "HGP",
        "microsatellites_present": features.microsatellites,
        "neurotropism_present": features.neurotropism,
        "tumor_border": features.tumor_boarder,
        "maturation_pattern": features.maturation_pattern,
        "mitotic_position": features.mitotic_figure_position,
    }

    til_scores = {"absent": 0, "non-brisk": 1, "brisk": 2}
    extracted["til_score"] = til_scores.get(features.tumor_infiltrating_lymphocytes, 0)

    staging_implications = []
    if features.tumor_infiltrating_lymphocytes == "brisk":
        staging_implications.append("Brisk TILs: favorable prognostic indicator")
    elif features.tumor_infiltrating_lymphocytes == "absent":
        staging_implications.append("Absent TILs: may indicate immune evasion")

    if features.regression:
        staging_implications.append("Regression present: consider Breslow depth uncertainty")

    if features.microsatellites:
        staging_implications.append("Microsatellites: upstage to III if confirmed")

    if features.neurotropism:
        staging_implications.append("Neurotropism: increased local recurrence risk")

    if features.maturation_pattern == "absent":
        staging_implications.append("Absent maturation: higher grade tumor")

    return {
        "extracted_features": extracted,
        "staging_implications": staging_implications,
        "feature_count": len([v for v in extracted.values() if v and v not in (0, "absent", False)]),
    }


class MicroscopicStagingAgent:
    """Sub-agent for microscopic staging features."""

    def __init__(self):
        self.agent_name = "MicroscopicStagingAgent"

    def evaluate(self, features: MicroscopicFeatures) -> Dict[str, Any]:
        """Evaluate microscopic features."""
        result = extract_microscopic_features(features)
        alerts = []

        if features.microsatellites:
            alerts.append({
                "type": "MICROSATELLITES_DETECTED",
                "severity": "WARNING",
                "message": "Microsatellites identified. May upstage disease.",
                "recommendation": "Confirm microsatellites on deeper sections. Stage accordingly."
            })

        if features.neurotropism:
            alerts.append({
                "type": "NEUROTROPISM",
                "severity": "WARNING",
                "message": "Neurotropism detected. Increased local recurrence risk.",
                "recommendation": "Consider wider local excision margins. Close follow-up."
            })

        if features.tumor_infiltrating_lymphocytes == "absent":
            alerts.append({
                "type": "NO_TILS",
                "severity": "ADVISORY",
                "message": "No tumor-infiltrating lymphocytes detected.",
                "recommendation": "May indicate immune evasion. Consider immunotherapy evaluation."
            })

        return {"microscopic_result": result, "alerts": alerts}
