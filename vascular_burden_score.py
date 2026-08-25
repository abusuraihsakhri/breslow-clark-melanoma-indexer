#!/usr/bin/env python3
"""
Vascular Burden Score Calculator for Breslow-Clark Melanoma Indexer.
Calculates VBS incorporating tumor vascularity, lymphovascular invasion,
and angiogenesis markers for enhanced melanoma prognostication.
"""

from typing import Dict, Any, Optional
from dataclasses import dataclass


@dataclass
class VascularFeatures:
    """Vascular features from pathology."""
    microvessel_density: float = 0.0
    lymphovascular_invasion: bool = False
    venous_invasion: bool = False
    tumor_infiltrating_lymphocytes_vascular: bool = False
    necrosis_vascular: bool = False
    vessel_cuffing: bool = False


def calculate_vbs(features: VascularFeatures, breslow_depth_mm: float,
                  clark_level: int) -> Dict[str, Any]:
    """Calculate Vascular Burden Score."""
    score = 0.0
    components = []

    if features.microvessel_density > 20:
        score += 1.0
        components.append(f"MVD >20/mm2: +1.0")
    elif features.microvessel_density > 10:
        score += 0.5
        components.append(f"MVD 10-20/mm2: +0.5")

    if features.lymphovascular_invasion:
        score += 2.0
        components.append("LVI present: +2.0")

    if features.venous_invasion:
        score += 1.5
        components.append("Venous invasion: +1.5")

    if features.necrosis_vascular:
        score += 1.0
        components.append("Necrosis near vessels: +1.0")

    if features.vessel_cuffing:
        score += 0.5
        components.append("Vessel cuffing: +0.5")

    depth_bonus = 0.0
    if breslow_depth_mm > 4.0:
        depth_bonus = 1.0
        components.append(f"Depth >4mm bonus: +1.0")
    elif breslow_depth_mm > 2.0:
        depth_bonus = 0.5
        components.append(f"Depth 2-4mm bonus: +0.5")

    if clark_level >= 4:
        depth_bonus += 0.5
        components.append(f"Clark >= IV bonus: +0.5")

    total_score = score + depth_bonus

    if total_score >= 4.0:
        risk_category = "HIGH"
        prognosis = "Poor prognosis. Consider adjuvant therapy and sentinel lymph node biopsy."
    elif total_score >= 2.0:
        risk_category = "MODERATE"
        prognosis = "Intermediate risk. Sentinel lymph node biopsy recommended."
    elif total_score >= 1.0:
        risk_category = "LOW"
        prognosis = "Low vascular burden. Standard management."
    else:
        risk_category = "MINIMAL"
        prognosis = "Minimal vascular burden."

    return {
        "vbs_score": round(total_score, 2),
        "risk_category": risk_category,
        "prognosis": prognosis,
        "components": components,
        "base_vascular_score": round(score, 2),
        "depth_bonus": round(depth_bonus, 2),
        "breslow_depth_mm": breslow_depth_mm,
        "clark_level": clark_level,
    }


class VascularBurdenAgent:
    """Sub-agent for VBS calculation."""

    def __init__(self):
        self.agent_name = "VascularBurdenAgent"

    def evaluate(self, features: VascularFeatures, breslow_depth_mm: float,
                 clark_level: int) -> Dict[str, Any]:
        """Evaluate vascular burden."""
        result = calculate_vbs(features, breslow_depth_mm, clark_level)
        alerts = []

        if result["risk_category"] == "HIGH":
            alerts.append({
                "type": "HIGH_VASCULAR_BURDEN",
                "severity": "WARNING",
                "message": f"High VBS ({result['vbs_score']:.2f}). "
                           f"LVI: {features.lymphovascular_invasion}.",
                "recommendation": "Consider adjuvant therapy. Strongly consider SLNB and imaging."
            })

        if features.lymphovascular_invasion:
            alerts.append({
                "type": "LYMPHOVASCULAR_INVASION",
                "severity": "WARNING",
                "message": "Lymphovascular invasion detected. Elevated metastatic risk.",
                "recommendation": "Sentinel lymph node biopsy mandatory. Consider systemic staging."
            })

        return {"vbs_result": result, "alerts": alerts}
