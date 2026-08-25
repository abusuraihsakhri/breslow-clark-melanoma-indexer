#!/usr/bin/env python3
"""
Sentinel Lymph Node Prediction for Breslow-Clark Melanoma Indexer.
Predicts sentinel lymph node positivity risk using Breslow depth,
Clark level, ulceration, and vascular burden score.
"""

from typing import Dict, Any, Optional
from dataclasses import dataclass


def predict_slnd_risk(breslow_depth_mm: float, clark_level: int, ulceration: bool,
                      vbs_score: float = 0.0, age_years: float = 50.0,
                      mitotic_rate: float = 0.0) -> Dict[str, Any]:
    """Predict SLN positivity risk using clinical-pathological features."""
    risk_score = 0.0
    risk_factors = []

    if breslow_depth_mm > 4.0:
        risk_score += 3.0
        risk_factors.append(f"Breslow >4mm ({breslow_depth_mm:.1f}mm): +3.0")
    elif breslow_depth_mm > 2.0:
        risk_score += 2.0
        risk_factors.append(f"Breslow 2-4mm ({breslow_depth_mm:.1f}mm): +2.0")
    elif breslow_depth_mm > 1.0:
        risk_score += 1.0
        risk_factors.append(f"Breslow 1-2mm ({breslow_depth_mm:.1f}mm): +1.0")

    if clark_level >= 4:
        risk_score += 1.5
        risk_factors.append(f"Clark >= IV ({clark_level}): +1.5")
    elif clark_level >= 3:
        risk_score += 0.5
        risk_factors.append(f"Clark III: +0.5")

    if ulceration:
        risk_score += 1.5
        risk_factors.append("Ulceration present: +1.5")

    if vbs_score >= 4.0:
        risk_score += 2.0
        risk_factors.append(f"High VBS ({vbs_score:.2f}): +2.0")
    elif vbs_score >= 2.0:
        risk_score += 1.0
        risk_factors.append(f"Moderate VBS ({vbs_score:.2f}): +1.0")

    if mitotic_rate > 5:
        risk_score += 1.0
        risk_factors.append(f"High mitotic rate ({mitotic_rate}/mm2): +1.0")
    elif mitotic_rate > 0:
        risk_score += 0.5
        risk_factors.append(f"Mitotic rate {mitotic_rate}/mm2: +0.5")

    if age_years > 60:
        risk_score += 0.5
        risk_factors.append(f"Age >60 ({age_years:.0f}): +0.5")

    sln_positive_pct = min(95.0, 5.0 + (risk_score * 8.0))

    if risk_score >= 6.0:
        recommendation = "Strongly recommended. High SLN positivity risk."
        urgency = "URGENT"
    elif risk_score >= 3.0:
        recommendation = "Recommended. Moderate SLN positivity risk."
        urgency = "ROUTINE"
    elif risk_score >= 1.5:
        recommendation = "Consider. Low-moderate risk. Discuss with patient."
        urgency = "ELECTIVE"
    else:
        recommendation = "May omit. Low risk. Individualize decision."
        urgency = "OPTIONAL"

    return {
        "predicted_sln_positive_pct": round(sln_positive_pct, 1),
        "risk_score": round(risk_score, 2),
        "recommendation": recommendation,
        "urgency": urgency,
        "risk_factors": risk_factors,
        "breslow_depth_mm": breslow_depth_mm,
        "clark_level": clark_level,
        "ulceration": ulceration,
    }


class SLNPredictionAgent:
    """Sub-agent for SLN prediction."""

    def __init__(self):
        self.agent_name = "SLNPredictionAgent"

    def evaluate(self, breslow_depth_mm: float, clark_level: int, ulceration: bool,
                 vbs_score: float = 0.0, age_years: float = 50.0,
                 mitotic_rate: float = 0.0) -> Dict[str, Any]:
        """Evaluate SLN prediction."""
        result = predict_slnd_risk(breslow_depth_mm, clark_level, ulceration,
                                    vbs_score, age_years, mitotic_rate)
        alerts = []

        if result["predicted_sln_positive_pct"] > 20:
            alerts.append({
                "type": "HIGH_SLN_RISK",
                "severity": "WARNING",
                "message": f"Predicted SLN positivity {result['predicted_sln_positive_pct']:.1f}%.",
                "recommendation": result["recommendation"]
            })

        return {"sln_result": result, "alerts": alerts}
