#!/usr/bin/env python3
"""
Sentinel lymph node biopsy positivity risk model for cutaneous melanoma.

Logistic regression calibrated against the ranges published by the MSKCC
melanoma nomogram and AJCC 8 staging cohorts:

    logit(p) = b0 + b1*Breslow(mm) + b2*ulceration + b3*(age - 50)/10
               + b4*LVI + b5*head_neck_site

Coefficients are an approximation of published nomograms (transparent,
auditable weights rather than a black box); probabilities are intended for
SLNB decision support, not prognosis.
"""

import math
from dataclasses import dataclass


INTERCEPT = -2.80
COEF_THICKNESS_PER_MM = 0.55
COEF_ULCERATION = 0.90
COEF_AGE_PER_DECADE_OVER_50 = -0.30
COEF_LYMPHOVASCULAR_INVASION = 0.80
COEF_HEAD_NECK_SITE = 0.30

SITE_ADJUSTMENTS = {
    "head_neck": COEF_HEAD_NECK_SITE,
    "trunk": 0.10,
    "extremity": 0.0,
    "acral": 0.20,
}


@dataclass
class SlnbCandidate:
    breslow_mm: float
    age_years: float = 50.0
    ulcerated: bool = False
    lymphovascular_invasion: bool = False
    site: str = "trunk"


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def predict_slnb_positivity(c: SlnbCandidate) -> dict:
    """Return probability of sentinel-node metastasis with component breakdown."""
    site = c.site.lower()
    if site not in SITE_ADJUSTMENTS:
        raise ValueError(f"site must be one of {sorted(SITE_ADJUSTMENTS)}")
    terms = {
        "thickness": COEF_THICKNESS_PER_MM * c.breslow_mm,
        "ulceration": COEF_ULCERATION * int(c.ulcerated),
        "age_decades_over_50": COEF_AGE_PER_DECADE_OVER_50 * ((c.age_years - 50.0) / 10.0),
        "lymphovascular_invasion": COEF_LYMPHOVASCULAR_INVASION * int(c.lymphovascular_invasion),
        "site_adjustment": SITE_ADJUSTMENTS[site],
    }
    logit = INTERCEPT + sum(terms.values())
    p = _sigmoid(logit)

    if p < 0.05:
        band, action = "very low", "SLNB generally not indicated (<5%)"
    elif p < 0.10:
        band, action = "low", "Discuss SLNB; many guidelines offer observation"
    elif p < 0.25:
        band, action = "moderate", "SLNB recommended"
    else:
        band, action = "high", "SLNB strongly recommended; complete node dissection if positive"

    return {
        "probability_slnb_positive": round(p, 4),
        "risk_band": band,
        "recommendation": action,
        "logit": round(logit, 3),
        "contributions": {k: round(v, 3) for k, v in terms.items()},
        "model": "logistic approximation of MSKCC-style melanoma SLNB nomogram",
    }


if __name__ == "__main__":
    print("SLNB positivity predictions")
    print("-" * 74)
    cases = [
        ("young thin trunk", SlnbCandidate(0.7, 38)),
        ("1 mm extremity", SlnbCandidate(1.0, 62)),
        ("2 mm ulcerated head/neck", SlnbCandidate(2.1, 70, True, site="head_neck")),
        ("4 mm LVI+ acral", SlnbCandidate(4.2, 55, True, True, "acral")),
        ("elderly thin LVI-", SlnbCandidate(0.85, 81)),
    ]
    for name, c in cases:
        r = predict_slnb_positivity(c)
        print(f"{name:26s} p(SLNBN+)={r['probability_slnb_positive']:.3f} "
              f"[{r['risk_band']:9s}] -> {r['recommendation']}")
