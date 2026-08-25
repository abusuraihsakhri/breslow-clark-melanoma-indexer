#!/usr/bin/env python3
"""
AJCC 8th edition melanoma T-category assignment from microstaging.

T category is defined by Breslow thickness (mm) with ulceration status as
the a/b modifier (mitotic rate is REPORTED but no longer T-defining since
AJCC 8; Clark level is used for historical comparison only):

    <= 1.0 mm   -> T1  (T1a: <0.8 mm non-ulcerated)
                        (T1b: <0.8 mm ulcerated, OR 0.8-1.0 mm any ulceration)
    >1.0-2.0    -> T2a/T2b
    >2.0-4.0    -> T3a/T3b
    >4.0        -> T4a/T4b

Also flags sentinel lymph node biopsy discussion per NCCN guidance
(consider SLNB when >= T1b, i.e. thickness >= 0.8 mm or ulcerated).
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional


@dataclass
class MelanomaSpecimen:
    breslow_mm: float
    ulcerated: bool = False
    mitotic_rate_per_mm2: Optional[float] = None
    clark_level: Optional[int] = None          # historical, I-V
    lymphovascular_invasion: bool = False
    microsatellites: bool = False              # upgrades to N1c regardless of T


def stage_t(specimen: MelanomaSpecimen) -> Dict[str, Any]:
    t = specimen.breslow_mm
    if t <= 0:
        raise ValueError("Breslow thickness must be positive")
    suffix = "b" if specimen.ulcerated else "a"

    if t <= 1.0:
        if t < 0.8:
            category = f"T1{'b' if specimen.ulcerated else 'a'}"
            basis = "<0.8 mm"
        else:
            category = "T1b"  # 0.8-1.0 mm is T1b regardless of ulceration
            basis = "0.8-1.0 mm qualifies as T1b with or without ulceration"
    elif t <= 2.0:
        category, basis = f"T2{suffix}", ">1.0-2.0 mm"
    elif t <= 4.0:
        category, basis = f"T3{suffix}", ">2.0-4.0 mm"
    else:
        category, basis = f"T4{suffix}", ">4.0 mm"

    slnb_discussion = specimen.microsatellites or not category.endswith("a") \
        or t >= 0.8

    notes = []
    if specimen.mitotic_rate_per_mm2 is not None:
        notes.append(f"Mitotic rate {specimen.mitotic_rate_per_mm2}/mm2 recorded "
                     "(report only in AJCC 8; not T-category defining)")
    if specimen.clark_level is not None:
        notes.append(f"Clark level {['I','II','III','IV','V'][specimen.clark_level-1]} "
                     "(historical descriptor; dropped from AJCC 8 staging)")
    if specimen.microsatellites:
        notes.append("Microsatellitess present -> at least N1c; T category still reported")
    if specimen.lymphovascular_invasion:
        notes.append("Lymphovascular invasion documented")

    return {
        "breslow_mm": round(t, 3),
        "ulcerated": specimen.ulcerated,
        "t_category": category,
        "category_basis": basis,
        "slnb_discussion_indicated": bool(slnb_discussion),
        "ajcc_edition": 8,
        "notes": notes,
    }


def stage_group_hint(t_category: str, node_status: str = "cN0") -> str:
    """Coarse clinical stage-group hint from T + nodal status."""
    if "N+" in node_status or node_status.startswith("pN"):
        return "III (regional nodes involved)"
    thin = {"T1a": "IA", "T1b": "IB", "T2a": "IB"}
    if t_category in thin and node_status == "cN0":
        return f"Stage {thin[t_category]}"
    if node_status == "cN0":
        return "Stage IIA-IIC (T2b/T3/T4, node negative)"
    return "Undetermined"


if __name__ == "__main__":
    cases = [
        ("thin non-ulcerated", MelanomaSpecimen(0.5)),
        ("0.9 mm borderline", MelanomaSpecimen(0.9, mitotic_rate_per_mm2=4)),
        ("ulcerated thin", MelanomaSpecimen(0.6, True)),
        ("intermediate", MelanomaSpecimen(1.8, True, 7, clark_level=4)),
        ("thick", MelanomaSpecimen(5.5, True, 12, lymphovascular_invasion=True)),
    ]
    for name, spec in cases:
        r = stage_t(spec)
        print(f"{name:22s} {r['breslow_mm']:>5}mm ulcer={r['ulcerated']!s:5s} "
              f"-> {r['t_category']:4s} | SLNB discuss={r['slnb_discussion_indicated']}")
    print("\nStage hints:", stage_group_hint("T1a"), "|",
          stage_group_hint("T4b"), "|", stage_group_hint("T2a", "pN1"))
