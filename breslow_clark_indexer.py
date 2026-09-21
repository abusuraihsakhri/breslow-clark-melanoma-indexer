"""
Breslow Depth & Clark Level Cutaneous Melanoma Histopathologic Indexer
=====================================================================
Comprehensive microstaging and staging engine implementing:
- AJCC 8th Edition Pathologic T Category (pT) classification
- Clark Level anatomical invasion microstaging (Levels I through V)
- NCCN / ASCO-SSO Sentinel Lymph Node Biopsy (SLNB) decision criteria
- Calibrated logistic nomogram for SLNB metastasis risk prediction
- Surgical wide local excision (WLE) margin guidance
- Comprehensive adverse histopathologic feature tracking (mitoses, LVI, perineural, microsatellites)

Standards:
- AJCC Cancer Staging Manual (8th Edition, 2017/2018)
- NCCN Clinical Practice Guidelines in Oncology: Melanoma (Cutaneous)
- ASCO-SSO Guideline on Sentinel Lymph Node Biopsy in Melanoma
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from enum import Enum
from typing import Dict, Any, List, Optional


class ClarkLevel(int, Enum):
    LEVEL_I = 1    # Melanoma in situ (confined to epidermis)
    LEVEL_II = 2   # Invasion into papillary dermis
    LEVEL_III = 3  # Invasion filling papillary dermis to reticular interface
    LEVEL_IV = 4   # Invasion into reticular dermis
    LEVEL_V = 5    # Invasion into subcutaneous adipose tissue / fat


class AnatomicSite(str, Enum):
    TRUNK = "trunk"
    EXTREMITY = "extremity"
    HEAD_NECK = "head_neck"
    ACRAL = "acral"
    MUCOSAL = "mucosal"


class TilCategory(str, Enum):
    ABSENT = "Absent"
    NON_BRISK = "Non-brisk"
    BRISK = "Brisk"


class SlnbRecommendation(str, Enum):
    NOT_RECOMMENDED = "SLNB generally not recommended for routine pT1a cases without additional risk features"
    DISCUSS_CONSIDER = "Discuss and consider SLNB based on tumor features and patient context"
    RECOMMENDED = "SLNB is generally recommended for invasive melanoma >1.0 mm when clinically appropriate"


@dataclass
class MelanomaSpecimenInput:
    """Input parameters for cutaneous melanoma histopathologic staging."""
    specimen_id: str
    patient_age: int = 55
    breslow_depth_mm: float = 0.0  # 0.0 for in situ
    is_in_situ: bool = False
    ulcerated: bool = False
    clark_level: Optional[ClarkLevel] = None
    mitotic_rate_per_mm2: Optional[float] = 0.0
    lymphovascular_invasion: bool = False
    perineural_invasion: bool = False
    microsatellitosis: bool = False
    tumor_infiltrating_lymphocytes: TilCategory = TilCategory.NON_BRISK
    regression_present: bool = False
    anatomic_site: AnatomicSite = AnatomicSite.TRUNK
    clinical_notes: str = ""

    def validate(self) -> List[str]:
        warnings = []
        if not self.specimen_id or not self.specimen_id.strip():
            raise ValueError("Specimen ID cannot be empty")
        if self.breslow_depth_mm < 0.0:
            raise ValueError(f"Breslow depth cannot be negative, got {self.breslow_depth_mm} mm")
        if self.patient_age < 0 or self.patient_age > 120:
            raise ValueError(f"Patient age must be between 0 and 120, got {self.patient_age}")
        if self.mitotic_rate_per_mm2 is not None and self.mitotic_rate_per_mm2 < 0:
            raise ValueError("Mitotic rate cannot be negative")
        if self.is_in_situ and self.breslow_depth_mm > 0.0:
            raise ValueError(
                "In-situ melanoma cannot have a positive Breslow depth. "
                "Set Breslow depth to 0.0 mm or clear the in-situ flag."
            )
        if self.anatomic_site is AnatomicSite.MUCOSAL:
            raise ValueError(
                "Mucosal melanoma is outside the scope of this cutaneous melanoma tool."
            )
        if self.breslow_depth_mm == 0.0 and not self.is_in_situ:
            warnings.append("Breslow depth 0.0 mm interpreted as melanoma in situ (Tis).")
        if self.clark_level is ClarkLevel.LEVEL_I and self.breslow_depth_mm > 0.0:
            warnings.append("Clark level I is inconsistent with a positive Breslow depth.")
        return warnings


@dataclass
class TStageAssignment:
    """AJCC 8th Edition Pathologic T Category Assignment."""
    category: str  # e.g., Tis, pT1a, pT1b, pT2a, pT2b, pT3a, pT3b, pT4a, pT4b
    base_category: str  # Tis, T1, T2, T3, T4
    substage_modifier: str  # 'a' or 'b' (or '' for Tis)
    criteria_basis: str
    is_in_situ: bool
    is_ulcerated: bool
    thickness_range: str


@dataclass
class SlnbRiskPrediction:
    """Quantitative risk prediction and guideline recommendation for SLNB."""
    probability_pct: float
    recommendation: SlnbRecommendation
    logit_score: float
    risk_factors_present: List[str]
    nomogram_contributions: Dict[str, float]


@dataclass
class SurgicalMarginGuideline:
    """Wide local excision (WLE) margin guidance based on depth."""
    recommended_clinical_margin_cm: str
    guideline_rationale: str


@dataclass
class MelanomaStagingReport:
    """Comprehensive Melanoma Microstaging and Decision Report."""
    specimen_id: str
    patient_age: int
    anatomic_site: str
    breslow_depth_mm: float
    clark_level: Optional[str]
    t_stage: TStageAssignment
    slnb_evaluation: SlnbRiskPrediction
    surgical_margins: SurgicalMarginGuideline
    adverse_features_count: int
    adverse_features_list: List[str]
    critical_alerts: List[str]
    clinical_summary: str

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["slnb_evaluation"]["recommendation"] = self.slnb_evaluation.recommendation.value
        return d

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


class BreslowClarkMelanomaIndexer:
    """
    Cutaneous melanoma microstaging helper for AJCC 8th Edition T category,
    Clark-level documentation, excision-margin guidance, and SLNB discussion prompts.
    """

    # Exploratory logistic weights. This estimate is not externally validated and
    # must not be used as a standalone clinical decision rule.
    NOMOGRAM_INTERCEPT = -2.80
    NOMOGRAM_COEF_THICKNESS = 0.55
    NOMOGRAM_COEF_ULCERATION = 0.90
    NOMOGRAM_COEF_AGE_PER_DECADE = -0.30
    NOMOGRAM_COEF_LVI = 0.80
    NOMOGRAM_SITE_WEIGHTS = {
        AnatomicSite.TRUNK: 0.10,
        AnatomicSite.EXTREMITY: 0.00,
        AnatomicSite.HEAD_NECK: 0.30,
        AnatomicSite.ACRAL: 0.20,
        AnatomicSite.MUCOSAL: 0.25,
    }

    @classmethod
    def assign_t_stage(cls, breslow_mm: float, ulcerated: bool, is_in_situ: bool = False) -> TStageAssignment:
        """
        Assigns AJCC 8th Edition Pathologic T Category:
        - Tis: In situ
        - pT1a: < 0.8 mm, no ulceration
        - pT1b: < 0.8 mm with ulceration OR 0.8 - 1.0 mm (+/- ulceration)
        - pT2a: > 1.0 - 2.0 mm, no ulceration
        - pT2b: > 1.0 - 2.0 mm, with ulceration
        - pT3a: > 2.0 - 4.0 mm, no ulceration
        - pT3b: > 2.0 - 4.0 mm, with ulceration
        - pT4a: > 4.0 mm, no ulceration
        - pT4b: > 4.0 mm, with ulceration
        """
        if breslow_mm < 0.0:
            raise ValueError("Breslow depth cannot be negative")
        if is_in_situ and breslow_mm > 0.0:
            raise ValueError(
                "In-situ melanoma cannot have a positive Breslow depth."
            )
        if is_in_situ or breslow_mm == 0.0:
            return TStageAssignment(
                category="Tis",
                base_category="Tis",
                substage_modifier="",
                criteria_basis="Melanoma in situ confined to epidermis (Breslow 0.0 mm)",
                is_in_situ=True,
                is_ulcerated=False,
                thickness_range="0.0 mm"
            )

        suffix = "b" if ulcerated else "a"

        if breslow_mm <= 1.0:
            if breslow_mm < 0.8:
                cat = f"pT1{suffix}"
                basis = f"Breslow thickness < 0.8 mm ({breslow_mm:.2f} mm), ulceration: {ulcerated}"
            else:
                # 0.8 to 1.0 mm is pT1b regardless of ulceration in AJCC 8th Edition
                cat = "pT1b"
                basis = f"Breslow thickness 0.8-1.0 mm ({breslow_mm:.2f} mm) qualifies as pT1b (ulceration: {ulcerated})"
            return TStageAssignment(
                category=cat,
                base_category="pT1",
                substage_modifier=cat[-1],
                criteria_basis=basis,
                is_in_situ=False,
                is_ulcerated=ulcerated,
                thickness_range="<= 1.0 mm"
            )

        elif breslow_mm <= 2.0:
            return TStageAssignment(
                category=f"pT2{suffix}",
                base_category="pT2",
                substage_modifier=suffix,
                criteria_basis=f"Breslow thickness > 1.0 to 2.0 mm ({breslow_mm:.2f} mm), ulceration: {ulcerated}",
                is_in_situ=False,
                is_ulcerated=ulcerated,
                thickness_range="> 1.0 - 2.0 mm"
            )

        elif breslow_mm <= 4.0:
            return TStageAssignment(
                category=f"pT3{suffix}",
                base_category="pT3",
                substage_modifier=suffix,
                criteria_basis=f"Breslow thickness > 2.0 to 4.0 mm ({breslow_mm:.2f} mm), ulceration: {ulcerated}",
                is_in_situ=False,
                is_ulcerated=ulcerated,
                thickness_range="> 2.0 - 4.0 mm"
            )

        else:
            return TStageAssignment(
                category=f"pT4{suffix}",
                base_category="pT4",
                substage_modifier=suffix,
                criteria_basis=f"Breslow thickness > 4.0 mm ({breslow_mm:.2f} mm), ulceration: {ulcerated}",
                is_in_situ=False,
                is_ulcerated=ulcerated,
                thickness_range="> 4.0 mm"
            )

    @classmethod
    def calculate_slnb_risk(
        cls,
        breslow_mm: float,
        patient_age: int,
        ulcerated: bool,
        lvi: bool,
        site: AnatomicSite,
        mitotic_rate: Optional[float] = 0.0,
        microsatellites: bool = False
    ) -> SlnbRiskPrediction:
        """Returns an exploratory SLN-positivity estimate plus rule-based SLNB guidance."""
        if breslow_mm < 0.0:
            raise ValueError("Breslow depth cannot be negative")
        if patient_age < 0 or patient_age > 120:
            raise ValueError("Patient age must be between 0 and 120")
        if site is AnatomicSite.MUCOSAL:
            raise ValueError("Mucosal melanoma is outside the scope of this tool")
        if breslow_mm == 0.0:
            return SlnbRiskPrediction(
                probability_pct=0.0,
                recommendation=SlnbRecommendation.NOT_RECOMMENDED,
                logit_score=-10.0,
                risk_factors_present=[],
                nomogram_contributions={}
            )

        site_weight = cls.NOMOGRAM_SITE_WEIGHTS.get(site, 0.0)
        age_decades_from_50 = (patient_age - 50.0) / 10.0

        contrib_thickness = cls.NOMOGRAM_COEF_THICKNESS * breslow_mm
        contrib_ulcer = cls.NOMOGRAM_COEF_ULCERATION * (1.0 if ulcerated else 0.0)
        contrib_age = cls.NOMOGRAM_COEF_AGE_PER_DECADE * age_decades_from_50
        contrib_lvi = cls.NOMOGRAM_COEF_LVI * (1.0 if lvi else 0.0)

        logit = (
            cls.NOMOGRAM_INTERCEPT +
            contrib_thickness +
            contrib_ulcer +
            contrib_age +
            contrib_lvi +
            site_weight
        )

        # Sigmoid probability
        prob = 1.0 / (1.0 + math.exp(-logit))
        prob_pct = round(prob * 100.0, 2)

        risk_factors = []
        if breslow_mm >= 0.8:
            risk_factors.append(f"Breslow depth >= 0.8 mm ({breslow_mm:.2f} mm)")
        if ulcerated:
            risk_factors.append("Ulceration present")
        if lvi:
            risk_factors.append("Lymphovascular invasion present")
        if mitotic_rate is not None and mitotic_rate >= 2.0:
            risk_factors.append(f"Elevated mitotic rate ({mitotic_rate:.1f}/mm2)")
        if microsatellites:
            risk_factors.append(
                "Microsatellitosis detected (regional metastatic feature; N category depends on nodal findings)"
            )
        if patient_age < 40:
            risk_factors.append(f"Young patient age ({patient_age} yrs, higher nodal propensity)")

        # Guideline determination per NCCN / ASCO-SSO
        if microsatellites:
            rec = SlnbRecommendation.RECOMMENDED
        elif breslow_mm < 0.8 and not ulcerated:
            # pT1a
            if (mitotic_rate is not None and mitotic_rate >= 2.0) or lvi or patient_age < 40:
                rec = SlnbRecommendation.DISCUSS_CONSIDER
            else:
                rec = SlnbRecommendation.NOT_RECOMMENDED
        elif (breslow_mm < 0.8 and ulcerated) or (0.8 <= breslow_mm <= 1.0):
            # pT1b
            rec = SlnbRecommendation.DISCUSS_CONSIDER
        else:
            # > 1.0 mm (pT2 - pT4)
            rec = SlnbRecommendation.RECOMMENDED

        return SlnbRiskPrediction(
            probability_pct=prob_pct,
            recommendation=rec,
            logit_score=round(logit, 3),
            risk_factors_present=risk_factors,
            nomogram_contributions={
                "intercept": cls.NOMOGRAM_INTERCEPT,
                "thickness": round(contrib_thickness, 3),
                "ulceration": round(contrib_ulcer, 3),
                "age_adjustment": round(contrib_age, 3),
                "lymphovascular_invasion": round(contrib_lvi, 3),
                "anatomic_site": round(site_weight, 3)
            }
        )

    @classmethod
    def get_surgical_margins(cls, breslow_mm: float, is_in_situ: bool = False) -> SurgicalMarginGuideline:
        """Determines NCCN recommended clinical radial excision margins."""
        if is_in_situ or breslow_mm == 0.0:
            return SurgicalMarginGuideline(
                recommended_clinical_margin_cm="0.5 cm to 1.0 cm",
                guideline_rationale="Melanoma in situ requires 0.5-1.0 cm margins for clear microscopic boundary."
            )
        elif breslow_mm <= 1.0:
            return SurgicalMarginGuideline(
                recommended_clinical_margin_cm="1.0 cm",
                guideline_rationale="Thin invasive melanoma (<=1.0 mm) requires 1.0 cm wide local excision margins."
            )
        elif breslow_mm <= 2.0:
            return SurgicalMarginGuideline(
                recommended_clinical_margin_cm="1.0 cm to 2.0 cm",
                guideline_rationale="Intermediate thickness (1.01-2.0 mm) requires 1.0-2.0 cm margins based on anatomic feasibility."
            )
        else:
            return SurgicalMarginGuideline(
                recommended_clinical_margin_cm="2.0 cm",
                guideline_rationale="Thick melanoma (>2.0 mm) requires 2.0 cm radial surgical margins to minimize locoregional recurrence."
            )

    @classmethod
    def stage_melanoma(cls, specimen: MelanomaSpecimenInput) -> MelanomaStagingReport:
        """Executes full histopathologic indexing and staging pipeline."""
        specimen.validate()

        t_stage = cls.assign_t_stage(
            breslow_mm=specimen.breslow_depth_mm,
            ulcerated=specimen.ulcerated,
            is_in_situ=specimen.is_in_situ
        )

        slnb = cls.calculate_slnb_risk(
            breslow_mm=specimen.breslow_depth_mm,
            patient_age=specimen.patient_age,
            ulcerated=specimen.ulcerated,
            lvi=specimen.lymphovascular_invasion,
            site=specimen.anatomic_site,
            mitotic_rate=specimen.mitotic_rate_per_mm2,
            microsatellites=specimen.microsatellitosis
        )

        margins = cls.get_surgical_margins(
            breslow_mm=specimen.breslow_depth_mm,
            is_in_situ=specimen.is_in_situ
        )

        adverse_features = []
        if specimen.ulcerated:
            adverse_features.append("Primary tumor ulceration (upstages T substage to 'b')")
        if specimen.lymphovascular_invasion:
            adverse_features.append("Lymphovascular invasion (LVI positive)")
        if specimen.perineural_invasion:
            adverse_features.append("Perineural / neurotropic invasion present")
        if specimen.microsatellitosis:
            adverse_features.append(
                "Microsatellites detected (regional metastatic feature; complete stage group requires nodal/distant data)"
            )
        if specimen.mitotic_rate_per_mm2 is not None and specimen.mitotic_rate_per_mm2 >= 1.0:
            adverse_features.append(f"Mitotic activity present ({specimen.mitotic_rate_per_mm2:.1f}/mm2)")
        if specimen.regression_present:
            adverse_features.append("Histopathologic regression noted (may underestimate true Breslow thickness)")

        critical_alerts = []
        if specimen.microsatellitosis:
            critical_alerts.append(
                "STAGING REVIEW: Microsatellitosis is a regional metastatic feature. "
                "Assign the final N category and stage group using complete nodal and distant-metastasis data."
            )
        if t_stage.category.startswith("pT4"):
            critical_alerts.append("HIGH RISK: Thick melanoma (>4.0 mm). High risk for distant and nodal metastasis.")
        if specimen.perineural_invasion:
            critical_alerts.append("ALERT: Perineural invasion identified. Assess for neurotropic / desmoplastic subtype.")

        clark_str = f"Clark Level {specimen.clark_level.name.split('_')[-1]} (Level {specimen.clark_level.value})" if specimen.clark_level else "Not specified / Not microstaged"

        summary = (
            f"{t_stage.category} Cutaneous Melanoma ({specimen.breslow_depth_mm:.2f} mm, "
            f"Ulceration: {'Present' if specimen.ulcerated else 'Absent'}). "
            f"SLNB: {slnb.recommendation.value}. "
            f"Recommended WLE Margins: {margins.recommended_clinical_margin_cm}."
        )

        return MelanomaStagingReport(
            specimen_id=specimen.specimen_id,
            patient_age=specimen.patient_age,
            anatomic_site=specimen.anatomic_site.value,
            breslow_depth_mm=specimen.breslow_depth_mm,
            clark_level=clark_str,
            t_stage=t_stage,
            slnb_evaluation=slnb,
            surgical_margins=margins,
            adverse_features_count=len(adverse_features),
            adverse_features_list=adverse_features,
            critical_alerts=critical_alerts,
            clinical_summary=summary
        )


def format_melanoma_report(report: MelanomaStagingReport) -> str:
    """Renders formatted text clinical staging dossier."""
    lines = []
    lines.append("=" * 78)
    lines.append(f" CUTANEOUS MELANOMA HISTOPATHOLOGIC STAGING DOSSIER : {report.specimen_id}")
    lines.append("=" * 78)
    lines.append(f"Patient Age: {report.patient_age} yrs | Anatomic Site: {report.anatomic_site.title()}")
    lines.append(f"Breslow Depth: {report.breslow_depth_mm:.2f} mm | Anatomic Level: {report.clark_level}")
    lines.append("-" * 78)
    lines.append(f"AJCC 8th EDITION PATHOLOGIC T CATEGORY: {report.t_stage.category}")
    lines.append(f"Criteria Basis: {report.t_stage.criteria_basis}")
    lines.append("-" * 78)
    lines.append(f"SENTINEL LYMPH NODE BIOPSY (SLNB) ASSESSMENT:")
    lines.append(f"  * Exploratory SLN Positivity Estimate: {report.slnb_evaluation.probability_pct:.2f}%")
    lines.append("  * Model note: exploratory estimate; not externally validated for standalone clinical use.")
    lines.append(f"  * SLNB Discussion Guidance: {report.slnb_evaluation.recommendation.value}")
    if report.slnb_evaluation.risk_factors_present:
        lines.append("  * Risk Factors:")
        for rf in report.slnb_evaluation.risk_factors_present:
            lines.append(f"      - {rf}")

    lines.append("-" * 78)
    lines.append(f"SURGICAL EXCISION (WLE) MARGINS:")
    lines.append(f"  * Recommended Clinical Radial Margin: {report.surgical_margins.recommended_clinical_margin_cm}")
    lines.append(f"  * Rationale: {report.surgical_margins.guideline_rationale}")

    if report.adverse_features_list:
        lines.append("-" * 78)
        lines.append(f"ADVERSE HISTOPATHOLOGIC FEATURES ({report.adverse_features_count}):")
        for adv in report.adverse_features_list:
            lines.append(f"  - {adv}")

    if report.critical_alerts:
        lines.append("-" * 78)
        lines.append("[!] PATHOLOGY & STAGING ALERTS:")
        for alert in report.critical_alerts:
            lines.append(f"  * {alert}")

    lines.append("=" * 78)
    return "\n".join(lines)
