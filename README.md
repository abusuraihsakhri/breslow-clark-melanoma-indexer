# Breslow Depth & Clark Level Cutaneous Melanoma Indexer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![Standards](https://img.shields.io/badge/Standards-AJCC%208th%20Ed%20%7C%20NCCN%20v2.2024-brightgreen.svg)
![Tests](https://img.shields.io/badge/Tests-Pytest%20Passing-success.svg)

A clinical-grade computational pathology engine for **cutaneous melanoma microstaging**, pathologic T-category assignment (AJCC 8th Edition), anatomical invasion stratification (Clark Levels I–V), quantitative Sentinel Lymph Node Biopsy (SLNB) risk modeling, and Wide Local Excision (WLE) margin guidance.

---

## 🔬 Dermatopathology & Staging Formulations

### 1. AJCC 8th Edition Pathologic T Category (pT)

Cutaneous melanoma primary tumor staging is defined primarily by **Breslow maximal thickness (mm)**, with histologic **ulceration status** serving as the primary sub-staging dichotomy (`a` for non-ulcerated, `b` for ulcerated).

$$\text{pT Category} = f(\text{Breslow Depth } t \text{ [mm]}, \text{Ulceration } u)$$

| AJCC 8th Edition Category | Breslow Thickness ($t$) | Ulceration Status ($u$) | Microstaging Criteria & Historical Notes |
| :--- | :--- | :--- | :--- |
| **Tis** | $0.0\text{ mm}$ (confined to epidermis) | Absent | Melanoma in situ, basement membrane intact |
| **pT1a** | $< 0.8\text{ mm}$ | Absent | Non-ulcerated thin melanoma without dermal ulcer defect |
| **pT1b** | $< 0.8\text{ mm}$ with ulceration **OR** $0.8\text{ mm} \le t \le 1.0\text{ mm}$ | Any ($\pm$ ulceration) | Under AJCC 8th Ed, $0.8-1.0\text{ mm}$ is pT1b regardless of ulceration; mitotic rate removed from T-defining criteria |
| **pT2a** | $1.01\text{ mm} - 2.00\text{ mm}$ | Absent | Intermediate thickness, non-ulcerated |
| **pT2b** | $1.01\text{ mm} - 2.00\text{ mm}$ | Present | Intermediate thickness, ulcerated |
| **pT3a** | $2.01\text{ mm} - 4.00\text{ mm}$ | Absent | Moderately thick, intact epidermis |
| **pT3b** | $2.01\text{ mm} - 4.00\text{ mm}$ | Present | Moderately thick, ulcerated |
| **pT4a** | $> 4.00\text{ mm}$ | Absent | Thick invasive melanoma, non-ulcerated |
| **pT4b** | $> 4.00\text{ mm}$ | Present | Thick invasive melanoma, ulcerated |

> **Note on AJCC 8 Updates:** The mitotic rate (previously defining T1b at $\ge 1/\text{mm}^2$ in AJCC 7th Edition) is still required to be reported by pathology protocol but is no longer a formal T-category defining threshold. Clark level invasion is maintained for anatomical depth documentation and historical comparability.

---

### 2. Clark Level of Anatomical Invasion

Clark levels measure micro-anatomical cutaneous depth relative to skin microarchitecture:

- **Level I:** Intraepidermal / in situ lesions strictly confined above the basement membrane.
- **Level II:** Invasion extending past the dermal-epidermal junction into the loose collagen network of the **papillary dermis**.
- **Level III:** Tumor cells expand and fill the papillary dermis, accumulating at the **papillary-reticular dermis interface**.
- **Level IV:** Penetration of dense collagen fascicles within the **reticular dermis**.
- **Level V:** Transdermal invasion deep into the **subcutaneous adipose tissue (panniculus)**.

---

### 3. Sentinel Lymph Node Biopsy (SLNB) Indications (NCCN / ASCO-SSO)

SLNB evaluation is based on predicted nodal metastasis risk:

1. **pT1a ($<0.8\text{ mm}$, non-ulcerated):** Predicted risk $<5\%$. SLNB is **generally not recommended**, except in high-risk subsets (e.g., age $<40$ years, lymphovascular invasion [LVI], significant mitoses $\ge 2/\text{mm}^2$, or indeterminate margins).
2. **pT1b ($<0.8\text{ mm}$ ulcerated, or $0.8-1.0\text{ mm}$):** Predicted risk $5-10\%$. Clinicians should **discuss and consider** SLNB.
3. **pT2a to pT4b ($>1.0\text{ mm}$):** Predicted risk $>10\%$. SLNB is **routinely recommended** for surgical nodal staging.
4. **Adverse Microstaging Features:** Presence of microscopic satellitosis directly upstages disease to **Stage III (N1c equivalent)**.

$$\text{logit}(p) = \beta_0 + \beta_1 t + \beta_2 u + \beta_3 \left(\frac{\text{Age}-50}{10}\right) + \beta_4 \text{LVI} + \beta_{\text{site}}$$

$$p(\text{SLN Positive}) = \frac{1}{1 + e^{-\text{logit}(p)}}$$

---

### 4. Surgical Wide Local Excision (WLE) Radial Margins

Guideline radial clinical excision margins according to primary tumor thickness:

- **Melanoma in situ (Tis):** $0.5\text{ cm} - 1.0\text{ cm}$ margin.
- **Thin invasive ($\le 1.0\text{ mm}$):** $1.0\text{ cm}$ radial margin.
- **Intermediate thickness ($1.01\text{ mm} - 2.00\text{ mm}$):** $1.0\text{ cm} - 2.0\text{ cm}$ radial margin.
- **Thick melanoma ($> 2.00\text{ mm}$):** $2.0\text{ cm}$ radial margin.

---

## 💻 CLI Quickstart & Usage

The CLI supports interactive single-case entry, benchmark demonstration cases, parameter flags, and batch processing of pathology CSVs.

### Batch Processing CSV
Process a cohort of melanoma pathology records:
```bash
python cli.py batch -i sample.csv -o results.csv
```

Or using standard arguments:
```bash
python cli.py --batch-csv sample.csv --output results.csv
```

### Interactive Mode
Walk through an interactive microstaging consultation:
```bash
python cli.py --interactive
```

### Direct CLI Evaluation
Stage a primary cutaneous melanoma via CLI arguments:
```bash
python cli.py --specimen-id "MEL-2026-09" --depth 1.45 --ulcerated --clark-level 4 --lvi --age 62 --site trunk
```

### Benchmark Demo Scenarios
Run pre-configured clinical benchmark cases (`in_situ`, `t1a`, `t1b`, `t2b`, `t4b`, or `all`):
```bash
python cli.py --demo all
```

---

## 🐍 Python API Quickstart

```python
from breslow_clark_indexer import (
    BreslowClarkMelanomaIndexer,
    MelanomaSpecimenInput,
    ClarkLevel,
    AnatomicSite,
    format_melanoma_report,
)

# 1. Instantiate pathology specimen
specimen = MelanomaSpecimenInput(
    specimen_id="SPEC-PATH-001",
    patient_age=58,
    breslow_depth_mm=1.45,
    ulcerated=True,
    clark_level=ClarkLevel.LEVEL_IV,
    mitotic_rate_per_mm2=3.0,
    lymphovascular_invasion=True,
    anatomic_site=AnatomicSite.HEAD_NECK,
)

# 2. Compute microstaging & clinical recommendations
report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)

# 3. Access structured results
print(f"Pathologic T Stage: {report.t_stage.category}")
print(f"SLNB Recommendation: {report.slnb_evaluation.recommendation.value}")
print(f"SLNB Positivity Risk: {report.slnb_evaluation.probability_pct}%")
print(f"WLE Margin: {report.surgical_margins.recommended_clinical_margin_cm}")

# 4. Render human-readable pathology dossier
print(format_melanoma_report(report))
```

---

## 🧪 Testing

Run the test suite with pytest:
```bash
python -m pytest -p no:zarr
```

All tests execute in isolated environments validating AJCC 8 boundary conditions, Clark level mapping, SLNB risk probabilities, and CLI batch output formatting.
