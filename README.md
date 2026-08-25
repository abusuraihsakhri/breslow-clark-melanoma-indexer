# Breslow Depth & Clark Level Cutaneous Melanoma Indexer

A deterministic, guideline-calibrated computational pathology and clinical decision-support engine implementing the **AJCC 8th Edition Pathologic T Category (pT)** rules, **Clark Level** anatomical invasion indexing, **NCCN / ASCO-SSO** Sentinel Lymph Node Biopsy (SLNB) indications, and Wide Local Excision (WLE) surgical margin standards.

---

## Clinical Staging Architecture

Histopathologic microstaging of primary cutaneous melanoma provides essential prognostic stratification and determines surgical and nodal management.

### AJCC 8th Edition Pathologic T Category (pT) Criteria

| Category | Breslow Thickness (mm) | Ulceration Status |
| :--- | :--- | :--- |
| **Tis** | In situ ($0.0\text{ mm}$) | Not applicable |
| **pT1a** | $< 0.8\text{ mm}$ | Without ulceration |
| **pT1b** | $< 0.8\text{ mm}$ with ulceration, OR $0.8 - 1.0\text{ mm}$ | With or without ulceration |
| **pT2a** | $> 1.0 - 2.0\text{ mm}$ | Without ulceration |
| **pT2b** | $> 1.0 - 2.0\text{ mm}$ | With ulceration |
| **pT3a** | $> 2.0 - 4.0\text{ mm}$ | Without ulceration |
| **pT3b** | $> 2.0 - 4.0\text{ mm}$ | With ulceration |
| **pT4a** | $> 4.0\text{ mm}$ | Without ulceration |
| **pT4b** | $> 4.0\text{ mm}$ | With ulceration |

### Clark Level of Anatomical Invasion

- **Level I**: Intraepidermal (Melanoma in situ, above basement membrane)
- **Level II**: Invasive into the papillary dermis
- **Level III**: Fills and expands the papillary dermis to the papillary-reticular dermal junction
- **Level IV**: Invades into the reticular dermis
- **Level V**: Invades into the subcutaneous adipose tissue (fat)

### Sentinel Lymph Node Biopsy (SLNB) Guidance & Nomogram

- **$\text{pT1a } (<0.8\text{ mm non-ulcerated})$**: SLNB generally not recommended ($<5\%$ positive risk) unless high-risk adverse features are present (elevated mitoses $\ge 2/\text{mm}^2$, LVI, age $<40$).
- **$\text{pT1b } (<0.8\text{ mm ulcerated or } 0.8-1.0\text{ mm})$**: Discuss and consider SLNB ($5-10\%$ risk).
- **$\text{pT2 - pT4 } (>1.0\text{ mm})$**: SLNB recommended ($>10\%$ risk).
- **Calibrated Logistic Risk Nomogram**:
  $$\text{Logit}(p) = -2.80 + 0.55 \times \text{Thickness} + 0.90 \times \text{Ulceration} - 0.30 \times \left(\frac{\text{Age}-50}{10}\right) + 0.80 \times \text{LVI} + \text{Site Factor}$$
  $$P(\text{SLN}^+) = \frac{1}{1 + e^{-\text{Logit}}}$$

### Recommended Wide Local Excision (WLE) Margins

- **Melanoma in situ (Tis)**: $0.5\text{ cm} - 1.0\text{ cm}$
- **$\le 1.0\text{ mm}$**: $1.0\text{ cm}$
- **$1.01 - 2.0\text{ mm}$**: $1.0\text{ cm} - 2.0\text{ cm}$
- **$> 2.0\text{ mm}$**: $2.0\text{ cm}$

---

## Installation & Setup

Requires **Python 3.9+** (zero external dependencies, pure standard library).

```bash
git clone https://github.com/abusuraihsakhri/breslow-clark-melanoma-indexer.git
cd breslow-clark-melanoma-indexer
```

---

## CLI Usage Examples

### 1. Pre-Configured Benchmark Scenarios

```bash
python cli.py --demo t1a
python cli.py --demo t1b
python cli.py --demo t2b
python cli.py --demo t4b
```

### 2. Direct Argument Staging with JSON Output

```bash
python cli.py --specimen-id MEL-2026-092 --age 62 --depth 1.85 \
  --ulcerated --clark-level 4 --lvi --site head_neck --json
```

### 3. Batch CSV Evaluation

```bash
python cli.py --batch-csv input_cases.csv --output staged_results.csv
```

### 4. Interactive Staging Mode

```bash
python cli.py --interactive
```

---

## Python API Integration

```python
from breslow_clark_indexer import (
    MelanomaSpecimenInput,
    ClarkLevel,
    AnatomicSite,
    BreslowClarkMelanomaIndexer,
    format_melanoma_report,
)

specimen = MelanomaSpecimenInput(
    specimen_id="PAT-MEL-402",
    patient_age=58,
    breslow_depth_mm=1.45,
    ulcerated=True,
    clark_level=ClarkLevel.LEVEL_IV,
    lymphovascular_invasion=False,
    anatomic_site=AnatomicSite.TRUNK
)

report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
print(format_melanoma_report(report))
```

---

## Unit Testing

Run the automated test suite with 27 unit test cases:

```bash
python -m unittest test_breslow_clark_indexer.py -v
```

---

## License

MIT License. Authored and maintained by Dr. Abu Suraih Sakhri.
