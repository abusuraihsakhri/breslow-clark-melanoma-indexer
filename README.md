# Breslow Clark Melanoma Indexer

> **Domain:** Medical Oncology & Cancer Staging Systems  
> **Reference Guidelines & Standards:** `AJCC Cancer Staging Manual & NCCN Clinical Practice Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

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

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`MelanomaSpecimen`** — dedicated module for melanoma specimen evaluation and state verification.
- **`ClarkLevel`** — dedicated module for clark level evaluation and state verification.
- **`AnatomicSite`** — dedicated module for anatomic site evaluation and state verification.
- **`TilCategory`** — dedicated module for til category evaluation and state verification.
- **`SlnbRecommendation`** — dedicated module for slnb recommendation evaluation and state verification.
- **`MelanomaSpecimenInput`**: Input parameters for cutaneous melanoma histopathologic staging.

---

## 📐 Mathematical Formulation & Logic

```text
  NON_BRISK = "Non-brisk"
  BRISK = "Brisk"
  slnb = cls.calculate_slnb_risk(
  risk_score = 0.0
  Calculates VBS incorporating tumor vascularity, lymphovascular invasion,
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --interactive <value> --demo <value> --specimen-id <value> --age <value>
```

### Parameter Reference
- `--interactive`: Specifies input measurement or parameter value.
- `--demo`: Specifies input measurement or parameter value.
- `--specimen-id`: Specifies input measurement or parameter value.
- `--age`: Specifies input measurement or parameter value.
- `--depth`: Specifies input measurement or parameter value.
- `--in-situ`: Specifies input measurement or parameter value.
- `--ulcerated`: Specifies input measurement or parameter value.
- `--clark-level`: Specifies input measurement or parameter value.
- `--mitoses`: Specifies input measurement or parameter value.
- `--lvi`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t breslow-clark-melanoma-indexer .
docker run -p 8000:8000 breslow-clark-melanoma-indexer
```
