# Breslow–Clark Melanoma Indexer

### [Open the Live Application →](https://abusuraihsakhri.github.io/breslow-clark-melanoma-indexer/)

A compact cutaneous melanoma microstaging tool with a Python CLI/API and a browser interface for GitHub Pages.

## What it does

The maintained workflow assigns an AJCC 8th Edition pathologic T category from Breslow thickness and ulceration, records Clark level as a historical/anatomic descriptor, provides general wide-local-excision margin guidance, and returns an SLNB discussion prompt. A numeric SLN-positivity estimate is retained for backwards compatibility, but it is an exploratory model and is not externally validated for standalone clinical use.

This project does **not** determine a complete melanoma stage group. Nodal status, distant metastasis, subtype, anatomic site, pathology context, and current clinical guidelines still need independent assessment.

## Browser application

The repository includes a static, client-side interface designed for GitHub Pages. All calculations run locally in the browser; entered values are not transmitted to a server by the application.

The browser interface intentionally excludes mucosal melanoma because this repository implements cutaneous melanoma staging logic.

## Python usage

Run an example case:

```bash
python cli.py --specimen-id MEL-001 --depth 1.45 --ulcerated --clark-level 4 --age 62 --site trunk
```

Run the interactive mode:

```bash
python cli.py --interactive
```

Process a CSV cohort:

```bash
python cli.py batch -i sample.csv -o results.csv
```

Use the Python API:

```python
from breslow_clark_indexer import (
    AnatomicSite,
    BreslowClarkMelanomaIndexer,
    MelanomaSpecimenInput,
)

specimen = MelanomaSpecimenInput(
    specimen_id="SPEC-001",
    patient_age=58,
    breslow_depth_mm=1.45,
    ulcerated=True,
    anatomic_site=AnatomicSite.TRUNK,
)

report = BreslowClarkMelanomaIndexer.stage_melanoma(specimen)
print(report.t_stage.category)
print(report.surgical_margins.recommended_clinical_margin_cm)
```

## Validation and scope

- Negative Breslow depth and negative mitotic rate are rejected.
- An in-situ flag cannot be combined with a positive Breslow depth.
- Mucosal melanoma is rejected because it is outside the cutaneous staging scope.
- Empty batch CSV files return an error instead of raising an indexing exception.
- Clark level is documented but is not used to define AJCC 8 T category.
- The SLN probability value is an exploratory estimate, not a validated nomogram implementation.

## Testing

The GitHub Actions workflow runs the Python test suite on Python 3.10, 3.11, and 3.12, performs a CLI batch smoke test, checks Python compilation, and runs browser-model tests with Node.js.

Locally:

```bash
python -m unittest discover -s tests -v
python cli.py batch -i sample.csv -o out_smoke.csv
node tests/test_web_model.mjs
```

## Technology

- Python standard library for the CLI/API
- Vanilla HTML, CSS, and JavaScript for the browser application
- GitHub Actions for continuous integration and Pages deployment
- No runtime third-party Python dependencies

## Privacy

The browser app performs calculations entirely client-side and does not send form data to an application backend. The repository is public, so do not commit patient data, credentials, or other sensitive information.

## License

MIT License. See [LICENSE](LICENSE).
