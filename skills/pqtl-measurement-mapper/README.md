# pQTL Measurement Mapper

Local mapping tool for analyte inputs (proteins/metabolites) and trait inputs (disease/phenotype).

## Quick Start

Set up the local Uvicorn app from the repository root:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m uvicorn --app-dir skills/pqtl-measurement-mapper/web app:app --reload --host 127.0.0.1 --port 8000
```

Then open [http://127.0.0.1:8000](http://127.0.0.1:8000).

If your checkout uses `skills-src` instead of `skills`, run:

```bash
.venv/bin/python -m uvicorn --app-dir skills-src/pqtl-measurement-mapper/web app:app --reload --host 127.0.0.1 --port 8000
```

## Modes

- `map` is for analytes and writes TSV outputs.
- `trait-map` is for disease/phenotype traits and writes TSV outputs plus a Catalog bulk-add TSV when final mapped terms are still missing from the local Catalog export.

## Notes

- Pyodide browser-only mode is currently hidden in the UI while parity and performance fixes are in progress.
