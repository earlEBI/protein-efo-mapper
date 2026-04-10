# pQTL Measurement Mapper

Local mapping tool for analyte inputs (proteins/metabolites) and trait inputs (disease/phenotype).

## Quick Start

Set up the local Uvicorn app from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -U pip
python3 -m pip install -r requirements.txt
python3 -m pip install -e .
.venv/bin/analyte-efo-mapper setup-bundled-caches
.venv/bin/python -m uvicorn --app-dir skills/pqtl-measurement-mapper/web app:app --reload --host 127.0.0.1 --port 8000
```

Then open [http://127.0.0.1:8000](http://127.0.0.1:8000).

## Modes

- `map` is for analytes and writes TSV outputs.
- `trait-map` writes TSV outputs only. The repo bundles `efo.obo`, UniProt/metabolite alias caches, trait caches, and the Catalog export snapshot. `setup-bundled-caches` builds the local `measurement_index.json`, and trait-map checks mapped terms against the bundled Catalog export snapshot, prefers exact existing Catalog terms when possible, and writes a Catalog bulk-add TSV when final mapped terms are still missing. A fresher local `catalog_trait_export.tsv` can override the bundled snapshot.

## Notes

- Pyodide browser-only mode is currently hidden in the UI while parity and performance fixes are in progress.
