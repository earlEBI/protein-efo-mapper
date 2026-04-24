# pQTL Measurement Mapper

This Codex skill maps analyte inputs (proteins or metabolites) and trait inputs
(disease or phenotype) to validated EFO or OBA terms for downstream GWAS
curation.

Use the bundled mapper when you are already inside an `analyte-efo-mapper`
checkout. If this skill is installed elsewhere, locate or bootstrap an
`analyte-efo-mapper` installation first and use its bundled resources.

## Expected Setup

Before first use, ensure the mapper resources are prepared:

```bash
.venv/bin/analyte-efo-mapper setup-bundled-caches
```

In GWAS curation workflows, `gwas-curation/scripts/preflight_check.py` and
`gwas-curation/scripts/ensure_pqtl_mapper.py` can provision or locate the
dependency automatically.

## Quick Start

If you want the local app, run it from the `analyte-efo-mapper` repository
root:

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
- `trait-map` is for disease or phenotype inputs and writes TSV outputs. It
  uses the bundled Catalog export snapshot to check whether mapped terms
  already exist in the Catalog DB, prefers exact existing Catalog terms when
  possible, and writes a Catalog bulk-add TSV when final mapped terms are still
  missing. A fresher local `catalog_trait_export.tsv` can override the bundled
  snapshot.

## Notes

- Pyodide browser-only mode is currently hidden in the UI while parity and performance fixes are in progress.
