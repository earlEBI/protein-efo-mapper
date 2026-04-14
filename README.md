# analyte-efo-mapper

Local mapper for:

- analyte inputs (proteins and metabolites) to measurement ontology terms
- trait inputs (disease, phenotype, ICD10, PheCode, UKB field context) to EFO, MONDO, HP, and OBA terms

## Quick Start: Uvicorn App

From the repository root:

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

The web app supports:

- `map` for analyte mapping
- `trait-map` for disease and phenotype mapping

Trait runs write TSV outputs only. The bundled core caches include `efo.obo`, UniProt aliases, HMDB/metabolite aliases, trait caches, UKB/ICD10 side caches, the ICD10 alias cache, the MONDO SSSOM cache, and the Catalog export snapshot. `setup-bundled-caches` mainly builds the local `measurement_index.json` from those bundled files and verifies that the quality-preserving ICD10 caches are present. Trait runs write TSV outputs only, check mapped terms against the bundled Catalog export snapshot, prefer exact existing Catalog terms when possible, and write a Catalog bulk-add TSV when final mapped terms are still missing. A fresher local `catalog_trait_export.tsv` can override the bundled snapshot.

## More Docs

- Full curator guide: [docs/README.md](docs/README.md)
- Skill-local notes: [skills/pqtl-measurement-mapper/README.md](skills/pqtl-measurement-mapper/README.md)
