# Pyodide Trait Mapper

This directory contains a browser-only trait-mapping app:

- `index.html`: UI
- `worker.js`: Web Worker that loads Pyodide
- `mapper_core.py`: fast cache-mode Python mapper logic
- `cache_sample.json`: small starter cache bundle (fast mode)

Modes:

- `fast`: lightweight cache JSON + simplified logic (quick startup, trait mapping only)
- `cli_compat`: runs the same `map_measurement_efo.py` code path in-browser using a runtime ZIP bundle
  - supports `trait-map`
  - supports `map` (analyte) when runtime bundle includes analyte assets

## Build a larger cache bundle

From repo root:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/build_pyodide_trait_bundle.py \
  --trait-cache skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv \
  --output skills/pqtl-measurement-mapper/web/pyodide/cache_bundle.json \
  --limit 0
```

Then in the prototype UI, choose `cache_bundle.json` as the optional override bundle.

## Build CLI-Compatible Runtime Bundle

From repo root:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/build_pyodide_trait_runtime_bundle.py \
  --repo-root . \
  --output skills/pqtl-measurement-mapper/web/pyodide/trait_runtime_bundle.zip
```

This ZIP includes:

- `map_measurement_efo.py`
- `trait_mapping_cache.tsv`
- `efo.obo`
- UKB and ICD10 supplementary caches used by trait-map
- manual curator seed/override caches

The UI will auto-load `trait_runtime_bundle.zip` when `cli_compat` mode is selected (or you can upload a custom ZIP).

To include analyte map assets (`measurement_index.json`) for protein/metabolite mapping:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/build_pyodide_trait_runtime_bundle.py \
  --repo-root . \
  --include-analyte-assets \
  --output skills/pqtl-measurement-mapper/web/pyodide/analyte_runtime_bundle.zip
```

Note: analyte-enabled runtime bundles are substantially larger.

## Parity QC Against UKB Master Mapping

Use the QC helper to compare manual mappings vs CLI trait-map vs Pyodide mapper:

```bash
.venv/bin/python qc/ukb_master_pyodide_cli_parity.py \
  --manual-input ./path/to/UK_Biobank_master_file.tsv \
  --run-cli \
  --pyodide-engine cli_compat
```

Outputs default to `./tmp/ukb_master_cli_pyodide_manual_{compare,focus,summary}.*`.
