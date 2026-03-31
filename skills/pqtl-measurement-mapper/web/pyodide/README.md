# Pyodide Trait Mapper Prototype

This directory contains a browser-only trait-mapping prototype:

- `index.html`: UI
- `worker.js`: Web Worker that loads Pyodide
- `mapper_core.py`: browser-safe Python mapper logic
- `cache_sample.json`: small starter cache bundle

## Build a larger cache bundle

From repo root:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/build_pyodide_trait_bundle.py \
  --trait-cache skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv \
  --output skills/pqtl-measurement-mapper/web/pyodide/cache_bundle.json \
  --limit 0
```

Then in the prototype UI, choose `cache_bundle.json` as the optional override bundle.
