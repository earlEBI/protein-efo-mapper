# Pyodide Web Mapper Migration Plan

## Goal

Run trait mapping directly in the browser (no server-side subprocess) so the web UI can be deployed as a static app when needed.

## Current State (March 30, 2026)

- Web UI is served by FastAPI at `skills/pqtl-measurement-mapper/web/app.py`.
- Mapping runs by spawning:
  - `map_measurement_efo.py map`
  - `map_measurement_efo.py trait-map`
- This means deployment currently requires Python runtime + local cache files on the server.

## Target State

1. Browser app loads:
- Pyodide runtime
- Trait-mapping core Python module
- Prebuilt cache bundle (JSON)

2. Browser worker performs mapping locally:
- Input parsing
- Cache lookup + lexical ranking
- QC restrictions (`EFO`, `MONDO`, `HP`, `OBA` only)
- TSV output generation

3. Optional server mode remains available:
- For very large runs
- For writeback / heavy cache refresh operations

## Implementation Phases

1. POC scaffold (completed)
- Add `web/pyodide/index.html`, `worker.js`, `mapper_core.py`
- Add cache-bundle builder script:
  `skills/pqtl-measurement-mapper/scripts/build_pyodide_trait_bundle.py`
- Add route in FastAPI app:
  `/pyodide`

2. CLI-compatible runtime mode (completed)
- Add runtime ZIP builder:
  `skills/pqtl-measurement-mapper/scripts/build_pyodide_trait_runtime_bundle.py`
- Add `cli_compat` engine to worker/UI, executing the same
  `map_measurement_efo.py trait-map` script in-browser using bundled references.
- Extend `cli_compat` to analyte `map` mode when runtime bundles include
  `measurement_index.json` (`--include-analyte-assets`).

3. Shared logic extraction (next)
- Move trait scoring / normalization into importable pure module
- Reuse same module in both CLI and Pyodide worker
- Reduce drift between web and CLI behavior

4. Full parity and performance (next)
- Build full browser bundle from `trait_mapping_cache.tsv`
- Chunk/cache by prefix and first-token shard if payload is too large
- Add small regression fixtures and compare Pyodide vs CLI outputs

5. Static deploy option (next)
- Ship static files only (any CDN/object storage)
- Keep server only for optional advanced workflows

## Performance Notes

- Pyodide startup has fixed cost (runtime load).
- Mapping throughput is good once runtime is warm, but initial load is cache-size dependent.
- Best practical path:
  - keep worker alive for session reuse
  - use compact JSON bundle
  - avoid full efo.obo parsing in browser for first iteration

## Validation Strategy

1. Use a targeted trait fixture with known difficult rows.
2. Compare output to current CLI results and curator mappings.
3. Track:
- mapped rate
- clear-worse count
- potentially-worse review count
- ontology-prefix violations (must stay zero outside EFO/MONDO/HP/OBA)

## Known Constraints

- Very large 30k+ rows can still run in browser, but UX may be better with chunked execution.
- Browser version should avoid network lookups (e.g., Zooma) by default for reproducibility.
- Server mode remains useful for long jobs and curated-cache update workflows.
