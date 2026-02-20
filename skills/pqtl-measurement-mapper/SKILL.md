---
name: pqtl-measurement-mapper
description: Map protein or metabolite identifiers (including UniProt IDs, metabolite IDs, names, and gene symbols) to measurement trait terms in EFO, including OBA terms imported into EFO. Use when Codex must normalize pQTL/mQTL analyte inputs to controlled measurement EFO IDs and labels with machine-checkable validation status.
---
# pQTL Measurement Mapper

Use this skill to convert analyte-like inputs into measurement ontology terms suitable for GWAS/pQTL metadata curation.

## Runtime

- Use `.venv/bin/python` for bundled scripts.
- The mapper uses only Python standard library modules.
- Use local index mode for scale: build once, map many.

## Input Contract

Provide an input file with one analyte per row:
- `.txt` or `.list`: one value per line
- `.csv` or `.tsv`: include one query column and optionally an input-type column.
  - Query column (one of):
  - `query`, `id`, `identifier`, `protein_id`, `metabolite_id`, `metabolite_name`, `gene`, `gene_symbol`, `symbol`, `term`
- Optional type column (one of): `input_type`, `query_type`, `id_type`, `type`
  - Allowed values: `accession`, `gene_symbol`, `gene_id`, `protein_name`, `metabolite_id`, `metabolite_name`, `auto`
- If no known column exists, the first column is used.

Recommended TSV template:

```tsv
query	input_type
Q9ULI3	accession
LIPC	gene_symbol
Beta-glucosidase 2	protein_name
ENSG00000172137	gene_id
CHEBI:17234	metabolite_id
glucose	metabolite_name
```

Examples of acceptable input values:
- `P02649`
- `APOE`
- `HMDB0000122`
- `apolipoprotein E`

## Output Contract

Write a TSV with one or more candidate mappings per input:
- `input_query`
- `input_type`
- `mapped_efo_id` (underscore format, for example `EFO_0802947`)
- `mapped_label`
- `confidence`
- `matched_via` (`direct-id`, `uniprot`, `synonym`, `none`)
- `validation` (`validated`, `not_validated`, `not_mapped`)
- `evidence`

Preferred downstream use: keep rows with `validation=validated`, then manually review medium-confidence rows.

Cache files used by default:
- `references/analyte_to_efo_cache.tsv` (exact analyte -> EFO mappings)
- `references/efo_measurement_terms_cache.tsv` (full measurement branch list from EFO/OBA imports: id/label/synonyms)
- `references/uniprot_aliases.tsv` (optional accession -> alias expansion table)
- `references/metabolite_aliases.tsv` (optional metabolite concept aliases and IDs: HMDB/ChEBI/KEGG)
- Preferred matrix ordering for tie-breaks: `plasma,blood,serum`

## Procedure

1. Build index once
- Create a JSON index from local cache/reference TSV files.
- Include term exact-match and token indexes for fast retrieval.
- Include analyte cache and optional UniProt accession alias index.
- Include optional metabolite concept alias index (HMDB/ChEBI/KEGG + synonyms).
- Build term cache from full EFO measurement branch (recommended), not sparse bootstrap queries.
- Generate `uniprot_aliases.tsv` from UniProt export when needed.
- For external users, preload all human UniProt aliases once to avoid alias-missing failures on new accessions.

2. Map queries in bulk
- Read the input list and map each query from the local index.
- Run exact analyte cache hits first.
- Expand UniProt accessions using local alias table when available.
- Resolve metabolite IDs/names using local metabolite alias concept index when available.
- In `--name-mode strict` (default), free-text names/gene symbols are mapped only after exact alias resolution to UniProt accessions.
- For mixed inputs, keep `--entity-type auto`; for single-entity batches, use `--entity-type protein` or `--entity-type metabolite`.
- Apply exact term and token-retrieval lexical scoring.
- Support compound ID cells (for example `Q9NZ08;Q9NZ08-2`) and isoform canonicalization.
- For new IDs, optionally auto-enrich UniProt aliases before mapping.
- Enforce matrix/tissue context (`blood` by default); explicit non-blood terms (for example cerebrospinal fluid) are rejected unless context is changed.
- In default `blood` context, serum-specific terms are excluded unless explicitly added with `--additional-contexts serum`.
- Optional free-text context add-on is available via `--additional-context-keywords` (for example `aorta`).

3. Validate and review
- Mark `validated` when mapped ID exists in local term index.
- For UniProt accession queries, enforce accession identity consistency (gene/protein aliases and numeric-family consistency) before accepting a candidate.
- Review medium/low-confidence rows.

4. Write back high-confidence hits (optional)
- Append validated mappings to analyte cache for future exact matches.

## Bundled Script

Build local index:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/build_efo_measurement_cache.py \
  --output skills/pqtl-measurement-mapper/references/efo_measurement_terms_cache.tsv
```

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/build_uniprot_aliases.py \
  --download-human \
  --download-to /tmp/uniprot_human_export.tsv \
  --output skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv \
  --source uniprot-human-live
```

Alternative from a local UniProt export file:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/build_uniprot_aliases.py \
  --input data/uniprot_export.tsv \
  --output skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv \
  --source uniprot-export-20260219
```

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  index-build \
  --term-cache skills/pqtl-measurement-mapper/references/efo_measurement_terms_cache.tsv \
  --analyte-cache skills/pqtl-measurement-mapper/references/analyte_to_efo_cache.tsv \
  --uniprot-aliases skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv \
  --metabolite-aliases skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --output-index skills/pqtl-measurement-mapper/references/measurement_index.json
```

Build metabolite alias table from HMDB/ChEBI/KEGG resources (optional, recommended for metabolite mapping):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  metabolite-alias-build \
  --output skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --download-common-sources
```

Optional HMDB download in the same command:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  metabolite-alias-build \
  --output skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --download-common-sources \
  --download-hmdb
```

Refresh EFO measurement cache and rebuild index (one command):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  refresh-efo-cache \
  --download-url https://github.com/EBISPOT/efo/releases/latest/download/efo.obo \
  --term-cache skills/pqtl-measurement-mapper/references/efo_measurement_terms_cache.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --analyte-cache skills/pqtl-measurement-mapper/references/analyte_to_efo_cache.tsv \
  --uniprot-aliases skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv \
  --metabolite-aliases skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv
```

By default this refresh command downloads EFO OBO from the URL above when `--efo-obo` is not provided.  
Index rebuild is enabled by default and should usually be left on.

Map in bulk:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  map \
  --input data/analytes.tsv \
  --output data/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --entity-type auto \
  --top-k 2 \
  --min-score 0.55 \
  --workers 8 \
  --name-mode strict \
  --progress \
  --measurement-context blood \
  --matrix-priority plasma,blood,serum \
  --cache-writeback
```

Allow serum alongside blood (optional):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  map \
  --input data/analytes.tsv \
  --output data/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --measurement-context blood \
  --additional-contexts serum
```

Free-text context keyword filter (optional):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  map \
  --input data/analytes.tsv \
  --output data/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --measurement-context auto \
  --additional-context-keywords aorta
```

Performance tips for large batches:
- Keep `--auto-enrich-uniprot` off during the main map pass unless you specifically need live alias fetches.
- Use `--parallel-mode process --workers 8` on machines where multiprocessing is allowed.
- If process mode is restricted in your environment, the script falls back to thread mode automatically; in that case `--workers 1` is often fastest for CPU-bound mapping.
- Keep `--name-mode strict` for production to avoid semantic drift on free-text protein names. Use `--name-mode fuzzy` only as a fallback review mode.

Strict high-precision mode for new terms (recommended):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --name-mode strict \
  --auto-enrich-uniprot \
  --uniprot-aliases skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv \
  --term-cache skills/pqtl-measurement-mapper/references/efo_measurement_terms_cache.tsv \
  --unmapped-output final_output/analytes_unmapped.tsv \
  --measurement-context blood \
  --matrix-priority plasma,blood,serum \
  --progress \
  --workers 8
```

Optional fallback mode (returns a row for every input, but unresolved rows are `not_validated`):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py \
  map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --auto-enrich-uniprot \
  --force-map-best \
  --fallback-efo-id EFO:0001444
```

LLM-assisted fallback mode (for unresolved rows only, OpenAI-compatible API):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/llm_fallback_mapper.py \
  --mapped-input final_output/analytes_efo.tsv \
  --output final_output/analytes_efo_llm.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --measurement-context blood \
  --api-base https://api.openai.com/v1 \
  --model gpt-4.1-mini \
  --api-key-env OPENAI_API_KEY \
  --max-candidates 15 \
  --workers 4
```

Web app mode (upload + progress + downloads):

```bash
.venv/bin/pip install -r requirements.txt
```

```bash
.venv/bin/python -m uvicorn \
  --app-dir skills/pqtl-measurement-mapper/web \
  app:app \
  --host 127.0.0.1 \
  --port 8000
```

Then open `http://127.0.0.1:8000` in a browser.

Notes:
- Works with any platform exposing an OpenAI-compatible `chat/completions` API.
- LLM output is accepted only when selected ID is from local shortlist and passes strict identity post-validation.

## References

- Mapping and adjudication rubric: `references/mapping-guidelines.md`
