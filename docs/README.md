# Analyte EFO Mapper Docs

This folder is a static GitHub Pages site.

## For New Users

The website is documentation and templates only.  
Actual mapping runs locally with Python.

Quick local setup:

```bash
git clone https://github.com/earlEBI/analyte-efo-mapper.git
cd analyte-efo-mapper
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install -e .
cp docs/input_template.tsv data/analytes.tsv
analyte-efo-mapper setup-bundled-caches
```

By default, setup now builds and uses a lightweight UniProt cache for indexing
(reviewed-like accessions, gene symbol required) at:
`skills/pqtl-measurement-mapper/references/uniprot_aliases_light.tsv`

To force setup to use the full bundled UniProt cache instead:
`analyte-efo-mapper setup-bundled-caches --uniprot-profile full`

Bundled offline caches used by setup:
- `skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv` (protein alias cache)
- `skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv` (metabolite HMDB-derived alias cache)
- `skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv` (disease/phenotype trait cache)
- `references/ukb/fieldsum.txt` (UKB data-field title catalog)
- `references/ukb/field.txt` (UKB field metadata including main category IDs)
- `references/ukb/category.txt` (UKB category labels)
- `references/ukb/catbrowse.txt` (UKB category parent/child tree)

These are local repository files; setup does not download these caches from the internet.
Trait mapping also needs `skills/pqtl-measurement-mapper/references/efo.obo`. If missing, `setup-bundled-caches`
and `trait-map` auto-provision it from the bundled pinned URL (default release asset:
`efo.obo.gz` in this repo).
If you maintain your own fork/release asset, override with:
`--efo-obo-bundled-url https://github.com/<owner>/<repo>/releases/latest/download/efo.obo.gz`

Cache provenance and refresh/build path:
- `uniprot_aliases.tsv`: bundled source alias cache in repo; setup builds `uniprot_aliases_light.tsv` by default for faster mapping.
- `uniprot_aliases_light.tsv`: generated locally by `setup-bundled-caches` (or `uniprot-alias-build-light`), optional online enrichment via `uniprot-alias-enrich`.
- `metabolite_aliases.tsv`: bundled HMDB-derived alias cache in repo; can be rebuilt from your pinned HMDB SDF via `metabolite-alias-build`.
- `trait_mapping_cache.tsv`: bundled curated disease/phenotype cache in repo.
- `references/ukb/fieldsum.txt`: bundled UKB field title catalog used by `trait-map`.
- `references/ukb/field.txt`: bundled UKB field metadata used by `trait-map` for category-aware matching.
- `references/ukb/category.txt`: bundled UKB category label catalog used by `trait-map` when category IDs/titles appear in input.
- `references/ukb/catbrowse.txt`: bundled UKB category tree used to resolve category-path context in `trait-map`.
- `efo_measurement_terms_cache.tsv`: bundled measurement-term cache in repo; can be rebuilt from OBO via `refresh-efo-cache`.
- `measurement_index.json`: local compiled index built from caches by `setup-bundled-caches` or `index-build`.
- `efo.obo`: local ontology file used by `trait-map` for fallback and obsolete-ID checks.
- `efo_measurement_terms_cache.tsv` and `efo.obo` are intentionally separate:
  - measurement cache is optimized for protein/metabolite mapping runtime
  - full `efo.obo` is used for trait fallback + obsolete-term validation

Optional (recommended if you plan to use `input_type=gene_id` heavily): backfill UniProt gene IDs into the lightweight cache during setup:

```bash
analyte-efo-mapper setup-bundled-caches --uniprot-light-enrich-gene-ids
```

Or run enrichment explicitly:

```bash
analyte-efo-mapper uniprot-alias-enrich \
  --uniprot-aliases skills/pqtl-measurement-mapper/references/uniprot_aliases_light.tsv \
  --workers 8
```

If you want to rebuild the lightweight cache manually:

```bash
analyte-efo-mapper uniprot-alias-build-light \
  --source skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv \
  --output skills/pqtl-measurement-mapper/references/uniprot_aliases_light.tsv
```

CLI command and upgrade:

- Run commands via:
  - `analyte-efo-mapper ...`
- Upgrade after new GitHub changes:
  - `git pull --rebase origin main`
  - `python -m pip install -e .`
  - `analyte-efo-mapper setup-bundled-caches`

Quick Start (5 lines):

```bash
cp docs/input_template.tsv data/analytes.tsv
analyte-efo-mapper map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --withheld-triage-output final_output/withheld_for_review_triage.tsv
```

This writes your main results to `final_output/analytes_efo.tsv` and triaged withheld candidates to `final_output/withheld_for_review_triage.tsv`.
`setup-bundled-caches` is offline/local by default and validates bundled caches, including bundled UKB field/category dictionaries.
If you add `--uniprot-light-enrich-gene-ids`, setup will also perform optional online UniProt backfill.

Disease/phenotype trait mapping mode (optional):

```bash
analyte-efo-mapper trait-map \
  --input data/traits.tsv \
  --output final_output/traits_mapped.tsv \
  --trait-cache skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv \
  --ukb-field-catalog references/ukb/fieldsum.txt \
  --efo-obo skills/pqtl-measurement-mapper/references/efo.obo \
  --min-score 0.82 \
  --review-output final_output/traits_review.tsv \
  --stream-output \
  --flush-every 10 \
  --memoize-queries \
  --progress
```

Notes:
- This mode is cache-first (curated trait cache), then falls back to `efo.obo`.
- Input can include free-text trait, ICD10, and/or PheCode columns.
- UKB-aware trait context uses bundled `references/ukb/fieldsum.txt`, `field.txt`, `category.txt`, and `catbrowse.txt` when available.
- Output includes a per-row `provenance` field suitable for curation notes tabs.
- Before mapping, trait cache IDs are checked against `efo.obo`; obsolete IDs are remapped via `replaced_by` (or single `consider`) where possible, and unresolved IDs are warned and excluded from auto output.

How mapping works (order and method):

- Protein/metabolite `map` order:
  - Normalize query + infer/resolve input type.
  - Try deterministic identity resolution first:
    - proteins: UniProt accession, gene symbol/gene ID/mnemonic -> accession
    - metabolites: HMDB/ChEBI/KEGG IDs, concept aliases
  - Retrieve candidate EFO/OBA measurement terms (exact label/synonym and indexed lexical/token retrieval).
  - Rerank candidates by lexical score + identity/context signals.
  - Validate candidate identity/context gates; emit validated rows.
  - If best candidate is plausible but blocked by validation, write to withheld/review outputs (if enabled).

- Trait `trait-map` order:
  - Exact cache hit by ICD10.
  - Exact cache hit by PheCode.
  - Exact cache hit by reported trait text.
  - Cache fuzzy text (only when no ontology-exact match exists).
  - Exact `efo.obo` label/synonym match.
  - Fuzzy `efo.obo` fallback.
  - Before all cache-based trait mapping, cache ontology IDs are checked against `efo.obo`; obsolete IDs are remapped (via `replaced_by` or single `consider`) or excluded.

Validation:

- Identity validation:
  - proteins: accession/alias subject checks prevent near-name drift (including number mismatches).
  - metabolites: HMDB/ChEBI/KEGG concept consistency checks prevent concept collisions.
- Context validation:
  - default `--measurement-context blood`: allows blood/plasma/unlabeled; excludes serum unless explicitly added.
  - `--additional-contexts` and `--additional-context-keywords` can expand allowed matrices/tissues.
- Structural safeguards:
  - ratio terms are blocked unless the query indicates ratio/composite intent.
  - lower-confidence lexical-only candidates are withheld from auto-validated output.

Withheld triage flow:

- Enable with `--withheld-triage-output final_output/withheld_for_review_triage.tsv`.
- Mapper writes non-auto-validated top candidates there with triage status:
  - `reject_high`: likely wrong target.
  - `review_needed`: unresolved identity.
  - `accept_medium`: moderate support.
  - `accept_high`: strong support.
- Triage includes query label + suggested mapped label/subject to support fast manual QC.

Practical expectation:

- Accession/ID inputs are most reliable.
- Gene symbol/gene ID are reliable when resolution to a unique accession succeeds.
- Free-text names are useful but should be treated as review-first for edge cases.

Run mapping (mixed protein + metabolite inputs):

```bash
analyte-efo-mapper map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --entity-type auto \
  --workers 8 \
  --parallel-mode process \
  --progress
```

You can run this directly after clone/install.  
No local HMDB file is required for normal mapping because the repo ships a metabolite alias cache.

Optional (recommended for metabolite mapping): pin your HMDB file version, then build cache once:

```bash
mkdir -p skills/pqtl-measurement-mapper/references/hmdb_source
cp /path/to/your/structures.sdf skills/pqtl-measurement-mapper/references/hmdb_source/structures.sdf
```

Then build aliases (auto-detects local HMDB file):

```bash
analyte-efo-mapper metabolite-alias-build \
  --output skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --no-merge-existing
```

If you already downloaded HMDB locally (for example `structures.sdf`), build from local HMDB only:

```bash
analyte-efo-mapper metabolite-alias-build \
  --output skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --hmdb-xml /path/to/structures.sdf \
  --no-merge-existing
```

HMDB download/cache notes:
- The command auto-detects HMDB files in `skills/pqtl-measurement-mapper/references/hmdb_source/` and `skills/pqtl-measurement-mapper/references/metabolite_downloads/`.
- If no local HMDB file is found, it auto-downloads from `--hmdb-url` (default: this repo's GitHub release asset).
- You can point `--hmdb-url` to your own private/public GitHub release bundle.
- Maintainer setup: upload your pinned `structures.sdf.gz` as a GitHub Release asset so users can fetch the exact same HMDB version.
- Use `--no-merge-existing` to regenerate a cleaner alias cache (normalized IDs and fewer noisy aliases).
- Practical model:
  - End users: no HMDB step needed unless rebuilding aliases.
  - Rebuild step: pulls your pinned bundle from GitHub release automatically when local HMDB is absent.
  - Fully offline installs: include `structures.sdf` or `structures.sdf.gz` in `references/hmdb_source/`.

Example with explicit pinned bundle URL:

```bash
analyte-efo-mapper metabolite-alias-build \
  --output skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --hmdb-url https://github.com/<owner>/<repo>/releases/download/<tag>/structures.sdf.gz \
  --no-merge-existing
```

Optional (recommended after upgrading): refresh EFO measurement cache so metabolite xrefs
(CHEBI/HMDB/KEGG) are available for direct ID-to-term matching:

```bash
analyte-efo-mapper refresh-efo-cache \
  --term-cache skills/pqtl-measurement-mapper/references/efo_measurement_terms_cache.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json
```

Then map with it:

```bash
analyte-efo-mapper map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --metabolite-aliases skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --entity-type auto
```

Optional triaged withheld output:

```bash
analyte-efo-mapper map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --withheld-triage-output final_output/withheld_for_review_triage.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --entity-type auto \
  --workers 8 \
  --parallel-mode process \
  --progress
```

`final_output/withheld_for_review_triage.tsv` classifies top withheld candidates as:

- `reject_high`: likely wrong target.
- `review_needed`: unresolved identity, send to manual review.
- `accept_medium`: accession mismatch but same primary symbol (often alias/isoform-equivalent).
- `accept_high`: candidate subject resolves to same accession.
- Includes `query_primary_label` and `suggested_subject` so you can QC full-name alignment directly.
- For metabolite rows, triage uses metabolite ID/xref overlap (HMDB/ChEBI/KEGG) and emits statuses such as `supports_query_metabolite_id`, `supports_query_metabolite_concept`, and `conflicts_with_query_metabolite`.

Optional broader review queue output:

```bash
analyte-efo-mapper map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --review-queue-output final_output/review_queue.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --entity-type auto \
  --workers 8 \
  --parallel-mode process \
  --progress
```

`final_output/review_queue.tsv` is for manual review of closest unresolved candidates; it is useful but can contain probable mismatches.

When `--review-output` is not set, triage mode generates an internal temporary withheld-review file automatically.

Refresh EFO/OBA measurement cache and rebuild index:

```bash
analyte-efo-mapper refresh-efo-cache \
  --download-url https://github.com/EBISPOT/efo/releases/latest/download/efo.obo \
  --term-cache skills/pqtl-measurement-mapper/references/efo_measurement_terms_cache.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --analyte-cache skills/pqtl-measurement-mapper/references/analyte_to_efo_cache.tsv \
  --uniprot-aliases skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv \
  --metabolite-aliases skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv
```

Notes:

- You do not need a local OBO file for this command.
- Index rebuild is enabled by default and should usually be kept on.

Defaults and options:

- Default `--measurement-context blood` excludes serum-specific terms.
- Add serum explicitly with `--additional-contexts serum`.
Measurement context usage:

- Blood default (blood + plasma + unlabeled, excludes serum):
  - `--measurement-context blood`
- Plasma-only:
  - `--measurement-context plasma`
- Blood plus serum:
  - `--measurement-context blood --additional-contexts serum`
- CSF-only:
  - `--measurement-context cerebrospinal_fluid`
- Free-text context add-on:
  - `--additional-context-keywords aorta`
  - `--additional-context-keywords "adipose tissue"`
- Free-text only filter mode:
  - `--measurement-context auto --additional-context-keywords aorta`

Useful optional arguments:

- `--top-k` number of candidates to keep per input.
- `--min-score` minimum score threshold.
- `--workers` parallel workers.
- `--name-mode strict|fuzzy` handling for name-like inputs.
- `--entity-type auto|protein|metabolite` routing/validation mode.
- `--auto-enrich-uniprot` checks UniProt for accession-like queries missing from your local alias file, appends returned aliases to the file set by `--uniprot-aliases`, and rebuilds index before mapping.
  - It is incremental for current input queries, not a full UniProt rebuild.
  - If you are using the lightweight setup profile, pass `--uniprot-aliases skills/pqtl-measurement-mapper/references/uniprot_aliases_light.tsv`.
- `setup-bundled-caches` default UniProt profile is `light` (reviewed-like + requires gene symbol) for faster index build.
- `setup-bundled-caches --uniprot-profile full` uses the full bundled UniProt alias cache.
- `uniprot-alias-build-light` creates/refreshes a lightweight UniProt cache from a source alias TSV.
- `uniprot-alias-enrich` backfills missing UniProt gene IDs (and aliases) across a chosen UniProt alias cache file (for example `uniprot_aliases_light.tsv`).
  - Use `--refresh-all` if you want to refresh all rows, not only rows with empty `gene_ids`.
- `--metabolite-aliases` local metabolite concept aliases for HMDB/ChEBI/KEGG/name resolution.
- `--unmapped-output` write unresolved rows to a separate TSV.
- `--review-queue-output` write optional broader manual-review suggestions (`review_queue.tsv`).
- `--withheld-triage-output` triage top withheld candidates into accept/reject/review buckets.

The published site URL is:

`https://earlebi.github.io/analyte-efo-mapper/`
