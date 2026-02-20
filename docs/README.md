# GitHub Pages for pQTL Mapper

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
cp docs/input_template.tsv data/analytes.tsv
```

Run mapping (mixed protein + metabolite inputs):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --entity-type auto \
  --workers 8 \
  --parallel-mode process \
  --progress
```

Optional: build a metabolite alias resource from online sources first (no local source files required):

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py metabolite-alias-build \
  --output skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --download-common-sources
```

Optional: also attempt HMDB download in the same command:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py metabolite-alias-build \
  --output skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --download-common-sources \
  --download-hmdb
```

If HMDB download fails, rerun without `--download-hmdb` (ChEBI + KEGG still build), or provide a local `--hmdb-xml` file.

Then map with it:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --metabolite-aliases skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv \
  --entity-type auto
```

Optional review queue output:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --review-output final_output/review_queue.tsv \
  --index skills/pqtl-measurement-mapper/references/measurement_index.json \
  --entity-type auto \
  --workers 8 \
  --parallel-mode process \
  --progress
```

`final_output/review_queue.tsv` is for manual review of closest unresolved candidates; it is useful but can contain probable mismatches.

Refresh EFO/OBA measurement cache and rebuild index:

```bash
.venv/bin/python skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py refresh-efo-cache \
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
- Output IDs are written with underscores (for example `EFO_0802947`).

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
- `--auto-enrich-uniprot` checks UniProt for accession-like queries missing from your local alias file, appends returned aliases to `uniprot_aliases.tsv`, and rebuilds index before mapping.
  - It is incremental for current input queries, not a full UniProt rebuild.
- `--metabolite-aliases` local metabolite concept aliases for HMDB/ChEBI/KEGG/name resolution.
- `--unmapped-output` write unresolved rows to a separate TSV.
- `--review-output` write optional manual-review suggestions (`review_queue.tsv`).

The published site URL is:

`https://earlebi.github.io/analyte-efo-mapper/`
