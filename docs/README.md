# Analyte Mapper (Curator Guide)

This tool maps analytes, traits, ICD10 codes, PheCodes, and UKB data-field references to ontology terms used in GWAS curation.

It runs locally and uses bundled caches so most workflows are offline after setup.

## 1) Quick Start

```bash
mkdir -p ~/src
cd ~/src

git clone https://github.com/earlEBI/analyte-efo-mapper.git
cd analyte-efo-mapper

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install -e .

analyte-efo-mapper setup-bundled-caches
analyte-efo-mapper cache-status --strict --output-json final_output/analyte_mapper_cache_status.json
```

`setup-bundled-caches` reuses the bundled UKB and ICD10 side caches when present, and only rebuilds or downloads missing pieces.  
If local UKB metadata files are missing, setup downloads the official UKB metadata tables automatically.  
If live MONDO refresh fails (for example offline), setup falls back to the local MONDO cache file.  
The compiled `measurement_index.json` is generated locally during setup; it is not required in Git history.

### Fresh Machine Install

For a clean machine starting from GitHub:

```bash
mkdir -p ~/src
cd ~/src

git clone https://github.com/earlEBI/analyte-efo-mapper.git
cd analyte-efo-mapper

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install -e .

analyte-efo-mapper setup-bundled-caches
analyte-efo-mapper cache-status --strict
```

What setup does on a fresh clone:

- uses the bundled core caches already tracked in the repo
- uses the bundled UKB field context, UKB supplement cache, and ICD10 side caches when present
- downloads missing official UKB metadata files when needed
- refreshes or reuses the MONDO ICD10 cache only when the ICD10 supplement must be rebuilt
- builds the ICD10 label cache when missing
- generates the local `skills/pqtl-measurement-mapper/references/measurement_index.json`

After that, normal `map` and `trait-map` runs are local/offline unless you explicitly request online refresh behavior.

Map analytes:

```bash
analyte-efo-mapper map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv
```

Map disease/phenotype traits:

```bash
analyte-efo-mapper trait-map \
  --input data/traits.tsv \
  --output final_output/traits_mapped.tsv \
  --review-output final_output/traits_review.tsv
```

## 2) Input Formats

### Analyte mode (`map`)
Use `.txt`, `.csv`, or `.tsv`.

- Query columns accepted: `query`, `id`, `identifier`, `protein_id`, `metabolite_id`, `metabolite_name`, `gene`, `gene_symbol`, `symbol`, `term`
- Optional type columns: `input_type`, `query_type`, `id_type`, `type`
- Allowed `input_type`: `accession`, `gene_symbol`, `protein_name`, `metabolite_id`, `metabolite_name`, `auto`
- Advanced/conditional: `gene_id` routing exists in the code, but the bundled offline UniProt alias caches do not include enough `gene_ids` to recommend it by default. Only use `gene_id` after explicitly enriching or rebuilding the UniProt alias cache with gene IDs.

Recommended analyte columns:

| Column | Required | Purpose |
|---|---|---|
| `query` | yes | Input analyte identifier or name |
| `input_type` | optional | Explicit routing (`accession`, `gene_symbol`, `protein_name`, `metabolite_id`, `metabolite_name`, `auto`) |

Recommended analyte TSV template:

```tsv
query	input_type
Q9ULI3	accession
P14625	accession
LIPC	gene_symbol
Beta-glucosidase 2	protein_name
HMDB0000122	metabolite_id
CHEBI:38553	metabolite_id
glucose	metabolite_name
```

Bundled example file: `docs/analyte_input_template.tsv`

Create a starter analyte input file:

```bash
mkdir -p data
cat > data/analytes.tsv <<'EOF'
query	input_type
Q9ULI3	accession
LIPC	gene_symbol
glucose	metabolite_name
EOF
```

Analyte context options (`map`):
- `--measurement-context` controls primary matrix/tissue (`blood`, `plasma`, `serum`, `cerebrospinal_fluid`, `urine`, `saliva`, `tissue`, `auto`)
- `--additional-contexts` adds allowed contexts
- `--additional-context-keywords` adds free-text context filters
- `--matrix-priority` sets tie-break order (default: `plasma,blood,serum`)

Example:

```bash
analyte-efo-mapper map \
  --input data/analytes.tsv \
  --output final_output/analytes_efo.tsv \
  --measurement-context blood \
  --additional-contexts serum \
  --additional-context-keywords "aorta,adipose tissue" \
  --matrix-priority plasma,blood,serum
```

### Trait mode (`trait-map`)
Use `.txt`, `.csv`, or `.tsv`.

- Query columns accepted: `query`, `trait`, `reported_trait`, `disease_trait`, `phenotype`, `term`
- Recommended curation columns: `query`, `trait_scale`, `code`, `data_type`
- Optional columns: `icd10`, `phecode`, `input_type`, `trait_scale`, `code`, `data_type`
- Allowed `input_type`: `trait_text`, `icd10`, `phecode`, `ukb_field`, `auto`
- Allowed `trait_scale`: `binary`, `quantitative`, `auto`

Recommended trait columns:

| Column | Required | Purpose |
|---|---|---|
| `query` | yes | Trait text or coded text (for example UKB/ICD10 union-style) |
| `trait_scale` | recommended | `binary` or `quantitative` branch preference |
| `code` | recommended | Source code (`20002`, `Union`, etc.) |
| `data_type` | recommended | Source type (`UKB Data Field`, `ICD10`, etc.) |
| `icd10` | optional | Explicit ICD10 code when available |
| `phecode` | optional | Explicit PheCode when available |
| `input_type` | optional | Manual route (`trait_text`, `icd10`, `phecode`, `ukb_field`, `auto`) |

Recommended trait TSV template:

```tsv
query	trait_scale	code	data_type
20002#1587#aortic regurgitation | incompetence	binary	20002	UKB Data Field
Union#M2571#M25.71 Osteophyte (Shoulder region)	binary	Union	ICD10
41202#M1997#M19.97 Arthrosis | unspecified (Ankle and foot)	binary	41202	UKB Data Field
Hand grip strength (left)	quantitative	46	UKB Data Field
```

Bundled example file: `docs/trait_input_template.tsv`

Create a starter trait input file:

```bash
mkdir -p data
cat > data/traits.tsv <<'EOF'
query	trait_scale	code	data_type
20002#1587#aortic regurgitation | incompetence	binary	20002	UKB Data Field
Hand grip strength (left)	quantitative	46	UKB Data Field
EOF
```

`trait-map` auto-routing now distinguishes:
- ICD10 (including union-like strings such as `Union#A071#A07.1`)
- UKB data field strings (for example `UKB data field 22435`)
- PheCode-like values
- Free text

## 3) Mapping Priority (What Runs First)

### `map` (analytes)
1. Exact analyte cache
2. UniProt / metabolite alias identity resolution
3. Exact ontology label/synonym
4. Lexical/token fallback
5. Review queue for non-validated rows

### `trait-map` (traits)
1. Cache exact ICD10
2. Cache exact ICD10 supplemental
3. Cache exact PheCode
4. Cache exact UKB field ID
5. Cache exact UKB field supplemental
6. Cache exact text (trait / UKB field title)
7. Cache fuzzy text
8. `efo.obo` exact ICD10 xref
9. `efo.obo` exact label/synonym
10. `efo.obo` fuzzy
11. UKB category fallback

## 4) Output (QC-Friendly Columns)

Trait mapping outputs now include:

- Input trace: `input_row_id`, `input_query`, `negation`, `input_type`, `input_trait_scale`, `input_source_code`, `input_source_data_type`, `input_icd10`, `input_icd10_label`, `input_phecode`
- UKB context: `input_ukb_field_id`, `input_ukb_field_label`, `input_ukb_category_id`, `input_ukb_category_label`, `input_ukb_category_path`
- Mapping: `mapped_trait_id`, `mapped_trait_label`, `confidence`, `matched_via`, `matched_on`, `source_file`
- QC flags: `validation`, `qc_review_flag`, `term_not_in_efo`, `mondo_missing_ids`
- Audit text: `evidence`, `provenance`

`term_not_in_efo` is set to `term not in EFO` when a MONDO mapping is not imported in EFO and no better active EFO-family term is available.
`negation` is `yes` when traits contain cues like `none`, `no`, `not`, `without`, or `none of the above`.
Rows with `negation=yes` are always `review_required` if mapped; non-informative options like `none of the above` remain `not_mapped`.
`input_trait_scale` influences branch preference: `quantitative` biases toward measurement mappings; `binary` biases toward non-measurement mappings.

## 5) Cache Refresh From New GWAS Studies TSV

New command:

```bash
analyte-efo-mapper trait-cache-refresh \
  --studies-tsv /path/to/gwas-catalog-studies.tsv
```

Default outputs:
- `final_output/Catalog-mapped-traits_qc.tsv`
- `final_output/Catalog-mapped-traits_high_confidence.tsv`
- `final_output/Catalog-mapped-traits_high_confidence_to_add.tsv`
- `final_output/Catalog-mapped-traits_qc_summary.json`
- `final_output/Catalog-mapped-traits_refresh_state.json` (last processed studies hash)
- Updated trait cache (in place unless `--output-cache` is set)

If the studies TSV hash is unchanged, refresh is skipped automatically (unless `--force-refresh` is set).

Policy used:
- high confidence ontology resolution
- unambiguous mapping per normalized trait
- specific label overlap filter
- non-conflicting with existing curated cache

## 6) Bundled Cache Files

Setup validates and uses:

- `skills/pqtl-measurement-mapper/references/efo_measurement_terms_cache.tsv`
- `skills/pqtl-measurement-mapper/references/analyte_to_efo_cache.tsv`
- `skills/pqtl-measurement-mapper/references/uniprot_aliases.tsv`
- `skills/pqtl-measurement-mapper/references/uniprot_aliases_light.tsv` (generated if `--uniprot-profile light`)
- `skills/pqtl-measurement-mapper/references/metabolite_aliases.tsv` (HMDB-derived alias cache)
- `skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv`
- `skills/pqtl-measurement-mapper/references/efo.obo`
- `references/ukb/fieldsum.txt`
- `references/ukb/field.txt`
- `references/ukb/category.txt`
- `references/ukb/catbrowse.txt`
- `references/ukb/field_trait_supplement_cache.tsv` (generated)
- `references/icd10/mondo.sssom.tsv` (downloaded/refreshed in setup)
- `references/icd10/icd10_label_index.tsv` (generated/downloaded from the official CMS ICD10 release when missing)
- `references/icd10/icd10_trait_supplement_cache.tsv` (generated from EFO ICD10 xrefs plus local MONDO cache)

## 7) What The Final Index Contains

`skills/pqtl-measurement-mapper/references/measurement_index.json` is a generated local artifact produced by `setup-bundled-caches` or `index-build`.

It contains the compiled lookup structures:

- `term_meta`: ontology ID -> label/synonyms/source
- `exact_term_index`: normalized label/synonym -> ontology IDs
- `token_index`: token -> ontology IDs
- `analyte_index`: normalized analyte key -> cached mapping rows
- `accession_alias_index`: UniProt accession -> alias terms
- `symbol_accession_index`: symbol-like alias -> accessions
- `alias_accession_index`: normalized alias phrase -> accessions
- `metabolite_alias_concept_index`: metabolite alias/id -> concept IDs
- `metabolite_concept_index`: concept ID -> labels/aliases/xrefs

The exact JSON bytes do not need to be versioned in Git for a fresh install. A new user can recreate the functional index from the bundled caches plus the setup-generated UKB/ICD10 side caches.

## 8) Setup Audit Manifest

`setup-bundled-caches` now writes a setup inventory JSON (default):

- `final_output/analyte_mapper_cache_manifest.json`

This lists cache file presence/size, cache role, schema/header previews, and runtime index key counts.

You can run a post-setup cache audit any time:

```bash
analyte-efo-mapper cache-status
```

For CI/validation gates:

```bash
analyte-efo-mapper cache-status --strict --output-json final_output/analyte_mapper_cache_status.json
```

## 9) Curator QC Routine

1. Run `trait-map` with `--review-output`.
2. Start with `validation=review_required` rows.
3. Prioritize rows where `term_not_in_efo` is set.
4. Use `matched_on`, `source_file`, `evidence`, and `provenance` for decisions.
5. Keep final curated mappings in trait cache and rerun setup to refresh supplements/index.
