# pQTL/mQTL Measurement Mapping Guidelines

Use this reference when adjudicating automatic mappings from IDs/gene symbols to EFO measurement terms.

## Scope

Map each input to a *measurement* trait term in EFO (including OBA terms imported into EFO), not to:
- disease labels
- pathways
- tissues
- generic proteins/metabolites without a measurement context

## Preferred Inputs

- UniProt accessions (e.g., `P02649`)
- HGNC gene symbols (e.g., `APOE`)
- metabolite IDs (e.g., HMDB/CHEBI style IDs)
- protein/metabolite names

## Resolution Strategy

1. Keep the raw input as a candidate synonym.
2. Expand UniProt accessions using local alias table (`uniprot_aliases.tsv`) when available.
3. Match against local EFO term cache (`efo_measurement_terms_cache.tsv`) using:
- exact normalized synonym/label matches
- token retrieval + lexical reranking
4. Prefer terms where label/synonyms contain measurement language:
- `measurement`, `level`, `concentration`, `abundance`, `circulating`, `plasma`, `serum`
5. Apply matrix priority tie-break with this default order:
- `plasma` > `blood` > `serum`
6. Enforce measurement context:
- default context is `blood`
- reject candidates explicitly tagged to conflicting matrices (for example `cerebrospinal fluid`) unless context is explicitly changed

## Hybrid Cache Files

- `analyte_to_efo_cache.tsv`: curated exact-key mappings used before lexical search.
- `efo_measurement_terms_cache.tsv`: full EFO measurement branch list (`id`, `label`, `synonyms`) used for lexical ranking and validation.
- `uniprot_aliases.tsv`: local accession alias expansion table (recommend preloading all human UniProt entries).
- `measurement_index.json`: compiled index produced by `index-build`.

Build `efo_measurement_terms_cache.tsv` from full EFO OBO measurement subtree with:
- `scripts/build_efo_measurement_cache.py --output references/efo_measurement_terms_cache.tsv`

Keep both files tab-delimited with headers.

## Validation Rules

A mapping is `validated` only when:
- the candidate EFO/OBA ID exists in local term index (`measurement_index.json`).

If no indexed ID match exists, mark as `not_validated` even if lexical similarity is high.

## Optional LLM Fallback

Use LLM fallback only for rows still marked `not_mapped` after deterministic mapping.

Safety gates:
- LLM must choose from a precomputed local shortlist.
- Chosen ID must pass strict accession identity post-validation.
- If either gate fails, keep row as `not_mapped`.

## Confidence Interpretation

- `>= 0.80`: high confidence, usually acceptable if validated
- `0.60-0.79`: medium confidence, review top 2 candidates
- `< 0.60`: low confidence, manual review required

## Manual Review Checklist

- Confirm the mapped term is a measurement trait relevant to pQTL/mQTL analyses.
- Confirm analyte identity (protein/metabolite) is preserved.
- Reject mappings that collapse to disease endpoints.
- If multiple biological matrices exist (plasma/serum/CSF), keep the matrix-matched term when available.
