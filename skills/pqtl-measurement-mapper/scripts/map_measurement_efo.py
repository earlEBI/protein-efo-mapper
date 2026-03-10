#!/usr/bin/env python3
"""Bulk/offline mapper for analyte/trait queries -> ontology IDs.

Architecture:
- setup-bundled-caches: validate/use bundled local caches and build local index (offline-first)
- cache-status: audit presence/missing state for setup caches/index
- icd10-label-cache-build: build an official ICD10 label cache from CMS tabular-order releases
- index-build: build a local JSON index from cache/reference files
- map: map protein/metabolite query files using local index only (fast, scalable)
- trait-map: map disease/phenotype traits via curated trait cache + efo.obo fallback
"""

from __future__ import annotations

import argparse
import csv
import difflib
import gzip
import hashlib
import html
import json
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.error
import urllib.request
import time
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, UTC
from functools import lru_cache
from pathlib import Path
from typing import Any

SEARCHABLE_COLUMNS = [
    "query",
    "id",
    "identifier",
    "protein_id",
    "metabolite_id",
    "metabolite_name",
    "gene",
    "gene_symbol",
    "symbol",
    "term",
]
INPUT_TYPE_COLUMNS = [
    "input_type",
    "query_type",
    "id_type",
    "type",
]
INPUT_TYPE_ALIASES = {
    "": "auto",
    "auto": "auto",
    "any": "auto",
    "accession": "accession",
    "uniprot": "accession",
    "uniprot_id": "accession",
    "uniprot_accession": "accession",
    "protein_id": "accession",
    "gene_symbol": "gene_symbol",
    "symbol": "gene_symbol",
    "gene": "gene_symbol",
    "hgnc_symbol": "gene_symbol",
    "gene_id": "gene_id",
    "ensembl_gene_id": "gene_id",
    "entrez_gene_id": "gene_id",
    "protein_name": "protein_name",
    "name": "protein_name",
    "protein": "protein_name",
    "analyte_name": "protein_name",
    "metabolite_id": "metabolite_id",
    "hmdb_id": "metabolite_id",
    "chebi_id": "metabolite_id",
    "kegg_id": "metabolite_id",
    "metabolite_name": "metabolite_name",
    "compound_name": "metabolite_name",
    "compound": "metabolite_name",
}
ENTITY_TYPE_ALIASES = {
    "": "auto",
    "auto": "auto",
    "any": "auto",
    "both": "auto",
    "protein": "protein",
    "proteins": "protein",
    "metabolite": "metabolite",
    "metabolites": "metabolite",
    "compound": "metabolite",
    "compounds": "metabolite",
}

MEASUREMENT_KEYWORDS = {
    "measurement",
    "level",
    "levels",
    "concentration",
    "abundance",
    "amount",
    "assay",
    "protein",
    "metabolite",
    "circulating",
    "blood",
    "plasma",
    "serum",
    "csf",
}
GENERIC_TOKENS = MEASUREMENT_KEYWORDS | {"trait", "quantification", "quantitative"}
GENERIC_NAME_QUERY_TOKENS = {
    "protein",
    "subunit",
    "member",
    "family",
    "class",
    "type",
    "candidate",
    "readthrough",
    "alpha",
    "beta",
    "gamma",
    "small",
    "open",
    "reading",
    "frame",
    "chromosome",
    "chain",
    "receptor",
    "nmd",
}
IDENTITY_STOP_TOKENS = GENERIC_TOKENS | {
    "protein",
    "proteins",
    "candidate",
    "readthrough",
    "nmd",
    "open",
    "reading",
    "frame",
    "chromosome",
}
IDENTITY_AMBIGUOUS_TOKENS = {
    "and",
    "or",
    "of",
    "in",
    "protein",
    "proteins",
    "molecule",
    "group",
    "family",
    "member",
    "subunit",
    "chain",
    "domain",
    "type",
    "factor",
    "receptor",
    "binding",
    "related",
    "activated",
    "kinase",
    "dependent",
    "associated",
    "like",
    "antigen",
    "glycoprotein",
    "containing",
    "isozyme",
    "atp",
    "putative",
    "probable",
    "homolog",
    "isoform",
    "component",
    "alpha",
    "beta",
    "gamma",
    "delta",
    "epsilon",
    "kappa",
    "lambda",
}
IDENTITY_RELATION_TOKENS = {
    "binding",
    "interacting",
    "receptor",
    "inhibitor",
    "activator",
    "associated",
    "ligand",
}
TYPE_SPECIFIER_TOKENS = {
    "platelet",
    "muscle",
    "liver",
    "hepatic",
    "cardiac",
    "heart",
    "brain",
    "neuronal",
    "kidney",
    "renal",
    "pancreatic",
    "endothelial",
    "epithelial",
    "erythroid",
    "myeloid",
}
SUBJECT_PHRASE_EQUIVALENTS: tuple[tuple[str, str], ...] = (
    # Common immunology shorthand used in UniProt aliases vs EFO/OBA labels.
    (r"\bmbl[- ]associated\b", "mannan binding lectin"),
    (r"\bmbl\b", "mannan binding lectin"),
    # Common enzyme/protein shorthand.
    (r"\bpaf\b", "platelet activating factor"),
)
LIPID_PANEL_PREFIX_RE = re.compile(
    r"^\s*(?:main parameters|lipoprotein main fractions|(?:vldl|ldl|hdl)\s+subfractions)\s*,\s*",
    re.IGNORECASE,
)
LIPID_PANEL_REPLACEMENTS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"\bapo[-\s]?a1\b", re.IGNORECASE), "apolipoprotein a 1"),
    (re.compile(r"\bapo[-\s]?a2\b", re.IGNORECASE), "apolipoprotein a 2"),
    # Handle common OCR/Unicode confusables for ApoB (B/8/beta/Cyrillic Ve).
    (re.compile(r"\bapo[-\s]?(?:8|β|Β|в|В|ß)\b", re.IGNORECASE), "apolipoprotein b"),
    (re.compile(r"\bapolipoprotein\s+a[-\s]?1\b", re.IGNORECASE), "apolipoprotein a-i"),
    (re.compile(r"\bapolipoprotein\s+a[-\s]?2\b", re.IGNORECASE), "apolipoprotein a-ii"),
    (re.compile(r"\bapo[-\s]?b100\b", re.IGNORECASE), "apolipoprotein b-100"),
    (re.compile(r"\bapo[-\s]?b\b", re.IGNORECASE), "apolipoprotein b"),
    (re.compile(r"\bapob\b", re.IGNORECASE), "apolipoprotein b"),
    (
        re.compile(
            r"\b(apolipoprotein\s+(?:a[-\s]?(?:i{1,3}|[12])|b(?:[-\s]?100)?))\s+in\s+(?:all\s+)?(?:(?:very\s+small|small|medium|large|very\s+large|extremely\s+large)\s+)?(?:vldl|ldl|hdl|idl)\b",
            re.IGNORECASE,
        ),
        r"\1",
    ),
    (
        re.compile(
            r"\b(?:total\s+)?apolipoprotein\s+b(?:[-\s]?100)?\s+particle\s+(?:number|concentration|count|no\.?)\b",
            re.IGNORECASE,
        ),
        "apolipoprotein b measurement",
    ),
    (re.compile(r"\bcholesterol esters\b", re.IGNORECASE), "cholesteryl esters"),
    (re.compile(r"\bce\b", re.IGNORECASE), "cholesteryl esters"),
    (re.compile(r"\bch\b", re.IGNORECASE), "cholesterol"),
    (re.compile(r"\btriglycerides\b", re.IGNORECASE), "triglyceride"),
    (re.compile(r"\bmufa\b", re.IGNORECASE), "monounsaturated fatty acids"),
    (re.compile(r"\bpufa\b", re.IGNORECASE), "polyunsaturated fatty acids"),
    (re.compile(r"\bufa\b", re.IGNORECASE), "unsaturated fatty acids"),
    (re.compile(r"\bsfa\b", re.IGNORECASE), "saturated fatty acids"),
    (
        re.compile(r"\bmonounsaturated fatty acid\s*/\s*total fatty acid(?:s)?(?:\s*\(%\))?\b", re.IGNORECASE),
        "monounsaturated fatty acids to total fatty acids percentage",
    ),
    (
        re.compile(r"\bpolyunsaturated fatty acid\s*/\s*total fatty acid(?:s)?(?:\s*\(%\))?\b", re.IGNORECASE),
        "polyunsaturated fatty acids to total fatty acids percentage",
    ),
    (
        re.compile(r"\bsaturated fatty acid\s*/\s*total fatty acid(?:s)?(?:\s*\(%\))?\b", re.IGNORECASE),
        "saturated fatty acids to total fatty acids percentage",
    ),
    (
        re.compile(r"\bunsaturated fatty acid\s*/\s*total fatty acid(?:s)?(?:\s*\(%\))?\b", re.IGNORECASE),
        "unsaturated fatty acids to total fatty acids percentage",
    ),
    (
        re.compile(r"\bmufa\s*/\s*pufa\b", re.IGNORECASE),
        "monounsaturated fatty acids to polyunsaturated fatty acids ratio",
    ),
    (
        re.compile(r"\bufa\s*/\s*pufa\b", re.IGNORECASE),
        "unsaturated fatty acids to polyunsaturated fatty acids ratio",
    ),
    (re.compile(r"\bacetoacetic acid\b", re.IGNORECASE), "acetoacetate"),
    (re.compile(r"\bformic acid\b", re.IGNORECASE), "formate"),
    (re.compile(r"\bcitric acid\b", re.IGNORECASE), "citric acid/isocitric acid"),
    (re.compile(r"\bdihydrothymine\b", re.IGNORECASE), "5,6-dihydrothymine"),
    (re.compile(r"\bacetic acid\b", re.IGNORECASE), "acetate"),
    (re.compile(r"\bin all (vldl|ldl|hdl|idl)\b", re.IGNORECASE), r"in \1"),
    (re.compile(r"\bextremely large (vldl|ldl|hdl)\b", re.IGNORECASE), r"very large \1"),
    (re.compile(r"\bparticle number\b", re.IGNORECASE), "particle concentration"),
    (re.compile(r"\ball all\b", re.IGNORECASE), "all"),
    (re.compile(r"\bn-acetylglycoproteins[12]\b", re.IGNORECASE), "n-acetylglycoproteins"),
    (re.compile(r"\bmeasurement levels?\b", re.IGNORECASE), "measurement"),
    (re.compile(r"\blevels?\s+measurement\b", re.IGNORECASE), "measurement"),
    (re.compile(r"\bmeasurements\b", re.IGNORECASE), "measurement"),
)
LIPID_MEASUREMENT_VARIANT_RE = re.compile(
    r"^(?:"
    r"apolipoprotein a 1|apolipoprotein a 2|apolipoprotein b(?:-100)?|"
    r"cholesterol|free cholesterol|cholesteryl ester(?:s)?|phospholipid(?:s)?|triglyceride(?:s)?|"
    r"polyunsaturated fatty acid(?:s)?|monounsaturated fatty acid(?:s)?|unsaturated fatty acid(?:s)?|saturated fatty acid(?:s)?"
    r")(?: in (?:vldl|ldl|hdl|idl))?$",
    re.IGNORECASE,
)
PROTEIN_CUE_RE = re.compile(
    r"\b(?:apolipoprotein|apo[-\s]?[ab][0-9]{0,3}|protein|receptor|subunit|chain)\b",
    re.IGNORECASE,
)
METABOLITE_CUE_RE = re.compile(
    r"\b(?:"
    r"cholesterol|cholesteryl|triglyceride(?:s)?|phospholipid(?:s)?|"
    r"fatty acid(?:s)?|acetoacetate|acetoacetic|formate|formic|citric|dihydrothymine|"
    r"hdl|ldl|vldl|idl|chylomicron(?:s)?|mufa|pufa|ufa|sfa|"
    r"glucose|lactate|pyruvate|ketone"
    r")\b",
    re.IGNORECASE,
)
METABOLITE_CONCEPT_ALIAS_EXCLUDE = {
    "cholesteryl ester",
    "cholesteryl esters",
    "cholesteryl ester measurement",
    "cholesteryl esters measurement",
}
ROMAN_NUMERAL_TOKENS = frozenset(
    {
        "i",
        "ii",
        "iii",
        "iv",
        "v",
        "vi",
        "vii",
        "viii",
        "ix",
        "x",
        "xi",
        "xii",
        "xiii",
        "xiv",
        "xv",
    }
)

# UniProt accession formats:
# - 6-char: [OPQ][0-9][A-Z0-9]{3}[0-9]
# - 6/10-char: [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
# In UniProt, reviewed (Swiss-Prot) primaries are typically in the OPQ 6-char
# namespace. We use this as a local heuristic because the bundled alias cache
# does not carry an explicit reviewed/unreviewed flag per accession.
UNIPROT_REVIEWED_LIKE_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
ENSEMBL_GENE_ID_RE = re.compile(r"^ENSG\d{6,}(?:\.\d+)?$", re.IGNORECASE)
ENTREZ_GENE_ID_RE = re.compile(r"^(?:NCBIGENE:|GENEID:|ENTREZ:)?(\d{2,9})$", re.IGNORECASE)
HMDB_ID_RE = re.compile(r"^HMDB(\d{3,})$", re.IGNORECASE)
CHEBI_ID_RE = re.compile(r"^CHEBI[:_](\d{1,7})$", re.IGNORECASE)
KEGG_ID_RE = re.compile(r"^(?:KEGG[:_])?(?:CPD:|DR:)?([CD]\d{5})$", re.IGNORECASE)
EFO_RE = re.compile(r"^(EFO|OBA)[:_]\d+$", re.IGNORECASE)
ALLOWED_ONTOLOGY_ID_RE = re.compile(r"^(EFO|OBA)[:_]\d+$", re.IGNORECASE)
MONDO_CURIE_RE = re.compile(r"MONDO[:_]\d+", re.IGNORECASE)
ICD10_CODE_HEAD_RE = r"[A-TV-Z][0-9][0-9A-Z]"
ICD10_XREF_RE = re.compile(
    rf"(?:ICD10(?:CM|WHO)?)[_:]({ICD10_CODE_HEAD_RE}(?:\.[0-9A-Z]{{1,4}})?)",
    re.IGNORECASE,
)
ICD10_XREF_URL_RE = re.compile(
    rf"/ICD10(?:CM|WHO)?/({ICD10_CODE_HEAD_RE}(?:\.[0-9A-Z]{{1,4}})?)",
    re.IGNORECASE,
)
MULTI_ID_SEP_RE = re.compile(r"[;|]")
HTML_TAG_RE = re.compile(r"<[^>]+>")
METABOLITE_NOISE_RE = re.compile(r"^[\[\](){},.:;+*/\\|=_-]+$")
METABOLITE_CHEBI_INLINE_RE = re.compile(r"CHEBI[:_](\\d{1,7})", re.IGNORECASE)
METABOLITE_HMDB_INLINE_RE = re.compile(r"HMDB[:_]?([0-9]{3,})", re.IGNORECASE)
METABOLITE_KEGG_INLINE_RE = re.compile(r"(?:CPD:|COMPOUND:|KEGG[:_])([CD]\\d{5})", re.IGNORECASE)
SCRIPT_RUNTIME_PATH = Path(globals().get("SCRIPT_PATH", __file__)).resolve()
SKILL_DIR = SCRIPT_RUNTIME_PATH.parents[1]
REPO_ROOT = SCRIPT_RUNTIME_PATH.parents[3]

DEFAULT_ANALYTE_CACHE = SKILL_DIR / "references" / "analyte_to_efo_cache.tsv"
DEFAULT_TERM_CACHE = SKILL_DIR / "references" / "efo_measurement_terms_cache.tsv"
DEFAULT_INDEX = SKILL_DIR / "references" / "measurement_index.json"
DEFAULT_UNIPROT_ALIASES = SKILL_DIR / "references" / "uniprot_aliases.tsv"
DEFAULT_UNIPROT_ALIASES_LIGHT = SKILL_DIR / "references" / "uniprot_aliases_light.tsv"
DEFAULT_METABOLITE_ALIASES = SKILL_DIR / "references" / "metabolite_aliases.tsv"
DEFAULT_METABOLITE_DOWNLOAD_DIR = SKILL_DIR / "references" / "metabolite_downloads"
DEFAULT_HMDB_LOCAL_DIR = SKILL_DIR / "references" / "hmdb_source"
DEFAULT_HMDB_LOCAL_FILENAMES = (
    "structures.sdf",
    "structures.sdf.gz",
    "hmdb_metabolites.xml",
    "hmdb_metabolites.xml.gz",
)
DEFAULT_HMDB_XML_URL = "https://github.com/earlEBI/analyte-efo-mapper/releases/latest/download/structures.sdf.gz"
DEFAULT_CHEBI_NAMES_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/names.tsv.gz"
DEFAULT_KEGG_CONV_URL = "https://rest.kegg.jp/conv/chebi/compound"
DEFAULT_EFO_OBO_URL = "https://github.com/EBISPOT/efo/releases/latest/download/efo.obo"
DEFAULT_EFO_OBO_LOCAL = SKILL_DIR / "references" / "efo.obo"
DEFAULT_EFO_OBO_BUNDLED_URL = DEFAULT_EFO_OBO_URL
DEFAULT_TRAIT_CACHE = SKILL_DIR / "references" / "trait_mapping_cache.tsv"
LEGACY_TRAIT_CACHE = REPO_ROOT / "final_output" / "code-EFO-mappings_final_mapping_cache.tsv"
DEFAULT_UKB_FIELD_CATALOG = REPO_ROOT / "references" / "ukb" / "fieldsum.txt"
DEFAULT_UKB_FIELD_METADATA = REPO_ROOT / "references" / "ukb" / "field.txt"
DEFAULT_UKB_CATEGORY_CATALOG = REPO_ROOT / "references" / "ukb" / "category.txt"
DEFAULT_UKB_CATEGORY_TREE = REPO_ROOT / "references" / "ukb" / "catbrowse.txt"
DEFAULT_UKB_FIELD_CATALOG_URL = "https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=16"
DEFAULT_UKB_FIELD_METADATA_URL = "https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=1"
DEFAULT_UKB_CATEGORY_CATALOG_URL = "https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=3"
DEFAULT_UKB_CATEGORY_TREE_URL = "https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=13"
DEFAULT_UKB_FIELD_SUPPLEMENT_CACHE = REPO_ROOT / "references" / "ukb" / "field_trait_supplement_cache.tsv"
DEFAULT_ICD10_SUPPLEMENT_CACHE = REPO_ROOT / "references" / "icd10" / "icd10_trait_supplement_cache.tsv"
DEFAULT_ICD10_LABEL_CACHE = REPO_ROOT / "references" / "icd10" / "icd10_label_index.tsv"
DEFAULT_MONDO_SSSOM_CACHE = REPO_ROOT / "references" / "icd10" / "mondo.sssom.tsv"
DEFAULT_ICD10_LABEL_SOURCE_URL = "https://www.cms.gov/files/zip/2026-code-descriptions-tabular-order.zip"
DEFAULT_MONDO_SSSOM_URL = "https://purl.obolibrary.org/obo/mondo/mappings/mondo.sssom.tsv"
DEFAULT_SETUP_CACHE_MANIFEST = REPO_ROOT / "final_output" / "analyte_mapper_cache_manifest.json"
DEFAULT_SETUP_CACHE_STATUS = REPO_ROOT / "final_output" / "analyte_mapper_cache_status.json"
DEFAULT_CATALOG_QC_OUTPUT = REPO_ROOT / "final_output" / "Catalog-mapped-traits_qc.tsv"
DEFAULT_CATALOG_HIGH_CONF_OUTPUT = REPO_ROOT / "final_output" / "Catalog-mapped-traits_high_confidence.tsv"
DEFAULT_CATALOG_TO_ADD_OUTPUT = REPO_ROOT / "final_output" / "Catalog-mapped-traits_high_confidence_to_add.tsv"
DEFAULT_CATALOG_SUMMARY_OUTPUT = REPO_ROOT / "final_output" / "Catalog-mapped-traits_qc_summary.json"
DEFAULT_CATALOG_REFRESH_STATE_OUTPUT = REPO_ROOT / "final_output" / "Catalog-mapped-traits_refresh_state.json"
DEFAULT_TRAIT_QC_RISK_SUFFIX = "_qc_risk.tsv"
DEFAULT_MEASUREMENT_ROOT = "EFO:0001444"
DEFAULT_MATRIX_PRIORITY = ["plasma", "blood", "serum"]
DEFAULT_MEASUREMENT_CONTEXT = "blood"
CACHE_REQUIREMENT_BY_NAME = {
    "index": "required",
    "term_cache": "required",
    "analyte_cache": "required",
    "uniprot_aliases_source": "required",
    "uniprot_aliases_index": "required",
    "metabolite_aliases": "required",
    "trait_cache": "required",
    "ukb_field_catalog": "required",
    "ukb_field_metadata": "required",
    "ukb_category_catalog": "required",
    "ukb_category_tree": "required",
    "ukb_field_supplement_cache": "required",
    "icd10_supplement_cache": "required",
    "icd10_label_cache": "recommended",
    "efo_obo": "required",
    "mondo_sssom_cache": "recommended",
    "hmdb_source_dir": "optional",
    "metabolite_download_dir": "optional",
}
CACHE_ROLE_BY_NAME = {
    "index": "runtime_index",
    "term_cache": "runtime_source",
    "analyte_cache": "runtime_source",
    "uniprot_aliases_source": "rebuild_source",
    "uniprot_aliases_index": "derived_runtime_source",
    "metabolite_aliases": "runtime_source",
    "trait_cache": "runtime_source",
    "ukb_field_catalog": "supporting_metadata",
    "ukb_field_metadata": "supporting_metadata",
    "ukb_category_catalog": "supporting_metadata",
    "ukb_category_tree": "supporting_metadata",
    "ukb_field_supplement_cache": "derived_supporting_cache",
    "icd10_supplement_cache": "derived_supporting_cache",
    "icd10_label_cache": "supporting_cache",
    "efo_obo": "ontology_source",
    "mondo_sssom_cache": "refresh_cache",
    "hmdb_source_dir": "source_workspace",
    "metabolite_download_dir": "download_workspace",
}
CACHE_DESCRIPTION_BY_NAME = {
    "index": "Primary runtime JSON used by mapper lookups.",
    "term_cache": "Measurement term TSV used to build the runtime index.",
    "analyte_cache": "Curated exact analyte-to-EFO TSV mappings.",
    "uniprot_aliases_source": "Full UniProt alias TSV kept for reproducible rebuilds.",
    "uniprot_aliases_index": "Lightweight UniProt alias TSV used during index build.",
    "metabolite_aliases": "Metabolite concept and alias TSV used during index build.",
    "trait_cache": "Trait mapping TSV used by disease/phenotype workflows.",
    "ukb_field_catalog": "UK Biobank field catalog listing available fields.",
    "ukb_field_metadata": "UK Biobank field metadata with titles and categories.",
    "ukb_category_catalog": "UK Biobank category catalog.",
    "ukb_category_tree": "UK Biobank category hierarchy.",
    "ukb_field_supplement_cache": "Derived UKB trait supplement cache built during setup.",
    "icd10_supplement_cache": "Derived ICD10 trait supplement cache built during setup.",
    "icd10_label_cache": "Optional ICD10 label lookup cache.",
    "efo_obo": "Local EFO ontology fallback used for trait resolution.",
    "mondo_sssom_cache": "Optional MONDO SSSOM cache used when refreshing ICD10 mappings.",
    "hmdb_source_dir": "Optional HMDB source workspace for metabolite refreshes.",
    "metabolite_download_dir": "Optional download workspace for metabolite refreshes.",
}
# Token-retrieved candidates are noisier than exact/synonym matches. Keep them
# out of auto-validated output unless confidence remains reasonably high.
TOKEN_SUBJECT_EXACT_AUTO_VALIDATE_MIN_SCORE = 0.70
TOKEN_SEMANTIC_EXACT_AUTO_VALIDATE_MIN_SCORE = 0.60
NEEDS_NEW_TERM_MATCH_TAG = "needs-new-term"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
HTTP_USER_AGENT = "pqtl-measurement-mapper/1.0"
MEASUREMENT_CONTEXT_PATTERNS = {
    "blood": ["blood", "circulating"],
    "plasma": ["plasma", "blood plasma"],
    "serum": ["serum", "blood serum"],
    "cerebrospinal_fluid": ["cerebrospinal fluid", "csf"],
    "urine": ["urine", "urinary"],
    "saliva": ["saliva", "salivary"],
    "tissue": ["tissue"],
}
MEASUREMENT_CONTEXT_ALIASES = {
    "auto": "auto",
    "any": "auto",
    "none": "auto",
    "blood": "blood",
    "whole blood": "blood",
    "plasma": "plasma",
    "blood plasma": "plasma",
    "serum": "serum",
    "blood serum": "serum",
    "circulating": "blood",
    "csf": "cerebrospinal_fluid",
    "cerebrospinal fluid": "cerebrospinal_fluid",
    "cerebrospinal_fluid": "cerebrospinal_fluid",
    "urine": "urine",
    "urinary": "urine",
    "saliva": "saliva",
    "salivary": "saliva",
    "tissue": "tissue",
}

TRAIT_QUERY_COLUMNS = [
    "query",
    "trait",
    "reported_trait",
    "disease_trait",
    "phenotype",
    "term",
]
TRAIT_ICD10_COLUMNS = [
    "icd10",
    "icd_10",
    "icd",
    "diagnosis_code",
]
TRAIT_PHECODE_COLUMNS = [
    "phecode",
    "phe_code",
    "phe",
]
TRAIT_INPUT_TYPE_COLUMNS = [
    "input_type",
    "query_type",
    "id_type",
    "type",
]
TRAIT_CODE_COLUMNS = [
    "code",
    "trait_code",
    "field_code",
    "field_id",
    "ukb_field_id",
    "icd10_code",
]
TRAIT_DATA_TYPE_COLUMNS = [
    "data_type",
    "code_type",
    "source_type",
]
TRAIT_SCALE_COLUMNS = [
    "trait_scale",
    "scale",
    "binary_quantitative",
    "binary/quantitative",
    "binary / quantitative",
    "outcome_scale",
    "outcome_type",
]
TRAIT_INPUT_TYPE_ALIASES = {
    "": "auto",
    "auto": "auto",
    "any": "auto",
    "text": "trait_text",
    "trait": "trait_text",
    "trait_text": "trait_text",
    "reported_trait": "trait_text",
    "disease_trait": "trait_text",
    "phenotype": "trait_text",
    "icd10": "icd10",
    "icd": "icd10",
    "phecode": "phecode",
    "phe": "phecode",
    "ukb_field": "ukb_field",
    "ukb data field": "ukb_field",
    "ukb-data-field": "ukb_field",
}
TRAIT_SCALE_ALIASES = {
    "": "auto",
    "auto": "auto",
    "any": "auto",
    "unspecified": "auto",
    "unknown": "auto",
    "binary": "binary",
    "case_control": "binary",
    "case-control": "binary",
    "casecontrol": "binary",
    "disease": "binary",
    "qualitative": "binary",
    "categorical": "binary",
    "quantitative": "quantitative",
    "continuous": "quantitative",
    "measurement": "quantitative",
    "numeric": "quantitative",
    "quant": "quantitative",
}
TRAIT_DATA_TYPE_ALIASES = {
    "": "auto",
    "auto": "auto",
    "ukb_data_field": "ukb_field",
    "ukb_field": "ukb_field",
    "ukb data field": "ukb_field",
    "icd10": "icd10",
    "icd_10": "icd10",
    "icd-10": "icd10",
    "phecode": "phecode",
    "phe_code": "phecode",
}
TRAIT_STOP_TOKENS = {
    "disease",
    "disorder",
    "syndrome",
    "trait",
    "traits",
    "phenotype",
    "phenotypes",
    "and",
    "or",
    "of",
    "in",
    "with",
    "without",
    "unspecified",
    "other",
    "type",
    "status",
    "measurement",
    "level",
    "levels",
    "ukb",
    "data",
    "field",
}
TRAIT_SYNONYM_SCOPE_RANK = {
    "LABEL": 4,
    "EXACT": 3,
    "NARROW": 2,
    "RELATED": 1,
    "BROAD": 0,
}
TRAIT_OVERSPECIFICITY_HINT_TOKENS = {
    "carcinoma",
    "adenocarcinoma",
    "sarcoma",
    "lymphoma",
    "leukemia",
    "leukaemia",
    "melanoma",
    "neuroblastoma",
    "glioma",
    "small",
    "large",
    "cell",
    "squamous",
    "invasive",
    "metastatic",
    "familial",
    "hereditary",
    "primary",
    "secondary",
    "acute",
    "chronic",
    "stage",
    "grade",
}
TRAIT_NON_DISEASE_ENTITY_PREFIXES = {
    "UBERON",
    "CHEBI",
    "PR",
    "GO",
    "CL",
    "HANCESTRO",
    "GSSO",
    "PATO",
    "OBA",
}
TRAIT_LOCATION_NONSPECIFIC_TOKENS = {
    "unspecified",
    "unknown",
    "other",
    "nos",
    "arthrosis",
    "osteoarthritis",
    "arthritis",
    "disorder",
    "disease",
    "syndrome",
    "condition",
    "diagnosis",
}
TRAIT_WEAK_FUZZY_TOKENS = {
    "incompetence",
    "insufficiency",
    "insufficient",
    "unspecified",
    "unknown",
    "other",
    "nos",
    "condition",
    "conditions",
    "disorder",
    "disease",
    "syndrome",
    "phenotype",
    "trait",
    "traits",
    "status",
    "case",
    "cases",
}
TRAIT_PIPE_RIGHT_WEAK_TOKENS = TRAIT_WEAK_FUZZY_TOKENS | {
    "inconclusive",
}
TRAIT_NEGATION_CUES = {
    "non",
    "without",
    "never",
    "no",
    "not",
    "none",
    "minus",
}
TRAIT_MULTI_CONNECTOR_RE = re.compile(r"(?:\b(?:and|or|vs|versus)\b|/|\\|,)")
UKB_DATA_FIELD_RE = re.compile(r"\bUKB(?:\s*data)?\s*field\s*(\d{1,7})\b", re.IGNORECASE)
UKB_CATEGORY_RE = re.compile(r"\bcategory\s*[:#\-\(\)\[\] ]*(\d{1,7})\b", re.IGNORECASE)
UKB_HASHED_TRAIT_RE = re.compile(r"^\s*(\d{1,7})#(.+)$")
STRICT_ICD10_CODE_RE = re.compile(rf"^{ICD10_CODE_HEAD_RE}(?:\.[0-9A-Z]{{1,4}})?$", re.IGNORECASE)
TRAIT_MULTI_HINT_TOKENS = {
    "combined",
    "pleiotropy",
    "multitrait",
    "multi",
    "comorbidity",
    "comorbidities",
}
TRAIT_CANONICAL_SIMILARITY_SUBS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"\bnon\s+high\s+density\s+lipoprotein\b"), "non hdl"),
    (re.compile(r"\bhigh\s+density\s+lipoprotein\b"), "hdl"),
    (re.compile(r"\bnon\s+low\s+density\s+lipoprotein\b"), "non ldl"),
    (re.compile(r"\blow\s+density\s+lipoprotein\b"), "ldl"),
    (re.compile(r"\bnon\s+very\s+low\s+density\s+lipoprotein\b"), "non vldl"),
    (re.compile(r"\bvery\s+low\s+density\s+lipoprotein\b"), "vldl"),
    (re.compile(r"\bintermediate\s+density\s+lipoprotein\b"), "idl"),
    (re.compile(r"\bhdl\s*c\b"), "hdl"),
    (re.compile(r"\bldl\s*c\b"), "ldl"),
    (re.compile(r"\bvldl\s*c\b"), "vldl"),
    (re.compile(r"\bidl\s*c\b"), "idl"),
    (re.compile(r"\bconcentration\s+of\s+idl\s+particles?\b"), "idl measurement"),
    (re.compile(r"\b(?:average|avg)\s+heart\s+rate\b"), "heart rate"),
    (re.compile(r"\b(?:average|avg)\s+pulse\s+rate\b"), "heart rate"),
)
TRAIT_STRICT_MODIFIER_TOKENS = {"free"}
TRAIT_QUERY_NOISE_SUBS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"\(\s*qc(?:['`’]?d)?\s*\)", re.IGNORECASE), " "),
    (re.compile(r"\(\s*gene\s*[- ]?based\s+burden\s*\)", re.IGNORECASE), " "),
    (re.compile(r"\bgene\s*[- ]?based\s+burden\b", re.IGNORECASE), " "),
    (re.compile(r"\(\s*unknown\s+primary\s*\)", re.IGNORECASE), " "),
    (re.compile(r"\+\s*/\s*-?\s*", re.IGNORECASE), " with or without "),
    (re.compile(r"\+\s*\|\s*-\s*", re.IGNORECASE), " with or without "),
    (re.compile(r"\bqc(?:['`’]?d)\b", re.IGNORECASE), " "),
    (re.compile(r"^\s*nmr\s+", re.IGNORECASE), " "),
)
TRAIT_NEGATION_FLAG_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\bnone\s+of\s+the\s+above\b", re.IGNORECASE),
    re.compile(r"\bnone\b", re.IGNORECASE),
    re.compile(r"\bno\b", re.IGNORECASE),
    re.compile(r"\bnot\b", re.IGNORECASE),
    re.compile(r"\bwithout\b", re.IGNORECASE),
    re.compile(r"\babsence\s+of\b", re.IGNORECASE),
)
TRAIT_NONINFORMATIVE_NEGATION_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\bnone\s+of\s+the\s+above\b", re.IGNORECASE),
    re.compile(r"^\s*none\b", re.IGNORECASE),
    re.compile(r"^\s*no\s+(?:known\s+|reported\s+)?(?:medical\s+)?(?:condition|disease|diagnosis)\b", re.IGNORECASE),
)
UKB_FIELD_SUPPLEMENT_PUBLICATION = "UKBFieldSupplement_v1"
UKB_FIELD_SUPPLEMENT_GENERIC_DISEASE_ID = "EFO_0000408"
UKB_FIELD_SUPPLEMENT_GENERIC_DISEASE_LABEL = "disease"
UKB_FIELD_SUPPLEMENT_GENERIC_MEASUREMENT_ID = "EFO_0001444"
UKB_FIELD_SUPPLEMENT_GENERIC_MEASUREMENT_LABEL = "measurement"
ICD10_SUPPLEMENT_PUBLICATION = "ICD10Supplement_v1"
TRAIT_CACHE_COLUMNS = [
    "UKBB data field",
    "ICD10",
    "PheCode",
    "Lookup text",
    "Curated reported trait",
    "main ontology URI(s)",
    "main ontology label(s)",
    "Publication",
    "User submitted trait",
    "qc_status",
    "qc_resolution_modes",
    "qc_label_match_modes",
    "qc_unresolved_ids",
    "final_update_applied",
    "final_update_reason",
]
DISEASE_HINT_TOKENS = {
    "disease",
    "disorder",
    "syndrome",
    "phenotype",
    "abnormality",
    "injury",
    "infection",
    "cancer",
    "carcinoma",
    "tumor",
    "tumour",
    "neoplasm",
    "arthritis",
    "myelitis",
    "pneumonia",
    "urticaria",
    "diabetes",
    "asthma",
    "eczema",
    "allergy",
    "allergic",
    "hypersensitivity",
    "anaphylaxis",
    "adenoma",
    "metastatic",
    "neoplasm",
    "depression",
    "anxiety",
    "phobia",
    "paralysis",
    "paresis",
    "palsy",
}
UKB_DISEASE_CATEGORY_PHRASES = (
    "medical condition",
    "medical conditions",
    "hospital",
    "cause of death",
    "cancer",
    "disease",
    "disorder",
    "diagnosis",
)
DISEASE_SUFFIX_HINTS = (
    "itis",
    "osis",
    "oma",
    "emia",
    "uria",
    "pathy",
    "algia",
    "phobia",
    "plegia",
    "paresis",
    "plasia",
)
MEASUREMENT_HINT_TOKENS = {
    "measurement",
    "amount",
    "concentration",
    "level",
    "levels",
    "ratio",
    "percentage",
    "percent",
    "count",
    "size",
    "diameter",
    "fraction",
    "abundance",
}
GENERIC_TRAIT_LABEL_KEYS = {
    "disease",
    "disorder",
    "condition",
    "abnormality",
    "phenotype",
    "syndrome",
    "medical procedure",
    "measurement",
}
TRAIT_QC_RISK_FIELDS = [
    "input_row_id",
    "input_query",
    "input_source_data_type",
    "input_type",
    "mapped_trait_id",
    "mapped_trait_label",
    "validation",
    "matched_via",
    "confidence",
    "specificity_loss_flag",
    "specificity_loss_notes",
    "risk_reasons",
    "evidence",
]
ICD10_PARENT_FALLBACK_METHODS = {
    "cache_parent_icd10",
    "cache_parent_icd10_supplement",
    "efo_obo_parent_icd10_xref",
}
ICD10_PARENT_FALLBACK_METHOD_KEYS = {
    "cache_parent_icd10",
    "cache_parent_icd10_supplement",
    "efo_obo_parent_icd10_xref",
}
TRAIT_PREFERRED_PREFIX_ORDER = {
    "EFO": 0,
    "MONDO": 1,
    "HP": 2,
    "OBA": 3,
    "NCIT": 4,
    "PATO": 5,
    "GO": 6,
    "MP": 7,
    "UBERON": 8,
    "HANCESTRO": 9,
    "GSSO": 10,
}


@dataclass
class Candidate:
    efo_id: str
    label: str
    score: float
    matched_via: str
    evidence: str
    is_validated: bool
    roundtrip_status: str = "not_applicable"
    roundtrip_accessions: tuple[str, ...] = ()
    roundtrip_symbols: tuple[str, ...] = ()


@dataclass(frozen=True)
class AccessionProfile:
    alias_terms: tuple[str, ...]
    expected_symbol_nums: frozenset[str]
    subject_terms: frozenset[str]
    symbol_terms: tuple[str, ...]


@dataclass(frozen=True)
class TraitCacheRecord:
    rownum: int
    reported_trait: str
    lookup_text: str
    icd10: str
    phecode: str
    mapped_ids: tuple[str, ...]
    mapped_labels: tuple[str, ...]
    publication: str
    is_ukb_field_supplement: bool
    is_icd10_supplement: bool


@dataclass(frozen=True)
class TraitOntologyTerm:
    term_id: str
    label: str
    synonyms: tuple[str, ...]
    label_key: str = ""
    exact_synonym_keys: frozenset[str] = frozenset()
    narrow_synonym_keys: frozenset[str] = frozenset()
    related_synonym_keys: frozenset[str] = frozenset()
    broad_synonym_keys: frozenset[str] = frozenset()


@dataclass(frozen=True)
class TraitObsoleteTerm:
    term_id: str
    label: str
    replaced_by: tuple[str, ...]
    consider: tuple[str, ...]


class ProgressReporter:
    def __init__(self, total: int, enabled: bool = True) -> None:
        self.total = max(0, total)
        self.enabled = enabled and self.total > 0
        self.done = 0
        self.start = time.time()
        self.last_tty_draw = 0.0
        self.is_tty = sys.stderr.isatty()
        self.next_non_tty_log = 0.05

        if self.enabled and not self.is_tty:
            print(f"[progress] 0/{self.total} (0.0%)", file=sys.stderr, flush=True)

    def _render_tty(self) -> None:
        now = time.time()
        if self.done < self.total and (now - self.last_tty_draw) < 0.08:
            return
        self.last_tty_draw = now

        ratio = self.done / self.total if self.total else 1.0
        width = 30
        filled = int(width * ratio)
        bar = "#" * filled + "-" * (width - filled)
        elapsed = max(0.0, now - self.start)
        rate = self.done / elapsed if elapsed > 0 else 0.0
        eta = (self.total - self.done) / rate if rate > 0 else 0.0
        msg = f"\r[{bar}] {self.done}/{self.total} {ratio*100:5.1f}%  elapsed {elapsed:5.1f}s  eta {eta:5.1f}s"
        print(msg, end="", file=sys.stderr, flush=True)
        if self.done >= self.total:
            print(file=sys.stderr, flush=True)

    def _render_non_tty(self) -> None:
        ratio = self.done / self.total if self.total else 1.0
        if self.done >= self.total or ratio >= self.next_non_tty_log:
            print(f"[progress] {self.done}/{self.total} ({ratio*100:.1f}%)", file=sys.stderr, flush=True)
            self.next_non_tty_log += 0.05

    def update(self, inc: int = 1) -> None:
        if not self.enabled:
            return
        self.done = min(self.total, self.done + max(0, inc))
        if self.is_tty:
            self._render_tty()
        else:
            self._render_non_tty()


@lru_cache(maxsize=200_000)
def normalize(text: str | None) -> str:
    if text is None:
        return ""
    stripped = text.strip()
    if not stripped:
        return ""
    if "\n" in stripped or "\t" in stripped or "\r" in stripped or "  " in stripped:
        return re.sub(r"\s+", " ", stripped)
    return stripped


@lru_cache(maxsize=200_000)
def normalize_id(value: str) -> str:
    return normalize(value).upper().replace("_", ":")


@lru_cache(maxsize=200_000)
def format_ontology_id_for_output(value: str) -> str:
    efo_id = normalize_id(value)
    if not efo_id:
        return ""
    return efo_id.replace(":", "_")


@lru_cache(maxsize=200_000)
def norm_key(value: str) -> str:
    return normalize(value).lower()


@lru_cache(maxsize=50_000)
def normalize_input_type(value: str) -> str:
    key = norm_key(value).replace("-", "_").replace(" ", "_")
    return INPUT_TYPE_ALIASES.get(key, "auto")


@lru_cache(maxsize=50_000)
def normalize_entity_type(value: str) -> str:
    key = norm_key(value).replace("-", "_").replace(" ", "_")
    return ENTITY_TYPE_ALIASES.get(key, "auto")


def effective_name_mode(default_name_mode: str, input_type: str) -> str:
    kind = normalize_input_type(input_type)
    if kind == "protein_name":
        # Protein-name inputs are often not uniquely resolvable via strict alias
        # lookup; enable fuzzy candidate retrieval but rely on validation gates.
        return "fuzzy"
    if kind == "metabolite_name":
        # Metabolite free-text names are commonly synonym-heavy and benefit from
        # lexical candidate retrieval before downstream validation gates.
        return "fuzzy"
    if kind == "gene_symbol":
        return "strict"
    if kind == "gene_id":
        return "strict"
    if kind == "accession":
        return "strict"
    if kind == "metabolite_id":
        return "strict"
    return default_name_mode


def tokenize(text: str) -> set[str]:
    return {tok for tok in re.split(r"[^a-z0-9]+", text.lower()) if tok}


@lru_cache(maxsize=200_000)
def canonical_accession(value: str) -> str:
    token = normalize(value).upper()
    if UNIPROT_ACCESSION_RE.match(token):
        return token
    if "-" in token:
        base = token.split("-", 1)[0]
        if UNIPROT_ACCESSION_RE.match(base):
            return base
    return ""


@lru_cache(maxsize=200_000)
def gene_id_alias_variants(value: str, include_unprefixed_numeric: bool = False) -> tuple[str, ...]:
    token = normalize(value).upper().replace(" ", "")
    if not token:
        return ()
    variants: set[str] = set()

    ensembl_match = ENSEMBL_GENE_ID_RE.match(token)
    if ensembl_match:
        variants.add(token.split(".", 1)[0])

    entrez_match = ENTREZ_GENE_ID_RE.match(token)
    if entrez_match:
        digits = entrez_match.group(1)
        variants.add(f"NCBIGENE:{digits}")
        variants.add(f"GENEID:{digits}")
        variants.add(f"ENTREZ:{digits}")
        if include_unprefixed_numeric:
            variants.add(digits)

    return tuple(sorted(variants))


@lru_cache(maxsize=200_000)
def canonical_metabolite_id(value: str) -> str:
    token = normalize(value).upper().replace(" ", "")
    if not token:
        return ""
    token = token.replace("_", ":")

    hmdb_match = HMDB_ID_RE.match(token)
    if hmdb_match:
        digits = hmdb_match.group(1)
        # HMDB accessions are typically zero-padded to 7+ digits.
        return f"HMDB{digits.zfill(7)}"

    chebi_match = CHEBI_ID_RE.match(token)
    if chebi_match:
        return f"CHEBI:{chebi_match.group(1)}"

    kegg_match = KEGG_ID_RE.match(token)
    if kegg_match:
        return kegg_match.group(1).upper()

    return ""


@lru_cache(maxsize=200_000)
def normalize_lipid_panel_phrase(value: str) -> str:
    phrase = normalize(value)
    if not phrase:
        return ""
    phrase = LIPID_PANEL_PREFIX_RE.sub("", phrase)
    phrase = normalize(re.sub(r"\s*\([^)]*\)", " ", phrase))
    phrase = normalize(phrase)
    for pattern, replacement in LIPID_PANEL_REPLACEMENTS:
        phrase = pattern.sub(replacement, phrase)
    phrase = normalize(re.sub(r"\s+", " ", phrase))
    return phrase


@lru_cache(maxsize=200_000)
def lipid_panel_query_variants(value: str) -> tuple[str, ...]:
    query = normalize(value)
    if not query:
        return ()
    variants: set[str] = {query}
    normalized = normalize_lipid_panel_phrase(query)
    if normalized:
        variants.add(normalized)

    for variant in list(variants):
        base = normalize(variant)
        if not base:
            continue
        singular_forms = {
            normalize(base.replace("cholesteryl esters", "cholesteryl ester")),
            normalize(base.replace("phospholipids", "phospholipid")),
            normalize(base.replace("triglycerides", "triglyceride")),
        }
        for form in singular_forms:
            if form:
                variants.add(form)
        if "/" in base:
            variants.add(normalize(base.replace("/", " / ")))
        for form in {base, *singular_forms}:
            if form and LIPID_MEASUREMENT_VARIANT_RE.match(norm_key(form)):
                variants.add(normalize(f"{form} measurement"))
                variants.add(normalize(f"{form} level"))

    return tuple(v for v in variants if v)


def query_has_protein_cue(value: str) -> bool:
    normalized = normalize(value)
    if PROTEIN_CUE_RE.search(normalized):
        return True
    folded = normalize_lipid_panel_phrase(normalized)
    return bool(folded and PROTEIN_CUE_RE.search(folded))


def query_has_metabolite_cue(value: str) -> bool:
    normalized = normalize(value)
    if METABOLITE_CUE_RE.search(normalized):
        return True
    folded = normalize_lipid_panel_phrase(normalized)
    return bool(folded and METABOLITE_CUE_RE.search(folded))


@lru_cache(maxsize=300_000)
def metabolite_descriptor_tag(value: str) -> str:
    text = norm_key(value)
    if not text:
        return ""
    if "free cholesterol" in text:
        return "free_cholesterol"
    if "cholesteryl ester" in text or "cholesterol ester" in text or "esterified cholesterol" in text:
        return "cholesteryl_ester"
    if re.search(r"\bcholesterol\b", text):
        return "cholesterol"
    if re.search(r"\bphospholipid(?:s)?\b", text):
        return "phospholipid"
    if re.search(r"\btriglyceride(?:s)?\b", text):
        return "triglyceride"
    return ""


@lru_cache(maxsize=300_000)
def has_lipid_chain_signature(value: str) -> bool:
    text = norm_key(value)
    if not text:
        return False
    return re.search(r"\b(?:c)?\d{1,2}:\d{1,2}\b", text) is not None


@lru_cache(maxsize=300_000)
def lipoprotein_class_tags(value: str) -> frozenset[str]:
    text = norm_key(value)
    if not text:
        return frozenset()
    tags: set[str] = set()
    for token in ("vldl", "ldl", "hdl", "idl"):
        if re.search(rf"\b{token}\b", text):
            tags.add(token)
    return frozenset(tags)


@lru_cache(maxsize=300_000)
def lipoprotein_size_tags(value: str) -> frozenset[str]:
    text = norm_key(value)
    if not text:
        return frozenset()
    tags: set[str] = set()
    if re.search(r"\bvery small\b", text):
        tags.add("very_small")
    if re.search(r"\bvery large\b", text) or re.search(r"\bextremely large\b", text):
        tags.add("very_large")
    # Add base sizes only when the "very" qualifiers are absent.
    if "very small" not in text and re.search(r"\bsmall\b", text):
        tags.add("small")
    if "very large" not in text and "extremely large" not in text and re.search(r"\blarge\b", text):
        tags.add("large")
    if re.search(r"\bmedium\b", text):
        tags.add("medium")
    return frozenset(tags)


@lru_cache(maxsize=300_000)
def apolipoprotein_subtype_tags(value: str) -> frozenset[str]:
    text = norm_key(value)
    if not text:
        return frozenset()
    tags: set[str] = set()
    if re.search(r"\b(?:apo[-\s]?a1|apoa1|apolipoprotein\s+a[-\s]?(?:1|i))\b", text):
        tags.add("a1")
    if re.search(r"\b(?:apo[-\s]?a2|apoa2|apolipoprotein\s+a[-\s]?(?:2|ii))\b", text):
        tags.add("a2")
    if re.search(
        r"\b(?:apo[-\s]?(?:b(?:100)?|8|β|Β|в|В|ß)|apob|apolipoprotein\s+b(?:[-\s]?100)?)\b",
        text,
    ):
        tags.add("b")
    return frozenset(tags)


@lru_cache(maxsize=300_000)
def protein_discriminator_tags(value: str) -> frozenset[str]:
    text = norm_key(value)
    if not text:
        return frozenset()
    tags: set[str] = set()
    for tok in subject_core_tokens(text):
        if any(ch.isdigit() for ch in tok):
            tags.add(tok)
            continue
        if tok in ROMAN_NUMERAL_TOKENS:
            tags.add(tok)
            continue
        if split_symbol_tail(tok) is not None:
            tags.add(tok)
    return frozenset(tags)


@lru_cache(maxsize=200_000)
def has_all_lipoprotein_scope(value: str) -> bool:
    text = norm_key(value)
    if not text:
        return False
    return re.search(r"\bin all (?:vldl|ldl|hdl|idl)\b", text) is not None


@lru_cache(maxsize=300_000)
def measurement_state_tags(value: str) -> frozenset[str]:
    text = norm_key(value)
    if not text:
        return frozenset()
    tags: set[str] = set()
    if re.search(r"\bfasting\b|\bduring fasting\b", text):
        tags.add("fasting")
    if re.search(r"\bpostprandial\b|\bafter meal\b|\bnon[- ]fasting\b", text):
        tags.add("postprandial")
    if re.search(r"\bgestational\b|\bgravid\b|\bpregnan", text):
        tags.add("gestational")
    return frozenset(tags)


def condition_mismatch_penalty(query: str, label: str, synonyms: list[str]) -> float:
    query_tags = measurement_state_tags(query)
    candidate_tags = measurement_state_tags(" ".join([label] + synonyms))
    if not candidate_tags:
        return 0.0
    if query_tags:
        if query_tags & candidate_tags:
            return 0.0
        # Query requested a specific state, candidate does not match it.
        return 0.25
    # Query is unqualified; avoid selecting state-specific terms by default.
    return 0.12


@lru_cache(maxsize=400_000)
def normalize_metabolite_alias_text(value: str, ascii_fallback: str = "") -> str:
    raw = normalize(value)
    if not raw:
        raw = normalize(ascii_fallback)
    if not raw:
        return ""
    # ChEBI names.tsv frequently contains HTML formatting tags.
    clean = html.unescape(raw)
    clean = HTML_TAG_RE.sub("", clean)
    clean = normalize(clean)
    if not clean:
        clean = normalize(ascii_fallback)
    if not clean:
        return ""
    # Normalize common minus glyph variants and excessive punctuation spacing.
    clean = clean.replace("\u2212", "-")
    clean = clean.replace("<", " ").replace(">", " ")
    clean = re.sub(r"\s*([,/;:+])\s*", r"\1 ", clean)
    clean = re.sub(r"\s+", " ", clean).strip(" _")
    clean = normalize(clean)
    canonical = canonical_metabolite_id(clean)
    if canonical:
        return canonical
    return clean


@lru_cache(maxsize=400_000)
def metabolite_alias_quality_rank(alias: str) -> tuple[int, int, int, int, str]:
    text = normalize_metabolite_alias_text(alias)
    if not text:
        return (9, 9, 9, 999, "")
    lower = text.lower()
    toks = tokenize(lower)
    letters = sum(ch.isalpha() for ch in text)
    digits = sum(ch.isdigit() for ch in text)
    punct = sum((not ch.isalnum()) and (not ch.isspace()) for ch in text)
    length_penalty = 0
    if len(text) < 3:
        length_penalty = 3
    elif len(text) > 120:
        length_penalty = 3
    elif len(text) > 80:
        length_penalty = 2
    elif len(text) > 50:
        length_penalty = 1

    noise_penalty = 0
    if METABOLITE_NOISE_RE.match(text):
        noise_penalty += 3
    if letters == 0 and not canonical_metabolite_id(text):
        noise_penalty += 3
    if digits > letters and letters > 0:
        noise_penalty += 1
    if punct > max(2, len(text) // 3):
        noise_penalty += 1
    if len(toks) == 1 and len(next(iter(toks), "")) <= 2 and not canonical_metabolite_id(text):
        noise_penalty += 2

    # Prefer commonly used names over synthetic markup-like labels.
    commonness_penalty = 0
    if re.search(r"\b(?:iupac|brand name)\b", lower):
        commonness_penalty += 1
    return (noise_penalty, length_penalty, commonness_penalty, len(text), lower)


@lru_cache(maxsize=400_000)
def keep_metabolite_alias(alias: str) -> bool:
    text = normalize_metabolite_alias_text(alias)
    if not text:
        return False
    if canonical_metabolite_id(text):
        return True
    if len(text) < 2:
        return False
    if METABOLITE_NOISE_RE.match(text):
        return False
    letters = sum(ch.isalpha() for ch in text)
    if letters == 0 and sum(ch.isdigit() for ch in text) > 0:
        return False
    return True


def prune_metabolite_aliases(label: str, aliases: list[str], ids: list[str], max_aliases: int = 64) -> list[str]:
    label_clean = normalize_metabolite_alias_text(label)
    cleaned: set[str] = set()
    for alias in aliases:
        alias_clean = normalize_metabolite_alias_text(alias)
        if alias_clean and keep_metabolite_alias(alias_clean):
            cleaned.add(alias_clean)
    if label_clean:
        cleaned.add(label_clean)
    for met_id in ids:
        canonical = canonical_metabolite_id(met_id)
        if canonical:
            cleaned.add(canonical)

    ids_only = sorted([a for a in cleaned if canonical_metabolite_id(a)])
    text_aliases = [a for a in cleaned if not canonical_metabolite_id(a)]
    text_aliases.sort(key=metabolite_alias_quality_rank)

    if len(text_aliases) > max_aliases:
        text_aliases = text_aliases[:max_aliases]
    return sorted(set(ids_only + text_aliases))


def extract_metabolite_ids_from_text(value: str) -> list[str]:
    text = normalize(value)
    if not text:
        return []
    ids: set[str] = set()
    # Try direct canonical parse first.
    direct = canonical_metabolite_id(text)
    if direct:
        ids.add(direct)

    for raw in re.split(r"[|;, \t]+", text):
        token = normalize(raw)
        if not token:
            continue
        canonical = canonical_metabolite_id(token)
        if canonical:
            ids.add(canonical)

    for match in METABOLITE_CHEBI_INLINE_RE.findall(text):
        ids.add(f"CHEBI:{match}")
    for match in METABOLITE_HMDB_INLINE_RE.findall(text):
        ids.add(f"HMDB{match.zfill(7)}")
    for match in METABOLITE_KEGG_INLINE_RE.findall(text):
        ids.add(match.upper())
    return sorted(ids)


def choose_chebi_primary_label(name_rows: list[tuple[str, str, str]]) -> str:
    if not name_rows:
        return ""

    type_rank = {
        "INN": 0,
        "UNIPROT NAME": 1,
        "SYNONYM": 2,
        "IUPAC NAME": 3,
        "BRAND NAME": 4,
        "": 5,
    }
    status_rank = {
        "1": 0,
        "3": 1,
        "9": 2,
        "": 3,
    }

    best = ""
    best_key: tuple[int, int, int, int, str] | None = None
    for name, row_type, status_id in name_rows:
        clean = normalize_metabolite_alias_text(name)
        if not clean or not keep_metabolite_alias(clean):
            continue
        type_key = type_rank.get(normalize(row_type).upper(), 6)
        status_key = status_rank.get(normalize(status_id), 4)
        quality = metabolite_alias_quality_rank(clean)
        key = (type_key, status_key, *quality)
        if best_key is None or key < best_key:
            best_key = key
            best = clean
    return best


def metabolite_id_bucket(met_id: str) -> str:
    mid = canonical_metabolite_id(met_id)
    if mid.startswith("HMDB"):
        return "hmdb"
    if mid.startswith("CHEBI:"):
        return "chebi"
    if re.match(r"^[CD]\d{5}$", mid):
        return "kegg"
    return "other"


@lru_cache(maxsize=100_000)
def query_variants(query: str) -> tuple[str, ...]:
    q = normalize(query)
    variants: set[str] = {q}
    parts = [normalize(p) for p in MULTI_ID_SEP_RE.split(q) if normalize(p)]
    variants.update(parts)
    for part in list(variants):
        variants.update(lipid_panel_query_variants(part))
    for part in list(variants):
        canonical = canonical_accession(part)
        if canonical:
            variants.add(canonical)
        metabolite_id = canonical_metabolite_id(part)
        if metabolite_id:
            variants.add(metabolite_id)
    return tuple(v for v in variants if v)


def is_measurement_like(label: str, synonyms: list[str]) -> bool:
    bag = tokenize(label)
    for syn in synonyms:
        bag |= tokenize(syn)
    return bool(bag & MEASUREMENT_KEYWORDS)


def similarity_score_prepared(query_lower: str, query_toks: set[str], label: str, synonyms: list[str]) -> float:
    label_toks = tokenize(label)
    overlap = len(query_toks & label_toks) / max(len(query_toks), 1)

    seq = difflib.SequenceMatcher(None, query_lower, label.lower()).ratio()
    syn_best = 0.0
    for syn in synonyms:
        syn_best = max(syn_best, difflib.SequenceMatcher(None, query_lower, syn.lower()).ratio())

    keyword_bonus = 0.2 if label_toks & MEASUREMENT_KEYWORDS else 0.0
    score = 0.45 * seq + 0.35 * syn_best + 0.2 * overlap + keyword_bonus
    return min(1.0, max(0.0, score))


def similarity_score(query: str, label: str, synonyms: list[str]) -> float:
    query_lower = query.lower()
    query_toks = tokenize(query)
    return similarity_score_prepared(query_lower, query_toks, label, synonyms)


def load_query_inputs(path: Path) -> list[tuple[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".txt", ".list"}:
        return [(normalize(line), "auto") for line in path.read_text(encoding="utf-8").splitlines() if normalize(line)]

    delimiter = "\t" if suffix in {".tsv", ".tab"} else ","
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter=delimiter))

    if not rows:
        return []

    # Support single-column headerless TSV/CSV lists where every line is a query.
    # If the first value is a known header name, skip it; otherwise include all rows.
    if rows and len(rows[0]) == 1 and all(len(r) <= 1 for r in rows):
        first = normalize(rows[0][0] if rows[0] else "")
        start_idx = 1 if first.lower() in SEARCHABLE_COLUMNS else 0
        values: list[tuple[str, str]] = []
        for r in rows[start_idx:]:
            raw = r[0] if r else ""
            if raw and normalize(raw):
                values.append((normalize(raw), "auto"))
        return values

    # Multi-column inputs are treated as headered tabular files.
    fieldnames = rows[0]
    field_lookup = {normalize(f).lower(): f for f in fieldnames}
    source_field = None
    for name in SEARCHABLE_COLUMNS:
        if name in field_lookup:
            source_field = field_lookup[name]
            break
    if source_field is None:
        source_field = fieldnames[0]
    type_field = None
    for name in INPUT_TYPE_COLUMNS:
        if name in field_lookup:
            type_field = field_lookup[name]
            break

    values: list[tuple[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            raw = row.get(source_field, "")
            if raw and normalize(raw):
                raw_type = row.get(type_field, "") if type_field else ""
                values.append((normalize(raw), normalize_input_type(raw_type)))
    return values


def load_queries(path: Path) -> list[str]:
    return [query for query, _input_type in load_query_inputs(path)]


@lru_cache(maxsize=50_000)
def normalize_trait_input_type(value: str) -> str:
    key = norm_key(value).replace("-", "_").replace(" ", "_")
    return TRAIT_INPUT_TYPE_ALIASES.get(key, "auto")


@lru_cache(maxsize=50_000)
def normalize_trait_scale(value: str) -> str:
    key = norm_key(value).replace("-", "_").replace(" ", "_")
    return TRAIT_SCALE_ALIASES.get(key, "auto")


@lru_cache(maxsize=50_000)
def normalize_trait_data_type(value: str) -> str:
    key = norm_key(value).replace("-", "_").replace(" ", "_")
    return TRAIT_DATA_TYPE_ALIASES.get(key, "auto")


@lru_cache(maxsize=200_000)
def normalize_icd10_code(value: str) -> str:
    token = normalize(value).upper()
    if not token:
        return ""
    token = token.replace("ICD-10", "ICD10").replace("ICD 10", "ICD10")

    def compact_to_strict(raw: str) -> str:
        compact = normalize(raw).upper().replace("ICD10", "").replace(" ", "")
        compact = compact.strip("()[]{}")
        if STRICT_ICD10_CODE_RE.fullmatch(compact):
            return compact
        # Common compact form in union-style strings, for example A071 -> A07.1
        compact_match = re.fullmatch(rf"({ICD10_CODE_HEAD_RE})([0-9A-Z]{{1,4}})", compact)
        if compact_match:
            return f"{compact_match.group(1)}.{compact_match.group(2)}"
        return ""

    # Parse tokenized segments first so mixed strings like Union#A071#A07.1
    # resolve to the most specific strict code.
    segmented_candidates: list[str] = []
    for raw_part in re.split(r"[\s,;|#/\\]+", token):
        code = compact_to_strict(raw_part)
        if code:
            segmented_candidates.append(code)
    if segmented_candidates:
        segmented_candidates.sort(key=lambda item: ("." not in item, -len(item), item))
        return segmented_candidates[0]

    strict_match = re.search(rf"\b({ICD10_CODE_HEAD_RE}(?:\.[0-9A-Z]{{1,4}})?)\b", token)
    if strict_match:
        return strict_match.group(1)

    compact_match = re.search(rf"\b({ICD10_CODE_HEAD_RE}[0-9A-Z]{{1,4}})\b", token)
    if compact_match:
        compact_code = compact_to_strict(compact_match.group(1))
        if compact_code:
            return compact_code

    if ":" in token:
        token = token.split(":", 1)[0]
    token = token.replace("ICD10", "").strip().replace(" ", "")
    return token


def is_strict_icd10_code(value: str) -> bool:
    token = normalize(value).upper()
    return bool(token and STRICT_ICD10_CODE_RE.fullmatch(token))


@lru_cache(maxsize=500_000)
def icd10_parent_fallback_codes(value: str) -> tuple[str, ...]:
    token = normalize_icd10_code(value)
    if not is_strict_icd10_code(token):
        return ()

    variants: list[str] = []
    if "." in token:
        head, tail = token.split(".", 1)
        for width in range(len(tail) - 1, 0, -1):
            parent = f"{head}.{tail[:width]}"
            if parent != token and is_strict_icd10_code(parent) and parent not in variants:
                variants.append(parent)
        if head != token and is_strict_icd10_code(head) and head not in variants:
            variants.append(head)
    return tuple(variants)


@lru_cache(maxsize=500_000)
def extract_icd10_code_from_xref(value: str) -> str:
    raw = normalize(value)
    if not raw:
        return ""
    decoded = urllib.parse.unquote(raw)
    for pattern in (ICD10_XREF_RE, ICD10_XREF_URL_RE):
        match = pattern.search(decoded)
        if match:
            code = normalize_icd10_code(match.group(1))
            if code:
                return code
    return ""


@lru_cache(maxsize=500_000)
def extract_mondo_id(value: str) -> str:
    raw = normalize(value)
    if not raw:
        return ""
    direct = canonicalize_trait_ontology_id(raw)
    if trait_ontology_prefix(direct) == "MONDO":
        return direct
    decoded = urllib.parse.unquote(raw)
    match = MONDO_CURIE_RE.search(decoded)
    if not match:
        return ""
    token = canonicalize_trait_ontology_id(match.group(0))
    return token if trait_ontology_prefix(token) == "MONDO" else ""


@lru_cache(maxsize=200_000)
def normalize_phecode(value: str) -> str:
    token = normalize(value)
    if not token:
        return ""
    match = re.search(r"\b([0-9]{1,4}(?:\.[0-9A-Z]{1,4})?)\b", token.upper())
    if match:
        return match.group(1)
    return token.replace(" ", "")


@lru_cache(maxsize=200_000)
def looks_like_phecode_token(value: str) -> bool:
    raw = normalize(value).upper()
    if not raw:
        return False
    raw = raw.replace("PHECODE", "").replace("PHE", "").strip(" :#-")
    return bool(re.fullmatch(r"[0-9]{1,4}(?:\.[0-9A-Z]{1,4})?", raw))


@lru_cache(maxsize=200_000)
def canonicalize_trait_ontology_id(value: str) -> str:
    raw = normalize(value)
    if not raw:
        return ""
    token = raw.split()[0]
    if token.lower().startswith("http"):
        token = token.rsplit("/", 1)[-1]
    if ":" in token:
        parts = token.split(":")
        last = parts[-1]
        if re.match(r"^[A-Za-z]+_[0-9A-Za-z]+$", last):
            token = last
        elif (
            len(parts) >= 2
            and re.match(r"^[A-Za-z]+$", parts[-2])
            and re.match(r"^[0-9A-Za-z_]+$", last)
        ):
            token = f"{parts[-2].upper()}_{last}"
        else:
            token = last
    token = token.replace(":", "_")
    if "_" in token:
        prefix, local = token.split("_", 1)
        return f"{prefix.upper()}_{local}"
    return token.upper()


def trait_ontology_prefix(ontology_id: str) -> str:
    oid = canonicalize_trait_ontology_id(ontology_id)
    if "_" not in oid:
        return oid
    return oid.split("_", 1)[0]


def split_multi_ids(value: str) -> list[str]:
    return [normalize(part) for part in re.split(r"[|,;]+", str(value)) if normalize(part)]


def split_multi_labels(value: str, expected_n: int) -> list[str]:
    raw = normalize(value)
    if not raw:
        return []
    if "|" in raw:
        return [normalize(part) for part in raw.split("|") if normalize(part)]
    if expected_n > 1 and "," in raw:
        return [normalize(part) for part in raw.split(",") if normalize(part)]
    return [raw]


def informative_trait_tokens(text: str) -> set[str]:
    toks = {tok for tok in tokenize(text) if len(tok) >= 3}
    filtered = {tok for tok in toks if tok not in TRAIT_STOP_TOKENS}
    return filtered or toks


def trait_nonweak_tokens(tokens: set[str]) -> set[str]:
    return {
        tok
        for tok in tokens
        if tok
        and tok not in TRAIT_STOP_TOKENS
        and tok not in TRAIT_WEAK_FUZZY_TOKENS
    }


def trait_anchor_tokens(tokens: set[str]) -> set[str]:
    anchors: set[str] = set()
    for tok in trait_nonweak_tokens(tokens):
        if len(tok) < 4:
            continue
        if tok in DISEASE_HINT_TOKENS:
            continue
        if tok in TRAIT_LOCATION_NONSPECIFIC_TOKENS:
            continue
        anchors.add(tok)
    return anchors


@lru_cache(maxsize=300_000)
def trim_trailing_weak_trait_tokens(text: str) -> str:
    key = norm_key(text)
    if not key:
        return ""
    toks = [tok for tok in re.findall(r"[a-z0-9]+", key) if tok]
    while toks and toks[-1] in TRAIT_PIPE_RIGHT_WEAK_TOKENS:
        toks.pop()
    return " ".join(toks)


@lru_cache(maxsize=300_000)
def strip_trait_query_noise(text: str) -> str:
    clean = normalize(text)
    if not clean:
        return ""
    for pattern, replacement in TRAIT_QUERY_NOISE_SUBS:
        clean = pattern.sub(replacement, clean)
    clean = re.sub(r"\s+", " ", clean).strip()
    return normalize(clean)


@lru_cache(maxsize=300_000)
def trait_query_exact_norm_variants(text: str) -> tuple[str, ...]:
    variants: list[str] = []
    for candidate in (norm_key(text), norm_key(strip_trait_query_noise(text))):
        if candidate and candidate not in variants:
            variants.append(candidate)
        paren_stripped = normalize(re.sub(r"\([^)]*\)", " ", candidate))
        paren_norm = norm_key(paren_stripped)
        if paren_norm and paren_norm not in variants:
            variants.append(paren_norm)
        trimmed = trim_trailing_weak_trait_tokens(candidate)
        if trimmed and trimmed not in variants:
            variants.append(trimmed)
        if "|" in candidate:
            for part in candidate.split("|"):
                part_norm = norm_key(part)
                if part_norm and part_norm not in variants:
                    variants.append(part_norm)
                part_trimmed = trim_trailing_weak_trait_tokens(part_norm)
                if part_trimmed and part_trimmed not in variants:
                    variants.append(part_trimmed)
    return tuple(variants)


@lru_cache(maxsize=300_000)
def normalize_trait_text_for_similarity(text: str) -> str:
    key = norm_key(strip_trait_query_noise(text))
    if not key:
        return ""
    for pattern, replacement in TRAIT_CANONICAL_SIMILARITY_SUBS:
        key = pattern.sub(replacement, key)
    return re.sub(r"\s+", " ", key).strip()


@lru_cache(maxsize=300_000)
def trait_has_negation(text: str) -> bool:
    raw = normalize(text)
    if not raw:
        return False
    clean = normalize_trait_text_for_similarity(raw) or norm_key(raw)
    if not clean:
        return False
    # "with or without X" is an inclusive clinical expression, not a strict
    # negation cue that should block validation.
    clean = re.sub(r"\bwith\s+or\s+without\b", " ", clean, flags=re.IGNORECASE)
    return any(pattern.search(clean) for pattern in TRAIT_NEGATION_FLAG_PATTERNS)


@lru_cache(maxsize=300_000)
def trait_is_noninformative_negation(text: str) -> bool:
    raw = normalize(text)
    if not raw:
        return False
    return any(pattern.search(raw) for pattern in TRAIT_NONINFORMATIVE_NEGATION_PATTERNS)


@lru_cache(maxsize=200_000)
def extract_ukb_field_ids(text: str) -> tuple[str, ...]:
    key = normalize(text)
    if not key:
        return ()
    field_ids: list[str] = []
    for match in UKB_DATA_FIELD_RE.findall(key):
        field_id = normalize(match)
        if field_id and field_id not in field_ids:
            field_ids.append(field_id)
    return tuple(field_ids)


@lru_cache(maxsize=200_000)
def parse_ukb_hashed_trait_text(text: str) -> tuple[str, str, str]:
    """Parse common UKB coded trait strings like `20002#1425#cerebral aneurysm`.

    Returns `(field_id, code, label)` with empty strings when a component is not
    present or the pattern does not match.
    """
    raw = normalize(text)
    if not raw:
        return ("", "", "")
    match = UKB_HASHED_TRAIT_RE.match(raw)
    if not match:
        return ("", "", "")
    field_id = normalize(match.group(1))
    remainder = normalize(match.group(2))
    if not field_id or not remainder:
        return ("", "", "")

    parts = [normalize(part) for part in remainder.split("#")]

    def clean_label(raw_label: str, parsed_field_id: str) -> str:
        label = normalize(raw_label)
        if not label:
            return ""
        label = re.sub(r"\+\s*/\s*-?\s*", " with or without ", label, flags=re.IGNORECASE)
        label = re.sub(r"\+\s*\|\s*-\s*", " with or without ", label, flags=re.IGNORECASE)
        label = re.sub(
            rf"\s*\(UKB(?:\s*data)?\s*field\s*(?:{re.escape(parsed_field_id)}|[0-9]{{1,7}})\)\s*$",
            "",
            label,
            flags=re.IGNORECASE,
        )
        label = re.sub(rf"^\s*{ICD10_CODE_HEAD_RE}(?:\.[0-9A-Z]{{1,4}})?\s+", "", label, flags=re.IGNORECASE)
        label = re.sub(rf"^\s*{ICD10_CODE_HEAD_RE}[0-9A-Z]{{1,4}}\s+", "", label, flags=re.IGNORECASE)
        label = re.sub(r"\bdvt\b", "deep venous thrombosis", label, flags=re.IGNORECASE)
        if "|" in label:
            pipe_parts = [normalize(part) for part in label.split("|") if normalize(part)]
            if pipe_parts:
                # UKB coded disease labels often use `|` as alias-ish separators
                # (for example "aortic regurgitation | incompetence"). Keep the
                # clinical anchor phrase and only drop a weak trailing segment.
                right_tokens = set(tokenize(norm_key(pipe_parts[-1])))
                right_is_weak = bool(right_tokens) and right_tokens.issubset(TRAIT_PIPE_RIGHT_WEAK_TOKENS)
                if right_is_weak and len(pipe_parts) > 1:
                    kept = [part for part in pipe_parts[:-1] if part]
                    label = " ".join(kept) if kept else pipe_parts[0]
                else:
                    label = " ".join(pipe_parts)
            else:
                label = label.replace("|", " ")
        return normalize(label)

    if len(parts) == 1:
        return (field_id, "", clean_label(parts[0], field_id))

    code = normalize(parts[0])
    label = clean_label("#".join(parts[1:]), field_id)
    return (field_id, code, label)


@lru_cache(maxsize=300_000)
def parse_icd10_labeled_trait_text(text: str) -> tuple[str, str]:
    """Parse text snippets that embed ICD10 code + label.

    Examples:
    - `Union#M2571#M25.71 Osteophyte (Shoulder region)`
    - `40001#I269#I26.9 Pulmonary embolism without mention of acute cor pulmonale`
    """
    raw = normalize(text)
    if not raw:
        return ("", "")

    candidates: list[str] = [raw]
    if "#" in raw:
        last_segment = normalize(raw.rsplit("#", 1)[-1])
        if last_segment:
            candidates.insert(0, last_segment)

    for segment in candidates:
        clean = normalize(segment)
        if not clean:
            continue
        strict_match = re.match(rf"^\s*({ICD10_CODE_HEAD_RE}(?:\.[0-9A-Z]{{1,4}})?)\s+(.+)$", clean, re.IGNORECASE)
        if strict_match:
            code = normalize_icd10_code(strict_match.group(1))
            label = normalize(strict_match.group(2))
            if is_strict_icd10_code(code) and label:
                return code, label
            continue
        compact_match = re.match(rf"^\s*({ICD10_CODE_HEAD_RE}[0-9A-Z]{{1,4}})\s+(.+)$", clean, re.IGNORECASE)
        if compact_match:
            code = normalize_icd10_code(compact_match.group(1))
            label = normalize(compact_match.group(2))
            if is_strict_icd10_code(code) and label:
                return code, label

    return ("", "")


def ukb_title_query_variants(field_id: str, title: str) -> tuple[str, ...]:
    fid = normalize(field_id)
    clean_title = normalize(title)
    if not fid or not clean_title:
        return ()
    variants = [
        clean_title,
        f"{clean_title} (UKB data field {fid})",
        f"{clean_title} (UKB data field {fid}) (Gene-based burden)",
    ]
    deduped: list[str] = []
    seen: set[str] = set()
    for variant in variants:
        key = normalize(variant)
        if key and key not in seen:
            seen.add(key)
            deduped.append(key)
    return tuple(deduped)


@lru_cache(maxsize=200_000)
def extract_ukb_category_ids(text: str) -> tuple[str, ...]:
    key = normalize(text)
    if not key:
        return ()
    category_ids: list[str] = []
    for match in UKB_CATEGORY_RE.findall(key):
        category_id = normalize(match)
        if category_id and category_id not in category_ids:
            category_ids.append(category_id)
    return tuple(category_ids)


def ukb_category_query_variants(category_id: str, title: str, path_titles: tuple[str, ...]) -> tuple[str, ...]:
    cid = normalize(category_id)
    clean_title = normalize(title)
    clean_path = tuple(normalize(item) for item in path_titles if normalize(item))

    variants: list[str] = []
    if clean_title:
        variants.append(clean_title)
    if clean_path:
        variants.append(" ".join(clean_path))
        variants.append(" > ".join(clean_path))
    if cid and clean_title:
        variants.append(f"Category {cid} {clean_title}")
        variants.append(f"Category {cid}: {clean_title}")

    deduped: list[str] = []
    seen: set[str] = set()
    for variant in variants:
        key = normalize(variant)
        if key and key not in seen:
            seen.add(key)
            deduped.append(key)
    return tuple(deduped)


def ukb_category_path_titles(
    category_id: str,
    *,
    category_titles: dict[str, str],
    category_parents: dict[str, tuple[str, ...]],
    max_depth: int = 12,
) -> tuple[str, ...]:
    cid = normalize(category_id)
    if not cid:
        return ()
    seen: set[str] = set()
    chain: list[str] = []
    current = cid
    depth = 0
    while current and depth < max_depth and current not in seen:
        seen.add(current)
        chain.append(current)
        parents = category_parents.get(current, ())
        current = parents[0] if parents else ""
        depth += 1
    chain.reverse()
    return tuple(normalize(category_titles.get(cat_id, "")) for cat_id in chain if normalize(category_titles.get(cat_id, "")))


def load_ukb_field_category_index(path: Path | None) -> tuple[dict[str, str], dict[str, tuple[str, ...]]]:
    if path is None or not path.exists():
        return {}, {}

    field_to_category: dict[str, str] = {}
    category_to_fields: dict[str, list[str]] = {}
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                field_id = normalize(row.get("field_id", ""))
                category_id = normalize(row.get("main_category", ""))
                if not field_id or not category_id:
                    continue
                field_to_category[field_id] = category_id
                if field_id not in category_to_fields.setdefault(category_id, []):
                    category_to_fields[category_id].append(field_id)
    except OSError:
        return {}, {}

    frozen_category_to_fields = {cat: tuple(field_ids) for cat, field_ids in category_to_fields.items()}
    return field_to_category, frozen_category_to_fields


def load_ukb_category_title_index(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}

    titles: dict[str, str] = {}
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                category_id = normalize(row.get("category_id", ""))
                title = normalize(row.get("title", ""))
                if not category_id or not title:
                    continue
                titles[category_id] = title
    except OSError:
        return {}
    return titles


def load_ukb_category_parent_index(path: Path | None) -> dict[str, tuple[str, ...]]:
    if path is None or not path.exists():
        return {}

    parents: dict[str, list[str]] = {}
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                parent_id = normalize(row.get("parent_id", ""))
                child_id = normalize(row.get("child_id", ""))
                if not parent_id or not child_id:
                    continue
                if parent_id not in parents.setdefault(child_id, []):
                    parents[child_id].append(parent_id)
    except OSError:
        return {}

    return {child_id: tuple(parent_ids) for child_id, parent_ids in parents.items()}


def normalize_negation_target(token: str) -> str:
    tok = norm_key(token)
    if len(tok) <= 3:
        return tok
    if tok.endswith("ies"):
        return tok[:-3] + "y"
    if tok.endswith("s") and not tok.endswith("ss"):
        return tok[:-1]
    return tok


def extract_negated_trait_targets(text: str) -> set[str]:
    key = normalize_trait_text_for_similarity(text)
    if not key:
        return set()
    toks = [tok for tok in re.split(r"[^a-z0-9]+", key) if tok]
    targets: set[str] = set()
    for idx, tok in enumerate(toks[:-1]):
        if tok not in TRAIT_NEGATION_CUES:
            continue
        target = normalize_negation_target(toks[idx + 1])
        if not target or target in TRAIT_STOP_TOKENS:
            continue
        targets.add(target)
    return targets


def trait_negation_mismatch(query_text: str, candidate_text: str) -> bool:
    qnorm = normalize_trait_text_for_similarity(query_text)
    cnorm = normalize_trait_text_for_similarity(candidate_text)
    if not qnorm or not cnorm:
        return False
    q_neg = extract_negated_trait_targets(qnorm)
    c_neg = extract_negated_trait_targets(cnorm)
    if not q_neg and not c_neg:
        return False
    if q_neg == c_neg:
        return False
    q_toks = {normalize_negation_target(tok) for tok in tokenize(qnorm)}
    c_toks = {normalize_negation_target(tok) for tok in tokenize(cnorm)}
    return bool((q_neg & c_toks) or (c_neg & q_toks))


def trait_looks_multi_concept(text: str) -> bool:
    key = normalize_trait_text_for_similarity(text)
    if not key:
        return False
    if TRAIT_MULTI_CONNECTOR_RE.search(key):
        return True
    toks = set(tokenize(key))
    return bool(toks & TRAIT_MULTI_HINT_TOKENS)


def skip_multi_mapped_cache_fuzzy_candidate(
    *,
    query_text: str,
    query_tokens: set[str],
    record: TraitCacheRecord,
    mapped_id_count: int | None = None,
) -> bool:
    # For single-trait queries, avoid fuzzy promotion of composite cache rows
    # (e.g. "Crohn's disease or Leprosy") that can inject unrelated terms.
    count = mapped_id_count if mapped_id_count is not None else len(record.mapped_ids)
    if count <= 1:
        return False
    if not query_tokens:
        return False
    if trait_looks_multi_concept(query_text):
        return False

    candidate_text = record.lookup_text or record.reported_trait
    if not trait_looks_multi_concept(candidate_text):
        return False

    candidate_norm = normalize_trait_text_for_similarity(candidate_text)
    candidate_tokens = informative_trait_tokens(candidate_norm)
    if not candidate_tokens:
        return False
    return query_tokens.issubset(candidate_tokens)


def has_disease_hint(text: str) -> bool:
    toks = tokenize(text)
    if any(tok in DISEASE_HINT_TOKENS for tok in toks):
        return True
    return any(tok.endswith(DISEASE_SUFFIX_HINTS) for tok in toks)


def has_ukb_disease_context(text: str) -> bool:
    key = norm_key(text)
    if not key:
        return False
    if has_disease_hint(key):
        return True
    return any(phrase in key for phrase in UKB_DISEASE_CATEGORY_PHRASES)


def is_generic_trait_label(label: str) -> bool:
    key = norm_key(label)
    if not key:
        return True
    if key in GENERIC_TRAIT_LABEL_KEYS:
        return True
    toks = [tok for tok in tokenize(key) if tok not in TRAIT_STOP_TOKENS]
    if len(toks) <= 1 and any(tok in GENERIC_TRAIT_LABEL_KEYS for tok in toks):
        return True
    return False


def trait_query_label_overlap_ratio(query_text: str, label_text: str) -> float:
    query_norm = normalize_trait_text_for_similarity(query_text)
    label_norm = normalize_trait_text_for_similarity(label_text)
    query_tokens = informative_trait_tokens(query_norm)
    label_tokens = informative_trait_tokens(label_norm)
    if not query_tokens or not label_tokens:
        return 0.0
    return len(query_tokens & label_tokens) / max(1, len(query_tokens))


def trait_nonweak_overlap(query_text: str, label_text: str) -> set[str]:
    query_norm = normalize_trait_text_for_similarity(query_text)
    label_norm = normalize_trait_text_for_similarity(label_text)
    query_tokens = informative_trait_tokens(query_norm)
    label_tokens = informative_trait_tokens(label_norm)
    if not query_tokens or not label_tokens:
        return set()
    overlap = query_tokens & label_tokens
    return trait_nonweak_tokens(overlap)


def trait_has_anchor_overlap(query_text: str, label_text: str) -> bool:
    query_norm = normalize_trait_text_for_similarity(query_text)
    label_norm = normalize_trait_text_for_similarity(label_text)
    query_tokens = informative_trait_tokens(query_norm)
    label_tokens = informative_trait_tokens(label_norm)
    if not query_tokens or not label_tokens:
        return False
    query_anchors = trait_anchor_tokens(query_tokens)
    label_anchors = trait_anchor_tokens(label_tokens)
    if not query_anchors or not label_anchors:
        return False
    return not query_anchors.isdisjoint(label_anchors)


def trait_location_hint_tokens(text: str) -> set[str]:
    norm = normalize_trait_text_for_similarity(text)
    if not norm:
        return set()
    tokens = informative_trait_tokens(norm)
    hints: set[str] = set()
    for tok in tokens:
        if len(tok) < 3:
            continue
        if tok.isdigit():
            continue
        if re.fullmatch(r"[a-tv-z][0-9]{2,4}", tok):
            continue
        if tok in TRAIT_STOP_TOKENS:
            continue
        if tok in DISEASE_HINT_TOKENS:
            continue
        if tok in TRAIT_LOCATION_NONSPECIFIC_TOKENS:
            continue
        hints.add(tok)
    return hints


def trait_specificity_loss_flag(
    *,
    query_text: str,
    icd10: str,
    matched_via: str,
    mapped_id_text: str,
    mapped_label_text: str,
    ontology_terms: dict[str, TraitOntologyTerm],
) -> tuple[str, str]:
    method_key = norm_key(matched_via)
    if method_key not in ICD10_PARENT_FALLBACK_METHOD_KEYS:
        return "no", ""

    strict_icd10 = normalize_icd10_code(icd10)
    if not strict_icd10 or "." not in strict_icd10:
        return "no", ""

    input_location_tokens = trait_location_hint_tokens(query_text)
    if not input_location_tokens:
        return "no", ""

    mapped_location_tokens: set[str] = set()
    mapped_ids = [canonicalize_trait_ontology_id(part) for part in split_multi_ids(mapped_id_text)]
    mapped_ids = [item for item in mapped_ids if item]
    mapped_labels = split_multi_labels(mapped_label_text, expected_n=len(mapped_ids))

    if not mapped_ids:
        mapped_location_tokens |= trait_location_hint_tokens(mapped_label_text)
    else:
        for idx, mapped_id in enumerate(mapped_ids):
            label_text = normalize(mapped_labels[idx] if idx < len(mapped_labels) else "")
            term = ontology_terms.get(mapped_id)
            if label_text:
                mapped_location_tokens |= trait_location_hint_tokens(label_text)
            if term is not None and term.label:
                mapped_location_tokens |= trait_location_hint_tokens(term.label)

    if input_location_tokens.isdisjoint(mapped_location_tokens):
        token_preview = "|".join(sorted(input_location_tokens)[:4])
        return "yes", f"location detail in input not represented by mapped term ({token_preview})"
    return "no", ""


def load_ukb_field_title_index(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}

    field_titles: dict[str, str] = {}
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                field_id = normalize(row.get("field_id", ""))
                title = normalize(row.get("title", ""))
                if not field_id or not title:
                    continue
                field_titles[field_id] = title
    except OSError:
        return {}
    return field_titles


def load_trait_query_inputs(path: Path, *, default_trait_scale: str = "auto") -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    default_scale_norm = normalize_trait_scale(default_trait_scale)

    suffix = path.suffix.lower()
    if suffix in {".txt", ".list"}:
        rows: list[dict[str, str]] = []
        for idx, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            query = normalize(line)
            if not query:
                continue
            rows.append(
                {
                    "row_id": str(idx),
                    "query": query,
                    "input_type": "auto",
                    "trait_scale": default_scale_norm,
                    "source_code": "",
                    "source_data_type": "",
                    "icd10": "",
                    "phecode": "",
                }
            )
        return rows

    delimiter = "\t" if suffix in {".tsv", ".tab"} else ","
    with path.open("r", encoding="utf-8", newline="") as handle:
        parsed_rows = list(csv.reader(handle, delimiter=delimiter))
    if not parsed_rows:
        return []

    if parsed_rows and len(parsed_rows[0]) == 1 and all(len(r) <= 1 for r in parsed_rows):
        first = normalize(parsed_rows[0][0] if parsed_rows[0] else "")
        start_idx = 1 if first.lower() in TRAIT_QUERY_COLUMNS else 0
        rows: list[dict[str, str]] = []
        for idx, row in enumerate(parsed_rows[start_idx:], start=1):
            query = normalize(row[0] if row else "")
            if not query:
                continue
            rows.append(
                {
                    "row_id": str(idx),
                    "query": query,
                    "input_type": "auto",
                    "trait_scale": default_scale_norm,
                    "source_code": "",
                    "source_data_type": "",
                    "icd10": "",
                    "phecode": "",
                }
            )
        return rows

    fieldnames = parsed_rows[0]
    field_lookup = {normalize(field).lower(): field for field in fieldnames}

    query_field = next((field_lookup[name] for name in TRAIT_QUERY_COLUMNS if name in field_lookup), fieldnames[0])
    icd_field = next((field_lookup[name] for name in TRAIT_ICD10_COLUMNS if name in field_lookup), None)
    phe_field = next((field_lookup[name] for name in TRAIT_PHECODE_COLUMNS if name in field_lookup), None)
    type_field = next((field_lookup[name] for name in TRAIT_INPUT_TYPE_COLUMNS if name in field_lookup), None)
    scale_field = next((field_lookup[name] for name in TRAIT_SCALE_COLUMNS if name in field_lookup), None)
    code_field = next((field_lookup[name] for name in TRAIT_CODE_COLUMNS if name in field_lookup), None)
    data_type_field = next((field_lookup[name] for name in TRAIT_DATA_TYPE_COLUMNS if name in field_lookup), None)

    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for idx, row in enumerate(reader, start=1):
            raw_query = normalize(row.get(query_field, ""))
            raw_icd = normalize(row.get(icd_field, "")) if icd_field else ""
            raw_phe = normalize(row.get(phe_field, "")) if phe_field else ""
            input_type = normalize_trait_input_type(row.get(type_field, "")) if type_field else "auto"
            trait_scale = normalize_trait_scale(row.get(scale_field, "")) if scale_field else default_scale_norm
            raw_code = normalize(row.get(code_field, "")) if code_field else ""
            raw_data_type = normalize(row.get(data_type_field, "")) if data_type_field else ""
            inferred_type_from_data_type = normalize_trait_data_type(raw_data_type)

            if input_type == "auto" and inferred_type_from_data_type != "auto":
                input_type = inferred_type_from_data_type

            normalized_row_icd = normalize_icd10_code(raw_icd) if raw_icd else ""
            if input_type == "icd10":
                if is_strict_icd10_code(normalized_row_icd):
                    raw_icd = normalized_row_icd
                else:
                    raw_icd = ""
                    icd_from_code = normalize_icd10_code(raw_code)
                    if is_strict_icd10_code(icd_from_code):
                        raw_icd = icd_from_code
                    else:
                        icd_from_query = normalize_icd10_code(raw_query)
                        if is_strict_icd10_code(icd_from_query):
                            raw_icd = icd_from_query
            elif normalized_row_icd:
                raw_icd = normalized_row_icd
            if input_type == "phecode" and not raw_phe:
                raw_phe = normalize_phecode(raw_code or raw_query)
            if input_type == "trait_text" and not raw_query:
                raw_query = raw_icd or raw_phe

            if not raw_query and not raw_icd and not raw_phe:
                continue

            if input_type == "auto":
                normalized_icd10 = normalize_icd10_code(raw_icd or raw_query)
                if is_strict_icd10_code(normalized_icd10):
                    input_type = "icd10"
                    if not raw_icd:
                        raw_icd = normalized_icd10
                elif extract_ukb_field_ids(raw_query):
                    input_type = "ukb_field"
                elif looks_like_phecode_token(raw_phe or raw_query):
                    input_type = "phecode"
                    if not raw_phe:
                        raw_phe = normalize_phecode(raw_query)
                else:
                    input_type = "trait_text"

            if input_type == "ukb_field":
                field_code = normalize(raw_code)
                if not field_code:
                    code_match = re.match(r"^\s*([0-9]{1,7})\b", raw_query)
                    field_code = normalize(code_match.group(1)) if code_match else ""
                if field_code and re.fullmatch(r"[0-9]{1,7}", field_code):
                    ukb_token = f"ukb data field {field_code}"
                    if raw_query:
                        if ukb_token not in norm_key(raw_query):
                            raw_query = normalize(f"{raw_query} (UKB data field {field_code})")
                    else:
                        raw_query = normalize(f"UKB data field {field_code}")

            rows.append(
                {
                    "row_id": str(idx),
                    "query": raw_query,
                    "input_type": input_type,
                    "trait_scale": trait_scale,
                    "source_code": raw_code,
                    "source_data_type": raw_data_type,
                    "icd10": raw_icd,
                    "phecode": raw_phe,
                }
            )
    return rows


def is_ukb_field_supplement_publication(value: str) -> bool:
    return normalize(value).lower().startswith(normalize(UKB_FIELD_SUPPLEMENT_PUBLICATION).lower())


def is_icd10_supplement_publication(value: str) -> bool:
    return normalize(value).lower().startswith(normalize(ICD10_SUPPLEMENT_PUBLICATION).lower())


def load_trait_cache_index(path: Path, extra_paths: list[Path] | None = None) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Trait cache not found: {path}")

    cache_paths: list[Path] = [path]
    for extra in extra_paths or []:
        if extra in cache_paths:
            continue
        if extra.exists():
            cache_paths.append(extra)

    records: list[TraitCacheRecord] = []
    icd10_index: dict[str, list[int]] = {}
    icd10_label_index: dict[str, set[str]] = {}
    phecode_index: dict[str, list[int]] = {}
    ukb_field_index: dict[str, list[int]] = {}
    text_exact_index: dict[str, list[int]] = {}
    token_index: dict[str, set[int]] = {}
    token_freq: Counter[str] = Counter()
    rownum = 0

    for cache_path in cache_paths:
        with cache_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                rownum += 1
                reported_trait = normalize(row.get("Curated reported trait") or row.get("DISEASE/TRAIT") or "")
                lookup_text = normalize(row.get("Lookup text") or row.get("DISEASE/TRAIT") or reported_trait)
                icd10 = normalize_icd10_code(row.get("ICD10") or "")
                phecode = normalize_phecode(row.get("PheCode") or "")
                publication = normalize(row.get("Publication") or "")
                is_ukb_field_supplement = is_ukb_field_supplement_publication(publication)
                is_icd10_supplement = is_icd10_supplement_publication(publication)

                raw_ids = split_multi_ids(row.get("main ontology URI(s)") or row.get("MAPPED_TRAIT_URI") or "")
                if not raw_ids:
                    continue
                canon_ids: list[str] = []
                for item in raw_ids:
                    oid = canonicalize_trait_ontology_id(item)
                    if not oid:
                        continue
                    if oid not in canon_ids:
                        canon_ids.append(oid)
                if not canon_ids:
                    continue

                raw_labels = split_multi_labels(
                    row.get("main ontology label(s)") or row.get("MAPPED_TRAIT") or "",
                    expected_n=len(canon_ids),
                )
                labels: list[str] = []
                for idx, oid in enumerate(canon_ids):
                    mapped_label = raw_labels[idx] if idx < len(raw_labels) else ""
                    labels.append(normalize(mapped_label))

                record = TraitCacheRecord(
                    rownum=rownum,
                    reported_trait=reported_trait,
                    lookup_text=lookup_text,
                    icd10=icd10,
                    phecode=phecode,
                    mapped_ids=tuple(canon_ids),
                    mapped_labels=tuple(labels),
                    publication=publication,
                    is_ukb_field_supplement=is_ukb_field_supplement,
                    is_icd10_supplement=is_icd10_supplement,
                )
                rec_idx = len(records)
                records.append(record)

                if not is_ukb_field_supplement:
                    if icd10:
                        icd10_index.setdefault(icd10, []).append(rec_idx)
                        label_hints: set[str] = set()
                        for text_value in (reported_trait, lookup_text):
                            parsed_code, parsed_label = parse_icd10_labeled_trait_text(text_value)
                            if parsed_code == icd10 and parsed_label:
                                label_hints.add(parsed_label)
                        if lookup_text and lookup_text != icd10:
                            label_hints.add(lookup_text)
                        for label_hint in label_hints:
                            if label_hint:
                                icd10_label_index.setdefault(icd10, set()).add(label_hint)
                    if phecode:
                        phecode_index.setdefault(phecode, []).append(rec_idx)

                ukb_ids = set(extract_ukb_field_ids(reported_trait)) | set(extract_ukb_field_ids(lookup_text))
                explicit_ukb_field = normalize(row.get("UKBB data field") or row.get("UKB data field") or "")
                if explicit_ukb_field and explicit_ukb_field not in {"-", "NA", "N/A"}:
                    if re.fullmatch(r"\d{1,7}", explicit_ukb_field):
                        ukb_ids.add(explicit_ukb_field)
                    else:
                        ukb_ids |= set(extract_ukb_field_ids(explicit_ukb_field))
                for ukb_field_id in ukb_ids:
                    ukb_field_index.setdefault(ukb_field_id, []).append(rec_idx)

                if not is_ukb_field_supplement and not is_icd10_supplement:
                    for text_value in {reported_trait, lookup_text}:
                        text_norm = norm_key(text_value)
                        if not text_norm:
                            continue
                        text_exact_index.setdefault(text_norm, []).append(rec_idx)
                        for tok in informative_trait_tokens(text_norm):
                            token_index.setdefault(tok, set()).add(rec_idx)
                            token_freq[tok] += 1

    return {
        "records": records,
        "icd10_index": icd10_index,
        "icd10_label_index": {k: tuple(sorted(v)) for k, v in icd10_label_index.items()},
        "phecode_index": phecode_index,
        "ukb_field_index": ukb_field_index,
        "text_exact_index": text_exact_index,
        "token_index": token_index,
        "token_freq": token_freq,
    }


def load_trait_ontology_index(obo_path: Path) -> dict[str, Any]:
    if not obo_path.exists():
        raise FileNotFoundError(f"EFO OBO not found: {obo_path}")

    terms: dict[str, TraitOntologyTerm] = {}
    obsolete_terms: dict[str, TraitObsoleteTerm] = {}
    exact_index: dict[str, set[str]] = {}
    token_index: dict[str, set[str]] = {}
    token_freq: Counter[str] = Counter()
    icd10_xref_index: dict[str, set[str]] = {}
    xref_label_hints: dict[str, str] = {}
    term_parents: dict[str, tuple[str, ...]] = {}

    current_id = ""
    current_label = ""
    current_synonyms: list[tuple[str, str]] = []
    current_obsolete = False
    current_replaced_by: list[str] = []
    current_consider: list[str] = []
    current_icd10_xrefs: set[str] = set()
    current_preferred_xrefs: set[str] = set()
    current_parents: list[str] = []
    in_term = False

    def flush_term() -> None:
        nonlocal current_id, current_label, current_synonyms, current_obsolete
        nonlocal current_replaced_by, current_consider
        nonlocal current_icd10_xrefs
        nonlocal current_preferred_xrefs
        nonlocal current_parents

        if not current_id:
            current_id = ""
            current_label = ""
            current_synonyms = []
            current_obsolete = False
            current_replaced_by = []
            current_consider = []
            current_icd10_xrefs = set()
            current_preferred_xrefs = set()
            current_parents = []
            return

        prefix = trait_ontology_prefix(current_id)
        prefix_allowed = prefix in TRAIT_PREFERRED_PREFIX_ORDER

        if current_obsolete:
            if prefix_allowed:
                replaced_by: list[str] = []
                for item in current_replaced_by:
                    oid = canonicalize_trait_ontology_id(item)
                    if not oid:
                        continue
                    if trait_ontology_prefix(oid) not in TRAIT_PREFERRED_PREFIX_ORDER:
                        continue
                    if oid not in replaced_by:
                        replaced_by.append(oid)

                consider: list[str] = []
                for item in current_consider:
                    oid = canonicalize_trait_ontology_id(item)
                    if not oid:
                        continue
                    if trait_ontology_prefix(oid) not in TRAIT_PREFERRED_PREFIX_ORDER:
                        continue
                    if oid not in consider:
                        consider.append(oid)

                obsolete_terms[current_id] = TraitObsoleteTerm(
                    term_id=current_id,
                    label=normalize(current_label),
                    replaced_by=tuple(replaced_by),
                    consider=tuple(consider),
                )
            current_id = ""
            current_label = ""
            current_synonyms = []
            current_obsolete = False
            current_replaced_by = []
            current_consider = []
            current_icd10_xrefs = set()
            current_preferred_xrefs = set()
            current_parents = []
            return

        if not current_label:
            current_id = ""
            current_label = ""
            current_synonyms = []
            current_obsolete = False
            current_replaced_by = []
            current_consider = []
            current_icd10_xrefs = set()
            current_preferred_xrefs = set()
            current_parents = []
            return

        if prefix_allowed:
            normalized_label = normalize(current_label)
            label_key = norm_key(normalized_label)
            synonym_values: list[str] = []
            seen_synonym_values: set[str] = set()
            exact_synonym_keys: set[str] = set()
            narrow_synonym_keys: set[str] = set()
            related_synonym_keys: set[str] = set()
            broad_synonym_keys: set[str] = set()
            for raw_synonym, raw_scope in current_synonyms:
                syn_value = normalize(raw_synonym)
                if not syn_value:
                    continue
                if syn_value not in seen_synonym_values:
                    seen_synonym_values.add(syn_value)
                    synonym_values.append(syn_value)
                syn_key = norm_key(syn_value)
                if not syn_key:
                    continue
                if raw_scope == "EXACT":
                    exact_synonym_keys.add(syn_key)
                elif raw_scope == "NARROW":
                    narrow_synonym_keys.add(syn_key)
                elif raw_scope == "RELATED":
                    related_synonym_keys.add(syn_key)
                elif raw_scope == "BROAD":
                    broad_synonym_keys.add(syn_key)
            term = TraitOntologyTerm(
                term_id=current_id,
                label=normalized_label,
                synonyms=tuple(synonym_values),
                label_key=label_key,
                exact_synonym_keys=frozenset(exact_synonym_keys),
                narrow_synonym_keys=frozenset(narrow_synonym_keys),
                related_synonym_keys=frozenset(related_synonym_keys),
                broad_synonym_keys=frozenset(broad_synonym_keys),
            )
            terms[current_id] = term
            parents: list[str] = []
            for raw_parent in current_parents:
                parent_id = canonicalize_trait_ontology_id(raw_parent)
                if not parent_id:
                    continue
                if parent_id not in parents:
                    parents.append(parent_id)
            term_parents[current_id] = tuple(parents)

            phrases = {term.label, *term.synonyms}
            for phrase in phrases:
                key = norm_key(phrase)
                if not key:
                    continue
                exact_index.setdefault(key, set()).add(current_id)
                for tok in informative_trait_tokens(key):
                    token_index.setdefault(tok, set()).add(current_id)
                    token_freq[tok] += 1

        icd10_targets = set(current_preferred_xrefs)
        if prefix_allowed:
            icd10_targets.add(current_id)
        if icd10_targets:
            for code in sorted(current_icd10_xrefs):
                icd10_xref_index.setdefault(code, set()).update(icd10_targets)
        if current_label:
            normalized_label = normalize(current_label)
            for xref_id in current_preferred_xrefs:
                if xref_id not in terms and normalized_label:
                    xref_label_hints.setdefault(xref_id, normalized_label)

        current_id = ""
        current_label = ""
        current_synonyms = []
        current_obsolete = False
        current_replaced_by = []
        current_consider = []
        current_icd10_xrefs = set()
        current_preferred_xrefs = set()
        current_parents = []

    with obo_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                flush_term()
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                flush_term()
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                current_id = canonicalize_trait_ontology_id(line[4:])
            elif line.startswith("name: "):
                current_label = line[6:].strip()
            elif line.startswith("synonym: "):
                match = re.match(r'^synonym:\s+"(.+?)"\s+(EXACT|BROAD|NARROW|RELATED)\b', line)
                if match:
                    current_synonyms.append((match.group(1).strip(), match.group(2).upper()))
            elif line.startswith("replaced_by: "):
                current_replaced_by.append(line.split(":", 1)[1].strip())
            elif line.startswith("consider: "):
                current_consider.append(line.split(":", 1)[1].strip())
            elif line.startswith("is_obsolete: "):
                current_obsolete = line.split(":", 1)[1].strip().lower() == "true"
            elif line.startswith("xref: "):
                xref_text = line[6:].strip()
                xref_token = xref_text.split()[0]
                xref_id = canonicalize_trait_ontology_id(xref_token)
                if trait_ontology_prefix(xref_id) in TRAIT_PREFERRED_PREFIX_ORDER:
                    current_preferred_xrefs.add(xref_id)
                icd10_code = extract_icd10_code_from_xref(xref_text)
                if icd10_code:
                    current_icd10_xrefs.add(icd10_code)
            elif line.startswith("property_value: "):
                property_text = line
                mondo_id = extract_mondo_id(property_text)
                if mondo_id and trait_ontology_prefix(mondo_id) in TRAIT_PREFERRED_PREFIX_ORDER:
                    current_preferred_xrefs.add(mondo_id)
                icd10_code = extract_icd10_code_from_xref(property_text)
                if icd10_code:
                    current_icd10_xrefs.add(icd10_code)
            elif line.startswith("is_a: "):
                parent_part = line[6:].split("!", 1)[0].strip()
                if parent_part:
                    current_parents.append(parent_part)

    flush_term()
    for term_id, label in xref_label_hints.items():
        if term_id in terms or not label:
            continue
        terms[term_id] = TraitOntologyTerm(
            term_id=term_id,
            label=label,
            synonyms=(),
            label_key=norm_key(label),
        )

    measurement_root = canonicalize_trait_ontology_id(DEFAULT_MEASUREMENT_ROOT)
    measurement_cache: dict[str, bool] = {}
    visiting: set[str] = set()

    def is_measurement_descendant(term_id: str) -> bool:
        if term_id == measurement_root:
            return True
        if term_id in measurement_cache:
            return measurement_cache[term_id]
        if term_id in visiting:
            return False
        visiting.add(term_id)
        parents = term_parents.get(term_id, ())
        result = any(is_measurement_descendant(parent_id) for parent_id in parents)
        visiting.remove(term_id)
        measurement_cache[term_id] = result
        return result

    measurement_term_ids = {
        term_id
        for term_id in terms
        if is_measurement_descendant(term_id)
    }
    return {
        "terms": terms,
        "obsolete_terms": obsolete_terms,
        "exact_index": exact_index,
        "token_index": token_index,
        "token_freq": token_freq,
        "icd10_xref_index": icd10_xref_index,
        "measurement_term_ids": measurement_term_ids,
    }


def build_trait_cache_term_resolution(
    trait_cache_index: dict[str, Any],
    ontology_index: dict[str, Any],
    *,
    issue_sample_limit: int = 8,
) -> tuple[dict[int, tuple[tuple[str, ...], tuple[str, ...]]], dict[str, int], list[str]]:
    records: list[TraitCacheRecord] = trait_cache_index.get("records", [])
    ontology_terms: dict[str, TraitOntologyTerm] = ontology_index.get("terms", {})
    obsolete_terms: dict[str, TraitObsoleteTerm] = ontology_index.get("obsolete_terms", {})

    resolution: dict[int, tuple[tuple[str, ...], tuple[str, ...]]] = {}
    stats: dict[str, int] = {
        "records_total": len(records),
        "rows_with_obsolete_ids": 0,
        "rows_with_unresolved_ids": 0,
        "rows_without_active_terms": 0,
        "ids_total": 0,
        "ids_active": 0,
        "ids_obsolete": 0,
        "ids_obsolete_replaced_by": 0,
        "ids_obsolete_consider_single": 0,
        "ids_obsolete_unresolved": 0,
        "ids_unknown": 0,
    }
    issues: list[str] = []

    for rec in records:
        resolved_ids: list[str] = []
        unresolved_mondo_ids: list[str] = []
        row_has_obsolete = False
        row_has_unresolved = False

        for raw_id in rec.mapped_ids:
            oid = canonicalize_trait_ontology_id(raw_id)
            if not oid:
                continue
            stats["ids_total"] += 1
            if oid in ontology_terms:
                stats["ids_active"] += 1
                if oid not in resolved_ids:
                    resolved_ids.append(oid)
                continue

            obsolete_meta = obsolete_terms.get(oid)
            if obsolete_meta is not None:
                row_has_obsolete = True
                stats["ids_obsolete"] += 1

                replacements = [rid for rid in obsolete_meta.replaced_by if rid in ontology_terms]
                if replacements:
                    stats["ids_obsolete_replaced_by"] += 1
                    for rid in replacements:
                        if rid not in resolved_ids:
                            resolved_ids.append(rid)
                    continue

                consider = [cid for cid in obsolete_meta.consider if cid in ontology_terms]
                if len(consider) == 1:
                    stats["ids_obsolete_consider_single"] += 1
                    cid = consider[0]
                    if cid not in resolved_ids:
                        resolved_ids.append(cid)
                    continue

                stats["ids_obsolete_unresolved"] += 1
                row_has_unresolved = True
                if len(issues) < issue_sample_limit:
                    issues.append(
                        f"row {rec.rownum}: unresolved obsolete id {oid} "
                        f"for '{rec.reported_trait or rec.lookup_text}'"
                    )
                continue

            stats["ids_unknown"] += 1
            row_has_unresolved = True
            if trait_ontology_prefix(oid) == "MONDO":
                unresolved_mondo_ids.append(oid)
            if len(issues) < issue_sample_limit:
                issues.append(
                    f"row {rec.rownum}: id {oid} not found in ontology "
                    f"for '{rec.reported_trait or rec.lookup_text}'"
                )

        if row_has_obsolete:
            stats["rows_with_obsolete_ids"] += 1
        if row_has_unresolved:
            stats["rows_with_unresolved_ids"] += 1

        if not resolved_ids:
            # Keep unresolved MONDO-only rows as explicit review-required
            # candidates so output can carry `term not in EFO` QC flags.
            if unresolved_mondo_ids:
                fallback_labels: list[str] = []
                for unresolved_id in unresolved_mondo_ids:
                    fallback_label = ""
                    for idx, mapped_id in enumerate(rec.mapped_ids):
                        if canonicalize_trait_ontology_id(mapped_id) == unresolved_id:
                            if idx < len(rec.mapped_labels):
                                fallback_label = normalize(rec.mapped_labels[idx])
                            break
                    fallback_labels.append(fallback_label or unresolved_id)
                resolution[rec.rownum] = (tuple(unresolved_mondo_ids), tuple(fallback_labels))
                continue
            stats["rows_without_active_terms"] += 1
            continue

        resolved_labels = tuple(ontology_terms[oid].label for oid in resolved_ids)
        resolution[rec.rownum] = (tuple(resolved_ids), resolved_labels)

    return resolution, stats, issues


def pick_seed_tokens(tokens: set[str], token_freq: Counter[str], max_tokens: int = 6) -> list[str]:
    if not tokens:
        return []
    items = list(tokens)
    items.sort(key=lambda tok: (token_freq.get(tok, 10**9), len(tok), tok))
    return items[:max_tokens]


@lru_cache(maxsize=300_000)
def _trait_similarity_candidate_profile(candidate_text: str) -> tuple[str, frozenset[str]]:
    cnorm = normalize_trait_text_for_similarity(candidate_text)
    if not cnorm:
        return "", frozenset()
    return cnorm, frozenset(informative_trait_tokens(cnorm))


@lru_cache(maxsize=1_000_000)
def _sequence_ratio(a: str, b: str) -> float:
    return difflib.SequenceMatcher(None, a, b).ratio()


def trait_similarity_score(query: str, query_toks: set[str], candidate_text: str) -> float:
    qnorm = normalize_trait_text_for_similarity(query)
    cnorm, c_toks_frozen = _trait_similarity_candidate_profile(candidate_text)
    if not qnorm or not cnorm:
        return 0.0
    negation_mismatch = trait_negation_mismatch(qnorm, cnorm)
    c_toks = set(c_toks_frozen)
    overlap_tokens = query_toks & c_toks
    overlap = len(overlap_tokens) / max(len(query_toks), 1)
    nonweak_overlap_tokens = trait_nonweak_tokens(overlap_tokens)
    query_anchor_tokens = trait_anchor_tokens(query_toks)
    candidate_anchor_tokens = trait_anchor_tokens(c_toks)
    seq = _sequence_ratio(qnorm, cnorm)
    substring_bonus = 0.1 if len(cnorm) >= 6 and (cnorm in qnorm or qnorm in cnorm) else 0.0
    strict_modifier_mismatches = sum(
        1 for tok in TRAIT_STRICT_MODIFIER_TOKENS if ((tok in query_toks) != (tok in c_toks))
    )
    modifier_penalty = 0.25 * strict_modifier_mismatches
    negation_penalty = 0.12 if negation_mismatch else 0.0
    weak_overlap_penalty = 0.18 if overlap_tokens and not nonweak_overlap_tokens else 0.0
    anchor_penalty = 0.20 if query_anchor_tokens and query_anchor_tokens.isdisjoint(candidate_anchor_tokens) else 0.0
    raw_score = (
        0.6 * overlap
        + 0.4 * seq
        + substring_bonus
        - modifier_penalty
        - negation_penalty
        - weak_overlap_penalty
        - anchor_penalty
    )
    return max(0.0, min(1.0, raw_score))


def trait_modifier_label_mismatch(query_toks: set[str], candidate_label: str) -> bool:
    if not query_toks or not candidate_label:
        return False
    label_norm = normalize_trait_text_for_similarity(candidate_label)
    if not label_norm:
        return False
    label_toks = informative_trait_tokens(label_norm)
    return any((tok in query_toks) != (tok in label_toks) for tok in TRAIT_STRICT_MODIFIER_TOKENS)


def trait_match_scope_rank(
    *,
    query_exact_keys: set[str],
    term: TraitOntologyTerm | None,
    label_text: str,
) -> int:
    if not query_exact_keys:
        return -1
    label_key = norm_key(label_text)
    if label_key and label_key in query_exact_keys:
        return TRAIT_SYNONYM_SCOPE_RANK["LABEL"]
    if term is None:
        return -1
    if query_exact_keys & term.exact_synonym_keys:
        return TRAIT_SYNONYM_SCOPE_RANK["EXACT"]
    if query_exact_keys & term.narrow_synonym_keys:
        return TRAIT_SYNONYM_SCOPE_RANK["NARROW"]
    if query_exact_keys & term.related_synonym_keys:
        return TRAIT_SYNONYM_SCOPE_RANK["RELATED"]
    if query_exact_keys & term.broad_synonym_keys:
        return TRAIT_SYNONYM_SCOPE_RANK["BROAD"]
    return -1


def trait_overspecificity_penalties(
    *,
    query_tokens: set[str],
    label_text: str,
) -> tuple[int, int, int]:
    if not label_text:
        return 0, 0, 0
    label_norm = normalize_trait_text_for_similarity(label_text) or norm_key(label_text)
    label_tokens = informative_trait_tokens(label_norm)
    if not label_tokens:
        return 0, 0, 0
    extra_tokens = label_tokens - query_tokens
    subtype_penalty = len(extra_tokens & TRAIT_OVERSPECIFICITY_HINT_TOKENS)
    extra_penalty = len([tok for tok in extra_tokens if tok not in TRAIT_STOP_TOKENS and len(tok) >= 4])
    overlap = len(label_tokens & query_tokens)
    return subtype_penalty, extra_penalty, overlap


def trait_candidate_lexical_priority(
    candidate: Candidate,
    *,
    query_exact_keys: set[str],
    query_tokens: set[str],
    ontology_terms: dict[str, TraitOntologyTerm],
) -> tuple[int, int, int, int]:
    ids = [canonicalize_trait_ontology_id(item) for item in split_multi_ids(candidate.efo_id)]
    ids = [item for item in ids if item]
    labels = split_multi_labels(candidate.label, expected_n=len(ids))
    if not ids:
        scope_rank = trait_match_scope_rank(
            query_exact_keys=query_exact_keys,
            term=None,
            label_text=normalize(candidate.label),
        )
        subtype_penalty, extra_penalty, overlap = trait_overspecificity_penalties(
            query_tokens=query_tokens,
            label_text=normalize(candidate.label),
        )
        return scope_rank, -subtype_penalty, -extra_penalty, overlap

    best_rank = -1
    best_subtype_penalty = 10**6
    best_extra_penalty = 10**6
    best_overlap = -1
    for idx, term_id in enumerate(ids):
        term = ontology_terms.get(term_id)
        label_text = normalize(labels[idx] if idx < len(labels) else "")
        if not label_text and term is not None:
            label_text = term.label
        scope_rank = trait_match_scope_rank(
            query_exact_keys=query_exact_keys,
            term=term,
            label_text=label_text,
        )
        subtype_penalty, extra_penalty, overlap = trait_overspecificity_penalties(
            query_tokens=query_tokens,
            label_text=label_text,
        )
        rank_tuple = (scope_rank, -subtype_penalty, -extra_penalty, overlap)
        best_tuple = (best_rank, -best_subtype_penalty, -best_extra_penalty, best_overlap)
        if rank_tuple > best_tuple:
            best_rank = scope_rank
            best_subtype_penalty = subtype_penalty
            best_extra_penalty = extra_penalty
            best_overlap = overlap
    return best_rank, -best_subtype_penalty, -best_extra_penalty, best_overlap


def candidate_is_non_disease_entity(
    candidate: Candidate,
    *,
    ontology_terms: dict[str, TraitOntologyTerm],
) -> bool:
    ids = [canonicalize_trait_ontology_id(item) for item in split_multi_ids(candidate.efo_id)]
    ids = [item for item in ids if item]
    if not ids:
        return False
    prefixes = {trait_ontology_prefix(item) for item in ids if item}
    if not prefixes or not prefixes.issubset(TRAIT_NON_DISEASE_ENTITY_PREFIXES):
        return False
    if has_disease_hint(candidate.label):
        return False
    for term_id in ids:
        term = ontology_terms.get(term_id)
        if term is not None and has_disease_hint(term.label):
            return False
    return True


def collapse_trait_cache_candidates(
    records: list[TraitCacheRecord],
    *,
    method: str,
    score: float,
    matched_on: str,
    validated_if_unique: bool,
    record_resolver: Any | None = None,
) -> list[Candidate]:
    grouped: dict[tuple[str, str], list[TraitCacheRecord]] = {}
    for record in records:
        if record_resolver is not None:
            resolved = record_resolver(record)
            if not resolved:
                continue
            resolved_ids, resolved_labels = resolved
            mapped_ids = "|".join(resolved_ids)
            mapped_labels = "|".join(resolved_labels)
        else:
            mapped_ids = "|".join(record.mapped_ids)
            mapped_labels = "|".join(record.mapped_labels)
        grouped.setdefault((mapped_ids, mapped_labels), []).append(record)

    unique_count = len(grouped)
    if unique_count == 0:
        return []
    candidates: list[Candidate] = []
    for (mapped_ids, mapped_labels), rows in grouped.items():
        source_rows = ",".join(str(row.rownum) for row in rows[:10])
        publication_set = sorted({row.publication for row in rows if row.publication})
        publication_hint = "|".join(publication_set[:3])
        evidence = f"{matched_on} cache rows={source_rows}"
        if publication_hint:
            evidence += f"; publications={publication_hint}"
        if unique_count > 1:
            evidence += "; multiple cache mappings for key"
        candidates.append(
            Candidate(
                efo_id=mapped_ids,
                label=mapped_labels,
                score=max(0.0, min(1.0, score)),
                matched_via=method,
                evidence=evidence,
                is_validated=validated_if_unique and unique_count == 1,
            )
        )
    return candidates


def pick_preferred_mapping_pair(mapped_ids: str, mapped_labels: str) -> tuple[str, str] | None:
    ids = [canonicalize_trait_ontology_id(item) for item in split_multi_ids(mapped_ids)]
    ids = [item for item in ids if item]
    if not ids:
        return None
    labels = split_multi_labels(mapped_labels, expected_n=len(ids))
    pairs: list[tuple[str, str]] = []
    for idx, item in enumerate(ids):
        label = normalize(labels[idx] if idx < len(labels) else "")
        pairs.append((item, label))
    pairs.sort(
        key=lambda pair: (
            TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(pair[0]), 999),
            pair[0],
        )
    )
    return pairs[0]


def pick_preferred_ontology_term(
    term_ids: set[str],
    *,
    ontology_terms: dict[str, TraitOntologyTerm],
) -> tuple[str, str] | None:
    if not term_ids:
        return None
    ordered = sorted(
        (term_id for term_id in term_ids if term_id in ontology_terms),
        key=lambda term_id: (
            TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(term_id), 999),
            term_id,
        ),
    )
    if not ordered:
        return None
    term_id = ordered[0]
    return term_id, ontology_terms[term_id].label


def is_generic_ukb_measurement_mapping(mapped_id: str, mapped_label: str) -> bool:
    oid = canonicalize_trait_ontology_id(mapped_id)
    label = norm_key(mapped_label)
    return oid == UKB_FIELD_SUPPLEMENT_GENERIC_MEASUREMENT_ID or label == "measurement"


def candidate_has_generic_measurement_mapping(mapped_ids: str, mapped_labels: str) -> bool:
    ids = [canonicalize_trait_ontology_id(item) for item in split_multi_ids(mapped_ids)]
    ids = [item for item in ids if item]
    if not ids:
        return is_generic_ukb_measurement_mapping("", mapped_labels)
    labels = split_multi_labels(mapped_labels, expected_n=len(ids))
    for idx, term_id in enumerate(ids):
        label = normalize(labels[idx] if idx < len(labels) else "")
        if is_generic_ukb_measurement_mapping(term_id, label):
            return True
    return False


def trait_qc_risk_reasons(row: dict[str, str]) -> list[str]:
    reasons: list[str] = []
    validation = norm_key(row.get("validation", ""))
    matched_via = norm_key(row.get("matched_via", ""))
    mapped_label = normalize(row.get("mapped_trait_label", ""))

    try:
        confidence = float(row.get("confidence") or "0")
    except ValueError:
        confidence = 0.0

    if row.get("term_not_in_efo"):
        reasons.append("term_not_in_efo")
    if norm_key(row.get("negation", "")) == "yes" and validation != "not_mapped":
        reasons.append("negation_with_mapping")
    if validation == "validated" and matched_via in {"cache_fuzzy_text", "efo_obo_fuzzy"}:
        reasons.append("validated_fuzzy_match")
    if validation == "validated" and confidence < 0.90:
        reasons.append("validated_low_confidence")
    if norm_key(row.get("specificity_loss_flag", "")) == "yes":
        reasons.append("specificity_loss")
    if is_generic_trait_label(mapped_label) and validation in {"validated", "review_required"}:
        reasons.append("generic_mapped_label")
    if matched_via in ICD10_PARENT_FALLBACK_METHOD_KEYS:
        reasons.append("icd10_parent_fallback")
    if matched_via == "cache_exact_ukb_field_supplement" and is_generic_trait_label(mapped_label):
        reasons.append("ukb_supplement_generic")
    return reasons


def write_trait_qc_risk_tsv(mapped_output: Path, risk_output: Path) -> int:
    risk_rows: list[dict[str, str]] = []
    with mapped_output.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            reasons = trait_qc_risk_reasons(row)
            if not reasons:
                continue
            risk_rows.append(
                {
                    "input_row_id": normalize(row.get("input_row_id", "")),
                    "input_query": normalize(row.get("input_query", "")),
                    "input_source_data_type": normalize(row.get("input_source_data_type", "")),
                    "input_type": normalize(row.get("input_type", "")),
                    "mapped_trait_id": normalize(row.get("mapped_trait_id", "")),
                    "mapped_trait_label": normalize(row.get("mapped_trait_label", "")),
                    "validation": normalize(row.get("validation", "")),
                    "matched_via": normalize(row.get("matched_via", "")),
                    "confidence": normalize(row.get("confidence", "")),
                    "specificity_loss_flag": normalize(row.get("specificity_loss_flag", "")),
                    "specificity_loss_notes": normalize(row.get("specificity_loss_notes", "")),
                    "risk_reasons": "|".join(sorted(set(reasons))),
                    "evidence": normalize(row.get("evidence", "")),
                }
            )

    risk_output.parent.mkdir(parents=True, exist_ok=True)
    with risk_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRAIT_QC_RISK_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(risk_rows)
    return len(risk_rows)


def parse_mondo_sssom_icd10_mappings(path: Path) -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = {}
    if not path.exists():
        return mapping
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            subject_id = normalize(row.get("subject_id", ""))
            object_id = normalize(row.get("object_id", ""))
            if not subject_id or not object_id:
                continue
            subject_icd10 = extract_icd10_code_from_xref(subject_id)
            object_icd10 = extract_icd10_code_from_xref(object_id)
            subject_mondo = extract_mondo_id(subject_id)
            object_mondo = extract_mondo_id(object_id)

            if subject_icd10 and object_mondo:
                mapping.setdefault(subject_icd10, set()).add(object_mondo)
            if object_icd10 and subject_mondo:
                mapping.setdefault(object_icd10, set()).add(subject_mondo)
    return mapping


def parse_cms_icd10_label_rows(zip_path: Path, *, source_url: str) -> tuple[list[dict[str, str]], str]:
    with zipfile.ZipFile(zip_path) as zf:
        member_names = [
            name
            for name in zf.namelist()
            if name.lower().endswith(".txt")
            and "order" in Path(name).name.lower()
            and "addenda" not in Path(name).name.lower()
        ]
        if not member_names:
            raise ValueError(f"No CMS ICD10 tabular-order text file found in {zip_path}")
        member_names.sort(key=lambda name: ("icd10cm_order" not in Path(name).name.lower(), len(name), name))
        member_name = member_names[0]
        raw_lines = zf.read(member_name).decode("utf-8", errors="replace").splitlines()

    rows: list[dict[str, str]] = []
    seen_codes: set[str] = set()
    for raw_line in raw_lines:
        line = raw_line.rstrip("\r\n")
        if len(line) < 17:
            continue
        raw_code = normalize(line[6:13])
        if not raw_code:
            continue
        label = normalize(line[77:] if len(line) > 77 else "") or normalize(line[16:76])
        code = normalize_icd10_code(raw_code)
        if not code or not is_strict_icd10_code(code) or not label or code in seen_codes:
            continue
        seen_codes.add(code)
        rows.append(
            {
                "ICD10": code,
                "label": label,
                "source": "CMS ICD-10-CM tabular order",
                "source_url": source_url,
            }
        )

    if not rows:
        raise ValueError(f"No ICD10 label rows parsed from {zip_path}")
    return rows, member_name


def build_icd10_label_cache(
    *,
    output_path: Path,
    source_url: str,
    timeout: float,
) -> dict[str, Any]:
    url = normalize(source_url)
    if not url:
        raise ValueError("ICD10 label source URL is required")

    with tempfile.TemporaryDirectory(prefix="icd10-label-cache-") as tmpdir:
        archive_path = Path(tmpdir) / Path(urllib.parse.urlparse(url).path).name
        if not archive_path.name:
            archive_path = Path(tmpdir) / "icd10_labels_source.zip"
        download_to_file(url, archive_path, timeout=timeout)
        archive_size_bytes = archive_path.stat().st_size if archive_path.exists() else 0
        rows, member_name = parse_cms_icd10_label_rows(archive_path, source_url=url)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ICD10", "label", "source", "source_url"],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    return {
        "codes": len(rows),
        "output_path": str(output_path),
        "source_url": url,
        "source_member": member_name,
        "archive_size_bytes": archive_size_bytes,
    }


def load_icd10_label_cache(path: Path) -> dict[str, tuple[str, ...]]:
    label_index: dict[str, set[str]] = {}
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            code = normalize_icd10_code(row.get("ICD10") or row.get("icd10") or "")
            label = normalize(row.get("label") or row.get("icd10_label") or row.get("description") or "")
            if not code or not is_strict_icd10_code(code) or not label:
                continue
            label_index.setdefault(code, set()).add(label)
    return {code: tuple(sorted(labels)) for code, labels in label_index.items()}


def build_icd10_supplement_cache(
    *,
    output_path: Path,
    trait_cache_path: Path,
    efo_obo_path: Path,
    mondo_sssom_cache_path: Path,
    mondo_sssom_url: str,
    mondo_download_timeout: float,
) -> dict[str, Any]:
    trait_cache_index = load_trait_cache_index(trait_cache_path)
    base_icd10_codes = set(trait_cache_index.get("icd10_index", {}).keys())
    ontology_index = load_trait_ontology_index(efo_obo_path)
    ontology_terms: dict[str, TraitOntologyTerm] = ontology_index.get("terms", {})
    ontology_icd10_xref_index: dict[str, set[str]] = ontology_index.get("icd10_xref_index", {})

    code_to_term_ids: dict[str, set[str]] = {}
    code_sources: dict[str, set[str]] = {}
    for raw_code, term_ids in ontology_icd10_xref_index.items():
        code = normalize_icd10_code(raw_code)
        if not code:
            continue
        for term_id in term_ids:
            oid = canonicalize_trait_ontology_id(term_id)
            if oid in ontology_terms:
                code_to_term_ids.setdefault(code, set()).add(oid)
                code_sources.setdefault(code, set()).add("efo_obo_xref")

    mondo_codes_added = 0
    mondo_terms_not_in_efo = 0
    mondo_downloaded = False
    mondo_download_error = ""
    mondo_cache_source = "none"
    mondo_cache_used = False
    mondo_mapping_codes = 0
    mondo_url = normalize(mondo_sssom_url)
    if mondo_url:
        try:
            mondo_sssom_cache_path.parent.mkdir(parents=True, exist_ok=True)
            download_to_file(mondo_url, mondo_sssom_cache_path, timeout=mondo_download_timeout)
            mondo_downloaded = True
            mondo_cache_source = "downloaded"
        except Exception as exc:  # noqa: BLE001 - setup should degrade gracefully offline.
            mondo_download_error = str(exc)
    if mondo_sssom_cache_path.exists() and mondo_sssom_cache_path.stat().st_size > 0:
        try:
            mondo_mapping = parse_mondo_sssom_icd10_mappings(mondo_sssom_cache_path)
            mondo_mapping_codes = len(mondo_mapping)
            if mondo_mapping:
                mondo_cache_used = True
                if mondo_cache_source == "none":
                    mondo_cache_source = "cached"
            for code, term_ids in mondo_mapping.items():
                valid_terms: set[str] = set()
                for term_id in term_ids:
                    canonical = canonicalize_trait_ontology_id(term_id)
                    if not canonical:
                        continue
                    if canonical in ontology_terms:
                        valid_terms.add(canonical)
                    elif trait_ontology_prefix(canonical) == "MONDO":
                        mondo_terms_not_in_efo += 1
                if not valid_terms:
                    continue
                before = len(code_to_term_ids.get(code, set()))
                code_to_term_ids.setdefault(code, set()).update(valid_terms)
                after = len(code_to_term_ids.get(code, set()))
                if after > before:
                    mondo_codes_added += 1
                code_sources.setdefault(code, set()).add("mondo_sssom")
        except Exception as exc:  # noqa: BLE001 - setup should degrade gracefully offline.
            if not mondo_download_error:
                mondo_download_error = str(exc)

    method_counts: Counter[str] = Counter()
    supplement_rows: list[dict[str, str]] = []
    total_codes = 0
    base_covered_codes = 0

    def method_for_sources(sources: set[str]) -> str:
        if "mondo_sssom" in sources and "efo_obo_xref" in sources:
            return "icd10_supplement_mondo_sssom_efo_obo_xref"
        if "mondo_sssom" in sources:
            return "icd10_supplement_mondo_sssom"
        return "icd10_supplement_efo_obo_xref"

    sortable_codes = sorted(
        code_to_term_ids.keys(),
        key=lambda code: (code[:1], int(re.sub(r"\D", "", code[:3]) or 0), code),
    )
    for code in sortable_codes:
        total_codes += 1
        if code in base_icd10_codes:
            base_covered_codes += 1
            continue
        selected = pick_preferred_ontology_term(code_to_term_ids.get(code, set()), ontology_terms=ontology_terms)
        if selected is None:
            continue
        mapped_id, mapped_label = selected
        method = method_for_sources(code_sources.get(code, set()))
        method_counts[method] += 1
        supplement_rows.append(
            {
                "UKBB data field": "-",
                "ICD10": code,
                "PheCode": "-",
                "Lookup text": f"ICD10 {code} {mapped_label}",
                "Curated reported trait": f"ICD10 {code}: {mapped_label}",
                "main ontology URI(s)": mapped_id,
                "main ontology label(s)": mapped_label,
                "Publication": ICD10_SUPPLEMENT_PUBLICATION,
                "User submitted trait": "-",
                "qc_status": "review_required",
                "qc_resolution_modes": method,
                "qc_label_match_modes": method,
                "qc_unresolved_ids": "",
                "final_update_applied": "yes",
                "final_update_reason": "icd10_supplement_cache_generated",
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TRAIT_CACHE_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(supplement_rows)

    return {
        "total_codes": total_codes,
        "base_covered_codes": base_covered_codes,
        "supplement_rows": len(supplement_rows),
        "method_efo_obo_xref": method_counts.get("icd10_supplement_efo_obo_xref", 0),
        "method_mondo_sssom": method_counts.get("icd10_supplement_mondo_sssom", 0),
        "method_mondo_sssom_efo_obo_xref": method_counts.get(
            "icd10_supplement_mondo_sssom_efo_obo_xref", 0
        ),
        "mondo_url_used": mondo_url,
        "mondo_cache_path": str(mondo_sssom_cache_path),
        "mondo_cache_source": mondo_cache_source,
        "mondo_cache_used": mondo_cache_used,
        "mondo_mapping_codes": mondo_mapping_codes,
        "mondo_downloaded": mondo_downloaded,
        "mondo_codes_added": mondo_codes_added,
        "mondo_terms_not_in_efo": mondo_terms_not_in_efo,
        "mondo_download_error": mondo_download_error,
    }


def build_ukb_field_supplement_cache(
    *,
    output_path: Path,
    trait_cache_path: Path,
    efo_obo_path: Path,
    ukb_field_catalog_path: Path | None,
    ukb_field_metadata_path: Path | None,
    ukb_category_catalog_path: Path | None,
    ukb_category_tree_path: Path | None,
) -> dict[str, int]:
    if ukb_field_metadata_path is None or not ukb_field_metadata_path.exists():
        raise FileNotFoundError(
            "UKB field metadata is required to build supplemental field cache: "
            f"{ukb_field_metadata_path}"
        )

    trait_cache_index = load_trait_cache_index(trait_cache_path)
    ontology_index = load_trait_ontology_index(efo_obo_path)
    trait_cache_resolution, _, _ = build_trait_cache_term_resolution(trait_cache_index, ontology_index)
    records: list[TraitCacheRecord] = trait_cache_index.get("records", [])
    ukb_field_index: dict[str, list[int]] = trait_cache_index.get("ukb_field_index", {})
    text_exact_index: dict[str, list[int]] = trait_cache_index.get("text_exact_index", {})
    ontology_terms: dict[str, TraitOntologyTerm] = ontology_index.get("terms", {})
    ontology_exact_index: dict[str, set[str]] = ontology_index.get("exact_index", {})
    ontology_token_index: dict[str, set[str]] = ontology_index.get("token_index", {})
    ontology_token_freq: Counter[str] = ontology_index.get("token_freq", Counter())

    ukb_field_titles = load_ukb_field_title_index(ukb_field_catalog_path)
    metadata_field_titles = load_ukb_field_title_index(ukb_field_metadata_path)
    for field_id, title in metadata_field_titles.items():
        if field_id in ukb_field_titles:
            continue
        ukb_field_titles[field_id] = title

    ukb_field_to_category, ukb_category_to_fields = load_ukb_field_category_index(ukb_field_metadata_path)
    ukb_category_titles = load_ukb_category_title_index(ukb_category_catalog_path)
    ukb_category_parents = load_ukb_category_parent_index(ukb_category_tree_path)

    def resolve_record(record: TraitCacheRecord) -> tuple[tuple[str, ...], tuple[str, ...]] | None:
        resolved = trait_cache_resolution.get(record.rownum)
        if resolved is not None:
            return resolved
        if record.mapped_ids:
            return record.mapped_ids, record.mapped_labels
        return None

    def choose_mapping_from_records(
        matched_records: list[TraitCacheRecord],
        *,
        method: str,
        matched_on: str,
    ) -> tuple[str, str, str] | None:
        candidates = collapse_trait_cache_candidates(
            matched_records,
            method=method,
            score=0.9,
            matched_on=matched_on,
            validated_if_unique=False,
            record_resolver=resolve_record,
        )
        if len(candidates) != 1:
            return None
        pair = pick_preferred_mapping_pair(candidates[0].efo_id, candidates[0].label)
        if pair is None:
            return None
        mapped_id, mapped_label = pair
        return mapped_id, mapped_label, method

    def choose_mapping_from_ontology_fuzzy(
        text: str,
        *,
        method: str,
        min_score: float = 0.88,
        min_margin: float = 0.04,
    ) -> tuple[str, str, str] | None:
        query_norm = normalize_trait_text_for_similarity(text)
        query_tokens = informative_trait_tokens(query_norm)
        if not query_norm or not query_tokens:
            return None
        seed_tokens = pick_seed_tokens(query_tokens, ontology_token_freq, max_tokens=7)
        candidate_term_ids: set[str] = set()
        for tok in seed_tokens:
            candidate_term_ids |= ontology_token_index.get(tok, set())
        if not candidate_term_ids:
            return None

        scored_terms: list[tuple[float, str]] = []
        for term_id in candidate_term_ids:
            term = ontology_terms.get(term_id)
            if term is None:
                continue
            best_score = trait_similarity_score(query_norm, query_tokens, term.label)
            for syn in term.synonyms[:10]:
                best_score = max(best_score, trait_similarity_score(query_norm, query_tokens, syn))
            if best_score <= 0.0:
                continue
            scored_terms.append((best_score, term_id))
        if not scored_terms:
            return None

        scored_terms.sort(
            key=lambda item: (
                item[0],
                -TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(item[1]), 999),
                item[1],
            ),
            reverse=True,
        )
        top_score, top_id = scored_terms[0]
        second_score = scored_terms[1][0] if len(scored_terms) > 1 else 0.0
        if top_score < min_score:
            return None
        if (top_score - second_score) < min_margin:
            return None
        top_term = ontology_terms.get(top_id)
        if top_term is None:
            return None
        return top_term.term_id, top_term.label, method

    def category_context(category_id: str) -> tuple[str, tuple[str, ...]]:
        category_title = normalize(ukb_category_titles.get(category_id, ""))
        path_titles = ukb_category_path_titles(
            category_id,
            category_titles=ukb_category_titles,
            category_parents=ukb_category_parents,
        )
        return category_title, path_titles

    category_default_mapping: dict[str, tuple[str, str, str]] = {}
    for category_id in sorted({normalize(cat) for cat in ukb_field_to_category.values() if normalize(cat)}):
        category_title, path_titles = category_context(category_id)
        category_field_ids = tuple(normalize(item) for item in ukb_category_to_fields.get(category_id, ()) if normalize(item))
        matched_record_idxs: list[int] = []
        seen_rec_idxs: set[int] = set()
        for field_id in category_field_ids:
            for rec_idx in ukb_field_index.get(field_id, []):
                record = records[rec_idx]
                if record.is_ukb_field_supplement:
                    continue
                if rec_idx in seen_rec_idxs:
                    continue
                seen_rec_idxs.add(rec_idx)
                matched_record_idxs.append(rec_idx)
        if matched_record_idxs:
            matched_records = [records[idx] for idx in matched_record_idxs]
            chosen = choose_mapping_from_records(
                matched_records,
                method="ukb_field_supplement_category_cache_unique",
                matched_on="UKB category field set",
            )
            if chosen is not None and is_generic_ukb_measurement_mapping(chosen[0], chosen[1]):
                chosen = None
            if chosen is not None:
                category_default_mapping[category_id] = chosen
                continue

        category_exact_keys: list[str] = []
        for variant in ukb_category_query_variants(category_id, category_title, path_titles):
            key = norm_key(variant)
            if key and key not in category_exact_keys:
                category_exact_keys.append(key)
        for title in path_titles:
            key = norm_key(title)
            if key and key not in category_exact_keys:
                category_exact_keys.append(key)

        matched_record_idxs = []
        seen_rec_idxs = set()
        for key in category_exact_keys:
            for rec_idx in text_exact_index.get(key, []):
                record = records[rec_idx]
                if record.is_ukb_field_supplement:
                    continue
                if rec_idx in seen_rec_idxs:
                    continue
                seen_rec_idxs.add(rec_idx)
                matched_record_idxs.append(rec_idx)
        if matched_record_idxs:
            matched_records = [records[idx] for idx in matched_record_idxs]
            chosen = choose_mapping_from_records(
                matched_records,
                method="ukb_field_supplement_category_text_exact",
                matched_on="UKB category title/path text",
            )
            if chosen is not None and is_generic_ukb_measurement_mapping(chosen[0], chosen[1]):
                chosen = None
            if chosen is not None:
                category_default_mapping[category_id] = chosen
                continue

        exact_term_ids: set[str] = set()
        for key in category_exact_keys:
            exact_term_ids |= ontology_exact_index.get(key, set())
        selected_term = pick_preferred_ontology_term(exact_term_ids, ontology_terms=ontology_terms)
        if selected_term is not None:
            term_id, term_label = selected_term
            category_default_mapping[category_id] = (
                term_id,
                term_label,
                "ukb_field_supplement_category_ontology_exact",
            )
            continue

        category_text = " ".join([category_title, *path_titles]).lower()
        looks_disease = (
            has_disease_hint(category_text)
            or any(phrase in category_text for phrase in UKB_DISEASE_CATEGORY_PHRASES)
        )
        if looks_disease:
            category_default_mapping[category_id] = (
                UKB_FIELD_SUPPLEMENT_GENERIC_DISEASE_ID,
                UKB_FIELD_SUPPLEMENT_GENERIC_DISEASE_LABEL,
                "ukb_field_supplement_generic_disease",
            )
        else:
            category_default_mapping[category_id] = (
                UKB_FIELD_SUPPLEMENT_GENERIC_MEASUREMENT_ID,
                UKB_FIELD_SUPPLEMENT_GENERIC_MEASUREMENT_LABEL,
                "ukb_field_supplement_generic_measurement",
            )

    base_covered_fields = set(ukb_field_index.keys())
    supplement_rows: list[dict[str, str]] = []
    method_counts: Counter[str] = Counter()
    total_fields = 0

    def make_row(field_id: str, field_title: str, mapped_id: str, mapped_label: str, method: str) -> dict[str, str]:
        clean_title = normalize(field_title) or f"UKB field {field_id}"
        return {
            "UKBB data field": field_id,
            "ICD10": "-",
            "PheCode": "-",
            "Lookup text": clean_title,
            "Curated reported trait": f"{clean_title} (UKB data field {field_id})",
            "main ontology URI(s)": mapped_id,
            "main ontology label(s)": mapped_label,
            "Publication": UKB_FIELD_SUPPLEMENT_PUBLICATION,
            "User submitted trait": "-",
            "qc_status": "review_required",
            "qc_resolution_modes": method,
            "qc_label_match_modes": method,
            "qc_unresolved_ids": "",
            "final_update_applied": "yes",
            "final_update_reason": "ukb_field_supplement_cache_generated",
        }

    sortable_fields = sorted(
        ukb_field_titles.items(),
        key=lambda item: (int(item[0]) if item[0].isdigit() else 10**9, item[0]),
    )
    for field_id, field_title in sortable_fields:
        if not field_id:
            continue
        total_fields += 1
        if field_id in base_covered_fields:
            continue

        chosen_mapping: tuple[str, str, str] | None = None
        field_exact_keys: list[str] = []
        for variant in ukb_title_query_variants(field_id, field_title):
            key = norm_key(variant)
            if key and key not in field_exact_keys:
                field_exact_keys.append(key)

        matched_record_idxs: list[int] = []
        seen_rec_idxs: set[int] = set()
        for key in field_exact_keys:
            for rec_idx in text_exact_index.get(key, []):
                record = records[rec_idx]
                if record.is_ukb_field_supplement:
                    continue
                if rec_idx in seen_rec_idxs:
                    continue
                seen_rec_idxs.add(rec_idx)
                matched_record_idxs.append(rec_idx)
        if matched_record_idxs:
            matched_records = [records[idx] for idx in matched_record_idxs]
            chosen_mapping = choose_mapping_from_records(
                matched_records,
                method="ukb_field_supplement_field_title_cache_exact",
                matched_on="UKB field title",
            )
            if (
                chosen_mapping is not None
                and is_generic_ukb_measurement_mapping(chosen_mapping[0], chosen_mapping[1])
            ):
                chosen_mapping = None

        if chosen_mapping is None:
            exact_term_ids: set[str] = set()
            for key in field_exact_keys:
                exact_term_ids |= ontology_exact_index.get(key, set())
            selected_term = pick_preferred_ontology_term(exact_term_ids, ontology_terms=ontology_terms)
            if selected_term is not None:
                term_id, term_label = selected_term
                chosen_mapping = (
                    term_id,
                    term_label,
                    "ukb_field_supplement_field_title_ontology_exact",
                )
                if is_generic_ukb_measurement_mapping(chosen_mapping[0], chosen_mapping[1]):
                    chosen_mapping = None

        if chosen_mapping is None:
            chosen_mapping = choose_mapping_from_ontology_fuzzy(
                field_title,
                method="ukb_field_supplement_field_title_ontology_fuzzy",
                min_score=0.95,
                min_margin=0.02,
            )

        if chosen_mapping is None:
            category_id = normalize(ukb_field_to_category.get(field_id, ""))
            chosen_mapping = category_default_mapping.get(category_id)

        if chosen_mapping is None:
            chosen_mapping = (
                UKB_FIELD_SUPPLEMENT_GENERIC_MEASUREMENT_ID,
                UKB_FIELD_SUPPLEMENT_GENERIC_MEASUREMENT_LABEL,
                "ukb_field_supplement_generic_measurement",
            )

        mapped_id, mapped_label, method = chosen_mapping
        method_counts[method] += 1
        supplement_rows.append(make_row(field_id, field_title, mapped_id, mapped_label, method))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TRAIT_CACHE_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(supplement_rows)

    return {
        "total_fields": total_fields,
        "base_covered_fields": len(base_covered_fields & set(ukb_field_titles)),
        "supplement_rows": len(supplement_rows),
        "method_category_cache_unique": method_counts.get("ukb_field_supplement_category_cache_unique", 0),
        "method_category_text_exact": method_counts.get("ukb_field_supplement_category_text_exact", 0),
        "method_category_ontology_exact": method_counts.get("ukb_field_supplement_category_ontology_exact", 0),
        "method_field_title_cache_exact": method_counts.get("ukb_field_supplement_field_title_cache_exact", 0),
        "method_field_title_ontology_exact": method_counts.get("ukb_field_supplement_field_title_ontology_exact", 0),
        "method_field_title_ontology_fuzzy": method_counts.get("ukb_field_supplement_field_title_ontology_fuzzy", 0),
        "method_generic_disease": method_counts.get("ukb_field_supplement_generic_disease", 0),
        "method_generic_measurement": method_counts.get("ukb_field_supplement_generic_measurement", 0),
    }


def confidence_bucket(score: float) -> str:
    if score >= 0.93:
        return "high"
    if score >= 0.82:
        return "medium"
    return "low"


def make_trait_provenance(
    *,
    method: str,
    input_type: str,
    trait_scale: str,
    input_value: str,
    matched_on: str,
    mapped_ids: str,
    mapped_labels: str,
    source_file: str,
    score: float,
    validated: bool,
) -> str:
    confidence = confidence_bucket(score)
    needs_review = "no" if validated else "yes"
    normalized_input = norm_key(input_value)
    return (
        "mapping_provenance=v1;"
        f"method={method};"
        f"input_type={input_type};"
        f"trait_scale={trait_scale};"
        f"input_value={input_value};"
        f"normalized_input={normalized_input};"
        f"matched_on={matched_on};"
        f"mapped_ids={mapped_ids};"
        f"mapped_labels={mapped_labels};"
        f"source_file={source_file};"
        f"confidence={confidence};"
        f"needs_review={needs_review}"
    )


def trait_candidate_prefix_rank(candidate_id: str) -> int:
    parts = split_multi_ids(candidate_id)
    if not parts:
        return 999
    best = 999
    for part in parts:
        rank = TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(part), 999)
        if rank < best:
            best = rank
    return best


def mondo_import_qc(
    mapped_id_text: str,
    *,
    ontology_terms: dict[str, TraitOntologyTerm],
) -> tuple[str, str]:
    mapped_ids = [canonicalize_trait_ontology_id(part) for part in split_multi_ids(mapped_id_text)]
    mapped_ids = [item for item in mapped_ids if item]
    if not mapped_ids:
        return "", ""

    missing_mondo_ids = [
        item
        for item in mapped_ids
        if trait_ontology_prefix(item) == "MONDO" and item not in ontology_terms
    ]
    if not missing_mondo_ids:
        return "", ""

    has_active_non_mondo = any(
        item in ontology_terms and trait_ontology_prefix(item) in {"EFO", "OBA", "HP", "NCIT"}
        for item in mapped_ids
    )
    if has_active_non_mondo:
        return "", ""

    return "term not in EFO", "|".join(missing_mondo_ids)


def trait_term_branch(
    *,
    term_id: str,
    label: str,
    ontology_terms: dict[str, TraitOntologyTerm],
    measurement_term_ids: set[str],
) -> str:
    oid = canonicalize_trait_ontology_id(term_id)
    term = ontology_terms.get(oid)
    text = norm_key(label or (term.label if term is not None else ""))
    toks = tokenize(text)
    prefix = trait_ontology_prefix(oid)

    if oid in measurement_term_ids:
        return "measurement"
    if prefix == "OBA":
        return "measurement"
    if toks & MEASUREMENT_HINT_TOKENS:
        return "measurement"
    if has_disease_hint(text):
        return "non_measurement"
    if prefix in {"MONDO", "HP"}:
        return "non_measurement"
    if text:
        # Default unresolved trait labels to non-measurement so quantitative
        # inputs only auto-validate when measurement evidence is explicit.
        return "non_measurement"
    return "unknown"


def trait_candidate_branch(
    *,
    mapped_id_text: str,
    mapped_label_text: str,
    ontology_terms: dict[str, TraitOntologyTerm],
    measurement_term_ids: set[str],
) -> str:
    mapped_ids = [canonicalize_trait_ontology_id(part) for part in split_multi_ids(mapped_id_text)]
    mapped_ids = [item for item in mapped_ids if item]
    if not mapped_ids:
        return "unknown"

    labels = split_multi_labels(mapped_label_text, expected_n=len(mapped_ids))
    measurement_count = 0
    non_measurement_count = 0
    unknown_count = 0
    for idx, term_id in enumerate(mapped_ids):
        label = normalize(labels[idx] if idx < len(labels) else "")
        branch = trait_term_branch(
            term_id=term_id,
            label=label,
            ontology_terms=ontology_terms,
            measurement_term_ids=measurement_term_ids,
        )
        if branch == "measurement":
            measurement_count += 1
        elif branch == "non_measurement":
            non_measurement_count += 1
        else:
            unknown_count += 1

    if measurement_count and non_measurement_count:
        return "mixed"
    if measurement_count:
        return "measurement"
    if non_measurement_count:
        return "non_measurement"
    if unknown_count:
        return "unknown"
    return "unknown"


def trait_scale_mismatch(trait_scale: str, candidate_branch: str) -> bool:
    scale = normalize_trait_scale(trait_scale)
    if scale == "auto":
        return False
    if scale == "quantitative":
        return candidate_branch in {"non_measurement", "mixed"}
    if scale == "binary":
        return candidate_branch in {"measurement", "mixed"}
    return False


def trait_scale_score_adjustment(trait_scale: str, candidate_branch: str) -> float:
    scale = normalize_trait_scale(trait_scale)
    if scale == "auto" or candidate_branch == "unknown":
        return 0.0
    if trait_scale_mismatch(scale, candidate_branch):
        return -0.20
    return 0.06


TRAIT_OUTPUT_FIELDS = [
    "input_row_id",
    "input_query",
    "negation",
    "input_type",
    "input_trait_scale",
    "input_source_code",
    "input_source_data_type",
    "input_icd10",
    "input_icd10_label",
    "input_phecode",
    "input_ukb_field_id",
    "input_ukb_field_label",
    "input_ukb_category_id",
    "input_ukb_category_label",
    "input_ukb_category_path",
    "mapped_trait_id",
    "mapped_trait_label",
    "confidence",
    "matched_via",
    "matched_on",
    "source_file",
    "validation",
    "qc_review_flag",
    "specificity_loss_flag",
    "specificity_loss_notes",
    "term_not_in_efo",
    "mondo_missing_ids",
    "evidence",
    "provenance",
]


def map_trait_queries(
    query_inputs: list[dict[str, str]],
    *,
    trait_cache_index: dict[str, Any],
    ontology_index: dict[str, Any],
    trait_cache_resolution: dict[int, tuple[tuple[str, ...], tuple[str, ...]]] | None,
    trait_cache_path: Path,
    efo_obo_path: Path,
    ukb_field_titles: dict[str, str] | None,
    ukb_category_titles: dict[str, str] | None,
    ukb_category_parents: dict[str, tuple[str, ...]] | None,
    ukb_category_to_fields: dict[str, tuple[str, ...]] | None,
    ukb_field_to_category: dict[str, str] | None,
    top_k: int,
    min_score: float,
    force_map_best: bool,
    show_progress: bool,
    row_callback: Any | None = None,
    collect_rows: bool = True,
    memoize_queries: bool = True,
    query_cache_max_entries: int = 200_000,
) -> list[dict[str, str]]:
    records: list[TraitCacheRecord] = trait_cache_index["records"]
    icd10_index: dict[str, list[int]] = trait_cache_index["icd10_index"]
    icd10_label_index: dict[str, tuple[str, ...]] = trait_cache_index.get("icd10_label_index", {})
    official_icd10_label_index: dict[str, tuple[str, ...]] = trait_cache_index.get("official_icd10_label_index", {})
    phecode_index: dict[str, list[int]] = trait_cache_index["phecode_index"]
    ukb_field_index: dict[str, list[int]] = trait_cache_index.get("ukb_field_index", {})
    text_exact_index: dict[str, list[int]] = trait_cache_index["text_exact_index"]
    cache_token_index: dict[str, set[int]] = trait_cache_index["token_index"]
    cache_token_freq: Counter[str] = trait_cache_index["token_freq"]

    ontology_terms: dict[str, TraitOntologyTerm] = ontology_index["terms"]
    measurement_term_ids: set[str] = set(ontology_index.get("measurement_term_ids", set()))
    ontology_exact_index: dict[str, set[str]] = ontology_index["exact_index"]
    ontology_token_index: dict[str, set[str]] = ontology_index["token_index"]
    ontology_token_freq: Counter[str] = ontology_index["token_freq"]
    ontology_icd10_xref_index: dict[str, set[str]] = ontology_index.get("icd10_xref_index", {})
    ukb_field_titles = ukb_field_titles or {}
    ukb_category_titles = ukb_category_titles or {}
    ukb_category_parents = ukb_category_parents or {}
    ukb_category_to_fields = ukb_category_to_fields or {}
    ukb_field_to_category = ukb_field_to_category or {}
    ukb_category_title_norms = {
        category_id: norm_key(title)
        for category_id, title in ukb_category_titles.items()
        if " " in norm_key(title) and len(norm_key(title)) >= 12
    }
    cache_resolution = trait_cache_resolution or {}
    category_mapping_only_after_other_methods = True

    def resolve_record(record: TraitCacheRecord) -> tuple[tuple[str, ...], tuple[str, ...]] | None:
        if cache_resolution:
            return cache_resolution.get(record.rownum)
        if not record.mapped_ids:
            return None
        return record.mapped_ids, record.mapped_labels

    rows: list[dict[str, str]] = []
    query_result_cache: dict[tuple[str, str, str, str], tuple[dict[str, str], ...]] = {}
    ukb_category_fallback_cache: dict[tuple[str, ...], tuple[Candidate, ...]] = {}
    ukb_field_fallback_cache: dict[str, tuple[Candidate, ...]] = {}
    progress = ProgressReporter(total=len(query_inputs), enabled=show_progress)

    via_rank = {
        "cache_exact_icd10": 0,
        "cache_exact_icd10_supplement": 0,
        "cache_parent_icd10": 1,
        "cache_parent_icd10_supplement": 1,
        "cache_exact_phecode": 1,
        "efo_obo_exact_icd10_xref": 1,
        "efo_obo_parent_icd10_xref": 2,
        "cache_exact_text": 2,
        "cache_exact_ukb_field": 2,
        "cache_exact_ukb_field_supplement": 2,
        "cache_exact_ukb_field_title": 2,
        "cache_exact_ukb_category": 2,
        "cache_exact_ukb_category_title": 2,
        "cache_exact_ukb_category_field_title": 2,
        "cache_fuzzy_text": 3,
        "efo_obo_exact": 4,
        "efo_obo_exact_ukb_field_title": 4,
        "efo_obo_exact_ukb_category": 4,
        "efo_obo_fuzzy": 5,
    }

    def build_ukb_category_fallback_candidates(fallback_key: tuple[str, ...]) -> tuple[Candidate, ...]:
        if not fallback_key:
            return ()
        cached_fallback = ukb_category_fallback_cache.get(fallback_key)
        if cached_fallback is not None:
            return cached_fallback

        fallback_candidates: list[Candidate] = []
        fallback_text_keys: list[str] = []
        fallback_category_field_ids: set[str] = set()

        for category_id in fallback_key:
            category_title = normalize(ukb_category_titles.get(category_id, ""))
            path_titles = ukb_category_path_titles(
                category_id,
                category_titles=ukb_category_titles,
                category_parents=ukb_category_parents,
            )
            for variant in ukb_category_query_variants(category_id, category_title, path_titles):
                key = norm_key(variant)
                if key and key not in fallback_text_keys:
                    fallback_text_keys.append(key)
            # Include direct parent-path titles (for example "physical activity measurement")
            # so technical UKB child fields can still map at the broader phenotype level.
            for title in path_titles:
                key = norm_key(title)
                if key and key not in fallback_text_keys:
                    fallback_text_keys.append(key)
            for field_id in ukb_category_to_fields.get(category_id, ()):
                field_id_norm = normalize(field_id)
                if field_id_norm:
                    fallback_category_field_ids.add(field_id_norm)

        if fallback_category_field_ids:
            matched_record_idxs: list[int] = []
            seen_rec_idxs: set[int] = set()
            for field_id in sorted(fallback_category_field_ids):
                for rec_idx in ukb_field_index.get(field_id, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                ukb_category_candidates = collapse_trait_cache_candidates(
                    matched_records,
                    method="cache_exact_ukb_category",
                    score=0.962,
                    matched_on="UKB category field set (fallback)",
                    validated_if_unique=True,
                    record_resolver=resolve_record,
                )
                if len(ukb_category_candidates) == 1:
                    fallback_candidates.extend(ukb_category_candidates)

        if not fallback_candidates and fallback_text_keys:
            matched_record_idxs = []
            seen_rec_idxs: set[int] = set()
            for key in fallback_text_keys:
                for rec_idx in text_exact_index.get(key, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                ukb_category_title_candidates = collapse_trait_cache_candidates(
                    matched_records,
                    method="cache_exact_ukb_category_title",
                    score=0.955,
                    matched_on="UKB category title/path text (fallback)",
                    validated_if_unique=False,
                    record_resolver=resolve_record,
                )
                if len(ukb_category_title_candidates) == 1:
                    fallback_candidates.extend(ukb_category_title_candidates)

        if not fallback_candidates and fallback_text_keys:
            exact_term_ids: set[str] = set()
            for key in fallback_text_keys:
                exact_term_ids |= ontology_exact_index.get(key, set())
            if exact_term_ids:
                ordered = sorted(
                    exact_term_ids,
                    key=lambda tid: (
                        TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(tid), 999),
                        tid,
                    ),
                )
                for term_id in ordered[: max(top_k, 5)]:
                    term = ontology_terms[term_id]
                    fallback_candidates.append(
                        Candidate(
                            efo_id=term.term_id,
                            label=term.label,
                            score=0.90 if len(ordered) == 1 else 0.86,
                            matched_via="efo_obo_exact_ukb_category",
                            evidence="exact label/synonym match in efo.obo (via UKB category fallback)",
                            is_validated=False,
                        )
                    )

        cached = tuple(fallback_candidates)
        ukb_category_fallback_cache[fallback_key] = cached
        return cached

    should_precompute_ukb_field_fallback = bool(ukb_field_to_category) and any(
        (
            extract_ukb_field_ids(normalize(item.get("query", "")))
            or extract_ukb_category_ids(normalize(item.get("query", "")))
            or re.fullmatch(r"(?:f)?[0-9]{1,7}", norm_key(normalize(item.get("query", "")))) is not None
        )
        for item in query_inputs
    )
    if should_precompute_ukb_field_fallback:
        precompute_category_ids = sorted(
            {
                category_id
                for category_id in (normalize(raw) for raw in ukb_field_to_category.values())
                if category_id
            }
        )
        for category_id in precompute_category_ids:
            build_ukb_category_fallback_candidates((category_id,))
        for field_id, category_id in ukb_field_to_category.items():
            fid = normalize(field_id)
            cid = normalize(category_id)
            if not fid or not cid:
                continue
            precomputed = ukb_category_fallback_cache.get((cid,), ())
            if precomputed:
                ukb_field_fallback_cache[fid] = precomputed

    for item in query_inputs:
        input_row_id = normalize(item.get("row_id", ""))
        query = normalize(item.get("query", ""))
        input_type = normalize_trait_input_type(item.get("input_type", "auto"))
        trait_scale = normalize_trait_scale(item.get("trait_scale", "auto"))
        source_code = normalize(item.get("source_code", ""))
        source_data_type = normalize(item.get("source_data_type", ""))
        icd10 = normalize_icd10_code(item.get("icd10", ""))
        icd10_label = ""
        icd10_match_label = ""
        phecode = normalize_phecode(item.get("phecode", ""))
        parsed_query_icd10, parsed_query_icd10_label = parse_icd10_labeled_trait_text(query)
        if parsed_query_icd10_label and (not icd10 or parsed_query_icd10 == icd10):
            icd10_label = parsed_query_icd10_label
            icd10_match_label = parsed_query_icd10_label

        if icd10 and not is_strict_icd10_code(icd10):
            icd10 = ""
        if input_type == "icd10" and not icd10:
            source_code_icd10 = normalize_icd10_code(source_code)
            if is_strict_icd10_code(source_code_icd10):
                icd10 = source_code_icd10
        if input_type == "icd10" and not icd10 and is_strict_icd10_code(parsed_query_icd10):
            icd10 = parsed_query_icd10
        if input_type == "icd10" and not icd10:
            query_icd10 = normalize_icd10_code(query)
            if is_strict_icd10_code(query_icd10):
                icd10 = query_icd10
        if icd10 and not icd10_label and parsed_query_icd10_label and parsed_query_icd10 == icd10:
            icd10_label = parsed_query_icd10_label
            icd10_match_label = parsed_query_icd10_label
        if icd10 and not icd10_label:
            official_icd10_labels = tuple(label for label in official_icd10_label_index.get(icd10, ()) if label)
            if official_icd10_labels:
                icd10_label = "|".join(official_icd10_labels)
                if len(official_icd10_labels) == 1:
                    icd10_match_label = official_icd10_labels[0]
        if icd10 and not icd10_label:
            cached_icd10_labels = tuple(label for label in icd10_label_index.get(icd10, ()) if label)
            if cached_icd10_labels:
                icd10_label = "|".join(cached_icd10_labels)
                if len(cached_icd10_labels) == 1:
                    icd10_match_label = cached_icd10_labels[0]
        if input_type == "phecode" and not phecode:
            phecode = normalize_phecode(query)

        query_norm_hint = norm_key(query)
        ukb_field_ids = extract_ukb_field_ids(query)
        if not ukb_field_ids:
            numeric_field_match = re.fullmatch(r"(?:f)?([0-9]{1,7})", query_norm_hint)
            if numeric_field_match:
                maybe_field_id = normalize(numeric_field_match.group(1))
                if maybe_field_id and (maybe_field_id in ukb_field_titles or maybe_field_id in ukb_field_to_category):
                    ukb_field_ids = (maybe_field_id,)

        if input_type == "auto":
            if icd10:
                input_type = "icd10"
            elif extract_ukb_field_ids(query):
                input_type = "ukb_field"
            elif phecode or looks_like_phecode_token(query):
                input_type = "phecode"
                if not phecode:
                    phecode = normalize_phecode(query)
            else:
                input_type = "trait_text"
        elif input_type == "ukb_field" and not ukb_field_ids and query:
            numeric_field_match = re.fullmatch(r"(?:f)?([0-9]{1,7})", query_norm_hint)
            if numeric_field_match:
                ukb_field_ids = (normalize(numeric_field_match.group(1)),)
        if not query:
            query = icd10 or phecode

        query_for_matching = query
        prefer_icd10_from_ukb = False
        parsed_ukb_field_id = ""
        if input_type == "ukb_field" and query:
            parsed_field_id, parsed_code, parsed_label = parse_ukb_hashed_trait_text(query)
            parsed_ukb_field_id = parsed_field_id
            if parsed_field_id and parsed_field_id not in ukb_field_ids:
                ukb_field_ids = (*ukb_field_ids, parsed_field_id)
            if parsed_label:
                query_for_matching = parsed_label
            elif parsed_code:
                query_for_matching = parsed_code

            # For UKB rows that carry explicit ICD10 values inside coded text,
            # route through the ICD10 exact pathway when the field metadata
            # explicitly indicates ICD10 coding.
            if parsed_code and parsed_field_id:
                field_title_key = norm_key(ukb_field_titles.get(parsed_field_id, ""))
                if "icd10" in field_title_key:
                    prefer_icd10_from_ukb = True
                    if not icd10:
                        code_icd10 = normalize_icd10_code(parsed_code)
                        if is_strict_icd10_code(code_icd10):
                            icd10 = code_icd10

        if not icd10 and query:
            allow_query_icd10_inference = (
                input_type != "ukb_field" or not parsed_ukb_field_id or prefer_icd10_from_ukb
            )
            if allow_query_icd10_inference:
                inferred_icd10 = normalize_icd10_code(query)
                if is_strict_icd10_code(inferred_icd10):
                    icd10 = inferred_icd10
        if input_type == "icd10" and icd10_match_label:
            query_for_matching = icd10_match_label

        if not query_for_matching:
            query_for_matching = query

        query_norm_hint = norm_key(query)
        ukb_category_ids = extract_ukb_category_ids(query)
        if not ukb_category_ids and query_norm_hint and ukb_category_title_norms:
            inferred_category_ids = tuple(
                category_id
                for category_id, title_norm in ukb_category_title_norms.items()
                if title_norm and title_norm in query_norm_hint
            )
            if inferred_category_ids:
                ukb_category_ids = inferred_category_ids

        resolved_ukb_category_ids: list[str] = []
        seen_category_ids: set[str] = set()
        for category_id in ukb_category_ids:
            cid = normalize(category_id)
            if cid and cid not in seen_category_ids:
                seen_category_ids.add(cid)
                resolved_ukb_category_ids.append(cid)
        for field_id in ukb_field_ids:
            category_id = normalize(ukb_field_to_category.get(field_id, ""))
            if category_id and category_id not in seen_category_ids:
                seen_category_ids.add(category_id)
                resolved_ukb_category_ids.append(category_id)

        ukb_exact_text_keys: list[str] = []
        ukb_primary_title = ""
        for ukb_field_id in ukb_field_ids:
            title = normalize(ukb_field_titles.get(ukb_field_id, ""))
            if not title:
                continue
            if not ukb_primary_title:
                ukb_primary_title = title
            for variant in ukb_title_query_variants(ukb_field_id, title):
                key = norm_key(variant)
                if key and key not in ukb_exact_text_keys:
                    ukb_exact_text_keys.append(key)

        ukb_category_exact_text_keys: list[str] = []
        ukb_category_field_exact_text_keys: list[str] = []
        ukb_category_field_ids: set[str] = set()
        for category_id in ukb_category_ids:
            category_title = normalize(ukb_category_titles.get(category_id, ""))
            path_titles = ukb_category_path_titles(
                category_id,
                category_titles=ukb_category_titles,
                category_parents=ukb_category_parents,
            )
            for variant in ukb_category_query_variants(category_id, category_title, path_titles):
                key = norm_key(variant)
                if key and key not in ukb_category_exact_text_keys:
                    ukb_category_exact_text_keys.append(key)
            for field_id in ukb_category_to_fields.get(category_id, ()):
                field_id_norm = normalize(field_id)
                if not field_id_norm:
                    continue
                ukb_category_field_ids.add(field_id_norm)
                field_title = normalize(ukb_field_titles.get(field_id_norm, ""))
                if not field_title:
                    continue
                for variant in ukb_title_query_variants(field_id_norm, field_title):
                    key = norm_key(variant)
                    if key and key not in ukb_category_field_exact_text_keys:
                        ukb_category_field_exact_text_keys.append(key)

        ukb_field_label_values = [
            normalize(ukb_field_titles.get(field_id, ""))
            for field_id in ukb_field_ids
            if normalize(ukb_field_titles.get(field_id, ""))
        ]
        ukb_category_label_values = [
            normalize(ukb_category_titles.get(category_id, ""))
            for category_id in resolved_ukb_category_ids
            if normalize(ukb_category_titles.get(category_id, ""))
        ]
        ukb_category_path_values = [
            " > ".join(
                title
                for title in ukb_category_path_titles(
                    category_id,
                    category_titles=ukb_category_titles,
                    category_parents=ukb_category_parents,
                )
                if title
            )
            for category_id in resolved_ukb_category_ids
        ]
        ukb_category_path_values = [value for value in ukb_category_path_values if value]
        ukb_field_id_text = "|".join(ukb_field_ids)
        ukb_field_label_text = "|".join(ukb_field_label_values)
        ukb_category_id_text = "|".join(resolved_ukb_category_ids)
        ukb_category_label_text = "|".join(ukb_category_label_values)
        ukb_category_path_text = "|".join(ukb_category_path_values)

        cache_key = (query, input_type, trait_scale, icd10, phecode)
        if memoize_queries and cache_key in query_result_cache:
            cached_rows = query_result_cache[cache_key]
            for cached in cached_rows:
                row = dict(cached)
                row["input_row_id"] = input_row_id
                row["input_query"] = query
                row["input_trait_scale"] = trait_scale
                row["input_source_code"] = source_code
                row["input_source_data_type"] = source_data_type
                if row_callback is not None:
                    row_callback(row)
                if collect_rows:
                    rows.append(row)
            progress.update()
            continue

        query_exact_norms = trait_query_exact_norm_variants(query_for_matching)
        query_exact_key_set = {key for key in query_exact_norms if key}
        query_norm = query_exact_norms[0] if query_exact_norms else ""
        query_match_norm = normalize_trait_text_for_similarity(query_for_matching)
        query_tokens = informative_trait_tokens(query_match_norm or query_norm)
        query_has_negation = trait_has_negation(query) or trait_has_negation(query_for_matching)
        has_ontology_exact = any(bool(ontology_exact_index.get(key)) for key in query_exact_norms)
        if not has_ontology_exact and ukb_exact_text_keys and not prefer_icd10_from_ukb:
            has_ontology_exact = any(bool(ontology_exact_index.get(key)) for key in ukb_exact_text_keys)

        fuzzy_query_norm = query_match_norm
        fuzzy_query_tokens = query_tokens

        candidates: list[Candidate] = []
        best_any: list[Candidate] = []
        generated_rows: list[dict[str, str]] = []

        if input_type == "trait_text" and trait_is_noninformative_negation(query_for_matching):
            generated_rows.append(
                {
                    "input_row_id": input_row_id,
                    "input_query": query,
                    "negation": "yes",
                    "input_type": input_type,
                    "input_trait_scale": trait_scale,
                    "input_source_code": source_code,
                    "input_source_data_type": source_data_type,
                    "input_icd10": icd10,
                    "input_icd10_label": icd10_label,
                    "input_phecode": phecode,
                    "input_ukb_field_id": ukb_field_id_text,
                    "input_ukb_field_label": ukb_field_label_text,
                    "input_ukb_category_id": ukb_category_id_text,
                    "input_ukb_category_label": ukb_category_label_text,
                    "input_ukb_category_path": ukb_category_path_text,
                    "mapped_trait_id": "",
                    "mapped_trait_label": "",
                    "confidence": "0.000",
                    "matched_via": "none",
                    "matched_on": "",
                    "source_file": "",
                    "validation": "not_mapped",
                    "qc_review_flag": "yes",
                    "specificity_loss_flag": "no",
                    "specificity_loss_notes": "",
                    "term_not_in_efo": "",
                    "mondo_missing_ids": "",
                    "evidence": "non-informative negation trait option",
                    "provenance": "",
                }
            )
            if memoize_queries and len(query_result_cache) < max(0, query_cache_max_entries):
                query_result_cache[cache_key] = tuple(dict(row) for row in generated_rows)
            for row in generated_rows:
                if row_callback is not None:
                    row_callback(dict(row))
                if collect_rows:
                    rows.append(dict(row))
            progress.update()
            continue

        if icd10 and icd10 in icd10_index:
            matched_records = [records[idx] for idx in icd10_index[icd10]]
            curated_records = [record for record in matched_records if not record.is_icd10_supplement]
            supplemental_records = [record for record in matched_records if record.is_icd10_supplement]
            if curated_records:
                candidates.extend(
                    collapse_trait_cache_candidates(
                        curated_records,
                        method="cache_exact_icd10",
                        score=1.0,
                        matched_on="ICD10",
                        validated_if_unique=True,
                        record_resolver=resolve_record,
                    )
                )
            elif supplemental_records:
                icd10_supplement_candidates = collapse_trait_cache_candidates(
                    supplemental_records,
                    method="cache_exact_icd10_supplement",
                    score=0.93,
                    matched_on="ICD10 supplemental cache",
                    validated_if_unique=False,
                    record_resolver=resolve_record,
                )
                if len(icd10_supplement_candidates) == 1:
                    candidates.extend(icd10_supplement_candidates)

        if not candidates and icd10:
            for parent_icd10 in icd10_parent_fallback_codes(icd10):
                if parent_icd10 not in icd10_index:
                    continue
                matched_records = [records[idx] for idx in icd10_index[parent_icd10]]
                curated_records = [record for record in matched_records if not record.is_icd10_supplement]
                supplemental_records = [record for record in matched_records if record.is_icd10_supplement]
                if curated_records:
                    parent_candidates = collapse_trait_cache_candidates(
                        curated_records,
                        method="cache_parent_icd10",
                        score=0.90,
                        matched_on=f"ICD10 parent fallback ({parent_icd10})",
                        validated_if_unique=False,
                        record_resolver=resolve_record,
                    )
                    if parent_candidates:
                        candidates.extend(parent_candidates)
                        break
                if supplemental_records:
                    parent_supplement_candidates = collapse_trait_cache_candidates(
                        supplemental_records,
                        method="cache_parent_icd10_supplement",
                        score=0.86,
                        matched_on=f"ICD10 parent supplemental fallback ({parent_icd10})",
                        validated_if_unique=False,
                        record_resolver=resolve_record,
                    )
                    if len(parent_supplement_candidates) == 1:
                        candidates.extend(parent_supplement_candidates)
                        break

        if not candidates and icd10:
            xref_term_ids = ontology_icd10_xref_index.get(icd10, set())
            if xref_term_ids:
                ordered = sorted(
                    (term_id for term_id in xref_term_ids if term_id in ontology_terms),
                    key=lambda tid: (
                        TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(tid), 999),
                        tid,
                    ),
                )
                if ordered:
                    for term_id in ordered[: max(top_k, 5)]:
                        term = ontology_terms[term_id]
                        candidates.append(
                            Candidate(
                                efo_id=term.term_id,
                                label=term.label,
                                score=0.94 if len(ordered) == 1 else 0.88,
                                matched_via="efo_obo_exact_icd10_xref",
                                evidence="exact ICD10 xref match in efo.obo",
                                is_validated=len(ordered) == 1,
                            )
                        )

        if not candidates and icd10:
            for parent_icd10 in icd10_parent_fallback_codes(icd10):
                xref_term_ids = ontology_icd10_xref_index.get(parent_icd10, set())
                if not xref_term_ids:
                    continue
                ordered = sorted(
                    (term_id for term_id in xref_term_ids if term_id in ontology_terms),
                    key=lambda tid: (
                        TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(tid), 999),
                        tid,
                    ),
                )
                if not ordered:
                    continue
                for term_id in ordered[: max(top_k, 5)]:
                    term = ontology_terms[term_id]
                    candidates.append(
                        Candidate(
                            efo_id=term.term_id,
                            label=term.label,
                            score=0.86 if len(ordered) == 1 else 0.82,
                            matched_via="efo_obo_parent_icd10_xref",
                            evidence=f"ICD10 parent fallback xref match in efo.obo ({parent_icd10})",
                            is_validated=False,
                        )
                    )
                if candidates:
                    break

        if phecode and phecode in phecode_index:
            matched_records = [records[idx] for idx in phecode_index[phecode]]
            candidates.extend(
                collapse_trait_cache_candidates(
                    matched_records,
                    method="cache_exact_phecode",
                    score=0.99,
                    matched_on="PheCode",
                    validated_if_unique=True,
                    record_resolver=resolve_record,
                )
            )

        if ukb_field_ids and not prefer_icd10_from_ukb:
            matched_record_idxs: list[int] = []
            seen_rec_idxs: set[int] = set()
            for ukb_field_id in ukb_field_ids:
                for rec_idx in ukb_field_index.get(ukb_field_id, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                curated_records = [record for record in matched_records if not record.is_ukb_field_supplement]
                supplemental_records = [record for record in matched_records if record.is_ukb_field_supplement]

                if curated_records:
                    ukb_field_candidates = collapse_trait_cache_candidates(
                        curated_records,
                        method="cache_exact_ukb_field",
                        score=0.985,
                        matched_on="UKB data field code",
                        validated_if_unique=True,
                        record_resolver=resolve_record,
                    )
                    # Field IDs may fan out to multiple ontology mappings in cache;
                    # only use this fast path when the mapping is unambiguous.
                    if len(ukb_field_candidates) == 1:
                        candidates.extend(ukb_field_candidates)
                elif supplemental_records:
                    supplemental_candidates = collapse_trait_cache_candidates(
                        supplemental_records,
                        method="cache_exact_ukb_field_supplement",
                        score=0.90,
                        matched_on="UKB data field supplemental cache",
                        validated_if_unique=False,
                        record_resolver=resolve_record,
                    )
                    if len(supplemental_candidates) == 1:
                        candidates.extend(supplemental_candidates)

        if (not category_mapping_only_after_other_methods) and ukb_category_field_ids and not prefer_icd10_from_ukb:
            matched_record_idxs = []
            seen_rec_idxs: set[int] = set()
            for ukb_field_id in sorted(ukb_category_field_ids):
                for rec_idx in ukb_field_index.get(ukb_field_id, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                ukb_category_candidates = collapse_trait_cache_candidates(
                    matched_records,
                    method="cache_exact_ukb_category",
                    score=0.982,
                    matched_on="UKB category field set",
                    validated_if_unique=True,
                    record_resolver=resolve_record,
                )
                # Category IDs can include multiple child fields with different
                # trait mappings; only auto-accept the unambiguous case.
                if len(ukb_category_candidates) == 1:
                    candidates.extend(ukb_category_candidates)

        if query_exact_norms:
            matched_record_idxs: list[int] = []
            seen_rec_idxs: set[int] = set()
            for key in query_exact_norms:
                for rec_idx in text_exact_index.get(key, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                candidates.extend(
                    collapse_trait_cache_candidates(
                        matched_records,
                        method="cache_exact_text",
                        score=0.97,
                        matched_on="Reported trait text",
                        validated_if_unique=True,
                        record_resolver=resolve_record,
                    )
                )

        if ukb_exact_text_keys and not prefer_icd10_from_ukb:
            matched_record_idxs: list[int] = []
            seen_rec_idxs: set[int] = set()
            for key in ukb_exact_text_keys:
                for rec_idx in text_exact_index.get(key, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                ukb_title_candidates = collapse_trait_cache_candidates(
                    matched_records,
                    method="cache_exact_ukb_field_title",
                    score=0.975,
                    matched_on="UKB field title text",
                    validated_if_unique=True,
                    record_resolver=resolve_record,
                )
                if len(ukb_title_candidates) == 1:
                    candidates.extend(ukb_title_candidates)

        if (not category_mapping_only_after_other_methods) and ukb_category_exact_text_keys and not prefer_icd10_from_ukb:
            matched_record_idxs = []
            seen_rec_idxs: set[int] = set()
            for key in ukb_category_exact_text_keys:
                for rec_idx in text_exact_index.get(key, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                ukb_category_title_candidates = collapse_trait_cache_candidates(
                    matched_records,
                    method="cache_exact_ukb_category_title",
                    score=0.972,
                    matched_on="UKB category title text",
                    validated_if_unique=True,
                    record_resolver=resolve_record,
                )
                if len(ukb_category_title_candidates) == 1:
                    candidates.extend(ukb_category_title_candidates)

        if (not category_mapping_only_after_other_methods) and ukb_category_field_exact_text_keys and not prefer_icd10_from_ukb:
            matched_record_idxs = []
            seen_rec_idxs: set[int] = set()
            for key in ukb_category_field_exact_text_keys:
                for rec_idx in text_exact_index.get(key, []):
                    if rec_idx in seen_rec_idxs:
                        continue
                    seen_rec_idxs.add(rec_idx)
                    matched_record_idxs.append(rec_idx)
            if matched_record_idxs:
                matched_records = [records[idx] for idx in matched_record_idxs]
                ukb_category_field_title_candidates = collapse_trait_cache_candidates(
                    matched_records,
                    method="cache_exact_ukb_category_field_title",
                    score=0.974,
                    matched_on="UKB category child field title text",
                    validated_if_unique=True,
                    record_resolver=resolve_record,
                )
                if len(ukb_category_field_title_candidates) == 1:
                    candidates.extend(ukb_category_field_title_candidates)

        if not candidates and fuzzy_query_tokens and not has_ontology_exact:
            seed_tokens = pick_seed_tokens(fuzzy_query_tokens, cache_token_freq, max_tokens=6)
            candidate_idxs: set[int] = set()
            for tok in seed_tokens:
                candidate_idxs |= cache_token_index.get(tok, set())
            scored_cache: list[tuple[float, TraitCacheRecord]] = []
            for rec_idx in candidate_idxs:
                rec = records[rec_idx]
                resolved_mapping = resolve_record(rec)
                if not resolved_mapping:
                    continue
                resolved_ids, _resolved_labels = resolved_mapping
                if skip_multi_mapped_cache_fuzzy_candidate(
                    query_text=fuzzy_query_norm,
                    query_tokens=fuzzy_query_tokens,
                    record=rec,
                    mapped_id_count=len(resolved_ids),
                ):
                    continue
                score_a = trait_similarity_score(fuzzy_query_norm, fuzzy_query_tokens, rec.reported_trait)
                score_b = trait_similarity_score(fuzzy_query_norm, fuzzy_query_tokens, rec.lookup_text)
                score = max(score_a, score_b)
                if score <= 0.0:
                    continue
                scored_cache.append((score, rec))
            scored_cache.sort(key=lambda item: item[0], reverse=True)
            best_any.extend(
                collapse_trait_cache_candidates(
                    [rec for _score, rec in scored_cache[: max(3, top_k * 4)]],
                    method="cache_fuzzy_text",
                    score=max(scored_cache[0][0], min_score) if scored_cache else 0.0,
                    matched_on="cache lexical candidate set",
                    validated_if_unique=False,
                    record_resolver=resolve_record,
                )
            )
            for score, rec in scored_cache[: max(3, top_k * 4)]:
                resolved_mapping = resolve_record(rec)
                if not resolved_mapping:
                    continue
                resolved_ids, resolved_labels = resolved_mapping
                mapped_labels = "|".join(resolved_labels)
                modifier_mismatch = trait_modifier_label_mismatch(fuzzy_query_tokens, mapped_labels)
                adjusted_score = max(0.0, score - (0.25 if modifier_mismatch else 0.0))
                if adjusted_score < min_score:
                    continue
                mapped_ids = "|".join(resolved_ids)
                candidates.append(
                    Candidate(
                        efo_id=mapped_ids,
                        label=mapped_labels,
                        score=adjusted_score,
                        matched_via="cache_fuzzy_text",
                        evidence=f"cache row={rec.rownum}; matched reported trait '{rec.reported_trait}'",
                        is_validated=adjusted_score >= 0.93 and not modifier_mismatch,
                    )
                )

        if (not candidates or (input_type == "icd10" and bool(icd10_label))) and (
            query_exact_norms or (ukb_exact_text_keys and not prefer_icd10_from_ukb)
        ):
            exact_term_ids: set[str] = set()
            matched_via_ukb_title = False
            for key in query_exact_norms:
                exact_term_ids |= ontology_exact_index.get(key, set())
            if not prefer_icd10_from_ukb:
                for key in ukb_exact_text_keys:
                    term_ids = ontology_exact_index.get(key, set())
                    if term_ids:
                        matched_via_ukb_title = True
                        exact_term_ids |= term_ids
            if exact_term_ids:
                ordered = sorted(
                    exact_term_ids,
                    key=lambda tid: (
                        TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(tid), 999),
                        tid,
                    ),
                )
                for term_id in ordered[: max(top_k, 5)]:
                    term = ontology_terms[term_id]
                    candidates.append(
                        Candidate(
                            efo_id=term.term_id,
                            label=term.label,
                            score=0.93 if len(ordered) == 1 else 0.88,
                            matched_via=(
                                "efo_obo_exact_ukb_field_title"
                                if matched_via_ukb_title
                                else "efo_obo_exact"
                            ),
                            evidence=(
                                "exact label/synonym match in efo.obo (via UKB field title)"
                                if matched_via_ukb_title
                                else "exact label/synonym match in efo.obo"
                            ),
                            is_validated=len(ordered) == 1,
                        )
                    )

        if not candidates and fuzzy_query_tokens:
            seed_tokens = pick_seed_tokens(fuzzy_query_tokens, ontology_token_freq, max_tokens=7)
            candidate_term_ids: set[str] = set()
            for tok in seed_tokens:
                candidate_term_ids |= ontology_token_index.get(tok, set())
            query_anchor_tokens = trait_anchor_tokens(fuzzy_query_tokens)

            scored_terms: list[tuple[float, str]] = []
            for term_id in candidate_term_ids:
                term = ontology_terms[term_id]
                best = trait_similarity_score(fuzzy_query_norm, fuzzy_query_tokens, term.label)
                for syn in term.synonyms[:10]:
                    best = max(best, trait_similarity_score(fuzzy_query_norm, fuzzy_query_tokens, syn))
                if best <= 0.0:
                    continue
                scored_terms.append((best, term_id))

            scored_terms.sort(
                key=lambda item: (
                    item[0],
                    -TRAIT_PREFERRED_PREFIX_ORDER.get(trait_ontology_prefix(item[1]), 999),
                ),
                reverse=True,
            )

            if scored_terms:
                top_score = scored_terms[0][0]
                second_score = scored_terms[1][0] if len(scored_terms) > 1 else 0.0
                for score, term_id in scored_terms[: max(5, top_k * 6)]:
                    term = ontology_terms[term_id]
                    modifier_mismatch = trait_modifier_label_mismatch(fuzzy_query_tokens, term.label)
                    label_tokens = informative_trait_tokens(normalize_trait_text_for_similarity(term.label))
                    overlap_tokens = fuzzy_query_tokens & label_tokens
                    nonweak_overlap_tokens = trait_nonweak_tokens(overlap_tokens)
                    label_anchor_tokens = trait_anchor_tokens(label_tokens)
                    weak_only_overlap = bool(overlap_tokens) and not nonweak_overlap_tokens
                    anchor_mismatch = bool(query_anchor_tokens) and query_anchor_tokens.isdisjoint(label_anchor_tokens)
                    scope_rank = trait_match_scope_rank(
                        query_exact_keys=query_exact_key_set,
                        term=term,
                        label_text=term.label,
                    )
                    exact_scope_bonus = 0.08 if scope_rank >= TRAIT_SYNONYM_SCOPE_RANK["EXACT"] else 0.0
                    adjusted_score = (
                        score
                        - (0.25 if modifier_mismatch else 0.0)
                        - (0.14 if weak_only_overlap else 0.0)
                        - (0.18 if anchor_mismatch else 0.0)
                        + exact_scope_bonus
                    )
                    adjusted_score = max(0.0, min(1.0, adjusted_score))
                    if adjusted_score < min_score:
                        continue
                    evidence_bits = ["token lexical fallback in efo.obo"]
                    if exact_scope_bonus > 0.0:
                        evidence_bits.append("exact label/synonym scope bonus")
                    if weak_only_overlap:
                        evidence_bits.append("weak-token-only overlap penalty")
                    if anchor_mismatch:
                        evidence_bits.append("anchor-token mismatch penalty")
                    candidates.append(
                        Candidate(
                            efo_id=term.term_id,
                            label=term.label,
                            score=adjusted_score,
                            matched_via="efo_obo_fuzzy",
                            evidence="; ".join(evidence_bits),
                            is_validated=(
                                adjusted_score >= 0.93
                                and (top_score - second_score) >= 0.08
                                and not modifier_mismatch
                                and not weak_only_overlap
                                and not anchor_mismatch
                            ),
                        )
                    )
                top_term = ontology_terms[scored_terms[0][1]]
                top_label_tokens = informative_trait_tokens(normalize_trait_text_for_similarity(top_term.label))
                top_overlap_tokens = fuzzy_query_tokens & top_label_tokens
                top_nonweak_overlap = trait_nonweak_tokens(top_overlap_tokens)
                top_anchor_tokens = trait_anchor_tokens(top_label_tokens)
                top_anchor_mismatch = bool(query_anchor_tokens) and query_anchor_tokens.isdisjoint(top_anchor_tokens)
                top_weak_only_overlap = bool(top_overlap_tokens) and not top_nonweak_overlap
                if not (top_anchor_mismatch and top_weak_only_overlap):
                    top_score_adjusted = max(
                        0.0,
                        min(
                            1.0,
                            top_score
                            - (0.12 if top_weak_only_overlap else 0.0)
                            - (0.16 if top_anchor_mismatch else 0.0),
                        ),
                    )
                    top_evidence = "best lexical candidate from efo.obo"
                    if top_weak_only_overlap:
                        top_evidence = f"{top_evidence}; weak-token-only overlap"
                    if top_anchor_mismatch:
                        top_evidence = f"{top_evidence}; anchor-token mismatch"
                    best_any.append(
                        Candidate(
                            efo_id=top_term.term_id,
                            label=top_term.label,
                            score=top_score_adjusted,
                            matched_via="efo_obo_fuzzy",
                            evidence=top_evidence,
                            is_validated=False,
                        )
                    )

        if not candidates and (ukb_field_ids or ukb_category_ids) and not prefer_icd10_from_ukb:
            fallback_candidates: list[Candidate] = []
            seen_fallback_keys: set[tuple[str, str, str]] = set()
            for field_id in ukb_field_ids:
                for candidate in ukb_field_fallback_cache.get(field_id, ()):
                    dedup_key = (candidate.efo_id, candidate.label, candidate.matched_via)
                    if dedup_key in seen_fallback_keys:
                        continue
                    seen_fallback_keys.add(dedup_key)
                    fallback_candidates.append(candidate)

            if not fallback_candidates:
                fallback_category_ids: list[str] = []
                seen_category_ids: set[str] = set()
                for category_id in ukb_category_ids:
                    cid = normalize(category_id)
                    if cid and cid not in seen_category_ids:
                        seen_category_ids.add(cid)
                        fallback_category_ids.append(cid)
                for field_id in ukb_field_ids:
                    category_id = normalize(ukb_field_to_category.get(field_id, ""))
                    if category_id and category_id not in seen_category_ids:
                        seen_category_ids.add(category_id)
                        fallback_category_ids.append(category_id)

                fallback_key = tuple(sorted(fallback_category_ids))
                if fallback_key:
                    fallback_candidates.extend(build_ukb_category_fallback_candidates(fallback_key))

            if fallback_candidates:
                candidates.extend(fallback_candidates)

        def apply_trait_scale_bias(target_candidates: list[Candidate]) -> None:
            scale_norm = normalize_trait_scale(trait_scale)
            if scale_norm == "auto":
                return
            for candidate in target_candidates:
                branch = trait_candidate_branch(
                    mapped_id_text=candidate.efo_id,
                    mapped_label_text=candidate.label,
                    ontology_terms=ontology_terms,
                    measurement_term_ids=measurement_term_ids,
                )
                candidate.score = max(
                    0.0,
                    min(1.0, candidate.score + trait_scale_score_adjustment(scale_norm, branch)),
                )
                if trait_scale_mismatch(scale_norm, branch):
                    candidate.is_validated = False
                    candidate.evidence = (
                        f"{candidate.evidence}; trait_scale={scale_norm} branch mismatch ({branch})"
                    )
                elif branch != "unknown":
                    candidate.evidence = f"{candidate.evidence}; trait_scale={scale_norm} aligned ({branch})"

        apply_trait_scale_bias(candidates)
        apply_trait_scale_bias(best_any)
        blocked_trait_scale_mismatch_count = 0
        blocked_non_disease_entity_count = 0
        blocked_generic_measurement_found = False
        query_has_disease_intent = (
            has_disease_hint(query_for_matching)
            or has_disease_hint(query)
            or has_ukb_disease_context(ukb_field_label_text)
            or has_ukb_disease_context(ukb_category_label_text)
            or has_ukb_disease_context(ukb_category_path_text)
            or bool(icd10)
            or input_type == "icd10"
        )

        def filter_trait_scale_mismatch(target_candidates: list[Candidate]) -> list[Candidate]:
            nonlocal blocked_trait_scale_mismatch_count
            scale_norm = normalize_trait_scale(trait_scale)
            if scale_norm == "auto":
                return target_candidates
            kept: list[Candidate] = []
            for candidate in target_candidates:
                branch = trait_candidate_branch(
                    mapped_id_text=candidate.efo_id,
                    mapped_label_text=candidate.label,
                    ontology_terms=ontology_terms,
                    measurement_term_ids=measurement_term_ids,
                )
                mismatch = (
                    (scale_norm == "binary" and branch in {"measurement", "mixed"})
                    or (scale_norm == "quantitative" and branch in {"non_measurement", "mixed"})
                )
                if mismatch:
                    blocked_trait_scale_mismatch_count += 1
                    continue
                kept.append(candidate)
            return kept

        candidates = filter_trait_scale_mismatch(candidates)
        best_any = filter_trait_scale_mismatch(best_any)

        def filter_non_disease_entity_candidates(target_candidates: list[Candidate]) -> list[Candidate]:
            nonlocal blocked_non_disease_entity_count
            if not query_has_disease_intent:
                return target_candidates
            kept: list[Candidate] = []
            for candidate in target_candidates:
                if candidate_is_non_disease_entity(candidate, ontology_terms=ontology_terms):
                    blocked_non_disease_entity_count += 1
                    continue
                kept.append(candidate)
            return kept

        candidates = filter_non_disease_entity_candidates(candidates)
        best_any = filter_non_disease_entity_candidates(best_any)

        def apply_parent_icd10_label_demotions(target_candidates: list[Candidate]) -> None:
            if input_type != "icd10" or not icd10_label:
                return
            for candidate in target_candidates:
                if candidate.matched_via not in ICD10_PARENT_FALLBACK_METHODS:
                    continue
                overlap = trait_query_label_overlap_ratio(icd10_label, candidate.label)
                if overlap < 0.45:
                    candidate.is_validated = False
                    candidate.score = min(candidate.score, 0.88)
                    candidate.evidence = (
                        f"{candidate.evidence}; ICD10 parent fallback demoted by code-label mismatch"
                    )

        apply_parent_icd10_label_demotions(candidates)
        apply_parent_icd10_label_demotions(best_any)
        query_anchor_tokens_for_fuzzy = trait_anchor_tokens(
            informative_trait_tokens(normalize_trait_text_for_similarity(query_for_matching))
        )

        def apply_high_risk_fuzzy_demotions(target_candidates: list[Candidate]) -> None:
            for candidate in target_candidates:
                if candidate.matched_via != "efo_obo_fuzzy":
                    continue
                overlap = trait_query_label_overlap_ratio(query_for_matching, candidate.label)
                nonweak_overlap = trait_nonweak_overlap(query_for_matching, candidate.label)
                anchor_overlap = trait_has_anchor_overlap(query_for_matching, candidate.label)
                weak_only_overlap = overlap > 0.0 and not nonweak_overlap
                anchor_mismatch = bool(query_anchor_tokens_for_fuzzy) and not anchor_overlap
                high_risk = (
                    is_generic_trait_label(candidate.label)
                    or overlap < 0.45
                    or weak_only_overlap
                    or (query_has_disease_intent and anchor_mismatch)
                )
                if high_risk:
                    candidate.is_validated = False
                    max_score = 0.89
                    reason_parts = ["high-risk ontology fuzzy match demoted"]
                    if weak_only_overlap:
                        max_score = min(max_score, 0.85)
                        reason_parts.append("weak-token-only overlap")
                    if query_has_disease_intent and anchor_mismatch:
                        max_score = min(max_score, 0.84)
                        reason_parts.append("missing anchor-token agreement")
                    candidate.score = min(candidate.score, max_score)
                    candidate.evidence = f"{candidate.evidence}; {'; '.join(reason_parts)}"

        apply_high_risk_fuzzy_demotions(candidates)
        apply_high_risk_fuzzy_demotions(best_any)

        filtered_candidates: list[Candidate] = []
        for candidate in candidates:
            if candidate_has_generic_measurement_mapping(candidate.efo_id, candidate.label):
                blocked_generic_measurement_found = True
                continue
            filtered_candidates.append(candidate)
        candidates = filtered_candidates
        best_any = [
            candidate
            for candidate in best_any
            if not candidate_has_generic_measurement_mapping(candidate.efo_id, candidate.label)
        ]

        deduped: dict[tuple[str, str, str], Candidate] = {}
        for candidate in candidates:
            key = (candidate.efo_id, candidate.label, candidate.matched_via)
            prior = deduped.get(key)
            if prior is None or candidate.score > prior.score:
                deduped[key] = candidate
        candidates = list(deduped.values())
        query_lexical_norm = normalize_trait_text_for_similarity(query_for_matching) or (
            query_exact_norms[0] if query_exact_norms else norm_key(query_for_matching)
        )
        query_lexical_tokens = informative_trait_tokens(query_lexical_norm)
        candidates.sort(
            key=lambda candidate: (
                candidate.score,
                *trait_candidate_lexical_priority(
                    candidate,
                    query_exact_keys=query_exact_key_set,
                    query_tokens=query_lexical_tokens,
                    ontology_terms=ontology_terms,
                ),
                -via_rank.get(candidate.matched_via, 999),
                -trait_candidate_prefix_rank(candidate.efo_id),
                candidate.efo_id,
            ),
            reverse=True,
        )

        emitted = candidates[: max(1, top_k)]
        used_best_any_fallback = False
        auto_force_map_best = (
            input_type == "ukb_field"
            and bool(parsed_ukb_field_id)
            and "#" in query
            and query_has_disease_intent
        )
        if not emitted and (force_map_best or auto_force_map_best) and best_any:
            emitted = [best_any[0]]
            used_best_any_fallback = True

        if not emitted:
            generated_rows.append(
                {
                    "input_row_id": input_row_id,
                    "input_query": query,
                    "negation": "yes" if query_has_negation else "no",
                    "input_type": input_type,
                    "input_trait_scale": trait_scale,
                    "input_source_code": source_code,
                    "input_source_data_type": source_data_type,
                    "input_icd10": icd10,
                    "input_icd10_label": icd10_label,
                    "input_phecode": phecode,
                    "input_ukb_field_id": ukb_field_id_text,
                    "input_ukb_field_label": ukb_field_label_text,
                    "input_ukb_category_id": ukb_category_id_text,
                    "input_ukb_category_label": ukb_category_label_text,
                    "input_ukb_category_path": ukb_category_path_text,
                    "mapped_trait_id": "",
                    "mapped_trait_label": "",
                    "confidence": "0.000",
                    "matched_via": "none",
                    "matched_on": "",
                    "source_file": "",
                    "validation": "not_mapped",
                    "qc_review_flag": "yes",
                    "specificity_loss_flag": "no",
                    "specificity_loss_notes": "",
                    "term_not_in_efo": "",
                    "mondo_missing_ids": "",
                    "evidence": (
                        "no cache/efo candidate above threshold"
                        if not query_has_negation
                        else "no cache/efo candidate above threshold; negation cue in input query"
                    ),
                    "provenance": "",
                }
            )
            if blocked_generic_measurement_found:
                generated_rows[-1]["evidence"] = (
                    f"{generated_rows[-1]['evidence']}; generic measurement term candidate blocked"
                )
            if blocked_trait_scale_mismatch_count > 0:
                generated_rows[-1]["evidence"] = (
                    f"{generated_rows[-1]['evidence']}; blocked trait-scale mismatched candidates="
                    f"{blocked_trait_scale_mismatch_count}"
                )
            if blocked_non_disease_entity_count > 0:
                generated_rows[-1]["evidence"] = (
                    f"{generated_rows[-1]['evidence']}; blocked non-disease entity candidates="
                    f"{blocked_non_disease_entity_count}"
                )
        else:
            for candidate in emitted:
                source_file = (
                    DEFAULT_UKB_FIELD_SUPPLEMENT_CACHE.name
                    if candidate.matched_via == "cache_exact_ukb_field_supplement"
                    else DEFAULT_ICD10_SUPPLEMENT_CACHE.name
                    if candidate.matched_via == "cache_exact_icd10_supplement"
                    else trait_cache_path.name
                    if candidate.matched_via.startswith("cache_")
                    else efo_obo_path.name
                )
                mapped_ids = "|".join(
                    format_ontology_id_for_output(part)
                    for part in split_multi_ids(candidate.efo_id)
                )
                if not mapped_ids:
                    mapped_ids = format_ontology_id_for_output(candidate.efo_id)
                mapped_labels = normalize(candidate.label)
                specificity_loss_flag, specificity_loss_notes = trait_specificity_loss_flag(
                    query_text=query_for_matching,
                    icd10=icd10,
                    matched_via=candidate.matched_via,
                    mapped_id_text=mapped_ids,
                    mapped_label_text=mapped_labels,
                    ontology_terms=ontology_terms,
                )
                term_not_in_efo, mondo_missing_ids = mondo_import_qc(
                    candidate.efo_id,
                    ontology_terms=ontology_terms,
                )
                validated = (
                    candidate.is_validated
                    and candidate.score >= min_score
                    and not term_not_in_efo
                    and not query_has_negation
                )
                validation = "validated" if validated else "review_required"
                matched_on = (
                    "ICD10"
                    if candidate.matched_via == "cache_exact_icd10"
                    else "ICD10 (supplemental cache)"
                    if candidate.matched_via == "cache_exact_icd10_supplement"
                    else "ICD10 parent fallback"
                    if candidate.matched_via == "cache_parent_icd10"
                    else "ICD10 parent fallback (supplemental cache)"
                    if candidate.matched_via == "cache_parent_icd10_supplement"
                    else "PheCode"
                    if candidate.matched_via == "cache_exact_phecode"
                    else "ICD10 xref in efo.obo"
                    if candidate.matched_via == "efo_obo_exact_icd10_xref"
                    else "ICD10 parent xref in efo.obo"
                    if candidate.matched_via == "efo_obo_parent_icd10_xref"
                    else "UKB data field"
                    if candidate.matched_via == "cache_exact_ukb_field"
                    else "UKB data field (supplemental cache)"
                    if candidate.matched_via == "cache_exact_ukb_field_supplement"
                    else "UKB field title"
                    if candidate.matched_via == "cache_exact_ukb_field_title"
                    else "UKB category"
                    if candidate.matched_via == "cache_exact_ukb_category"
                    else "UKB category title"
                    if candidate.matched_via == "cache_exact_ukb_category_title"
                    else "UKB category child field title"
                    if candidate.matched_via == "cache_exact_ukb_category_field_title"
                    else "UKB category fallback"
                    if candidate.matched_via == "efo_obo_exact_ukb_category"
                    else "Reported trait"
                    if candidate.matched_via.startswith("cache_")
                    else "EFO/OBO label or synonym"
                )
                effective_input_type = input_type
                if "ukb_" in candidate.matched_via:
                    effective_input_type = "ukb_field"
                primary_input_value = (
                    icd10
                    if effective_input_type == "icd10"
                    else phecode
                    if effective_input_type == "phecode"
                    else ukb_field_id_text
                    if effective_input_type == "ukb_field" and ukb_field_id_text
                    else query
                )
                provenance = make_trait_provenance(
                    method=candidate.matched_via,
                    input_type=effective_input_type,
                    trait_scale=trait_scale,
                    input_value=primary_input_value,
                    matched_on=matched_on,
                    mapped_ids=mapped_ids,
                    mapped_labels=mapped_labels,
                    source_file=source_file,
                    score=candidate.score,
                    validated=validated,
                )
                evidence = candidate.evidence
                if term_not_in_efo:
                    evidence = f"{evidence}; MONDO term missing from EFO import"
                if query_has_negation:
                    evidence = f"{evidence}; negation cue in input query"
                if specificity_loss_flag == "yes" and specificity_loss_notes:
                    evidence = f"{evidence}; {specificity_loss_notes}"
                if used_best_any_fallback and not force_map_best:
                    evidence = f"{evidence}; auto best-candidate fallback for coded UKB trait"
                generated_rows.append(
                    {
                        "input_row_id": input_row_id,
                        "input_query": query,
                        "negation": "yes" if query_has_negation else "no",
                        "input_type": effective_input_type,
                        "input_trait_scale": trait_scale,
                        "input_source_code": source_code,
                        "input_source_data_type": source_data_type,
                        "input_icd10": icd10,
                        "input_icd10_label": icd10_label,
                        "input_phecode": phecode,
                        "input_ukb_field_id": ukb_field_id_text,
                        "input_ukb_field_label": ukb_field_label_text,
                        "input_ukb_category_id": ukb_category_id_text,
                        "input_ukb_category_label": ukb_category_label_text,
                        "input_ukb_category_path": ukb_category_path_text,
                        "mapped_trait_id": mapped_ids,
                        "mapped_trait_label": mapped_labels,
                        "confidence": f"{max(0.0, min(1.0, candidate.score)):.3f}",
                        "matched_via": candidate.matched_via,
                        "matched_on": matched_on,
                        "source_file": source_file,
                        "validation": validation,
                        "qc_review_flag": "no" if validation == "validated" else "yes",
                        "specificity_loss_flag": specificity_loss_flag,
                        "specificity_loss_notes": specificity_loss_notes,
                        "term_not_in_efo": term_not_in_efo,
                        "mondo_missing_ids": mondo_missing_ids,
                        "evidence": evidence,
                        "provenance": provenance,
                    }
                )

        if memoize_queries and len(query_result_cache) < max(0, query_cache_max_entries):
            query_result_cache[cache_key] = tuple(dict(row) for row in generated_rows)

        for row in generated_rows:
            if row_callback is not None:
                row_callback(dict(row))
            if collect_rows:
                rows.append(dict(row))
        progress.update()

    return rows


def write_trait_tsv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRAIT_OUTPUT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_trait_review_tsv(rows: list[dict[str, str]], path: Path) -> int:
    review_rows = [row for row in rows if row.get("validation") == "review_required"]
    write_trait_tsv(review_rows, path)
    return len(review_rows)


class TraitMapStreamingWriter:
    """Write trait-map rows incrementally so partial output survives interruptions."""

    def __init__(self, *, output_path: Path, review_path: Path | None, flush_every: int) -> None:
        self.output_path = output_path
        self.review_path = review_path
        self.flush_every = max(1, int(flush_every))
        self.output_handle: Any | None = None
        self.review_handle: Any | None = None
        self.output_writer: Any | None = None
        self.review_writer: Any | None = None
        self.rows_written = 0
        self.review_rows_written = 0
        self._pending_writes = 0

    def __enter__(self) -> "TraitMapStreamingWriter":
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.output_handle = self.output_path.open("w", encoding="utf-8", newline="")
        self.output_writer = csv.DictWriter(
            self.output_handle,
            fieldnames=TRAIT_OUTPUT_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
        )
        self.output_writer.writeheader()

        if self.review_path is not None:
            self.review_path.parent.mkdir(parents=True, exist_ok=True)
            self.review_handle = self.review_path.open("w", encoding="utf-8", newline="")
            self.review_writer = csv.DictWriter(
                self.review_handle,
                fieldnames=TRAIT_OUTPUT_FIELDS,
                delimiter="\t",
                extrasaction="ignore",
            )
            self.review_writer.writeheader()
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self.close()

    def write_row(self, row: dict[str, str]) -> None:
        if self.output_writer is None:
            raise RuntimeError("TraitMapStreamingWriter is not open")
        self.output_writer.writerow(row)
        self.rows_written += 1

        if self.review_writer is not None and row.get("validation") == "review_required":
            self.review_writer.writerow(row)
            self.review_rows_written += 1

        self._pending_writes += 1
        if self._pending_writes >= self.flush_every:
            self.flush()

    def flush(self) -> None:
        if self.output_handle is not None:
            self.output_handle.flush()
        if self.review_handle is not None:
            self.review_handle.flush()
        self._pending_writes = 0

    def close(self) -> None:
        self.flush()
        if self.review_handle is not None:
            self.review_handle.close()
            self.review_handle = None
            self.review_writer = None
        if self.output_handle is not None:
            self.output_handle.close()
            self.output_handle = None
            self.output_writer = None


def load_term_cache(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []

    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            efo_id = normalize_id(row.get("efo_id") or "")
            label = normalize(row.get("label") or "")
            if not efo_id or not label:
                continue
            if not ALLOWED_ONTOLOGY_ID_RE.match(efo_id):
                continue
            syn_field = row.get("synonyms") or ""
            syns = [normalize(x) for x in syn_field.split("|") if normalize(x)]
            metabolite_ids: set[str] = set()
            for col in (
                "metabolite_ids",
                "xrefs",
                "chebi_ids",
                "hmdb_ids",
                "kegg_ids",
                "chebi_id",
                "hmdb_id",
                "kegg_id",
            ):
                raw_val = row.get(col) or ""
                if not raw_val:
                    continue
                for part in parse_pipe_list(raw_val):
                    for met_id in extract_metabolite_ids_from_text(part):
                        metabolite_ids.add(met_id)
            rows.append(
                {
                    "efo_id": efo_id,
                    "label": label,
                    "synonyms": syns,
                    "metabolite_ids": sorted(metabolite_ids),
                }
            )
    return rows


def analyte_cache_entry_is_trusted(*, source: str, evidence: str) -> bool:
    source_key = norm_key(source)
    evidence_key = norm_key(evidence)

    if "manual" in source_key or "curated" in source_key:
        return True

    # Reject known weak/autogenerated lexical signals from older mapper runs.
    unsafe_markers = (
        "token retrieval + lexical score",
        "source=local-term-token",
        "forced-best",
        "withheld",
        "needs_new_term",
        "semantic-exact-needs-new-term",
    )
    if any(marker in evidence_key for marker in unsafe_markers):
        return False

    safe_source_prefixes = (
        "local-term-exact",
        "local-term-metabolite-xref",
        "direct-id",
        "local-term-cache",
    )
    if any(source_key.startswith(prefix) for prefix in safe_source_prefixes):
        return True

    safe_evidence_prefixes = (
        "exact normalized term match:",
        "metabolite external-id match:",
        "direct input efo/oba id",
    )
    if any(evidence_key.startswith(prefix) for prefix in safe_evidence_prefixes):
        return True

    # Legacy writeback rows often embed the original source in the evidence blob.
    if "source=local-term-cache" in evidence_key and "source=local-term-token" not in evidence_key:
        return True

    return False


def load_analyte_cache(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []

    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            analyte = norm_key(row.get("analyte_key") or "")
            efo_id = normalize_id(row.get("efo_id") or "")
            if not analyte or not efo_id:
                continue
            try:
                confidence = float(row.get("confidence") or "0.9")
            except ValueError:
                confidence = 0.9
            confidence = min(1.0, max(0.0, confidence))
            source = normalize(row.get("source") or "cache")
            evidence = normalize(row.get("evidence") or "cached mapping")
            if not analyte_cache_entry_is_trusted(source=source, evidence=evidence):
                continue
            rows.append(
                {
                    "analyte_key": analyte,
                    "efo_id": efo_id,
                    "label": normalize(row.get("label") or ""),
                    "confidence": confidence,
                    "validation": norm_key(row.get("validation") or ""),
                    "source": source,
                    "evidence": evidence,
                }
            )
    return rows


def load_uniprot_aliases(path: Path) -> dict[str, list[str]]:
    if not path.exists():
        return {}

    result: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            acc = normalize((row.get("accession") or "")).upper()
            if not acc:
                continue

            alias_values: list[str] = []
            aliases_col = row.get("aliases") or ""
            if aliases_col:
                alias_values.extend([normalize(x) for x in aliases_col.split("|") if normalize(x)])

            single_alias = row.get("alias") or ""
            if single_alias and normalize(single_alias):
                alias_values.append(normalize(single_alias))

            gene = row.get("gene") or ""
            if gene and normalize(gene):
                alias_values.append(normalize(gene))

            gene_ids_col = row.get("gene_ids") or row.get("gene_id") or ""
            if gene_ids_col:
                for gid_raw in gene_ids_col.split("|"):
                    gid = normalize(gid_raw)
                    if not gid:
                        continue
                    for gid_variant in gene_id_alias_variants(gid, include_unprefixed_numeric=False):
                        alias_values.append(gid_variant)

            if not alias_values:
                continue
            result.setdefault(acc, set()).update(alias_values)

    return {k: sorted(v) for k, v in result.items()}


def load_uniprot_alias_records(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    if not path.exists():
        return records

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            acc = normalize((row.get("accession") or "")).upper()
            if not acc:
                continue
            aliases = [normalize(x) for x in (row.get("aliases") or "").split("|") if normalize(x)]
            gene = normalize(row.get("gene") or "")
            gene_ids_raw = [
                normalize(x)
                for x in (row.get("gene_ids") or row.get("gene_id") or "").split("|")
                if normalize(x)
            ]
            gene_ids: set[str] = set()
            for gid in gene_ids_raw:
                gene_ids.update(gene_id_alias_variants(gid, include_unprefixed_numeric=False))
            source = normalize(row.get("source") or "")
            # Merge duplicate rows by accession if present.
            if acc in records:
                prior_aliases = [normalize(x) for x in records[acc].get("aliases", "").split("|") if normalize(x)]
                merged = sorted(set(prior_aliases) | set(aliases))
                records[acc]["aliases"] = "|".join(merged)
                prior_gene_ids = {
                    normalize(x)
                    for x in records[acc].get("gene_ids", "").split("|")
                    if normalize(x)
                }
                merged_gene_ids = sorted(prior_gene_ids | gene_ids)
                records[acc]["gene_ids"] = "|".join(merged_gene_ids)
                if not records[acc].get("gene") and gene:
                    records[acc]["gene"] = gene
                continue
            records[acc] = {
                "accession": acc,
                "aliases": "|".join(sorted(set(aliases))),
                "gene": gene,
                "gene_ids": "|".join(sorted(gene_ids)),
                "source": source,
            }
    return records


def write_uniprot_alias_records(path: Path, records: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["accession", "aliases", "gene", "gene_ids", "source"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for acc in sorted(records):
            row = records[acc]
            aliases = [normalize(x) for x in (row.get("aliases") or "").split("|") if normalize(x)]
            gene_ids_raw = [normalize(x) for x in (row.get("gene_ids") or "").split("|") if normalize(x)]
            gene_ids: set[str] = set()
            for gid in gene_ids_raw:
                gene_ids.update(gene_id_alias_variants(gid, include_unprefixed_numeric=False))
            writer.writerow(
                {
                    "accession": acc,
                    "aliases": "|".join(sorted(set(aliases))),
                    "gene": normalize(row.get("gene") or ""),
                    "gene_ids": "|".join(sorted(gene_ids)),
                    "source": normalize(row.get("source") or ""),
                }
            )


def build_lightweight_uniprot_alias_cache(
    *,
    source_path: Path,
    output_path: Path,
    reviewed_only: bool = True,
    require_gene: bool = True,
    add_source_tag: str = "lightweight-cache-build",
) -> tuple[int, int]:
    records = load_uniprot_alias_records(source_path)
    if not records:
        raise ValueError(f"No UniProt alias rows found in source: {source_path}")

    lightweight: dict[str, dict[str, str]] = {}
    for acc in sorted(records):
        row = records[acc]
        if reviewed_only and not is_reviewed_like_accession(acc):
            continue
        gene = normalize(row.get("gene") or "")
        if require_gene and not gene:
            continue

        aliases = set(parse_pipe_list(row.get("aliases") or ""))
        if gene:
            aliases.add(gene)
        clean_aliases = sorted({normalize(a) for a in aliases if normalize(a)})

        gene_ids: set[str] = set()
        for raw in parse_pipe_list(row.get("gene_ids") or row.get("gene_id") or ""):
            gene_ids.update(gene_id_alias_variants(raw, include_unprefixed_numeric=False))

        source_values = set(parse_pipe_list(row.get("source") or ""))
        if add_source_tag:
            source_values.add(add_source_tag)

        lightweight[acc] = {
            "accession": acc,
            "aliases": "|".join(clean_aliases),
            "gene": gene,
            "gene_ids": "|".join(sorted(gene_ids)),
            "source": "|".join(sorted(source_values)),
        }

    write_uniprot_alias_records(output_path, lightweight)
    return len(records), len(lightweight)


def parse_pipe_list(value: str) -> list[str]:
    return [normalize(x) for x in (value or "").split("|") if normalize(x)]


def normalize_metabolite_ids(values: list[str]) -> list[str]:
    ids: list[str] = []
    for raw in values:
        mid = canonical_metabolite_id(raw)
        if mid:
            ids.append(mid)
    return sorted(set(ids))


def default_metabolite_concept_id(label: str, ids: list[str], idx: int = 0) -> str:
    if ids:
        return ids[0]
    key = re.sub(r"[^a-z0-9]+", "_", norm_key(label)).strip("_")
    if not key:
        key = f"concept_{idx + 1}"
    return f"META:{key}"


def load_metabolite_alias_records(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    if not path.exists():
        return records

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_i, row in enumerate(reader):
            concept_raw = (
                row.get("concept_id")
                or row.get("metabolite_concept_id")
                or row.get("metabolite_id")
                or row.get("id")
                or ""
            )
            label = normalize_metabolite_alias_text(row.get("primary_label") or row.get("label") or row.get("name") or "")
            aliases = [
                normalize_metabolite_alias_text(alias)
                for alias in parse_pipe_list(row.get("aliases") or row.get("synonyms") or "")
            ]
            aliases = [alias for alias in aliases if alias and keep_metabolite_alias(alias)]

            hmdb_ids = parse_pipe_list(row.get("hmdb_ids") or row.get("hmdb_id") or "")
            chebi_ids = [
                f"CHEBI:{val}" if re.fullmatch(r"\d{1,7}", val) else val
                for val in parse_pipe_list(row.get("chebi_ids") or row.get("chebi_id") or "")
            ]
            kegg_ids = parse_pipe_list(row.get("kegg_ids") or row.get("kegg_id") or "")
            generic_ids = parse_pipe_list(row.get("ids") or row.get("metabolite_ids") or "")
            concept_seed = (
                f"CHEBI:{concept_raw}"
                if re.fullmatch(r"\d{1,7}", normalize(concept_raw))
                else concept_raw
            )
            all_ids = normalize_metabolite_ids(
                [concept_seed, *hmdb_ids, *chebi_ids, *kegg_ids, *generic_ids]
            )

            concept_id = normalize(concept_seed)
            concept_id_canonical = canonical_metabolite_id(concept_id) if concept_id else ""
            if concept_id_canonical:
                concept_id = concept_id_canonical
            if not concept_id:
                concept_id = default_metabolite_concept_id(label, all_ids, idx=row_i)

            source = normalize(row.get("source") or "metabolite-aliases")
            merged_aliases = prune_metabolite_aliases(
                label=label,
                aliases=aliases,
                ids=all_ids,
            )
            record = records.get(concept_id)
            if record is None:
                records[concept_id] = {
                    "concept_id": concept_id,
                    "primary_label": label,
                    "aliases": "|".join(merged_aliases),
                    "hmdb_ids": "|".join([mid for mid in all_ids if metabolite_id_bucket(mid) == "hmdb"]),
                    "chebi_ids": "|".join([mid for mid in all_ids if metabolite_id_bucket(mid) == "chebi"]),
                    "kegg_ids": "|".join([mid for mid in all_ids if metabolite_id_bucket(mid) == "kegg"]),
                    "source": source,
                }
                continue

            existing_aliases = set(parse_pipe_list(record.get("aliases") or ""))
            existing_aliases.update(merged_aliases)
            record["aliases"] = "|".join(
                prune_metabolite_aliases(
                    label=normalize_metabolite_alias_text(record.get("primary_label") or ""),
                    aliases=sorted(existing_aliases),
                    ids=normalize_metabolite_ids(
                        [
                            *parse_pipe_list(record.get("hmdb_ids") or ""),
                            *parse_pipe_list(record.get("chebi_ids") or ""),
                            *parse_pipe_list(record.get("kegg_ids") or ""),
                        ]
                    ),
                )
            )
            if not normalize(record.get("primary_label") or "") and label:
                record["primary_label"] = label

            existing_ids = normalize_metabolite_ids(
                [
                    *parse_pipe_list(record.get("hmdb_ids") or ""),
                    *parse_pipe_list(record.get("chebi_ids") or ""),
                    *parse_pipe_list(record.get("kegg_ids") or ""),
                    *all_ids,
                ]
            )
            record["hmdb_ids"] = "|".join([mid for mid in existing_ids if metabolite_id_bucket(mid) == "hmdb"])
            record["chebi_ids"] = "|".join([mid for mid in existing_ids if metabolite_id_bucket(mid) == "chebi"])
            record["kegg_ids"] = "|".join([mid for mid in existing_ids if metabolite_id_bucket(mid) == "kegg"])
            if source and source not in parse_pipe_list(record.get("source") or ""):
                record["source"] = "|".join(sorted(set(parse_pipe_list(record.get("source") or "")) | {source}))

    return records


def write_metabolite_alias_records(path: Path, records: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["concept_id", "primary_label", "aliases", "hmdb_ids", "chebi_ids", "kegg_ids", "source"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for concept_id in sorted(records):
            row = records[concept_id]
            primary_label = normalize_metabolite_alias_text(row.get("primary_label") or "")
            ids = normalize_metabolite_ids(
                [
                    *parse_pipe_list(row.get("hmdb_ids") or ""),
                    *parse_pipe_list(row.get("chebi_ids") or ""),
                    *parse_pipe_list(row.get("kegg_ids") or ""),
                ]
            )
            aliases = prune_metabolite_aliases(
                label=primary_label,
                aliases=parse_pipe_list(row.get("aliases") or ""),
                ids=ids,
            )
            writer.writerow(
                {
                    "concept_id": concept_id,
                    "primary_label": primary_label,
                    "aliases": "|".join(aliases),
                    "hmdb_ids": "|".join([mid for mid in ids if metabolite_id_bucket(mid) == "hmdb"]),
                    "chebi_ids": "|".join([mid for mid in ids if metabolite_id_bucket(mid) == "chebi"]),
                    "kegg_ids": "|".join([mid for mid in ids if metabolite_id_bucket(mid) == "kegg"]),
                    "source": normalize(row.get("source") or ""),
                }
            )


def download_to_file(url: str, destination: Path, timeout: float) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": HTTP_USER_AGENT,
            "Accept": "*/*",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response, destination.open("wb") as out:  # nosec B310
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    return destination


def filename_from_url(url: str, fallback: str) -> str:
    path = urllib.parse.urlparse(url).path
    name = Path(path).name
    return name or fallback


def decompress_gzip_file(source_gz: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source_gz, "rb") as src, destination.open("wb") as dst:
        while True:
            chunk = src.read(1024 * 1024)
            if not chunk:
                break
            dst.write(chunk)
    return destination


def looks_like_obo(path: Path, max_lines: int = 20000) -> bool:
    if not path.exists() or path.stat().st_size <= 0:
        return False
    seen_format = False
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for idx, raw_line in enumerate(handle):
                line = raw_line.strip()
                if line.startswith("format-version:"):
                    seen_format = True
                if line == "[Term]":
                    return True
                if idx >= max_lines:
                    break
    except OSError:
        return False
    return seen_format


def ensure_efo_obo_available(
    *,
    efo_obo_path: Path,
    bundled_url: str,
    timeout: float,
) -> tuple[Path, str]:
    if looks_like_obo(efo_obo_path):
        return efo_obo_path, "local"

    local_gz = Path(str(efo_obo_path) + ".gz")
    if local_gz.exists() and local_gz.stat().st_size > 0:
        decompress_gzip_file(local_gz, efo_obo_path)
        if not looks_like_obo(efo_obo_path):
            raise ValueError(
                f"Local gzip source exists but did not produce a valid OBO file: {local_gz}"
            )
        return efo_obo_path, "local-gz"

    # If caller passed a missing/relative custom path, prefer the repository's
    # bundled local OBO before attempting any network download.
    if efo_obo_path != DEFAULT_EFO_OBO_LOCAL and looks_like_obo(DEFAULT_EFO_OBO_LOCAL):
        return DEFAULT_EFO_OBO_LOCAL, "local-default"

    url = normalize(bundled_url)
    if not url:
        raise FileNotFoundError(
            f"EFO OBO not found at {efo_obo_path}, and no bundled URL was provided."
        )

    fallback_name = f"{efo_obo_path.name}.gz"
    download_name = filename_from_url(url, fallback_name)
    download_path = efo_obo_path.parent / download_name
    try:
        download_to_file(url, download_path, timeout=timeout)
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            "Failed to download bundled EFO OBO "
            f"(url={url}, http_status={exc.code}). "
            f"Provide a valid local file via --efo-obo (current: {efo_obo_path}) "
            "or override --efo-obo-bundled-url."
        ) from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(
            "Failed to download bundled EFO OBO "
            f"(url={url}, reason={exc.reason}). "
            f"Provide a valid local file via --efo-obo (current: {efo_obo_path}) "
            "or override --efo-obo-bundled-url."
        ) from exc
    if not download_path.exists() or download_path.stat().st_size <= 0:
        raise ValueError(f"Downloaded EFO OBO bundle is empty: {download_path}")

    if download_path.suffix.lower() == ".gz":
        decompress_gzip_file(download_path, efo_obo_path)
        if not looks_like_obo(efo_obo_path):
            raise ValueError(
                f"Downloaded gzip did not produce a valid OBO file: {download_path}"
            )
        return efo_obo_path, "downloaded-gz"

    if download_path != efo_obo_path:
        shutil.copy2(download_path, efo_obo_path)
    if not looks_like_obo(efo_obo_path):
        raise ValueError(f"Downloaded file is not a valid OBO: {download_path}")
    return efo_obo_path, "downloaded"


def looks_like_tsv_with_columns(path: Path, required_columns: tuple[str, ...]) -> bool:
    if not path.exists() or path.stat().st_size <= 0:
        return False
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
    except (OSError, StopIteration):
        return False
    if not header:
        return False
    normalized_header = {normalize(col) for col in header if normalize(col)}
    return all(normalize(col) in normalized_header for col in required_columns)


def ensure_tsv_available(
    *,
    path: Path,
    download_url: str,
    timeout: float,
    required_columns: tuple[str, ...],
    label: str,
) -> tuple[Path, str]:
    if looks_like_tsv_with_columns(path, required_columns):
        return path, "local"

    url = normalize(download_url)
    if not url:
        raise FileNotFoundError(f"{label} not found at {path}, and no download URL was provided.")

    try:
        download_to_file(url, path, timeout=timeout)
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"Failed to download {label} (url={url}, http_status={exc.code}). "
            f"Provide a valid local file via {path} or override the download URL."
        ) from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"Failed to download {label} (url={url}, reason={exc.reason}). "
            f"Provide a valid local file via {path} or override the download URL."
        ) from exc

    if not looks_like_tsv_with_columns(path, required_columns):
        raise ValueError(f"Downloaded {label} did not contain expected TSV columns: {path}")
    return path, "downloaded"


def ensure_ukb_metadata_available(
    *,
    ukb_field_catalog_path: Path,
    ukb_field_metadata_path: Path,
    ukb_category_catalog_path: Path,
    ukb_category_tree_path: Path,
    field_catalog_url: str,
    field_metadata_url: str,
    category_catalog_url: str,
    category_tree_url: str,
    timeout: float,
) -> dict[str, tuple[Path, str]]:
    return {
        "ukb_field_catalog": ensure_tsv_available(
            path=ukb_field_catalog_path,
            download_url=field_catalog_url,
            timeout=timeout,
            required_columns=("field_id", "title"),
            label="UKB field catalog",
        ),
        "ukb_field_metadata": ensure_tsv_available(
            path=ukb_field_metadata_path,
            download_url=field_metadata_url,
            timeout=timeout,
            required_columns=("field_id", "title", "main_category"),
            label="UKB field metadata",
        ),
        "ukb_category_catalog": ensure_tsv_available(
            path=ukb_category_catalog_path,
            download_url=category_catalog_url,
            timeout=timeout,
            required_columns=("category_id", "title"),
            label="UKB category catalog",
        ),
        "ukb_category_tree": ensure_tsv_available(
            path=ukb_category_tree_path,
            download_url=category_tree_url,
            timeout=timeout,
            required_columns=("parent_id", "child_id"),
            label="UKB category tree",
        ),
    }


def extract_first_xml_from_zip(zip_path: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        xml_members = [name for name in zf.namelist() if name.lower().endswith(".xml")]
        if not xml_members:
            raise ValueError(f"No XML file found in archive: {zip_path}")
        # Prefer obvious HMDB metabolite XML names when present.
        xml_members_sorted = sorted(
            xml_members,
            key=lambda name: (0 if "metabolite" in name.lower() else 1, len(name), name),
        )
        member = xml_members_sorted[0]
        # Zip members can contain nested directories; flatten to filename in output_dir.
        out_path = output_dir / Path(member).name
        with zf.open(member) as src, out_path.open("wb") as dst:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
    return out_path


def detect_local_hmdb_source(
    *,
    preferred_path: Path | None = None,
    download_dir: Path | None = None,
    local_dir: Path | None = None,
) -> Path | None:
    candidates: list[Path] = []
    if preferred_path is not None:
        candidates.append(preferred_path)
    for base in [local_dir, download_dir]:
        if base is None:
            continue
        for filename in DEFAULT_HMDB_LOCAL_FILENAMES:
            candidates.append(base / filename)
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate.resolve()) if candidate.exists() else str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.exists() and candidate.is_file():
            return candidate
    return None


def prepare_metabolite_source_paths(
    *,
    hmdb_xml_path: Path | None,
    chebi_names_path: Path | None,
    kegg_conv_path: Path | None,
    download_common_sources: bool,
    download_hmdb: bool,
    download_dir: Path,
    hmdb_url: str,
    chebi_names_url: str,
    kegg_conv_url: str,
    timeout: float,
) -> tuple[Path | None, Path | None, Path | None, list[str]]:
    notes: list[str] = []
    hmdb_path = hmdb_xml_path
    chebi_path = chebi_names_path
    kegg_path = kegg_conv_path

    download_chebi = download_common_sources and chebi_path is None
    download_kegg = download_common_sources and kegg_path is None
    should_download_hmdb = download_hmdb and hmdb_path is None

    if should_download_hmdb:
        hmdb_download = download_dir / filename_from_url(hmdb_url, "hmdb_metabolites.zip")
        try:
            download_to_file(hmdb_url, hmdb_download, timeout=timeout)
            if hmdb_download.suffix.lower() == ".zip":
                try:
                    hmdb_path = extract_first_xml_from_zip(hmdb_download, download_dir)
                except zipfile.BadZipFile as exc:
                    raise RuntimeError(
                        f"HMDB download from {hmdb_url} did not produce a valid ZIP archive "
                        "(common when a login/anti-bot HTML page is returned). "
                        "Provide a valid local HMDB file via --hmdb-xml or update --hmdb-url "
                        "to your pinned bundle."
                    ) from exc
            else:
                hmdb_path = hmdb_download
            notes.append(f"downloaded HMDB -> {hmdb_path}")
        except Exception as exc:
            raise RuntimeError(
                f"HMDB download failed from {hmdb_url}. "
                f"Reason: {exc.__class__.__name__}: {exc}. "
                "Provide --hmdb-xml <local-file> or set --hmdb-url to your pinned bundle."
            ) from exc

    if download_chebi:
        chebi_download = download_dir / filename_from_url(chebi_names_url, "chebi_names.tsv.gz")
        try:
            download_to_file(chebi_names_url, chebi_download, timeout=timeout)
            chebi_path = chebi_download
            notes.append(f"downloaded ChEBI names -> {chebi_path}")
        except Exception as exc:
            raise RuntimeError(
                f"ChEBI names download failed from {chebi_names_url}. "
                f"Reason: {exc.__class__.__name__}: {exc}. "
                "Provide --chebi-names <local-file>."
            ) from exc

    if download_kegg:
        kegg_download = download_dir / filename_from_url(kegg_conv_url, "kegg_conv_chebi_compound.tsv")
        try:
            download_to_file(kegg_conv_url, kegg_download, timeout=timeout)
            kegg_path = kegg_download
            notes.append(f"downloaded KEGG conv -> {kegg_path}")
        except Exception as exc:
            raise RuntimeError(
                f"KEGG conv download failed from {kegg_conv_url}. "
                f"Reason: {exc.__class__.__name__}: {exc}. "
                "Provide --kegg-conv <local-file>."
            ) from exc

    return hmdb_path, chebi_path, kegg_path, notes


def open_text_with_optional_gzip(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def xml_local_name(tag: str) -> str:
    if "}" in tag:
        return tag.split("}", 1)[1]
    return tag


def upsert_metabolite_builder_record(
    records: dict[str, dict[str, str]],
    id_to_concept: dict[str, str],
    *,
    concept_id: str,
    label: str,
    aliases: list[str],
    ids: list[str],
    source: str,
) -> None:
    clean_label = normalize_metabolite_alias_text(label)
    clean_aliases = [normalize_metabolite_alias_text(a) for a in aliases]
    clean_aliases = [a for a in clean_aliases if a and keep_metabolite_alias(a)]
    canonical_ids = normalize_metabolite_ids([concept_id, *ids])
    cid = canonical_metabolite_id(concept_id) if concept_id else ""
    if cid:
        concept_id = cid
    merged_concept_id = ""
    for mid in canonical_ids:
        existing = id_to_concept.get(mid, "")
        if existing:
            merged_concept_id = existing
            break
    concept_id = normalize(merged_concept_id or concept_id) or default_metabolite_concept_id(clean_label, canonical_ids)
    source_clean = normalize(source) or "metabolite-build"

    rec = records.get(concept_id)
    if rec is None:
        rec = {
            "concept_id": concept_id,
            "primary_label": clean_label,
            "aliases": "",
            "hmdb_ids": "",
            "chebi_ids": "",
            "kegg_ids": "",
            "source": source_clean,
        }
        records[concept_id] = rec

    if not normalize(rec.get("primary_label") or "") and clean_label:
        rec["primary_label"] = clean_label

    merged_aliases = set(parse_pipe_list(rec.get("aliases") or ""))
    merged_aliases.update(clean_aliases)
    if clean_label:
        merged_aliases.add(clean_label)
    for mid in canonical_ids:
        merged_aliases.add(mid)
    rec["aliases"] = "|".join(
        prune_metabolite_aliases(
            label=normalize_metabolite_alias_text(rec.get("primary_label") or ""),
            aliases=sorted(merged_aliases),
            ids=canonical_ids,
        )
    )

    merged_ids = normalize_metabolite_ids(
        [
            *parse_pipe_list(rec.get("hmdb_ids") or ""),
            *parse_pipe_list(rec.get("chebi_ids") or ""),
            *parse_pipe_list(rec.get("kegg_ids") or ""),
            *canonical_ids,
        ]
    )
    rec["hmdb_ids"] = "|".join([mid for mid in merged_ids if metabolite_id_bucket(mid) == "hmdb"])
    rec["chebi_ids"] = "|".join([mid for mid in merged_ids if metabolite_id_bucket(mid) == "chebi"])
    rec["kegg_ids"] = "|".join([mid for mid in merged_ids if metabolite_id_bucket(mid) == "kegg"])
    for mid in merged_ids:
        id_to_concept[mid] = concept_id
    if source_clean:
        rec["source"] = "|".join(sorted(set(parse_pipe_list(rec.get("source") or "")) | {source_clean}))


def parse_hmdb_xml_into_records(path: Path, records: dict[str, dict[str, str]], id_to_concept: dict[str, str]) -> int:
    count = 0
    with open_text_with_optional_gzip(path) as handle:
        context = ET.iterparse(handle, events=("end",))
        for _event, elem in context:
            if xml_local_name(elem.tag) != "metabolite":
                continue

            accession = ""
            name = ""
            aliases: list[str] = []
            ids: list[str] = []
            for child in list(elem):
                child_name = xml_local_name(child.tag)
                child_text = normalize(child.text or "")
                if child_name == "accession":
                    accession = child_text
                    if accession:
                        ids.append(accession)
                    continue
                if child_name == "name":
                    name = child_text
                    if name:
                        aliases.append(name)
                    continue
                if child_name in {"chebi_id", "kegg_id"} and child_text:
                    ids.append(child_text)
                    continue
                if child_name == "synonyms":
                    for syn in list(child):
                        if xml_local_name(syn.tag) != "synonym":
                            continue
                        syn_text = normalize(syn.text or "")
                        if syn_text:
                            aliases.append(syn_text)

            if accession or name or aliases:
                upsert_metabolite_builder_record(
                    records,
                    id_to_concept,
                    concept_id=accession,
                    label=name,
                    aliases=aliases,
                    ids=ids,
                    source="hmdb",
                )
                count += 1
            elem.clear()
    return count


def parse_hmdb_sdf_into_records(path: Path, records: dict[str, dict[str, str]], id_to_concept: dict[str, str]) -> int:
    count = 0
    fields: dict[str, str] = {}
    current_field = ""
    buffer: list[str] = []

    def flush_field() -> None:
        nonlocal current_field, buffer
        if not current_field:
            return
        value = normalize(" ".join(buffer))
        if value:
            prior = normalize(fields.get(current_field, ""))
            if prior:
                fields[current_field] = normalize(prior + "; " + value)
            else:
                fields[current_field] = value
        current_field = ""
        buffer = []

    def flush_record() -> None:
        nonlocal count
        flush_field()
        if not fields:
            return

        hmdb_id = canonical_metabolite_id(fields.get("HMDB_ID", "") or fields.get("DATABASE_ID", ""))
        name = normalize_metabolite_alias_text(fields.get("GENERIC_NAME", ""))
        synonyms_raw = fields.get("SYNONYMS", "")
        aliases: list[str] = []
        if synonyms_raw:
            for part in re.split(r"[;|]", synonyms_raw):
                alias = normalize_metabolite_alias_text(part)
                if alias:
                    aliases.append(alias)

        ids: list[str] = []
        if hmdb_id:
            ids.append(hmdb_id)
        # Some HMDB SDF exports include extra cross-reference fields.
        for key in (
            "CHEBI_ID",
            "CHEBI",
            "KEGG_ID",
            "KEGG",
            "KEGG_COMPOUND_ID",
            "COMPOUND_ID",
        ):
            raw_val = fields.get(key, "")
            if not raw_val:
                continue
            for token in re.split(r"[;|,\\s]+", raw_val):
                canonical = canonical_metabolite_id(token)
                if canonical:
                    ids.append(canonical)

        if hmdb_id or name or aliases:
            upsert_metabolite_builder_record(
                records,
                id_to_concept,
                concept_id=hmdb_id or "",
                label=name,
                aliases=aliases,
                ids=ids,
                source="hmdb-sdf",
            )
            count += 1
        fields.clear()

    with open_text_with_optional_gzip(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if line == "$$$$":
                flush_record()
                continue
            if line.startswith("> <") and line.endswith(">"):
                flush_field()
                current_field = normalize(line[3:-1].strip("<>")).upper()
                buffer = []
                continue
            if current_field:
                if line.strip() == "":
                    flush_field()
                else:
                    buffer.append(line)

    # Last record may not end with $$$$.
    flush_record()
    return count


def parse_hmdb_source_into_records(path: Path, records: dict[str, dict[str, str]], id_to_concept: dict[str, str]) -> int:
    name = path.name.lower()
    if name.endswith(".sdf") or name.endswith(".sdf.gz"):
        return parse_hmdb_sdf_into_records(path, records, id_to_concept)
    return parse_hmdb_xml_into_records(path, records, id_to_concept)


def parse_chebi_names_into_records(path: Path, records: dict[str, dict[str, str]], id_to_concept: dict[str, str]) -> int:
    count = 0
    with open_text_with_optional_gzip(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return 0
        field_map = {norm_key(name): name for name in reader.fieldnames if name}
        id_field = (
            field_map.get("chebi_id")
            or field_map.get("compound_id")
            or field_map.get("id")
        )
        name_field = field_map.get("name")
        type_field = field_map.get("type")
        status_field = field_map.get("status_id") or field_map.get("status")
        ascii_name_field = field_map.get("ascii_name")
        if not id_field or not name_field:
            return 0
        grouped: dict[str, list[tuple[str, str, str]]] = {}
        for row in reader:
            chebi_raw = normalize(row.get(id_field) or "")
            if not chebi_raw:
                continue
            chebi_id = canonical_metabolite_id(f"CHEBI:{chebi_raw}")
            if not chebi_id:
                continue
            name = normalize_metabolite_alias_text(
                row.get(name_field) or "",
                ascii_fallback=row.get(ascii_name_field) if ascii_name_field else "",
            )
            if not name:
                continue
            row_type = normalize(row.get(type_field) or "") if type_field else ""
            status_id = normalize(row.get(status_field) or "") if status_field else ""
            grouped.setdefault(chebi_id, []).append((name, row_type, status_id))
            count += 1

        for chebi_id, names in grouped.items():
            label = choose_chebi_primary_label(names)
            aliases = sorted({name for name, _row_type, _status_id in names if keep_metabolite_alias(name)})
            upsert_metabolite_builder_record(
                records,
                id_to_concept,
                concept_id=chebi_id,
                label=label,
                aliases=aliases,
                ids=[chebi_id],
                source="chebi",
            )
    return count


def parse_kegg_conv_into_records(path: Path, records: dict[str, dict[str, str]], id_to_concept: dict[str, str]) -> int:
    count = 0
    with open_text_with_optional_gzip(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line or "\t" not in line:
                continue
            left, right = line.split("\t", 1)
            left = normalize(left)
            right = normalize(right)
            left_id = canonical_metabolite_id(left)
            right_id = canonical_metabolite_id(right)
            if not left_id and not right_id:
                continue
            concept_id = left_id or right_id
            ids = [mid for mid in [left_id, right_id] if mid]
            upsert_metabolite_builder_record(
                records,
                id_to_concept,
                concept_id=concept_id,
                label="",
                aliases=[],
                ids=ids,
                source="kegg-conv",
            )
            count += 1
    return count


def build_metabolite_alias_records(
    *,
    output_path: Path,
    hmdb_xml_path: Path | None,
    chebi_names_path: Path | None,
    kegg_conv_path: Path | None,
    merge_existing: bool,
) -> tuple[int, int, int, int]:
    records = load_metabolite_alias_records(output_path) if merge_existing else {}
    id_to_concept: dict[str, str] = {}
    for concept_id, row in records.items():
        for mid in normalize_metabolite_ids(
            [
                concept_id,
                *parse_pipe_list(row.get("hmdb_ids") or ""),
                *parse_pipe_list(row.get("chebi_ids") or ""),
                *parse_pipe_list(row.get("kegg_ids") or ""),
            ]
        ):
            id_to_concept[mid] = concept_id
    hmdb_rows = 0
    chebi_rows = 0
    kegg_rows = 0

    if hmdb_xml_path is not None:
        hmdb_rows = parse_hmdb_source_into_records(hmdb_xml_path, records, id_to_concept)
    if chebi_names_path is not None:
        chebi_rows = parse_chebi_names_into_records(chebi_names_path, records, id_to_concept)
    if kegg_conv_path is not None:
        kegg_rows = parse_kegg_conv_into_records(kegg_conv_path, records, id_to_concept)

    write_metabolite_alias_records(output_path, records)
    return hmdb_rows, chebi_rows, kegg_rows, len(records)


def fetch_json(url: str, timeout: float) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": HTTP_USER_AGENT,
            "Accept": "application/json",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:  # nosec B310
        return json.loads(response.read().decode("utf-8"))


def parse_uniprot_entry_aliases(entry: dict[str, Any]) -> tuple[list[str], str, list[str]]:
    aliases: set[str] = set()
    gene_primary = ""
    gene_ids: set[str] = set()

    mnemonic = normalize(entry.get("uniProtkbId") or "")
    if mnemonic:
        aliases.add(mnemonic)

    def add_name_block(name_block: Any) -> None:
        if not isinstance(name_block, dict):
            return
        full = name_block.get("fullName") or {}
        if isinstance(full, dict):
            value = normalize(full.get("value") or "")
            if value:
                aliases.add(value)
        short_names = name_block.get("shortNames") or []
        if isinstance(short_names, list):
            for short in short_names:
                if isinstance(short, dict):
                    value = normalize(short.get("value") or "")
                    if value:
                        aliases.add(value)
        ec_numbers = name_block.get("ecNumbers") or []
        if isinstance(ec_numbers, list):
            for ec in ec_numbers:
                if isinstance(ec, dict):
                    value = normalize(ec.get("value") or "")
                    if value:
                        aliases.add(value)

    genes = entry.get("genes") or []
    if isinstance(genes, list):
        for gene in genes:
            if not isinstance(gene, dict):
                continue
            gene_name = ((gene.get("geneName") or {}).get("value")) if isinstance(gene.get("geneName"), dict) else ""
            if gene_name:
                clean = normalize(gene_name)
                aliases.add(clean)
                if not gene_primary:
                    gene_primary = clean
            syns = gene.get("synonyms") or []
            if isinstance(syns, list):
                for syn in syns:
                    if isinstance(syn, dict):
                        value = normalize(syn.get("value") or "")
                        if value:
                            aliases.add(value)
            orf_names = gene.get("orfNames") or []
            if isinstance(orf_names, list):
                for orf_name in orf_names:
                    if isinstance(orf_name, dict):
                        value = normalize(orf_name.get("value") or "")
                        if value:
                            aliases.add(value)
            locus_names = gene.get("orderedLocusNames") or []
            if isinstance(locus_names, list):
                for locus_name in locus_names:
                    if isinstance(locus_name, dict):
                        value = normalize(locus_name.get("value") or "")
                        if value:
                            aliases.add(value)

    cross_refs = entry.get("uniProtKBCrossReferences") or entry.get("dbReferences") or []
    if isinstance(cross_refs, list):
        for ref in cross_refs:
            if not isinstance(ref, dict):
                continue
            db = norm_key(ref.get("database") or ref.get("type") or "")
            id_values: list[str] = []
            ref_id = normalize(ref.get("id") or "")
            if ref_id:
                id_values.append(ref_id)
            properties = ref.get("properties") or []
            if isinstance(properties, list):
                for prop in properties:
                    if not isinstance(prop, dict):
                        continue
                    value = normalize(prop.get("value") or "")
                    if value:
                        id_values.append(value)
            if db == "ensembl":
                for raw in id_values:
                    for gid in gene_id_alias_variants(raw, include_unprefixed_numeric=False):
                        if gid.startswith("ENSG"):
                            gene_ids.add(gid)
            elif db in {"geneid", "entrez gene", "ncbigene"}:
                for raw in id_values:
                    for gid in gene_id_alias_variants(raw, include_unprefixed_numeric=False):
                        if gid.startswith(("NCBIGENE:", "GENEID:", "ENTREZ:")):
                            gene_ids.add(gid)

    protein_desc = entry.get("proteinDescription") or {}
    if isinstance(protein_desc, dict):
        add_name_block(protein_desc.get("recommendedName") or {})
        alt_names = protein_desc.get("alternativeNames") or []
        if isinstance(alt_names, list):
            for alt in alt_names:
                add_name_block(alt)
        sub_names = protein_desc.get("submissionNames") or []
        if isinstance(sub_names, list):
            for sub in sub_names:
                add_name_block(sub)
        for section_key in ("contains", "includes"):
            section_rows = protein_desc.get(section_key) or []
            if not isinstance(section_rows, list):
                continue
            for section in section_rows:
                if not isinstance(section, dict):
                    continue
                add_name_block(section.get("recommendedName") or {})
                section_alt_names = section.get("alternativeNames") or []
                if isinstance(section_alt_names, list):
                    for alt in section_alt_names:
                        add_name_block(alt)

    return sorted(aliases), gene_primary, sorted(gene_ids)


def fetch_uniprot_alias_payloads(
    accessions: list[str],
    *,
    timeout: float,
    workers: int,
) -> list[tuple[str, list[str], str, list[str], bool]]:
    def fetch_one(acc: str) -> tuple[str, list[str], str, list[str], bool]:
        try:
            entry = fetch_json(UNIPROT_ENTRY_URL.format(acc=acc), timeout=timeout)
            aliases, gene, gene_ids = parse_uniprot_entry_aliases(entry)
            if not aliases and not gene and not gene_ids:
                return acc, [], "", [], False
            return acc, aliases, gene, gene_ids, True
        except Exception:
            return acc, [], "", [], False

    if workers <= 1:
        return [fetch_one(acc) for acc in accessions]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(fetch_one, accessions))


def enrich_uniprot_aliases_for_queries(
    queries: list[str],
    alias_path: Path,
    timeout: float,
    workers: int,
    source_tag: str,
) -> tuple[int, int]:
    records = load_uniprot_alias_records(alias_path)
    accession_set: set[str] = set()
    for query in queries:
        for variant in query_variants(query):
            canonical = canonical_accession(variant)
            if canonical:
                accession_set.add(canonical)
    accessions = sorted(accession_set)
    missing = [acc for acc in accessions if acc not in records]
    if not missing:
        return 0, 0

    added = 0
    failed = 0
    for acc, aliases, gene, gene_ids, ok in fetch_uniprot_alias_payloads(missing, timeout=timeout, workers=workers):
        if not ok:
            failed += 1
            continue
        records[acc] = {
            "accession": acc,
            "aliases": "|".join(sorted(set(aliases))),
            "gene": normalize(gene),
            "gene_ids": "|".join(
                sorted(
                    {
                        gid
                        for raw in gene_ids
                        for gid in gene_id_alias_variants(raw, include_unprefixed_numeric=False)
                    }
                )
            ),
            "source": source_tag,
        }
        added += 1

    if added > 0:
        write_uniprot_alias_records(alias_path, records)
    return added, failed


def backfill_uniprot_gene_ids(
    *,
    alias_path: Path,
    timeout: float,
    workers: int,
    source_tag: str,
    refresh_all: bool = False,
) -> tuple[int, int, int]:
    records = load_uniprot_alias_records(alias_path)
    if not records:
        return 0, 0, 0

    targets: list[str] = []
    for acc in sorted(records):
        has_gene_ids = bool(parse_pipe_list(records[acc].get("gene_ids") or ""))
        if refresh_all or not has_gene_ids:
            targets.append(acc)
    if not targets:
        return 0, 0, 0

    updated = 0
    failed = 0
    for acc, aliases, gene, gene_ids, ok in fetch_uniprot_alias_payloads(targets, timeout=timeout, workers=workers):
        if not ok:
            failed += 1
            continue
        record = records.get(acc) or {"accession": acc, "aliases": "", "gene": "", "gene_ids": "", "source": ""}
        before = (
            normalize(record.get("aliases") or ""),
            normalize(record.get("gene") or ""),
            normalize(record.get("gene_ids") or ""),
            normalize(record.get("source") or ""),
        )
        merged_aliases = sorted(
            set(parse_pipe_list(record.get("aliases") or "")) | {normalize(a) for a in aliases if normalize(a)}
        )
        merged_gene_ids = sorted(
            set(parse_pipe_list(record.get("gene_ids") or ""))
            | {
                gid
                for raw in gene_ids
                for gid in gene_id_alias_variants(raw, include_unprefixed_numeric=False)
            }
        )
        merged_gene = normalize(record.get("gene") or "") or normalize(gene)
        merged_sources = sorted(set(parse_pipe_list(record.get("source") or "")) | ({source_tag} if source_tag else set()))

        records[acc] = {
            "accession": acc,
            "aliases": "|".join(merged_aliases),
            "gene": merged_gene,
            "gene_ids": "|".join(merged_gene_ids),
            "source": "|".join(merged_sources),
        }
        after = (
            normalize(records[acc].get("aliases") or ""),
            normalize(records[acc].get("gene") or ""),
            normalize(records[acc].get("gene_ids") or ""),
            normalize(records[acc].get("source") or ""),
        )
        if after != before:
            updated += 1

    if updated > 0:
        write_uniprot_alias_records(alias_path, records)
    return updated, failed, len(targets)


def build_index(
    term_rows: list[dict[str, Any]],
    analyte_rows: list[dict[str, Any]],
    accession_aliases: dict[str, list[str]],
    metabolite_alias_records: dict[str, dict[str, str]] | None = None,
) -> dict[str, Any]:
    term_meta: dict[str, dict[str, Any]] = {}
    exact_term_index: dict[str, set[str]] = {}
    token_index: dict[str, set[str]] = {}
    term_metabolite_id_index: dict[str, set[str]] = {}

    for row in term_rows:
        efo_id = row["efo_id"]
        label = row["label"]
        syns = row["synonyms"]
        metabolite_ids = normalize_metabolite_ids(row.get("metabolite_ids") or [])

        term_meta[efo_id] = {
            "label": label,
            "synonyms": syns,
            "subject": normalize_measurement_subject(label),
            "metabolite_ids": metabolite_ids,
        }
        for met_id in metabolite_ids:
            term_metabolite_id_index.setdefault(met_id, set()).add(efo_id)

        searchable = [label] + syns
        expanded_searchable: set[str] = set()
        for phrase in searchable:
            if not phrase:
                continue
            expanded_searchable.add(phrase)
            # Add normalized phrase variants so EFO/OBA synonym wording can be
            # matched against UniProt aliases with different shorthand.
            for variant in subject_phrase_variants(phrase):
                if variant:
                    expanded_searchable.add(variant)
        for phrase in expanded_searchable:
            key = norm_key(phrase)
            if key:
                exact_term_index.setdefault(key, set()).add(efo_id)
            for tok in tokenize(phrase):
                token_index.setdefault(tok, set()).add(efo_id)
    token_freq: dict[str, int] = {tok: len(efo_ids) for tok, efo_ids in token_index.items()}

    analyte_index: dict[str, list[dict[str, Any]]] = {}
    for row in analyte_rows:
        analyte_index.setdefault(row["analyte_key"], []).append(row)

    alias_accession_index: dict[str, set[str]] = {}
    alias_loose_index: dict[str, set[str]] = {}
    symbol_accession_index: dict[str, set[str]] = {}
    gene_id_accession_index: dict[str, set[str]] = {}
    for acc, aliases in accession_aliases.items():
        alias_accession_index.setdefault(norm_key(acc), set()).add(acc)
        acc_loose = loose_key(acc)
        if acc_loose:
            alias_loose_index.setdefault(acc_loose, set()).add(acc)
        for alias in aliases:
            alias_clean = normalize(alias)
            if not alias_clean:
                continue
            for alias_variant in alias_variants(alias_clean):
                key = norm_key(alias_variant)
                if key:
                    alias_accession_index.setdefault(key, set()).add(acc)
                loose = loose_key(alias_variant)
                if loose:
                    alias_loose_index.setdefault(loose, set()).add(acc)
                symbol_key = symbol_query_key(alias_variant)
                if symbol_key:
                    symbol_accession_index.setdefault(symbol_key, set()).add(acc)
                for gid in gene_id_alias_variants(alias_variant, include_unprefixed_numeric=False):
                    gene_id_accession_index.setdefault(norm_key(gid), set()).add(acc)

    metabolite_concept_index: dict[str, dict[str, Any]] = {}
    metabolite_id_concept_index: dict[str, set[str]] = {}
    metabolite_alias_concept_index: dict[str, set[str]] = {}
    metabolite_alias_loose_index: dict[str, set[str]] = {}
    metabolite_alias_records = metabolite_alias_records or {}
    for concept_id, row in metabolite_alias_records.items():
        concept = normalize(concept_id)
        if not concept:
            continue
        label = normalize_metabolite_alias_text(row.get("primary_label") or row.get("label") or "")
        aliases = [normalize_metabolite_alias_text(x) for x in parse_pipe_list(row.get("aliases") or "")]
        aliases = [x for x in aliases if x and keep_metabolite_alias(x)]
        ids = normalize_metabolite_ids(
            [
                concept,
                *parse_pipe_list(row.get("hmdb_ids") or ""),
                *parse_pipe_list(row.get("chebi_ids") or ""),
                *parse_pipe_list(row.get("kegg_ids") or ""),
                *parse_pipe_list(row.get("ids") or ""),
            ]
        )
        aliases = prune_metabolite_aliases(label=label, aliases=aliases, ids=ids)
        searchable_aliases = sorted(set([label, *aliases, *ids]) - {""})
        # Keep identity profile compact; aliases are already quality-ranked.
        identity_aliases = aliases[:12]
        subject_terms = sorted(extract_subject_terms(label, identity_aliases))
        metabolite_concept_index[concept] = {
            "label": label,
            "aliases": aliases,
            "ids": ids,
            "subject_terms": subject_terms,
            "source": normalize(row.get("source") or ""),
        }
        for met_id in ids:
            metabolite_id_concept_index.setdefault(met_id, set()).add(concept)
        for alias in searchable_aliases:
            for alias_variant in alias_variants(alias):
                key = norm_key(alias_variant)
                if key:
                    metabolite_alias_concept_index.setdefault(key, set()).add(concept)
                loose = loose_key(alias_variant)
                if loose:
                    metabolite_alias_loose_index.setdefault(loose, set()).add(concept)

    return {
        "version": "5",
        "built_at": datetime.now(UTC).isoformat(),
        "term_meta": term_meta,
        "exact_term_index": {k: sorted(v) for k, v in exact_term_index.items()},
        "token_index": {k: sorted(v) for k, v in token_index.items()},
        "token_freq": token_freq,
        "term_metabolite_id_index": {k: sorted(v) for k, v in term_metabolite_id_index.items()},
        "analyte_index": analyte_index,
        "accession_alias_index": accession_aliases,
        "alias_accession_index": {k: sorted(v) for k, v in alias_accession_index.items()},
        "alias_loose_index": {k: sorted(v) for k, v in alias_loose_index.items()},
        "symbol_accession_index": {k: sorted(v) for k, v in symbol_accession_index.items()},
        "gene_id_accession_index": {k: sorted(v) for k, v in gene_id_accession_index.items()},
        "metabolite_concept_index": metabolite_concept_index,
        "metabolite_id_concept_index": {k: sorted(v) for k, v in metabolite_id_concept_index.items()},
        "metabolite_alias_concept_index": {k: sorted(v) for k, v in metabolite_alias_concept_index.items()},
        "metabolite_alias_loose_index": {k: sorted(v) for k, v in metabolite_alias_loose_index.items()},
    }


def save_index(index: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(index, ensure_ascii=True, sort_keys=True), encoding="utf-8")


def load_index(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Index file not found: {path}")
    index = json.loads(path.read_text(encoding="utf-8"))
    term_meta = index.get("term_meta") or {}
    if isinstance(term_meta, dict):
        for meta in term_meta.values():
            if not isinstance(meta, dict):
                continue
            if "subject" not in meta:
                meta["subject"] = normalize_measurement_subject(normalize(meta.get("label") or ""))
            metabolite_ids = normalize_metabolite_ids(meta.get("metabolite_ids") or [])
            meta["metabolite_ids"] = metabolite_ids

    if "term_metabolite_id_index" not in index:
        rebuilt_term_metabolite_id_index: dict[str, set[str]] = {}
        if isinstance(term_meta, dict):
            for efo_id, meta in term_meta.items():
                if not isinstance(meta, dict):
                    continue
                for met_id in normalize_metabolite_ids(meta.get("metabolite_ids") or []):
                    rebuilt_term_metabolite_id_index.setdefault(met_id, set()).add(normalize_id(efo_id))
        index["term_metabolite_id_index"] = {k: sorted(v) for k, v in rebuilt_term_metabolite_id_index.items()}
    if "token_freq" not in index:
        token_index = index.get("token_index") or {}
        if isinstance(token_index, dict):
            index["token_freq"] = {
                tok: len(efo_ids) if isinstance(efo_ids, list) else 0
                for tok, efo_ids in token_index.items()
            }
        else:
            index["token_freq"] = {}

    if "alias_accession_index" not in index:
        alias_accession_index: dict[str, set[str]] = {}
        alias_loose_index: dict[str, set[str]] = {}
        symbol_accession_index: dict[str, set[str]] = {}
        gene_id_accession_index: dict[str, set[str]] = {}
        accession_aliases = index.get("accession_alias_index") or {}
        if isinstance(accession_aliases, dict):
            for acc, aliases in accession_aliases.items():
                acc_clean = normalize(acc).upper()
                if not acc_clean:
                    continue
                alias_accession_index.setdefault(norm_key(acc_clean), set()).add(acc_clean)
                acc_loose = loose_key(acc_clean)
                if acc_loose:
                    alias_loose_index.setdefault(acc_loose, set()).add(acc_clean)
                if isinstance(aliases, list):
                    for alias in aliases:
                        alias_clean = normalize(alias)
                        if not alias_clean:
                            continue
                        for alias_variant in alias_variants(alias_clean):
                            key = norm_key(alias_variant)
                            if key:
                                alias_accession_index.setdefault(key, set()).add(acc_clean)
                            loose = loose_key(alias_variant)
                            if loose:
                                alias_loose_index.setdefault(loose, set()).add(acc_clean)
                            symbol_key = symbol_query_key(alias_variant)
                            if symbol_key:
                                symbol_accession_index.setdefault(symbol_key, set()).add(acc_clean)
                            for gid in gene_id_alias_variants(alias_variant, include_unprefixed_numeric=False):
                                gene_id_accession_index.setdefault(norm_key(gid), set()).add(acc_clean)
        index["alias_accession_index"] = {k: sorted(v) for k, v in alias_accession_index.items()}
        index["alias_loose_index"] = {k: sorted(v) for k, v in alias_loose_index.items()}
        index["symbol_accession_index"] = {k: sorted(v) for k, v in symbol_accession_index.items()}
        index["gene_id_accession_index"] = {k: sorted(v) for k, v in gene_id_accession_index.items()}
    elif (
        "symbol_accession_index" not in index
        or "alias_loose_index" not in index
        or "gene_id_accession_index" not in index
    ):
        alias_loose_index: dict[str, set[str]] = {}
        symbol_accession_index: dict[str, set[str]] = {}
        gene_id_accession_index: dict[str, set[str]] = {}
        accession_aliases = index.get("accession_alias_index") or {}
        if isinstance(accession_aliases, dict):
            for acc, aliases in accession_aliases.items():
                acc_clean = normalize(acc).upper()
                if not acc_clean:
                    continue
                acc_loose = loose_key(acc_clean)
                if acc_loose:
                    alias_loose_index.setdefault(acc_loose, set()).add(acc_clean)
                if isinstance(aliases, list):
                    for alias in aliases:
                        alias_clean = normalize(alias)
                        if not alias_clean:
                            continue
                        for alias_variant in alias_variants(alias_clean):
                            loose = loose_key(alias_variant)
                            if loose:
                                alias_loose_index.setdefault(loose, set()).add(acc_clean)
                            symbol_key = symbol_query_key(alias_variant)
                            if symbol_key:
                                symbol_accession_index.setdefault(symbol_key, set()).add(acc_clean)
                            for gid in gene_id_alias_variants(alias_variant, include_unprefixed_numeric=False):
                                gene_id_accession_index.setdefault(norm_key(gid), set()).add(acc_clean)
        index["alias_loose_index"] = {k: sorted(v) for k, v in alias_loose_index.items()}
        index["symbol_accession_index"] = {k: sorted(v) for k, v in symbol_accession_index.items()}
        index["gene_id_accession_index"] = {k: sorted(v) for k, v in gene_id_accession_index.items()}

    if "metabolite_concept_index" not in index:
        index["metabolite_concept_index"] = {}
    if (
        "metabolite_id_concept_index" not in index
        or "metabolite_alias_concept_index" not in index
        or "metabolite_alias_loose_index" not in index
    ):
        concept_index = index.get("metabolite_concept_index") or {}
        metabolite_id_concept_index: dict[str, set[str]] = {}
        metabolite_alias_concept_index: dict[str, set[str]] = {}
        metabolite_alias_loose_index: dict[str, set[str]] = {}
        if isinstance(concept_index, dict):
            for concept_id, concept_meta in concept_index.items():
                concept = normalize(concept_id)
                if not concept or not isinstance(concept_meta, dict):
                    continue
                label = normalize_metabolite_alias_text(concept_meta.get("label") or "")
                aliases = [
                    normalize_metabolite_alias_text(x)
                    for x in (concept_meta.get("aliases") or [])
                    if normalize_metabolite_alias_text(x)
                ]
                ids = normalize_metabolite_ids([*aliases, *(concept_meta.get("ids") or []), concept])
                aliases = prune_metabolite_aliases(label=label, aliases=aliases, ids=ids)
                searchable = sorted(set([label, *aliases, *ids]) - {""})
                subject_terms = extract_subject_terms(label, aliases[:12])
                concept_meta["subject_terms"] = sorted(subject_terms)
                concept_meta["ids"] = ids
                concept_meta["aliases"] = aliases
                for mid in ids:
                    metabolite_id_concept_index.setdefault(mid, set()).add(concept)
                for alias in searchable:
                    for alias_variant in alias_variants(alias):
                        key = norm_key(alias_variant)
                        if key:
                            metabolite_alias_concept_index.setdefault(key, set()).add(concept)
                        loose = loose_key(alias_variant)
                        if loose:
                            metabolite_alias_loose_index.setdefault(loose, set()).add(concept)
        index["metabolite_id_concept_index"] = {k: sorted(v) for k, v in metabolite_id_concept_index.items()}
        index["metabolite_alias_concept_index"] = {k: sorted(v) for k, v in metabolite_alias_concept_index.items()}
        index["metabolite_alias_loose_index"] = {k: sorted(v) for k, v in metabolite_alias_loose_index.items()}
    return index


PROCESS_INDEX: dict[str, Any] | None = None
PROCESS_MAP_KW: dict[str, Any] = {}


def init_process_mapper(
    index_path: str,
    top_k: int,
    min_score: float,
    matrix_priority: list[str],
    measurement_context: str,
    additional_contexts: list[str],
    additional_context_keywords: list[str],
    name_mode: str,
    force_map_best: bool,
    fallback_efo_id: str,
    entity_type: str,
) -> None:
    global PROCESS_INDEX, PROCESS_MAP_KW
    PROCESS_INDEX = load_index(Path(index_path))
    PROCESS_MAP_KW = {
        "top_k": top_k,
        "min_score": min_score,
        "matrix_priority": matrix_priority,
        "measurement_context": measurement_context,
        "additional_contexts": additional_contexts,
        "additional_context_keywords": additional_context_keywords,
        "name_mode": name_mode,
        "force_map_best": force_map_best,
        "fallback_efo_id": fallback_efo_id,
        "entity_type": entity_type,
    }


def process_map_one(query_item: tuple[str, str]) -> list[dict[str, str]]:
    if PROCESS_INDEX is None:
        raise RuntimeError("Process mapper index is not initialized")
    query, input_type = query_item
    return map_one(
        query,
        PROCESS_INDEX,
        top_k=int(PROCESS_MAP_KW["top_k"]),
        min_score=float(PROCESS_MAP_KW["min_score"]),
        matrix_priority=list(PROCESS_MAP_KW["matrix_priority"]),
        measurement_context=str(PROCESS_MAP_KW["measurement_context"]),
        additional_contexts=list(PROCESS_MAP_KW["additional_contexts"]),
        additional_context_keywords=list(PROCESS_MAP_KW["additional_context_keywords"]),
        name_mode=effective_name_mode(str(PROCESS_MAP_KW["name_mode"]), input_type),
        force_map_best=bool(PROCESS_MAP_KW["force_map_best"]),
        fallback_efo_id=str(PROCESS_MAP_KW["fallback_efo_id"]),
        input_type=input_type,
        entity_type=str(PROCESS_MAP_KW.get("entity_type", "auto")),
    )


def merge_candidate(store: dict[str, Candidate], cand: Candidate) -> None:
    current = store.get(cand.efo_id)
    if current is None:
        store[cand.efo_id] = cand
        return
    if cand.is_validated and not current.is_validated:
        store[cand.efo_id] = cand
        return
    if cand.is_validated == current.is_validated:
        current_is_token = current.matched_via.startswith("local-term-token")
        cand_is_token = cand.matched_via.startswith("local-term-token")
        # Keep explicit/exact match evidence over token-retrieved evidence for
        # the same ontology ID. Token paths are noisier and can otherwise
        # downgrade a stable exact hit into a review-only candidate.
        if current_is_token and not cand_is_token:
            store[cand.efo_id] = cand
            return
        if cand_is_token and not current_is_token:
            return
    if cand.score > current.score:
        store[cand.efo_id] = cand
        return
    if cand.score == current.score:
        cand_exact = "subject-exact" in cand.matched_via
        cur_exact = "subject-exact" in current.matched_via
        if cand_exact and not cur_exact:
            store[cand.efo_id] = cand


def parse_matrix_priority(raw: str) -> list[str]:
    parts = [norm_key(x) for x in raw.split(",")]
    return [p for p in parts if p]


def parse_measurement_context(raw: str) -> str:
    key = norm_key(raw).replace("_", " ")
    return MEASUREMENT_CONTEXT_ALIASES.get(key, DEFAULT_MEASUREMENT_CONTEXT)


def parse_additional_contexts(raw: str) -> list[str]:
    contexts: list[str] = []
    seen: set[str] = set()
    for token in raw.split(","):
        key = norm_key(token).replace("_", " ")
        if not key:
            continue
        ctx = MEASUREMENT_CONTEXT_ALIASES.get(key, "")
        if not ctx or ctx == "auto":
            continue
        if ctx in seen:
            continue
        seen.add(ctx)
        contexts.append(ctx)
    return contexts


def parse_additional_context_keywords(raw: str) -> list[str]:
    keywords: list[str] = []
    seen: set[str] = set()
    for token in raw.split(","):
        key = normalize(token)
        if not key:
            continue
        lowered = norm_key(key)
        if lowered in seen:
            continue
        seen.add(lowered)
        keywords.append(key)
    return keywords


def contains_context_token(text: str, token: str) -> bool:
    token = norm_key(token)
    if not token:
        return False
    if " " in token:
        return token in text
    pattern = r"(?<![a-z0-9])" + re.escape(token) + r"(?![a-z0-9])"
    return re.search(pattern, text) is not None


def extract_measurement_context_tags(label: str, synonyms: list[str]) -> set[str]:
    tags: set[str] = set()
    text = norm_key(" ".join([label] + synonyms))
    for context, tokens in MEASUREMENT_CONTEXT_PATTERNS.items():
        for token in tokens:
            if contains_context_token(text, token):
                tags.add(context)
                break
    return tags


def extract_explicit_matrix_tags(label: str, synonyms: list[str]) -> set[str]:
    text = norm_key(" ".join([label] + synonyms))
    tags: set[str] = set()
    if re.search(r"\bin (?:blood serum|serum)\b", text):
        tags.add("serum")
    if re.search(r"\bin (?:blood plasma|plasma)\b", text):
        tags.add("plasma")
    if re.search(r"\bin blood\b(?!\s+(?:serum|plasma))", text):
        tags.add("blood")
    if re.search(r"\bin cerebrospinal fluid\b|\bcsf\b", text):
        tags.add("cerebrospinal_fluid")
    if re.search(r"\bin (?:urine|urinary)\b", text):
        tags.add("urine")
    if re.search(r"\bin (?:saliva|salivary)\b", text):
        tags.add("saliva")
    if re.search(r"\bin tissue\b", text):
        tags.add("tissue")
    return tags


def expand_context_scope(context: str) -> set[str]:
    if context == "blood":
        # Default blood context allows blood/plasma/circulating but excludes serum.
        return {"blood", "plasma", "circulating"}
    if context in {"plasma", "serum", "cerebrospinal_fluid", "urine", "saliva", "tissue"}:
        return {context}
    return set()


def allowed_context_tags(measurement_context: str, additional_contexts: list[str]) -> set[str]:
    allowed = expand_context_scope(measurement_context)
    for ctx in additional_contexts:
        allowed.update(expand_context_scope(ctx))
    return allowed


def has_additional_context_keyword_match(label: str, synonyms: list[str], context_keywords: list[str]) -> bool:
    if not context_keywords:
        return False
    text = norm_key(" ".join([label] + synonyms))
    for keyword in context_keywords:
        if contains_context_token(text, keyword):
            return True
    return False


def context_compatible(
    label: str,
    synonyms: list[str],
    measurement_context: str,
    additional_contexts: list[str] | tuple[str, ...] = (),
    additional_context_keywords: list[str] | tuple[str, ...] = (),
) -> bool:
    keyword_list = list(additional_context_keywords)
    if keyword_list and measurement_context == "auto":
        return has_additional_context_keyword_match(label, synonyms, keyword_list)
    if measurement_context == "auto":
        return True

    if keyword_list and has_additional_context_keyword_match(label, synonyms, keyword_list):
        return True

    allowed = allowed_context_tags(measurement_context, list(additional_contexts))
    if not allowed:
        return True
    explicit_tags = extract_explicit_matrix_tags(label, synonyms)
    if explicit_tags:
        return bool(explicit_tags & allowed)
    tags = extract_measurement_context_tags(label, synonyms)
    if not tags:
        return True
    return bool(tags & allowed)


def matrix_bonus(label: str, synonyms: list[str], matrix_priority: list[str]) -> float:
    if not matrix_priority:
        return 0.0
    matrix = measurement_matrix_bucket(label)
    if not matrix:
        return 0.0
    # Serum-specific measurements are only a fallback when alternatives exist.
    if matrix == "serum":
        return -0.25
    for idx, preferred in enumerate(matrix_priority):
        if preferred == matrix:
            # Ordered preference: first entry gets the largest bonus.
            return max(0.0, 0.06 - 0.02 * idx)
    return 0.0


def matrix_preference_rank(label: str, matrix_priority: list[str]) -> int:
    if not matrix_priority:
        return -1
    matrix = measurement_matrix_bucket(label)
    if not matrix:
        # Prefer matrix-agnostic terms over serum-only terms.
        return 0
    if matrix == "serum":
        return -2
    for idx, preferred in enumerate(matrix_priority):
        if preferred == matrix:
            return len(matrix_priority) - idx
    return 0


def measurement_matrix_bucket(label: str) -> str:
    text = norm_key(label)
    if not text:
        return ""
    if contains_context_token(text, "plasma"):
        return "plasma"
    if contains_context_token(text, "serum"):
        return "serum"
    if contains_context_token(text, "blood") or contains_context_token(text, "circulating"):
        return "blood"
    return ""


def measurement_label_rank(label: str) -> int:
    text = norm_key(label)
    if text.startswith("level of "):
        return 3
    if text.endswith(" measurement") or " measurement " in text:
        return 2
    if text.startswith("amount of "):
        return 1
    return 0


def is_ratio_trait(label: str, synonyms: list[str]) -> bool:
    text = norm_key(" ".join([label] + synonyms))
    if "ratio" in text or "quotient" in text:
        return True
    # Treat slash as ratio only for explicit "A / B" style patterns.
    # This avoids false positives for names like "paraoxonase/lactonase".
    if re.search(r"\b[a-z0-9][a-z0-9\-]*\s+/\s+[a-z0-9][a-z0-9\-]*\b", text):
        return True
    return False


def query_requests_ratio(query: str) -> bool:
    q = norm_key(query)
    if "ratio" in q or "quotient" in q:
        return True
    slash_match = re.search(r"\b([a-z0-9][a-z0-9\-]*)\s*/\s*([a-z0-9][a-z0-9\-]*)\b", q)
    if not slash_match:
        return False
    left, right = slash_match.group(1), slash_match.group(2)
    # Permit compact biochemical shorthand (for example MUFA/PUFA), while
    # avoiding broad slash-name false positives such as paraoxonase/lactonase.
    if len(left) <= 6 and len(right) <= 6:
        return True
    if any(token in {left, right} for token in {"total", "hdl", "ldl", "vldl", "idl"}):
        return True
    if "%" in q:
        return True
    return False


@lru_cache(maxsize=100_000)
def alias_boundary_pattern(alias_lower: str) -> re.Pattern[str]:
    return re.compile(r"(?<![A-Za-z0-9])" + re.escape(alias_lower) + r"(?![A-Za-z0-9])")


def contains_term(candidate_text_lower: str, alias: str) -> bool:
    # Keep strict token boundaries to avoid partial symbol matches.
    alias_lower = alias.lower()
    if not alias_lower:
        return False
    return alias_boundary_pattern(alias_lower).search(candidate_text_lower) is not None


@lru_cache(maxsize=200_000)
def alias_variants(alias: str) -> tuple[str, ...]:
    base = normalize(alias)
    if not base:
        return ()

    variants: set[str] = {base}
    # UniProt names sometimes include cleavage annotations in square brackets.
    # Keep the pre-cleavage canonical name as a clean alias.
    pre_cleaved = normalize(re.sub(r"\s*\[cleaved into:.*$", "", base, flags=re.IGNORECASE))
    if pre_cleaved:
        variants.add(pre_cleaved)
    # Strip UniProt annotation tails (for example "[Includes: ...]"), including
    # malformed/unclosed tails that occasionally appear in bulk exports.
    de_annotated = normalize(
        re.sub(r"\s*\[(?:includes|contains|cleaved into):.*$", "", pre_cleaved or base, flags=re.IGNORECASE)
    )
    if de_annotated:
        variants.add(de_annotated)
    # UniProt aliases sometimes append context qualifiers after a comma
    # (for example ", muscle-specific form"). Keep the core analyte phrase too.
    for phrase in [base, pre_cleaved, de_annotated]:
        if not phrase or "," not in phrase:
            continue
        lead = normalize(phrase.split(",", 1)[0])
        if lead and len(tokenize(norm_key(lead))) >= 3:
            variants.add(lead)

    # Add parenthetical expansions and stripped base form to support UniProt names
    # like "Collectin-10 (Collectin liver protein 1) (CL-L1)". Handle nested
    # parenthetical fragments iteratively (seen in some bulk UniProt aliases).
    stripped_source = normalize(
        re.sub(
            r"\s*\[[^\]]*\]",
            "",
            de_annotated or pre_cleaved or base,
        )
    )
    stripped = stripped_source
    while True:
        reduced = normalize(re.sub(r"\([^()]*\)", " ", stripped))
        if reduced == stripped:
            break
        stripped = reduced
    stripped = normalize(re.sub(r"[()]", " ", stripped))
    if stripped:
        variants.add(stripped)
    for inner in re.findall(r"\(([^)]+)\)", base):
        inner_clean = normalize(inner)
        # Ignore single-character parenthetical fragments (for example "G(I)/G(S)/G(T)")
        # because they create extremely noisy pseudo-symbols.
        if inner_clean and len(inner_clean) >= 2:
            # Ignore standalone Roman numerals (for example "VI"), which are
            # chain-family annotations and not analyte identity aliases.
            if re.fullmatch(r"[IVXLCDM]+", inner_clean.upper()):
                continue
            variants.add(inner_clean)
    return tuple(v for v in variants if v)


@lru_cache(maxsize=300_000)
def numeric_tokens(text: str) -> frozenset[str]:
    return frozenset(re.findall(r"\b\d+\b", text))


@lru_cache(maxsize=300_000)
def symbol_suffix_number(text: str) -> frozenset[str]:
    tokens = [tok for tok in re.split(r"[^a-z0-9]+", norm_key(text)) if tok]
    if not tokens:
        return frozenset()
    nums: set[str] = set()
    for tok in tokens:
        match = re.match(r"^[a-z\-]+(\d+)$", tok)
        if match:
            nums.add(match.group(1))
    if tokens[-1].isdigit() and any(any(ch.isalpha() for ch in tok) for tok in tokens[:-1]):
        nums.add(tokens[-1])
    return frozenset(nums)


@lru_cache(maxsize=300_000)
def normalize_measurement_subject(text: str) -> str:
    # Normalize labels/synonyms/aliases to the analyte phrase so identity checks
    # compare protein names rather than measurement wrappers.
    s = norm_key(text)
    if not s:
        return ""

    s = re.sub(r"^\s*(?:level|levels|amount|concentration)\s+of\s+", "", s)
    s = re.sub(
        r"\s+in\s+(?:blood serum|blood plasma|blood|plasma|serum|cerebrospinal fluid|csf|urine|saliva|tissue)\b.*$",
        "",
        s,
    )
    s = re.sub(
        r"^\s*(?:blood serum|blood plasma|blood|plasma|serum|cerebrospinal fluid|csf|urine|saliva|tissue)\s+",
        "",
        s,
    )
    s = re.sub(r"\s+(?:measurement|level|levels|amount|concentration)\b.*$", "", s)
    s = re.sub(r"\((?:human|mouse|rat)\)", "", s)
    s = re.sub(r"\s+", " ", s).strip(" -_/")
    return s


@lru_cache(maxsize=400_000)
def subject_phrase_variants(text: str) -> tuple[str, ...]:
    base = normalize_measurement_subject(text)
    if not base:
        return ()
    variants: set[str] = {base}
    variants.add(normalize(base.replace("-", " ")))
    variants.add(normalize(base.replace("/", " ")))

    # Apply controlled phrase equivalences to bridge UniProt/EFO wording gaps.
    for phrase in list(variants):
        if not phrase:
            continue
        mapped = phrase
        for pattern, replacement in SUBJECT_PHRASE_EQUIVALENTS:
            mapped = re.sub(pattern, replacement, mapped)
        mapped = normalize(mapped)
        if mapped:
            variants.add(mapped)

    # Canonical spacing form for robust token matching.
    canonicalized: set[str] = set()
    for phrase in variants:
        if not phrase:
            continue
        canonicalized.add(re.sub(r"\s+", " ", phrase).strip(" -_/"))
    return tuple(v for v in canonicalized if v)


def extract_subject_terms(label: str, synonyms: list[str]) -> frozenset[str]:
    terms: set[str] = set()
    for raw in [label] + synonyms:
        for subj in subject_phrase_variants(raw):
            if subj:
                terms.add(subj)
    return frozenset(terms)


@lru_cache(maxsize=800_000)
def subject_tokens(text: str) -> frozenset[str]:
    toks = [tok for tok in re.split(r"[^a-z0-9]+", norm_key(text)) if tok]
    keep = {
        tok
        for tok in toks
        if tok not in IDENTITY_STOP_TOKENS
    }
    return frozenset(keep)


@lru_cache(maxsize=800_000)
def subject_core_tokens(text: str) -> frozenset[str]:
    toks = subject_tokens(text)
    return frozenset(tok for tok in toks if len(tok) >= 2 or any(ch.isdigit() for ch in tok))


@lru_cache(maxsize=800_000)
def subject_suffix_token(text: str) -> str:
    toks = [tok for tok in re.split(r"[^a-z0-9]+", norm_key(text)) if tok]
    if not toks:
        return ""
    filtered = [tok for tok in toks if tok not in GENERIC_NAME_QUERY_TOKENS and tok not in GENERIC_TOKENS]
    if not filtered:
        return ""
    if (
        len(filtered) >= 2
        and len(filtered[-1]) == 1
        and filtered[-1].isalpha()
        and len(filtered[-2]) >= 2
    ):
        # Preserve member-letter identity in tails like "Ral-A"/"Ral-B".
        return filtered[-2] + filtered[-1]
    # Prefer a substantive trailing token over a single-letter tail token.
    for tok in reversed(filtered):
        if len(tok) >= 2 or any(ch.isdigit() for ch in tok):
            return tok
    return filtered[-1]


@lru_cache(maxsize=400_000)
def is_identifier_like_suffix(token: str) -> bool:
    if not token:
        return False
    if token in ROMAN_NUMERAL_TOKENS:
        return True
    if re.match(r"^[ivxlcdm]{1,4}[a-z]{1,2}$", token):
        # Family/member codes such as IA/IB/IC or IIA.
        return True
    if token.isdigit():
        return True
    if len(token) == 1 and token.isalpha():
        return True
    if re.match(r"^[a-z]+\d+$", token):
        return True
    if re.match(r"^\d+[a-z]+$", token):
        return True
    return False


@lru_cache(maxsize=400_000)
def split_symbol_tail(token: str) -> tuple[str, str] | None:
    t = norm_key(token)
    if not t:
        return None
    # Symbols with numeric tail (for example LTBP3, DEFA1, RHO6).
    m_num = re.match(r"^([a-z]{2,})(\d+)$", t)
    if m_num:
        return m_num.group(1), m_num.group(2)
    # Symbols with short letter tail (for example RHOF, RHOA).
    m_alpha = re.match(r"^([a-z]{2,})([a-z])$", t)
    if m_alpha and len(t) <= 6:
        return m_alpha.group(1), m_alpha.group(2)
    return None


def specific_identity_tokens(tokens: set[str]) -> set[str]:
    specific: set[str] = set()
    for tok in tokens:
        if tok in IDENTITY_AMBIGUOUS_TOKENS:
            continue
        if len(tok) >= 3 or any(ch.isdigit() for ch in tok):
            specific.add(tok)
    return specific


def subject_terms_match(
    profile_terms: frozenset[str],
    candidate_terms: frozenset[str],
    *,
    strict: bool = False,
) -> bool:
    if not profile_terms or not candidate_terms:
        return False
    # Fast prefilter before expensive pairwise term comparisons.
    profile_union: set[str] = set()
    for profile_term in profile_terms:
        profile_union.update(subject_core_tokens(profile_term))
    candidate_union: set[str] = set()
    for candidate_term in candidate_terms:
        candidate_union.update(subject_core_tokens(candidate_term))
    if not profile_union or not candidate_union:
        return False
    if profile_union.isdisjoint(candidate_union):
        return False
    profile_specific = specific_identity_tokens(profile_union)
    candidate_specific = specific_identity_tokens(candidate_union)
    if profile_specific and candidate_specific and profile_specific.isdisjoint(candidate_specific):
        return False
    for profile_term in profile_terms:
        profile_norm = norm_key(profile_term)
        profile_toks = subject_core_tokens(profile_term)
        if not profile_toks:
            continue
        for candidate_term in candidate_terms:
            candidate_norm = norm_key(candidate_term)
            if profile_norm == candidate_norm:
                return True
            profile_suffix = subject_suffix_token(profile_term)
            candidate_suffix = subject_suffix_token(candidate_term)
            profile_tail = split_symbol_tail(profile_suffix)
            candidate_tail = split_symbol_tail(candidate_suffix)
            if (
                profile_tail is not None
                and candidate_tail is not None
                and profile_tail[0] == candidate_tail[0]
                and profile_tail[1] != candidate_tail[1]
            ):
                # Reject family-member/symbol tail mismatches (for example RhoF vs Rho6).
                continue
            if (
                is_identifier_like_suffix(profile_suffix)
                and is_identifier_like_suffix(candidate_suffix)
                and profile_suffix != candidate_suffix
            ):
                continue
            candidate_toks = subject_core_tokens(candidate_term)
            if not candidate_toks:
                continue
            profile_type_toks = profile_toks & TYPE_SPECIFIER_TOKENS
            candidate_type_toks = candidate_toks & TYPE_SPECIFIER_TOKENS
            if profile_type_toks and candidate_type_toks and profile_type_toks.isdisjoint(candidate_type_toks):
                continue
            profile_num_toks = {tok for tok in profile_toks if any(ch.isdigit() for ch in tok)}
            candidate_num_toks = {tok for tok in candidate_toks if any(ch.isdigit() for ch in tok)}
            if profile_num_toks and candidate_num_toks and profile_num_toks.isdisjoint(candidate_num_toks):
                continue
            # A shared bare "1" is often non-specific in enzyme naming patterns
            # (for example "1,2-" vs "1,3-1,6-"), so compare nontrivial numeric
            # signatures as a second guard.
            profile_nontrivial_num_toks = {tok for tok in profile_num_toks if tok != "1"}
            candidate_nontrivial_num_toks = {tok for tok in candidate_num_toks if tok != "1"}
            if (
                profile_nontrivial_num_toks
                and candidate_nontrivial_num_toks
                and profile_nontrivial_num_toks.isdisjoint(candidate_nontrivial_num_toks)
            ):
                continue
            profile_roman_toks = {tok for tok in profile_toks if tok in ROMAN_NUMERAL_TOKENS}
            candidate_roman_toks = {tok for tok in candidate_toks if tok in ROMAN_NUMERAL_TOKENS}
            if profile_roman_toks and candidate_roman_toks and profile_roman_toks.isdisjoint(candidate_roman_toks):
                continue
            profile_id_toks = profile_num_toks | profile_roman_toks
            if profile_id_toks and not (profile_id_toks & candidate_toks):
                continue
            overlap = profile_toks & candidate_toks
            if not overlap:
                continue
            non_id_overlap = {tok for tok in overlap if not is_identifier_like_suffix(tok)}
            if not specific_identity_tokens(non_id_overlap):
                continue
            if strict:
                profile_non_id = {tok for tok in profile_toks if not is_identifier_like_suffix(tok)}
                candidate_non_id = {tok for tok in candidate_toks if not is_identifier_like_suffix(tok)}
                if not profile_non_id or not candidate_non_id:
                    continue
                shared_non_id = profile_non_id & candidate_non_id
                profile_specific = specific_identity_tokens(profile_non_id)
                candidate_specific = specific_identity_tokens(candidate_non_id)
                if profile_specific and candidate_specific and not (profile_specific & candidate_specific):
                    continue
                if profile_non_id == candidate_non_id:
                    return True
                if profile_non_id.issubset(candidate_non_id):
                    extra = candidate_non_id - profile_non_id
                    if extra & IDENTITY_RELATION_TOKENS:
                        continue
                    if specific_identity_tokens(extra):
                        # Do not accept when one side drops/adds a specific identity token.
                        continue
                    if len(extra) <= 1 and len(shared_non_id) >= 2:
                        return True
                if candidate_non_id.issubset(profile_non_id):
                    extra = profile_non_id - candidate_non_id
                    if extra & IDENTITY_RELATION_TOKENS:
                        continue
                    if specific_identity_tokens(extra):
                        continue
                    if len(extra) <= 1 and len(shared_non_id) >= 2:
                        return True
                continue
            min_len = min(len(profile_toks), len(candidate_toks))
            if profile_toks.issubset(candidate_toks):
                extra = candidate_toks - profile_toks
                if any(is_identifier_like_suffix(tok) for tok in extra):
                    continue
                if extra & IDENTITY_RELATION_TOKENS:
                    continue
                if len(extra) <= 2:
                    return True
            if candidate_toks.issubset(profile_toks):
                extra = profile_toks - candidate_toks
                if any(is_identifier_like_suffix(tok) for tok in extra):
                    continue
                if extra & IDENTITY_RELATION_TOKENS:
                    continue
                if len(extra) <= 2:
                    return True
            if len(overlap) >= 2 and (len(overlap) / max(1, min_len)) >= 0.6:
                if len(non_id_overlap) >= 2:
                    return True
    return False


def is_symbol_like_alias(alias: str) -> bool:
    raw = normalize(alias)
    if not raw:
        return False
    if len(raw) < 2:
        return False
    if len(raw) > 16 or " " in raw:
        return False
    if raw != raw.upper():
        return False
    if re.fullmatch(r"[IVXLCDM]+", raw):
        return False
    return re.match(r"^[A-Z0-9][A-Z0-9\-]*$", raw) is not None


def is_weak_symbol_alias(alias: str) -> bool:
    raw = normalize(alias).upper()
    if not is_symbol_like_alias(raw):
        return False
    # Very short symbols (for example C2) are highly ambiguous with metabolite shorthands.
    return len(raw) <= 3


def subject_suffix_numbers(subject_terms: frozenset[str]) -> frozenset[str]:
    nums: set[str] = set()
    for term in subject_terms:
        nums |= symbol_suffix_number(term)
    return frozenset(nums)


@lru_cache(maxsize=300_000)
def informative_subject_tokens(subject_terms: frozenset[str]) -> frozenset[str]:
    informative: set[str] = set()
    for term in subject_terms:
        for tok in subject_core_tokens(term):
            if tok in IDENTITY_STOP_TOKENS or tok in IDENTITY_AMBIGUOUS_TOKENS:
                continue
            if is_identifier_like_suffix(tok):
                continue
            if len(tok) >= 3:
                informative.add(tok)
    return frozenset(informative)


def symbol_query_key(query: str) -> str:
    raw = normalize(query).upper()
    if is_symbol_like_alias(raw):
        return raw
    return ""


def coag_factor_symbol_key(query: str, index: dict[str, Any]) -> str:
    raw = normalize(query).upper()
    match = re.match(r"^FA([0-9]{1,2})$", raw)
    if not match:
        return ""
    candidate = f"F{match.group(1)}"
    symbol_index = index.get("symbol_accession_index", {})
    hits = symbol_index.get(candidate, [])
    if not hits:
        return ""

    accession_index = index.get("accession_alias_index", {})
    # Require coagulation-factor evidence on at least one accession to avoid
    # treating arbitrary FA<digits> tokens as gene symbols.
    for hit in hits:
        aliases = [normalize(a) for a in accession_index.get(normalize(hit).upper(), []) if normalize(a)]
        joined = " ".join(aliases).lower()
        if "coagulation factor" in joined or "hageman factor" in joined:
            return candidate
    return ""


def defensin_alpha_shorthand_symbol_key(query: str, index: dict[str, Any]) -> str:
    raw = normalize(query).upper()
    match = re.match(r"^DEF([0-9]{1,3})$", raw)
    if not match:
        return ""
    candidate = f"DEFA{match.group(1)}"
    symbol_index = index.get("symbol_accession_index", {})
    hits = symbol_index.get(candidate, [])
    if not hits:
        return ""

    accession_index = index.get("accession_alias_index", {})
    for hit in hits:
        aliases = [normalize(a) for a in accession_index.get(normalize(hit).upper(), []) if normalize(a)]
        joined = " ".join(aliases).lower()
        if "defensin, alpha" in joined or "neutrophil defensin" in joined:
            return candidate
    return ""


def uniprot_mnemonic_symbol_keys(query: str, index: dict[str, Any]) -> list[str]:
    raw = normalize(query).upper()
    match = re.match(r"^([A-Z0-9]+)_HUMAN$", raw)
    if not match:
        return []

    stem = match.group(1)
    symbol_index = index.get("symbol_accession_index", {})
    keys: list[str] = []
    if stem in symbol_index and symbol_index.get(stem):
        keys.append(stem)

    factor_key = coag_factor_symbol_key(stem, index)
    if factor_key and factor_key not in keys:
        keys.append(factor_key)
    return keys


def loose_key(text: str) -> str:
    toks = re.split(r"[^a-z0-9]+", norm_key(text))
    return " ".join(tok for tok in toks if tok)


def informative_lookup_tokens(text: str) -> set[str]:
    toks = tokenize(norm_key(text))
    keep: set[str] = set()
    for tok in toks:
        if tok in GENERIC_NAME_QUERY_TOKENS or tok in GENERIC_TOKENS:
            continue
        if tok.isdigit():
            continue
        keep.add(tok)
    return keep


def is_informative_lookup_variant(text: str) -> bool:
    if not text:
        return False
    if canonical_accession(text) or is_symbol_like_alias(text):
        return True
    toks = informative_lookup_tokens(text)
    if not toks:
        return False
    return any(any(ch.isalpha() for ch in tok) and len(tok) >= 3 for tok in toks)


def rank_accession_hits(hits: set[str], index: dict[str, Any]) -> list[str]:
    alias_index = index.get("accession_alias_index", {})
    return sorted(
        hits,
        key=lambda acc: (
            acc.startswith("A0A"),
            -len(alias_index.get(acc, [])),
            acc,
        ),
    )


@lru_cache(maxsize=200_000)
def is_reviewed_like_accession(accession: str) -> bool:
    acc = normalize(accession).upper()
    return bool(UNIPROT_REVIEWED_LIKE_RE.match(acc))


def accession_descriptive_alias_length(accession: str, index: dict[str, Any]) -> int:
    acc = normalize(accession).upper()
    if not acc:
        return 0
    accession_index = index.get("accession_alias_index", {})
    aliases = accession_index.get(acc, []) if isinstance(accession_index, dict) else []
    best = 0
    for alias in aliases:
        alias_clean = normalize(alias)
        if not alias_clean:
            continue
        if canonical_accession(alias_clean):
            continue
        if is_symbol_like_alias(alias_clean):
            continue
        alias_norm = normalize(alias_clean)
        if not alias_norm:
            continue
        best = max(best, len(alias_norm))
    return best


def prefer_reviewed_accessions(accessions: list[str], index: dict[str, Any]) -> list[str]:
    # Accessions come from tied top alias/symbol scores. For symbol/name input,
    # prefer reviewed-like primaries and then prefer the accession carrying the
    # richest descriptive alias text.
    normalized = sorted({normalize(acc).upper() for acc in accessions if normalize(acc)})
    if len(normalized) <= 1:
        return normalized

    reviewed = [acc for acc in normalized if is_reviewed_like_accession(acc)]
    pool = reviewed if reviewed else normalized
    if len(pool) <= 1:
        return pool

    accession_index = index.get("accession_alias_index", {})
    scored: list[tuple[str, tuple[int, int, int, str]]] = []
    for acc in pool:
        aliases = accession_index.get(acc, []) if isinstance(accession_index, dict) else []
        descriptive_len = accession_descriptive_alias_length(acc, index)
        alias_count = len(aliases)
        # Prefer non-A0A accessions as a weak tie-breaker after descriptive evidence.
        non_a0a = 1 if not acc.startswith("A0A") else 0
        scored.append((acc, (descriptive_len, alias_count, non_a0a, acc)))

    best_score = max(score for _acc, score in scored)
    winners = sorted([acc for acc, score in scored if score == best_score])
    if len(winners) == 1:
        return winners
    return winners


def name_lookup_variants(query: str) -> list[str]:
    base = normalize(query)
    variants: set[str] = set(query_variants(query))
    variants.add(base)

    stripped = base
    while True:
        reduced = normalize(re.sub(r"\([^()]*\)", " ", stripped))
        if reduced == stripped:
            break
        stripped = reduced
    if stripped:
        variants.add(stripped)
    before_paren = normalize(base.split("(", 1)[0])
    if before_paren:
        variants.add(before_paren)

    for inner in re.findall(r"\(([^)]+)\)", base):
        inner_clean = normalize(inner)
        if inner_clean:
            variants.add(inner_clean)

    expanded: set[str] = set()
    for v in list(variants):
        expanded.add(v)
        punct_norm = normalize(re.sub(r"[^A-Za-z0-9\- ]+", " ", v))
        if punct_norm:
            expanded.add(punct_norm)
            sf = re.sub(r"\bsub[- ]*family member (\d+)\b", r"-\1", punct_norm, flags=re.IGNORECASE)
            sf = normalize(sf)
            if sf:
                expanded.add(sf)

        for part in re.split(r"[:/]|,\s+(?=[A-Za-z])", v):
            part_clean = normalize(part)
            if part_clean:
                expanded.add(part_clean)

    return [v for v in expanded if v]


def strict_query_tokens(query: str) -> frozenset[str]:
    subject = normalize_measurement_subject(query)
    base = subject if subject else norm_key(query)
    toks = {
        tok
        for tok in tokenize(base)
        if len(tok) >= 3 and tok not in GENERIC_NAME_QUERY_TOKENS and tok not in GENERIC_TOKENS
    }
    return frozenset(toks)


def resolve_query_accessions(query: str, index: dict[str, Any], input_type: str = "auto") -> tuple[list[str], bool, bool]:
    # Returns: (resolved accessions, used_alias_resolution, ambiguous_alias_resolution)
    # Composite ratio-like analyte strings (for example MUFA/PUFA, ApoB/ApoA1 ratio)
    # should not collapse to a single protein accession via alias expansion.
    if query_requests_ratio(query):
        return [], False, False

    input_kind = normalize_input_type(input_type)
    accession_index = index.get("alias_accession_index", {})
    symbol_index = index.get("symbol_accession_index", {})
    gene_id_index = index.get("gene_id_accession_index", {})
    loose_index = index.get("alias_loose_index", {})
    accession_scores: dict[str, int] = {}
    used_alias_resolution = False
    ambiguous_alias_resolution = False
    variants = name_lookup_variants(query)

    for variant in variants:
        acc = canonical_accession(variant)
        if acc:
            return [acc], False, False

    gene_id_query_variants: set[str] = set()
    for variant in variants:
        for gid in gene_id_alias_variants(
            variant,
            include_unprefixed_numeric=(input_kind == "gene_id"),
        ):
            gene_id_query_variants.add(gid)

    if gene_id_query_variants:
        for gid in sorted(gene_id_query_variants):
            hits = gene_id_index.get(norm_key(gid), [])
            if not hits:
                continue
            used_alias_resolution = True
            for hit in hits:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + 10
    elif input_kind == "gene_id":
        # For explicit gene-id input, avoid falling back to symbol/name heuristics.
        return [], False, False

    factor_symbol = coag_factor_symbol_key(query, index)
    if factor_symbol:
        hits = symbol_index.get(factor_symbol, [])
        if hits:
            used_alias_resolution = True
            for hit in hits:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + 6

    for mnemonic_symbol in uniprot_mnemonic_symbol_keys(query, index):
        hits = symbol_index.get(mnemonic_symbol, [])
        if hits:
            used_alias_resolution = True
            for hit in hits:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + 6

    defensin_symbol = defensin_alpha_shorthand_symbol_key(query, index)
    if defensin_symbol:
        hits = symbol_index.get(defensin_symbol, [])
        if hits:
            used_alias_resolution = True
            for hit in hits:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + 5

    symbol_key = symbol_query_key(query)
    if symbol_key:
        hits = symbol_index.get(symbol_key, [])
        if hits:
            used_alias_resolution = True
            for hit in hits:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + 5
    for variant in variants:
        variant_symbol = symbol_query_key(variant)
        if not variant_symbol:
            continue
        hits = symbol_index.get(variant_symbol, [])
        if hits:
            used_alias_resolution = True
            for hit in hits:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + 5

    for variant in variants:
        if not is_informative_lookup_variant(variant):
            continue
        variant_specificity = len(informative_lookup_tokens(variant))
        key = norm_key(variant)
        if not key:
            continue
        exact_hits = set(accession_index.get(key, []))
        if len(exact_hits) == 1:
            used_alias_resolution = True
            only = normalize(next(iter(exact_hits))).upper()
            if only:
                accession_scores[only] = accession_scores.get(only, 0) + 3 + min(2, variant_specificity)
            continue
        if len(exact_hits) > 1:
            if variant_specificity < 2 and len(exact_hits) > 8:
                continue
            ambiguous_alias_resolution = True
            used_alias_resolution = True
            hit_weight = max(1, min(4, variant_specificity))
            for hit in rank_accession_hits(exact_hits, index)[:25]:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + hit_weight
            continue

        loose = loose_key(variant)
        if not loose:
            continue
        loose_hits = set(loose_index.get(loose, []))
        if len(loose_hits) == 1:
            used_alias_resolution = True
            only = normalize(next(iter(loose_hits))).upper()
            if only:
                accession_scores[only] = accession_scores.get(only, 0) + 1 + min(2, variant_specificity)
        elif len(loose_hits) > 1:
            if variant_specificity < 2 and len(loose_hits) > 8:
                continue
            ambiguous_alias_resolution = True
            used_alias_resolution = True
            hit_weight = max(1, min(3, variant_specificity))
            for hit in rank_accession_hits(loose_hits, index)[:15]:
                hit_norm = normalize(hit).upper()
                if hit_norm:
                    accession_scores[hit_norm] = accession_scores.get(hit_norm, 0) + hit_weight

    if accession_scores:
        top = max(accession_scores.values())
        chosen = sorted(acc for acc, score in accession_scores.items() if score == top)
        # For non-accession inputs, collapse tied symbol/alias hits to a
        # reviewed-prioritized subset when possible.
        if not any(canonical_accession(v) for v in variants):
            chosen = prefer_reviewed_accessions(chosen, index)
        return chosen, used_alias_resolution, ambiguous_alias_resolution
    return [], used_alias_resolution, ambiguous_alias_resolution


def accession_profiles_for_accessions(accessions: list[str], index: dict[str, Any]) -> list[list[str]]:
    profiles: list[list[str]] = []
    accession_index = index.get("accession_alias_index", {})
    for acc in sorted(set(accessions)):
        aliases = [acc] + [normalize(a) for a in accession_index.get(acc, []) if normalize(a)]
        # Deduplicate while preserving stable order.
        seen: set[str] = set()
        ordered: list[str] = []
        for alias in aliases:
            key = norm_key(alias)
            if key in seen:
                continue
            seen.add(key)
            ordered.append(alias)
        profiles.append(ordered)
    return profiles


def compile_accession_profiles(profiles: list[list[str]]) -> tuple[AccessionProfile, ...]:
    compiled: list[AccessionProfile] = []
    for aliases in profiles:
        alias_terms: list[str] = []
        seen_terms: set[str] = set()
        expected_symbol_nums: set[str] = set()
        subject_terms: set[str] = set()
        symbol_terms: set[str] = set()
        for alias in aliases:
            alias_clean = normalize(alias)
            if not alias_clean:
                continue
            for alias_variant in alias_variants(alias_clean):
                alias_key = norm_key(alias_variant)
                if alias_key in seen_terms:
                    continue
                seen_terms.add(alias_key)
                alias_terms.append(alias_variant)
                for subject in subject_phrase_variants(alias_variant):
                    if subject:
                        subject_terms.add(subject)
                # Do not treat accession identifiers as gene/symbol aliases;
                # otherwise numeric mismatch checks can incorrectly reject
                # valid candidates (for example P26447 vs S100A4).
                if is_symbol_like_alias(alias_variant) and not canonical_accession(alias_variant):
                    symbol_terms.add(alias_variant.lower())
                    expected_symbol_nums |= symbol_suffix_number(alias_variant)
        if alias_terms:
            compiled.append(
                AccessionProfile(
                    alias_terms=tuple(alias_terms),
                    expected_symbol_nums=frozenset(expected_symbol_nums),
                    subject_terms=frozenset(subject_terms),
                    symbol_terms=tuple(sorted(symbol_terms)),
                )
            )
    return tuple(compiled)


def protein_term_roundtrip_queries(
    *,
    label: str,
    synonyms: list[str],
    subject_terms: frozenset[str] | None = None,
    meta_subject: str = "",
) -> list[str]:
    queries: list[str] = []
    seen: set[str] = set()

    def add_query(value: str) -> None:
        cleaned = normalize(value)
        if not cleaned:
            return
        key = norm_key(cleaned)
        if not key or key in seen:
            return
        seen.add(key)
        queries.append(cleaned)

    primary_subject = normalize(meta_subject or normalize_measurement_subject(label))
    add_query(primary_subject)

    terms = subject_terms if subject_terms is not None else extract_subject_terms(label, synonyms)
    ranked_terms = sorted(
        (normalize(term) for term in terms if normalize(term)),
        key=lambda term: (
            len(informative_lookup_tokens(term)),
            len(term),
            term,
        ),
        reverse=True,
    )
    for term in ranked_terms:
        add_query(term)
    for synonym in synonyms[:8]:
        add_query(normalize_measurement_subject(synonym))
    add_query(normalize_measurement_subject(label))
    return queries


def resolve_term_protein_accessions(
    efo_id: str,
    index: dict[str, Any],
    *,
    label: str = "",
    synonyms: list[str] | None = None,
    subject_terms: frozenset[str] | None = None,
) -> tuple[str, ...]:
    term_meta = index.get("term_meta", {})
    meta = term_meta.get(efo_id, {}) if isinstance(term_meta, dict) else {}
    if isinstance(meta, dict) and "protein_accessions" in meta:
        cached = tuple(
            sorted(
                {
                    normalize(accession).upper()
                    for accession in (meta.get("protein_accessions") or [])
                    if canonical_accession(accession)
                }
            )
        )
        meta["protein_accessions"] = list(cached)
        return cached

    label_clean = normalize(label or (meta.get("label") if isinstance(meta, dict) else "") or "")
    syns = (
        [normalize(x) for x in synonyms if normalize(x)]
        if synonyms is not None
        else [
            normalize(x)
            for x in ((meta.get("synonyms") if isinstance(meta, dict) else []) or [])
            if normalize(x)
        ]
    )
    meta_subject = normalize((meta.get("subject") if isinstance(meta, dict) else "") or "")
    queries = protein_term_roundtrip_queries(
        label=label_clean,
        synonyms=syns,
        subject_terms=subject_terms,
        meta_subject=meta_subject,
    )

    best_hits: tuple[str, ...] = ()
    best_rank: tuple[int, int, int, int, str] | None = None
    primary_subject_key = norm_key(meta_subject or normalize_measurement_subject(label_clean))
    for subject_query in queries:
        hits, _used_alias_resolution, _ambiguous_alias_resolution = resolve_query_accessions(
            subject_query,
            index,
            input_type="protein_name",
        )
        normalized_hits = tuple(
            sorted(
                {
                    normalize(accession).upper()
                    for accession in hits
                    if canonical_accession(accession)
                }
            )
        )
        if not normalized_hits:
            continue
        rank = (
            1 if len(normalized_hits) == 1 else 0,
            1 if norm_key(subject_query) == primary_subject_key and primary_subject_key else 0,
            len(informative_lookup_tokens(subject_query)),
            -len(normalized_hits),
            subject_query,
        )
        if best_rank is None or rank > best_rank:
            best_rank = rank
            best_hits = normalized_hits

    if isinstance(meta, dict):
        meta["protein_accessions"] = list(best_hits)
    return best_hits


def resolve_term_protein_symbols(
    efo_id: str,
    index: dict[str, Any],
    *,
    label: str = "",
    synonyms: list[str] | None = None,
    subject_terms: frozenset[str] | None = None,
) -> tuple[str, ...]:
    term_meta = index.get("term_meta", {})
    meta = term_meta.get(efo_id, {}) if isinstance(term_meta, dict) else {}
    if isinstance(meta, dict) and "protein_symbols" in meta:
        cached = tuple(
            sorted(
                {
                    normalize(symbol).upper()
                    for symbol in (meta.get("protein_symbols") or [])
                    if is_symbol_like_alias(symbol)
                }
            )
        )
        meta["protein_symbols"] = list(cached)
        return cached

    candidate_accessions = resolve_term_protein_accessions(
        efo_id,
        index,
        label=label,
        synonyms=synonyms,
        subject_terms=subject_terms,
    )
    symbols: set[str] = set()
    for accession in candidate_accessions:
        symbols |= accession_symbol_aliases(accession, index)

    label_clean = normalize(label or (meta.get("label") if isinstance(meta, dict) else "") or "")
    syns = (
        [normalize(x) for x in synonyms if normalize(x)]
        if synonyms is not None
        else [
            normalize(x)
            for x in ((meta.get("synonyms") if isinstance(meta, dict) else []) or [])
            if normalize(x)
        ]
    )
    meta_subject = normalize((meta.get("subject") if isinstance(meta, dict) else "") or "")
    for subject_query in protein_term_roundtrip_queries(
        label=label_clean,
        synonyms=syns,
        subject_terms=subject_terms,
        meta_subject=meta_subject,
    ):
        symbol = symbol_query_key(subject_query)
        if symbol:
            symbols.add(symbol)
        for mnemonic_symbol in uniprot_mnemonic_symbol_keys(subject_query, index):
            symbols.add(mnemonic_symbol)
        factor_symbol = coag_factor_symbol_key(subject_query, index)
        if factor_symbol:
            symbols.add(factor_symbol)
        defensin_symbol = defensin_alpha_shorthand_symbol_key(subject_query, index)
        if defensin_symbol:
            symbols.add(defensin_symbol)

    normalized = tuple(sorted(symbol for symbol in symbols if is_symbol_like_alias(symbol)))
    if isinstance(meta, dict):
        meta["protein_symbols"] = list(normalized)
    return normalized


def query_protein_roundtrip_symbols(
    query: str,
    resolved_accessions: list[str] | tuple[str, ...],
    index: dict[str, Any],
) -> tuple[str, ...]:
    symbols: set[str] = set()
    for accession in resolved_accessions:
        symbols |= accession_symbol_aliases(accession, index)
    query_symbol = symbol_query_key(query)
    if query_symbol:
        symbols.add(query_symbol)
    return tuple(sorted(symbol for symbol in symbols if is_symbol_like_alias(symbol)))


def protein_roundtrip_status(
    query_accessions: list[str] | tuple[str, ...],
    query_symbols: list[str] | tuple[str, ...],
    candidate_accessions: list[str] | tuple[str, ...],
    candidate_symbols: list[str] | tuple[str, ...],
) -> str:
    query_set = {
        normalize(accession).upper()
        for accession in query_accessions
        if canonical_accession(accession)
    }
    query_symbol_set = {
        normalize(symbol).upper()
        for symbol in query_symbols
        if is_symbol_like_alias(symbol)
    }
    if not query_set and not query_symbol_set:
        return "not_applicable"
    candidate_set = {
        normalize(accession).upper()
        for accession in candidate_accessions
        if canonical_accession(accession)
    }
    if candidate_set & query_set:
        return "accession_match"
    candidate_symbol_set = {
        normalize(symbol).upper()
        for symbol in candidate_symbols
        if is_symbol_like_alias(symbol)
    }
    if candidate_symbol_set & query_symbol_set:
        return "symbol_match"
    if not candidate_set and not candidate_symbol_set:
        return "unresolved"
    return "mismatch"


def protein_accession_roundtrip_evidence(
    roundtrip_status: str,
    candidate_accessions: list[str] | tuple[str, ...],
    candidate_symbols: list[str] | tuple[str, ...],
) -> str:
    if roundtrip_status == "not_applicable":
        return ""
    accessions = [
        normalize(accession).upper()
        for accession in candidate_accessions
        if canonical_accession(accession)
    ]
    symbols = [
        normalize(symbol).upper()
        for symbol in candidate_symbols
        if is_symbol_like_alias(symbol)
    ]
    accession_text = "|".join(accessions) if accessions else "none"
    symbol_text = "|".join(symbols) if symbols else "none"
    return (
        f"protein_roundtrip={roundtrip_status}; "
        f"term_accessions={accession_text}; term_symbols={symbol_text}"
    )


def accession_identity_match(
    query: str,
    label: str,
    synonyms: list[str],
    index: dict[str, Any],
    precomputed_profiles: tuple[AccessionProfile, ...] | None = None,
    candidate_text: str | None = None,
    candidate_nums: frozenset[str] | None = None,
    candidate_subject_terms: frozenset[str] | None = None,
) -> bool:
    profiles = (
        precomputed_profiles
        if precomputed_profiles is not None
        else compile_accession_profiles(
            accession_profiles_for_accessions(resolve_query_accessions(query, index)[0], index)
        )
    )
    # If query is not accession-based, keep legacy behavior.
    if not profiles:
        return True

    merged = candidate_text if candidate_text is not None else normalize(" ".join([label] + synonyms))
    merged_lower = merged.lower()
    subject_terms = candidate_subject_terms if candidate_subject_terms is not None else extract_subject_terms(label, synonyms)
    candidate_symbol_nums = subject_suffix_numbers(subject_terms)
    candidate_informative_tokens = informative_subject_tokens(subject_terms)
    candidate_type_tokens: set[str] = set()
    for term in subject_terms:
        candidate_type_tokens |= set(subject_core_tokens(term)) & TYPE_SPECIFIER_TOKENS
    query_is_accession_input = any(canonical_accession(v) for v in query_variants(query))
    symbol_query = (not query_is_accession_input) and bool(
        symbol_query_key(query) or uniprot_mnemonic_symbol_keys(query, index)
    )

    for profile in profiles:
        profile_informative_tokens = informative_subject_tokens(profile.subject_terms)
        profile_type_tokens: set[str] = set()
        for term in profile.subject_terms:
            profile_type_tokens |= set(subject_core_tokens(term)) & TYPE_SPECIFIER_TOKENS
        if profile_type_tokens and candidate_type_tokens and profile_type_tokens.isdisjoint(candidate_type_tokens):
            continue
        strong_symbols = [sym for sym in profile.symbol_terms if not is_weak_symbol_alias(sym)]
        weak_symbols = [sym for sym in profile.symbol_terms if is_weak_symbol_alias(sym)]
        symbol_match_strong = any(contains_term(merged_lower, sym) for sym in strong_symbols)
        symbol_match_weak = any(contains_term(merged_lower, sym) for sym in weak_symbols)
        symbol_match = symbol_match_strong or symbol_match_weak
        subject_match = subject_terms_match(profile.subject_terms, subject_terms)
        if not symbol_match and not subject_match:
            continue
        if symbol_match_weak and not symbol_match_strong and profile_informative_tokens:
            # Weak short-symbol hits (for example "C2") are only accepted when
            # there is informative subject-token agreement.
            if not (profile_informative_tokens & candidate_informative_tokens):
                continue
            if not symbol_query and not subject_terms_match(profile.subject_terms, subject_terms, strict=True):
                continue
        # When we only have a loose subject phrase match (no symbol hit), require
        # strict subject compatibility to avoid family-member drift.
        if not symbol_match and subject_match:
            if not subject_terms_match(profile.subject_terms, subject_terms, strict=True):
                continue
        if (
            profile.expected_symbol_nums
            and candidate_symbol_nums
            and profile.expected_symbol_nums.isdisjoint(candidate_symbol_nums)
        ):
            # Reject number-family mismatch (e.g., protein 4 vs protein 1).
            continue
        return True
    return False


def profile_subject_exact_match_bonus(
    profiles: tuple[AccessionProfile, ...],
    candidate_subject_terms: frozenset[str],
) -> int:
    if not profiles or not candidate_subject_terms:
        return 0
    candidate_keys = {norm_key(t) for t in candidate_subject_terms if norm_key(t)}
    if not candidate_keys:
        return 0
    for profile in profiles:
        profile_keys = {norm_key(t) for t in profile.subject_terms if norm_key(t)}
        if candidate_keys & profile_keys:
            return 1
    return 0


def profile_subject_exact_match_rank(
    profiles: tuple[AccessionProfile, ...],
    candidate_subject_terms: frozenset[str],
) -> int:
    if not profiles or not candidate_subject_terms:
        return 0
    candidate_keys = {norm_key(t) for t in candidate_subject_terms if norm_key(t)}
    if not candidate_keys:
        return 0
    for profile in profiles:
        profile_keys = {norm_key(t) for t in profile.subject_terms if norm_key(t)}
        if candidate_keys & profile_keys:
            return 1
    return 0


def is_metabolite_like_input_type(input_type: str) -> bool:
    kind = normalize_input_type(input_type)
    return kind in {"metabolite_id", "metabolite_name"}


def is_protein_like_input_type(input_type: str) -> bool:
    kind = normalize_input_type(input_type)
    return kind in {"accession", "gene_symbol", "gene_id", "protein_name"}


def resolve_query_metabolite_concepts(query: str, index: dict[str, Any]) -> tuple[list[str], bool, bool]:
    id_index = index.get("metabolite_id_concept_index", {})
    alias_index = index.get("metabolite_alias_concept_index", {})
    loose_index = index.get("metabolite_alias_loose_index", {})
    concept_scores: dict[str, int] = {}
    used_alias_resolution = False
    ambiguous_alias_resolution = False

    variants = name_lookup_variants(query)
    for variant in variants:
        met_id = canonical_metabolite_id(variant)
        if met_id:
            hits = id_index.get(met_id, [])
            if len(hits) == 1:
                return [normalize(hits[0])], False, False
            if hits:
                for hit in hits:
                    concept = normalize(hit)
                    if concept:
                        concept_scores[concept] = concept_scores.get(concept, 0) + 8

        key = norm_key(variant)
        if key:
            if (
                key in METABOLITE_CONCEPT_ALIAS_EXCLUDE
                or key.startswith("cholesteryl ester")
                or key.startswith("cholesteryl esters")
            ):
                continue
            hits = alias_index.get(key, [])
            if len(hits) == 1:
                used_alias_resolution = True
                concept = normalize(hits[0])
                if concept:
                    concept_scores[concept] = concept_scores.get(concept, 0) + 6
            elif len(hits) > 1:
                used_alias_resolution = True
                ambiguous_alias_resolution = True
                for hit in sorted({normalize(h) for h in hits if normalize(h)})[:25]:
                    concept_scores[hit] = concept_scores.get(hit, 0) + 4

        loose = loose_key(variant)
        if loose:
            hits = loose_index.get(loose, [])
            if len(hits) == 1:
                used_alias_resolution = True
                concept = normalize(hits[0])
                if concept:
                    concept_scores[concept] = concept_scores.get(concept, 0) + 3
            elif len(hits) > 1:
                used_alias_resolution = True
                ambiguous_alias_resolution = True
                for hit in sorted({normalize(h) for h in hits if normalize(h)})[:20]:
                    concept_scores[hit] = concept_scores.get(hit, 0) + 2

    if concept_scores:
        ranked = sorted(concept_scores.items(), key=lambda item: (-item[1], item[0]))
        top_score = ranked[0][1]
        chosen = sorted(c for c, score in ranked if score == top_score)
        # Avoid forcing unstable concept identity from noisy alias collisions.
        if len(chosen) > 1:
            return [], used_alias_resolution, True
        if (
            ambiguous_alias_resolution
            and not any(canonical_metabolite_id(v) for v in query_variants(query))
            and top_score < 7
        ):
            return [], used_alias_resolution, True
        return chosen, used_alias_resolution, ambiguous_alias_resolution
    return [], used_alias_resolution, ambiguous_alias_resolution


def detect_entity_route(
    *,
    query: str,
    input_type: str,
    forced_entity_type: str,
    resolved_accessions: list[str],
    resolved_metabolites: list[str],
) -> str:
    forced = normalize_entity_type(forced_entity_type)
    if forced in {"protein", "metabolite"}:
        return forced
    # Keep mixed NMR-style files usable even when a row is marked metabolite
    # but clearly names a protein analyte (for example Apo-A1/Apo-B100).
    if query_has_protein_cue(query) and resolved_accessions:
        return "protein"
    if query_has_protein_cue(query):
        return "protein"
    if is_protein_like_input_type(input_type):
        return "protein"
    if is_metabolite_like_input_type(input_type):
        return "metabolite"
    if any(canonical_accession(v) for v in query_variants(query)):
        return "protein"
    if any(canonical_metabolite_id(v) for v in query_variants(query)):
        return "metabolite"
    if resolved_accessions:
        return "protein"
    if resolved_metabolites:
        return "metabolite"
    if query_has_metabolite_cue(query):
        return "metabolite"
    return "protein"


def has_useful_external_identity(
    *,
    query: str,
    route: str,
    resolved_accessions: list[str],
    resolved_metabolites: list[str],
    index: dict[str, Any],
) -> bool:
    if route == "protein":
        if any(canonical_accession(v) for v in query_variants(query)):
            return True
        return bool(resolved_accessions)

    if route == "metabolite":
        if any(canonical_metabolite_id(v) for v in query_variants(query)):
            return True
        concept_index = index.get("metabolite_concept_index", {})
        for concept_id in resolved_metabolites:
            concept_meta = concept_index.get(concept_id) or {}
            if not isinstance(concept_meta, dict):
                continue
            for met_id in concept_meta.get("ids") or []:
                if canonical_metabolite_id(met_id):
                    return True
        return False

    return False


def concept_subject_terms(concept_ids: list[str], index: dict[str, Any]) -> tuple[frozenset[str], ...]:
    concept_index = index.get("metabolite_concept_index", {})
    profiles: list[frozenset[str]] = []
    for concept_id in concept_ids:
        meta = concept_index.get(concept_id) or {}
        if not isinstance(meta, dict):
            continue
        subject_terms = frozenset(normalize(x) for x in (meta.get("subject_terms") or []) if normalize(x))
        if subject_terms:
            profiles.append(subject_terms)
    return tuple(profiles)


def metabolite_identity_match(
    query: str,
    candidate_subject_terms: frozenset[str],
    candidate_nums: frozenset[str],
    concept_subject_profiles: tuple[frozenset[str], ...],
) -> bool:
    if not concept_subject_profiles:
        return True
    query_nums = numeric_tokens(normalize_measurement_subject(query) or norm_key(query))
    candidate_core: set[str] = set()
    for candidate_term in candidate_subject_terms:
        candidate_core.update(subject_core_tokens(candidate_term))
    candidate_core = specific_identity_tokens(candidate_core)
    candidate_num_tokens = {tok for tok in candidate_core if any(ch.isdigit() for ch in tok)}
    for profile_terms in concept_subject_profiles:
        if not profile_terms:
            continue
        profile_core: set[str] = set()
        for profile_term in profile_terms:
            profile_core.update(subject_core_tokens(profile_term))
        profile_core = specific_identity_tokens(profile_core)
        if profile_core and candidate_core and profile_core.isdisjoint(candidate_core):
            continue
        profile_num_tokens = {tok for tok in profile_core if any(ch.isdigit() for ch in tok)}
        if profile_num_tokens and candidate_num_tokens and profile_num_tokens.isdisjoint(candidate_num_tokens):
            continue
        # Metabolite concept validation must stay strict to avoid near-name drift
        # (for example creatinine aliases like "creatine anhydride" incorrectly
        # validating plain creatine measurements).
        if subject_terms_match(profile_terms, candidate_subject_terms, strict=True):
            if query_nums and candidate_nums and query_nums.isdisjoint(candidate_nums):
                continue
            return True
    return False


def build_metabolite_query_aliases(
    query: str,
    index: dict[str, Any],
    resolved_concepts: list[str],
    name_mode: str,
    is_metabolite_id_input: bool,
) -> list[str]:
    aliases: set[str] = set()
    concept_index = index.get("metabolite_concept_index", {})

    if name_mode == "strict" and not is_metabolite_id_input and not resolved_concepts:
        for variant in name_lookup_variants(query):
            aliases.add(variant)
        return [x for x in aliases if x]

    for variant in name_lookup_variants(query):
        aliases.add(variant)
        met_id = canonical_metabolite_id(variant) if is_metabolite_id_input else ""
        if met_id:
            aliases.add(met_id)
        aliases.add(variant.replace("_", " "))
        aliases.add(variant.replace("-", " "))

    for concept_id in resolved_concepts:
        concept = concept_index.get(concept_id) or {}
        if not isinstance(concept, dict):
            continue
        label = normalize(concept.get("label") or "")
        if label:
            aliases.add(label)
        for alias in (concept.get("aliases") or [])[:24]:
            alias_clean = normalize(alias)
            if alias_clean:
                aliases.add(alias_clean)
        for met_id in concept.get("ids") or []:
            met_id_clean = canonical_metabolite_id(met_id)
            if met_id_clean:
                aliases.add(met_id_clean)

    expanded_aliases: set[str] = set()
    for alias in aliases:
        if not alias:
            continue
        expanded_aliases.add(alias)
        for variant in subject_phrase_variants(alias):
            if variant:
                expanded_aliases.add(variant)
    return sorted(set(x for x in expanded_aliases if x))


def build_query_aliases(
    query: str,
    index: dict[str, Any],
    resolved_accessions: list[str],
    name_mode: str,
    is_accession_input: bool,
) -> list[str]:
    aliases: set[str] = set()
    accession_index = index.get("accession_alias_index", {})

    # In strict mode, unresolved free-text names are not fuzzily expanded.
    if name_mode == "strict" and not is_accession_input and not resolved_accessions:
        for variant in name_lookup_variants(query):
            aliases.add(variant)
        return [x for x in aliases if x]

    for variant in name_lookup_variants(query):
        aliases.add(variant)
        canonical = canonical_accession(variant) if is_accession_input else ""
        if canonical:
            aliases.add(canonical)
            for alias in accession_index.get(canonical, []):
                for alias_variant in alias_variants(alias):
                    if alias_variant:
                        aliases.add(alias_variant)
            continue
        aliases.add(variant.replace("_", " "))
        aliases.add(variant.replace("-", " "))

    for acc in resolved_accessions:
        aliases.add(acc)
        for alias in accession_index.get(acc, []):
            for alias_variant in alias_variants(alias):
                if alias_variant:
                    aliases.add(alias_variant)
    expanded_aliases: set[str] = set()
    for alias in aliases:
        if not alias:
            continue
        expanded_aliases.add(alias)
        for variant in subject_phrase_variants(alias):
            if variant:
                expanded_aliases.add(variant)
    cleaned_aliases = [x for x in expanded_aliases if x]

    # If we already have informative alias phrases, ignore weak short symbols
    # (for example C2) unless the query explicitly asked for symbol/mnemonic mode.
    query_has_symbol_intent = bool(symbol_query_key(query) or uniprot_mnemonic_symbol_keys(query, index))
    if not query_has_symbol_intent:
        has_informative_alias = any(
            (not canonical_accession(alias))
            and (not is_symbol_like_alias(alias))
            and bool(informative_lookup_tokens(alias))
            for alias in cleaned_aliases
        )
        if has_informative_alias:
            cleaned_aliases = [alias for alias in cleaned_aliases if not is_weak_symbol_alias(alias)]

    return sorted(set(cleaned_aliases))


def candidates_from_index(
    query: str,
    index: dict[str, Any],
    matrix_priority: list[str],
    measurement_context: str,
    additional_contexts: list[str],
    additional_context_keywords: list[str],
    name_mode: str,
    input_type: str = "auto",
    entity_type: str = "auto",
) -> list[Candidate]:
    term_meta: dict[str, dict[str, Any]] = index.get("term_meta", {})
    exact_term_index: dict[str, list[str]] = index.get("exact_term_index", {})
    token_index: dict[str, list[str]] = index.get("token_index", {})
    token_freq: dict[str, int] = index.get("token_freq", {})
    term_metabolite_id_index: dict[str, list[str]] = index.get("term_metabolite_id_index", {})
    analyte_index: dict[str, list[dict[str, Any]]] = index.get("analyte_index", {})

    candidates: dict[str, Candidate] = {}
    input_type_norm = normalize_input_type(input_type)
    resolved_accessions, used_alias_resolution, _ambiguous_alias_resolution = resolve_query_accessions(
        query,
        index,
        input_type=input_type_norm,
    )
    resolved_metabolites, used_met_alias_resolution, ambiguous_met_alias_resolution = resolve_query_metabolite_concepts(
        query, index
    )
    route = detect_entity_route(
        query=query,
        input_type=input_type_norm,
        forced_entity_type=entity_type,
        resolved_accessions=resolved_accessions,
        resolved_metabolites=resolved_metabolites,
    )
    if normalize_entity_type(entity_type) == "metabolite" and query_has_protein_cue(query):
        # Metabolite-only mode is intended for true metabolite rows; fail closed
        # on protein-like analytes so users can rerun in auto/protein mode.
        return []
    is_accession_input = any(canonical_accession(v) for v in query_variants(query))
    is_metabolite_id_input = any(canonical_metabolite_id(v) for v in query_variants(query))
    is_symbol_input = bool(symbol_query_key(query) or uniprot_mnemonic_symbol_keys(query, index))
    accession_profiles = (
        compile_accession_profiles(accession_profiles_for_accessions(resolved_accessions, index))
        if route == "protein"
        else tuple()
    )
    metabolite_subject_profiles = concept_subject_terms(resolved_metabolites, index) if route == "metabolite" else tuple()
    has_accession_profiles = bool(accession_profiles)
    has_metabolite_profiles = bool(metabolite_subject_profiles) and (
        is_metabolite_id_input or not ambiguous_met_alias_resolution
    )
    allow_ratio_terms = query_requests_ratio(query)
    query_apolipoprotein_subtypes = apolipoprotein_subtype_tags(normalize_lipid_panel_phrase(query) or query)
    query_protein_discriminator_tags = (
        protein_discriminator_tags(normalize_lipid_panel_phrase(query) or query)
        if route == "protein"
        else frozenset()
    )
    query_metabolite_descriptor = metabolite_descriptor_tag(query) if route == "metabolite" else ""
    query_has_chain_signature = has_lipid_chain_signature(query) if route == "metabolite" else False
    query_lipoprotein_classes = lipoprotein_class_tags(query) if route == "metabolite" else frozenset()
    query_lipoprotein_sizes = lipoprotein_size_tags(query) if route == "metabolite" else frozenset()
    query_all_lipoprotein_scope = has_all_lipoprotein_scope(query) if route == "metabolite" else False
    query_roundtrip_symbols = list(query_protein_roundtrip_symbols(query, resolved_accessions, index))
    strict_unresolved_name = (
        name_mode == "strict"
        and ((route == "protein" and not is_accession_input) or (route == "metabolite" and not is_metabolite_id_input))
        and not EFO_RE.match(query)
        and not (used_alias_resolution if route == "protein" else used_met_alias_resolution)
    )
    strict_tokens = (
        strict_query_tokens(query)
        if (
            name_mode == "strict"
            and input_type_norm != "gene_id"
            and (
                (route == "protein" and not is_accession_input and not is_symbol_input)
                or (route == "metabolite" and not is_metabolite_id_input)
            )
        )
        else frozenset()
    )
    term_cache: dict[str, tuple[str, list[str], str, frozenset[str], frozenset[str]]] = {}
    identity_cache: dict[str, bool] = {}
    context_cache: dict[tuple[str, str], bool] = {}
    token_gate_cache: dict[str, bool] = {}
    candidate_eval_cache: dict[str, tuple[bool, str, list[str], str, frozenset[str], frozenset[str], int, int, bool]] = {}
    protein_roundtrip_cache: dict[str, tuple[str, tuple[str, ...], tuple[str, ...]]] = {}

    def get_term_bundle(efo_id: str) -> tuple[str, list[str], str, frozenset[str], frozenset[str]]:
        cached_bundle = term_cache.get(efo_id)
        if cached_bundle is not None:
            return cached_bundle
        meta = term_meta.get(efo_id, {})
        label = normalize(meta.get("label") or "")
        syns = [normalize(x) for x in (meta.get("synonyms") or []) if normalize(x)]
        merged = normalize(" ".join([label] + syns)) if label or syns else ""
        nums = numeric_tokens(merged)
        # Use full EFO/OBA synonym-derived subject variants for identity checks.
        subject_terms_set: set[str] = set(extract_subject_terms(label, syns))
        meta_subject = normalize(meta.get("subject") or "")
        for subject in subject_phrase_variants(meta_subject):
            if subject:
                subject_terms_set.add(subject)
        if not subject_terms_set and label:
            subject_terms_set.update(subject_phrase_variants(label))
        subject_terms = frozenset(subject_terms_set)
        bundle = (label, syns, merged, nums, subject_terms)
        term_cache[efo_id] = bundle
        return bundle

    def identity_ok(
        efo_id: str,
        label: str,
        syns: list[str],
        merged: str,
        nums: frozenset[str],
        subject_terms: frozenset[str],
    ) -> bool:
        cached = identity_cache.get(efo_id)
        if cached is not None:
            return cached
        if route == "protein":
            candidate_text = " ".join([label] + syns[:10])
            if not has_accession_profiles:
                # Free-text protein rows without accession resolution are prone
                # to family-member drift in cache/token retrieval. Keep a
                # lightweight subtype/discriminator guard active.
                if query_apolipoprotein_subtypes:
                    candidate_apolipoprotein_subtypes = apolipoprotein_subtype_tags(candidate_text)
                    if (
                        not candidate_apolipoprotein_subtypes
                        or candidate_apolipoprotein_subtypes.isdisjoint(query_apolipoprotein_subtypes)
                    ):
                        identity_cache[efo_id] = False
                        return False
                if query_protein_discriminator_tags:
                    candidate_discriminator_tags = protein_discriminator_tags(candidate_text)
                    if (
                        candidate_discriminator_tags
                        and candidate_discriminator_tags.isdisjoint(query_protein_discriminator_tags)
                    ):
                        identity_cache[efo_id] = False
                        return False
                identity_cache[efo_id] = True
                return True

            ok = accession_identity_match(
                query,
                label,
                syns,
                index,
                precomputed_profiles=accession_profiles,
                candidate_text=merged,
                candidate_nums=nums,
                candidate_subject_terms=subject_terms,
            )
        else:
            if not has_metabolite_profiles:
                identity_cache[efo_id] = True
                return True
            ok = metabolite_identity_match(
                query=query,
                candidate_subject_terms=subject_terms,
                candidate_nums=nums,
                concept_subject_profiles=metabolite_subject_profiles,
            )
        identity_cache[efo_id] = ok
        return ok

    def profile_exact_bonus(subject_terms: frozenset[str]) -> int:
        if route != "protein" or not has_accession_profiles:
            return 0
        return profile_subject_exact_match_bonus(accession_profiles, subject_terms)

    def profile_exact_rank(subject_terms: frozenset[str]) -> int:
        if route != "protein" or not has_accession_profiles:
            return 0
        return profile_subject_exact_match_rank(accession_profiles, subject_terms)

    def context_ok(efo_id: str, label: str, syns: list[str]) -> bool:
        cache_key = (efo_id, label)
        cached = context_cache.get(cache_key)
        if cached is not None:
            return cached
        ok = context_compatible(
            label,
            syns,
            measurement_context,
            additional_contexts,
            additional_context_keywords,
        )
        context_cache[cache_key] = ok
        return ok

    def token_gate_ok(efo_id: str, label: str, syns: list[str]) -> bool:
        if not strict_tokens:
            return True
        cached = token_gate_cache.get(efo_id)
        if cached is not None:
            return cached
        candidate_tokens = tokenize(normalize_measurement_subject(label))
        ok = bool(candidate_tokens & strict_tokens)
        token_gate_cache[efo_id] = ok
        return ok

    def candidate_eval(
        efo_id: str,
    ) -> tuple[bool, str, list[str], str, frozenset[str], frozenset[str], int, int, bool]:
        cached_eval = candidate_eval_cache.get(efo_id)
        if cached_eval is not None:
            return cached_eval

        label, syns, merged, nums, subject_terms = get_term_bundle(efo_id)
        ok = bool(label)
        if ok and not allow_ratio_terms and is_ratio_trait(label, syns):
            ok = False
        if ok and not token_gate_ok(efo_id, label, syns):
            ok = False
        if ok and not context_ok(efo_id, label, syns):
            ok = False
        if ok and not identity_ok(efo_id, label, syns, merged, nums, subject_terms):
            ok = False
        if ok and query_apolipoprotein_subtypes:
            candidate_apolipoprotein_subtypes = apolipoprotein_subtype_tags(" ".join([label] + syns[:10]))
            if (
                not candidate_apolipoprotein_subtypes
                or candidate_apolipoprotein_subtypes.isdisjoint(query_apolipoprotein_subtypes)
            ):
                # Enforce subtype consistency for Apo-family queries regardless
                # of protein/metabolite routing.
                ok = False
        if ok and route == "metabolite":
            candidate_text = " ".join([label] + syns[:10])
            candidate_descriptor = metabolite_descriptor_tag(candidate_text)
            candidate_classes = lipoprotein_class_tags(candidate_text)
            candidate_sizes = lipoprotein_size_tags(candidate_text)
            if query_metabolite_descriptor:
                if query_metabolite_descriptor == "cholesterol":
                    if candidate_descriptor != "cholesterol":
                        ok = False
                elif candidate_descriptor != query_metabolite_descriptor:
                    ok = False
            if ok and not query_lipoprotein_classes and candidate_classes:
                # For unqualified metabolite queries (for example "cholesterol levels"),
                # avoid auto-promoting lipoprotein subclass terms by default.
                ok = False
            if ok and query_lipoprotein_classes:
                if not candidate_classes or candidate_classes.isdisjoint(query_lipoprotein_classes):
                    ok = False
            if ok and query_lipoprotein_classes:
                if query_all_lipoprotein_scope and candidate_sizes:
                    ok = False
                if query_lipoprotein_sizes:
                    # If query asks for a specific size, allow exact size matches.
                    # When no size exists on the candidate (class-level term), keep it
                    # as a fallback rather than dropping the class-consistent term.
                    if candidate_sizes and candidate_sizes.isdisjoint(query_lipoprotein_sizes):
                        ok = False
            if ok and has_metabolite_profiles and not query_has_chain_signature and has_lipid_chain_signature(merged):
                # Prevent generic analyte names from auto-mapping to specific
                # chain-composition terms (for example CE 22:6).
                ok = False

        exact_bonus = profile_exact_bonus(subject_terms) if ok else 0
        exact_rank = profile_exact_rank(subject_terms) if ok else 0
        measurement_like = is_measurement_like(label, syns) if ok else False
        result = (ok, label, syns, merged, nums, subject_terms, exact_bonus, exact_rank, measurement_like)
        candidate_eval_cache[efo_id] = result
        return result

    def candidate_protein_roundtrip(
        efo_id: str,
        label: str,
        syns: list[str],
        subject_terms: frozenset[str],
    ) -> tuple[str, tuple[str, ...], tuple[str, ...]]:
        if route != "protein" or len(resolved_accessions) != 1:
            return ("not_applicable", (), ())
        cached = protein_roundtrip_cache.get(efo_id)
        if cached is not None:
            return cached
        candidate_accessions = resolve_term_protein_accessions(
            efo_id,
            index,
            label=label,
            synonyms=syns,
            subject_terms=subject_terms,
        )
        candidate_symbols = resolve_term_protein_symbols(
            efo_id,
            index,
            label=label,
            synonyms=syns,
            subject_terms=subject_terms,
        )
        roundtrip_status = protein_roundtrip_status(
            resolved_accessions,
            query_roundtrip_symbols,
            candidate_accessions,
            candidate_symbols,
        )
        result = (roundtrip_status, candidate_accessions, candidate_symbols)
        protein_roundtrip_cache[efo_id] = result
        return result

    def finalize_candidate(
        candidate: Candidate,
        *,
        label: str,
        syns: list[str],
        subject_terms: frozenset[str],
    ) -> Candidate:
        roundtrip_status, roundtrip_accessions, roundtrip_symbols = candidate_protein_roundtrip(
            candidate.efo_id,
            label,
            syns,
            subject_terms,
        )
        matched_via = candidate.matched_via
        evidence = candidate.evidence
        is_validated = candidate.is_validated
        if roundtrip_status != "not_applicable":
            matched_via = f"{matched_via}-roundtrip-{roundtrip_status.replace('_', '-')}"
            roundtrip_evidence = protein_accession_roundtrip_evidence(
                roundtrip_status,
                roundtrip_accessions,
                roundtrip_symbols,
            )
            if roundtrip_evidence:
                evidence = f"{evidence}; {roundtrip_evidence}" if evidence else roundtrip_evidence
            is_validated = is_validated and roundtrip_status in {"accession_match", "symbol_match"}
        return Candidate(
            efo_id=candidate.efo_id,
            label=candidate.label,
            score=candidate.score,
            matched_via=matched_via,
            evidence=evidence,
            is_validated=is_validated,
            roundtrip_status=roundtrip_status,
            roundtrip_accessions=roundtrip_accessions,
            roundtrip_symbols=roundtrip_symbols,
        )

    def select_retrieval_tokens(tokens: set[str], *, max_tokens: int = 6) -> list[str]:
        ranked = [tok for tok in tokens if tok in token_index]
        if not ranked:
            return []
        ranked.sort(key=lambda tok: (token_freq.get(tok, 10**9), len(tok), tok))
        # De-prioritize very high-frequency tokens when possible.
        selected = [tok for tok in ranked if token_freq.get(tok, 10**9) <= 2500][:max_tokens]
        if not selected:
            selected = ranked[:max_tokens]
        return selected

    if EFO_RE.match(query):
        efo_id = normalize_id(query)
        meta = term_meta.get(efo_id, {})
        label, syns, _, _, _ = get_term_bundle(efo_id)
        if label and not context_ok(efo_id, label, syns):
            return []
        direct = [
            finalize_candidate(
                Candidate(
                    efo_id=efo_id,
                    label=label,
                    score=1.0,
                    matched_via="direct-id",
                    evidence="direct input EFO/OBA id",
                    is_validated=bool(meta),
                ),
                label=label,
                syns=syns,
                subject_terms=extract_subject_terms(label, syns),
            )
        ]
        return sorted(
            direct,
            key=lambda c: (
                c.score + matrix_bonus(c.label, [], matrix_priority),
                matrix_preference_rank(c.label, matrix_priority),
                c.efo_id,
            ),
            reverse=True,
        )

    query_keys = {norm_key(v) for v in query_variants(query)}
    for query_key in query_keys:
        for cached in analyte_index.get(query_key, []):
            efo_id = normalize_id(cached.get("efo_id") or "")
            meta = term_meta.get(efo_id, {})
            ok, label, syns, _merged, _nums, subject_terms, _exact_bonus, _exact_rank, _measurement_like = candidate_eval(
                efo_id
            )
            if not ok:
                continue
            validated = cached.get("validation") == "validated" or bool(meta)
            cached_source = norm_key(cached.get("source") or "")
            try:
                cached_score = float(cached.get("confidence") or 0.8)
            except ValueError:
                cached_score = 0.8
            # Prevent stale cache entries from outranking fresh exact-term hits.
            if "manual" in cached_source or "curated" in cached_source:
                score = min(1.0, max(0.0, cached_score))
            else:
                score = min(0.89, max(0.0, cached_score))
            merge_candidate(
                candidates,
                finalize_candidate(
                    Candidate(
                        efo_id=efo_id,
                        label=label,
                        score=score,
                        matched_via="local-analyte-cache",
                        evidence=normalize(cached.get("evidence") or "cached mapping"),
                        is_validated=validated,
                    ),
                    label=label,
                    syns=syns,
                    subject_terms=subject_terms,
                ),
            )

    # For metabolite route, prefer direct external-ID joins (CHEBI/HMDB/KEGG)
    # from term cache metadata when available, then let lexical reranking handle
    # remaining unresolved cases.
    if route == "metabolite":
        direct_metabolite_ids: set[str] = set()
        for variant in query_variants(query):
            met_id = canonical_metabolite_id(variant)
            if met_id:
                direct_metabolite_ids.add(met_id)
        concept_index = index.get("metabolite_concept_index", {})
        for concept_id in resolved_metabolites:
            concept_meta = concept_index.get(concept_id) or {}
            if not isinstance(concept_meta, dict):
                continue
            for met_id in concept_meta.get("ids") or []:
                canonical = canonical_metabolite_id(met_id)
                if canonical:
                    direct_metabolite_ids.add(canonical)

        for met_id in sorted(direct_metabolite_ids):
            for efo_id in term_metabolite_id_index.get(met_id, []):
                label, syns, _merged, _nums, subject_terms = get_term_bundle(efo_id)
                if not label:
                    continue
                if not allow_ratio_terms and is_ratio_trait(label, syns):
                    continue
                if not context_ok(efo_id, label, syns):
                    continue
                exact_rank = profile_exact_rank(subject_terms)
                score = min(1.0, max(0.92, 0.97 + matrix_bonus(label, syns, matrix_priority)))
                merge_candidate(
                    candidates,
                    finalize_candidate(
                        Candidate(
                            efo_id=efo_id,
                            label=label,
                            score=score,
                            matched_via="local-term-metabolite-xref" + ("-subject-exact" if exact_rank else ""),
                            evidence=f"metabolite external-id match: {met_id}",
                            is_validated=True,
                        ),
                        label=label,
                        syns=syns,
                        subject_terms=subject_terms,
                    ),
                )

    if strict_unresolved_name:
        # For unresolved free-text names, avoid fuzzy lexical mapping to prevent
        # semantic drift (e.g., matching wrong family members by shared tokens).
        return sorted(
            candidates.values(),
            key=lambda x: (
                x.score,
                matrix_preference_rank(x.label, matrix_priority),
                x.efo_id,
            ),
            reverse=True,
        )

    if route == "metabolite":
        expanded = build_metabolite_query_aliases(
            query,
            index,
            resolved_concepts=resolved_metabolites,
            name_mode=name_mode,
            is_metabolite_id_input=is_metabolite_id_input,
        )
    else:
        expanded = build_query_aliases(
            query,
            index,
            resolved_accessions=resolved_accessions,
            name_mode=name_mode,
            is_accession_input=is_accession_input,
        )
    term_parts: list[tuple[str, str, str, set[str], set[str]]] = []
    for term in expanded:
        term_key = norm_key(term)
        term_lower = term.lower()
        term_toks = tokenize(term)
        filtered_tokens = {tok for tok in term_toks if tok not in GENERIC_TOKENS}
        term_parts.append((term, term_key, term_lower, term_toks, filtered_tokens))
    seen_ids: set[str] = set()

    for term, key, term_lower, term_toks, filtered_tokens in term_parts:
        # Exact normalized label/synonym hits
        for efo_id in exact_term_index.get(key, []):
            ok, label, syns, _merged, _nums, subject_terms, exact_bonus, exact_rank, _measurement_like = candidate_eval(
                efo_id
            )
            if not ok:
                continue
            score = max(0.9, similarity_score_prepared(term_lower, term_toks, label, syns))
            score = min(1.0, max(0.0, score + matrix_bonus(label, syns, matrix_priority)))
            score = min(1.0, max(0.0, score - condition_mismatch_penalty(query, label, syns)))
            score = min(1.0, score + 0.02 * exact_bonus)
            merge_candidate(
                candidates,
                finalize_candidate(
                    Candidate(
                        efo_id=efo_id,
                        label=label,
                        score=score,
                        matched_via="local-term-exact" + ("-subject-exact" if exact_rank else ""),
                        evidence=f"exact normalized term match: {term}",
                        is_validated=True,
                    ),
                    label=label,
                    syns=syns,
                    subject_terms=subject_terms,
                ),
            )
            seen_ids.add(efo_id)

        # Token-based retrieval + lexical rerank
        retrieval_tokens = select_retrieval_tokens(filtered_tokens)
        candidate_ids: set[str] = set()
        for tok in retrieval_tokens:
            candidate_ids.update(token_index.get(tok, []))
        if len(candidate_ids) > 4000 and len(retrieval_tokens) >= 2:
            counts: dict[str, int] = {}
            for tok in retrieval_tokens[: min(4, len(retrieval_tokens))]:
                for efo_id in token_index.get(tok, []):
                    counts[efo_id] = counts.get(efo_id, 0) + 1
            narrowed = {efo_id for efo_id, c in counts.items() if c >= 2}
            if narrowed:
                candidate_ids = narrowed
            elif counts:
                candidate_ids = set(
                    sorted(
                        counts,
                        key=lambda efo_id: (-counts[efo_id], efo_id),
                    )[:4000]
                )

        for efo_id in candidate_ids:
            if efo_id in seen_ids:
                continue
            ok, label, syns, _merged, _nums, subject_terms, exact_bonus, exact_rank, measurement_like = candidate_eval(
                efo_id
            )
            if not ok:
                continue
            score = similarity_score_prepared(term_lower, term_toks, label, syns)
            if not measurement_like:
                score -= 0.35
            score += matrix_bonus(label, syns, matrix_priority)
            score -= condition_mismatch_penalty(query, label, syns)
            score += 0.02 * exact_bonus
            score = min(1.0, max(0.0, score))
            if score < 0.30:
                continue
            semantic_exact = False
            needs_new_term = False
            if route == "metabolite":
                candidate_text = " ".join([label] + syns[:10])
                candidate_descriptor = metabolite_descriptor_tag(candidate_text)
                candidate_classes = lipoprotein_class_tags(candidate_text)
                candidate_sizes = lipoprotein_size_tags(candidate_text)

                descriptor_ok = True
                if query_metabolite_descriptor:
                    if query_metabolite_descriptor == "cholesterol":
                        descriptor_ok = candidate_descriptor == "cholesterol"
                    else:
                        descriptor_ok = candidate_descriptor == query_metabolite_descriptor

                if query_lipoprotein_classes:
                    class_ok = bool(candidate_classes and not candidate_classes.isdisjoint(query_lipoprotein_classes))
                else:
                    class_ok = not candidate_classes

                size_ok = True
                if query_all_lipoprotein_scope and candidate_sizes:
                    size_ok = False
                if query_lipoprotein_sizes:
                    # Accept class-level fallback when candidate has no explicit size.
                    size_ok = (not candidate_sizes) or bool(candidate_sizes & query_lipoprotein_sizes)

                if descriptor_ok and class_ok and size_ok:
                    size_fallback_without_size = bool(
                        query_lipoprotein_sizes and query_lipoprotein_classes and not candidate_sizes
                    )
                    # Keep semantic auto-validation conservative but usable for
                    # normalized metabolite/lipoprotein phrasing variants.
                    if size_fallback_without_size:
                        # Missing size-specific ontology term; allow class-level fallback
                        # when descriptor/class semantics are already aligned.
                        min_semantic_score = 0.60
                        needs_new_term = True
                    else:
                        min_semantic_score = 0.85 if query_lipoprotein_classes else 0.90
                    semantic_exact = score >= min_semantic_score
            merge_candidate(
                candidates,
                finalize_candidate(
                    Candidate(
                        efo_id=efo_id,
                        label=label,
                        score=score,
                        matched_via=(
                            "local-term-token"
                            + ("-subject-exact" if exact_rank else "")
                            + ("-concept-validated" if route == "metabolite" and has_metabolite_profiles else "")
                            + ("-semantic-exact" if semantic_exact else "")
                            + (f"-{NEEDS_NEW_TERM_MATCH_TAG}" if needs_new_term else "")
                        ),
                        evidence=f"token retrieval + lexical score using: {term}",
                        is_validated=True,
                    ),
                    label=label,
                    syns=syns,
                    subject_terms=subject_terms,
                ),
            )

    return sorted(
        candidates.values(),
        key=lambda x: (
            x.score,
            1 if "subject-exact" in x.matched_via else 0,
            matrix_preference_rank(x.label, matrix_priority),
            measurement_label_rank(x.label),
            x.efo_id,
        ),
        reverse=True,
    )


def map_one(
    query: str,
    index: dict[str, Any],
    top_k: int,
    min_score: float,
    matrix_priority: list[str],
    measurement_context: str,
    additional_contexts: list[str],
    additional_context_keywords: list[str],
    name_mode: str,
    force_map_best: bool,
    fallback_efo_id: str,
    input_type: str = "auto",
    entity_type: str = "auto",
) -> list[dict[str, str]]:
    query = normalize(query)
    query_apolipoprotein_subtypes = apolipoprotein_subtype_tags(normalize_lipid_panel_phrase(query) or query)
    input_type_norm = normalize_input_type(input_type)
    resolved_accessions, _used_alias_resolution, _ambiguous_alias_resolution = resolve_query_accessions(
        query,
        index,
        input_type=input_type_norm,
    )
    resolved_metabolites, _used_met_alias_resolution, _ambiguous_met_alias_resolution = resolve_query_metabolite_concepts(
        query, index
    )
    route = detect_entity_route(
        query=query,
        input_type=input_type_norm,
        forced_entity_type=entity_type,
        resolved_accessions=resolved_accessions,
        resolved_metabolites=resolved_metabolites,
    )
    has_external_identity = has_useful_external_identity(
        query=query,
        route=route,
        resolved_accessions=resolved_accessions,
        resolved_metabolites=resolved_metabolites,
        index=index,
    )
    candidates = candidates_from_index(
        query,
        index,
        matrix_priority=matrix_priority,
        measurement_context=measurement_context,
        additional_contexts=additional_contexts,
        additional_context_keywords=additional_context_keywords,
        name_mode=name_mode,
        input_type=input_type_norm,
        entity_type=entity_type,
    )
    eligible = [c for c in candidates if c.score >= min_score]

    def apolipoprotein_subtype_consistent(cand: Candidate) -> bool:
        if not query_apolipoprotein_subtypes:
            return True
        candidate_tags = apolipoprotein_subtype_tags(cand.label)
        return bool(candidate_tags and (candidate_tags & query_apolipoprotein_subtypes))

    def auto_validates(cand: Candidate) -> bool:
        if not cand.is_validated:
            return False
        if not apolipoprotein_subtype_consistent(cand):
            return False
        if NEEDS_NEW_TERM_MATCH_TAG in cand.matched_via:
            return False
        if cand.matched_via.startswith("local-term-token"):
            has_subject_exact = "subject-exact" in cand.matched_via
            has_concept_validation = "concept-validated" in cand.matched_via
            has_semantic_exact = "semantic-exact" in cand.matched_via
            if not has_subject_exact and not has_concept_validation and not has_semantic_exact:
                return False
            if has_semantic_exact:
                if cand.score < TOKEN_SEMANTIC_EXACT_AUTO_VALIDATE_MIN_SCORE:
                    return False
            elif cand.score < TOKEN_SUBJECT_EXACT_AUTO_VALIDATE_MIN_SCORE:
                return False
        return True

    selected = [c for c in eligible if auto_validates(c)][:top_k]
    withheld = [c for c in eligible if not auto_validates(c)]

    if not selected:
        if force_map_best and candidates:
            best = candidates[0]
            return [
                {
                    "input_query": query,
                    "mapped_efo_id": format_ontology_id_for_output(best.efo_id),
                    "mapped_label": best.label,
                    "confidence": f"{best.score:.3f}",
                    "matched_via": f"{best.matched_via}-forced",
                    "validation": "not_validated",
                    "input_type": input_type_norm,
                    "evidence": (
                        f"forced-best below threshold ({min_score:.2f}); "
                        f"{best.evidence}"
                    ),
                }
            ]
        if withheld:
            best = withheld[0]
            if NEEDS_NEW_TERM_MATCH_TAG in best.matched_via:
                if not has_external_identity:
                    return [
                        {
                            "input_query": query,
                            "mapped_efo_id": "",
                            "mapped_label": "",
                            "confidence": "0.0",
                            "matched_via": "withheld-broad-term-suggestion",
                            "validation": "not_mapped",
                            "input_type": input_type_norm,
                            "evidence": (
                                f"reason_code=suggest_broad_term_mapping; candidate "
                                f"{format_ontology_id_for_output(best.efo_id)} ({best.label}) "
                                "is the broad-term fallback, but no stable external identifier "
                                "was available to justify new-term request"
                            ),
                        }
                    ]
                return [
                    {
                        "input_query": query,
                        "mapped_efo_id": "",
                        "mapped_label": "",
                        "confidence": "0.0",
                        "matched_via": "withheld-needs-new-term",
                        "validation": "not_mapped",
                        "input_type": input_type_norm,
                        "evidence": (
                            f"reason_code=needs_new_term; candidate {format_ontology_id_for_output(best.efo_id)} "
                            f"is a class-level fallback for a size-qualified query (score={best.score:.3f}, "
                            f"matched_via={best.matched_via})"
                        ),
                    }
                ]
            return [
                {
                    "input_query": query,
                    "mapped_efo_id": "",
                    "mapped_label": "",
                    "confidence": "0.0",
                    "matched_via": "withheld-for-review",
                    "validation": "not_mapped",
                    "input_type": input_type_norm,
                    "evidence": (
                        f"candidate withheld from auto-validation: {format_ontology_id_for_output(best.efo_id)} "
                        f"(score={best.score:.3f}, matched_via={best.matched_via})"
                    ),
                }
            ]
        if normalize_entity_type(entity_type) == "metabolite" and query_has_protein_cue(query):
            return [
                {
                    "input_query": query,
                    "mapped_efo_id": "",
                    "mapped_label": "",
                    "confidence": "0.0",
                    "matched_via": "blocked-protein-like-query",
                    "validation": "not_mapped",
                    "input_type": input_type_norm,
                    "evidence": "metabolite-only mode blocked protein-like query; rerun with --entity-type auto",
                }
            ]
        if force_map_best and fallback_efo_id:
            fallback_id = normalize_id(fallback_efo_id)
            term_meta = index.get("term_meta", {})
            fallback_label = normalize((term_meta.get(fallback_id, {}) or {}).get("label") or "measurement")
            return [
                {
                    "input_query": query,
                    "mapped_efo_id": format_ontology_id_for_output(fallback_id),
                    "mapped_label": fallback_label,
                    "confidence": "0.000",
                    "matched_via": "fallback-generic",
                    "validation": "not_validated",
                    "input_type": input_type_norm,
                    "evidence": "no specific candidate found; generic fallback applied",
                }
            ]
        return [
            {
                "input_query": query,
                "mapped_efo_id": "",
                "mapped_label": "",
                "confidence": "0.0",
                "matched_via": "none",
                "validation": "not_mapped",
                "input_type": input_type_norm,
                "evidence": "no candidate above threshold from local index",
            }
        ]

    out: list[dict[str, str]] = []
    for cand in selected:
        out.append(
            {
                "input_query": query,
                "mapped_efo_id": format_ontology_id_for_output(cand.efo_id),
                "mapped_label": cand.label,
                "confidence": f"{cand.score:.3f}",
                "matched_via": cand.matched_via,
                "validation": "validated" if cand.is_validated else "not_validated",
                "input_type": input_type_norm,
                "evidence": cand.evidence,
            }
        )
    return out


def map_queries(
    query_inputs: list[tuple[str, str]],
    index: dict[str, Any],
    top_k: int,
    min_score: float,
    workers: int,
    matrix_priority: list[str],
    measurement_context: str,
    additional_contexts: list[str],
    additional_context_keywords: list[str],
    name_mode: str,
    force_map_best: bool,
    fallback_efo_id: str,
    entity_type: str,
    show_progress: bool,
    parallel_mode: str,
    index_path: Path | None = None,
) -> list[dict[str, str]]:
    normalized_inputs = [
        (normalize(query), normalize_input_type(input_type))
        for query, input_type in query_inputs
        if normalize(query)
    ]
    progress = ProgressReporter(total=len(normalized_inputs), enabled=show_progress)
    if not normalized_inputs:
        return []

    unique_inputs = list(dict.fromkeys(normalized_inputs))
    input_counts = Counter(normalized_inputs)

    def run_one(item: tuple[str, str]) -> list[dict[str, str]]:
        query, input_type = item
        return map_one(
            query,
            index,
            top_k=top_k,
            min_score=min_score,
            matrix_priority=matrix_priority,
            measurement_context=measurement_context,
            additional_contexts=additional_contexts,
            additional_context_keywords=additional_context_keywords,
            name_mode=effective_name_mode(name_mode, input_type),
            force_map_best=force_map_best,
            fallback_efo_id=fallback_efo_id,
            input_type=input_type,
            entity_type=entity_type,
        )

    if workers <= 1:
        mapped_by_input: dict[tuple[str, str], list[dict[str, str]]] = {}
        for item in unique_inputs:
            mapped_by_input[item] = run_one(item)
            progress.update(input_counts[item])
        rows: list[dict[str, str]] = []
        for item in normalized_inputs:
            rows.extend(mapped_by_input[item])
        return rows

    selected_mode = parallel_mode
    if selected_mode == "auto":
        # Process pools usually win for large batches; keep threads for smaller jobs.
        selected_mode = "process" if len(unique_inputs) >= 200 else "thread"

    if selected_mode == "process" and index_path is not None:
        try:
            mapped_by_input: dict[tuple[str, str], list[dict[str, str]]] = {}
            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=init_process_mapper,
                initargs=(
                    str(index_path),
                    top_k,
                    min_score,
                    matrix_priority,
                    measurement_context,
                    additional_contexts,
                    additional_context_keywords,
                    name_mode,
                    force_map_best,
                    fallback_efo_id,
                    entity_type,
                ),
            ) as pool:
                future_to_query = {pool.submit(process_map_one, item): item for item in unique_inputs}
                for future in as_completed(future_to_query):
                    item = future_to_query[future]
                    mapped_by_input[item] = future.result()
                    progress.update(input_counts[item])

            rows: list[dict[str, str]] = []
            for item in normalized_inputs:
                rows.extend(mapped_by_input[item])
            return rows
        except Exception as exc:
            print(
                f"[WARN] process parallel mode unavailable ({exc}); falling back to thread mode",
                file=sys.stderr,
            )

    mapped_by_input: dict[tuple[str, str], list[dict[str, str]]] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_to_input = {
            pool.submit(
                run_one,
                item,
            ): item
            for item in unique_inputs
        }
        for future in as_completed(future_to_input):
            item = future_to_input[future]
            mapped_by_input[item] = future.result()
            progress.update(input_counts[item])

    rows: list[dict[str, str]] = []
    for item in normalized_inputs:
        rows.extend(mapped_by_input[item])
    return rows


def lexical_probe_candidates(
    query: str,
    index: dict[str, Any],
    matrix_priority: list[str],
    measurement_context: str,
    additional_contexts: list[str],
    additional_context_keywords: list[str],
    name_mode: str,
    input_type: str = "auto",
    entity_type: str = "auto",
) -> list[Candidate]:
    term_meta: dict[str, dict[str, Any]] = index.get("term_meta", {})
    exact_term_index: dict[str, list[str]] = index.get("exact_term_index", {})
    token_index: dict[str, list[str]] = index.get("token_index", {})
    candidates: dict[str, Candidate] = {}
    resolved_accessions, _used_alias_resolution, _ambiguous_alias_resolution = resolve_query_accessions(
        query,
        index,
        input_type=input_type,
    )
    resolved_metabolites, _used_met_alias_resolution, _ambiguous_met_alias_resolution = resolve_query_metabolite_concepts(
        query, index
    )
    route = detect_entity_route(
        query=query,
        input_type=normalize_input_type(input_type),
        forced_entity_type=entity_type,
        resolved_accessions=resolved_accessions,
        resolved_metabolites=resolved_metabolites,
    )
    is_accession_input = any(canonical_accession(v) for v in query_variants(query))
    is_metabolite_id_input = any(canonical_metabolite_id(v) for v in query_variants(query))
    allow_ratio_terms = query_requests_ratio(query)
    if route == "metabolite":
        expanded = build_metabolite_query_aliases(
            query,
            index,
            resolved_concepts=resolved_metabolites,
            name_mode=name_mode,
            is_metabolite_id_input=is_metabolite_id_input,
        )
    else:
        expanded = build_query_aliases(
            query,
            index,
            resolved_accessions=resolved_accessions,
            name_mode=name_mode,
            is_accession_input=is_accession_input,
        )
    term_parts: list[tuple[str, str, set[str], set[str]]] = []
    for term in expanded:
        term_lower = term.lower()
        term_toks = tokenize(term)
        filtered_tokens = {tok for tok in term_toks if tok not in GENERIC_TOKENS}
        term_parts.append((term, term_lower, term_toks, filtered_tokens))

    for term, term_lower, term_toks, filtered_tokens in term_parts:
        key = norm_key(term)
        candidate_ids: set[str] = set(exact_term_index.get(key, []))
        for tok in filtered_tokens:
            candidate_ids.update(token_index.get(tok, []))

        for efo_id in candidate_ids:
            meta = term_meta.get(efo_id, {})
            label = normalize(meta.get("label") or "")
            syns = [normalize(x) for x in (meta.get("synonyms") or []) if normalize(x)]
            if not label:
                continue
            if not allow_ratio_terms and is_ratio_trait(label, syns):
                continue
            if not context_compatible(
                label,
                syns,
                measurement_context,
                additional_contexts,
                additional_context_keywords,
            ):
                continue
            score = similarity_score_prepared(term_lower, term_toks, label, syns)
            if not is_measurement_like(label, syns):
                score -= 0.35
            score += matrix_bonus(label, syns, matrix_priority)
            score = min(1.0, max(0.0, score))
            if score < 0.20:
                continue
            merge_candidate(
                candidates,
                Candidate(
                    efo_id=efo_id,
                    label=label,
                    score=score,
                    matched_via="lexical-probe",
                    evidence=f"lexical probe using: {term}",
                    is_validated=True,
                ),
            )

    return sorted(
        candidates.values(),
        key=lambda x: (
            x.score,
            matrix_preference_rank(x.label, matrix_priority),
            measurement_label_rank(x.label),
            x.efo_id,
        ),
        reverse=True,
    )


def infer_review_reason(
    query: str,
    index: dict[str, Any],
    min_score: float,
    matrix_priority: list[str],
    measurement_context: str,
    additional_contexts: list[str],
    additional_context_keywords: list[str],
    name_mode: str,
    input_type: str = "auto",
    entity_type: str = "auto",
) -> tuple[str, str, list[str], list[Candidate]]:
    query = normalize(query)
    effective_mode = effective_name_mode(name_mode, input_type)
    resolved_accessions, used_alias_resolution, ambiguous_alias_resolution = resolve_query_accessions(
        query,
        index,
        input_type=input_type,
    )
    resolved_metabolites, used_met_alias_resolution, ambiguous_met_alias_resolution = resolve_query_metabolite_concepts(
        query, index
    )
    route = detect_entity_route(
        query=query,
        input_type=normalize_input_type(input_type),
        forced_entity_type=entity_type,
        resolved_accessions=resolved_accessions,
        resolved_metabolites=resolved_metabolites,
    )
    if normalize_entity_type(entity_type) == "metabolite" and query_has_protein_cue(query):
        return (
            "metabolite_mode_blocked_protein_like_query",
            "metabolite-only mode blocked a protein-like analyte name; rerun in auto mode",
            resolved_accessions,
            [],
        )
    resolved_identity = resolved_accessions if route == "protein" else resolved_metabolites
    has_external_identity = has_useful_external_identity(
        query=query,
        route=route,
        resolved_accessions=resolved_accessions,
        resolved_metabolites=resolved_metabolites,
        index=index,
    )
    is_accession_input = any(canonical_accession(v) for v in query_variants(query))
    is_metabolite_input = any(canonical_metabolite_id(v) for v in query_variants(query))
    used_resolution = used_alias_resolution if route == "protein" else used_met_alias_resolution
    ambiguous_resolution = ambiguous_alias_resolution if route == "protein" else ambiguous_met_alias_resolution

    strict_candidates = candidates_from_index(
        query,
        index,
        matrix_priority=matrix_priority,
        measurement_context=measurement_context,
        additional_contexts=additional_contexts,
        additional_context_keywords=additional_context_keywords,
        name_mode=effective_mode,
        input_type=input_type,
        entity_type=route,
    )

    if strict_candidates:
        best = strict_candidates[0]
        if best.score >= min_score:
            if NEEDS_NEW_TERM_MATCH_TAG in best.matched_via:
                if not has_external_identity:
                    return (
                        "suggest_broad_term_mapping",
                        (
                            f"best candidate {format_ontology_id_for_output(best.efo_id)} scored {best.score:.3f} but "
                            "no stable external identifier was available; prefer broad-term mapping"
                        ),
                        resolved_identity,
                        strict_candidates,
                    )
                return (
                    "needs_new_term",
                    (
                        f"best candidate {format_ontology_id_for_output(best.efo_id)} scored {best.score:.3f} but "
                        "drops required size granularity; request a new ontology term"
                    ),
                    resolved_identity,
                    strict_candidates,
                )
            return (
                "manual_validation_required",
                (
                    f"best candidate {format_ontology_id_for_output(best.efo_id)} scored {best.score:.3f} but did not pass "
                    "auto-validation gate"
                ),
                resolved_identity,
                strict_candidates,
            )
        return (
            "below_score_threshold",
            (
                f"best candidate {format_ontology_id_for_output(best.efo_id)} "
                f"scored {best.score:.3f}, below min_score={min_score:.3f}"
            ),
            resolved_identity,
            strict_candidates,
        )

    if measurement_context != "auto":
        any_context_candidates = candidates_from_index(
            query,
            index,
            matrix_priority=matrix_priority,
            measurement_context="auto",
            additional_contexts=[],
            additional_context_keywords=[],
            name_mode=effective_mode,
            input_type=input_type,
            entity_type=route,
        )
        if any_context_candidates:
            best = any_context_candidates[0]
            return (
                "context_mismatch",
                (
                    f"candidate {format_ontology_id_for_output(best.efo_id)} exists only outside "
                    f"requested context '{measurement_context}'"
                ),
                resolved_identity,
                any_context_candidates,
            )

    if (
        effective_mode == "strict"
        and not EFO_RE.match(query)
        and (
            (route == "protein" and not is_accession_input)
            or (route == "metabolite" and not is_metabolite_input)
        )
        and not used_resolution
    ):
        unresolved_msg = (
            "query did not resolve to a unique accession in strict mode"
            if route == "protein"
            else "query did not resolve to a unique metabolite concept in strict mode"
        )
        return (
            "unresolved_name_strict_mode",
            unresolved_msg,
            resolved_identity,
            [],
        )

    probe_candidates = lexical_probe_candidates(
        query,
        index,
        matrix_priority=matrix_priority,
        measurement_context=measurement_context,
        additional_contexts=additional_contexts,
        additional_context_keywords=additional_context_keywords,
        name_mode=effective_mode,
        input_type=input_type,
        entity_type=route,
    )
    if route == "protein" and resolved_accessions and probe_candidates:
        best = probe_candidates[0]
        return (
            "identity_conflict_family_member",
            (
                f"lexical candidate {format_ontology_id_for_output(best.efo_id)} exists but "
                "failed accession identity validation"
            ),
            resolved_accessions,
            probe_candidates,
        )
    if route == "metabolite" and resolved_metabolites and probe_candidates:
        best = probe_candidates[0]
        return (
            "identity_conflict_metabolite",
            (
                f"lexical candidate {format_ontology_id_for_output(best.efo_id)} exists but "
                "failed metabolite concept identity validation"
            ),
            resolved_metabolites,
            probe_candidates,
        )

    if ambiguous_resolution and not resolved_identity:
        return (
            "ambiguous_alias_resolution" if route == "protein" else "ambiguous_metabolite_alias_resolution",
            (
                "query alias matched multiple accessions without a stable winner"
                if route == "protein"
                else "query alias matched multiple metabolite concepts without a stable winner"
            ),
            resolved_identity,
            probe_candidates,
        )

    if route == "protein" and resolved_accessions:
        return (
            "no_measurement_term_for_resolved_accession",
            "resolved accession(s) found but no valid measurement term passed filters",
            resolved_accessions,
            probe_candidates,
        )
    if route == "metabolite" and resolved_metabolites:
        return (
            "no_measurement_term_for_resolved_metabolite",
            "resolved metabolite concept(s) found but no valid measurement term passed filters",
            resolved_metabolites,
            probe_candidates,
        )

    return (
        "no_candidate_found",
        "no local candidates were retrieved from exact/token indexes",
        resolved_identity,
        probe_candidates,
    )


def review_action(reason_code: str) -> tuple[str, str]:
    mapping = {
        "identity_conflict_family_member": (
            "keep_not_mapped",
            "Top lexical matches conflict with accession identity; do not auto-map.",
        ),
        "identity_conflict_metabolite": (
            "keep_not_mapped",
            "Top lexical matches conflict with metabolite concept identity; do not auto-map.",
        ),
        "context_mismatch": (
            "check_context",
            "Candidate exists in a different matrix/tissue context; verify expected context.",
        ),
        "below_score_threshold": (
            "manual_review",
            "Candidate exists but confidence is below threshold; review top suggestions.",
        ),
        "manual_validation_required": (
            "manual_review",
            "Candidate exists but was withheld from auto-validation; confirm identity manually.",
        ),
        "needs_new_term": (
            "request_new_efo_term",
            "Closest candidate loses required granularity (for example lipoprotein size); request a new ontology term.",
        ),
        "suggest_broad_term_mapping": (
            "map_broad_term",
            "No stable external identifier (UniProt/HMDB/ChEBI/KEGG) was available; use broad-term mapping.",
        ),
        "metabolite_mode_blocked_protein_like_query": (
            "rerun_auto_mode",
            "Protein-like analyte was blocked in metabolite-only mode; rerun with entity-type=auto.",
        ),
        "unresolved_name_strict_mode": (
            "resolve_identifier",
            "Resolve query to a stable UniProt accession or metabolite identifier first.",
        ),
        "ambiguous_alias_resolution": (
            "resolve_ambiguity",
            "Alias maps to multiple accessions; disambiguate before mapping.",
        ),
        "ambiguous_metabolite_alias_resolution": (
            "resolve_ambiguity",
            "Alias maps to multiple metabolite concepts; disambiguate before mapping.",
        ),
        "no_measurement_term_for_resolved_accession": (
            "curate_new_term",
            "Resolved accession has no suitable local measurement term; consider cache curation.",
        ),
        "no_measurement_term_for_resolved_metabolite": (
            "curate_new_term",
            "Resolved metabolite concept has no suitable local measurement term; consider cache curation.",
        ),
        "no_candidate_found": (
            "expand_resources",
            "No local candidate found; consider enriching aliases/term cache.",
        ),
    }
    return mapping.get(
        reason_code,
        ("manual_review", "Manual review required."),
    )


def describe_query_identity(
    query: str,
    index: dict[str, Any],
    resolved_accessions: list[str],
) -> dict[str, str]:
    accessions = sorted({normalize(acc).upper() for acc in resolved_accessions if normalize(acc)})
    accession_index = index.get("accession_alias_index", {})

    if not accessions:
        return {
            "query_primary_accession": "",
            "query_primary_symbol": symbol_query_key(query),
            "query_primary_label": normalize(query),
            "query_subject_hints": normalize_measurement_subject(query),
        }

    primary_accession = accessions[0]
    aliases = [normalize(a) for a in accession_index.get(primary_accession, []) if normalize(a)]

    primary_symbol = ""
    primary_label = ""
    for alias in aliases:
        if not primary_symbol and is_symbol_like_alias(alias):
            primary_symbol = alias
        alias_key = norm_key(alias)
        if canonical_accession(alias):
            continue
        if is_symbol_like_alias(alias):
            continue
        if re.match(r"^ec\s+\d", alias_key):
            continue
        if not re.search(r"[a-z]", alias_key):
            continue
        if not primary_label:
            primary_label = alias

    profiles = compile_accession_profiles(accession_profiles_for_accessions(accessions, index))
    subject_terms: set[str] = set()
    for profile in profiles:
        for term in profile.subject_terms:
            term_clean = normalize(term)
            if not term_clean:
                continue
            if canonical_accession(term_clean):
                continue
            if is_symbol_like_alias(term_clean.upper()):
                continue
            subject_terms.add(term_clean)

    ordered_hints = sorted(subject_terms, key=lambda t: (-len(t), t))
    query_subject_hints = "|".join(ordered_hints[:3])
    if not query_subject_hints:
        query_subject_hints = normalize_measurement_subject(primary_label or query)

    return {
        "query_primary_accession": primary_accession,
        "query_primary_symbol": primary_symbol,
        "query_primary_label": primary_label,
        "query_subject_hints": query_subject_hints,
    }


def write_review_tsv(
    rows: list[dict[str, str]],
    path: Path,
    *,
    index: dict[str, Any],
    min_score: float,
    matrix_priority: list[str],
    measurement_context: str,
    additional_contexts: list[str],
    additional_context_keywords: list[str],
    name_mode: str,
    entity_type: str,
    top_n: int,
    reason_code_filter: set[str] | None = None,
) -> int:
    pending_rows = [row for row in rows if row.get("validation") != "validated"]
    allowed_reason_codes = {normalize(code) for code in (reason_code_filter or set()) if normalize(code)}
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "input_query",
        "input_type",
        "query_primary_accession",
        "query_primary_symbol",
        "query_primary_label",
        "query_subject_hints",
        "reason_code",
        "reason",
        "recommended_action",
        "guidance",
        "resolved_accessions",
        "candidate_count",
        "suggestion_rank",
        "suggested_efo_id",
        "suggested_label",
        "suggested_subject",
        "suggested_score",
        "suggested_matrix",
        "candidate_roundtrip_status",
        "candidate_resolved_accessions",
        "candidate_resolved_symbols",
        "source_validation",
        "source_evidence",
    ]
    count = 0
    reason_cache: dict[tuple[str, str], tuple[str, str, list[str], list[Candidate]]] = {}
    identity_cache: dict[tuple[str, tuple[str, ...]], dict[str, str]] = {}
    roundtrip_cache: dict[tuple[tuple[str, ...], str], tuple[str, tuple[str, ...], tuple[str, ...]]] = {}
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in pending_rows:
            query = normalize(row.get("input_query") or "")
            if not query:
                continue
            input_type = normalize_input_type(row.get("input_type") or "")
            reason_key = (query, input_type)
            cached_reason = reason_cache.get(reason_key)
            if cached_reason is None:
                cached_reason = infer_review_reason(
                    query,
                    index=index,
                    min_score=min_score,
                    matrix_priority=matrix_priority,
                    measurement_context=measurement_context,
                    additional_contexts=additional_contexts,
                    additional_context_keywords=additional_context_keywords,
                    name_mode=name_mode,
                    input_type=input_type,
                    entity_type=entity_type,
                )
                reason_cache[reason_key] = cached_reason
            reason_code, reason, resolved_accessions, suggestions = cached_reason
            if allowed_reason_codes and normalize(reason_code) not in allowed_reason_codes:
                continue
            action, guidance = review_action(reason_code)
            identity_key = (query, tuple(resolved_accessions))
            identity = identity_cache.get(identity_key)
            if identity is None:
                identity = describe_query_identity(query, index, resolved_accessions)
                identity_cache[identity_key] = identity
            query_symbols = query_protein_roundtrip_symbols(query, resolved_accessions, index)
            top = suggestions[: max(1, top_n)]
            if not top:
                writer.writerow(
                    {
                        "input_query": query,
                        "input_type": input_type,
                        "query_primary_accession": identity["query_primary_accession"],
                        "query_primary_symbol": identity["query_primary_symbol"],
                        "query_primary_label": identity["query_primary_label"],
                        "query_subject_hints": identity["query_subject_hints"],
                        "reason_code": reason_code,
                        "reason": reason,
                        "recommended_action": action,
                        "guidance": guidance,
                        "resolved_accessions": "|".join(resolved_accessions),
                        "candidate_count": "0",
                        "suggestion_rank": "",
                        "suggested_efo_id": "",
                        "suggested_label": "",
                        "suggested_subject": "",
                        "suggested_score": "",
                        "suggested_matrix": "",
                        "candidate_roundtrip_status": "",
                        "candidate_resolved_accessions": "",
                        "candidate_resolved_symbols": "",
                        "source_validation": row.get("validation", ""),
                        "source_evidence": row.get("evidence", ""),
                    }
                )
                count += 1
                continue

            for rank, cand in enumerate(top, start=1):
                roundtrip_status = cand.roundtrip_status
                roundtrip_accessions = cand.roundtrip_accessions
                roundtrip_symbols = cand.roundtrip_symbols
                if resolved_accessions or query_symbols:
                    roundtrip_key = (tuple(sorted([*resolved_accessions, *query_symbols])), cand.efo_id)
                    cached_roundtrip = roundtrip_cache.get(roundtrip_key)
                    if cached_roundtrip is None:
                        candidate_accessions = resolve_term_protein_accessions(cand.efo_id, index)
                        candidate_symbols = resolve_term_protein_symbols(cand.efo_id, index)
                        cached_roundtrip = (
                            protein_roundtrip_status(
                                resolved_accessions,
                                query_symbols,
                                candidate_accessions,
                                candidate_symbols,
                            ),
                            candidate_accessions,
                            candidate_symbols,
                        )
                        roundtrip_cache[roundtrip_key] = cached_roundtrip
                    roundtrip_status, roundtrip_accessions, roundtrip_symbols = cached_roundtrip
                writer.writerow(
                    {
                        "input_query": query,
                        "input_type": input_type,
                        "query_primary_accession": identity["query_primary_accession"],
                        "query_primary_symbol": identity["query_primary_symbol"],
                        "query_primary_label": identity["query_primary_label"],
                        "query_subject_hints": identity["query_subject_hints"],
                        "reason_code": reason_code,
                        "reason": reason,
                        "recommended_action": action,
                        "guidance": guidance,
                        "resolved_accessions": "|".join(resolved_accessions),
                        "candidate_count": str(len(top)),
                        "suggestion_rank": str(rank),
                        "suggested_efo_id": format_ontology_id_for_output(cand.efo_id),
                        "suggested_label": cand.label,
                        "suggested_subject": normalize_measurement_subject(cand.label),
                        "suggested_score": f"{cand.score:.3f}",
                        "suggested_matrix": measurement_matrix_bucket(cand.label) or "none",
                        "candidate_roundtrip_status": roundtrip_status,
                        "candidate_resolved_accessions": "|".join(roundtrip_accessions),
                        "candidate_resolved_symbols": "|".join(roundtrip_symbols),
                        "source_validation": row.get("validation", ""),
                        "source_evidence": row.get("evidence", ""),
                    }
                )
                count += 1
    return count


def accession_symbol_aliases(accession: str, index: dict[str, Any]) -> set[str]:
    acc = normalize(accession).upper()
    if not acc:
        return set()
    accession_index = index.get("accession_alias_index", {})
    aliases = accession_index.get(acc, [])
    symbols: set[str] = set()
    for alias in aliases:
        alias_clean = normalize(alias)
        if not alias_clean:
            continue
        if canonical_accession(alias_clean):
            continue
        if is_symbol_like_alias(alias_clean):
            symbols.add(alias_clean.upper())
    return symbols


def write_withheld_triage_tsv(
    review_path: Path,
    triage_path: Path,
    *,
    index: dict[str, Any],
) -> tuple[int, dict[str, int]]:
    if not review_path.exists():
        raise FileNotFoundError(
            f"Review file not found: {review_path}. "
            "Run map with --review-output first."
        )

    input_rows: list[dict[str, str]] = []
    with review_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if normalize(row.get("suggestion_rank") or "") != "1":
                continue
            reason_code = normalize(row.get("reason_code") or "")
            if reason_code in {"needs_new_term", "suggest_broad_term_mapping"}:
                continue
            input_rows.append(row)

    fields = [
        "input_query",
        "query_primary_label",
        "query_subject_hints",
        "query_primary_accession",
        "query_primary_symbol",
        "suggested_efo_id",
        "suggested_label",
        "suggested_subject",
        "suggested_score",
        "accession_qc_status",
        "candidate_resolved_accessions",
        "candidate_resolved_symbols",
        "same_symbol",
        "triage_decision",
        "triage_rationale",
    ]

    decision_order = {
        "reject_high": 0,
        "review_needed": 1,
        "accept_medium": 2,
        "accept_high": 3,
    }
    decision_counts: dict[str, int] = {}
    triage_rows: list[dict[str, str]] = []
    term_meta: dict[str, dict[str, Any]] = index.get("term_meta", {})
    metabolite_concept_index: dict[str, dict[str, Any]] = index.get("metabolite_concept_index", {})

    def query_metabolite_identity(
        query_value: str,
        query_label: str,
        query_hints: str,
    ) -> tuple[set[str], set[str], set[str]]:
        query_met_ids: set[str] = set()
        for source in [query_value, query_label, query_hints]:
            source_clean = normalize(source)
            if not source_clean:
                continue
            for variant in query_variants(source_clean):
                met_id = canonical_metabolite_id(variant)
                if met_id:
                    query_met_ids.add(met_id)
            for met_id in extract_metabolite_ids_from_text(source_clean):
                query_met_ids.add(met_id)

        query_concepts: set[str] = set()
        for source in [query_value, query_label]:
            source_clean = normalize(source)
            if not source_clean:
                continue
            concepts, _used_alias_resolution, _ambiguous_alias_resolution = resolve_query_metabolite_concepts(
                source_clean,
                index,
            )
            for concept in concepts:
                concept_clean = normalize(concept)
                if concept_clean:
                    query_concepts.add(concept_clean)

        concept_ids: set[str] = set()
        for concept in query_concepts:
            concept_meta = metabolite_concept_index.get(concept) or {}
            if not isinstance(concept_meta, dict):
                continue
            for met_id in normalize_metabolite_ids(concept_meta.get("ids") or []):
                concept_ids.add(met_id)

        return query_met_ids, query_concepts, concept_ids

    for row in input_rows:
        input_type = normalize_input_type(row.get("input_type") or "")
        query_text = normalize(row.get("input_query") or "")
        query_acc = normalize(row.get("query_primary_accession") or "").upper()
        query_symbol = normalize(row.get("query_primary_symbol") or "").upper()
        query_label = normalize(row.get("query_primary_label") or query_text)
        query_hints = normalize(row.get("query_subject_hints") or "")
        candidate_subject = normalize(row.get("suggested_subject") or row.get("suggested_label") or "")

        resolved_accessions: list[str] = []
        resolved_symbols: set[str] = set()
        resolved_metabolite_ids: list[str] = []
        resolved_metabolite_concepts: list[str] = []

        is_metabolite_row = input_type in {"metabolite_id", "metabolite_name"} or any(
            canonical_metabolite_id(v) for v in query_variants(query_text)
        )

        if is_metabolite_row:
            query_met_ids, query_concepts, concept_ids = query_metabolite_identity(
                query_text,
                query_label,
                query_hints,
            )
            candidate_efo_id = normalize_id(row.get("suggested_efo_id") or "")
            candidate_meta = term_meta.get(candidate_efo_id, {}) if candidate_efo_id else {}
            candidate_met_ids = set(normalize_metabolite_ids(candidate_meta.get("metabolite_ids") or []))
            if not candidate_met_ids:
                candidate_met_ids = set(
                    normalize_metabolite_ids(
                        extract_metabolite_ids_from_text(row.get("suggested_label") or "")
                        + extract_metabolite_ids_from_text(candidate_subject)
                    )
                )

            id_overlap = query_met_ids & candidate_met_ids
            concept_overlap = concept_ids & candidate_met_ids

            resolved_metabolite_ids = sorted(candidate_met_ids)
            resolved_metabolite_concepts = sorted(query_concepts)
            same_symbol = bool(id_overlap or concept_overlap)

            if id_overlap:
                qc_status = "supports_query_metabolite_id"
                triage_decision = "accept_high"
                triage_rationale = "candidate term metabolite xref matches query metabolite ID"
            elif concept_overlap:
                qc_status = "supports_query_metabolite_concept"
                triage_decision = "accept_high"
                triage_rationale = "candidate term metabolite xref matches resolved query metabolite concept ID"
            elif (query_met_ids or query_concepts or concept_ids) and candidate_met_ids:
                qc_status = "conflicts_with_query_metabolite"
                triage_decision = "reject_high"
                triage_rationale = "candidate term metabolite xrefs conflict with resolved query metabolite identity"
            else:
                qc_status = "unresolved"
                triage_decision = "review_needed"
                triage_rationale = "metabolite identity could not be resolved to stable IDs in local alias/xref index"
        else:
            query_accessions = sorted(
                {
                    normalize(accession).upper()
                    for accession in parse_pipe_list(row.get("resolved_accessions") or "")
                    if canonical_accession(accession)
                }
            )
            if query_acc and query_acc not in query_accessions:
                query_accessions.append(query_acc)
                query_accessions.sort()
            query_symbols = set(query_protein_roundtrip_symbols(query_text, query_accessions, index))
            if query_symbol:
                query_symbols.add(query_symbol)

            candidate_efo_id = normalize_id(row.get("suggested_efo_id") or "")
            resolved_accessions = list(
                resolve_term_protein_accessions(
                    candidate_efo_id,
                    index,
                    label=row.get("suggested_label") or candidate_subject,
                )
            )
            resolved_symbols = set(
                resolve_term_protein_symbols(
                    candidate_efo_id,
                    index,
                    label=row.get("suggested_label") or candidate_subject,
                )
            )
            roundtrip_status = protein_roundtrip_status(
                query_accessions,
                sorted(query_symbols),
                resolved_accessions,
                sorted(resolved_symbols),
            )
            same_symbol = bool(query_symbols & resolved_symbols)

            if roundtrip_status == "accession_match":
                qc_status = "supports_query_accession"
                triage_decision = "accept_high"
                triage_rationale = "candidate term back-maps to the same UniProt accession as the query"
            elif roundtrip_status == "symbol_match":
                qc_status = "supports_query_symbol"
                triage_decision = "accept_medium"
                triage_rationale = "candidate term back-maps to the same gene symbol as the query"
            elif roundtrip_status == "mismatch":
                qc_status = "conflicts_with_query_identity"
                triage_decision = "reject_high"
                triage_rationale = "candidate term back-maps to a different accession/gene symbol than the query"
            else:
                triage_decision = "review_needed"
                qc_status = "unresolved"
                triage_rationale = "candidate term did not back-map cleanly to an accession or gene symbol"

        decision_counts[triage_decision] = decision_counts.get(triage_decision, 0) + 1
        triage_rows.append(
            {
                "input_query": row.get("input_query", ""),
                "query_primary_label": row.get("query_primary_label", ""),
                "query_subject_hints": row.get("query_subject_hints", ""),
                "query_primary_accession": query_acc,
                "query_primary_symbol": query_symbol,
                "suggested_efo_id": row.get("suggested_efo_id", ""),
                "suggested_label": row.get("suggested_label", ""),
                "suggested_subject": row.get("suggested_subject", ""),
                "suggested_score": row.get("suggested_score", ""),
                "accession_qc_status": qc_status,
                "candidate_resolved_accessions": "|".join(
                    resolved_accessions if not is_metabolite_row else resolved_metabolite_ids
                ),
                "candidate_resolved_symbols": "|".join(
                    sorted(resolved_symbols) if not is_metabolite_row else resolved_metabolite_concepts
                ),
                "same_symbol": "yes" if same_symbol else "no",
                "triage_decision": triage_decision,
                "triage_rationale": triage_rationale,
            }
        )

    triage_rows.sort(
        key=lambda row: (
            decision_order.get(row.get("triage_decision", ""), 9),
            row.get("input_query", ""),
            row.get("suggested_efo_id", ""),
        )
    )

    triage_path.parent.mkdir(parents=True, exist_ok=True)
    with triage_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(triage_rows)

    return len(triage_rows), decision_counts


def write_tsv(rows: list[dict[str, str]], path: Path) -> None:
    fields = [
        "input_query",
        "input_type",
        "mapped_efo_id",
        "mapped_label",
        "confidence",
        "matched_via",
        "validation",
        "evidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_unmapped_tsv(rows: list[dict[str, str]], path: Path) -> int:
    unmapped = [row for row in rows if row.get("validation") == "not_mapped"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["input_query", "input_type", "validation", "evidence"],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in unmapped:
            writer.writerow(
                {
                    "input_query": row.get("input_query", ""),
                    "input_type": row.get("input_type", ""),
                    "validation": row.get("validation", ""),
                    "evidence": row.get("evidence", ""),
                }
            )
    return len(unmapped)


def write_back_cache(rows: list[dict[str, str]], path: Path, min_confidence: float) -> int:
    existing: dict[tuple[str, str], dict[str, str]] = {}
    if path.exists():
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                analyte = norm_key(row.get("analyte_key") or "")
                efo_id = normalize_id(row.get("efo_id") or "")
                if analyte and efo_id:
                    existing[(analyte, efo_id)] = row

    fields = ["analyte_key", "efo_id", "label", "confidence", "validation", "source", "evidence"]
    new_rows = 0
    for row in rows:
        if row["validation"] != "validated":
            continue
        matched_via = normalize(row.get("matched_via") or "")
        evidence = normalize(row.get("evidence") or "")
        # Only persist cache rows that are stable enough to reuse later.
        if not analyte_cache_entry_is_trusted(source=matched_via, evidence=evidence):
            continue
        try:
            confidence = float(row["confidence"])
        except ValueError:
            continue
        if confidence < min_confidence:
            continue

        analyte = norm_key(row["input_query"])
        efo_id = normalize_id(row["mapped_efo_id"])
        if not analyte or not efo_id:
            continue

        key = (analyte, efo_id)
        if key in existing:
            continue

        existing[key] = {
            "analyte_key": analyte,
            "efo_id": efo_id,
            "label": normalize(row["mapped_label"]),
            "confidence": f"{min(1.0, max(0.0, confidence)):.3f}",
            "validation": "validated",
            "source": matched_via or "auto-writeback",
            "evidence": evidence,
        }
        new_rows += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(existing):
            writer.writerow(existing[key])
    return new_rows


def run_efo_cache_refresh(
    *,
    efo_obo: str,
    download_url: str,
    download_to: str,
    measurement_root: str,
    term_cache: str,
    source: str,
    timeout: float,
) -> None:
    builder_script = Path(__file__).resolve().parent / "build_efo_measurement_cache.py"
    if not builder_script.exists():
        raise FileNotFoundError(f"EFO cache builder script not found: {builder_script}")

    cmd: list[str] = [
        sys.executable,
        str(builder_script),
        "--output",
        term_cache,
        "--download-url",
        download_url,
        "--download-to",
        download_to,
        "--measurement-root",
        measurement_root,
        "--source",
        source,
        "--timeout",
        f"{timeout}",
    ]
    if normalize(efo_obo):
        cmd.extend(["--efo-obo", efo_obo])

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603
    if proc.stdout:
        print(proc.stdout.strip())
    if proc.returncode != 0:
        stderr_msg = proc.stderr.strip() if proc.stderr else "unknown error"
        raise RuntimeError(f"EFO cache refresh failed: {stderr_msg}")


def seed_default_trait_cache_if_needed(trait_cache_path: Path) -> tuple[Path, bool]:
    """Populate default bundled trait cache from legacy location when available."""
    if trait_cache_path.exists():
        return trait_cache_path, False
    if trait_cache_path == DEFAULT_TRAIT_CACHE and LEGACY_TRAIT_CACHE.exists():
        trait_cache_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(LEGACY_TRAIT_CACHE, trait_cache_path)
        return trait_cache_path, True
    return trait_cache_path, False


def write_setup_cache_manifest(
    *,
    output_path: Path,
    index_path: Path,
    index: dict[str, Any],
    file_entries: list[tuple[str, Path, bool]],
    generated_paths: list[Path],
) -> None:
    def file_entry(name: str, path: Path, indexed: bool) -> dict[str, Any]:
        exists = path.exists()
        format_label, schema_preview = describe_cache_file_contents(path) if exists else ("", [])
        return {
            "name": name,
            "path": str(path),
            "exists": exists,
            "size_bytes": path.stat().st_size if exists else 0,
            "indexed": indexed,
            "role": cache_role(name, indexed),
            "description": cache_description(name),
            "content_format": format_label,
            "schema_preview": schema_preview,
        }

    manifest = {
        "generated_at_utc": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "index_path": str(index_path),
        "index_exists": index_path.exists(),
        "index_size_bytes": index_path.stat().st_size if index_path.exists() else 0,
        "index_summary": {
            "term_meta": len(index.get("term_meta", {})),
            "exact_term_index": len(index.get("exact_term_index", {})),
            "token_index": len(index.get("token_index", {})),
            "accession_alias_index": len(index.get("accession_alias_index", {})),
            "symbol_accession_index": len(index.get("symbol_accession_index", {})),
            "alias_accession_index": len(index.get("alias_accession_index", {})),
            "metabolite_alias_concept_index": len(index.get("metabolite_alias_concept_index", {})),
            "metabolite_concept_index": len(index.get("metabolite_concept_index", {})),
        },
        "files": [file_entry(name, path, indexed) for name, path, indexed in file_entries],
        "generated_files": [str(path) for path in generated_paths if path.exists()],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


@dataclass(frozen=True)
class CacheStatusSpec:
    name: str
    path: Path
    indexed: bool
    requirement: str


def cache_role(name: str, indexed: bool) -> str:
    role = CACHE_ROLE_BY_NAME.get(name, "")
    if role:
        return role
    return "runtime_source" if indexed else "supporting_file"


def cache_description(name: str) -> str:
    return CACHE_DESCRIPTION_BY_NAME.get(name, "")


def detect_delimited_header(path: Path) -> tuple[str, list[str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            first_line = handle.readline().strip()
    except Exception:
        return "", []
    if not first_line:
        return "", []
    delimiter = "\t" if "\t" in first_line else "," if "," in first_line else ""
    if not delimiter:
        return "", []
    header = [normalize(part) for part in first_line.split(delimiter)]
    header = [part for part in header if part]
    if not header:
        return "", []
    return delimiter, header


def describe_cache_file_contents(path: Path) -> tuple[str, list[str]]:
    if not path.exists():
        return "", []
    if path.is_dir():
        try:
            children = sorted(child.name for child in path.iterdir())[:5]
        except Exception:
            children = []
        preview = f"entries={len(children)}" if children else "directory"
        return preview, children

    suffix = path.suffix.lower()
    if suffix == ".json":
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return "json", []
        if isinstance(payload, dict):
            keys = [normalize(str(key)) for key in list(payload.keys())[:8] if normalize(str(key))]
            return "json_keys", keys
        if isinstance(payload, list):
            return "json_list", []
        return "json", []

    delimiter, header = detect_delimited_header(path)
    if header:
        format_label = "tsv_columns" if delimiter == "\t" else "csv_columns"
        return format_label, header[:8]

    return f"{suffix.lstrip('.')}_file" if suffix else "file", []


def summarize_runtime_index(index_path: Path) -> dict[str, int]:
    if not index_path.exists() or index_path.is_dir():
        return {}
    try:
        payload = json.loads(index_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    if not isinstance(payload, dict):
        return {}
    return {
        "term_meta": len(payload.get("term_meta", {})),
        "exact_term_index": len(payload.get("exact_term_index", {})),
        "token_index": len(payload.get("token_index", {})),
        "accession_alias_index": len(payload.get("accession_alias_index", {})),
        "symbol_accession_index": len(payload.get("symbol_accession_index", {})),
        "alias_accession_index": len(payload.get("alias_accession_index", {})),
        "metabolite_alias_concept_index": len(payload.get("metabolite_alias_concept_index", {})),
        "metabolite_concept_index": len(payload.get("metabolite_concept_index", {})),
    }


def setup_ukb_paths(ukb_field_catalog_path: Path | None) -> tuple[Path, Path, Path, Path]:
    catalog_path = ukb_field_catalog_path or DEFAULT_UKB_FIELD_CATALOG
    roots: list[Path] = [catalog_path.parent]
    default_root = DEFAULT_UKB_FIELD_CATALOG.parent
    if default_root not in roots:
        roots.append(default_root)

    def pick(filename: str, fallback: Path) -> Path:
        for root in roots:
            candidate = root / filename
            if candidate.exists():
                return candidate
        return roots[0] / filename if roots else fallback

    return (
        catalog_path,
        pick("field.txt", DEFAULT_UKB_FIELD_METADATA),
        pick("category.txt", DEFAULT_UKB_CATEGORY_CATALOG),
        pick("catbrowse.txt", DEFAULT_UKB_CATEGORY_TREE),
    )


def collect_default_cache_status_specs(
    *,
    term_cache_path: Path,
    analyte_cache_path: Path,
    uniprot_alias_source_path: Path,
    uniprot_alias_index_path: Path,
    metabolite_alias_path: Path,
    trait_cache_path: Path,
    ukb_field_catalog_path: Path,
    ukb_field_metadata_path: Path,
    ukb_category_catalog_path: Path,
    ukb_category_tree_path: Path,
    ukb_field_supplement_cache_path: Path,
    icd10_supplement_cache_path: Path,
    icd10_label_cache_path: Path,
    mondo_sssom_cache_path: Path,
    efo_obo_path: Path,
    index_path: Path,
) -> list[CacheStatusSpec]:
    specs: list[CacheStatusSpec] = [
        CacheStatusSpec("index", index_path, True, CACHE_REQUIREMENT_BY_NAME.get("index", "required")),
        CacheStatusSpec("term_cache", term_cache_path, True, CACHE_REQUIREMENT_BY_NAME.get("term_cache", "required")),
        CacheStatusSpec("analyte_cache", analyte_cache_path, True, CACHE_REQUIREMENT_BY_NAME.get("analyte_cache", "required")),
        CacheStatusSpec(
            "uniprot_aliases_source",
            uniprot_alias_source_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("uniprot_aliases_source", "required"),
        ),
        CacheStatusSpec(
            "uniprot_aliases_index",
            uniprot_alias_index_path,
            True,
            CACHE_REQUIREMENT_BY_NAME.get("uniprot_aliases_index", "required"),
        ),
        CacheStatusSpec(
            "metabolite_aliases",
            metabolite_alias_path,
            True,
            CACHE_REQUIREMENT_BY_NAME.get("metabolite_aliases", "required"),
        ),
        CacheStatusSpec("hmdb_source_dir", DEFAULT_HMDB_LOCAL_DIR, False, CACHE_REQUIREMENT_BY_NAME.get("hmdb_source_dir", "optional")),
        CacheStatusSpec(
            "metabolite_download_dir",
            DEFAULT_METABOLITE_DOWNLOAD_DIR,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("metabolite_download_dir", "optional"),
        ),
        CacheStatusSpec("trait_cache", trait_cache_path, False, CACHE_REQUIREMENT_BY_NAME.get("trait_cache", "required")),
        CacheStatusSpec(
            "ukb_field_catalog",
            ukb_field_catalog_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("ukb_field_catalog", "required"),
        ),
        CacheStatusSpec(
            "ukb_field_metadata",
            ukb_field_metadata_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("ukb_field_metadata", "required"),
        ),
        CacheStatusSpec(
            "ukb_category_catalog",
            ukb_category_catalog_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("ukb_category_catalog", "required"),
        ),
        CacheStatusSpec(
            "ukb_category_tree",
            ukb_category_tree_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("ukb_category_tree", "required"),
        ),
        CacheStatusSpec(
            "ukb_field_supplement_cache",
            ukb_field_supplement_cache_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("ukb_field_supplement_cache", "required"),
        ),
        CacheStatusSpec(
            "icd10_supplement_cache",
            icd10_supplement_cache_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("icd10_supplement_cache", "required"),
        ),
        CacheStatusSpec(
            "icd10_label_cache",
            icd10_label_cache_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("icd10_label_cache", "recommended"),
        ),
        CacheStatusSpec(
            "mondo_sssom_cache",
            mondo_sssom_cache_path,
            False,
            CACHE_REQUIREMENT_BY_NAME.get("mondo_sssom_cache", "recommended"),
        ),
        CacheStatusSpec("efo_obo", efo_obo_path, False, CACHE_REQUIREMENT_BY_NAME.get("efo_obo", "required")),
    ]
    return specs


def load_cache_status_specs_from_manifest(manifest_path: Path) -> tuple[list[CacheStatusSpec], dict[str, Any]]:
    manifest_meta: dict[str, Any] = {
        "path": str(manifest_path),
        "exists": manifest_path.exists(),
        "loaded": False,
        "generated_at_utc": "",
        "error": "",
        "index_summary": {},
    }
    if not manifest_path.exists():
        return [], manifest_meta

    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - reporting should keep going with defaults.
        manifest_meta["error"] = str(exc)
        return [], manifest_meta

    manifest_meta["loaded"] = True
    manifest_meta["generated_at_utc"] = normalize(str(payload.get("generated_at_utc", "")))
    manifest_index_summary = payload.get("index_summary", {})
    manifest_meta["index_summary"] = manifest_index_summary if isinstance(manifest_index_summary, dict) else {}

    specs: list[CacheStatusSpec] = []
    seen_keys: set[tuple[str, str]] = set()

    index_path = normalize(str(payload.get("index_path", "")))
    if index_path:
        fresh_index_summary = summarize_runtime_index(Path(index_path))
        if fresh_index_summary:
            manifest_meta["index_summary"] = fresh_index_summary
        index_spec = CacheStatusSpec(
            name="index",
            path=Path(index_path),
            indexed=True,
            requirement=CACHE_REQUIREMENT_BY_NAME.get("index", "required"),
        )
        specs.append(index_spec)
        seen_keys.add((index_spec.name, str(index_spec.path)))

    for entry in payload.get("files", []):
        if not isinstance(entry, dict):
            continue
        name = normalize(str(entry.get("name", "")))
        path_str = normalize(str(entry.get("path", "")))
        if not name or not path_str:
            continue
        key = (name, path_str)
        if key in seen_keys:
            continue
        seen_keys.add(key)
        specs.append(
            CacheStatusSpec(
                name=name,
                path=Path(path_str),
                indexed=bool(entry.get("indexed", False)),
                requirement=CACHE_REQUIREMENT_BY_NAME.get(name, "optional"),
            )
        )
    return specs, manifest_meta


def evaluate_cache_status_specs(
    specs: list[CacheStatusSpec],
    *,
    manifest_meta: dict[str, Any],
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    summary: Counter[str] = Counter()

    for spec in specs:
        exists = spec.path.exists()
        size_bytes = spec.path.stat().st_size if exists else 0
        mtime_utc = ""
        content_format = ""
        schema_preview: list[str] = []
        if exists:
            mtime_utc = datetime.fromtimestamp(spec.path.stat().st_mtime, UTC).isoformat().replace("+00:00", "Z")
            content_format, schema_preview = describe_cache_file_contents(spec.path)
        path_type = "dir" if exists and spec.path.is_dir() else "file"
        status = "present" if exists else "missing"
        summary[f"{spec.requirement}_{status}"] += 1
        rows.append(
            {
                "name": spec.name,
                "path": str(spec.path),
                "exists": exists,
                "indexed": spec.indexed,
                "role": cache_role(spec.name, spec.indexed),
                "description": cache_description(spec.name),
                "requirement": spec.requirement,
                "status": status,
                "path_type": path_type,
                "size_bytes": size_bytes,
                "mtime_utc": mtime_utc,
                "content_format": content_format,
                "schema_preview": schema_preview,
            }
        )

    required_total = summary["required_present"] + summary["required_missing"]
    recommended_total = summary["recommended_present"] + summary["recommended_missing"]
    optional_total = summary["optional_present"] + summary["optional_missing"]
    status_payload = {
        "generated_at_utc": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "manifest": manifest_meta,
        "summary": {
            "files_total": len(rows),
            "required_total": required_total,
            "required_present": summary["required_present"],
            "required_missing": summary["required_missing"],
            "recommended_total": recommended_total,
            "recommended_present": summary["recommended_present"],
            "recommended_missing": summary["recommended_missing"],
            "optional_total": optional_total,
            "optional_present": summary["optional_present"],
            "optional_missing": summary["optional_missing"],
        },
        "files": rows,
    }
    return status_payload


def print_cache_status_report(report: dict[str, Any]) -> None:
    manifest = report.get("manifest", {})
    manifest_path = manifest.get("path", "")
    if manifest.get("loaded"):
        generated_at = normalize(str(manifest.get("generated_at_utc", "")))
        generated_suffix = f", generated_at={generated_at}" if generated_at else ""
        print(f"[OK] cache manifest loaded: {manifest_path}{generated_suffix}")
    elif manifest.get("exists"):
        manifest_error = normalize(str(manifest.get("error", "")))
        if manifest_error:
            print(f"[WARN] cache manifest could not be parsed: {manifest_path}; reason={manifest_error}")
        else:
            print(f"[OK] cache manifest present but not used: {manifest_path}")
    else:
        print(f"[WARN] cache manifest not found: {manifest_path} (falling back to default cache paths)")

    summary = report.get("summary", {})
    print(
        "[OK] cache status summary "
        f"(required_present={summary.get('required_present', 0)}/{summary.get('required_total', 0)}, "
        f"required_missing={summary.get('required_missing', 0)}, "
        f"recommended_missing={summary.get('recommended_missing', 0)}, "
        f"optional_missing={summary.get('optional_missing', 0)})"
    )
    index_summary = manifest.get("index_summary", {})
    if isinstance(index_summary, dict) and index_summary:
        summary_bits = [
            f"term_meta={int(index_summary.get('term_meta', 0))}",
            f"exact_term_index={int(index_summary.get('exact_term_index', 0))}",
            f"token_index={int(index_summary.get('token_index', 0))}",
            f"accession_alias_index={int(index_summary.get('accession_alias_index', 0))}",
            f"symbol_accession_index={int(index_summary.get('symbol_accession_index', 0))}",
            f"alias_accession_index={int(index_summary.get('alias_accession_index', 0))}",
            f"metabolite_alias_concept_index={int(index_summary.get('metabolite_alias_concept_index', 0))}",
            f"metabolite_concept_index={int(index_summary.get('metabolite_concept_index', 0))}",
        ]
        print(f"[OK] runtime index summary ({', '.join(summary_bits)})")
    for row in report.get("files", []):
        status_flag = "OK" if row.get("exists") else "MISSING"
        role = normalize(str(row.get("role", ""))).replace("_", " ")
        description = normalize(str(row.get("description", "")))
        format_label = normalize(str(row.get("content_format", "")))
        schema_preview = row.get("schema_preview", [])
        preview_suffix = ""
        if isinstance(schema_preview, list) and schema_preview:
            preview_suffix = f" preview={','.join(str(item) for item in schema_preview[:6])}"
        print(
            f"[{status_flag}] requirement={row.get('requirement', 'optional')} "
            f"name={row.get('name', '')} "
            f"role={role or 'supporting file'} "
            f"indexed={bool(row.get('indexed', False))} "
            f"size_bytes={int(row.get('size_bytes', 0))} "
            f"format={format_label or row.get('path_type', 'file')} "
            f"path={row.get('path', '')}"
        )
        if description:
            print(f"[INFO] {row.get('name', '')}: {description}{preview_suffix}")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def trait_catalog_tokenize(text: str) -> set[str]:
    normalized_text = normalize_trait_text_for_similarity(text)
    return informative_trait_tokens(normalized_text)


def is_generic_trait_label(label: str) -> bool:
    key = norm_key(label)
    if not key:
        return True
    generic = {
        "disease",
        "disorder",
        "phenotype",
        "trait",
        "measurement",
        "clinical measurement",
        "symptom",
        "abnormality",
    }
    return key in generic


def parse_catalog_mapped_pairs(
    mapped_trait: str,
    mapped_trait_uri: str,
) -> tuple[list[tuple[str, str]], list[str]]:
    raw_ids = [canonicalize_trait_ontology_id(part) for part in split_multi_ids(mapped_trait_uri)]
    raw_ids = [item for item in raw_ids if item]
    labels = split_multi_labels(mapped_trait, expected_n=len(raw_ids))
    pairs: list[tuple[str, str]] = []
    issues: list[str] = []
    for idx, mapped_id in enumerate(raw_ids):
        if trait_ontology_prefix(mapped_id) not in TRAIT_PREFERRED_PREFIX_ORDER:
            continue
        label = normalize(labels[idx] if idx < len(labels) else "")
        pairs.append((mapped_id, label))
    if not pairs and normalize(mapped_trait_uri):
        issues.append("no_supported_mapped_ids")
    return pairs, issues


def resolve_catalog_mapping_pairs(
    pairs: list[tuple[str, str]],
    *,
    ontology_terms: dict[str, TraitOntologyTerm],
    obsolete_terms: dict[str, TraitObsoleteTerm],
) -> tuple[list[tuple[str, str]], list[str], list[str], list[str]]:
    resolved_pairs: list[tuple[str, str]] = []
    unresolved_ids: list[str] = []
    obsolete_ids: list[str] = []
    resolution_modes: list[str] = []

    for mapped_id, mapped_label in pairs:
        if mapped_id in ontology_terms:
            label = mapped_label or ontology_terms[mapped_id].label
            if (mapped_id, label) not in resolved_pairs:
                resolved_pairs.append((mapped_id, label))
            resolution_modes.append("direct")
            continue

        obsolete_meta = obsolete_terms.get(mapped_id)
        if obsolete_meta is None:
            unresolved_ids.append(mapped_id)
            resolution_modes.append("missing")
            continue

        obsolete_ids.append(mapped_id)
        replacements = [rid for rid in obsolete_meta.replaced_by if rid in ontology_terms]
        if replacements:
            for rid in replacements:
                label = ontology_terms[rid].label
                if (rid, label) not in resolved_pairs:
                    resolved_pairs.append((rid, label))
            resolution_modes.append("replaced_by")
            continue

        consider = [rid for rid in obsolete_meta.consider if rid in ontology_terms]
        if len(consider) == 1:
            rid = consider[0]
            label = ontology_terms[rid].label
            if (rid, label) not in resolved_pairs:
                resolved_pairs.append((rid, label))
            resolution_modes.append("consider_single")
            continue

        unresolved_ids.append(mapped_id)
        resolution_modes.append("obsolete_unresolved")

    return resolved_pairs, unresolved_ids, obsolete_ids, resolution_modes


def build_catalog_trait_refresh(
    *,
    studies_tsv: Path,
    trait_cache_path: Path,
    ontology_index: dict[str, Any],
    publication_tag: str,
    studies_sha256: str = "",
) -> dict[str, Any]:
    if not studies_tsv.exists():
        raise FileNotFoundError(f"GWAS Catalog studies TSV not found: {studies_tsv}")
    if not trait_cache_path.exists():
        raise FileNotFoundError(f"Trait cache TSV not found: {trait_cache_path}")

    ontology_terms: dict[str, TraitOntologyTerm] = ontology_index.get("terms", {})
    obsolete_terms: dict[str, TraitObsoleteTerm] = ontology_index.get("obsolete_terms", {})

    cache_index = load_trait_cache_index(trait_cache_path)
    records: list[TraitCacheRecord] = cache_index.get("records", [])
    existing_trait_mappings: dict[str, set[str]] = {}
    existing_exact_pairs: set[tuple[str, str]] = set()
    for rec in records:
        if not rec.reported_trait:
            continue
        trait_key = normalize_trait_text_for_similarity(rec.reported_trait)
        if not trait_key:
            continue
        mapped_key = "|".join(rec.mapped_ids)
        if mapped_key:
            existing_trait_mappings.setdefault(trait_key, set()).add(mapped_key)
        existing_exact_pairs.add((trait_key, mapped_key))

    qc_rows: list[dict[str, str]] = []
    by_trait_key: dict[str, set[str]] = {}
    with studies_tsv.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            trait = normalize(row.get("DISEASE/TRAIT", ""))
            mapped_trait = normalize(row.get("MAPPED_TRAIT", ""))
            mapped_trait_uri = normalize(row.get("MAPPED_TRAIT_URI", ""))
            if not trait or not mapped_trait_uri:
                continue

            parsed_pairs, parse_issues = parse_catalog_mapped_pairs(mapped_trait, mapped_trait_uri)
            resolved_pairs, unresolved_ids, obsolete_ids, resolution_modes = resolve_catalog_mapping_pairs(
                parsed_pairs,
                ontology_terms=ontology_terms,
                obsolete_terms=obsolete_terms,
            )

            resolved_ids = [pair[0] for pair in resolved_pairs]
            resolved_labels = [pair[1] for pair in resolved_pairs]
            resolved_id_key = "|".join(resolved_ids)
            if resolved_id_key:
                by_trait_key.setdefault(normalize_trait_text_for_similarity(trait), set()).add(resolved_id_key)

            source_label_tokens = trait_catalog_tokenize(mapped_trait)
            resolved_label_tokens = set()
            for label in resolved_labels:
                resolved_label_tokens |= trait_catalog_tokenize(label)
            overlap_n = len(source_label_tokens & resolved_label_tokens)
            denom = max(1, len(source_label_tokens))
            label_match_ratio = overlap_n / denom
            if not source_label_tokens and resolved_labels:
                label_match_ratio = 1.0

            label_mode = "no_source_label"
            if mapped_trait and resolved_labels:
                if norm_key(mapped_trait) == norm_key("|".join(resolved_labels)):
                    label_mode = "full_exact"
                elif overlap_n > 0:
                    label_mode = "token_overlap"
                else:
                    label_mode = "low_overlap"

            high_confidence = bool(resolved_ids) and not unresolved_ids and label_match_ratio >= 0.50
            qc_rows.append(
                {
                    "DISEASE/TRAIT": trait,
                    "MAPPED_TRAIT": mapped_trait,
                    "MAPPED_TRAIT_URI": mapped_trait_uri,
                    "resolved_mapped_ids": "|".join(resolved_ids),
                    "resolved_mapped_labels": "|".join(resolved_labels),
                    "resolved_id_count": str(len(resolved_ids)),
                    "unresolved_ids": "|".join(unresolved_ids + parse_issues),
                    "obsolete_ids": "|".join(obsolete_ids),
                    "id_resolution_modes": "|".join(resolution_modes),
                    "label_match_ratio": f"{label_match_ratio:.3f}",
                    "label_match_modes": label_mode,
                    "confidence": "high" if high_confidence else "low",
                    "qc_status": "ok" if high_confidence else "fail",
                }
            )

    high_conf_initial = [row for row in qc_rows if row.get("confidence") == "high" and row.get("qc_status") == "ok"]
    high_conf_by_trait: list[dict[str, str]] = []
    for row in high_conf_initial:
        trait_key = normalize_trait_text_for_similarity(row.get("DISEASE/TRAIT", ""))
        if not trait_key:
            continue
        unique_mappings = by_trait_key.get(trait_key, set())
        if len(unique_mappings) == 1:
            high_conf_by_trait.append(row)

    high_conf_specific: list[dict[str, str]] = []
    for row in high_conf_by_trait:
        trait = normalize(row.get("DISEASE/TRAIT", ""))
        mapped_labels = split_multi_labels(row.get("resolved_mapped_labels", ""), expected_n=1)
        if not mapped_labels:
            continue
        first_label = mapped_labels[0]
        if is_generic_trait_label(first_label):
            continue
        trait_tokens = trait_catalog_tokenize(trait)
        label_tokens = trait_catalog_tokenize(first_label)
        if trait_tokens and label_tokens and not (trait_tokens & label_tokens):
            continue
        high_conf_specific.append(row)

    rows_to_add: list[dict[str, str]] = []
    dedup_new: set[tuple[str, str]] = set()
    conflict_count = 0
    for row in high_conf_specific:
        trait = normalize(row.get("DISEASE/TRAIT", ""))
        trait_key = normalize_trait_text_for_similarity(trait)
        mapped_ids = normalize(row.get("resolved_mapped_ids", ""))
        mapped_labels = normalize(row.get("resolved_mapped_labels", ""))
        if not trait_key or not mapped_ids:
            continue
        cache_key = (trait_key, mapped_ids)
        if cache_key in existing_exact_pairs:
            continue
        existing_for_trait = existing_trait_mappings.get(trait_key, set())
        if existing_for_trait and mapped_ids not in existing_for_trait:
            conflict_count += 1
            continue
        if cache_key in dedup_new:
            continue
        dedup_new.add(cache_key)
        rows_to_add.append(
            {
                "UKBB data field": "-",
                "ICD10": "-",
                "PheCode": "-",
                "Lookup text": trait,
                "Curated reported trait": trait,
                "main ontology URI(s)": mapped_ids,
                "main ontology label(s)": mapped_labels,
                "Publication": publication_tag,
                "User submitted trait": "-",
                "qc_status": "validated",
                "qc_resolution_modes": "catalog_import_high_conf",
                "qc_label_match_modes": row.get("label_match_modes", ""),
                "qc_unresolved_ids": "",
                "final_update_applied": "yes",
                "final_update_reason": "catalog_import_high_conf_unambiguous_specific_nonconflicting",
            }
        )

    merged_rows: list[dict[str, str]] = []
    with trait_cache_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            merged_rows.append({column: normalize(row.get(column, "")) for column in TRAIT_CACHE_COLUMNS})
    merged_rows.extend(rows_to_add)

    summary = {
        "input_rows": len(qc_rows),
        "high_conf_initial": len(high_conf_initial),
        "after_unambiguous_trait_filter": len(high_conf_by_trait),
        "after_specific_filter_candidate_rows": len(high_conf_specific),
        "cache_rows_before": len(merged_rows) - len(rows_to_add),
        "cache_rows_added": len(rows_to_add),
        "cache_rows_after": len(merged_rows),
        "conflicting_rows_skipped": conflict_count,
        "policy": (
            "high && unambiguous mapping per normalized trait && "
            "specific_label_overlap && nonconflicting with existing cache"
        ),
        "studies_tsv": str(studies_tsv),
        "studies_tsv_sha256": studies_sha256 or file_sha256(studies_tsv),
        "trait_cache": str(trait_cache_path),
        "publication_tag": publication_tag,
    }

    return {
        "qc_rows": qc_rows,
        "high_conf_rows": high_conf_specific,
        "rows_to_add": rows_to_add,
        "merged_rows": merged_rows,
        "summary": summary,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Bulk/offline mapper for measurement EFO terms")
    sub = parser.add_subparsers(dest="command")

    p_idx = sub.add_parser("index-build", help="Build local JSON index from cache/reference TSVs")
    p_idx.add_argument("--term-cache", default=str(DEFAULT_TERM_CACHE), help="Term cache TSV path")
    p_idx.add_argument("--analyte-cache", default=str(DEFAULT_ANALYTE_CACHE), help="Analyte->EFO cache TSV path")
    p_idx.add_argument("--uniprot-aliases", default=str(DEFAULT_UNIPROT_ALIASES), help="UniProt accession alias TSV")
    p_idx.add_argument(
        "--metabolite-aliases",
        default=str(DEFAULT_METABOLITE_ALIASES),
        help="Metabolite alias TSV (HMDB/ChEBI/KEGG/local synonyms)",
    )
    p_idx.add_argument("--output-index", default=str(DEFAULT_INDEX), help="Output JSON index path")

    p_refresh = sub.add_parser(
        "refresh-efo-cache",
        help="Refresh full EFO measurement term cache from OBO and optionally rebuild index",
    )
    p_refresh.add_argument(
        "--term-cache",
        default=str(DEFAULT_TERM_CACHE),
        help="Output term cache TSV path",
    )
    p_refresh.add_argument(
        "--efo-obo",
        default="",
        help="Optional local EFO OBO path. If omitted, downloads from --download-url",
    )
    p_refresh.add_argument(
        "--download-url",
        default=DEFAULT_EFO_OBO_URL,
        help=(
            "EFO OBO URL used when --efo-obo is not provided "
            f"(default: {DEFAULT_EFO_OBO_URL})"
        ),
    )
    p_refresh.add_argument(
        "--download-to",
        default=str(DEFAULT_EFO_OBO_LOCAL),
        help="Local path to save downloaded EFO OBO",
    )
    p_refresh.add_argument(
        "--measurement-root",
        default=DEFAULT_MEASUREMENT_ROOT,
        help="Measurement ontology root ID",
    )
    p_refresh.add_argument(
        "--source",
        default="efo-obo-measurement-branch",
        help="Source tag written into refreshed term cache rows",
    )
    p_refresh.add_argument(
        "--timeout",
        type=float,
        default=60.0,
        help="Timeout (seconds) for OBO download",
    )
    p_refresh.add_argument(
        "--rebuild-index",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Rebuild local JSON index after cache refresh (recommended, default: on). "
            "Use --no-rebuild-index only for staged/advanced workflows."
        ),
    )
    p_refresh.add_argument(
        "--index",
        default=str(DEFAULT_INDEX),
        help="JSON index path to build when --rebuild-index is enabled",
    )
    p_refresh.add_argument(
        "--analyte-cache",
        default=str(DEFAULT_ANALYTE_CACHE),
        help="Analyte cache TSV used if rebuilding index",
    )
    p_refresh.add_argument(
        "--uniprot-aliases",
        default=str(DEFAULT_UNIPROT_ALIASES),
        help="UniProt alias TSV used if rebuilding index",
    )
    p_refresh.add_argument(
        "--metabolite-aliases",
        default=str(DEFAULT_METABOLITE_ALIASES),
        help="Metabolite alias TSV used if rebuilding index",
    )

    p_setup = sub.add_parser(
        "setup-bundled-caches",
        help=(
            "Offline-first setup using bundled local caches "
            "(protein aliases, metabolite aliases, trait cache); auto-provisions missing "
            "efo.obo and UKB metadata, builds ICD10 side caches, then builds the local index"
        ),
    )
    p_setup.add_argument(
        "--term-cache",
        default=str(DEFAULT_TERM_CACHE),
        help="Bundled measurement term cache TSV path",
    )
    p_setup.add_argument(
        "--analyte-cache",
        default=str(DEFAULT_ANALYTE_CACHE),
        help="Bundled analyte cache TSV path",
    )
    p_setup.add_argument(
        "--uniprot-aliases",
        default=str(DEFAULT_UNIPROT_ALIASES),
        help="Bundled protein alias cache TSV path",
    )
    p_setup.add_argument(
        "--uniprot-profile",
        choices=["light", "full"],
        default="light",
        help=(
            "UniProt cache profile for setup index build. light: build reviewed+gene lightweight cache "
            "(faster, recommended). full: use bundled full alias cache directly."
        ),
    )
    p_setup.add_argument(
        "--uniprot-light-output",
        default=str(DEFAULT_UNIPROT_ALIASES_LIGHT),
        help="Output path for generated lightweight UniProt alias cache",
    )
    p_setup.add_argument(
        "--uniprot-light-reviewed-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="For lightweight profile: keep reviewed-like accessions only (default: on)",
    )
    p_setup.add_argument(
        "--uniprot-light-require-gene",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="For lightweight profile: require non-empty primary gene symbol (default: on)",
    )
    p_setup.add_argument(
        "--uniprot-light-enrich-gene-ids",
        action="store_true",
        help=(
            "For lightweight profile: backfill missing UniProt gene_ids online after cache build "
            "(optional, can take time)."
        ),
    )
    p_setup.add_argument(
        "--uniprot-timeout",
        type=float,
        default=20.0,
        help="Timeout per UniProt request used by optional enrichment (seconds)",
    )
    p_setup.add_argument(
        "--uniprot-workers",
        type=int,
        default=8,
        help="Parallel workers used by optional UniProt enrichment",
    )
    p_setup.add_argument(
        "--uniprot-source",
        default="uniprot-live-backfill",
        help="Source tag for UniProt enrichment rows",
    )
    p_setup.add_argument(
        "--metabolite-aliases",
        default=str(DEFAULT_METABOLITE_ALIASES),
        help="Bundled metabolite alias cache TSV path",
    )
    p_setup.add_argument(
        "--trait-cache",
        default=str(DEFAULT_TRAIT_CACHE),
        help="Bundled trait cache TSV path",
    )
    p_setup.add_argument(
        "--ukb-field-catalog",
        default=str(DEFAULT_UKB_FIELD_CATALOG),
        help="Bundled UKB field catalog TSV path (for example references/ukb/fieldsum.txt)",
    )
    p_setup.add_argument(
        "--ukb-field-catalog-url",
        default=DEFAULT_UKB_FIELD_CATALOG_URL,
        help="Official UKB field catalog TSV URL used when the local file is missing",
    )
    p_setup.add_argument(
        "--ukb-field-metadata-url",
        default=DEFAULT_UKB_FIELD_METADATA_URL,
        help="Official UKB field metadata TSV URL used when the local file is missing",
    )
    p_setup.add_argument(
        "--ukb-category-catalog-url",
        default=DEFAULT_UKB_CATEGORY_CATALOG_URL,
        help="Official UKB category catalog TSV URL used when the local file is missing",
    )
    p_setup.add_argument(
        "--ukb-category-tree-url",
        default=DEFAULT_UKB_CATEGORY_TREE_URL,
        help="Official UKB category tree TSV URL used when the local file is missing",
    )
    p_setup.add_argument(
        "--ukb-download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) when auto-downloading missing UKB metadata files",
    )
    p_setup.add_argument(
        "--ukb-field-supplement-cache",
        default=str(DEFAULT_UKB_FIELD_SUPPLEMENT_CACHE),
        help=(
            "Output UKB supplemental field trait cache TSV path "
            "(generated from UKB metadata + ontology during setup)"
        ),
    )
    p_setup.add_argument(
        "--build-ukb-field-supplement",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Build UKB supplemental field trait cache during setup (default: on)",
    )
    p_setup.add_argument(
        "--icd10-supplement-cache",
        default=str(DEFAULT_ICD10_SUPPLEMENT_CACHE),
        help=(
            "Output ICD10 supplemental trait cache TSV path "
            "(generated from ontology ICD10 xrefs plus MONDO SSSOM enrichment when cache is available)"
        ),
    )
    p_setup.add_argument(
        "--icd10-label-cache",
        default=str(DEFAULT_ICD10_LABEL_CACHE),
        help="Output ICD10 label cache TSV path",
    )
    p_setup.add_argument(
        "--icd10-mondo-sssom-cache",
        default=str(DEFAULT_MONDO_SSSOM_CACHE),
        help=(
            "Local MONDO SSSOM cache path used by setup for ICD10 enrichment. "
            "setup-bundled-caches refreshes this file from --icd10-mondo-sssom-url when reachable."
        ),
    )
    p_setup.add_argument(
        "--build-icd10-supplement",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Build ICD10 supplemental trait cache during setup (default: on)",
    )
    p_setup.add_argument(
        "--build-icd10-label-cache",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Build ICD10 label cache during setup when missing (default: on)",
    )
    p_setup.add_argument(
        "--icd10-label-source-url",
        default=DEFAULT_ICD10_LABEL_SOURCE_URL,
        help="Official ICD10 label source ZIP URL used when the local label cache is missing",
    )
    p_setup.add_argument(
        "--icd10-label-download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) for ICD10 label source download during setup",
    )
    p_setup.add_argument(
        "--icd10-mondo-sssom-url",
        default=DEFAULT_MONDO_SSSOM_URL,
        help=(
            "MONDO SSSOM mappings URL used by setup to refresh local MONDO ICD10 cache. "
            "When unavailable (for example offline), setup falls back to the existing local cache file."
        ),
    )
    p_setup.add_argument(
        "--icd10-mondo-download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) for MONDO SSSOM refresh during ICD10 supplemental build",
    )
    p_setup.add_argument(
        "--efo-obo",
        default=str(DEFAULT_EFO_OBO_LOCAL),
        help=(
            "Local EFO OBO path required by trait-map. If missing, setup will auto-provision "
            "from --efo-obo-bundled-url."
        ),
    )
    p_setup.add_argument(
        "--efo-obo-bundled-url",
        default=DEFAULT_EFO_OBO_BUNDLED_URL,
        help=(
            "Pinned bundled EFO OBO URL (recommended: your repo release asset, e.g. efo.obo.gz). "
            "Used only when --efo-obo is missing."
        ),
    )
    p_setup.add_argument(
        "--efo-obo-download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) when auto-fetching bundled EFO OBO",
    )
    p_setup.add_argument(
        "--seed-legacy-trait-cache",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "If bundled trait cache is missing, copy from legacy "
            "final_output/code-EFO-mappings_final_mapping_cache.tsv when available "
            "(default: on)"
        ),
    )
    p_setup.add_argument(
        "--index",
        default=str(DEFAULT_INDEX),
        help="JSON index output path",
    )
    p_setup.add_argument(
        "--rebuild-index",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Rebuild index from bundled caches (default: on)",
    )
    p_setup.add_argument(
        "--cache-manifest-output",
        default=str(DEFAULT_SETUP_CACHE_MANIFEST),
        help=(
            "Write setup cache/index inventory JSON for auditability "
            "(default: final_output/analyte_mapper_cache_manifest.json)"
        ),
    )

    p_cache_status = sub.add_parser(
        "cache-status",
        help=(
            "Audit cache/index file presence after setup. "
            "Prints present/missing status and can write a JSON status report."
        ),
    )
    p_cache_status.add_argument(
        "--manifest",
        default=str(DEFAULT_SETUP_CACHE_MANIFEST),
        help=(
            "Setup cache manifest JSON path (default: final_output/analyte_mapper_cache_manifest.json). "
            "When present, paths are read from the manifest."
        ),
    )
    p_cache_status.add_argument(
        "--prefer-manifest",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use paths listed in --manifest when available (default: on)",
    )
    p_cache_status.add_argument("--term-cache", default=str(DEFAULT_TERM_CACHE), help="Measurement term cache TSV path")
    p_cache_status.add_argument("--analyte-cache", default=str(DEFAULT_ANALYTE_CACHE), help="Analyte cache TSV path")
    p_cache_status.add_argument(
        "--uniprot-aliases",
        default=str(DEFAULT_UNIPROT_ALIASES),
        help="UniProt alias source TSV path",
    )
    p_cache_status.add_argument(
        "--uniprot-light-output",
        default=str(DEFAULT_UNIPROT_ALIASES_LIGHT),
        help="Lightweight UniProt alias TSV path",
    )
    p_cache_status.add_argument(
        "--metabolite-aliases",
        default=str(DEFAULT_METABOLITE_ALIASES),
        help="Metabolite alias TSV path",
    )
    p_cache_status.add_argument("--trait-cache", default=str(DEFAULT_TRAIT_CACHE), help="Trait cache TSV path")
    p_cache_status.add_argument(
        "--ukb-field-catalog",
        default=str(DEFAULT_UKB_FIELD_CATALOG),
        help="UKB field catalog TSV path (used to resolve UKB metadata sibling files)",
    )
    p_cache_status.add_argument(
        "--ukb-field-supplement-cache",
        default=str(DEFAULT_UKB_FIELD_SUPPLEMENT_CACHE),
        help="UKB supplemental trait cache TSV path",
    )
    p_cache_status.add_argument(
        "--icd10-supplement-cache",
        default=str(DEFAULT_ICD10_SUPPLEMENT_CACHE),
        help="ICD10 supplemental trait cache TSV path",
    )
    p_cache_status.add_argument(
        "--icd10-label-cache",
        default=str(DEFAULT_ICD10_LABEL_CACHE),
        help="ICD10 label cache TSV path",
    )
    p_cache_status.add_argument(
        "--icd10-mondo-sssom-cache",
        default=str(DEFAULT_MONDO_SSSOM_CACHE),
        help="MONDO SSSOM ICD10 cache TSV path",
    )
    p_cache_status.add_argument("--efo-obo", default=str(DEFAULT_EFO_OBO_LOCAL), help="EFO OBO path")
    p_cache_status.add_argument("--index", default=str(DEFAULT_INDEX), help="JSON index path")
    p_cache_status.add_argument(
        "--output-json",
        default="",
        help=(
            "Optional JSON status output path "
            f"(for example {DEFAULT_SETUP_CACHE_STATUS})"
        ),
    )
    p_cache_status.add_argument(
        "--strict",
        action="store_true",
        help="Exit with code 1 when required cache files are missing",
    )

    p_icd10_labels = sub.add_parser(
        "icd10-label-cache-build",
        help="Build ICD10 label cache from the official CMS ICD-10-CM tabular-order release",
    )
    p_icd10_labels.add_argument(
        "--output",
        default=str(DEFAULT_ICD10_LABEL_CACHE),
        help="Output ICD10 label cache TSV path",
    )
    p_icd10_labels.add_argument(
        "--source-url",
        default=DEFAULT_ICD10_LABEL_SOURCE_URL,
        help="Official ICD10 source ZIP URL (default: CMS 2026 tabular-order release)",
    )
    p_icd10_labels.add_argument(
        "--download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) for ICD10 label source download",
    )

    p_uniprot = sub.add_parser(
        "uniprot-alias-enrich",
        help="Backfill gene IDs (and aliases) in UniProt alias cache via UniProt REST",
    )
    p_uniprot.add_argument(
        "--uniprot-aliases",
        default=str(DEFAULT_UNIPROT_ALIASES),
        help="UniProt alias TSV path to enrich in place",
    )
    p_uniprot.add_argument(
        "--timeout",
        type=float,
        default=20.0,
        help="Timeout per UniProt request (seconds)",
    )
    p_uniprot.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Parallel workers for UniProt API fetches",
    )
    p_uniprot.add_argument(
        "--source",
        default="uniprot-live-backfill",
        help="Source tag written into updated UniProt alias rows",
    )
    p_uniprot.add_argument(
        "--refresh-all",
        action="store_true",
        help="Refresh all accessions, not only rows missing gene IDs",
    )

    p_uniprot_light = sub.add_parser(
        "uniprot-alias-build-light",
        help="Build lightweight UniProt alias cache (reviewed+gene focused) from source TSV",
    )
    p_uniprot_light.add_argument(
        "--source",
        default=str(DEFAULT_UNIPROT_ALIASES),
        help="Source UniProt alias TSV path",
    )
    p_uniprot_light.add_argument(
        "--output",
        default=str(DEFAULT_UNIPROT_ALIASES_LIGHT),
        help="Output lightweight UniProt alias TSV path",
    )
    p_uniprot_light.add_argument(
        "--reviewed-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Keep reviewed-like accessions only (default: on)",
    )
    p_uniprot_light.add_argument(
        "--require-gene",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Require non-empty primary gene symbol (default: on)",
    )

    p_met_alias = sub.add_parser(
        "metabolite-alias-build",
        help=(
            "Build/merge local metabolite alias TSV from a pinned local HMDB XML/SDF "
            "(plus optional ChEBI/KEGG resources)"
        ),
    )
    p_met_alias.add_argument(
        "--output",
        default=str(DEFAULT_METABOLITE_ALIASES),
        help="Output metabolite alias TSV path",
    )
    p_met_alias.add_argument(
        "--hmdb-xml",
        default="",
        help=(
            "HMDB metabolites XML(.gz) or SDF(.gz) path. If omitted, the command auto-detects "
            "local files in references/hmdb_source and references/metabolite_downloads."
        ),
    )
    p_met_alias.add_argument(
        "--chebi-names",
        default="",
        help="Optional ChEBI names TSV(.gz) path",
    )
    p_met_alias.add_argument(
        "--kegg-conv",
        default="",
        help="Optional KEGG conv mapping file path",
    )
    p_met_alias.add_argument(
        "--merge-existing",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Merge parsed aliases into existing output TSV (default: on)",
    )
    p_met_alias.add_argument(
        "--download-common-sources",
        action="store_true",
        help=(
            "Download optional ChEBI names and KEGG mappings when local paths are not provided."
        ),
    )
    p_met_alias.add_argument(
        "--download-hmdb",
        action="store_true",
        help=(
            "Force-download HMDB bundle URL. By default, if no local HMDB source is found "
            "and no other sources are provided, HMDB download is attempted automatically."
        ),
    )
    p_met_alias.add_argument(
        "--download-dir",
        default=str(DEFAULT_METABOLITE_DOWNLOAD_DIR),
        help="Directory for downloaded metabolite source files",
    )
    p_met_alias.add_argument(
        "--download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) for metabolite source downloads",
    )
    p_met_alias.add_argument(
        "--hmdb-url",
        default=DEFAULT_HMDB_XML_URL,
        help=(
            "Pinned HMDB bundle URL (for example a GitHub release asset). "
            "Used for automatic fallback download when no local HMDB source is found."
        ),
    )
    p_met_alias.add_argument(
        "--chebi-names-url",
        default=DEFAULT_CHEBI_NAMES_URL,
        help="ChEBI names TSV(.gz) URL used with --download-common-sources",
    )
    p_met_alias.add_argument(
        "--kegg-conv-url",
        default=DEFAULT_KEGG_CONV_URL,
        help="KEGG conv URL used with --download-common-sources",
    )

    p_map = sub.add_parser("map", help="Map queries using local JSON index")
    p_map.add_argument(
        "--input",
        required=True,
        help=(
            "Input txt/csv/tsv with analyte queries. Optional tabular input_type column "
            "(input_type/query_type/id_type/type): accession, gene_symbol, gene_id, protein_name, "
            "metabolite_id, metabolite_name, auto"
        ),
    )
    p_map.add_argument("--output", required=True, help="Output mapping TSV")
    p_map.add_argument("--index", default=str(DEFAULT_INDEX), help="JSON index path")
    p_map.add_argument("--top-k", type=int, default=1, help="Candidates to keep per input")
    p_map.add_argument("--min-score", type=float, default=0.55, help="Minimum score to emit")
    p_map.add_argument("--workers", type=int, default=1, help="Parallel workers for map stage")
    p_map.add_argument(
        "--parallel-mode",
        choices=["thread", "process", "auto"],
        default="auto",
        help=(
            "Parallel backend when workers>1: thread, process, or auto "
            "(auto uses process mode for larger batches)"
        ),
    )
    p_map.add_argument(
        "--progress",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show mapping progress bar/logs (default: on)",
    )
    p_map.add_argument(
        "--matrix-priority",
        default=",".join(DEFAULT_MATRIX_PRIORITY),
        help="Preferred sample matrices in order, comma-separated (default: plasma,blood,serum)",
    )
    p_map.add_argument(
        "--measurement-context",
        default=DEFAULT_MEASUREMENT_CONTEXT,
        help=(
            "Primary measurement context: blood, plasma, serum, cerebrospinal_fluid, urine, "
            "saliva, tissue, or auto (default: blood). Blood context includes blood + plasma "
            "and excludes serum unless --additional-contexts includes serum."
        ),
    )
    p_map.add_argument(
        "--additional-contexts",
        default="",
        help=(
            "Optional extra contexts allowed in addition to --measurement-context, comma-separated "
            "(for example: serum or cerebrospinal_fluid). Default is none."
        ),
    )
    p_map.add_argument(
        "--additional-context-keywords",
        default="",
        help=(
            "Optional free-text context keywords/phrases, comma-separated (for example: aorta, "
            "adipose tissue). Terms matching these phrases are allowed in addition to controlled "
            "context filters. If --measurement-context auto is used, these keywords become the "
            "primary context filter."
        ),
    )
    p_map.add_argument(
        "--name-mode",
        choices=["strict", "fuzzy"],
        default="strict",
        help=(
            "Handling for free-text names: strict requires unique alias resolution "
            "(UniProt for proteins, concept alias for metabolites); fuzzy enables lexical fallback"
        ),
    )
    p_map.add_argument(
        "--entity-type",
        choices=["auto", "protein", "metabolite"],
        default="auto",
        help=(
            "Entity mapping mode. auto routes per query; protein uses UniProt identity gating; "
            "metabolite uses metabolite-concept identity gating and blocks protein-like queries."
        ),
    )
    p_map.add_argument(
        "--force-map-best",
        action="store_true",
        help="When no candidate passes min-score, emit the single best candidate as not_validated",
    )
    p_map.add_argument(
        "--fallback-efo-id",
        default="",
        help="Optional generic fallback ID (for example EFO:0001444) used when force-map-best is set and no candidate exists",
    )
    p_map.add_argument(
        "--auto-enrich-uniprot",
        action="store_true",
        help="Fetch missing UniProt accession aliases for current input, update alias cache, and rebuild index",
    )
    p_map.add_argument(
        "--uniprot-aliases",
        default=str(DEFAULT_UNIPROT_ALIASES),
        help="UniProt alias TSV path used for enrichment and index rebuild",
    )
    p_map.add_argument(
        "--metabolite-aliases",
        default=str(DEFAULT_METABOLITE_ALIASES),
        help="Metabolite alias TSV used for metabolite concept resolution and index rebuild",
    )
    p_map.add_argument(
        "--term-cache",
        default=str(DEFAULT_TERM_CACHE),
        help="Term cache TSV used when rebuilding index after enrichment",
    )
    p_map.add_argument(
        "--uniprot-timeout",
        type=float,
        default=20.0,
        help="Timeout for UniProt alias enrichment requests",
    )
    p_map.add_argument(
        "--uniprot-workers",
        type=int,
        default=8,
        help="Parallel workers for UniProt alias enrichment",
    )
    p_map.add_argument(
        "--uniprot-source",
        default="uniprot-live-auto",
        help="Source tag for auto-enriched alias rows",
    )
    p_map.add_argument(
        "--unmapped-output",
        help="Optional TSV path to write unresolved not_mapped rows",
    )
    p_map.add_argument(
        "--review-output",
        help=(
            "Optional TSV path for candidates withheld from auto-validation "
            "(reason_code in {manual_validation_required, needs_new_term, "
            "suggest_broad_term_mapping, metabolite_mode_blocked_protein_like_query})"
        ),
    )
    p_map.add_argument(
        "--review-queue-output",
        help=(
            "Optional TSV path for broader non-validated review queue with reason codes "
            "and suggestions"
        ),
    )
    p_map.add_argument(
        "--llm-review-output",
        help=argparse.SUPPRESS,
    )
    p_map.add_argument(
        "--withheld-triage-output",
        help=(
            "Optional TSV path to triage top withheld candidates into accept/reject/review "
            "buckets using accession identity checks. If --review-output is not set, a "
            "temporary internal review file is generated automatically."
        ),
    )
    p_map.add_argument(
        "--review-top-n",
        type=int,
        default=3,
        help="How many suggested candidate terms to include per review row (default: 3)",
    )
    p_map.add_argument("--cache-writeback", action="store_true", help="Write validated hits into analyte cache")
    p_map.add_argument("--analyte-cache", default=str(DEFAULT_ANALYTE_CACHE), help="Analyte cache path for write-back")
    p_map.add_argument("--cache-min-confidence", type=float, default=0.80, help="Min confidence for write-back")

    p_trait_cache = sub.add_parser(
        "trait-cache-refresh",
        help=(
            "Refresh catalog-derived trait cache rows from an updated GWAS Catalog studies TSV "
            "using high-confidence, unambiguous, non-conflicting filters"
        ),
    )
    p_trait_cache.add_argument(
        "--studies-tsv",
        required=True,
        help="GWAS Catalog studies TSV (download export with DISEASE/TRAIT, MAPPED_TRAIT, MAPPED_TRAIT_URI columns)",
    )
    p_trait_cache.add_argument(
        "--trait-cache",
        default=str(DEFAULT_TRAIT_CACHE),
        help="Existing trait mapping cache TSV to update",
    )
    p_trait_cache.add_argument(
        "--output-cache",
        default="",
        help="Optional output cache path (defaults to --trait-cache for in-place update)",
    )
    p_trait_cache.add_argument(
        "--publication-tag",
        default="CatalogMappedTraitsImport",
        help="Publication tag for inserted cache rows",
    )
    p_trait_cache.add_argument(
        "--efo-obo",
        default=str(DEFAULT_EFO_OBO_LOCAL),
        help="EFO OBO path used for ID resolution and MONDO-import checks",
    )
    p_trait_cache.add_argument(
        "--efo-obo-bundled-url",
        default=DEFAULT_EFO_OBO_BUNDLED_URL,
        help="Bundled EFO OBO URL used only when --efo-obo is missing",
    )
    p_trait_cache.add_argument(
        "--efo-obo-download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) when auto-fetching bundled EFO OBO",
    )
    p_trait_cache.add_argument(
        "--qc-output",
        default=str(DEFAULT_CATALOG_QC_OUTPUT),
        help="Optional QC TSV output path (set empty to skip)",
    )
    p_trait_cache.add_argument(
        "--high-confidence-output",
        default=str(DEFAULT_CATALOG_HIGH_CONF_OUTPUT),
        help="Optional high-confidence TSV output path (set empty to skip)",
    )
    p_trait_cache.add_argument(
        "--to-add-output",
        default=str(DEFAULT_CATALOG_TO_ADD_OUTPUT),
        help="Optional cache-row TSV output path for rows selected to add (set empty to skip)",
    )
    p_trait_cache.add_argument(
        "--summary-output",
        default=str(DEFAULT_CATALOG_SUMMARY_OUTPUT),
        help="Optional JSON summary output path (set empty to skip)",
    )
    p_trait_cache.add_argument(
        "--state-output",
        default=str(DEFAULT_CATALOG_REFRESH_STATE_OUTPUT),
        help="State JSON path used to track last processed studies TSV hash (set empty to skip)",
    )
    p_trait_cache.add_argument(
        "--force-refresh",
        action="store_true",
        help="Run refresh even when studies TSV hash matches last recorded state",
    )

    p_trait = sub.add_parser(
        "trait-map",
        help="Map disease/phenotype traits using curated cache first, then efo.obo fallback",
    )
    p_trait.add_argument(
        "--input",
        required=True,
        help=(
            "Input txt/csv/tsv for disease/phenotype mapping. Supported columns: "
            "query/trait/reported_trait, optional icd10 and/or phecode, optional input_type, "
            "optional trait_scale (binary|quantitative), optional code and data_type "
            "(for example code+data_type=UKB Data Field or ICD10). "
            "Auto-routing handles ICD10 strings (including union-style strings like 'Union#A071#A07.1'), "
            "UKB field text (for example 'UKB data field 22435'), PheCode-like values, and free text."
        ),
    )
    p_trait.add_argument("--output", required=True, help="Output TSV path")
    p_trait.add_argument(
        "--trait-cache",
        default=str(DEFAULT_TRAIT_CACHE),
        help=(
            "Curated trait mapping cache TSV (default: "
            "skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv)"
        ),
    )
    p_trait.add_argument(
        "--efo-obo",
        default=str(DEFAULT_EFO_OBO_LOCAL),
        help=(
            "EFO OBO path used for fallback label/synonym mapping. If missing, trait-map will "
            "auto-provision from --efo-obo-bundled-url."
        ),
    )
    p_trait.add_argument(
        "--efo-obo-bundled-url",
        default=DEFAULT_EFO_OBO_BUNDLED_URL,
        help=(
            "Pinned bundled EFO OBO URL (recommended: your repo release asset, e.g. efo.obo.gz). "
            "Used only when --efo-obo is missing."
        ),
    )
    p_trait.add_argument(
        "--efo-obo-download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) when auto-fetching bundled EFO OBO",
    )
    p_trait.add_argument("--top-k", type=int, default=1, help="Candidates to emit per input row")
    p_trait.add_argument("--min-score", type=float, default=0.82, help="Minimum confidence score to accept candidate")
    p_trait.add_argument(
        "--default-trait-scale",
        default="auto",
        choices=["auto", "binary", "quantitative"],
        help=(
            "Default trait scale when no trait_scale column is present. "
            "binary biases to non-measurement branch; quantitative biases to measurement branch."
        ),
    )
    p_trait.add_argument(
        "--force-map-best",
        action="store_true",
        help="Emit best candidate as review_required when no mapping reaches --min-score",
    )
    p_trait.add_argument(
        "--review-output",
        help="Optional TSV path for rows with validation=review_required",
    )
    p_trait.add_argument(
        "--qc-risk-output",
        help=(
            "Optional TSV path for curator-risk flags. "
            "If omitted, trait-map writes <output stem>_qc_risk.tsv."
        ),
    )
    p_trait.add_argument(
        "--ukb-field-catalog",
        default=str(DEFAULT_UKB_FIELD_CATALOG),
        help=(
            "Optional UKB field catalog TSV (for example references/ukb/fieldsum.txt) "
            "used to enrich UKB data-field queries by field title."
        ),
    )
    p_trait.add_argument(
        "--icd10-supplement-cache",
        default=str(DEFAULT_ICD10_SUPPLEMENT_CACHE),
        help=(
            "Optional ICD10 supplemental cache TSV (generated by setup-bundled-caches) "
            "used for ICD10 exact matches not present in the curated trait cache."
        ),
    )
    p_trait.add_argument(
        "--icd10-label-cache",
        default=str(DEFAULT_ICD10_LABEL_CACHE),
        help=(
            "Optional ICD10 code-label cache TSV (for example built with icd10-label-cache-build) "
            "used to populate input_icd10_label for bare ICD10 code queries."
        ),
    )
    p_trait.add_argument(
        "--stream-output",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write trait-map rows incrementally while processing (default: on)",
    )
    p_trait.add_argument(
        "--flush-every",
        type=int,
        default=50,
        help=(
            "When --stream-output is on, flush output files every N written rows "
            "(lower is safer for interruptions, default: 50)"
        ),
    )
    p_trait.add_argument(
        "--memoize-queries",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Reuse mapping results for duplicate normalized trait inputs (default: on)",
    )
    p_trait.add_argument(
        "--query-cache-max-entries",
        type=int,
        default=200000,
        help="Maximum duplicate-query memoization entries (default: 200000)",
    )
    p_trait.add_argument(
        "--progress",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show progress logs (default: on)",
    )

    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 2

    try:
        if args.command == "icd10-label-cache-build":
            output_path = Path(args.output)
            stats = build_icd10_label_cache(
                output_path=output_path,
                source_url=args.source_url,
                timeout=float(args.download_timeout),
            )
            print(
                f"[OK] built ICD10 label cache at {output_path} "
                f"(codes={stats.get('codes', 0)}, "
                f"source_member={stats.get('source_member', '')}, "
                f"archive_size_bytes={stats.get('archive_size_bytes', 0)}, "
                f"source_url={stats.get('source_url', '')})"
            )
            return 0

        if args.command == "setup-bundled-caches":
            term_cache_path = Path(args.term_cache)
            analyte_cache_path = Path(args.analyte_cache)
            uniprot_alias_source_path = Path(args.uniprot_aliases)
            uniprot_profile = normalize(args.uniprot_profile)
            uniprot_light_output_path = Path(args.uniprot_light_output)
            metabolite_alias_path = Path(args.metabolite_aliases)
            trait_cache_path = Path(args.trait_cache)
            ukb_setup_catalog_arg = normalize(args.ukb_field_catalog)
            ukb_setup_catalog_path = Path(args.ukb_field_catalog) if ukb_setup_catalog_arg else None
            ukb_field_supplement_cache_path = Path(args.ukb_field_supplement_cache)
            icd10_supplement_cache_path = Path(args.icd10_supplement_cache)
            icd10_label_cache_path = Path(args.icd10_label_cache)
            mondo_sssom_cache_path = Path(args.icd10_mondo_sssom_cache)
            efo_obo_setup_path = Path(args.efo_obo)
            index_path = Path(args.index)
            cache_manifest_arg = normalize(args.cache_manifest_output)
            cache_manifest_output_path = Path(args.cache_manifest_output) if cache_manifest_arg else None

            if args.seed_legacy_trait_cache:
                trait_cache_path, seeded = seed_default_trait_cache_if_needed(trait_cache_path)
                if seeded:
                    print(
                        f"[OK] seeded bundled trait cache at {trait_cache_path} "
                        f"from legacy file {LEGACY_TRAIT_CACHE}"
                    )

            efo_obo_setup_path, efo_obo_source = ensure_efo_obo_available(
                efo_obo_path=efo_obo_setup_path,
                bundled_url=args.efo_obo_bundled_url,
                timeout=float(args.efo_obo_download_timeout),
            )
            if efo_obo_source != "local":
                print(
                    f"[OK] provisioned EFO OBO at {efo_obo_setup_path} "
                    f"(source={efo_obo_source})"
                )

            required_paths = [
                ("term_cache", term_cache_path),
                ("analyte_cache", analyte_cache_path),
                ("uniprot_aliases", uniprot_alias_source_path),
                ("metabolite_aliases", metabolite_alias_path),
                ("trait_cache", trait_cache_path),
                ("efo_obo", efo_obo_setup_path),
            ]

            ukb_roots: list[Path] = []
            if ukb_setup_catalog_path is not None:
                ukb_roots.append(ukb_setup_catalog_path.parent)
            default_ukb_root = DEFAULT_UKB_FIELD_CATALOG.parent
            if default_ukb_root not in ukb_roots:
                ukb_roots.append(default_ukb_root)

            def pick_ukb_metadata_path(filename: str) -> Path | None:
                for root in ukb_roots:
                    candidate = root / filename
                    if candidate.exists():
                        return candidate
                if ukb_roots:
                    return ukb_roots[0] / filename
                return None

            ukb_field_catalog_path = ukb_setup_catalog_path or DEFAULT_UKB_FIELD_CATALOG
            ukb_field_metadata_path = pick_ukb_metadata_path("field.txt") or DEFAULT_UKB_FIELD_METADATA
            ukb_category_catalog_path = pick_ukb_metadata_path("category.txt") or DEFAULT_UKB_CATEGORY_CATALOG
            ukb_category_tree_path = pick_ukb_metadata_path("catbrowse.txt") or DEFAULT_UKB_CATEGORY_TREE
            if bool(args.build_ukb_field_supplement):
                ukb_provisioned = ensure_ukb_metadata_available(
                    ukb_field_catalog_path=ukb_field_catalog_path,
                    ukb_field_metadata_path=ukb_field_metadata_path,
                    ukb_category_catalog_path=ukb_category_catalog_path,
                    ukb_category_tree_path=ukb_category_tree_path,
                    field_catalog_url=args.ukb_field_catalog_url,
                    field_metadata_url=args.ukb_field_metadata_url,
                    category_catalog_url=args.ukb_category_catalog_url,
                    category_tree_url=args.ukb_category_tree_url,
                    timeout=float(args.ukb_download_timeout),
                )
                for cache_name, (cache_path, source) in ukb_provisioned.items():
                    if source != "local":
                        print(f"[OK] provisioned {cache_name} at {cache_path} (source={source})")
                required_paths.extend(
                    [
                        ("ukb_field_catalog", ukb_field_catalog_path),
                        ("ukb_field_metadata", ukb_field_metadata_path),
                        ("ukb_category_catalog", ukb_category_catalog_path),
                        ("ukb_category_tree", ukb_category_tree_path),
                    ]
                )
            missing = [f"{name}={path}" for name, path in required_paths if not path.exists()]
            if missing:
                raise FileNotFoundError(
                    "Missing setup inputs after bundled checks/provisioning: "
                    + ", ".join(missing)
                )

            term_rows = load_term_cache(term_cache_path)
            analyte_rows = load_analyte_cache(analyte_cache_path)
            uniprot_alias_index_path = uniprot_alias_source_path
            light_total = 0
            light_kept = 0
            if uniprot_profile == "light":
                light_total, light_kept = build_lightweight_uniprot_alias_cache(
                    source_path=uniprot_alias_source_path,
                    output_path=uniprot_light_output_path,
                    reviewed_only=bool(args.uniprot_light_reviewed_only),
                    require_gene=bool(args.uniprot_light_require_gene),
                )
                uniprot_alias_index_path = uniprot_light_output_path
                print(
                    f"[OK] built lightweight UniProt cache at {uniprot_alias_index_path} "
                    f"(source_rows={light_total}, kept_rows={light_kept}, "
                    f"reviewed_only={bool(args.uniprot_light_reviewed_only)}, "
                    f"require_gene={bool(args.uniprot_light_require_gene)})"
                )
                if args.uniprot_light_enrich_gene_ids:
                    updated, failed, targets = backfill_uniprot_gene_ids(
                        alias_path=uniprot_alias_index_path,
                        timeout=float(args.uniprot_timeout),
                        workers=max(1, int(args.uniprot_workers)),
                        source_tag=normalize(args.uniprot_source),
                        refresh_all=False,
                    )
                    print(
                        f"[OK] lightweight UniProt gene_id enrichment complete "
                        f"(targets={targets}, updated={updated}, failed={failed})"
                    )

            uniprot_aliases = load_uniprot_aliases(uniprot_alias_index_path)
            metabolite_alias_records = load_metabolite_alias_records(metabolite_alias_path)
            if bool(args.build_ukb_field_supplement):
                if (
                    ukb_field_supplement_cache_path.exists()
                    and ukb_field_supplement_cache_path.stat().st_size > 0
                ):
                    print(
                        f"[OK] using existing UKB supplemental trait cache at "
                        f"{ukb_field_supplement_cache_path}"
                    )
                else:
                    supplement_stats = build_ukb_field_supplement_cache(
                        output_path=ukb_field_supplement_cache_path,
                        trait_cache_path=trait_cache_path,
                        efo_obo_path=efo_obo_setup_path,
                        ukb_field_catalog_path=ukb_setup_catalog_path,
                        ukb_field_metadata_path=ukb_field_metadata_path,
                        ukb_category_catalog_path=ukb_category_catalog_path,
                        ukb_category_tree_path=ukb_category_tree_path,
                    )
                    print(
                        f"[OK] built UKB supplemental trait cache at {ukb_field_supplement_cache_path} "
                        f"(fields_total={supplement_stats.get('total_fields', 0)}, "
                        f"base_covered={supplement_stats.get('base_covered_fields', 0)}, "
                        f"supplement_rows={supplement_stats.get('supplement_rows', 0)}, "
                        f"field_title_cache_exact={supplement_stats.get('method_field_title_cache_exact', 0)}, "
                        f"field_title_ontology_exact={supplement_stats.get('method_field_title_ontology_exact', 0)}, "
                        f"field_title_ontology_fuzzy={supplement_stats.get('method_field_title_ontology_fuzzy', 0)}, "
                        f"category_cache_unique={supplement_stats.get('method_category_cache_unique', 0)}, "
                        f"category_text_exact={supplement_stats.get('method_category_text_exact', 0)}, "
                        f"category_ontology_exact={supplement_stats.get('method_category_ontology_exact', 0)}, "
                        f"generic_disease={supplement_stats.get('method_generic_disease', 0)}, "
                        f"generic_measurement={supplement_stats.get('method_generic_measurement', 0)})"
                    )

            if bool(args.build_icd10_supplement):
                if (
                    icd10_supplement_cache_path.exists()
                    and icd10_supplement_cache_path.stat().st_size > 0
                ):
                    print(
                        f"[OK] using existing ICD10 supplemental trait cache at "
                        f"{icd10_supplement_cache_path}"
                    )
                else:
                    icd10_supplement_stats = build_icd10_supplement_cache(
                        output_path=icd10_supplement_cache_path,
                        trait_cache_path=trait_cache_path,
                        efo_obo_path=efo_obo_setup_path,
                        mondo_sssom_cache_path=mondo_sssom_cache_path,
                        mondo_sssom_url=args.icd10_mondo_sssom_url,
                        mondo_download_timeout=float(args.icd10_mondo_download_timeout),
                    )
                    print(
                        f"[OK] built ICD10 supplemental trait cache at {icd10_supplement_cache_path} "
                        f"(codes_total={icd10_supplement_stats.get('total_codes', 0)}, "
                        f"base_covered={icd10_supplement_stats.get('base_covered_codes', 0)}, "
                        f"supplement_rows={icd10_supplement_stats.get('supplement_rows', 0)}, "
                        f"from_efo_obo_xref={icd10_supplement_stats.get('method_efo_obo_xref', 0)}, "
                        f"from_mondo_sssom={icd10_supplement_stats.get('method_mondo_sssom', 0)}, "
                        f"from_both={icd10_supplement_stats.get('method_mondo_sssom_efo_obo_xref', 0)}, "
                        f"mondo_cache_source={icd10_supplement_stats.get('mondo_cache_source', 'none')}, "
                        f"mondo_mapping_codes={icd10_supplement_stats.get('mondo_mapping_codes', 0)}, "
                        f"mondo_downloaded={icd10_supplement_stats.get('mondo_downloaded', False)}, "
                        f"mondo_codes_added={icd10_supplement_stats.get('mondo_codes_added', 0)}, "
                        f"mondo_terms_not_in_efo={icd10_supplement_stats.get('mondo_terms_not_in_efo', 0)})"
                    )
                    if icd10_supplement_stats.get("mondo_download_error"):
                        if icd10_supplement_stats.get("mondo_cache_used"):
                            print(
                                "[WARN] MONDO SSSOM live refresh failed; using local cached MONDO mappings. "
                                f"reason={icd10_supplement_stats.get('mondo_download_error')}"
                            )
                        else:
                            print(
                                "[WARN] MONDO SSSOM mappings unavailable; ICD10 supplement built from EFO ICD10 xrefs only. "
                                f"reason={icd10_supplement_stats.get('mondo_download_error')}"
                            )

            if bool(args.build_icd10_label_cache):
                existing_icd10_labels = load_icd10_label_cache(icd10_label_cache_path)
                if existing_icd10_labels:
                    print(
                        f"[OK] using existing ICD10 label cache at {icd10_label_cache_path} "
                        f"(codes={len(existing_icd10_labels)})"
                    )
                else:
                    icd10_label_stats = build_icd10_label_cache(
                        output_path=icd10_label_cache_path,
                        source_url=args.icd10_label_source_url,
                        timeout=float(args.icd10_label_download_timeout),
                    )
                    print(
                        f"[OK] built ICD10 label cache at {icd10_label_cache_path} "
                        f"(codes={icd10_label_stats.get('codes', 0)}, "
                        f"source_member={icd10_label_stats.get('source_member', '')})"
                    )

            extra_trait_cache_paths: list[Path] = []
            if ukb_field_supplement_cache_path.exists():
                extra_trait_cache_paths.append(ukb_field_supplement_cache_path)
            if icd10_supplement_cache_path.exists():
                extra_trait_cache_paths.append(icd10_supplement_cache_path)
            trait_cache_index = load_trait_cache_index(trait_cache_path, extra_paths=extra_trait_cache_paths)

            if args.rebuild_index or not index_path.exists():
                index = build_index(
                    term_rows=term_rows,
                    analyte_rows=analyte_rows,
                    accession_aliases=uniprot_aliases,
                    metabolite_alias_records=metabolite_alias_records,
                )
                save_index(index, index_path)
                print(
                    f"[OK] built index at {index_path} "
                    f"(terms={len(index.get('term_meta', {}))}, "
                    f"accessions={len(index.get('accession_alias_index', {}))}, "
                    f"gene_id_keys={len(index.get('gene_id_accession_index', {}))}, "
                    f"metabolite_concepts={len(index.get('metabolite_concept_index', {}))})"
                )
            else:
                index = load_index(index_path)
                print(
                    f"[OK] using existing index at {index_path} "
                    f"(terms={len(index.get('term_meta', {}))})"
                )

            print(
                "[OK] bundled cache setup complete "
                f"(uniprot_profile={uniprot_profile}, "
                f"uniprot_alias_path={uniprot_alias_index_path}, "
                f"protein_aliases={len(uniprot_aliases)}, "
                f"metabolite_concepts={len(metabolite_alias_records)}, "
                f"trait_records={len(trait_cache_index.get('records', []))}, "
                f"ukb_field_supplement={ukb_field_supplement_cache_path}, "
                f"icd10_supplement={icd10_supplement_cache_path}, "
                f"icd10_label_cache={icd10_label_cache_path}, "
                f"efo_obo={efo_obo_setup_path})"
            )
            if cache_manifest_output_path is not None:
                write_setup_cache_manifest(
                    output_path=cache_manifest_output_path,
                    index_path=index_path,
                    index=index,
                    file_entries=[
                        ("term_cache", term_cache_path, True),
                        ("analyte_cache", analyte_cache_path, True),
                        ("uniprot_aliases_source", uniprot_alias_source_path, False),
                        ("uniprot_aliases_index", uniprot_alias_index_path, True),
                        ("metabolite_aliases", metabolite_alias_path, True),
                        ("hmdb_source_dir", DEFAULT_HMDB_LOCAL_DIR, False),
                        ("metabolite_download_dir", DEFAULT_METABOLITE_DOWNLOAD_DIR, False),
                        ("trait_cache", trait_cache_path, False),
                        ("ukb_field_catalog", ukb_field_catalog_path, False),
                        ("ukb_field_metadata", ukb_field_metadata_path, False),
                        ("ukb_category_catalog", ukb_category_catalog_path, False),
                        ("ukb_category_tree", ukb_category_tree_path, False),
                        ("ukb_field_supplement_cache", ukb_field_supplement_cache_path, False),
                        ("icd10_supplement_cache", icd10_supplement_cache_path, False),
                        ("icd10_label_cache", icd10_label_cache_path, False),
                        ("mondo_sssom_cache", mondo_sssom_cache_path, False),
                        ("efo_obo", efo_obo_setup_path, False),
                    ],
                    generated_paths=[
                        index_path,
                        ukb_field_supplement_cache_path,
                        icd10_supplement_cache_path,
                        icd10_label_cache_path,
                        mondo_sssom_cache_path,
                    ],
                )
                print(f"[OK] wrote setup cache manifest to {cache_manifest_output_path}")
            return 0

        if args.command == "cache-status":
            manifest_path = Path(args.manifest)
            manifest_meta: dict[str, Any] = {
                "path": str(manifest_path),
                "exists": manifest_path.exists(),
                "loaded": False,
                "generated_at_utc": "",
                "error": "",
                "index_summary": {},
            }
            specs: list[CacheStatusSpec] = []
            if bool(args.prefer_manifest):
                specs, manifest_meta = load_cache_status_specs_from_manifest(manifest_path)
            if not specs:
                ukb_catalog_arg = normalize(args.ukb_field_catalog)
                ukb_catalog_path = Path(args.ukb_field_catalog) if ukb_catalog_arg else DEFAULT_UKB_FIELD_CATALOG
                (
                    ukb_catalog_path,
                    ukb_field_metadata_path,
                    ukb_category_catalog_path,
                    ukb_category_tree_path,
                ) = setup_ukb_paths(ukb_catalog_path)
                specs = collect_default_cache_status_specs(
                    term_cache_path=Path(args.term_cache),
                    analyte_cache_path=Path(args.analyte_cache),
                    uniprot_alias_source_path=Path(args.uniprot_aliases),
                    uniprot_alias_index_path=Path(args.uniprot_light_output),
                    metabolite_alias_path=Path(args.metabolite_aliases),
                    trait_cache_path=Path(args.trait_cache),
                    ukb_field_catalog_path=ukb_catalog_path,
                    ukb_field_metadata_path=ukb_field_metadata_path,
                    ukb_category_catalog_path=ukb_category_catalog_path,
                    ukb_category_tree_path=ukb_category_tree_path,
                    ukb_field_supplement_cache_path=Path(args.ukb_field_supplement_cache),
                    icd10_supplement_cache_path=Path(args.icd10_supplement_cache),
                    icd10_label_cache_path=Path(args.icd10_label_cache),
                    mondo_sssom_cache_path=Path(args.icd10_mondo_sssom_cache),
                    efo_obo_path=Path(args.efo_obo),
                    index_path=Path(args.index),
                )
            report = evaluate_cache_status_specs(specs, manifest_meta=manifest_meta)
            print_cache_status_report(report)

            output_json_arg = normalize(args.output_json)
            if output_json_arg:
                output_json_path = Path(args.output_json)
                output_json_path.parent.mkdir(parents=True, exist_ok=True)
                output_json_path.write_text(
                    json.dumps(report, indent=2, ensure_ascii=True) + "\n",
                    encoding="utf-8",
                )
                print(f"[OK] wrote cache status JSON to {output_json_path}")

            required_missing = int(report.get("summary", {}).get("required_missing", 0))
            if bool(args.strict) and required_missing > 0:
                print(
                    f"[ERROR] cache status check failed: required_missing={required_missing}",
                    file=sys.stderr,
                )
                return 1
            return 0

        if args.command == "uniprot-alias-build-light":
            source_path = Path(args.source)
            output_path = Path(args.output)
            total, kept = build_lightweight_uniprot_alias_cache(
                source_path=source_path,
                output_path=output_path,
                reviewed_only=bool(args.reviewed_only),
                require_gene=bool(args.require_gene),
            )
            print(
                f"[OK] wrote lightweight UniProt cache to {output_path} "
                f"(source_rows={total}, kept_rows={kept}, reviewed_only={bool(args.reviewed_only)}, "
                f"require_gene={bool(args.require_gene)})"
            )
            return 0

        if args.command == "uniprot-alias-enrich":
            updated, failed, targets = backfill_uniprot_gene_ids(
                alias_path=Path(args.uniprot_aliases),
                timeout=float(args.timeout),
                workers=max(1, int(args.workers)),
                source_tag=normalize(args.source),
                refresh_all=bool(args.refresh_all),
            )
            print(
                f"[OK] UniProt alias enrichment complete (targets={targets}, updated={updated}, failed={failed}) "
                f"-> {args.uniprot_aliases}"
            )
            return 0

        if args.command == "index-build":
            term_rows = load_term_cache(Path(args.term_cache))
            analyte_rows = load_analyte_cache(Path(args.analyte_cache))
            aliases = load_uniprot_aliases(Path(args.uniprot_aliases))
            metabolite_aliases = load_metabolite_alias_records(Path(args.metabolite_aliases))
            index = build_index(
                term_rows=term_rows,
                analyte_rows=analyte_rows,
                accession_aliases=aliases,
                metabolite_alias_records=metabolite_aliases,
            )
            save_index(index, Path(args.output_index))
            print(
                f"[OK] built index at {args.output_index} "
                f"(terms={len(index.get('term_meta', {}))}, "
                f"analyte_keys={len(index.get('analyte_index', {}))}, "
                f"accessions={len(index.get('accession_alias_index', {}))}, "
                f"gene_id_keys={len(index.get('gene_id_accession_index', {}))}, "
                f"metabolite_concepts={len(index.get('metabolite_concept_index', {}))})"
            )
            return 0

        if args.command == "refresh-efo-cache":
            run_efo_cache_refresh(
                efo_obo=args.efo_obo,
                download_url=args.download_url,
                download_to=args.download_to,
                measurement_root=args.measurement_root,
                term_cache=args.term_cache,
                source=args.source,
                timeout=args.timeout,
            )
            print(f"[OK] refreshed term cache at {args.term_cache}")

            if args.rebuild_index:
                term_rows = load_term_cache(Path(args.term_cache))
                analyte_rows = load_analyte_cache(Path(args.analyte_cache))
                aliases = load_uniprot_aliases(Path(args.uniprot_aliases))
                metabolite_aliases = load_metabolite_alias_records(Path(args.metabolite_aliases))
                index = build_index(
                    term_rows=term_rows,
                    analyte_rows=analyte_rows,
                    accession_aliases=aliases,
                    metabolite_alias_records=metabolite_aliases,
                )
                save_index(index, Path(args.index))
                print(
                    f"[OK] rebuilt index at {args.index} "
                    f"(terms={len(index.get('term_meta', {}))}, "
                    f"analyte_keys={len(index.get('analyte_index', {}))}, "
                    f"accessions={len(index.get('accession_alias_index', {}))}, "
                    f"gene_id_keys={len(index.get('gene_id_accession_index', {}))}, "
                    f"metabolite_concepts={len(index.get('metabolite_concept_index', {}))})"
                )
            return 0

        if args.command == "metabolite-alias-build":
            output_path = Path(args.output)
            hmdb_path = Path(args.hmdb_xml) if normalize(args.hmdb_xml) else None
            download_dir = Path(args.download_dir)
            if hmdb_path is None:
                hmdb_path = detect_local_hmdb_source(
                    download_dir=download_dir,
                    local_dir=DEFAULT_HMDB_LOCAL_DIR,
                )
                if hmdb_path is not None:
                    print(f"[OK] using local HMDB source -> {hmdb_path}")
            chebi_path = Path(args.chebi_names) if normalize(args.chebi_names) else None
            kegg_path = Path(args.kegg_conv) if normalize(args.kegg_conv) else None
            auto_download_hmdb = bool(args.download_hmdb) or (
                hmdb_path is None and chebi_path is None and kegg_path is None
            )
            hmdb_path, chebi_path, kegg_path, download_notes = prepare_metabolite_source_paths(
                hmdb_xml_path=hmdb_path,
                chebi_names_path=chebi_path,
                kegg_conv_path=kegg_path,
                download_common_sources=bool(args.download_common_sources),
                download_hmdb=auto_download_hmdb,
                download_dir=download_dir,
                hmdb_url=args.hmdb_url,
                chebi_names_url=args.chebi_names_url,
                kegg_conv_url=args.kegg_conv_url,
                timeout=float(args.download_timeout),
            )
            if hmdb_path is None and chebi_path is None and kegg_path is None:
                raise ValueError(
                    "No metabolite sources selected. Provide --hmdb-xml <local-HMDB-file> "
                    "(recommended, pinned HMDB version), or --chebi-names/--kegg-conv for optional "
                    "extra resources."
                )
            for source_path in [hmdb_path, chebi_path, kegg_path]:
                if source_path is not None and not source_path.exists():
                    raise FileNotFoundError(f"Metabolite source file not found: {source_path}")
            hmdb_rows, chebi_rows, kegg_rows, concept_count = build_metabolite_alias_records(
                output_path=output_path,
                hmdb_xml_path=hmdb_path,
                chebi_names_path=chebi_path,
                kegg_conv_path=kegg_path,
                merge_existing=bool(args.merge_existing),
            )
            for note in download_notes:
                print(f"[OK] {note}")
            print(
                f"[OK] wrote metabolite aliases to {output_path} "
                f"(hmdb_rows={hmdb_rows}, chebi_rows={chebi_rows}, kegg_rows={kegg_rows}, concepts={concept_count})"
            )
            return 0

        if args.command == "map":
            query_inputs = load_query_inputs(Path(args.input))
            if not query_inputs:
                raise ValueError("No usable queries found in input")
            queries = [query for query, _input_type in query_inputs]
            measurement_context = parse_measurement_context(args.measurement_context)
            additional_contexts = parse_additional_contexts(args.additional_contexts)
            additional_context_keywords = parse_additional_context_keywords(args.additional_context_keywords)
            index_path = Path(args.index)
            analyte_cache_path = Path(args.analyte_cache)
            uniprot_alias_path = Path(args.uniprot_aliases)
            metabolite_alias_path = Path(args.metabolite_aliases)
            term_cache_path = Path(args.term_cache)
            entity_type = normalize_entity_type(args.entity_type)
            matrix_priority = parse_matrix_priority(args.matrix_priority)

            enriched_added = 0
            enriched_failed = 0
            if args.auto_enrich_uniprot:
                enriched_added, enriched_failed = enrich_uniprot_aliases_for_queries(
                    queries,
                    alias_path=uniprot_alias_path,
                    timeout=args.uniprot_timeout,
                    workers=max(1, args.uniprot_workers),
                    source_tag=args.uniprot_source,
                )
                print(
                    f"[OK] UniProt auto-enrich added {enriched_added} aliases "
                    f"(failed={enriched_failed}) into {uniprot_alias_path}"
                )

            # Rebuild index when missing or after new alias enrichment.
            if enriched_added > 0 or not index_path.exists():
                term_rows = load_term_cache(term_cache_path)
                analyte_rows = load_analyte_cache(analyte_cache_path)
                aliases = load_uniprot_aliases(uniprot_alias_path)
                metabolite_aliases = load_metabolite_alias_records(metabolite_alias_path)
                index = build_index(
                    term_rows=term_rows,
                    analyte_rows=analyte_rows,
                    accession_aliases=aliases,
                    metabolite_alias_records=metabolite_aliases,
                )
                save_index(index, index_path)
                print(
                    f"[OK] rebuilt index at {index_path} "
                    f"(terms={len(index.get('term_meta', {}))}, "
                    f"analyte_keys={len(index.get('analyte_index', {}))}, "
                    f"accessions={len(index.get('accession_alias_index', {}))}, "
                    f"gene_id_keys={len(index.get('gene_id_accession_index', {}))}, "
                    f"metabolite_concepts={len(index.get('metabolite_concept_index', {}))})"
                )
            else:
                index = load_index(index_path)

            rows = map_queries(
                query_inputs,
                index=index,
                top_k=max(1, args.top_k),
                min_score=args.min_score,
                workers=max(1, args.workers),
                matrix_priority=matrix_priority,
                measurement_context=measurement_context,
                additional_contexts=additional_contexts,
                additional_context_keywords=additional_context_keywords,
                name_mode=args.name_mode,
                force_map_best=args.force_map_best,
                fallback_efo_id=args.fallback_efo_id,
                entity_type=entity_type,
                show_progress=args.progress,
                parallel_mode=args.parallel_mode,
                index_path=index_path,
            )
            write_tsv(rows, Path(args.output))
            if args.unmapped_output:
                unresolved = write_unmapped_tsv(rows, Path(args.unmapped_output))
                print(f"[OK] wrote {unresolved} unresolved rows to {args.unmapped_output}")
            review_output_arg = normalize(args.review_output or "")
            legacy_llm_review_output_arg = normalize(args.llm_review_output or "")
            review_output_path: Path | None = None
            if review_output_arg and legacy_llm_review_output_arg and review_output_arg != legacy_llm_review_output_arg:
                raise ValueError(
                    "--review-output and --llm-review-output were both set with different paths; "
                    "use one path only"
                )
            if review_output_arg:
                review_output_path = Path(args.review_output)
            elif legacy_llm_review_output_arg:
                review_output_path = Path(args.llm_review_output)
                print("[WARN] --llm-review-output is deprecated; use --review-output", file=sys.stderr)

            if review_output_path is not None:
                review_rows = write_review_tsv(
                    rows,
                    review_output_path,
                    index=index,
                    min_score=args.min_score,
                    matrix_priority=matrix_priority,
                    measurement_context=measurement_context,
                    additional_contexts=additional_contexts,
                    additional_context_keywords=additional_context_keywords,
                    name_mode=args.name_mode,
                    entity_type=entity_type,
                    top_n=max(1, args.review_top_n),
                    reason_code_filter={
                        "manual_validation_required",
                        "needs_new_term",
                        "suggest_broad_term_mapping",
                        "metabolite_mode_blocked_protein_like_query",
                    },
                )
                print(f"[OK] wrote {review_rows} withheld review rows to {review_output_path}")
            if args.review_queue_output:
                review_queue_rows = write_review_tsv(
                    rows,
                    Path(args.review_queue_output),
                    index=index,
                    min_score=args.min_score,
                    matrix_priority=matrix_priority,
                    measurement_context=measurement_context,
                    additional_contexts=additional_contexts,
                    additional_context_keywords=additional_context_keywords,
                    name_mode=args.name_mode,
                    entity_type=entity_type,
                    top_n=max(1, args.review_top_n),
                )
                print(f"[OK] wrote {review_queue_rows} review-queue rows to {args.review_queue_output}")
            if args.withheld_triage_output:
                triage_output_path = Path(args.withheld_triage_output)
                triage_input_path = review_output_path
                temp_review_path: Path | None = None
                try:
                    if triage_input_path is None:
                        triage_output_path.parent.mkdir(parents=True, exist_ok=True)
                        with tempfile.NamedTemporaryFile(
                            mode="w",
                            suffix=".tsv",
                            prefix="withheld_review_",
                            dir=str(triage_output_path.parent),
                            delete=False,
                        ) as tmp_handle:
                            temp_review_path = Path(tmp_handle.name)
                        temp_review_rows = write_review_tsv(
                            rows,
                            temp_review_path,
                            index=index,
                            min_score=args.min_score,
                            matrix_priority=matrix_priority,
                            measurement_context=measurement_context,
                            additional_contexts=additional_contexts,
                            additional_context_keywords=additional_context_keywords,
                            name_mode=args.name_mode,
                            entity_type=entity_type,
                            top_n=max(1, args.review_top_n),
                            reason_code_filter={"manual_validation_required"},
                        )
                        triage_input_path = temp_review_path
                        print(
                            f"[OK] generated {temp_review_rows} temporary withheld review rows "
                            f"for triage"
                        )

                    triage_rows, triage_counts = write_withheld_triage_tsv(
                        triage_input_path,
                        triage_output_path,
                        index=index,
                    )
                    triage_summary = ", ".join(
                        f"{key}={triage_counts.get(key, 0)}"
                        for key in ["reject_high", "review_needed", "accept_medium", "accept_high"]
                    )
                    print(
                        f"[OK] wrote {triage_rows} withheld triage rows to {args.withheld_triage_output} "
                        f"({triage_summary})"
                    )
                finally:
                    if temp_review_path is not None:
                        try:
                            temp_review_path.unlink(missing_ok=True)
                        except OSError:
                            pass
            if args.cache_writeback:
                added = write_back_cache(rows, Path(args.analyte_cache), args.cache_min_confidence)
                print(f"[OK] cache write-back added {added} entries to {args.analyte_cache}")
            print(f"[OK] wrote {len(rows)} rows to {args.output}")
            return 0

        if args.command == "trait-cache-refresh":
            trait_cache_path = Path(args.trait_cache)
            studies_tsv_path = Path(args.studies_tsv)
            output_cache_arg = normalize(args.output_cache)
            output_cache_path = Path(args.output_cache) if output_cache_arg else trait_cache_path
            state_output_arg = normalize(args.state_output)
            state_output_path = Path(args.state_output) if state_output_arg else None
            studies_sha256 = file_sha256(studies_tsv_path)

            if state_output_path is not None and state_output_path.exists() and not bool(args.force_refresh):
                try:
                    state = json.loads(state_output_path.read_text(encoding="utf-8"))
                except Exception:  # noqa: BLE001 - malformed state should not block refresh.
                    state = {}
                if (
                    normalize(state.get("studies_tsv_sha256", "")) == studies_sha256
                    and output_cache_path.exists()
                ):
                    print(
                        "[OK] trait cache refresh skipped; studies TSV hash unchanged "
                        f"({studies_sha256[:12]})"
                    )
                    print(f"[OK] using existing trait cache at {output_cache_path}")
                    return 0

            efo_obo_path, efo_obo_source = ensure_efo_obo_available(
                efo_obo_path=Path(args.efo_obo),
                bundled_url=args.efo_obo_bundled_url,
                timeout=float(args.efo_obo_download_timeout),
            )
            if efo_obo_source != "local":
                print(
                    f"[OK] provisioned EFO OBO at {efo_obo_path} "
                    f"(source={efo_obo_source})"
                )
            ontology_index = load_trait_ontology_index(efo_obo_path)

            refresh_result = build_catalog_trait_refresh(
                studies_tsv=studies_tsv_path,
                trait_cache_path=trait_cache_path,
                ontology_index=ontology_index,
                publication_tag=normalize(args.publication_tag) or "CatalogMappedTraitsImport",
                studies_sha256=studies_sha256,
            )

            def write_rows(path_arg: str, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
                if not normalize(path_arg):
                    return
                path = Path(path_arg)
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
                    writer.writeheader()
                    writer.writerows(rows)

            qc_fieldnames = [
                "DISEASE/TRAIT",
                "MAPPED_TRAIT",
                "MAPPED_TRAIT_URI",
                "resolved_mapped_ids",
                "resolved_mapped_labels",
                "resolved_id_count",
                "unresolved_ids",
                "obsolete_ids",
                "id_resolution_modes",
                "label_match_ratio",
                "label_match_modes",
                "confidence",
                "qc_status",
            ]

            write_rows(args.qc_output, refresh_result.get("qc_rows", []), qc_fieldnames)
            write_rows(
                args.high_confidence_output,
                refresh_result.get("high_conf_rows", []),
                qc_fieldnames,
            )
            write_rows(
                args.to_add_output,
                refresh_result.get("rows_to_add", []),
                TRAIT_CACHE_COLUMNS,
            )

            output_cache_path.parent.mkdir(parents=True, exist_ok=True)
            with output_cache_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=TRAIT_CACHE_COLUMNS,
                    delimiter="\t",
                    extrasaction="ignore",
                )
                writer.writeheader()
                writer.writerows(refresh_result.get("merged_rows", []))

            if normalize(args.summary_output):
                summary_path = Path(args.summary_output)
                summary_path.parent.mkdir(parents=True, exist_ok=True)
                summary_path.write_text(
                    json.dumps(refresh_result.get("summary", {}), indent=2, ensure_ascii=True) + "\n",
                    encoding="utf-8",
                )
                print(f"[OK] wrote trait-cache refresh summary to {summary_path}")

            if state_output_path is not None:
                state_output_path.parent.mkdir(parents=True, exist_ok=True)
                state_payload = {
                    "updated_at_utc": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
                    "studies_tsv": str(studies_tsv_path),
                    "studies_tsv_sha256": studies_sha256,
                    "trait_cache_output": str(output_cache_path),
                    "cache_rows_added": refresh_result.get("summary", {}).get("cache_rows_added", 0),
                }
                state_output_path.write_text(
                    json.dumps(state_payload, indent=2, ensure_ascii=True) + "\n",
                    encoding="utf-8",
                )
                print(f"[OK] wrote trait-cache refresh state to {state_output_path}")

            summary = refresh_result.get("summary", {})
            print(
                "[OK] refreshed trait cache from studies TSV "
                f"(input_rows={summary.get('input_rows', 0)}, "
                f"high_conf_initial={summary.get('high_conf_initial', 0)}, "
                f"after_unambiguous={summary.get('after_unambiguous_trait_filter', 0)}, "
                f"after_specific={summary.get('after_specific_filter_candidate_rows', 0)}, "
                f"added={summary.get('cache_rows_added', 0)}, "
                f"conflicts_skipped={summary.get('conflicting_rows_skipped', 0)})"
            )
            print(f"[OK] wrote updated trait cache to {output_cache_path}")
            return 0

        if args.command == "trait-map":
            query_inputs = load_trait_query_inputs(
                Path(args.input),
                default_trait_scale=args.default_trait_scale,
            )
            if not query_inputs:
                raise ValueError("No usable trait queries found in input")
            input_type_counts = Counter(normalize_trait_input_type(row.get("input_type", "auto")) for row in query_inputs)
            trait_scale_counts = Counter(normalize_trait_scale(row.get("trait_scale", "auto")) for row in query_inputs)
            print(
                "[OK] trait input rows loaded "
                f"(total={len(query_inputs)}, "
                f"auto={input_type_counts.get('auto', 0)}, "
                f"trait_text={input_type_counts.get('trait_text', 0)}, "
                f"icd10={input_type_counts.get('icd10', 0)}, "
                f"phecode={input_type_counts.get('phecode', 0)}, "
                f"ukb_field={input_type_counts.get('ukb_field', 0)}, "
                f"scale_auto={trait_scale_counts.get('auto', 0)}, "
                f"scale_binary={trait_scale_counts.get('binary', 0)}, "
                f"scale_quantitative={trait_scale_counts.get('quantitative', 0)})"
            )
            icd10_query_rows = sum(
                1
                for row in query_inputs
                if is_strict_icd10_code(
                    normalize_icd10_code(row.get("icd10", "") or row.get("query", ""))
                )
            )
            ukb_query_rows = sum(1 for row in query_inputs if extract_ukb_field_ids(row.get("query", "")))
            ukb_category_query_rows = sum(
                1 for row in query_inputs if extract_ukb_category_ids(row.get("query", ""))
            )

            trait_cache_path = Path(args.trait_cache)
            trait_cache_path, seeded = seed_default_trait_cache_if_needed(trait_cache_path)
            if seeded:
                print(
                    f"[OK] seeded bundled trait cache at {trait_cache_path} "
                    f"from legacy file {LEGACY_TRAIT_CACHE}"
                )
            efo_obo_path = Path(args.efo_obo)
            efo_obo_path, efo_obo_source = ensure_efo_obo_available(
                efo_obo_path=efo_obo_path,
                bundled_url=args.efo_obo_bundled_url,
                timeout=float(args.efo_obo_download_timeout),
            )
            if efo_obo_source != "local":
                print(
                    f"[OK] provisioned EFO OBO at {efo_obo_path} "
                    f"(source={efo_obo_source})"
                )
            ontology_index = load_trait_ontology_index(efo_obo_path)
            ukb_catalog_arg = normalize(args.ukb_field_catalog)
            ukb_catalog_path = Path(args.ukb_field_catalog) if ukb_catalog_arg else None
            ukb_field_titles = load_ukb_field_title_index(ukb_catalog_path)
            if ukb_catalog_path is not None and ukb_catalog_path.exists():
                print(
                    f"[OK] loaded UKB field catalog with {len(ukb_field_titles)} titles "
                    f"from {ukb_catalog_path}"
                )
            elif ukb_query_rows > 0:
                print(
                    "[WARN] input contains UKB data-field queries but UKB field catalog was not found; "
                    f"checked path: {ukb_catalog_path}"
                )

            ukb_roots: list[Path] = []
            if ukb_catalog_path is not None:
                ukb_roots.append(ukb_catalog_path.parent)
            default_ukb_root = DEFAULT_UKB_FIELD_CATALOG.parent
            if default_ukb_root not in ukb_roots:
                ukb_roots.append(default_ukb_root)

            def pick_ukb_metadata_path(filename: str) -> Path | None:
                for root in ukb_roots:
                    candidate = root / filename
                    if candidate.exists():
                        return candidate
                if ukb_roots:
                    return ukb_roots[0] / filename
                return None

            ukb_field_metadata_path = pick_ukb_metadata_path("field.txt")
            ukb_category_catalog_path = pick_ukb_metadata_path("category.txt")
            ukb_category_tree_path = pick_ukb_metadata_path("catbrowse.txt")
            ukb_field_supplement_cache_path = pick_ukb_metadata_path("field_trait_supplement_cache.tsv")
            icd10_supplement_cache_arg = normalize(args.icd10_supplement_cache)
            icd10_supplement_cache_path = (
                Path(args.icd10_supplement_cache)
                if icd10_supplement_cache_arg
                else DEFAULT_ICD10_SUPPLEMENT_CACHE
            )
            icd10_label_cache_arg = normalize(args.icd10_label_cache)
            icd10_label_cache_path = (
                Path(args.icd10_label_cache)
                if icd10_label_cache_arg
                else DEFAULT_ICD10_LABEL_CACHE
            )

            extra_trait_cache_paths: list[Path] = []
            if ukb_field_supplement_cache_path is not None and ukb_field_supplement_cache_path.exists():
                extra_trait_cache_paths.append(ukb_field_supplement_cache_path)
                print(
                    f"[OK] loaded UKB supplemental field trait cache from {ukb_field_supplement_cache_path}"
                )
            elif ukb_query_rows > 0:
                print(
                    "[WARN] UKB supplemental field trait cache not found; "
                    "setup can build it with setup-bundled-caches."
                )
            if icd10_supplement_cache_path.exists():
                extra_trait_cache_paths.append(icd10_supplement_cache_path)
                print(
                    f"[OK] loaded ICD10 supplemental trait cache from {icd10_supplement_cache_path}"
                )
            elif icd10_query_rows > 0:
                print(
                    "[WARN] ICD10 supplemental trait cache not found; "
                    "setup can build it with setup-bundled-caches."
                )
            trait_cache_index = load_trait_cache_index(trait_cache_path, extra_paths=extra_trait_cache_paths)
            official_icd10_label_index: dict[str, tuple[str, ...]] = {}
            if icd10_label_cache_path.exists():
                official_icd10_label_index = load_icd10_label_cache(icd10_label_cache_path)
                print(
                    f"[OK] loaded ICD10 label cache with {len(official_icd10_label_index)} codes "
                    f"from {icd10_label_cache_path}"
                )
            elif icd10_query_rows > 0:
                print(
                    "[WARN] ICD10 label cache not found; bare ICD10 code outputs may have blank "
                    "input_icd10_label. Build it with icd10-label-cache-build."
                )
            if official_icd10_label_index:
                trait_cache_index = dict(trait_cache_index)
                trait_cache_index["official_icd10_label_index"] = official_icd10_label_index

            ukb_field_to_category, ukb_category_to_fields = load_ukb_field_category_index(ukb_field_metadata_path)
            ukb_category_titles = load_ukb_category_title_index(ukb_category_catalog_path)
            ukb_category_parents = load_ukb_category_parent_index(ukb_category_tree_path)
            metadata_field_titles = load_ukb_field_title_index(ukb_field_metadata_path)
            merged_field_titles = 0
            for field_id, title in metadata_field_titles.items():
                if field_id in ukb_field_titles:
                    continue
                ukb_field_titles[field_id] = title
                merged_field_titles += 1

            if ukb_field_metadata_path is not None and ukb_field_metadata_path.exists():
                print(
                    f"[OK] loaded UKB field metadata with {len(ukb_category_to_fields)} category groups "
                    f"from {ukb_field_metadata_path}"
                )
                if merged_field_titles > 0:
                    print(
                        f"[OK] augmented UKB field title index with {merged_field_titles} extra titles "
                        f"from {ukb_field_metadata_path}"
                    )
            elif ukb_category_query_rows > 0:
                print(
                    "[WARN] input contains UKB category queries but UKB field metadata was not found; "
                    f"checked path: {ukb_field_metadata_path}"
                )

            if ukb_field_to_category:
                fields_with_titles = sum(1 for field_id in ukb_field_to_category if field_id in ukb_field_titles)
                print(
                    "[OK] UKB field title coverage: "
                    f"{fields_with_titles}/{len(ukb_field_to_category)} fields"
                )

            if ukb_category_titles:
                print(
                    f"[OK] loaded UKB category catalog with {len(ukb_category_titles)} titles "
                    f"from {ukb_category_catalog_path}"
                )
            elif ukb_category_query_rows > 0:
                print(
                    "[WARN] input contains UKB category queries but UKB category catalog was not found; "
                    f"checked path: {ukb_category_catalog_path}"
                )
                print(
                    "[WARN] download with: "
                    "curl -sS -L 'https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=3' "
                    f"-o {DEFAULT_UKB_CATEGORY_CATALOG}"
                )

            if ukb_category_parents:
                print(
                    f"[OK] loaded UKB category tree with {len(ukb_category_parents)} child links "
                    f"from {ukb_category_tree_path}"
                )
            elif ukb_category_query_rows > 0:
                print(
                    "[WARN] input contains UKB category queries but UKB category tree was not found; "
                    f"checked path: {ukb_category_tree_path}"
                )
                print(
                    "[WARN] download with: "
                    "curl -sS -L 'https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=13' "
                    f"-o {DEFAULT_UKB_CATEGORY_TREE}"
                )
            cache_resolution, cache_qc_stats, cache_qc_issues = build_trait_cache_term_resolution(
                trait_cache_index,
                ontology_index,
            )
            rows_with_obsolete_ids = cache_qc_stats.get("rows_with_obsolete_ids", 0)
            rows_with_unresolved_ids = cache_qc_stats.get("rows_with_unresolved_ids", 0)
            rows_without_active_terms = cache_qc_stats.get("rows_without_active_terms", 0)

            if rows_with_unresolved_ids > 0 or rows_without_active_terms > 0:
                print(
                    "[WARN] trait cache ontology QC: "
                    f"rows_with_obsolete_ids={rows_with_obsolete_ids}, "
                    f"rows_with_unresolved_ids={rows_with_unresolved_ids}, "
                    f"rows_without_active_terms={rows_without_active_terms}, "
                    f"ids_obsolete_replaced_by={cache_qc_stats.get('ids_obsolete_replaced_by', 0)}"
                )
                for issue in cache_qc_issues:
                    print(f"[WARN] {issue}")
            elif rows_with_obsolete_ids > 0:
                print(
                    "[OK] trait cache ontology QC: "
                    f"rows_checked={cache_qc_stats.get('records_total', 0)}; "
                    f"obsolete_rows_auto_remapped={rows_with_obsolete_ids}; "
                    f"ids_obsolete_replaced_by={cache_qc_stats.get('ids_obsolete_replaced_by', 0)}; "
                    "unresolved=0"
                )
            else:
                print(
                    "[OK] trait cache ontology QC: "
                    f"rows_checked={cache_qc_stats.get('records_total', 0)}; "
                    "no obsolete/unresolved ontology IDs in cache rows"
                )

            output_path = Path(args.output)
            review_output_path = Path(args.review_output) if args.review_output else None
            qc_risk_output_arg = normalize(args.qc_risk_output)
            qc_risk_output_path = (
                Path(args.qc_risk_output)
                if qc_risk_output_arg
                else output_path.with_name(f"{output_path.stem}{DEFAULT_TRAIT_QC_RISK_SUFFIX}")
            )
            flush_every = max(1, int(args.flush_every))
            memoize_queries = bool(args.memoize_queries)
            query_cache_max_entries = max(0, int(args.query_cache_max_entries))

            if bool(args.stream_output):
                with TraitMapStreamingWriter(
                    output_path=output_path,
                    review_path=review_output_path,
                    flush_every=flush_every,
                ) as streaming_writer:
                    map_trait_queries(
                        query_inputs,
                        trait_cache_index=trait_cache_index,
                        ontology_index=ontology_index,
                        trait_cache_resolution=cache_resolution,
                        trait_cache_path=trait_cache_path,
                        efo_obo_path=efo_obo_path,
                        ukb_field_titles=ukb_field_titles,
                        ukb_category_titles=ukb_category_titles,
                        ukb_category_parents=ukb_category_parents,
                        ukb_category_to_fields=ukb_category_to_fields,
                        ukb_field_to_category=ukb_field_to_category,
                        top_k=max(1, args.top_k),
                        min_score=max(0.0, min(1.0, args.min_score)),
                        force_map_best=bool(args.force_map_best),
                        show_progress=bool(args.progress),
                        row_callback=streaming_writer.write_row,
                        collect_rows=False,
                        memoize_queries=memoize_queries,
                        query_cache_max_entries=query_cache_max_entries,
                    )
                    total_rows = streaming_writer.rows_written
                    review_count = streaming_writer.review_rows_written
            else:
                rows = map_trait_queries(
                    query_inputs,
                    trait_cache_index=trait_cache_index,
                    ontology_index=ontology_index,
                    trait_cache_resolution=cache_resolution,
                    trait_cache_path=trait_cache_path,
                    efo_obo_path=efo_obo_path,
                    ukb_field_titles=ukb_field_titles,
                    ukb_category_titles=ukb_category_titles,
                    ukb_category_parents=ukb_category_parents,
                    ukb_category_to_fields=ukb_category_to_fields,
                    ukb_field_to_category=ukb_field_to_category,
                    top_k=max(1, args.top_k),
                    min_score=max(0.0, min(1.0, args.min_score)),
                    force_map_best=bool(args.force_map_best),
                    show_progress=bool(args.progress),
                    memoize_queries=memoize_queries,
                    query_cache_max_entries=query_cache_max_entries,
                )
                write_trait_tsv(rows, output_path)
                total_rows = len(rows)
                review_count = 0
                if review_output_path is not None:
                    review_count = write_trait_review_tsv(rows, review_output_path)

            if review_output_path is not None:
                print(f"[OK] wrote {review_count} review rows to {review_output_path}")
            risk_count = write_trait_qc_risk_tsv(output_path, qc_risk_output_path)
            print(f"[OK] wrote {risk_count} QC-risk rows to {qc_risk_output_path}")
            print(f"[OK] wrote {total_rows} rows to {output_path}")
            return 0

        raise ValueError(f"Unsupported command: {args.command}")
    except KeyboardInterrupt:
        print(
            "[WARN] Interrupted by user. Partial output files may already contain processed rows.",
            file=sys.stderr,
        )
        return 130
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
