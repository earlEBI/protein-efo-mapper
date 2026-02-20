#!/usr/bin/env python3
"""Bulk/offline mapper for protein/metabolite queries -> measurement EFO IDs.

Architecture:
- index-build: build a local JSON index from cache/reference files
- map: map large query files using local index only (fast, scalable)
"""

from __future__ import annotations

import argparse
import csv
import difflib
import json
import re
import subprocess
import sys
import urllib.request
import time
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
EFO_RE = re.compile(r"^(EFO|OBA)[:_]\d+$", re.IGNORECASE)
ALLOWED_ONTOLOGY_ID_RE = re.compile(r"^(EFO|OBA)[:_]\d+$", re.IGNORECASE)
MULTI_ID_SEP_RE = re.compile(r"[;|]")
DEFAULT_ANALYTE_CACHE = Path(__file__).resolve().parents[1] / "references" / "analyte_to_efo_cache.tsv"
DEFAULT_TERM_CACHE = Path(__file__).resolve().parents[1] / "references" / "efo_measurement_terms_cache.tsv"
DEFAULT_INDEX = Path(__file__).resolve().parents[1] / "references" / "measurement_index.json"
DEFAULT_UNIPROT_ALIASES = Path(__file__).resolve().parents[1] / "references" / "uniprot_aliases.tsv"
DEFAULT_EFO_OBO_URL = "https://github.com/EBISPOT/efo/releases/latest/download/efo.obo"
DEFAULT_EFO_OBO_LOCAL = Path(__file__).resolve().parents[1] / "references" / "efo.obo"
DEFAULT_MEASUREMENT_ROOT = "EFO:0001444"
DEFAULT_MATRIX_PRIORITY = ["plasma", "blood", "serum"]
DEFAULT_MEASUREMENT_CONTEXT = "blood"
# Token-retrieved candidates are noisier than exact/synonym matches. Keep them
# out of auto-validated output unless confidence remains reasonably high.
TOKEN_SUBJECT_EXACT_AUTO_VALIDATE_MIN_SCORE = 0.70
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


@dataclass
class Candidate:
    efo_id: str
    label: str
    score: float
    matched_via: str
    evidence: str
    is_validated: bool


@dataclass(frozen=True)
class AccessionProfile:
    alias_terms: tuple[str, ...]
    expected_symbol_nums: frozenset[str]
    subject_terms: frozenset[str]
    symbol_terms: tuple[str, ...]


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
def normalize(text: str) -> str:
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


def effective_name_mode(default_name_mode: str, input_type: str) -> str:
    kind = normalize_input_type(input_type)
    if kind == "protein_name":
        # Protein-name inputs are often not uniquely resolvable via strict alias
        # lookup; enable fuzzy candidate retrieval but rely on validation gates.
        return "fuzzy"
    if kind == "gene_symbol":
        return "strict"
    if kind == "gene_id":
        return "strict"
    if kind == "accession":
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


@lru_cache(maxsize=100_000)
def query_variants(query: str) -> tuple[str, ...]:
    q = normalize(query)
    variants: set[str] = {q}
    parts = [normalize(p) for p in MULTI_ID_SEP_RE.split(q) if normalize(p)]
    variants.update(parts)
    for part in list(variants):
        canonical = canonical_accession(part)
        if canonical:
            variants.add(canonical)
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
            rows.append({"efo_id": efo_id, "label": label, "synonyms": syns})
    return rows


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
            rows.append(
                {
                    "analyte_key": analyte,
                    "efo_id": efo_id,
                    "label": normalize(row.get("label") or ""),
                    "confidence": confidence,
                    "validation": norm_key(row.get("validation") or ""),
                    "source": normalize(row.get("source") or "cache"),
                    "evidence": normalize(row.get("evidence") or "cached mapping"),
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
            source = normalize(row.get("source") or "")
            # Merge duplicate rows by accession if present.
            if acc in records:
                prior_aliases = [normalize(x) for x in records[acc].get("aliases", "").split("|") if normalize(x)]
                merged = sorted(set(prior_aliases) | set(aliases))
                records[acc]["aliases"] = "|".join(merged)
                if not records[acc].get("gene") and gene:
                    records[acc]["gene"] = gene
                continue
            records[acc] = {
                "accession": acc,
                "aliases": "|".join(sorted(set(aliases))),
                "gene": gene,
                "source": source,
            }
    return records


def write_uniprot_alias_records(path: Path, records: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["accession", "aliases", "gene", "source"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for acc in sorted(records):
            row = records[acc]
            aliases = [normalize(x) for x in (row.get("aliases") or "").split("|") if normalize(x)]
            writer.writerow(
                {
                    "accession": acc,
                    "aliases": "|".join(sorted(set(aliases))),
                    "gene": normalize(row.get("gene") or ""),
                    "source": normalize(row.get("source") or ""),
                }
            )


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


def parse_uniprot_entry_aliases(entry: dict[str, Any]) -> tuple[list[str], str]:
    aliases: set[str] = set()
    gene_primary = ""

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

    return sorted(aliases), gene_primary


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

    def fetch_one(acc: str) -> tuple[str, list[str], str, bool]:
        try:
            entry = fetch_json(UNIPROT_ENTRY_URL.format(acc=acc), timeout=timeout)
            aliases, gene = parse_uniprot_entry_aliases(entry)
            if not aliases and not gene:
                return acc, [], "", False
            return acc, aliases, gene, True
        except Exception:
            return acc, [], "", False

    added = 0
    failed = 0
    if workers <= 1:
        results = [fetch_one(acc) for acc in missing]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(fetch_one, missing))

    for acc, aliases, gene, ok in results:
        if not ok:
            failed += 1
            continue
        records[acc] = {
            "accession": acc,
            "aliases": "|".join(sorted(set(aliases))),
            "gene": normalize(gene),
            "source": source_tag,
        }
        added += 1

    if added > 0:
        write_uniprot_alias_records(alias_path, records)
    return added, failed


def build_index(
    term_rows: list[dict[str, Any]],
    analyte_rows: list[dict[str, Any]],
    accession_aliases: dict[str, list[str]],
) -> dict[str, Any]:
    term_meta: dict[str, dict[str, Any]] = {}
    exact_term_index: dict[str, set[str]] = {}
    token_index: dict[str, set[str]] = {}

    for row in term_rows:
        efo_id = row["efo_id"]
        label = row["label"]
        syns = row["synonyms"]

        term_meta[efo_id] = {
            "label": label,
            "synonyms": syns,
            "subject": normalize_measurement_subject(label),
        }

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

    analyte_index: dict[str, list[dict[str, Any]]] = {}
    for row in analyte_rows:
        analyte_index.setdefault(row["analyte_key"], []).append(row)

    alias_accession_index: dict[str, set[str]] = {}
    alias_loose_index: dict[str, set[str]] = {}
    symbol_accession_index: dict[str, set[str]] = {}
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

    return {
        "version": "2",
        "built_at": datetime.now(UTC).isoformat(),
        "term_meta": term_meta,
        "exact_term_index": {k: sorted(v) for k, v in exact_term_index.items()},
        "token_index": {k: sorted(v) for k, v in token_index.items()},
        "analyte_index": analyte_index,
        "accession_alias_index": accession_aliases,
        "alias_accession_index": {k: sorted(v) for k, v in alias_accession_index.items()},
        "alias_loose_index": {k: sorted(v) for k, v in alias_loose_index.items()},
        "symbol_accession_index": {k: sorted(v) for k, v in symbol_accession_index.items()},
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

    if "alias_accession_index" not in index:
        alias_accession_index: dict[str, set[str]] = {}
        alias_loose_index: dict[str, set[str]] = {}
        symbol_accession_index: dict[str, set[str]] = {}
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
        index["alias_accession_index"] = {k: sorted(v) for k, v in alias_accession_index.items()}
        index["alias_loose_index"] = {k: sorted(v) for k, v in alias_loose_index.items()}
        index["symbol_accession_index"] = {k: sorted(v) for k, v in symbol_accession_index.items()}
    elif "symbol_accession_index" not in index or "alias_loose_index" not in index:
        alias_loose_index: dict[str, set[str]] = {}
        symbol_accession_index: dict[str, set[str]] = {}
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
        index["alias_loose_index"] = {k: sorted(v) for k, v in alias_loose_index.items()}
        index["symbol_accession_index"] = {k: sorted(v) for k, v in symbol_accession_index.items()}
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
    )


def merge_candidate(store: dict[str, Candidate], cand: Candidate) -> None:
    current = store.get(cand.efo_id)
    if current is None:
        store[cand.efo_id] = cand
        return
    if cand.is_validated and not current.is_validated:
        store[cand.efo_id] = cand
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
    return re.search(r"\b[a-z0-9][a-z0-9\-]*\s+/\s+[a-z0-9][a-z0-9\-]*\b", q) is not None


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


def resolve_query_accessions(query: str, index: dict[str, Any]) -> tuple[list[str], bool, bool]:
    # Returns: (resolved accessions, used_alias_resolution, ambiguous_alias_resolution)
    accession_index = index.get("alias_accession_index", {})
    symbol_index = index.get("symbol_accession_index", {})
    loose_index = index.get("alias_loose_index", {})
    accession_scores: dict[str, int] = {}
    used_alias_resolution = False
    ambiguous_alias_resolution = False
    variants = name_lookup_variants(query)

    for variant in variants:
        acc = canonical_accession(variant)
        if acc:
            return [acc], False, False

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
) -> list[Candidate]:
    term_meta: dict[str, dict[str, Any]] = index.get("term_meta", {})
    exact_term_index: dict[str, list[str]] = index.get("exact_term_index", {})
    token_index: dict[str, list[str]] = index.get("token_index", {})
    analyte_index: dict[str, list[dict[str, Any]]] = index.get("analyte_index", {})

    candidates: dict[str, Candidate] = {}
    resolved_accessions, used_alias_resolution, _ambiguous_alias_resolution = resolve_query_accessions(query, index)
    is_accession_input = any(canonical_accession(v) for v in query_variants(query))
    is_symbol_input = bool(symbol_query_key(query) or uniprot_mnemonic_symbol_keys(query, index))
    accession_profiles = compile_accession_profiles(accession_profiles_for_accessions(resolved_accessions, index))
    has_accession_profiles = bool(accession_profiles)
    allow_ratio_terms = query_requests_ratio(query)
    strict_unresolved_name = (
        name_mode == "strict"
        and not is_accession_input
        and not EFO_RE.match(query)
        and not used_alias_resolution
    )
    strict_tokens = (
        strict_query_tokens(query) if (name_mode == "strict" and not is_accession_input and not is_symbol_input) else frozenset()
    )
    term_cache: dict[str, tuple[str, list[str], str, frozenset[str], frozenset[str]]] = {}
    identity_cache: dict[str, bool] = {}
    context_cache: dict[tuple[str, str], bool] = {}
    token_gate_cache: dict[str, bool] = {}
    candidate_eval_cache: dict[str, tuple[bool, str, list[str], str, frozenset[str], frozenset[str], int, int, bool]] = {}

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
        if not has_accession_profiles:
            return True
        cached = identity_cache.get(efo_id)
        if cached is not None:
            return cached
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
        identity_cache[efo_id] = ok
        return ok

    def profile_exact_bonus(subject_terms: frozenset[str]) -> int:
        if not has_accession_profiles:
            return 0
        return profile_subject_exact_match_bonus(accession_profiles, subject_terms)

    def profile_exact_rank(subject_terms: frozenset[str]) -> int:
        if not has_accession_profiles:
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

        exact_bonus = profile_exact_bonus(subject_terms) if ok else 0
        exact_rank = profile_exact_rank(subject_terms) if ok else 0
        measurement_like = is_measurement_like(label, syns) if ok else False
        result = (ok, label, syns, merged, nums, subject_terms, exact_bonus, exact_rank, measurement_like)
        candidate_eval_cache[efo_id] = result
        return result

    if EFO_RE.match(query):
        efo_id = normalize_id(query)
        meta = term_meta.get(efo_id, {})
        label, syns, _, _, _ = get_term_bundle(efo_id)
        if label and not context_ok(efo_id, label, syns):
            return []
        direct = [
            Candidate(
                efo_id=efo_id,
                label=label,
                score=1.0,
                matched_via="direct-id",
                evidence="direct input EFO/OBA id",
                is_validated=bool(meta),
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
            cached_label = normalize(cached.get("label") or "")
            meta_label, syns, merged, nums, subject_terms = get_term_bundle(efo_id)
            label = cached_label or meta_label
            if not allow_ratio_terms and is_ratio_trait(label, syns):
                continue
            if not token_gate_ok(efo_id, label, syns):
                continue
            if not context_ok(efo_id, label, syns):
                continue
            if not identity_ok(efo_id, label, syns, merged, nums, subject_terms):
                continue
            validated = cached.get("validation") == "validated" or bool(meta)
            merge_candidate(
                candidates,
                Candidate(
                    efo_id=efo_id,
                    label=label,
                    score=float(cached.get("confidence") or 0.8),
                    matched_via="local-analyte-cache",
                    evidence=normalize(cached.get("evidence") or "cached mapping"),
                    is_validated=validated,
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
            ok, label, syns, _merged, _nums, _subject_terms, exact_bonus, exact_rank, _measurement_like = candidate_eval(efo_id)
            if not ok:
                continue
            score = max(0.9, similarity_score_prepared(term_lower, term_toks, label, syns))
            score = min(1.0, max(0.0, score + matrix_bonus(label, syns, matrix_priority)))
            score = min(1.0, score + 0.02 * exact_bonus)
            merge_candidate(
                candidates,
                Candidate(
                    efo_id=efo_id,
                    label=label,
                    score=score,
                    matched_via="local-term-exact" + ("-subject-exact" if exact_rank else ""),
                    evidence=f"exact normalized term match: {term}",
                    is_validated=True,
                ),
            )
            seen_ids.add(efo_id)

        # Token-based retrieval + lexical rerank
        candidate_ids: set[str] = set()
        for tok in filtered_tokens:
            candidate_ids.update(token_index.get(tok, []))

        for efo_id in candidate_ids:
            if efo_id in seen_ids:
                continue
            ok, label, syns, _merged, _nums, _subject_terms, exact_bonus, exact_rank, measurement_like = candidate_eval(efo_id)
            if not ok:
                continue
            score = similarity_score_prepared(term_lower, term_toks, label, syns)
            if not measurement_like:
                score -= 0.35
            score += matrix_bonus(label, syns, matrix_priority)
            score += 0.02 * exact_bonus
            score = min(1.0, max(0.0, score))
            if score < 0.30:
                continue
            merge_candidate(
                candidates,
                Candidate(
                    efo_id=efo_id,
                    label=label,
                    score=score,
                    matched_via="local-term-token" + ("-subject-exact" if exact_rank else ""),
                    evidence=f"token retrieval + lexical score using: {term}",
                    is_validated=True,
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
) -> list[dict[str, str]]:
    query = normalize(query)
    input_type_norm = normalize_input_type(input_type)
    candidates = candidates_from_index(
        query,
        index,
        matrix_priority=matrix_priority,
        measurement_context=measurement_context,
        additional_contexts=additional_contexts,
        additional_context_keywords=additional_context_keywords,
        name_mode=name_mode,
    )
    eligible = [c for c in candidates if c.score >= min_score]

    def auto_validates(cand: Candidate) -> bool:
        if not cand.is_validated:
            return False
        if cand.matched_via.startswith("local-term-token"):
            if "subject-exact" not in cand.matched_via:
                return False
            if cand.score < TOKEN_SUBJECT_EXACT_AUTO_VALIDATE_MIN_SCORE:
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
) -> list[Candidate]:
    term_meta: dict[str, dict[str, Any]] = index.get("term_meta", {})
    exact_term_index: dict[str, list[str]] = index.get("exact_term_index", {})
    token_index: dict[str, list[str]] = index.get("token_index", {})
    candidates: dict[str, Candidate] = {}
    resolved_accessions, _used_alias_resolution, _ambiguous_alias_resolution = resolve_query_accessions(query, index)
    is_accession_input = any(canonical_accession(v) for v in query_variants(query))
    allow_ratio_terms = query_requests_ratio(query)
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
) -> tuple[str, str, list[str], list[Candidate]]:
    query = normalize(query)
    effective_mode = effective_name_mode(name_mode, input_type)
    strict_candidates = candidates_from_index(
        query,
        index,
        matrix_priority=matrix_priority,
        measurement_context=measurement_context,
        additional_contexts=additional_contexts,
        additional_context_keywords=additional_context_keywords,
        name_mode=effective_mode,
    )
    resolved_accessions, used_alias_resolution, ambiguous_alias_resolution = resolve_query_accessions(query, index)
    is_accession_input = any(canonical_accession(v) for v in query_variants(query))

    if strict_candidates:
        best = strict_candidates[0]
        if best.score >= min_score:
            return (
                "manual_validation_required",
                (
                    f"best candidate {format_ontology_id_for_output(best.efo_id)} scored {best.score:.3f} but did not pass "
                    "auto-validation gate"
                ),
                resolved_accessions,
                strict_candidates,
            )
        return (
            "below_score_threshold",
            (
                f"best candidate {format_ontology_id_for_output(best.efo_id)} "
                f"scored {best.score:.3f}, below min_score={min_score:.3f}"
            ),
            resolved_accessions,
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
        )
        if any_context_candidates:
            best = any_context_candidates[0]
            return (
                "context_mismatch",
                (
                    f"candidate {format_ontology_id_for_output(best.efo_id)} exists only outside "
                    f"requested context '{measurement_context}'"
                ),
                resolved_accessions,
                any_context_candidates,
            )

    if effective_mode == "strict" and not is_accession_input and not EFO_RE.match(query) and not used_alias_resolution:
        return (
            "unresolved_name_strict_mode",
            "query did not resolve to a unique accession in strict mode",
            resolved_accessions,
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
    )
    if resolved_accessions and probe_candidates:
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

    if ambiguous_alias_resolution and not resolved_accessions:
        return (
            "ambiguous_alias_resolution",
            "query alias matched multiple accessions without a stable winner",
            resolved_accessions,
            probe_candidates,
        )

    if resolved_accessions:
        return (
            "no_measurement_term_for_resolved_accession",
            "resolved accession(s) found but no valid measurement term passed filters",
            resolved_accessions,
            probe_candidates,
        )

    return (
        "no_candidate_found",
        "no local candidates were retrieved from exact/token indexes",
        resolved_accessions,
        probe_candidates,
    )


def review_action(reason_code: str) -> tuple[str, str]:
    mapping = {
        "identity_conflict_family_member": (
            "keep_not_mapped",
            "Top lexical matches conflict with accession identity; do not auto-map.",
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
        "unresolved_name_strict_mode": (
            "resolve_identifier",
            "Resolve query to a UniProt accession or canonical gene symbol first.",
        ),
        "ambiguous_alias_resolution": (
            "resolve_ambiguity",
            "Alias maps to multiple accessions; disambiguate before mapping.",
        ),
        "no_measurement_term_for_resolved_accession": (
            "curate_new_term",
            "Resolved accession has no suitable local measurement term; consider cache curation.",
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
    top_n: int,
) -> int:
    pending_rows = [row for row in rows if row.get("validation") != "validated"]
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
        "source_validation",
        "source_evidence",
    ]
    count = 0
    reason_cache: dict[tuple[str, str], tuple[str, str, list[str], list[Candidate]]] = {}
    identity_cache: dict[tuple[str, tuple[str, ...]], dict[str, str]] = {}
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
                )
                reason_cache[reason_key] = cached_reason
            reason_code, reason, resolved_accessions, suggestions = cached_reason
            action, guidance = review_action(reason_code)
            identity_key = (query, tuple(resolved_accessions))
            identity = identity_cache.get(identity_key)
            if identity is None:
                identity = describe_query_identity(query, index, resolved_accessions)
                identity_cache[identity_key] = identity
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
                        "source_validation": row.get("validation", ""),
                        "source_evidence": row.get("evidence", ""),
                    }
                )
                count += 1
                continue

            for rank, cand in enumerate(top, start=1):
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
                        "source_validation": row.get("validation", ""),
                        "source_evidence": row.get("evidence", ""),
                    }
                )
                count += 1
    return count


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
            "source": "auto-writeback",
            "evidence": normalize(row["evidence"]),
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


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Bulk/offline mapper for measurement EFO terms")
    sub = parser.add_subparsers(dest="command")

    p_idx = sub.add_parser("index-build", help="Build local JSON index from cache/reference TSVs")
    p_idx.add_argument("--term-cache", default=str(DEFAULT_TERM_CACHE), help="Term cache TSV path")
    p_idx.add_argument("--analyte-cache", default=str(DEFAULT_ANALYTE_CACHE), help="Analyte->EFO cache TSV path")
    p_idx.add_argument("--uniprot-aliases", default=str(DEFAULT_UNIPROT_ALIASES), help="UniProt accession alias TSV")
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

    p_map = sub.add_parser("map", help="Map queries using local JSON index")
    p_map.add_argument(
        "--input",
        required=True,
        help=(
            "Input txt/csv/tsv with analyte queries. Optional tabular input_type column "
            "(input_type/query_type/id_type/type): accession, gene_symbol, gene_id, protein_name, auto"
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
            "Handling for non-accession free-text names/gene symbols: "
            "strict requires unique alias->UniProt resolution; fuzzy enables lexical fallback"
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
        help="Optional TSV path for non-validated rows with reason codes and review suggestions",
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

    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 2

    try:
        if args.command == "index-build":
            term_rows = load_term_cache(Path(args.term_cache))
            analyte_rows = load_analyte_cache(Path(args.analyte_cache))
            aliases = load_uniprot_aliases(Path(args.uniprot_aliases))
            index = build_index(term_rows=term_rows, analyte_rows=analyte_rows, accession_aliases=aliases)
            save_index(index, Path(args.output_index))
            print(
                f"[OK] built index at {args.output_index} "
                f"(terms={len(index.get('term_meta', {}))}, "
                f"analyte_keys={len(index.get('analyte_index', {}))}, "
                f"accessions={len(index.get('accession_alias_index', {}))})"
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
                index = build_index(
                    term_rows=term_rows,
                    analyte_rows=analyte_rows,
                    accession_aliases=aliases,
                )
                save_index(index, Path(args.index))
                print(
                    f"[OK] rebuilt index at {args.index} "
                    f"(terms={len(index.get('term_meta', {}))}, "
                    f"analyte_keys={len(index.get('analyte_index', {}))}, "
                    f"accessions={len(index.get('accession_alias_index', {}))})"
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
            term_cache_path = Path(args.term_cache)

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
                index = build_index(
                    term_rows=term_rows,
                    analyte_rows=analyte_rows,
                    accession_aliases=aliases,
                )
                save_index(index, index_path)
                print(
                    f"[OK] rebuilt index at {index_path} "
                    f"(terms={len(index.get('term_meta', {}))}, "
                    f"analyte_keys={len(index.get('analyte_index', {}))}, "
                    f"accessions={len(index.get('accession_alias_index', {}))})"
                )
            else:
                index = load_index(index_path)

            rows = map_queries(
                query_inputs,
                index=index,
                top_k=max(1, args.top_k),
                min_score=args.min_score,
                workers=max(1, args.workers),
                matrix_priority=parse_matrix_priority(args.matrix_priority),
                measurement_context=measurement_context,
                additional_contexts=additional_contexts,
                additional_context_keywords=additional_context_keywords,
                name_mode=args.name_mode,
                force_map_best=args.force_map_best,
                fallback_efo_id=args.fallback_efo_id,
                show_progress=args.progress,
                parallel_mode=args.parallel_mode,
                index_path=index_path,
            )
            write_tsv(rows, Path(args.output))
            if args.unmapped_output:
                unresolved = write_unmapped_tsv(rows, Path(args.unmapped_output))
                print(f"[OK] wrote {unresolved} unresolved rows to {args.unmapped_output}")
            if args.review_output:
                review_rows = write_review_tsv(
                    rows,
                    Path(args.review_output),
                    index=index,
                    min_score=args.min_score,
                    matrix_priority=parse_matrix_priority(args.matrix_priority),
                    measurement_context=measurement_context,
                    additional_contexts=additional_contexts,
                    additional_context_keywords=additional_context_keywords,
                    name_mode=args.name_mode,
                    top_n=max(1, args.review_top_n),
                )
                print(f"[OK] wrote {review_rows} review rows to {args.review_output}")
            if args.cache_writeback:
                added = write_back_cache(rows, Path(args.analyte_cache), args.cache_min_confidence)
                print(f"[OK] cache write-back added {added} entries to {args.analyte_cache}")
            print(f"[OK] wrote {len(rows)} rows to {args.output}")
            return 0

        raise ValueError(f"Unsupported command: {args.command}")
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
