"""Browser-safe trait mapper core for Pyodide worker execution."""

from __future__ import annotations

import csv
import io
import json
import re
from dataclasses import dataclass

ALLOWED_PREFIXES = ("EFO_", "MONDO_", "HP_", "OBA_")
QUERY_COL_CANDIDATES = ("query", "trait", "reported_trait", "term", "input_query")
SCALE_COL_CANDIDATES = ("trait_scale", "input_trait_scale", "scale")
STOPWORDS = {
    "a",
    "an",
    "and",
    "at",
    "by",
    "for",
    "from",
    "in",
    "of",
    "on",
    "or",
    "the",
    "to",
    "ukb",
    "data",
    "field",
    "gene",
    "based",
    "burden",
}
AGE_MARKERS = ("age at", "age of onset", "year when", "year of", "first diagnosed", "age diagnosed")
AGE_LABEL_MARKERS = ("age at", "age of onset", "age at diagnosis", "onset age")
NONSPECIFIC_LABELS = {
    "disease",
    "measurement",
    "area",
    "finding",
    "phenotype",
}


def _clean_text(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "").strip())


def _normalize_text(value: str) -> str:
    text = _clean_text(value).lower()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _tokenize(value: str) -> set[str]:
    norm = _normalize_text(value)
    tokens = {token for token in norm.split(" ") if token and token not in STOPWORDS}
    return tokens


def _split_multi(value: str) -> list[str]:
    if not value:
        return []
    parts = re.split(r"\s*\|\s*|\s*;\s*", value.strip())
    return [_clean_text(part) for part in parts if _clean_text(part)]


def _id_allowed(term_id: str) -> bool:
    return any(term_id.startswith(prefix) for prefix in ALLOWED_PREFIXES)


def _is_age_query(query: str) -> bool:
    qn = _normalize_text(query)
    return any(marker in qn for marker in AGE_MARKERS)


def _extract_source_code(query: str) -> str:
    # UKB data field style: "(UKB data field 20008)" or "(UKB data filed 20008)"
    match = re.search(r"ukb data fi(?:el|le)d\s+([0-9_#.-]+)", query.lower())
    if match:
        return match.group(1)
    # ICD10 range examples: "40001_C73-C75"
    match = re.search(r"\b([A-Z][0-9]{2}(?:\.[0-9]+)?(?:-[A-Z]?[0-9]{2}(?:\.[0-9]+)?)?)\b", query)
    if match:
        return match.group(1)
    return ""


@dataclass(frozen=True)
class CacheEntry:
    trait_id: str
    trait_label: str
    source_code: str
    lookup_text: str
    curated_trait: str
    search_keys: tuple[str, ...]
    tokens: frozenset[str]
    label_norm: str


class TraitMapper:
    def __init__(self, entries: list[CacheEntry]) -> None:
        self.entries = entries
        self.exact_index: dict[str, list[int]] = {}
        self.token_index: dict[str, list[int]] = {}
        for idx, entry in enumerate(entries):
            for key in entry.search_keys:
                self.exact_index.setdefault(key, []).append(idx)
            for token in entry.tokens:
                self.token_index.setdefault(token, []).append(idx)

    @classmethod
    def from_bundle_json(cls, payload_text: str) -> "TraitMapper":
        payload = json.loads(payload_text)
        raw_entries = payload.get("entries", [])
        entries: list[CacheEntry] = []
        for row in raw_entries:
            trait_id = _clean_text(str(row.get("id", "")))
            trait_label = _clean_text(str(row.get("label", "")))
            if not trait_id or not trait_label or not _id_allowed(trait_id):
                continue

            lookup_text = _clean_text(str(row.get("lookup_text", "")))
            curated_trait = _clean_text(str(row.get("curated_trait", "")))
            source_code = _clean_text(str(row.get("source_code", "")))
            synonyms = [_clean_text(str(v)) for v in row.get("synonyms", []) if _clean_text(str(v))]

            keys = []
            for candidate in (lookup_text, curated_trait, trait_label, *synonyms):
                norm = _normalize_text(candidate)
                if norm and norm not in keys:
                    keys.append(norm)

            tokens: set[str] = set()
            for candidate in (lookup_text, curated_trait, trait_label, *synonyms):
                tokens.update(_tokenize(candidate))

            if not keys and not tokens:
                continue

            entries.append(
                CacheEntry(
                    trait_id=trait_id,
                    trait_label=trait_label,
                    source_code=source_code,
                    lookup_text=lookup_text,
                    curated_trait=curated_trait,
                    search_keys=tuple(keys),
                    tokens=frozenset(tokens),
                    label_norm=_normalize_text(trait_label),
                )
            )
        return cls(entries)

    def _candidate_indexes(self, query: str, query_norm: str, query_tokens: set[str]) -> set[int]:
        candidates: set[int] = set()
        candidates.update(self.exact_index.get(query_norm, []))
        for token in query_tokens:
            for idx in self.token_index.get(token, []):
                candidates.add(idx)
        if candidates:
            return candidates
        # Fallback for very sparse tokenization.
        return set(range(len(self.entries)))

    def _score(
        self,
        entry: CacheEntry,
        query: str,
        query_norm: str,
        query_tokens: set[str],
        query_scale: str,
        is_age: bool,
        extracted_code: str,
    ) -> tuple[float, str]:
        score = 0.0
        reason = "lexical_tokens"

        if query_norm and query_norm in entry.search_keys:
            score = 1.0
            reason = "exact_text"
        else:
            for key in entry.search_keys:
                if query_norm and (query_norm in key or key in query_norm):
                    score = max(score, 0.90)
                    reason = "substring"
            if query_tokens and entry.tokens:
                overlap = len(query_tokens & entry.tokens)
                if overlap:
                    containment = overlap / max(1, len(query_tokens))
                    union = len(query_tokens | entry.tokens)
                    jaccard = overlap / max(1, union)
                    lexical = 0.30 + 0.45 * containment + 0.25 * jaccard
                    if lexical > score:
                        score = lexical
                        reason = "lexical_tokens"

        if extracted_code and entry.source_code and extracted_code.lower() == entry.source_code.lower():
            score = max(score, 0.98)
            reason = "source_code"

        if is_age:
            if any(marker in entry.label_norm for marker in AGE_LABEL_MARKERS):
                score += 0.10
                if reason == "lexical_tokens":
                    reason = "age_composite"
            else:
                score -= 0.10

        if query_scale == "quantitative":
            if {"measurement", "density", "mass", "level", "concentration"} & entry.tokens:
                score += 0.06
        elif query_scale == "binary":
            if "measurement" in entry.tokens and "disease" not in entry.tokens:
                score -= 0.06

        if entry.label_norm in NONSPECIFIC_LABELS:
            score -= 0.24

        if entry.label_norm in {"disease", "measurement", "area"}:
            score -= 0.08

        return max(0.0, min(1.0, score)), reason

    def map_query(self, query: str, query_scale: str = "auto") -> dict[str, str]:
        query_text = _clean_text(query)
        query_norm = _normalize_text(query_text)
        if not query_norm:
            return {
                "mapped_trait_id": "",
                "mapped_trait_label": "",
                "confidence": "0.000",
                "matched_via": "none",
                "evidence": "empty_query",
                "validation": "not_mapped",
            }

        query_tokens = _tokenize(query_text)
        is_age = _is_age_query(query_text)
        extracted_code = _extract_source_code(query_text)
        candidates = self._candidate_indexes(query_text, query_norm, query_tokens)

        best_idx = -1
        best_score = 0.0
        best_reason = "none"
        for idx in candidates:
            entry = self.entries[idx]
            score, reason = self._score(
                entry=entry,
                query=query_text,
                query_norm=query_norm,
                query_tokens=query_tokens,
                query_scale=(query_scale or "auto"),
                is_age=is_age,
                extracted_code=extracted_code,
            )
            if score > best_score:
                best_score = score
                best_idx = idx
                best_reason = reason

        if best_idx < 0:
            return {
                "mapped_trait_id": "",
                "mapped_trait_label": "",
                "confidence": "0.000",
                "matched_via": "none",
                "evidence": "no_candidates",
                "validation": "not_mapped",
            }

        best = self.entries[best_idx]
        return {
            "mapped_trait_id": best.trait_id,
            "mapped_trait_label": best.trait_label,
            "confidence": f"{best_score:.3f}",
            "matched_via": best_reason,
            "evidence": f"lookup={best.lookup_text or '-'};curated={best.curated_trait or '-'}",
            "validation": "validated",
        }


def _sniff_delimiter(text: str) -> str:
    header = text.splitlines()[0] if text.splitlines() else ""
    if "\t" in header:
        return "\t"
    if "," in header:
        return ","
    return "\t"


def _extract_query_from_row(row: dict[str, str], fallback_col: str) -> str:
    for key in QUERY_COL_CANDIDATES:
        if key in row and _clean_text(row.get(key, "")):
            return _clean_text(row.get(key, ""))
    return _clean_text(row.get(fallback_col, ""))


def _extract_scale_from_row(row: dict[str, str]) -> str:
    for key in SCALE_COL_CANDIDATES:
        val = _clean_text(row.get(key, "")).lower()
        if val in {"auto", "binary", "quantitative"}:
            return val
    return "auto"


def _parse_input_rows(input_tsv_text: str) -> list[dict[str, str]]:
    text = input_tsv_text or ""
    lines = [line for line in text.splitlines() if _clean_text(line)]
    if not lines:
        return []

    delimiter = _sniff_delimiter(text)
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    if not reader.fieldnames:
        return [{"query": _clean_text(line), "trait_scale": "auto"} for line in lines]

    fallback_col = reader.fieldnames[0]
    rows: list[dict[str, str]] = []
    for raw in reader:
        query = _extract_query_from_row(raw, fallback_col)
        if not query:
            continue
        rows.append({"query": query, "trait_scale": _extract_scale_from_row(raw)})

    if rows:
        return rows
    return [{"query": _clean_text(line), "trait_scale": "auto"} for line in lines]


def map_tsv(input_tsv_text: str, cache_json_text: str, min_score: float = 0.82) -> str:
    mapper = TraitMapper.from_bundle_json(cache_json_text)
    rows = _parse_input_rows(input_tsv_text)

    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "input_row_id",
            "input_query",
            "input_trait_scale",
            "mapped_trait_id",
            "mapped_trait_label",
            "confidence",
            "matched_via",
            "evidence",
            "validation",
        ]
    )

    threshold = max(0.0, min(1.0, float(min_score)))
    for idx, row in enumerate(rows, start=1):
        query = row.get("query", "")
        scale = row.get("trait_scale", "auto")
        mapped = mapper.map_query(query=query, query_scale=scale)
        score = float(mapped.get("confidence", "0.0") or "0.0")
        if score < threshold:
            mapped["mapped_trait_id"] = ""
            mapped["mapped_trait_label"] = ""
            mapped["validation"] = "not_mapped"
            mapped["matched_via"] = "below_threshold"
        writer.writerow(
            [
                idx,
                query,
                scale,
                mapped.get("mapped_trait_id", ""),
                mapped.get("mapped_trait_label", ""),
                mapped.get("confidence", "0.000"),
                mapped.get("matched_via", "none"),
                mapped.get("evidence", ""),
                mapped.get("validation", "not_mapped"),
            ]
        )
    return out.getvalue()
