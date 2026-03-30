#!/usr/bin/env python3
"""Build a browser-friendly trait cache bundle for the Pyodide mapper."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

ALLOWED_PREFIXES = ("EFO_", "MONDO_", "HP_", "OBA_")


def _clean(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "").strip())


def _split_multi(value: str) -> list[str]:
    if not value:
        return []
    parts = re.split(r"\s*\|\s*|\s*;\s*", value.strip())
    return [_clean(p) for p in parts if _clean(p)]


def _id_allowed(term_id: str) -> bool:
    return any(term_id.startswith(prefix) for prefix in ALLOWED_PREFIXES)


def _extract_source_code(row: dict[str, str]) -> str:
    field = _clean(row.get("UKBB data field", ""))
    if field and field != "-":
        return field
    for key in ("ICD10", "PheCode"):
        val = _clean(row.get(key, ""))
        if val and val != "-":
            return val
    return ""


def _parse_row(row: dict[str, str]) -> list[dict[str, object]]:
    ids = _split_multi(row.get("main ontology URI(s)", ""))
    labels = _split_multi(row.get("main ontology label(s)", ""))
    lookup_text = _clean(row.get("Lookup text", ""))
    curated_trait = _clean(row.get("Curated reported trait", ""))
    user_trait = _clean(row.get("User submitted trait", ""))
    source_code = _extract_source_code(row)

    entries: list[dict[str, object]] = []
    for idx, term_id in enumerate(ids):
        if not _id_allowed(term_id):
            continue
        label = labels[idx] if idx < len(labels) else (labels[-1] if labels else "")
        label = _clean(label)
        if not label:
            continue
        synonyms = [v for v in (user_trait,) if v and v != "-"]
        entries.append(
            {
                "id": term_id,
                "label": label,
                "lookup_text": lookup_text,
                "curated_trait": curated_trait,
                "source_code": source_code,
                "synonyms": synonyms,
            }
        )
    return entries


def _dedupe(entries: list[dict[str, object]]) -> list[dict[str, object]]:
    seen: set[tuple[str, str, str, str, str]] = set()
    out: list[dict[str, object]] = []
    for entry in entries:
        key = (
            str(entry["id"]),
            str(entry["label"]),
            str(entry["lookup_text"]),
            str(entry["curated_trait"]),
            str(entry["source_code"]),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(entry)
    return out


def _matches_keyword(entry: dict[str, object], keywords: list[str]) -> bool:
    if not keywords:
        return False
    haystack = " ".join(
        (
            str(entry.get("lookup_text", "")),
            str(entry.get("curated_trait", "")),
            str(entry.get("label", "")),
        )
    ).lower()
    return any(keyword in haystack for keyword in keywords)


def build_bundle(
    trait_cache_path: Path,
    output_path: Path,
    limit: int,
    keywords: list[str],
) -> tuple[int, int]:
    with trait_cache_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        parsed: list[dict[str, object]] = []
        for row in reader:
            parsed.extend(_parse_row(row))

    deduped = _dedupe(parsed)

    keywords_norm = [_clean(k).lower() for k in keywords if _clean(k)]
    prioritized: list[dict[str, object]] = []
    remainder: list[dict[str, object]] = []
    for entry in deduped:
        if _matches_keyword(entry, keywords_norm):
            prioritized.append(entry)
        else:
            remainder.append(entry)

    merged = prioritized + remainder
    if limit > 0:
        merged = merged[:limit]

    payload = {
        "version": "pyodide-trait-bundle-v1",
        "source": trait_cache_path.name,
        "allowed_prefixes": list(ALLOWED_PREFIXES),
        "entry_count": len(merged),
        "entries": merged,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"))
    return len(parsed), len(merged)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build browser trait cache bundle for Pyodide.")
    parser.add_argument(
        "--trait-cache",
        default="skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv",
        help="Path to trait mapping cache TSV",
    )
    parser.add_argument(
        "--output",
        default="skills/pqtl-measurement-mapper/web/pyodide/cache_bundle.json",
        help="Output JSON path",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional max entries (0 = all)",
    )
    parser.add_argument(
        "--keyword",
        action="append",
        default=[],
        help="Prioritize entries containing this keyword in lookup/trait/label",
    )
    args = parser.parse_args()

    parsed_count, kept_count = build_bundle(
        trait_cache_path=Path(args.trait_cache),
        output_path=Path(args.output),
        limit=max(0, int(args.limit)),
        keywords=list(args.keyword),
    )
    print(f"[ok] parsed={parsed_count} kept={kept_count} output={args.output}")


if __name__ == "__main__":
    main()
