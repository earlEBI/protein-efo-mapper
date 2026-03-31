#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

DEFAULT_BASELINE = Path("/Users/earl/Downloads/PMID34662886_studies_export.tsv")
DEFAULT_CURRENT = Path("/tmp/PMID34662886_reportedTrait_mapped_after_seed_lock_full_20260321c5k.tsv")
DEFAULT_EFO_OBO = Path("/Users/earl/src/analyte-efo-mapper/skills/pqtl-measurement-mapper/references/efo.obo")

DEFAULT_OUT_ALL = Path("/tmp/PMID34662886_worse_vs_curator_all_rows_20260321.tsv")
DEFAULT_OUT_SURE = Path("/tmp/PMID34662886_worse_vs_curator_sure_20260321.tsv")
DEFAULT_OUT_REVIEW = Path("/tmp/PMID34662886_worse_vs_curator_review_20260321.tsv")

ICD10_RE = re.compile(r"ICD10\s+([A-Z][0-9][0-9A-Z]*(?:\.[0-9A-Z]+)?)\s*:\s*([^\(]+)", re.IGNORECASE)
MULTI_SPLIT_RE = re.compile(r"\s*(?:\||,|;|/|\band\b|\bor\b)\s*", re.IGNORECASE)

GENERIC_LABELS = {
    "disease",
    "disorder",
    "medical procedure",
    "planned process",
    "test result",
    "wellbeing measurement",
    "drug use measurement",
    "clinical treatment",
    "treatment",
    "injury",
    "sign or symptom",
    "mental or behavioural disorder",
    "mental disorder",
}
WEAK_CURRENT_MATCH_VIA = {
    "efo_obo_fuzzy",
    "cache_fuzzy_text",
    "low_confidence_reasoned_fallback",
    "cache_exact_text_rescued",
}
TOKEN_STOP = {
    "a",
    "an",
    "and",
    "at",
    "by",
    "code",
    "disease",
    "disorder",
    "for",
    "gene",
    "based",
    "burden",
    "icd10",
    "in",
    "is",
    "not",
    "of",
    "or",
    "other",
    "reported",
    "self",
    "specified",
    "trait",
    "type",
    "ukb",
    "field",
    "with",
    "without",
}


@dataclass(frozen=True)
class Resolution:
    status: str
    resolved_id: str
    resolved_label: str


def normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", (text or "").strip())


def normalize_lookup(text: str) -> str:
    value = normalize_space(text).lower().replace('"', "").replace("'", "")
    value = re.sub(r"[^a-z0-9\s\-|]", " ", value)
    return normalize_space(value)


def normalize_term_id(raw: str) -> str:
    token = (raw or "").strip()
    if token.startswith("efo:EFO_"):
        return token.split(":", 1)[1]
    if token.startswith("EFO:"):
        return token.replace(":", "_")
    return token.replace(":", "_")


def normalize_pipe(text: str) -> str:
    value = normalize_lookup(text)
    if not value:
        return ""
    parts = [p.strip() for p in re.split(r"\s*\|\s*", value) if p.strip()]
    return "|".join(sorted(parts))


def split_current_ids(text: str) -> set[str]:
    value = normalize_space(text)
    if not value:
        return set()
    ids: set[str] = set()
    for raw in re.split(r"[|;,]", value):
        term_id = normalize_term_id(raw)
        if term_id:
            ids.add(term_id)
    return ids


def parse_icd10(reported_trait: str) -> tuple[str, str]:
    match = ICD10_RE.search(reported_trait or "")
    if not match:
        return "", ""
    return match.group(1).upper().strip(), normalize_space(match.group(2))


def strip_query_noise(text: str) -> str:
    value = normalize_space(text)
    value = re.sub(r"\([^)]*\)", " ", value)
    value = ICD10_RE.sub(r"\2", value)
    value = re.sub(r"\bUKB data field\s*[0-9]{1,7}\b", " ", value, flags=re.IGNORECASE)
    value = re.sub(r"\bGene-based burden\b", " ", value, flags=re.IGNORECASE)
    return normalize_space(value)


def tokenize(text: str) -> set[str]:
    tokens: set[str] = set()
    for tok in re.findall(r"[a-z0-9]+", normalize_lookup(text)):
        if tok in TOKEN_STOP:
            continue
        if len(tok) <= 2:
            continue
        tokens.add(tok)
    return tokens


def f1_similarity(lhs: str, rhs: str) -> float:
    lt = tokenize(lhs)
    rt = tokenize(rhs)
    if not lt or not rt:
        return 0.0
    inter = len(lt & rt)
    if inter == 0:
        return 0.0
    precision = inter / len(lt)
    recall = inter / len(rt)
    return 2.0 * precision * recall / (precision + recall)


def is_generic_label(label: str) -> bool:
    key = normalize_lookup(label)
    if not key:
        return True
    if key in GENERIC_LABELS:
        return True
    parts = [p for p in MULTI_SPLIT_RE.split(key) if p]
    if len(parts) == 1 and len(tokenize(parts[0])) <= 1:
        return True
    return False


def labels_semantically_equivalent(manual_label: str, current_label: str) -> bool:
    manual_key = normalize_lookup(manual_label)
    current_key = normalize_lookup(current_label)
    if not manual_key or not current_key:
        return False
    if manual_key == current_key:
        return True
    pair_aliases = (
        ("whooping cough", "pertussis"),
        ("barrett's esophagus", "barrett esophagus"),
        ("barrett oesophagus", "barrett esophagus"),
    )
    for left, right in pair_aliases:
        if (left in manual_key and right in current_key) or (right in manual_key and left in current_key):
            return True
    return False


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{(k or "").strip(): (v or "").strip() for k, v in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_efo_exact_index(efo_obo_path: Path) -> dict[str, set[tuple[str, str]]]:
    index: dict[str, set[tuple[str, str]]] = defaultdict(set)
    term_id = ""
    term_name = ""
    exact_synonyms: list[str] = []
    in_term = False

    def flush_term() -> None:
        nonlocal term_id, term_name, exact_synonyms
        if term_id:
            canonical_label = term_name or ""
            for text in [term_name] + exact_synonyms:
                norm = normalize_lookup(text)
                if norm:
                    index[norm].add((term_id, canonical_label))
        term_id = ""
        term_name = ""
        exact_synonyms = []

    with efo_obo_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                flush_term()
                in_term = True
                continue
            if not in_term:
                continue
            if not line:
                flush_term()
                in_term = False
                continue
            if line.startswith("id: "):
                term_id = normalize_term_id(line[4:].strip())
                continue
            if line.startswith("name: "):
                term_name = line[6:].strip()
                continue
            if line.startswith("synonym: "):
                match = re.match(r'synonym:\s+"(.+?)"\s+(EXACT|BROAD|NARROW|RELATED)', line)
                if match and match.group(2) == "EXACT":
                    exact_synonyms.append(match.group(1).strip())
    flush_term()
    return index


def resolve_manual_term(manual_efo: str, efo_exact_index: dict[str, set[tuple[str, str]]]) -> Resolution:
    text = normalize_space(manual_efo)
    if not text:
        return Resolution("blank", "", "")
    if MULTI_SPLIT_RE.search(text):
        return Resolution("multi", "", "")
    norm = normalize_lookup(text)
    matches = sorted(efo_exact_index.get(norm, set()))
    if len(matches) == 1:
        return Resolution("unique", matches[0][0], matches[0][1] or text)
    if len(matches) > 1:
        return Resolution("ambiguous", "", "")
    return Resolution("unresolved", "", "")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find rows where current mapping is worse than curator mapping (all rows).")
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--current", type=Path, default=DEFAULT_CURRENT)
    parser.add_argument("--efo-obo", type=Path, default=DEFAULT_EFO_OBO)
    parser.add_argument("--out-all", type=Path, default=DEFAULT_OUT_ALL)
    parser.add_argument("--out-sure", type=Path, default=DEFAULT_OUT_SURE)
    parser.add_argument("--out-review", type=Path, default=DEFAULT_OUT_REVIEW)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    baseline_rows = read_tsv(args.baseline)
    current_rows = read_tsv(args.current)
    if len(baseline_rows) != len(current_rows):
        raise ValueError(f"Row count mismatch: baseline={len(baseline_rows)} current={len(current_rows)}")

    efo_exact_index = build_efo_exact_index(args.efo_obo)

    out_rows: list[dict[str, str]] = []
    sure_rows: list[dict[str, str]] = []
    review_rows: list[dict[str, str]] = []
    decision_counts: Counter[str] = Counter()

    for rownum, (base, cur) in enumerate(zip(baseline_rows, current_rows), start=1):
        manual = normalize_space(base.get("efoTraits", ""))
        current_label = normalize_space(cur.get("mapped_trait_label", ""))
        current_id = normalize_space(cur.get("mapped_trait_id", ""))
        if normalize_pipe(manual) == normalize_pipe(current_label):
            continue

        reported_trait = normalize_space(base.get("reportedTrait", ""))
        query_core = strip_query_noise(reported_trait)
        icd10, icd10_label = parse_icd10(reported_trait)
        if icd10_label:
            query_core = strip_query_noise(icd10_label)

        sim_manual = f1_similarity(manual, query_core)
        sim_current = f1_similarity(current_label, query_core)
        margin = sim_manual - sim_current

        resolution = resolve_manual_term(manual, efo_exact_index)
        current_id_set = split_current_ids(current_id)
        manual_id_equivalent = (
            resolution.status == "unique"
            and bool(current_id_set)
            and normalize_term_id(resolution.resolved_id) in current_id_set
        )
        if manual_id_equivalent:
            continue
        if labels_semantically_equivalent(manual, current_label):
            continue
        current_generic = is_generic_label(current_label)
        manual_generic = is_generic_label(manual)
        matched_via = normalize_space(cur.get("matched_via", ""))
        matched_via_key = normalize_lookup(matched_via)
        validation = normalize_space(cur.get("validation", ""))

        decision = "not_worse_or_equivalent"
        reason = "no_strict_worse_signal"
        confidence = 0.55

        if manual and not current_label:
            decision = "worse_sure"
            reason = "current_unmapped_manual_present"
            confidence = 0.99
        elif (
            resolution.status == "unique"
            and current_id_set
            and normalize_term_id(resolution.resolved_id) not in current_id_set
        ):
            if sim_manual >= 0.45 and margin >= 0.22:
                decision = "worse_sure"
                reason = "manual_unique_id_with_strong_query_similarity_gain"
                confidence = min(0.99, 0.78 + margin * 0.35)
        if decision == "not_worse_or_equivalent":
            if current_generic and not manual_generic and sim_manual >= max(0.35, sim_current + 0.10):
                decision = "worse_sure"
                reason = "current_generic_manual_more_specific"
                confidence = min(0.97, 0.72 + max(margin, 0.0) * 0.3)
        if decision == "not_worse_or_equivalent":
            if matched_via_key in WEAK_CURRENT_MATCH_VIA and sim_current < 0.30 and sim_manual >= 0.40:
                decision = "worse_sure"
                reason = "weak_current_method_and_manual_has_stronger_query_overlap"
                confidence = min(0.95, 0.7 + max(margin, 0.0) * 0.35)
        if decision == "not_worse_or_equivalent":
            manual_tokens = tokenize(manual)
            current_tokens = tokenize(current_label)
            query_tokens = tokenize(query_core)
            manual_hits = len(manual_tokens & query_tokens)
            current_hits = len(current_tokens & query_tokens)
            if (
                not manual_generic
                and normalize_lookup(current_label) not in normalize_lookup(manual)
                and manual_hits >= current_hits + 2
                and manual_hits >= 2
            ):
                decision = "worse_sure"
                reason = "manual_anchor_token_match_exceeds_current"
                confidence = 0.78

        if decision == "not_worse_or_equivalent":
            if margin >= 0.10 and sim_manual >= 0.35:
                decision = "potential_worse_review"
                reason = "manual_overlap_higher_but_not_strict"
                confidence = 0.5
            elif current_generic and not manual_generic:
                decision = "potential_worse_review"
                reason = "possible_genericity_regression"
                confidence = 0.45

        row = {
            "rownum": str(rownum),
            "reportedTrait": reported_trait,
            "icd10": icd10,
            "manual_efoTraits": manual,
            "current_mapped_trait_id": current_id,
            "current_mapped_trait_label": current_label,
            "matched_via": matched_via,
            "validation": validation,
            "decision": decision,
            "decision_reason": reason,
            "confidence": f"{confidence:.3f}",
            "manual_query_similarity": f"{sim_manual:.3f}",
            "current_query_similarity": f"{sim_current:.3f}",
            "similarity_margin": f"{margin:.3f}",
            "manual_resolution": resolution.status,
            "manual_resolved_id": resolution.resolved_id,
            "manual_resolved_label": resolution.resolved_label,
            "current_is_generic": "yes" if current_generic else "no",
            "manual_is_generic": "yes" if manual_generic else "no",
        }
        out_rows.append(row)
        decision_counts[decision] += 1
        if decision == "worse_sure":
            sure_rows.append(row)
        elif decision == "potential_worse_review":
            review_rows.append(row)

    fieldnames = [
        "rownum",
        "reportedTrait",
        "icd10",
        "manual_efoTraits",
        "current_mapped_trait_id",
        "current_mapped_trait_label",
        "matched_via",
        "validation",
        "decision",
        "decision_reason",
        "confidence",
        "manual_query_similarity",
        "current_query_similarity",
        "similarity_margin",
        "manual_resolution",
        "manual_resolved_id",
        "manual_resolved_label",
        "current_is_generic",
        "manual_is_generic",
    ]
    write_tsv(args.out_all, out_rows, fieldnames)
    write_tsv(args.out_sure, sure_rows, fieldnames)
    write_tsv(args.out_review, review_rows, fieldnames)

    print("Worse-vs-curator audit complete")
    print(f"baseline_rows={len(baseline_rows)} current_rows={len(current_rows)}")
    print(f"changed_rows={len(out_rows)}")
    print("decisions=" + ",".join(f"{k}:{v}" for k, v in sorted(decision_counts.items(), key=lambda kv: kv[0])))
    print(f"worse_sure_rows={len(sure_rows)}")
    print(f"potential_worse_review_rows={len(review_rows)}")
    print(f"out_all={args.out_all}")
    print(f"out_sure={args.out_sure}")
    print(f"out_review={args.out_review}")


if __name__ == "__main__":
    main()
