#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

DEFAULT_BASELINE = Path("/Users/earl/Downloads/PMID34662886_studies_export.tsv")
DEFAULT_CURRENT = Path("/tmp/PMID34662886_reportedTrait_mapped_with_labelsyn_appendonly_v2.tsv")
DEFAULT_EFO_OBO = Path("/Users/earl/src/analyte-efo-mapper/skills/pqtl-measurement-mapper/references/efo.obo")
DEFAULT_CACHE_TEMPLATE = Path("/Users/earl/src/analyte-efo-mapper/references/icd10/icd10_trait_supplement_cache.tsv")

DEFAULT_OUT_ALL = Path("/tmp/PMID34662886_changed_all_rows_with_decision.tsv")
DEFAULT_OUT_MANUAL_BETTER = Path("/tmp/PMID34662886_icd10_manual_better_strict.tsv")
DEFAULT_OUT_AMBIGUOUS = Path("/tmp/PMID34662886_icd10_ambiguous_review.tsv")
DEFAULT_OUT_CACHE_PREVIEW = Path("/tmp/PMID34662886_icd10_cache_patch_preview.tsv")

ICD10_RE = re.compile(r"ICD10\s+([A-Z][0-9][0-9A-Z]*(?:\.[0-9A-Z]+)?)\s*:\s*([^\(]+)", re.IGNORECASE)
MULTI_TERM_RE = re.compile(r"\s*(?:\||,|/|;|\band\b|\bor\b)\s*", re.IGNORECASE)

STOPWORDS = {
    "a",
    "an",
    "and",
    "as",
    "at",
    "by",
    "classified",
    "condition",
    "conditions",
    "disease",
    "diseases",
    "disorder",
    "disorders",
    "due",
    "elsewhere",
    "for",
    "from",
    "history",
    "in",
    "injuries",
    "injury",
    "is",
    "it",
    "not",
    "of",
    "on",
    "or",
    "other",
    "otherwise",
    "personal",
    "problem",
    "problems",
    "related",
    "sequelae",
    "specified",
    "state",
    "system",
    "the",
    "to",
    "type",
    "unspecified",
    "use",
    "with",
    "without",
}

TOKEN_CANON = {
    "abdominal": "abdomen",
    "alcoholic": "alcohol",
    "arthralgia": "joint",
    "cancer": "neoplasm",
    "cancers": "neoplasm",
    "cardiac": "heart",
    "cervical": "cervix",
    "colonic": "colon",
    "diarrhoea": "diarrhea",
    "gastric": "stomach",
    "hepatic": "liver",
    "hypotensive": "hypotension",
    "intestinal": "intestine",
    "joints": "joint",
    "myopic": "myopia",
    "neoplasms": "neoplasm",
    "orthostatic": "orthostatic",
    "phlebitis": "phlebitis",
    "rectal": "rectum",
    "renal": "kidney",
    "respiratory": "respiratory",
    "thrombophlebitis": "phlebitis",
    "vascular": "vessel",
}

GENERIC_LABELS = {
    "abdominal symptom",
    "adverse effect",
    "arthropathy",
    "complication",
    "connective tissue disorder",
    "device complication",
    "disease",
    "digestive system disorder",
    "encounter with health service for adjustment and management of implanted device",
    "female reproductive system disorder",
    "follow-up",
    "hearing loss disorder",
    "head injury",
    "hypotensive disorder",
    "infectious disease",
    "inflammatory bowel disease 1",
    "liver disorder",
    "muscular disease",
    "nontoxic goiter",
    "pain",
    "pernicious anemia",
    "phlebitis",
    "planned process",
    "pleural disorder",
    "respiratory system disorder",
    "rhinitis",
    "skin disorder",
}


@dataclass(frozen=True)
class Resolution:
    status: str
    resolved_id: str
    resolved_label: str
    reason: str


def normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", (text or "").strip())


def normalize_pipe(text: str) -> str:
    value = normalize_space(text).lower().replace('"', "").replace("'", "")
    if not value:
        return ""
    parts = [p.strip() for p in re.split(r"\s*\|\s*", value) if p.strip()]
    return "|".join(sorted(parts))


def normalize_lookup(text: str) -> str:
    value = normalize_space(text).lower().replace('"', "").replace("'", "")
    value = re.sub(r"[^a-z0-9\s\-]", " ", value)
    return normalize_space(value)


def normalize_term_id(raw: str) -> str:
    token = (raw or "").strip()
    if token.startswith("efo:EFO_"):
        return token.split(":", 1)[1]
    if token.startswith("EFO:"):
        return token.replace(":", "_")
    return token.replace(":", "_")


def singularize(token: str) -> str:
    if token.endswith("ies") and len(token) > 4:
        return token[:-3] + "y"
    if token.endswith("es") and len(token) > 4:
        return token[:-2]
    if token.endswith("s") and len(token) > 3:
        return token[:-1]
    return token


def tokenize(text: str) -> set[str]:
    tokens: set[str] = set()
    for tok in re.findall(r"[a-z0-9]+", (text or "").lower()):
        if tok in STOPWORDS:
            continue
        if len(tok) <= 2:
            continue
        tok = singularize(tok)
        tok = TOKEN_CANON.get(tok, tok)
        if tok in STOPWORDS or len(tok) <= 2:
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


def contains_phrase(lhs: str, rhs: str) -> bool:
    l = normalize_lookup(lhs)
    r = normalize_lookup(rhs)
    if not l or not r:
        return False
    return l in r or r in l


def label_specificity(label: str) -> int:
    return len(tokenize(label))


def is_generic_label(label: str) -> bool:
    nl = normalize_lookup(label)
    if not nl:
        return True
    if nl in GENERIC_LABELS:
        return True
    tok = tokenize(nl)
    if not tok:
        return True
    if len(tok) <= 2 and any(t in {"disease", "disorder", "syndrome", "condition", "symptom"} for t in tok):
        return True
    return False


def parse_icd10_from_reported_trait(reported_trait: str) -> tuple[str, str]:
    match = ICD10_RE.search(reported_trait or "")
    if not match:
        return "", ""
    code = match.group(1).upper().strip()
    label = normalize_space(match.group(2))
    return code, label


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


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def build_efo_exact_index(efo_obo_path: Path) -> dict[str, set[tuple[str, str]]]:
    index: dict[str, set[tuple[str, str]]] = defaultdict(set)

    term_id = ""
    term_name = ""
    exact_synonyms: list[str] = []
    in_term = False

    def flush_term() -> None:
        nonlocal term_id, term_name, exact_synonyms
        if not term_id:
            term_id = ""
            term_name = ""
            exact_synonyms = []
            return
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
        return Resolution("unresolved", "", "", "manual_blank")

    if MULTI_TERM_RE.search(text):
        return Resolution("ambiguous", "", "", "manual_multi_component")

    norm = normalize_lookup(text)
    matches = sorted(efo_exact_index.get(norm, set()))
    if not matches:
        return Resolution("unresolved", "", "", "manual_not_in_efo_exact")
    if len(matches) > 1:
        ids = "|".join(m[0] for m in matches)
        return Resolution("ambiguous", "", "", f"manual_maps_to_multiple_ids:{ids}")

    resolved_id, resolved_label = matches[0]
    return Resolution("unique", resolved_id, resolved_label or text, "manual_exact_resolved")


def decide_changed_row(
    rownum: int,
    baseline_row: dict[str, str],
    current_row: dict[str, str],
    resolution: Resolution,
) -> dict[str, str]:
    reported_trait = baseline_row.get("reportedTrait", "")
    manual_efo = baseline_row.get("efoTraits", "")

    current_id = current_row.get("mapped_trait_id", "")
    current_label = current_row.get("mapped_trait_label", "")
    matched_via = current_row.get("matched_via", "")
    validation = current_row.get("validation", "")

    icd10, icd10_label = parse_icd10_from_reported_trait(reported_trait)

    manual_score = f1_similarity(manual_efo, icd10_label)
    current_score = f1_similarity(current_label, icd10_label)
    if contains_phrase(manual_efo, icd10_label):
        manual_score = min(1.0, manual_score + 0.2)
    if contains_phrase(current_label, icd10_label):
        current_score = min(1.0, current_score + 0.2)

    score_margin = manual_score - current_score
    current_generic = is_generic_label(current_label)
    equivalent_by_label = normalize_lookup(manual_efo) == normalize_lookup(current_label)
    equivalent_by_id = resolution.status == "unique" and resolution.resolved_id == current_id

    decision = "current_better_or_equivalent"
    reason = ""
    confidence = 0.6

    if not icd10:
        decision = "current_better_or_equivalent"
        reason = "non_icd10_row_out_of_scope_for_strict_icd10_override"
        confidence = 0.7
    elif equivalent_by_id or equivalent_by_label:
        decision = "current_better_or_equivalent"
        reason = "manual_and_current_equivalent"
        confidence = 0.99
    elif resolution.status != "unique":
        decision = "ambiguous_review"
        reason = resolution.reason
        confidence = 0.4
    elif not current_id:
        decision = "manual_better_sure"
        reason = "current_unmapped_manual_exact_resolved"
        confidence = 0.98
    else:
        manual_specific = label_specificity(manual_efo)
        current_specific = label_specificity(current_label)

        if manual_score >= 0.75 and score_margin >= 0.2 and manual_specific >= current_specific:
            decision = "manual_better_sure"
            reason = "strong_icd10_lexical_gain"
            confidence = min(0.99, 0.72 + (manual_score * 0.2) + max(score_margin, 0.0) * 0.2)
        elif current_generic and manual_score >= 0.55 and score_margin >= 0.1:
            decision = "manual_better_sure"
            reason = "current_label_generic_manual_more_specific"
            confidence = min(0.99, 0.7 + (manual_score * 0.2) + max(score_margin, 0.0) * 0.15)
        elif manual_score >= 0.55 and current_score < 0.35 and manual_specific > current_specific + 1:
            decision = "manual_better_sure"
            reason = "manual_specificity_gain"
            confidence = min(0.99, 0.7 + (manual_score * 0.2) + max(score_margin, 0.0) * 0.15)
        elif score_margin >= 0.1:
            decision = "ambiguous_review"
            reason = "manual_looks_better_but_margin_not_strict_enough"
            confidence = 0.55
        else:
            decision = "current_better_or_equivalent"
            reason = "no_strict_manual_advantage"
            confidence = 0.7

    return {
        "rownum": str(rownum),
        "reportedTrait": reported_trait,
        "icd10": icd10,
        "icd10_label": icd10_label,
        "manual_efoTraits": manual_efo,
        "current_mapped_trait_id": current_id,
        "current_mapped_trait_label": current_label,
        "matched_via": matched_via,
        "validation": validation,
        "decision": decision,
        "decision_reason": reason,
        "manual_resolved_id": resolution.resolved_id,
        "manual_resolved_label": resolution.resolved_label,
        "confidence": f"{confidence:.3f}",
        "manual_icd10_score": f"{manual_score:.3f}",
        "current_icd10_score": f"{current_score:.3f}",
        "score_margin": f"{score_margin:.3f}",
        "current_is_generic": "yes" if current_generic else "no",
    }


def enforce_icd10_consistency(rows: list[dict[str, str]]) -> None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["decision"] == "manual_better_sure" and row["icd10"]:
            grouped[row["icd10"]].append(row)

    for icd10, icd_rows in grouped.items():
        unique_targets = {(r["manual_resolved_id"], r["manual_resolved_label"]) for r in icd_rows}
        if len(unique_targets) <= 1:
            continue
        for row in icd_rows:
            row["decision"] = "ambiguous_review"
            row["decision_reason"] = f"{row['decision_reason']};conflicting_manual_target_for_icd10:{icd10}"
            row["confidence"] = "0.450"


def load_cache_header(cache_template: Path) -> list[str]:
    with cache_template.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in cache template: {cache_template}")
        return list(reader.fieldnames)


def build_cache_preview_rows(
    strict_rows: list[dict[str, str]],
    cache_header: list[str],
) -> tuple[list[dict[str, str]], list[str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in strict_rows:
        grouped[row["icd10"]].append(row)

    preview_rows: list[dict[str, str]] = []
    skipped_conflicts: list[str] = []

    for icd10 in sorted(grouped):
        rows = grouped[icd10]
        targets = {(r["manual_resolved_id"], r["manual_resolved_label"]) for r in rows}
        if len(targets) != 1:
            skipped_conflicts.append(icd10)
            continue

        # Deterministic pick: highest confidence then lowest rownum.
        sorted_rows = sorted(
            rows,
            key=lambda r: (-float(r["confidence"]), int(r["rownum"])),
        )
        best = sorted_rows[0]

        row = {key: "" for key in cache_header}
        row.update(
            {
                "UKBB data field": "-",
                "ICD10": icd10,
                "PheCode": "-",
                "Lookup text": f"ICD10 {icd10} {normalize_lookup(best['icd10_label'])}".strip(),
                "Curated reported trait": best["reportedTrait"],
                "main ontology URI(s)": best["manual_resolved_id"],
                "main ontology label(s)": best["manual_resolved_label"],
                "Publication": "PMID34662886_manual_qc_strict_preview",
                "User submitted trait": "-",
                "qc_status": "review_required",
                "qc_resolution_modes": "manual_curator_baseline_strict",
                "qc_label_match_modes": "manual_exact_efo_label_or_exact_synonym",
                "qc_unresolved_ids": "",
                "final_update_applied": "no",
                "final_update_reason": "preview_only_not_applied",
            }
        )
        preview_rows.append(row)

    return preview_rows, skipped_conflicts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Strict manual-vs-current audit for PMID34662886 (report-only)."
    )
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE, help="Curator baseline TSV")
    parser.add_argument("--current", type=Path, default=DEFAULT_CURRENT, help="Current mapper output TSV")
    parser.add_argument("--efo-obo", type=Path, default=DEFAULT_EFO_OBO, help="EFO OBO path")
    parser.add_argument(
        "--cache-template",
        type=Path,
        default=DEFAULT_CACHE_TEMPLATE,
        help="ICD10 cache template for preview header",
    )
    parser.add_argument("--out-all", type=Path, default=DEFAULT_OUT_ALL)
    parser.add_argument("--out-manual-better", type=Path, default=DEFAULT_OUT_MANUAL_BETTER)
    parser.add_argument("--out-ambiguous", type=Path, default=DEFAULT_OUT_AMBIGUOUS)
    parser.add_argument("--out-cache-preview", type=Path, default=DEFAULT_OUT_CACHE_PREVIEW)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    baseline_rows = read_tsv(args.baseline)
    current_rows = read_tsv(args.current)

    if len(baseline_rows) != len(current_rows):
        raise ValueError(
            f"Row count mismatch: baseline={len(baseline_rows)} current={len(current_rows)}"
        )

    efo_index = build_efo_exact_index(args.efo_obo)

    changed_rows: list[dict[str, str]] = []
    for i, (baseline_row, current_row) in enumerate(zip(baseline_rows, current_rows), start=1):
        manual_efo = baseline_row.get("efoTraits", "")
        current_label = current_row.get("mapped_trait_label", "")

        if normalize_pipe(manual_efo) == normalize_pipe(current_label):
            continue

        resolution = resolve_manual_term(manual_efo, efo_index)
        decided = decide_changed_row(i, baseline_row, current_row, resolution)
        changed_rows.append(decided)

    enforce_icd10_consistency(changed_rows)

    strict_manual_better = [
        row
        for row in changed_rows
        if row["icd10"] and row["decision"] == "manual_better_sure"
    ]
    ambiguous_icd10 = [
        row
        for row in changed_rows
        if row["icd10"] and row["decision"] == "ambiguous_review"
    ]

    all_fieldnames = [
        "rownum",
        "reportedTrait",
        "icd10",
        "icd10_label",
        "manual_efoTraits",
        "current_mapped_trait_id",
        "current_mapped_trait_label",
        "matched_via",
        "validation",
        "decision",
        "decision_reason",
        "manual_resolved_id",
        "manual_resolved_label",
        "confidence",
        "manual_icd10_score",
        "current_icd10_score",
        "score_margin",
        "current_is_generic",
    ]

    write_tsv(args.out_all, changed_rows, all_fieldnames)
    write_tsv(args.out_manual_better, strict_manual_better, all_fieldnames)
    write_tsv(args.out_ambiguous, ambiguous_icd10, all_fieldnames)

    cache_header = load_cache_header(args.cache_template)
    cache_preview_rows, skipped_conflicts = build_cache_preview_rows(strict_manual_better, cache_header)
    write_tsv(args.out_cache_preview, cache_preview_rows, cache_header)

    decision_counts = Counter(row["decision"] for row in changed_rows)
    changed_icd10_rows = [row for row in changed_rows if row["icd10"]]
    unique_preview_codes = {row.get("ICD10", "") for row in cache_preview_rows}

    print("Audit complete")
    print(f"baseline_rows={len(baseline_rows)} current_rows={len(current_rows)}")
    print(f"changed_rows={len(changed_rows)} changed_icd10_rows={len(changed_icd10_rows)}")
    print(
        "decisions="
        + ",".join(f"{k}:{v}" for k, v in sorted(decision_counts.items(), key=lambda kv: kv[0]))
    )
    print(
        f"manual_better_icd10_rows={len(strict_manual_better)} "
        f"manual_better_unique_icd10_codes={len({r['icd10'] for r in strict_manual_better})}"
    )
    print(
        f"ambiguous_icd10_rows={len(ambiguous_icd10)} "
        f"ambiguous_unique_icd10_codes={len({r['icd10'] for r in ambiguous_icd10})}"
    )
    print(
        f"cache_preview_rows={len(cache_preview_rows)} cache_preview_unique_icd10_codes={len(unique_preview_codes)}"
    )
    if skipped_conflicts:
        print(f"cache_preview_skipped_conflicts={len(skipped_conflicts)} codes={','.join(sorted(skipped_conflicts))}")

    print(f"out_all={args.out_all}")
    print(f"out_manual_better={args.out_manual_better}")
    print(f"out_ambiguous={args.out_ambiguous}")
    print(f"out_cache_preview={args.out_cache_preview}")

    print(f"sha256_all={file_sha256(args.out_all)}")
    print(f"sha256_manual_better={file_sha256(args.out_manual_better)}")
    print(f"sha256_ambiguous={file_sha256(args.out_ambiguous)}")
    print(f"sha256_cache_preview={file_sha256(args.out_cache_preview)}")


if __name__ == "__main__":
    main()
