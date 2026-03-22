#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_BASELINE = Path("/Users/earl/Downloads/PMID34662886_studies_export.tsv")
DEFAULT_CURRENT = Path("/tmp/PMID34662886_reportedTrait_mapped_after_seed_lock_full_20260321c5f.tsv")
DEFAULT_OUT_FULL = Path("/tmp/PMID34662886_drug_use_qc_full_20260321.tsv")
DEFAULT_OUT_PRIORITY = Path("/tmp/PMID34662886_drug_use_qc_priority_20260321.tsv")
DEFAULT_OUT_SUMMARY = Path("/tmp/PMID34662886_drug_use_qc_summary_20260321.tsv")

MEDICATION_CONTEXT_RE = re.compile(
    r"\b("
    r"treatment or medication use|treatment or medication code|medication use|drug use|"
    r"medication for|metformin|insulin|supplement|analgesic|opioid|antihypertensive|"
    r"cholesterol lowering medication|laxative|nsaid|ssri"
    r")\b",
    re.IGNORECASE,
)
GENERIC_DRUG_LABELS = {
    "drug use measurement",
    "treatment",
}
ALLOWED_NON_MEDICATION_LABEL_HINTS = (
    "drug allergy",
    "adverse effect",
    "response to xenobiotic",
    "medication adherence",
    "encounter with health service",
    "poison",
)
MEDICATION_SPECIFIC_LABEL_HINTS = (
    "use measurement",
    "medication",
    "drug",
    "supplement",
    "antihypertensive",
    "antithrombotic",
    "opioid",
    "ssri",
    "hmg coa reductase",
    "calcium channel blocker",
    "diuretic",
    "beta blocking",
)

TREATMENT_MEDICATION_TOKEN_RE = re.compile(r"treatment or medication use\s*-\s*([^()]+)", re.IGNORECASE)


def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", (text or "").strip())


def norm_key(text: str) -> str:
    return normalize(text).lower()


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


def extract_medication_token(query: str) -> str:
    text = normalize(query)
    if not text:
        return ""
    match = TREATMENT_MEDICATION_TOKEN_RE.search(text)
    if not match:
        return ""
    token = normalize(match.group(1))
    token = re.sub(r"\s+", " ", token).strip(" -")
    return token


def is_medication_context(*parts: str) -> bool:
    text = " ".join(part for part in parts if normalize(part))
    if not text:
        return False
    return bool(MEDICATION_CONTEXT_RE.search(text))


def mapped_label_looks_medication(label: str) -> bool:
    key = norm_key(label)
    if not key:
        return False
    if any(hint in key for hint in MEDICATION_SPECIFIC_LABEL_HINTS):
        return True
    return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Medication/drug-use QC helper for PMID34662886.")
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--current", type=Path, default=DEFAULT_CURRENT)
    parser.add_argument("--out-full", type=Path, default=DEFAULT_OUT_FULL)
    parser.add_argument("--out-priority", type=Path, default=DEFAULT_OUT_PRIORITY)
    parser.add_argument("--out-summary", type=Path, default=DEFAULT_OUT_SUMMARY)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    baseline_rows = read_tsv(args.baseline)
    current_rows = read_tsv(args.current)
    if len(baseline_rows) != len(current_rows):
        raise ValueError(f"Row count mismatch: baseline={len(baseline_rows)} current={len(current_rows)}")

    rows: list[dict[str, str]] = []
    token_to_labels: dict[str, set[str]] = defaultdict(set)
    token_to_rownums: dict[str, list[int]] = defaultdict(list)

    for rownum, (base, cur) in enumerate(zip(baseline_rows, current_rows), start=1):
        reported_trait = base.get("reportedTrait", "")
        manual_efo = base.get("efoTraits", "")
        input_query = cur.get("input_query", "")
        mapped_id = cur.get("mapped_trait_id", "")
        mapped_label = cur.get("mapped_trait_label", "")
        matched_via = cur.get("matched_via", "")
        validation = cur.get("validation", "")

        medication_context = is_medication_context(reported_trait, input_query, manual_efo)
        mapped_looks_medication = mapped_label_looks_medication(mapped_label)
        manual_looks_medication = is_medication_context(manual_efo)

        if not (medication_context or mapped_looks_medication or manual_looks_medication):
            continue

        token = extract_medication_token(input_query or reported_trait)
        if token:
            token_to_labels[token].add(norm_key(mapped_label))
            token_to_rownums[token].append(rownum)

        flags: list[str] = []
        mapped_label_key = norm_key(mapped_label)
        manual_label_key = norm_key(manual_efo)
        matched_via_key = norm_key(matched_via)
        query_key = norm_key(input_query or reported_trait)

        if (
            medication_context
            and mapped_label_key in GENERIC_DRUG_LABELS
            and manual_label_key
            and manual_label_key not in GENERIC_DRUG_LABELS
        ):
            flags.append("generic_current_vs_manual")
        if medication_context and matched_via_key in {"efo_obo_fuzzy", "cache_fuzzy_text"}:
            flags.append("fuzzy_drug_mapping")
        if medication_context and not mapped_looks_medication:
            if not any(hint in mapped_label_key for hint in ALLOWED_NON_MEDICATION_LABEL_HINTS):
                flags.append("non_medication_target_for_medication_query")
        if "ukb data field 6177" in query_key and "insulin" in query_key and mapped_id != "EFO_0009924":
            flags.append("ukb_6177_not_diabetes_drug_use_measurement")
        if "started insulin within one year" in query_key and mapped_id != "EFO_0009924":
            flags.append("started_insulin_not_diabetes_drug_use_measurement")

        rows.append(
            {
                "rownum": str(rownum),
                "reportedTrait": reported_trait,
                "input_query": input_query,
                "manual_efoTraits": manual_efo,
                "mapped_trait_id": mapped_id,
                "mapped_trait_label": mapped_label,
                "matched_via": matched_via,
                "validation": validation,
                "medication_token": token,
                "flags": "|".join(flags),
            }
        )

    # Post-pass inconsistency check per explicit medication token
    inconsistent_tokens = {
        token
        for token, labels in token_to_labels.items()
        if token and len({label for label in labels if label}) > 1
    }
    rows_by_rownum = {int(row["rownum"]): row for row in rows}
    for token in inconsistent_tokens:
        for rownum in token_to_rownums.get(token, []):
            row = rows_by_rownum.get(rownum)
            if row is None:
                continue
            existing = [flag for flag in row.get("flags", "").split("|") if flag]
            if "inconsistent_mapping_for_same_medication_token" not in existing:
                existing.append("inconsistent_mapping_for_same_medication_token")
            row["flags"] = "|".join(existing)

    priority_rows = [row for row in rows if row.get("flags")]

    fieldnames = [
        "rownum",
        "reportedTrait",
        "input_query",
        "manual_efoTraits",
        "mapped_trait_id",
        "mapped_trait_label",
        "matched_via",
        "validation",
        "medication_token",
        "flags",
    ]
    write_tsv(args.out_full, rows, fieldnames)
    write_tsv(args.out_priority, priority_rows, fieldnames)

    flag_counter: Counter[str] = Counter()
    mapped_counter: Counter[str] = Counter()
    for row in rows:
        mapped_counter[row["mapped_trait_label"]] += 1
        for flag in [f for f in row.get("flags", "").split("|") if f]:
            flag_counter[flag] += 1

    summary_rows: list[dict[str, str]] = []
    summary_rows.append({"metric": "drug_rows_total", "value": str(len(rows))})
    summary_rows.append({"metric": "priority_rows_total", "value": str(len(priority_rows))})
    for label, count in mapped_counter.most_common(25):
        summary_rows.append({"metric": f"mapped_label::{label}", "value": str(count)})
    for flag, count in flag_counter.most_common():
        summary_rows.append({"metric": f"flag::{flag}", "value": str(count)})
    write_tsv(args.out_summary, summary_rows, ["metric", "value"])

    print("Drug-use QC complete")
    print(f"rows_total={len(rows)} priority_rows={len(priority_rows)}")
    print(f"inconsistent_tokens={len(inconsistent_tokens)}")
    print(f"out_full={args.out_full}")
    print(f"out_priority={args.out_priority}")
    print(f"out_summary={args.out_summary}")


if __name__ == "__main__":
    main()
