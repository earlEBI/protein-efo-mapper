#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "traits-4k.csv"
OUT_ROOT = ROOT / "final_output" / "traits_4k_qc"
BATCH_SIZE = 50
PYTHON = ROOT / ".venv" / "bin" / "python"
MAPPER = ROOT / "skills" / "pqtl-measurement-mapper" / "scripts" / "map_measurement_efo.py"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        rows = []
        for row in r:
            fixed = {}
            for k, v in row.items():
                nk = (k or "").replace("\ufeff", "").strip()
                fixed[nk] = (v or "").strip()
            rows.append(fixed)
        return rows


def write_input_chunk(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["query", "trait_scale", "code", "data_type", "further_info"]
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def run_mapper(inp: Path, out: Path, review: Path, qc_risk: Path, qc_sum: Path) -> None:
    cmd = [
        str(PYTHON),
        str(MAPPER),
        "trait-map",
        "--input",
        str(inp),
        "--output",
        str(out),
        "--review-output",
        str(review),
        "--qc-risk-output",
        str(qc_risk),
        "--qc-summary-output",
        str(qc_sum),
        "--top-k",
        "1",
        "--min-score",
        "0.82",
    ]
    subprocess.run(cmd, check=True, cwd=str(ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def query_stem(query: str) -> str:
    q = (query or "").strip().lower()
    if q.startswith("comparative body size at age 10") or q.startswith("comparative height size at age 10"):
        return q.split(" - ", 1)[0]
    return q


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows = read_rows(INPUT)
    total = len(rows)

    all_output: list[dict[str, str]] = []
    all_flags: list[dict[str, str]] = []
    batch_summaries: list[dict[str, str]] = []

    for i in range(0, total, BATCH_SIZE):
        chunk = rows[i : i + BATCH_SIZE]
        batch_idx = i // BATCH_SIZE + 1
        bdir = OUT_ROOT / f"batch_{batch_idx:03d}"
        inp = bdir / "input.tsv"
        out = bdir / "mapped.tsv"
        review = bdir / "review.tsv"
        qc_risk = bdir / "qc_risk.tsv"
        qc_sum = bdir / "qc_summary.json"

        write_input_chunk(chunk, inp)
        run_mapper(inp, out, review, qc_risk, qc_sum)

        mapped = read_tsv(out)
        all_output.extend(mapped)

        val_counts = Counter(r.get("validation", "") for r in mapped)
        via_counts = Counter(r.get("matched_via", "") for r in mapped)

        flags = []
        stem_to_ids: dict[str, set[str]] = defaultdict(set)
        for r in mapped:
            stem_to_ids[query_stem(r.get("input_query", ""))].add(r.get("mapped_trait_id", ""))

            mapped_id = (r.get("mapped_trait_id", "") or "").strip()
            label = (r.get("mapped_trait_label", "") or "").strip()
            validation = (r.get("validation", "") or "").strip()
            via = (r.get("matched_via", "") or "").strip()

            reason = []
            if not mapped_id:
                reason.append("empty_mapping")
            if validation != "validated":
                reason.append(f"validation={validation or 'missing'}")
            if label in {"SR", "normal", "measurement", "disease"}:
                reason.append("generic_or_acronym_label")
            if via in {"efo_obo_fuzzy", "cache_fuzzy_text"}:
                reason.append(f"fuzzy={via}")

            if reason:
                fr = dict(r)
                fr["qc_reason"] = ";".join(reason)
                fr["batch"] = str(batch_idx)
                flags.append(fr)

        for stem, ids in stem_to_ids.items():
            clean_ids = {x for x in ids if x}
            if len(clean_ids) > 1:
                for r in mapped:
                    if query_stem(r.get("input_query", "")) == stem:
                        fr = dict(r)
                        fr["qc_reason"] = f"inconsistent_stem_ids={','.join(sorted(clean_ids))}"
                        fr["batch"] = str(batch_idx)
                        flags.append(fr)

        all_flags.extend(flags)
        batch_summaries.append(
            {
                "batch": str(batch_idx),
                "start_row": str(i + 1),
                "end_row": str(i + len(chunk)),
                "rows": str(len(chunk)),
                "validated": str(val_counts.get("validated", 0)),
                "review_required": str(val_counts.get("review_required", 0)),
                "not_mapped": str(val_counts.get("not_mapped", 0)),
                "flags": str(len(flags)),
                "top_via": ",".join(f"{k}:{v}" for k, v in via_counts.most_common(5)),
            }
        )
        print(f"batch {batch_idx:03d} rows {i+1}-{i+len(chunk)} flags={len(flags)}")

    if all_output:
        out_all = OUT_ROOT / "traits_4k_mapped_all.tsv"
        with out_all.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_output[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(all_output)

    if all_flags:
        flags_path = OUT_ROOT / "traits_4k_qc_flags.tsv"
        fields = list(all_flags[0].keys())
        with flags_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(all_flags)

    summary_path = OUT_ROOT / "traits_4k_batch_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(batch_summaries[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(batch_summaries)

    print("done", total, "rows")
    print("output:", OUT_ROOT)


if __name__ == "__main__":
    main()
