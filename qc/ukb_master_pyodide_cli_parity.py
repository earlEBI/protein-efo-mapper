#!/usr/bin/env python3
"""Compare UKB master manual mappings vs CLI trait-map vs Pyodide mapper output."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import runpy
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any

ID_PATTERN = re.compile(r"\b(?:EFO|MONDO|HP|OBA)_[0-9]{3,}\b")
DEFAULT_OUTPUT_DIR = Path("tmp")


def normalize(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "").strip())


def parse_ids(value: str) -> set[str]:
    if not value:
        return set()
    cleaned = value.replace(",", "|")
    return set(ID_PATTERN.findall(cleaned.upper()))


def compare_id_sets(tool_ids: set[str], manual_ids: set[str]) -> str:
    if not manual_ids and not tool_ids:
        return "both_unmapped"
    if not manual_ids and tool_ids:
        return "tool_only"
    if manual_ids and not tool_ids:
        return "manual_only"
    if tool_ids & manual_ids:
        if tool_ids == manual_ids:
            return "exact_match"
        return "partial_overlap"
    return "different_nonoverlap"


def load_rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return [dict(row) for row in reader]


def run_cli_trait_map(
    python_bin: Path,
    mapper_script: Path,
    input_path: Path,
    output_path: Path,
    review_path: Path,
    qc_risk_path: Path,
    qc_summary_path: Path,
    trait_cache_path: Path,
    efo_obo_path: Path,
    ukb_field_catalog: Path,
    icd10_supplement_cache: Path,
    icd10_label_cache: Path,
) -> None:
    cmd = [
        str(python_bin),
        str(mapper_script),
        "trait-map",
        "--input",
        str(input_path),
        "--output",
        str(output_path),
        "--review-output",
        str(review_path),
        "--qc-risk-output",
        str(qc_risk_path),
        "--qc-summary-output",
        str(qc_summary_path),
        "--trait-cache",
        str(trait_cache_path),
        "--efo-obo",
        str(efo_obo_path),
        "--ukb-field-catalog",
        str(ukb_field_catalog),
        "--icd10-supplement-cache",
        str(icd10_supplement_cache),
        "--icd10-label-cache",
        str(icd10_label_cache),
        "--progress",
        "--stream-output",
        "--memoize-queries",
    ]
    subprocess.run(cmd, check=True)


def load_mapper_core(mapper_core_path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("pyodide_mapper_core", str(mapper_core_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load mapper core module from {mapper_core_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run_pyodide_fast_map(
    mapper_core_path: Path,
    cache_bundle_path: Path,
    manual_rows: list[dict[str, str]],
    min_score: float,
    pyodide_out: Path,
) -> None:
    module = load_mapper_core(mapper_core_path)
    cache_text = cache_bundle_path.read_text(encoding="utf-8")

    input_lines = ["query\ttrait_scale"]
    for row in manual_rows:
        query = normalize(row.get("ZOOMA QUERY", ""))
        if query:
            input_lines.append(f"{query}\tauto")
    input_tsv = "\n".join(input_lines) + "\n"

    mapped_tsv = module.map_tsv(input_tsv, cache_text, float(min_score))
    pyodide_out.write_text(mapped_tsv, encoding="utf-8")


def run_pyodide_cli_compat_map(
    runtime_bundle_path: Path,
    manual_input_path: Path,
    min_score: float,
    pyodide_out: Path,
) -> None:
    if not runtime_bundle_path.exists():
        raise FileNotFoundError(f"Runtime bundle not found: {runtime_bundle_path}")

    temp_root = Path(tempfile.mkdtemp(prefix="pyodide_runtime_parity_"))
    runtime_repo = temp_root / "repo"
    try:
        with zipfile.ZipFile(runtime_bundle_path, "r") as archive:
            archive.extractall(runtime_repo)

        script_path = runtime_repo / "skills" / "pqtl-measurement-mapper" / "scripts" / "map_measurement_efo.py"
        if not script_path.exists():
            raise FileNotFoundError(f"Mapper script missing in runtime bundle: {script_path}")

        review_path = pyodide_out.with_name(f"{pyodide_out.stem}_review.tsv")
        qc_risk_path = pyodide_out.with_name(f"{pyodide_out.stem}_qc_risk.tsv")
        qc_summary_path = pyodide_out.with_name(f"{pyodide_out.stem}_qc_summary.json")
        trait_cache = runtime_repo / "skills" / "pqtl-measurement-mapper" / "references" / "trait_mapping_cache.tsv"
        efo_obo = runtime_repo / "skills" / "pqtl-measurement-mapper" / "references" / "efo.obo"
        ukb_fieldsum = runtime_repo / "references" / "ukb" / "fieldsum.txt"
        icd10_supp = runtime_repo / "references" / "icd10" / "icd10_trait_supplement_cache.tsv"
        icd10_labels = runtime_repo / "references" / "icd10" / "icd10_label_index.tsv"

        argv = [
            str(script_path),
            "trait-map",
            "--input",
            str(manual_input_path),
            "--output",
            str(pyodide_out),
            "--review-output",
            str(review_path),
            "--qc-risk-output",
            str(qc_risk_path),
            "--qc-summary-output",
            str(qc_summary_path),
            "--trait-cache",
            str(trait_cache),
            "--efo-obo",
            str(efo_obo),
            "--ukb-field-catalog",
            str(ukb_fieldsum),
            "--icd10-supplement-cache",
            str(icd10_supp),
            "--icd10-label-cache",
            str(icd10_labels),
            "--min-score",
            str(float(min_score)),
            "--stream-output",
            "--memoize-queries",
            "--no-progress",
        ]
        old_argv = sys.argv[:]
        try:
            sys.argv = argv
            try:
                runpy.run_path(str(script_path), run_name="__main__")
            except SystemExit as exc:
                code = exc.code if isinstance(exc.code, int) else 0
                if code not in (0, None):
                    raise RuntimeError(f"CLI-compat run exited with code {code}") from exc
        finally:
            sys.argv = old_argv
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def map_rows_by_query(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        query = normalize(row.get(key, ""))
        if query and query not in out:
            out[query] = row
    return out


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ratio(count: int, total: int) -> float:
    return round((count / total) if total else 0.0, 6)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manual-input", type=Path, required=True, help="UK_Biobank_master_file.tsv path")
    parser.add_argument("--run-cli", action="store_true", help="Run CLI trait-map before comparison")
    parser.add_argument("--python-bin", type=Path, default=Path(".venv/bin/python"))
    parser.add_argument(
        "--mapper-script",
        type=Path,
        default=Path("skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py"),
    )
    parser.add_argument(
        "--trait-cache",
        type=Path,
        default=Path("skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv"),
    )
    parser.add_argument("--efo-obo", type=Path, default=Path("skills/pqtl-measurement-mapper/references/efo.obo"))
    parser.add_argument("--ukb-field-catalog", type=Path, default=Path("references/ukb/fieldsum.txt"))
    parser.add_argument(
        "--icd10-supplement-cache",
        type=Path,
        default=Path("references/icd10/icd10_trait_supplement_cache.tsv"),
    )
    parser.add_argument(
        "--icd10-label-cache",
        type=Path,
        default=Path("references/icd10/icd10_label_index.tsv"),
    )
    parser.add_argument(
        "--mapper-core",
        type=Path,
        default=Path("skills/pqtl-measurement-mapper/web/pyodide/mapper_core.py"),
    )
    parser.add_argument(
        "--cache-bundle",
        type=Path,
        default=Path("skills/pqtl-measurement-mapper/web/pyodide/cache_bundle.json"),
    )
    parser.add_argument(
        "--runtime-bundle",
        type=Path,
        default=Path("skills/pqtl-measurement-mapper/web/pyodide/trait_runtime_bundle.zip"),
    )
    parser.add_argument(
        "--pyodide-engine",
        choices=["fast", "cli_compat"],
        default="fast",
        help="Pyodide mapping engine: fast (mapper_core) or cli_compat (runtime bundle + map_measurement_efo.py)",
    )
    parser.add_argument("--min-score", type=float, default=0.82)
    parser.add_argument("--cli-output", type=Path, default=DEFAULT_OUTPUT_DIR / "ukb_master_cli_mapped.tsv")
    parser.add_argument("--cli-review-output", type=Path, default=DEFAULT_OUTPUT_DIR / "ukb_master_cli_review.tsv")
    parser.add_argument("--cli-qc-risk-output", type=Path, default=DEFAULT_OUTPUT_DIR / "ukb_master_cli_qc_risk.tsv")
    parser.add_argument("--cli-qc-summary-output", type=Path, default=DEFAULT_OUTPUT_DIR / "ukb_master_cli_qc_summary.json")
    parser.add_argument("--pyodide-output", type=Path, default=DEFAULT_OUTPUT_DIR / "ukb_master_pyodide_mapped.tsv")
    parser.add_argument(
        "--compare-output",
        type=Path,
        default=DEFAULT_OUTPUT_DIR / "ukb_master_cli_pyodide_manual_compare.tsv",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_OUTPUT_DIR / "ukb_master_cli_pyodide_manual_summary.json",
    )
    parser.add_argument(
        "--focus-output",
        type=Path,
        default=DEFAULT_OUTPUT_DIR / "ukb_master_cli_pyodide_manual_focus.tsv",
        help=(
            "Focused rows where CLI is manual_only/different_nonoverlap and "
            "Pyodide is exact_match/partial_overlap."
        ),
    )
    args = parser.parse_args()

    manual_rows = [row for row in load_rows(args.manual_input) if normalize(row.get("ZOOMA QUERY", ""))]
    if not manual_rows:
        raise SystemExit("No usable manual rows found")

    if args.run_cli:
        run_cli_trait_map(
            python_bin=args.python_bin,
            mapper_script=args.mapper_script,
            input_path=args.manual_input,
            output_path=args.cli_output,
            review_path=args.cli_review_output,
            qc_risk_path=args.cli_qc_risk_output,
            qc_summary_path=args.cli_qc_summary_output,
            trait_cache_path=args.trait_cache,
            efo_obo_path=args.efo_obo,
            ukb_field_catalog=args.ukb_field_catalog,
            icd10_supplement_cache=args.icd10_supplement_cache,
            icd10_label_cache=args.icd10_label_cache,
        )

    if not args.cli_output.exists():
        raise SystemExit(f"CLI output not found: {args.cli_output}")

    if args.pyodide_engine == "cli_compat":
        run_pyodide_cli_compat_map(
            runtime_bundle_path=args.runtime_bundle,
            manual_input_path=args.manual_input,
            min_score=args.min_score,
            pyodide_out=args.pyodide_output,
        )
    else:
        run_pyodide_fast_map(
            mapper_core_path=args.mapper_core,
            cache_bundle_path=args.cache_bundle,
            manual_rows=manual_rows,
            min_score=args.min_score,
            pyodide_out=args.pyodide_output,
        )

    cli_rows = load_rows(args.cli_output)
    pyodide_rows = load_rows(args.pyodide_output)
    cli_by_query = map_rows_by_query(cli_rows, "input_query")
    pyodide_by_query = map_rows_by_query(pyodide_rows, "input_query")

    compare_rows: list[dict[str, str]] = []
    focus_rows: list[dict[str, str]] = []
    cli_vs_manual_counts: Counter[str] = Counter()
    pyodide_vs_manual_counts: Counter[str] = Counter()
    pyodide_vs_cli_counts: Counter[str] = Counter()

    for row in manual_rows:
        query = normalize(row.get("ZOOMA QUERY", ""))
        manual_ids = parse_ids(normalize(row.get("MAPPED_TERM_URI", "")))
        cli = cli_by_query.get(query, {})
        py = pyodide_by_query.get(query, {})

        cli_ids = parse_ids(cli.get("mapped_trait_id", ""))
        pyodide_ids = parse_ids(py.get("mapped_trait_id", ""))

        cli_vs_manual = compare_id_sets(cli_ids, manual_ids)
        pyodide_vs_manual = compare_id_sets(pyodide_ids, manual_ids)
        pyodide_vs_cli = compare_id_sets(pyodide_ids, cli_ids)

        cli_vs_manual_counts[cli_vs_manual] += 1
        pyodide_vs_manual_counts[pyodide_vs_manual] += 1
        pyodide_vs_cli_counts[pyodide_vs_cli] += 1

        out_row = {
            "query": query,
            "manual_mapping_type": normalize(row.get("MAPPING_TYPE", "")),
            "manual_ids": "|".join(sorted(manual_ids)),
            "manual_label": normalize(row.get("MAPPED_TERM_LABEL", "")),
            "cli_ids": "|".join(sorted(cli_ids)),
            "cli_label": normalize(cli.get("mapped_trait_label", "")),
            "cli_confidence": normalize(cli.get("confidence", "")),
            "cli_matched_via": normalize(cli.get("matched_via", "")),
            "pyodide_ids": "|".join(sorted(pyodide_ids)),
            "pyodide_label": normalize(py.get("mapped_trait_label", "")),
            "pyodide_confidence": normalize(py.get("confidence", "")),
            "pyodide_matched_via": normalize(py.get("matched_via", "")),
            "cli_vs_manual": cli_vs_manual,
            "pyodide_vs_manual": pyodide_vs_manual,
            "pyodide_vs_cli": pyodide_vs_cli,
        }
        compare_rows.append(out_row)

        cli_bad = cli_vs_manual in {"manual_only", "different_nonoverlap"}
        pyodide_good = pyodide_vs_manual in {"exact_match", "partial_overlap"}
        if cli_bad and pyodide_good:
            focus_rows.append(out_row)

    fieldnames = [
        "query",
        "manual_mapping_type",
        "manual_ids",
        "manual_label",
        "cli_ids",
        "cli_label",
        "cli_confidence",
        "cli_matched_via",
        "pyodide_ids",
        "pyodide_label",
        "pyodide_confidence",
        "pyodide_matched_via",
        "cli_vs_manual",
        "pyodide_vs_manual",
        "pyodide_vs_cli",
    ]
    write_tsv(args.compare_output, compare_rows, fieldnames)
    write_tsv(args.focus_output, focus_rows, fieldnames)

    total = len(compare_rows)
    summary = {
        "rows": total,
        "pyodide_engine": args.pyodide_engine,
        "cli_vs_manual": dict(cli_vs_manual_counts),
        "pyodide_vs_manual": dict(pyodide_vs_manual_counts),
        "pyodide_vs_cli": dict(pyodide_vs_cli_counts),
        "rates": {
            "cli_manual_exact_or_partial": ratio(
                cli_vs_manual_counts.get("exact_match", 0) + cli_vs_manual_counts.get("partial_overlap", 0),
                total,
            ),
            "pyodide_manual_exact_or_partial": ratio(
                pyodide_vs_manual_counts.get("exact_match", 0) + pyodide_vs_manual_counts.get("partial_overlap", 0),
                total,
            ),
            "cli_manual_only_or_nonoverlap": ratio(
                cli_vs_manual_counts.get("manual_only", 0) + cli_vs_manual_counts.get("different_nonoverlap", 0),
                total,
            ),
            "pyodide_manual_only_or_nonoverlap": ratio(
                pyodide_vs_manual_counts.get("manual_only", 0) + pyodide_vs_manual_counts.get("different_nonoverlap", 0),
                total,
            ),
            "focus_rows_share": ratio(len(focus_rows), total),
        },
        "focus_rows": len(focus_rows),
        "files": {
            "manual_input": str(args.manual_input),
            "cli_output": str(args.cli_output),
            "pyodide_output": str(args.pyodide_output),
            "runtime_bundle": str(args.runtime_bundle),
            "compare_output": str(args.compare_output),
            "focus_output": str(args.focus_output),
        },
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
