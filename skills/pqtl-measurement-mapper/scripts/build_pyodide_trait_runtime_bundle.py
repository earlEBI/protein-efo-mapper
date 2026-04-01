#!/usr/bin/env python3
"""Build a CLI-compatible trait-map runtime bundle for Pyodide execution."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import zipfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class BundleFile:
    relative_path: str
    required: bool = True


RUNTIME_FILES: tuple[BundleFile, ...] = (
    BundleFile("skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py"),
    BundleFile("skills/pqtl-measurement-mapper/references/trait_mapping_cache.tsv"),
    BundleFile("skills/pqtl-measurement-mapper/references/efo.obo"),
    BundleFile("skills/pqtl-measurement-mapper/references/manual_curator_seed_cache.tsv"),
    BundleFile("skills/pqtl-measurement-mapper/references/manual_curator_rescue_overrides.tsv"),
    BundleFile("references/ukb/fieldsum.txt"),
    BundleFile("references/ukb/field.txt"),
    BundleFile("references/ukb/category.txt"),
    BundleFile("references/ukb/catbrowse.txt"),
    BundleFile("references/ukb/field_trait_supplement_cache.tsv"),
    BundleFile("references/icd10/icd10_trait_supplement_cache.tsv"),
    BundleFile("references/icd10/icd10_label_index.tsv"),
    BundleFile("references/icd10/icd10_label_alias_index.tsv", required=False),
    BundleFile("references/icd10/mondo.sssom.tsv", required=False),
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def normalize_relative(path: str) -> str:
    return path.replace("\\", "/").lstrip("/")


def build_runtime_bundle(repo_root: Path, output_zip: Path, compresslevel: int = 6) -> dict[str, object]:
    selected_paths: list[Path] = []
    manifest_files: list[dict[str, object]] = []
    missing_required: list[str] = []

    for entry in RUNTIME_FILES:
        rel = normalize_relative(entry.relative_path)
        abs_path = (repo_root / rel).resolve()
        if not abs_path.exists():
            if entry.required:
                missing_required.append(rel)
            continue
        if not abs_path.is_file():
            if entry.required:
                missing_required.append(rel)
            continue
        selected_paths.append(abs_path)
        manifest_files.append(
            {
                "path": rel,
                "size_bytes": abs_path.stat().st_size,
                "sha256": file_sha256(abs_path),
                "required": entry.required,
            }
        )

    if missing_required:
        missing_text = ", ".join(sorted(missing_required))
        raise FileNotFoundError(f"Missing required runtime files: {missing_text}")

    output_zip.parent.mkdir(parents=True, exist_ok=True)
    compression = zipfile.ZIP_DEFLATED
    with zipfile.ZipFile(
        output_zip,
        mode="w",
        compression=compression,
        compresslevel=max(0, min(9, int(compresslevel))),
    ) as archive:
        for abs_path in selected_paths:
            rel = abs_path.relative_to(repo_root).as_posix()
            archive.write(abs_path, arcname=rel)

        manifest = {
            "version": "pyodide-trait-runtime-v1",
            "repo_root": str(repo_root),
            "file_count": len(selected_paths),
            "files": manifest_files,
        }
        archive.writestr(
            "runtime_manifest.json",
            json.dumps(manifest, ensure_ascii=False, indent=2),
        )

    archive_size = output_zip.stat().st_size
    input_size = sum(path.stat().st_size for path in selected_paths)
    ratio = (archive_size / input_size) if input_size else 0.0

    return {
        "output_zip": str(output_zip),
        "file_count": len(selected_paths),
        "input_size_bytes": input_size,
        "archive_size_bytes": archive_size,
        "compression_ratio": round(ratio, 4),
    }


def main() -> None:
    default_repo_root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=str(default_repo_root),
        help="Repository root containing skills/ and references/ directories",
    )
    parser.add_argument(
        "--output",
        default="skills/pqtl-measurement-mapper/web/pyodide/trait_runtime_bundle.zip",
        help="Output runtime ZIP path (relative to --repo-root if not absolute)",
    )
    parser.add_argument(
        "--compresslevel",
        type=int,
        default=6,
        help="ZIP deflate compression level (0-9, default: 6)",
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    output = Path(args.output)
    output_zip = output if output.is_absolute() else (repo_root / output).resolve()

    summary = build_runtime_bundle(
        repo_root=repo_root,
        output_zip=output_zip,
        compresslevel=int(args.compresslevel),
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
