#!/usr/bin/env python3
"""Build uniprot_aliases.tsv from a UniProt export TSV/CSV.

Expected source export contains at least an accession-like column and optionally
protein/gene name columns. The script normalizes and deduplicates aliases.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

CANONICAL_OUTPUT_FIELDS = ["accession", "aliases", "gene", "source"]
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/stream"
USER_AGENT = "pqtl-measurement-mapper/1.0"

ACCESSION_KEYS = ["entry", "accession", "uniprotkb", "uniprot_id"]
GENE_KEYS = ["gene names (primary)", "gene names", "gene", "genes"]
PROTEIN_KEYS = ["protein names", "protein name", "protein"]
ALT_NAME_KEYS = ["alternative products", "protein names", "gene names"]


def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def normalize_acc(text: str) -> str:
    return normalize(text).upper()


def split_multi_values(text: str) -> list[str]:
    if not text:
        return []
    raw = normalize(text)
    if ";" in raw or "|" in raw:
        parts = re.split(r"[;|]", raw)
    else:
        parts = [raw]
    return [normalize(part) for part in parts if normalize(part)]


def split_gene_values(text: str) -> list[str]:
    out: list[str] = []
    for part in split_multi_values(text):
        # Gene fields commonly look like: "HP HPA HPT".
        if re.fullmatch(r"[A-Za-z0-9\- ]{2,120}", part) and " " in part:
            tokens = [t for t in part.split(" ") if t]
            if len(tokens) <= 16 and all(len(t) <= 24 for t in tokens):
                out.extend(tokens)
                continue
        out.append(part)
    return out


def choose_field(fieldnames: list[str], candidates: list[str]) -> str | None:
    lookup = {f.lower(): f for f in fieldnames}
    for key in candidates:
        if key in lookup:
            return lookup[key]
    return None


def infer_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    return ","


def build_uniprot_human_stream_url(reviewed_only: bool) -> str:
    query = "organism_id:9606"
    if reviewed_only:
        query += " AND reviewed:true"
    params = urllib.parse.urlencode(
        {
            "format": "tsv",
            "fields": "accession,gene_primary,gene_names,protein_name",
            "query": query,
        }
    )
    return f"{UNIPROT_BASE_URL}?{params}"


def download_uniprot_export(url: str, out_path: Path, timeout: float) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "text/tab-separated-values",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:  # nosec B310
        out_path.write_bytes(response.read())


def build_aliases(input_path: Path, output_path: Path, source_tag: str) -> tuple[int, int]:
    delimiter = infer_delimiter(input_path)

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError("Input file has no header")

        accession_col = choose_field(reader.fieldnames, ACCESSION_KEYS)
        gene_col = choose_field(reader.fieldnames, GENE_KEYS)
        protein_col = choose_field(reader.fieldnames, PROTEIN_KEYS)
        alt_col = choose_field(reader.fieldnames, ALT_NAME_KEYS)

        if not accession_col:
            raise ValueError(
                "Could not find accession column. Expected one of: "
                + ", ".join(ACCESSION_KEYS)
            )

        merged: dict[str, dict[str, set[str] | str]] = {}
        input_rows = 0

        for row in reader:
            input_rows += 1
            accession = normalize_acc(row.get(accession_col, ""))
            if not accession:
                continue

            gene_primary = ""
            aliases: set[str] = set()

            if gene_col:
                gene_values = split_gene_values(row.get(gene_col, ""))
                if gene_values:
                    gene_primary = gene_values[0]
                    aliases.update(gene_values)

            if protein_col:
                aliases.update(split_multi_values(row.get(protein_col, "")))

            if alt_col and alt_col != protein_col and alt_col != gene_col:
                aliases.update(split_multi_values(row.get(alt_col, "")))

            # Remove empty and self-duplicate aliases.
            aliases = {normalize(a) for a in aliases if normalize(a)}

            bucket = merged.setdefault(
                accession,
                {"aliases": set(), "gene": ""},
            )
            alias_bucket = bucket["aliases"]
            if isinstance(alias_bucket, set):
                alias_bucket.update(aliases)

            if gene_primary and not bucket["gene"]:
                bucket["gene"] = gene_primary

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        for accession in sorted(merged):
            row = merged[accession]
            aliases = sorted(row["aliases"]) if isinstance(row["aliases"], set) else []
            writer.writerow(
                {
                    "accession": accession,
                    "aliases": "|".join(aliases),
                    "gene": normalize(str(row["gene"])),
                    "source": source_tag,
                }
            )

    return input_rows, len(merged)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build uniprot_aliases.tsv from UniProt export")
    parser.add_argument("--input", help="UniProt export file (TSV/CSV)")
    parser.add_argument(
        "--download-human",
        action="store_true",
        help="Download UniProt human proteome export first, then build alias TSV",
    )
    parser.add_argument(
        "--download-to",
        default="/tmp/uniprot_human_export.tsv",
        help="Path for downloaded UniProt export when --download-human is used",
    )
    parser.add_argument(
        "--reviewed-only",
        action="store_true",
        help="When downloading, keep only reviewed (Swiss-Prot) human entries",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=120.0,
        help="HTTP timeout seconds for download mode",
    )
    parser.add_argument("--output", required=True, help="Output uniprot_aliases.tsv path")
    parser.add_argument("--source", default="uniprot-export", help="Source tag written to output")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        input_path: Path
        source_tag = args.source

        if args.download_human:
            input_path = Path(args.download_to)
            url = build_uniprot_human_stream_url(reviewed_only=args.reviewed_only)
            download_uniprot_export(url, input_path, timeout=args.timeout)
            print(f"[OK] downloaded UniProt export to {input_path}")
            if source_tag == "uniprot-export":
                source_tag = "uniprot-human-reviewed-live" if args.reviewed_only else "uniprot-human-live"
        else:
            if not args.input:
                raise ValueError("--input is required unless --download-human is set")
            input_path = Path(args.input)

        input_rows, output_rows = build_aliases(input_path, Path(args.output), source_tag)
        print(f"[OK] processed {input_rows} rows and wrote {output_rows} accessions to {args.output}")
        return 0
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
