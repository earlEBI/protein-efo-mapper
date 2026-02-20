#!/usr/bin/env python3
"""Build full EFO measurement branch cache from EFO OBO.

Default behavior:
- Download EFO OBO if local file is not supplied.
- Parse terms and parent relationships.
- Extract descendants of measurement root (EFO:0001444).
- Write TSV cache: efo_id, label, synonyms, source.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.request
from collections import deque
from pathlib import Path

DEFAULT_EFO_OBO_URL = "https://github.com/EBISPOT/efo/releases/latest/download/efo.obo"
DEFAULT_MEASUREMENT_ROOT = "EFO:0001444"
ALLOWED_ID_RE = re.compile(r"^(EFO|OBA):\d+$", re.IGNORECASE)


def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def normalize_obo_id(text: str) -> str:
    raw = normalize(text)
    if ":" in raw:
        prefix, rest = raw.split(":", 1)
        if prefix.lower() in {"efo", "oba"} and "_" in rest:
            base, number = rest.split("_", 1)
            if base.upper() in {"EFO", "OBA"}:
                return f"{base.upper()}:{number}"
    if "_" in raw:
        base, number = raw.split("_", 1)
        if base.upper() in {"EFO", "OBA"}:
            return f"{base.upper()}:{number}"
    return raw.upper()


def parse_obo_terms(obo_path: Path) -> dict[str, dict[str, object]]:
    terms: dict[str, dict[str, object]] = {}

    current: dict[str, object] | None = None
    with obo_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")

            if line == "[Term]":
                if current and current.get("id"):
                    term_id = str(current["id"])
                    terms[term_id] = current
                current = {"id": "", "name": "", "synonyms": [], "parents": [], "obsolete": False}
                continue

            if not current:
                continue

            if not line:
                continue

            if line.startswith("id: "):
                current["id"] = normalize_obo_id(line[4:])
                continue

            if line.startswith("name: "):
                current["name"] = normalize(line[6:])
                continue

            if line.startswith("synonym: "):
                # synonym: "text" EXACT []
                match = re.search(r'"([^"]+)"', line)
                if match:
                    syn = normalize(match.group(1))
                    if syn:
                        cast = current["synonyms"]
                        if isinstance(cast, list):
                            cast.append(syn)
                continue

            if line.startswith("is_a: "):
                parent = normalize_obo_id(line[6:].split(" ! ", 1)[0])
                if parent:
                    cast = current["parents"]
                    if isinstance(cast, list):
                        cast.append(parent)
                continue

            if line.startswith("is_obsolete: true"):
                current["obsolete"] = True
                continue

    if current and current.get("id"):
        terms[str(current["id"])] = current

    return terms


def measurement_descendants(terms: dict[str, dict[str, object]], root_id: str) -> set[str]:
    children: dict[str, set[str]] = {}
    for term_id, data in terms.items():
        parents = data.get("parents")
        if not isinstance(parents, list):
            continue
        for parent in parents:
            if not isinstance(parent, str):
                continue
            children.setdefault(parent, set()).add(term_id)

    seen: set[str] = set()
    queue: deque[str] = deque([root_id])

    while queue:
        current = queue.popleft()
        if current in seen:
            continue
        seen.add(current)
        for child in children.get(current, set()):
            if child not in seen:
                queue.append(child)

    return seen


def write_cache(terms: dict[str, dict[str, object]], selected_ids: set[str], output_path: Path, source: str) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["efo_id", "label", "synonyms", "source"], delimiter="\t")
        writer.writeheader()

        for term_id in sorted(selected_ids):
            if not ALLOWED_ID_RE.match(term_id):
                continue

            data = terms.get(term_id)
            if not data:
                continue
            if data.get("obsolete") is True:
                continue

            name = normalize(str(data.get("name") or ""))
            if not name:
                continue

            syns_raw = data.get("synonyms")
            syns: list[str] = []
            if isinstance(syns_raw, list):
                syns = sorted({normalize(s) for s in syns_raw if isinstance(s, str) and normalize(s)})

            writer.writerow(
                {
                    "efo_id": term_id,
                    "label": name,
                    "synonyms": "|".join(syns),
                    "source": source,
                }
            )
            written += 1

    return written


def download_obo(url: str, out_path: Path, timeout: float) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "pqtl-measurement-cache-builder/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # nosec B310
        out_path.write_bytes(resp.read())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build full measurement branch cache from EFO OBO")
    parser.add_argument(
        "--efo-obo",
        help="Path to local EFO OBO file. If omitted, script downloads from --download-url",
    )
    parser.add_argument(
        "--download-url",
        default=DEFAULT_EFO_OBO_URL,
        help="EFO OBO URL used when --efo-obo is not provided",
    )
    parser.add_argument(
        "--download-to",
        default="/tmp/efo.obo",
        help="Where to save downloaded OBO before parsing",
    )
    parser.add_argument(
        "--measurement-root",
        default=DEFAULT_MEASUREMENT_ROOT,
        help="Root ontology ID for measurement branch extraction",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path (efo_id,label,synonyms,source)",
    )
    parser.add_argument("--source", default="efo-obo-measurement-branch", help="Source label written to output")
    parser.add_argument("--timeout", type=float, default=60.0, help="Download timeout seconds")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.efo_obo:
            obo_path = Path(args.efo_obo)
        else:
            obo_path = Path(args.download_to)
            download_obo(args.download_url, obo_path, timeout=args.timeout)
            print(f"[OK] downloaded EFO OBO to {obo_path}")

        if not obo_path.exists():
            raise FileNotFoundError(f"EFO OBO not found: {obo_path}")

        terms = parse_obo_terms(obo_path)
        if not terms:
            raise ValueError("No terms parsed from OBO")

        root_id = normalize_obo_id(args.measurement_root)
        selected_ids = measurement_descendants(terms, root_id)
        if root_id not in selected_ids:
            raise ValueError(f"Measurement root not found in ontology graph: {root_id}")

        count = write_cache(terms, selected_ids, Path(args.output), args.source)
        print(f"[OK] wrote {count} measurement terms to {args.output} (root={root_id})")
        return 0
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
