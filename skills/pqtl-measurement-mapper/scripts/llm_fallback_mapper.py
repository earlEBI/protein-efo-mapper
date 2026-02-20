#!/usr/bin/env python3
"""LLM-assisted fallback for unresolved analyte->measurement term mappings.

This script is provider-agnostic for any OpenAI-compatible chat completion API.
It reads an existing mapping TSV, finds `validation=not_mapped` rows, builds a
strict candidate shortlist from the local index, asks an LLM to choose one
candidate ID (or NONE), and applies strict post-validation before accepting.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import map_measurement_efo as core

ID_RE = re.compile(r"\b(?:EFO|OBA):\d+\b", re.IGNORECASE)


@dataclass
class ShortCandidate:
    efo_id: str
    label: str
    score: float
    synonyms: list[str]
    evidence: str


def normalize(text: str) -> str:
    return core.normalize(text)


def build_shortlist(
    query: str,
    index: dict,
    matrix_priority: list[str],
    measurement_context: str,
    max_candidates: int,
    min_score: float,
) -> list[ShortCandidate]:
    term_meta = index.get("term_meta", {})
    exact_term_index = index.get("exact_term_index", {})
    token_index = index.get("token_index", {})

    by_id: dict[str, ShortCandidate] = {}
    terms = core.build_query_aliases(query, index)

    def upsert(candidate: ShortCandidate) -> None:
        cur = by_id.get(candidate.efo_id)
        if cur is None or candidate.score > cur.score:
            by_id[candidate.efo_id] = candidate

    for term in terms:
        key = core.norm_key(term)
        for efo_id in exact_term_index.get(key, []):
            meta = term_meta.get(efo_id, {})
            label = normalize(meta.get("label") or "")
            syns = [normalize(x) for x in (meta.get("synonyms") or []) if normalize(x)]
            if not label:
                continue
            if not core.accession_identity_match(query, label, syns, index):
                continue
            if not core.context_compatible(label, syns, measurement_context):
                continue
            score = max(0.9, core.similarity_score(term, label, syns))
            score = min(1.0, max(0.0, score + core.matrix_bonus(label, syns, matrix_priority)))
            upsert(
                ShortCandidate(
                    efo_id=efo_id,
                    label=label,
                    score=score,
                    synonyms=syns,
                    evidence=f"exact: {term}",
                )
            )

        ids: set[str] = set()
        for tok in core.tokenize(term):
            if tok in core.GENERIC_TOKENS:
                continue
            ids.update(token_index.get(tok, []))

        for efo_id in ids:
            meta = term_meta.get(efo_id, {})
            label = normalize(meta.get("label") or "")
            syns = [normalize(x) for x in (meta.get("synonyms") or []) if normalize(x)]
            if not label:
                continue
            if not core.accession_identity_match(query, label, syns, index):
                continue
            if not core.context_compatible(label, syns, measurement_context):
                continue
            score = core.similarity_score(term, label, syns)
            if not core.is_measurement_like(label, syns):
                score -= 0.35
            score += core.matrix_bonus(label, syns, matrix_priority)
            score = min(1.0, max(0.0, score))
            if score < min_score:
                continue
            upsert(
                ShortCandidate(
                    efo_id=efo_id,
                    label=label,
                    score=score,
                    synonyms=syns,
                    evidence=f"token: {term}",
                )
            )

    ranked = sorted(
        by_id.values(),
        key=lambda c: (c.score, core.matrix_preference_rank(c.label, matrix_priority), c.efo_id),
        reverse=True,
    )
    return ranked[:max_candidates]


def chat_complete(
    api_base: str,
    api_key: str,
    model: str,
    system_prompt: str,
    user_prompt: str,
    timeout: float,
) -> str:
    url = api_base.rstrip("/") + "/chat/completions"
    payload = {
        "model": model,
        "temperature": 0,
        "max_tokens": 32,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
    }
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # nosec B310
        body = json.loads(resp.read().decode("utf-8"))

    choices = body.get("choices") or []
    if not choices:
        return ""
    msg = (choices[0] or {}).get("message") or {}
    content = msg.get("content")
    if isinstance(content, str):
        return content.strip()
    if isinstance(content, list):
        chunks: list[str] = []
        for part in content:
            if isinstance(part, dict):
                text = part.get("text")
                if isinstance(text, str):
                    chunks.append(text)
        return "\n".join(chunks).strip()
    return ""


def extract_choice_id(text: str) -> str:
    if not text:
        return ""
    if "NONE" in text.upper():
        return "NONE"
    m = ID_RE.search(text)
    if not m:
        return ""
    return normalize(m.group(0)).upper().replace("_", ":")


def llm_choose_candidate(
    query: str,
    candidates: list[ShortCandidate],
    api_base: str,
    api_key: str,
    model: str,
    timeout: float,
) -> str:
    system = (
        "Choose exactly one ontology ID from candidates or NONE. "
        "Do not invent IDs. Respect exact analyte identity and numeric family."
    )

    lines = [f"Query: {query}", "Candidates:"]
    for idx, c in enumerate(candidates, start=1):
        syn_preview = "; ".join(c.synonyms[:2])
        lines.append(f"{idx}. {c.efo_id}\t{c.label}\tscore={c.score:.3f}\tsyn={syn_preview}")
    lines.append("Return only: EFO:... or OBA:... or NONE")

    content = chat_complete(
        api_base=api_base,
        api_key=api_key,
        model=model,
        system_prompt=system,
        user_prompt="\n".join(lines),
        timeout=timeout,
    )
    return extract_choice_id(content)


def process_row(
    row: dict[str, str],
    index: dict,
    matrix_priority: list[str],
    measurement_context: str,
    max_candidates: int,
    shortlist_min_score: float,
    dry_run: bool,
    api_base: str,
    api_key: str,
    model: str,
    timeout: float,
) -> dict[str, str]:
    if row.get("validation") != "not_mapped":
        return row

    query = normalize(row.get("input_query") or "")
    if not query:
        return row

    shortlist = build_shortlist(
        query=query,
        index=index,
        matrix_priority=matrix_priority,
        measurement_context=measurement_context,
        max_candidates=max_candidates,
        min_score=shortlist_min_score,
    )
    if not shortlist:
        row["evidence"] = "llm-fallback: no shortlist candidates"
        return row

    if dry_run:
        row["evidence"] = f"llm-fallback dry-run: shortlist={len(shortlist)}"
        return row

    choice = llm_choose_candidate(
        query=query,
        candidates=shortlist,
        api_base=api_base,
        api_key=api_key,
        model=model,
        timeout=timeout,
    )
    if not choice or choice == "NONE":
        row["evidence"] = "llm-fallback: model returned NONE"
        return row

    chosen = None
    for c in shortlist:
        if c.efo_id.upper() == choice.upper():
            chosen = c
            break
    if chosen is None:
        row["evidence"] = "llm-fallback: model chose ID outside shortlist"
        return row

    # Strict post-validation gate.
    if not core.accession_identity_match(query, chosen.label, chosen.synonyms, index):
        row["evidence"] = "llm-fallback: rejected by accession identity validation"
        return row
    if not core.context_compatible(chosen.label, chosen.synonyms, measurement_context):
        row["evidence"] = "llm-fallback: rejected by matrix/tissue context validation"
        return row

    row["mapped_efo_id"] = chosen.efo_id
    row["mapped_label"] = chosen.label
    row["confidence"] = f"{chosen.score:.3f}"
    row["matched_via"] = "llm-fallback"
    row["validation"] = "validated"
    row["evidence"] = f"llm-fallback selected from shortlist ({len(shortlist)}): {chosen.evidence}"
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LLM-assisted fallback for unresolved mappings")
    parser.add_argument("--mapped-input", required=True, help="Existing mapping TSV (from map_measurement_efo.py)")
    parser.add_argument("--output", required=True, help="Output TSV with LLM fallback applied")
    parser.add_argument("--index", default=str(core.DEFAULT_INDEX), help="Measurement JSON index")
    parser.add_argument(
        "--matrix-priority",
        default=",".join(core.DEFAULT_MATRIX_PRIORITY),
        help="Preferred sample matrices in order, comma-separated",
    )
    parser.add_argument(
        "--measurement-context",
        default=core.DEFAULT_MEASUREMENT_CONTEXT,
        help="Expected measurement context: blood, cerebrospinal_fluid, urine, saliva, tissue, or auto",
    )
    parser.add_argument("--max-candidates", type=int, default=15, help="Shortlist size passed to LLM")
    parser.add_argument("--shortlist-min-score", type=float, default=0.20, help="Min score for shortlist candidates")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers for LLM fallback rows")

    parser.add_argument("--api-base", default="https://api.openai.com/v1", help="OpenAI-compatible API base URL")
    parser.add_argument("--model", default="gpt-4.1-mini", help="Model name for fallback")
    parser.add_argument("--api-key-env", default="OPENAI_API_KEY", help="Env var that contains API key")
    parser.add_argument("--timeout", type=float, default=45.0, help="HTTP timeout seconds")
    parser.add_argument("--dry-run", action="store_true", help="Do not call LLM; only annotate shortlist availability")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        index = core.load_index(Path(args.index))
        matrix_priority = core.parse_matrix_priority(args.matrix_priority)
        measurement_context = core.parse_measurement_context(args.measurement_context)

        with Path(args.mapped_input).open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            fieldnames = reader.fieldnames

        if not fieldnames:
            raise ValueError("Mapped input has no header")

        api_key = ""
        if not args.dry_run:
            api_key = os.environ.get(args.api_key_env, "").strip()
            if not api_key:
                raise ValueError(f"Missing API key in env var: {args.api_key_env}")

        def runner(row: dict[str, str]) -> dict[str, str]:
            return process_row(
                row=dict(row),
                index=index,
                matrix_priority=matrix_priority,
                measurement_context=measurement_context,
                max_candidates=max(1, args.max_candidates),
                shortlist_min_score=args.shortlist_min_score,
                dry_run=args.dry_run,
                api_base=args.api_base,
                api_key=api_key,
                model=args.model,
                timeout=args.timeout,
            )

        if args.workers <= 1:
            out_rows = [runner(r) for r in rows]
        else:
            with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
                out_rows = list(pool.map(runner, rows))

        unresolved_before = sum(1 for r in rows if r.get("validation") == "not_mapped")
        unresolved_after = sum(1 for r in out_rows if r.get("validation") == "not_mapped")
        accepted = sum(1 for r in out_rows if r.get("matched_via") == "llm-fallback")

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(out_rows)

        print(
            f"[OK] wrote {len(out_rows)} rows to {args.output}; "
            f"unresolved_before={unresolved_before}, unresolved_after={unresolved_after}, "
            f"llm_accepted={accepted}"
        )
        return 0
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
