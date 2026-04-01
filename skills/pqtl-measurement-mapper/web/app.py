#!/usr/bin/env python3
"""Minimal web app for pQTL measurement mapping.

Features:
- Upload analyte file (.txt/.tsv/.csv)
- Run mapper in background job
- Live progress/log polling
- Download mapped + unmapped TSV outputs
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse

ROOT_DIR = Path(__file__).resolve().parents[3]
SKILL_DIR = ROOT_DIR / "skills" / "pqtl-measurement-mapper"
MAPPER_SCRIPT = SKILL_DIR / "scripts" / "map_measurement_efo.py"
DEFAULT_INDEX = SKILL_DIR / "references" / "measurement_index.json"
DEFAULT_UNIPROT_ALIASES = SKILL_DIR / "references" / "uniprot_aliases.tsv"
DEFAULT_TERM_CACHE = SKILL_DIR / "references" / "efo_measurement_terms_cache.tsv"
DEFAULT_ANALYTE_CACHE = SKILL_DIR / "references" / "analyte_to_efo_cache.tsv"
DEFAULT_TRAIT_CACHE = SKILL_DIR / "references" / "trait_mapping_cache.tsv"
DEFAULT_EFO_OBO = SKILL_DIR / "references" / "efo.obo"
DEFAULT_UKB_FIELD_CATALOG = ROOT_DIR / "references" / "ukb" / "fieldsum.txt"
PYODIDE_WEB_DIR = SKILL_DIR / "web" / "pyodide"
PYTHON_BIN = ROOT_DIR / ".venv" / "bin" / "python"
JOB_ROOT = ROOT_DIR / "final_output" / "web_jobs"
EXPORT_ROOT = ROOT_DIR / "final_output" / "web_exports"
PROGRESS_RE = re.compile(r"\[progress\]\s+(\d+)/(\d+)\s+\((\d+\.?\d*)%\)")


@dataclass
class JobState:
    job_id: str
    status: str
    created_at: float
    updated_at: float
    progress_current: int = 0
    progress_total: int = 0
    progress_percent: float = 0.0
    message: str = "queued"
    error: str = ""
    logs: list[str] = field(default_factory=list)
    input_path: str = ""
    mapped_path: str = ""
    unmapped_path: str = ""
    review_path: str = ""
    qc_risk_path: str = ""
    qc_summary_path: str = ""
    local_export_dir: str = ""
    options: dict[str, Any] = field(default_factory=dict)


app = FastAPI(title="pQTL Measurement Mapper")
JOBS: dict[str, JobState] = {}
JOBS_LOCK = threading.Lock()


def add_log(job: JobState, line: str) -> None:
    line = line.rstrip("\n")
    if not line:
        return
    job.logs.append(line)
    if len(job.logs) > 500:
        job.logs = job.logs[-500:]


def set_progress_from_line(job: JobState, line: str) -> None:
    m = PROGRESS_RE.search(line)
    if not m:
        return
    cur = int(m.group(1))
    total = int(m.group(2))
    pct = float(m.group(3))
    job.progress_current = cur
    job.progress_total = total
    job.progress_percent = pct


def get_job(job_id: str) -> JobState:
    with JOBS_LOCK:
        job = JOBS.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


def write_csv_copy(tsv_path: Path) -> Path | None:
    if not tsv_path.exists() or tsv_path.suffix.lower() != ".tsv":
        return None
    csv_path = tsv_path.with_suffix(".csv")
    with tsv_path.open("r", encoding="utf-8", newline="") as src, csv_path.open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst)
        for row in reader:
            writer.writerow(row)
    return csv_path


def stage_job_exports(job: JobState) -> None:
    export_dir = Path(job.local_export_dir)
    export_dir.mkdir(parents=True, exist_ok=True)
    artifacts = [
        Path(job.input_path),
        Path(job.mapped_path),
        Path(job.unmapped_path),
        Path(job.review_path),
        Path(job.qc_risk_path),
        Path(job.qc_summary_path),
    ]
    for path in artifacts:
        if path.exists():
            shutil.copy2(path, export_dir / path.name)
            csv_path = write_csv_copy(path)
            if csv_path is not None and csv_path.exists():
                shutil.copy2(csv_path, export_dir / csv_path.name)


def run_mapping_job(job_id: str) -> None:
    job = get_job(job_id)

    try:
        mode = str(job.options.get("mode", "map")).strip().lower()
        if mode == "trait-map":
            cmd = [
                str(PYTHON_BIN),
                str(MAPPER_SCRIPT),
                "trait-map",
                "--input",
                job.input_path,
                "--output",
                job.mapped_path,
                "--trait-cache",
                str(job.options.get("trait_cache_path", DEFAULT_TRAIT_CACHE)),
                "--efo-obo",
                str(job.options.get("efo_obo_path", DEFAULT_EFO_OBO)),
                "--ukb-field-catalog",
                str(job.options.get("ukb_field_catalog_path", DEFAULT_UKB_FIELD_CATALOG)),
                "--top-k",
                str(job.options.get("top_k", 1)),
                "--min-score",
                str(job.options.get("min_score", 0.82)),
                "--default-trait-scale",
                str(job.options.get("default_trait_scale", "auto")),
                "--review-output",
                job.review_path,
                "--qc-risk-output",
                job.qc_risk_path,
                "--qc-summary-output",
                job.qc_summary_path,
                "--flush-every",
                str(job.options.get("flush_every", 50)),
                "--progress",
            ]
            if job.options.get("stream_output", True):
                cmd.append("--stream-output")
            else:
                cmd.append("--no-stream-output")
            if job.options.get("memoize_queries", True):
                cmd.append("--memoize-queries")
            else:
                cmd.append("--no-memoize-queries")
            if job.options.get("force_map_best", False):
                cmd.append("--force-map-best")
        else:
            cmd = [
                str(PYTHON_BIN),
                str(MAPPER_SCRIPT),
                "map",
                "--input",
                job.input_path,
                "--output",
                job.mapped_path,
                "--index",
                str(job.options.get("index_path", DEFAULT_INDEX)),
                "--unmapped-output",
                job.unmapped_path,
                "--top-k",
                str(job.options.get("top_k", 1)),
                "--min-score",
                str(job.options.get("min_score", 0.55)),
                "--workers",
                str(job.options.get("workers", 8)),
                "--measurement-context",
                str(job.options.get("measurement_context", "blood")),
                "--entity-type",
                str(job.options.get("entity_type", "auto")),
                "--matrix-priority",
                str(job.options.get("matrix_priority", "plasma,blood,serum")),
                "--progress",
            ]

            if job.options.get("auto_enrich", True):
                cmd.extend(
                    [
                        "--auto-enrich-uniprot",
                        "--uniprot-aliases",
                        str(job.options.get("uniprot_aliases_path", DEFAULT_UNIPROT_ALIASES)),
                        "--term-cache",
                        str(job.options.get("term_cache_path", DEFAULT_TERM_CACHE)),
                    ]
                )

            if job.options.get("cache_writeback", False):
                cmd.extend(
                    [
                        "--cache-writeback",
                        "--analyte-cache",
                        str(job.options.get("analyte_cache_path", DEFAULT_ANALYTE_CACHE)),
                        "--cache-min-confidence",
                        str(job.options.get("cache_min_confidence", 0.80)),
                    ]
                )

            if job.options.get("force_map_best", False):
                cmd.append("--force-map-best")
                fallback_id = str(job.options.get("fallback_efo_id", "")).strip()
                if fallback_id:
                    cmd.extend(["--fallback-efo-id", fallback_id])

        job.status = "running"
        job.message = "mapping" if mode == "map" else "trait-mapping"
        job.updated_at = time.time()

        proc = subprocess.Popen(
            cmd,
            cwd=str(ROOT_DIR),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert proc.stdout is not None
        for line in proc.stdout:
            with JOBS_LOCK:
                live = JOBS[job_id]
                add_log(live, line)
                set_progress_from_line(live, line)
                live.updated_at = time.time()

        ret = proc.wait()

        with JOBS_LOCK:
            live = JOBS[job_id]
            live.updated_at = time.time()
            if ret != 0:
                live.status = "failed"
                live.message = "mapper failed"
                if live.logs:
                    live.error = live.logs[-1]
                return

            live.status = "completed"
            live.message = "completed"
            live.progress_percent = 100.0
            if live.progress_total == 0:
                live.progress_total = 1
                live.progress_current = 1
            stage_job_exports(live)

    except Exception as exc:
        with JOBS_LOCK:
            live = JOBS[job_id]
            live.status = "failed"
            live.error = str(exc)
            live.message = "exception"
            live.updated_at = time.time()


@app.get("/", response_class=HTMLResponse)
def home(request: Request) -> str:
    server_fallback = ""
    raw_job_id = (request.query_params.get("job_id", "") or "").strip()
    if raw_job_id:
        safe_job_id = re.sub(r"[^a-zA-Z0-9_-]", "", raw_job_id)
        job_dir = JOB_ROOT / safe_job_id
        mode = "trait-map" if (job_dir / "traits_mapped.tsv").exists() else "map"
        links: list[str] = []
        if (job_dir / "mapped.tsv").exists() or (job_dir / "traits_mapped.tsv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/mapped">Download mapped</a>')
        if (job_dir / "mapped.csv").exists() or (job_dir / "traits_mapped.csv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/mapped_csv">Download mapped CSV</a>')
        if (job_dir / "unmapped.tsv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/unmapped">Download unmapped</a>')
        if (job_dir / "unmapped.csv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/unmapped_csv">Download unmapped CSV</a>')
        if (job_dir / "traits_review.tsv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/review">Download review</a>')
        if (job_dir / "traits_review.csv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/review_csv">Download review CSV</a>')
        if (job_dir / "traits_qc_risk.tsv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/qc_risk">Download QC risk</a>')
        if (job_dir / "traits_qc_risk.csv").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/qc_risk_csv">Download QC risk CSV</a>')
        if (job_dir / "traits_qc_summary.json").exists():
            links.append(f'<a href="/api/jobs/{safe_job_id}/download/qc_summary">Download QC summary</a>')
        if links:
            server_fallback = (
                f'<div class="status">job={safe_job_id} mode={mode} (server fallback)</div>'
                f'<div class="status">Local output folder: {str(EXPORT_ROOT / safe_job_id)}</div>'
                f'<div class="links">{"".join(links)}</div>'
            )
        else:
            server_fallback = (
                f'<div class="status">job={safe_job_id} queued/running (server fallback). '
                "Refresh this page in a few seconds.</div>"
            )

    html = """<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width,initial-scale=1\" />
  <title>pQTL Mapper</title>
  <style>
    :root {
      --bg: #f5f7f9;
      --card: #ffffff;
      --ink: #12222f;
      --muted: #4f6879;
      --line: #d3dde5;
      --accent: #0d6d63;
      --accent-2: #0f8f81;
      --danger: #b42318;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: ui-sans-serif, -apple-system, BlinkMacSystemFont, \"Segoe UI\", sans-serif;
      background: radial-gradient(1200px 600px at -10% -10%, #d8efe8 0%, var(--bg) 55%);
      color: var(--ink);
    }
    .wrap {
      max-width: 980px;
      margin: 24px auto;
      padding: 0 16px 28px;
    }
    .hero {
      background: linear-gradient(135deg, #e8faf4, #f6fbff);
      border: 1px solid #cfe9dd;
      border-radius: 16px;
      padding: 18px;
      margin-bottom: 16px;
    }
    h1 {
      margin: 0 0 6px;
      font-size: 1.45rem;
      letter-spacing: 0.01em;
    }
    .sub {
      margin: 0;
      color: var(--muted);
      font-size: 0.96rem;
    }
    .grid {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 12px;
    }
    @media (max-width: 860px) {
      .grid { grid-template-columns: 1fr; }
    }
    .card {
      background: var(--card);
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 14px;
    }
    label {
      display: block;
      font-size: 0.85rem;
      color: var(--muted);
      margin: 10px 0 6px;
      font-weight: 600;
    }
    input[type=file], input[type=text], input[type=number], select {
      width: 100%;
      border: 1px solid #bcd0dd;
      border-radius: 8px;
      padding: 9px 10px;
      font-size: 0.95rem;
      background: #fff;
    }
    .row {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 8px;
    }
    .check {
      display: flex;
      align-items: center;
      gap: 8px;
      margin-top: 10px;
      color: var(--muted);
      font-size: 0.9rem;
    }
    .btn {
      border: 0;
      border-radius: 10px;
      background: linear-gradient(120deg, var(--accent), var(--accent-2));
      color: #fff;
      font-weight: 700;
      letter-spacing: 0.02em;
      padding: 11px 14px;
      cursor: pointer;
      width: 100%;
      margin-top: 12px;
    }
    .btn:disabled {
      opacity: 0.6;
      cursor: not-allowed;
    }
    .status {
      font-size: 0.93rem;
      color: var(--muted);
      margin-top: 8px;
    }
    .bar {
      height: 12px;
      background: #e7eef3;
      border-radius: 999px;
      overflow: hidden;
      margin-top: 10px;
      border: 1px solid #d6e1ea;
    }
    .bar > div {
      height: 100%;
      width: 0%;
      background: linear-gradient(90deg, #0f8f81, #18b79f);
      transition: width 0.35s ease;
    }
    .links {
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
      margin-top: 12px;
    }
    .links a {
      font-size: 0.84rem;
      border: 1px solid #bfd3e0;
      border-radius: 8px;
      text-decoration: none;
      color: #10405e;
      padding: 7px 9px;
      background: #f8fcff;
    }
    #logs {
      margin: 0;
      background: #0c1f2b;
      color: #cde8f2;
      border-radius: 10px;
      border: 1px solid #19435a;
      min-height: 210px;
      max-height: 410px;
      overflow: auto;
      padding: 12px;
      font-size: 0.79rem;
      line-height: 1.45;
    }
    .help {
      margin-top: 12px;
    }
    .help p {
      margin: 8px 0;
      color: var(--muted);
      font-size: 0.92rem;
    }
    .help strong {
      color: var(--ink);
    }
    .cmd {
      margin: 6px 0 10px;
      background: #0c1f2b;
      color: #cde8f2;
      border-radius: 10px;
      border: 1px solid #19435a;
      overflow: auto;
      padding: 10px 12px;
      font-size: 0.82rem;
      line-height: 1.45;
      white-space: pre;
    }
    .cmd code {
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, \"Courier New\", monospace;
    }
    .err { color: var(--danger); }
  </style>
</head>
<body>
  <div class=\"wrap\">
    <section class=\"hero\">
      <h1>pQTL Measurement Mapper</h1>
      <p class=\"sub\">Run analyte map (protein/metabolite) or trait-map, track progress, and download result artifacts.</p>
      <p class=\"sub\" style=\"margin-top:8px;\">
        Browser-only prototype available:
        <a href=\"/pyodide\" style=\"font-weight:700;color:#0d6d63;text-decoration:none;\">Open Pyodide Trait Mapper</a>
      </p>
    </section>

    <div class=\"grid\">
      <section class=\"card\">
        <form id=\"jobForm\" method=\"post\" action=\"/api/jobs\" enctype=\"multipart/form-data\">
          <label>Input File (.txt/.tsv/.csv)</label>
          <input name=\"file\" type=\"file\" required />

          <label>Mode</label>
          <select name=\"mode\">
            <option value=\"map\" selected>map (analyte: protein/metabolite)</option>
            <option value=\"trait-map\">trait-map (disease/phenotype)</option>
          </select>

          <div class=\"row\">
            <div>
              <label>Measurement Context</label>
              <select name=\"measurement_context\">
                <option value=\"blood\" selected>blood</option>
                <option value=\"cerebrospinal_fluid\">cerebrospinal_fluid</option>
                <option value=\"urine\">urine</option>
                <option value=\"saliva\">saliva</option>
                <option value=\"tissue\">tissue</option>
                <option value=\"auto\">auto</option>
              </select>
            </div>
            <div>
              <label>Workers</label>
              <input name=\"workers\" type=\"number\" min=\"1\" value=\"8\" />
            </div>
          </div>

          <label>Entity Type</label>
          <select name=\"entity_type\">
            <option value=\"auto\" selected>auto (Recommended for mixed panels)</option>
            <option value=\"protein\">protein</option>
            <option value=\"metabolite\">metabolite (Advanced: metabolite-only rows; protein-like rows are blocked)</option>
          </select>

          <div class=\"row\">
            <div>
              <label>Top K</label>
              <input name=\"top_k\" type=\"number\" min=\"1\" value=\"1\" />
            </div>
            <div>
              <label>Min Score</label>
              <input name=\"min_score\" type=\"number\" min=\"0\" max=\"1\" step=\"0.01\" value=\"0.55\" />
            </div>
          </div>

          <label>Matrix Priority</label>
          <input name=\"matrix_priority\" type=\"text\" value=\"plasma,blood,serum\" />

          <label class=\"check\"><input name=\"auto_enrich\" type=\"checkbox\" checked />Auto-enrich missing UniProt aliases</label>
          <label class=\"check\"><input name=\"cache_writeback\" type=\"checkbox\" />Write high-confidence results to analyte cache</label>
          <label class=\"check\"><input name=\"force_map_best\" type=\"checkbox\" />Force fallback row when unresolved</label>

          <label>Fallback EFO ID (optional)</label>
          <input name=\"fallback_efo_id\" type=\"text\" placeholder=\"EFO:0001444\" />

          <label>Trait Default Scale</label>
          <select name=\"default_trait_scale\">
            <option value=\"auto\" selected>auto</option>
            <option value=\"binary\">binary</option>
            <option value=\"quantitative\">quantitative</option>
          </select>

          <div class=\"row\">
            <div>
              <label>Trait Flush Every</label>
              <input name=\"flush_every\" type=\"number\" min=\"1\" value=\"50\" />
            </div>
            <div></div>
          </div>

          <label class=\"check\"><input name=\"stream_output\" type=\"checkbox\" checked />Trait: stream output while running</label>
          <label class=\"check\"><input name=\"memoize_queries\" type=\"checkbox\" checked />Trait: memoize duplicate queries</label>

          <button class=\"btn\" id=\"runBtn\" type=\"submit\">Run Mapping</button>
          <noscript><div class=\"status err\">JavaScript is disabled. Submit will return JSON only; enable JavaScript for live progress and download links.</div></noscript>
          <div class=\"status\" id=\"status\">Idle</div>
          <div class=\"bar\"><div id=\"barFill\"></div></div>

          <div class=\"links\" id=\"links\"></div>
          __SERVER_FALLBACK__
        </form>
      </section>

      <section class=\"card\">
        <h3 style=\"margin:0 0 8px;font-size:1rem;\">Job Logs</h3>
        <pre id=\"logs\"></pre>
        <div class=\"status\" id=\"jobMeta\"></div>
      </section>
    </div>

    <section class=\"card help\">
      <h3 style=\"margin:0 0 8px;font-size:1rem;\">Update / Pull Fixes</h3>
      <p>Run mapper commands from the repo root with <code>analyte-efo-mapper ...</code>.</p>
      <p><strong>Upgrade / pull latest (external-user flow)</strong></p>
      <div class=\"cmd\"><code>git pull --rebase origin main
python -m pip install -e .
analyte-efo-mapper setup-bundled-caches</code></div>
      <p>Bundled UKB dictionaries used by trait mode:
      <code>references/ukb/fieldsum.txt</code>, <code>references/ukb/field.txt</code>,
      <code>references/ukb/category.txt</code>, <code>references/ukb/catbrowse.txt</code>.</p>
      <p><strong>If you have local uncommitted changes</strong></p>
      <div class=\"cmd\"><code>git stash
git pull --rebase origin main
python -m pip install -e .
analyte-efo-mapper setup-bundled-caches
git stash pop</code></div>
      <p><strong>After mapper/cache changes</strong></p>
      <div class=\"cmd\"><code>analyte-efo-mapper index-build</code></div>
      <p>Restart this web app after update so the new code and index are loaded.</p>
      <p>If you are not on <code>main</code>, replace it with your current branch.</p>
    </section>
  </div>

<script>
const form = document.getElementById('jobForm');
const statusEl = document.getElementById('status');
const barFill = document.getElementById('barFill');
const logsEl = document.getElementById('logs');
const linksEl = document.getElementById('links');
const runBtn = document.getElementById('runBtn');
const jobMeta = document.getElementById('jobMeta');

let timer = null;
let currentJobId = null;

function boolFromForm(fd, key) {
  return fd.get(key) !== null;
}

function setStatus(msg, isError=false) {
  statusEl.textContent = msg;
  statusEl.className = isError ? 'status err' : 'status';
}

function setProgress(percent) {
  const pct = Math.max(0, Math.min(100, percent || 0));
  barFill.style.width = pct.toFixed(1) + '%';
}

function renderLinks(downloads, ready) {
  linksEl.innerHTML = '';
  if (!ready || !downloads) return;
  const labels = {
    mapped: 'Download mapped TSV',
    mapped_csv: 'Download mapped CSV',
    unmapped: 'Download unmapped TSV',
    unmapped_csv: 'Download unmapped CSV',
    review: 'Download review TSV',
    review_csv: 'Download review CSV',
    qc_risk: 'Download QC risk TSV',
    qc_risk_csv: 'Download QC risk CSV',
    qc_summary: 'Download QC summary JSON',
  };
  for (const [key, href] of Object.entries(downloads)) {
    const a = document.createElement('a');
    a.href = href;
    a.textContent = labels[key] || `Download ${key}`;
    linksEl.appendChild(a);
  }
}

async function pollJob(jobId) {
  try {
    const res = await fetch(`/api/jobs/${jobId}`);
    if (!res.ok) throw new Error(`status ${res.status}`);
    const data = await res.json();

    setProgress(data.progress_percent || 0);
    setStatus(`${data.status}: ${data.message || ''}`.trim(), data.status === 'failed');
    logsEl.textContent = (data.logs || []).join('\n');
    logsEl.scrollTop = logsEl.scrollHeight;
    jobMeta.textContent = `job=${data.job_id} updated=${new Date(data.updated_at * 1000).toLocaleTimeString()}`;

    if (data.status === 'completed') {
      clearInterval(timer);
      timer = null;
      runBtn.disabled = false;
      renderLinks(data.downloads, true);
      return;
    }
    if (data.status === 'failed') {
      clearInterval(timer);
      timer = null;
      runBtn.disabled = false;
      renderLinks(null, false);
      return;
    }
  } catch (err) {
    setStatus(`Polling failed: ${err}`, true);
    clearInterval(timer);
    timer = null;
    runBtn.disabled = false;
  }
}

form.addEventListener('submit', async (e) => {
  e.preventDefault();
  if (timer) {
    clearInterval(timer);
    timer = null;
  }

  const fd = new FormData(form);
  fd.set('auto_enrich', boolFromForm(fd, 'auto_enrich'));
  fd.set('cache_writeback', boolFromForm(fd, 'cache_writeback'));
  fd.set('force_map_best', boolFromForm(fd, 'force_map_best'));
  fd.set('stream_output', boolFromForm(fd, 'stream_output'));
  fd.set('memoize_queries', boolFromForm(fd, 'memoize_queries'));

  runBtn.disabled = true;
  renderLinks(null, false);
  logsEl.textContent = '';
  setProgress(0);
  setStatus('Submitting job...');

  try {
    const res = await fetch('/api/jobs', {
      method: 'POST',
      headers: { 'X-Requested-With': 'fetch' },
      body: fd
    });
    if (!res.ok) {
      const txt = await res.text();
      throw new Error(txt || `status ${res.status}`);
    }
    const data = await res.json();
    currentJobId = data.job_id;
    setStatus(`queued: ${data.job_id}`);
    timer = setInterval(() => pollJob(currentJobId), 1000);
    pollJob(currentJobId);
  } catch (err) {
    runBtn.disabled = false;
    setStatus(`Submit failed: ${err}`, true);
  }
});

const qs = new URLSearchParams(window.location.search);
const restoreJobId = qs.get('job_id');
if (restoreJobId) {
  currentJobId = restoreJobId;
  runBtn.disabled = true;
  setStatus(`queued: ${restoreJobId}`);
  timer = setInterval(() => pollJob(currentJobId), 1000);
  pollJob(currentJobId);
}
</script>
</body>
</html>
"""
    return html.replace("__SERVER_FALLBACK__", server_fallback)


def _resolve_pyodide_asset(asset_path: str) -> Path:
    safe = (asset_path or "").lstrip("/")
    candidate = (PYODIDE_WEB_DIR / safe).resolve()
    root = PYODIDE_WEB_DIR.resolve()
    if root not in candidate.parents and candidate != root:
        raise HTTPException(status_code=404, detail="Asset not found")
    if not candidate.exists() or not candidate.is_file():
        raise HTTPException(status_code=404, detail="Asset not found")
    return candidate


@app.get("/pyodide")
def pyodide_home() -> FileResponse:
    index_path = PYODIDE_WEB_DIR / "index.html"
    if not index_path.exists():
        raise HTTPException(status_code=404, detail="Pyodide prototype not available")
    return FileResponse(path=str(index_path), media_type="text/html")


@app.get("/pyodide/{asset_path:path}")
def pyodide_asset(asset_path: str) -> FileResponse:
    path = _resolve_pyodide_asset(asset_path)
    media_type = None
    if path.suffix == ".js":
        media_type = "application/javascript"
    elif path.suffix == ".json":
        media_type = "application/json"
    elif path.suffix == ".zip":
        media_type = "application/zip"
    elif path.suffix == ".py":
        media_type = "text/x-python"
    elif path.suffix == ".css":
        media_type = "text/css"
    return FileResponse(path=str(path), media_type=media_type)


@app.post("/api/jobs")
def create_job(
    request: Request,
    file: UploadFile = File(...),
    mode: str = Form("map"),
    workers: int = Form(8),
    top_k: int = Form(1),
    min_score: float = Form(0.55),
    matrix_priority: str = Form("plasma,blood,serum"),
    measurement_context: str = Form("blood"),
    entity_type: str = Form("auto"),
    auto_enrich: bool = Form(True),
    cache_writeback: bool = Form(False),
    force_map_best: bool = Form(False),
    fallback_efo_id: str = Form(""),
    default_trait_scale: str = Form("auto"),
    stream_output: bool = Form(True),
    memoize_queries: bool = Form(True),
    flush_every: int = Form(50),
) -> JSONResponse:
    if not file.filename:
        raise HTTPException(status_code=400, detail="Missing input filename")

    ext = Path(file.filename).suffix.lower()
    if ext not in {".txt", ".tsv", ".csv", ".list", ".tab"}:
        raise HTTPException(status_code=400, detail="Unsupported file type")
    mode_clean = (mode or "").strip().lower() or "map"
    if mode_clean not in {"map", "trait-map"}:
        raise HTTPException(status_code=400, detail="Invalid mode")

    entity_type_clean = (entity_type or "").strip().lower() or "auto"
    if entity_type_clean not in {"auto", "protein", "metabolite"}:
        raise HTTPException(status_code=400, detail="Invalid entity_type")
    default_trait_scale_clean = (default_trait_scale or "").strip().lower() or "auto"
    if default_trait_scale_clean not in {"auto", "binary", "quantitative"}:
        raise HTTPException(status_code=400, detail="Invalid default_trait_scale")

    job_id = uuid.uuid4().hex[:12]
    job_dir = JOB_ROOT / job_id
    export_dir = EXPORT_ROOT / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    export_dir.mkdir(parents=True, exist_ok=True)

    input_path = job_dir / f"input{ext if ext else '.txt'}"
    mapped_path = job_dir / ("traits_mapped.tsv" if mode_clean == "trait-map" else "mapped.tsv")
    unmapped_path = job_dir / "unmapped.tsv"
    review_path = job_dir / "traits_review.tsv"
    qc_risk_path = job_dir / "traits_qc_risk.tsv"
    qc_summary_path = job_dir / "traits_qc_summary.json"

    content = file.file.read()
    input_path.write_bytes(content)

    now = time.time()
    job = JobState(
        job_id=job_id,
        status="queued",
        created_at=now,
        updated_at=now,
        input_path=str(input_path),
        mapped_path=str(mapped_path),
        unmapped_path=str(unmapped_path),
        review_path=str(review_path),
        qc_risk_path=str(qc_risk_path),
        qc_summary_path=str(qc_summary_path),
        local_export_dir=str(export_dir),
        options={
            "mode": mode_clean,
            "workers": max(1, int(workers)),
            "top_k": max(1, int(top_k)),
            "min_score": max(0.0, float(min_score)),
            "matrix_priority": matrix_priority.strip() or "plasma,blood,serum",
            "measurement_context": measurement_context.strip() or "blood",
            "entity_type": entity_type_clean,
            "auto_enrich": bool(auto_enrich),
            "cache_writeback": bool(cache_writeback),
            "force_map_best": bool(force_map_best),
            "fallback_efo_id": fallback_efo_id.strip(),
            "default_trait_scale": default_trait_scale_clean,
            "stream_output": bool(stream_output),
            "memoize_queries": bool(memoize_queries),
            "flush_every": max(1, int(flush_every)),
            "index_path": str(DEFAULT_INDEX),
            "uniprot_aliases_path": str(DEFAULT_UNIPROT_ALIASES),
            "term_cache_path": str(DEFAULT_TERM_CACHE),
            "analyte_cache_path": str(DEFAULT_ANALYTE_CACHE),
            "trait_cache_path": str(DEFAULT_TRAIT_CACHE),
            "efo_obo_path": str(DEFAULT_EFO_OBO),
            "ukb_field_catalog_path": str(DEFAULT_UKB_FIELD_CATALOG),
        },
    )

    with JOBS_LOCK:
        JOBS[job_id] = job

    thread = threading.Thread(target=run_mapping_job, args=(job_id,), daemon=True)
    thread.start()

    requested_with = (request.headers.get("x-requested-with", "") or "").strip().lower()
    accept = request.headers.get("accept", "")
    if requested_with != "fetch" and "text/html" in accept:
        return RedirectResponse(url=f"/?job_id={job_id}", status_code=303)

    return JSONResponse({"job_id": job_id, "status": "queued"})


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str) -> JSONResponse:
    job = get_job(job_id)
    with JOBS_LOCK:
        payload = {
            "job_id": job.job_id,
            "status": job.status,
            "created_at": job.created_at,
            "updated_at": job.updated_at,
            "progress_current": job.progress_current,
            "progress_total": job.progress_total,
            "progress_percent": job.progress_percent,
            "message": job.message,
            "error": job.error,
            "logs": job.logs,
            "local_export_dir": job.local_export_dir,
            "local_paths": {
                "input": str(Path(job.local_export_dir) / Path(job.input_path).name),
                "mapped": str(Path(job.local_export_dir) / Path(job.mapped_path).name),
                "unmapped": str(Path(job.local_export_dir) / Path(job.unmapped_path).name),
                "review": str(Path(job.local_export_dir) / Path(job.review_path).name),
                "qc_risk": str(Path(job.local_export_dir) / Path(job.qc_risk_path).name),
                "qc_summary": str(Path(job.local_export_dir) / Path(job.qc_summary_path).name),
            },
            "downloads": (
                {
                    "mapped": f"/api/jobs/{job_id}/download/mapped",
                    "mapped_csv": f"/api/jobs/{job_id}/download/mapped_csv",
                    "review": f"/api/jobs/{job_id}/download/review",
                    "review_csv": f"/api/jobs/{job_id}/download/review_csv",
                    "qc_risk": f"/api/jobs/{job_id}/download/qc_risk",
                    "qc_risk_csv": f"/api/jobs/{job_id}/download/qc_risk_csv",
                    "qc_summary": f"/api/jobs/{job_id}/download/qc_summary",
                }
                if job.options.get("mode") == "trait-map"
                else {
                    "mapped": f"/api/jobs/{job_id}/download/mapped",
                    "mapped_csv": f"/api/jobs/{job_id}/download/mapped_csv",
                    "unmapped": f"/api/jobs/{job_id}/download/unmapped",
                    "unmapped_csv": f"/api/jobs/{job_id}/download/unmapped_csv",
                }
            ),
        }
    return JSONResponse(payload)


@app.get("/api/jobs/{job_id}/download/{kind}")
def download_job_file(job_id: str, kind: str) -> FileResponse:
    job = get_job(job_id)
    media_type = "text/tab-separated-values"

    if kind == "mapped":
        path = Path(job.mapped_path)
    elif kind == "mapped_csv":
        path = Path(job.mapped_path).with_suffix(".csv")
        media_type = "text/csv"
    elif kind == "unmapped":
        path = Path(job.unmapped_path)
    elif kind == "unmapped_csv":
        path = Path(job.unmapped_path).with_suffix(".csv")
        media_type = "text/csv"
    elif kind == "review":
        path = Path(job.review_path)
    elif kind == "review_csv":
        path = Path(job.review_path).with_suffix(".csv")
        media_type = "text/csv"
    elif kind == "qc_risk":
        path = Path(job.qc_risk_path)
    elif kind == "qc_risk_csv":
        path = Path(job.qc_risk_path).with_suffix(".csv")
        media_type = "text/csv"
    elif kind == "qc_summary":
        path = Path(job.qc_summary_path)
        media_type = "application/json"
    else:
        raise HTTPException(status_code=404, detail="Unknown artifact")

    if not path.exists():
        raise HTTPException(status_code=404, detail="File not available yet")

    return FileResponse(path=str(path), filename=path.name, media_type=media_type)


@app.get("/healthz")
def health() -> JSONResponse:
    return JSONResponse({"ok": True})
