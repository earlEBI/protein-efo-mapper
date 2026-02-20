#!/usr/bin/env python3
"""Minimal web app for pQTL measurement mapping.

Features:
- Upload analyte file (.txt/.tsv/.csv)
- Run mapper in background job
- Live progress/log polling
- Download mapped + unmapped TSV outputs
"""

from __future__ import annotations

import re
import subprocess
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse

ROOT_DIR = Path(__file__).resolve().parents[3]
SKILL_DIR = ROOT_DIR / "skills" / "pqtl-measurement-mapper"
MAPPER_SCRIPT = SKILL_DIR / "scripts" / "map_measurement_efo.py"
DEFAULT_INDEX = SKILL_DIR / "references" / "measurement_index.json"
DEFAULT_UNIPROT_ALIASES = SKILL_DIR / "references" / "uniprot_aliases.tsv"
DEFAULT_TERM_CACHE = SKILL_DIR / "references" / "efo_measurement_terms_cache.tsv"
DEFAULT_ANALYTE_CACHE = SKILL_DIR / "references" / "analyte_to_efo_cache.tsv"
PYTHON_BIN = ROOT_DIR / ".venv" / "bin" / "python"
JOB_ROOT = ROOT_DIR / "final_output" / "web_jobs"
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


def run_mapping_job(job_id: str) -> None:
    job = get_job(job_id)

    try:
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
        job.message = "mapping"
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

    except Exception as exc:
        with JOBS_LOCK:
            live = JOBS[job_id]
            live.status = "failed"
            live.error = str(exc)
            live.message = "exception"
            live.updated_at = time.time()


@app.get("/", response_class=HTMLResponse)
def home() -> str:
    return """<!doctype html>
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
    pre {
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
    .err { color: var(--danger); }
  </style>
</head>
<body>
  <div class=\"wrap\">
    <section class=\"hero\">
      <h1>pQTL Measurement Mapper</h1>
      <p class=\"sub\">Upload analytes, run deterministic mapping, track progress, and download mapped/unmapped results.</p>
    </section>

    <div class=\"grid\">
      <section class=\"card\">
        <form id=\"jobForm\">
          <label>Input File (.txt/.tsv/.csv)</label>
          <input name=\"file\" type=\"file\" required />

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

          <button class=\"btn\" id=\"runBtn\" type=\"submit\">Run Mapping</button>
          <div class=\"status\" id=\"status\">Idle</div>
          <div class=\"bar\"><div id=\"barFill\"></div></div>

          <div class=\"links\" id=\"links\"></div>
        </form>
      </section>

      <section class=\"card\">
        <h3 style=\"margin:0 0 8px;font-size:1rem;\">Job Logs</h3>
        <pre id=\"logs\"></pre>
        <div class=\"status\" id=\"jobMeta\"></div>
      </section>
    </div>
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

function renderLinks(jobId, ready) {
  linksEl.innerHTML = '';
  if (!ready) return;
  const mapped = document.createElement('a');
  mapped.href = `/api/jobs/${jobId}/download/mapped`;
  mapped.textContent = 'Download mapped TSV';
  linksEl.appendChild(mapped);

  const unmapped = document.createElement('a');
  unmapped.href = `/api/jobs/${jobId}/download/unmapped`;
  unmapped.textContent = 'Download unmapped TSV';
  linksEl.appendChild(unmapped);
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
      renderLinks(jobId, true);
      return;
    }
    if (data.status === 'failed') {
      clearInterval(timer);
      timer = null;
      runBtn.disabled = false;
      renderLinks(jobId, false);
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

  runBtn.disabled = true;
  renderLinks('', false);
  logsEl.textContent = '';
  setProgress(0);
  setStatus('Submitting job...');

  try {
    const res = await fetch('/api/jobs', { method: 'POST', body: fd });
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
</script>
</body>
</html>
"""


@app.post("/api/jobs")
def create_job(
    file: UploadFile = File(...),
    workers: int = Form(8),
    top_k: int = Form(1),
    min_score: float = Form(0.55),
    matrix_priority: str = Form("plasma,blood,serum"),
    measurement_context: str = Form("blood"),
    auto_enrich: bool = Form(True),
    cache_writeback: bool = Form(False),
    force_map_best: bool = Form(False),
    fallback_efo_id: str = Form(""),
) -> JSONResponse:
    if not file.filename:
        raise HTTPException(status_code=400, detail="Missing input filename")

    ext = Path(file.filename).suffix.lower()
    if ext not in {".txt", ".tsv", ".csv", ".list", ".tab"}:
        raise HTTPException(status_code=400, detail="Unsupported file type")

    job_id = uuid.uuid4().hex[:12]
    job_dir = JOB_ROOT / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    input_path = job_dir / f"input{ext if ext else '.txt'}"
    mapped_path = job_dir / "mapped.tsv"
    unmapped_path = job_dir / "unmapped.tsv"

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
        options={
            "workers": max(1, int(workers)),
            "top_k": max(1, int(top_k)),
            "min_score": max(0.0, float(min_score)),
            "matrix_priority": matrix_priority.strip() or "plasma,blood,serum",
            "measurement_context": measurement_context.strip() or "blood",
            "auto_enrich": bool(auto_enrich),
            "cache_writeback": bool(cache_writeback),
            "force_map_best": bool(force_map_best),
            "fallback_efo_id": fallback_efo_id.strip(),
            "index_path": str(DEFAULT_INDEX),
            "uniprot_aliases_path": str(DEFAULT_UNIPROT_ALIASES),
            "term_cache_path": str(DEFAULT_TERM_CACHE),
            "analyte_cache_path": str(DEFAULT_ANALYTE_CACHE),
        },
    )

    with JOBS_LOCK:
        JOBS[job_id] = job

    thread = threading.Thread(target=run_mapping_job, args=(job_id,), daemon=True)
    thread.start()

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
            "downloads": {
                "mapped": f"/api/jobs/{job_id}/download/mapped",
                "unmapped": f"/api/jobs/{job_id}/download/unmapped",
            },
        }
    return JSONResponse(payload)


@app.get("/api/jobs/{job_id}/download/{kind}")
def download_job_file(job_id: str, kind: str) -> FileResponse:
    job = get_job(job_id)

    if kind == "mapped":
        path = Path(job.mapped_path)
    elif kind == "unmapped":
        path = Path(job.unmapped_path)
    else:
        raise HTTPException(status_code=404, detail="Unknown artifact")

    if not path.exists():
        raise HTTPException(status_code=404, detail="File not available yet")

    return FileResponse(path=str(path), filename=path.name, media_type="text/tab-separated-values")


@app.get("/healthz")
def health() -> JSONResponse:
    return JSONResponse({"ok": True})
