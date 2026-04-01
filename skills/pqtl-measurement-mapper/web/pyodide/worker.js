/* global loadPyodide */

const PYODIDE_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.3/full/";
const CORE_PATH = "./mapper_core.py";
const DEFAULT_CACHE_PATH = "./cache_sample.json";
const DEFAULT_RUNTIME_BUNDLE_PATH = "./trait_runtime_bundle.zip";
const RUNTIME_ZIP_FS_PATH = "/runtime/trait_runtime_bundle.zip";
const RUNTIME_REPO_ROOT = "/runtime/repo";
const RUNTIME_WORK_ROOT = "/runtime/work";
const SUPPORTED_ENGINES = new Set(["fast", "cli_compat"]);
const SUPPORTED_TASK_TYPES = new Set(["trait", "analyte"]);

let pyodideReadyPromise = null;
let pyodide = null;
let defaultCacheText = "";
let runtimeBundleReady = false;
let runtimeBundleKey = "";
let currentRunLogLines = [];

function resetRunLogs() {
  currentRunLogLines = [];
}

function pushRunLog(line) {
  currentRunLogLines.push(line);
  if (currentRunLogLines.length > 3000) {
    currentRunLogLines = currentRunLogLines.slice(-3000);
  }
}

function getRunLogsText() {
  return currentRunLogLines.join("\n");
}

function postStatus(message) {
  self.postMessage({ type: "status", payload: { message } });
}

function postLog(line, stream = "stdout") {
  self.postMessage({ type: "log", payload: { line, stream } });
}

function postProgress(current, total, percent, line = "") {
  self.postMessage({
    type: "progress",
    payload: {
      current: Number(current) || 0,
      total: Number(total) || 0,
      percent: Number(percent) || 0,
      line: String(line || ""),
    },
  });
}

function handleRuntimeText(rawText, stream = "stdout") {
  if (!rawText) {
    return;
  }
  const lines = String(rawText).split(/\r?\n/);
  for (const rawLine of lines) {
    const line = String(rawLine || "").trim();
    if (!line) {
      continue;
    }
    pushRunLog(line);
    postLog(line, stream);
    const match = line.match(/^\[progress\]\s+(\d+)\/(\d+)\s+\(([\d.]+)%\)/i);
    if (match) {
      postProgress(match[1], match[2], match[3], line);
    }
  }
}

async function loadPyodideRuntime() {
  if (pyodideReadyPromise) {
    return pyodideReadyPromise;
  }

  pyodideReadyPromise = (async () => {
    postStatus("Loading Pyodide runtime...");
    if (typeof loadPyodide === "undefined") {
      importScripts(`${PYODIDE_INDEX_URL}pyodide.js`);
    }
    pyodide = await loadPyodide({ indexURL: PYODIDE_INDEX_URL });
    pyodide.setStdout({ batched: (text) => handleRuntimeText(text, "stdout") });
    pyodide.setStderr({ batched: (text) => handleRuntimeText(text, "stderr") });

    postStatus("Loading fast mapper core...");
    const coreResponse = await fetch(CORE_PATH);
    if (!coreResponse.ok) {
      throw new Error(`Failed to load ${CORE_PATH}`);
    }
    const coreSource = await coreResponse.text();
    await pyodide.runPythonAsync(coreSource);

    const cacheResponse = await fetch(DEFAULT_CACHE_PATH);
    if (cacheResponse.ok) {
      defaultCacheText = await cacheResponse.text();
    } else {
      defaultCacheText = "";
      postStatus(`Default fast cache not found (${DEFAULT_CACHE_PATH}); upload a cache JSON when using fast mode.`);
    }

    postStatus("Pyodide mapper ready.");
  })();

  return pyodideReadyPromise;
}

function parseEngine(payload) {
  const raw = String(payload?.engine ?? "fast").trim().toLowerCase();
  if (SUPPORTED_ENGINES.has(raw)) {
    return raw;
  }
  return "fast";
}

function parseTaskType(payload) {
  const raw = String(payload?.taskType ?? "trait").trim().toLowerCase();
  if (SUPPORTED_TASK_TYPES.has(raw)) {
    return raw;
  }
  return "trait";
}

function parseEntityType(payload) {
  const raw = String(payload?.entityType ?? "auto").trim().toLowerCase();
  if (["auto", "protein", "metabolite"].includes(raw)) {
    return raw;
  }
  return "auto";
}

function parseMeasurementContext(payload) {
  const raw = String(payload?.measurementContext ?? "blood").trim().toLowerCase();
  if (
    [
      "blood",
      "plasma",
      "serum",
      "cerebrospinal_fluid",
      "urine",
      "saliva",
      "tissue",
      "auto",
    ].includes(raw)
  ) {
    return raw;
  }
  return "blood";
}

async function ensureRuntimeBundle(runtimeBundleBytes, runtimeBundleName) {
  await loadPyodideRuntime();

  let bundleBytes = null;
  let bundleKey = "";
  if (runtimeBundleBytes && runtimeBundleBytes.byteLength > 0) {
    bundleBytes = new Uint8Array(runtimeBundleBytes);
    bundleKey = `custom:${String(runtimeBundleName || "runtime_bundle.zip")}:${bundleBytes.byteLength}`;
  } else {
    bundleKey = "default:trait_runtime_bundle.zip";
    if (runtimeBundleReady && runtimeBundleKey === bundleKey) {
      return;
    }
    postStatus("Fetching CLI-compat runtime bundle...");
    const response = await fetch(DEFAULT_RUNTIME_BUNDLE_PATH);
    if (!response.ok) {
      throw new Error(
        "Missing default runtime bundle (trait_runtime_bundle.zip). Build it with "
        + "build_pyodide_trait_runtime_bundle.py or upload a custom runtime ZIP."
      );
    }
    bundleBytes = new Uint8Array(await response.arrayBuffer());
  }

  if (runtimeBundleReady && runtimeBundleKey === bundleKey) {
    return;
  }

  postStatus("Unpacking CLI-compat runtime bundle...");
  pyodide.FS.mkdirTree("/runtime");
  pyodide.FS.writeFile(RUNTIME_ZIP_FS_PATH, bundleBytes);
  pyodide.globals.set("__RUNTIME_ZIP_PATH__", RUNTIME_ZIP_FS_PATH);
  pyodide.globals.set("__RUNTIME_REPO_ROOT__", RUNTIME_REPO_ROOT);
  await pyodide.runPythonAsync(`
import pathlib
import shutil
import zipfile

zip_path = pathlib.Path(__RUNTIME_ZIP_PATH__)
repo_root = pathlib.Path(__RUNTIME_REPO_ROOT__)
if repo_root.exists():
    shutil.rmtree(repo_root)
repo_root.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(zip_path, "r") as archive:
    archive.extractall(repo_root)
`);
  pyodide.globals.delete("__RUNTIME_ZIP_PATH__");
  pyodide.globals.delete("__RUNTIME_REPO_ROOT__");

  runtimeBundleReady = true;
  runtimeBundleKey = bundleKey;
}

async function runFastMap(payload) {
  await loadPyodideRuntime();
  resetRunLogs();

  const taskType = parseTaskType(payload);
  if (taskType !== "trait") {
    throw new Error("Fast mode currently supports trait mapping only. Use cli_compat for analyte mapping.");
  }
  const inputText = String(payload?.inputText ?? "");
  const minScore = Number(payload?.minScore ?? 0.82);
  const cacheText = payload?.cacheText ? String(payload.cacheText) : defaultCacheText;
  if (!cacheText) {
    throw new Error("Fast mode requires a cache JSON. Upload one or provide cache_sample.json.");
  }

  pyodide.globals.set("__BROWSER_INPUT__", inputText);
  pyodide.globals.set("__BROWSER_CACHE__", cacheText);
  pyodide.globals.set("__BROWSER_MIN_SCORE__", minScore);

  try {
    postStatus("Running fast trait mapping in worker...");
    const resultText = await pyodide.runPythonAsync(
      "map_tsv(__BROWSER_INPUT__, __BROWSER_CACHE__, float(__BROWSER_MIN_SCORE__))"
    );
    self.postMessage({
      type: "mapped",
      payload: {
        outputText: resultText,
        logs: getRunLogsText(),
        engine: "fast",
        taskType: "trait",
      },
    });
  } finally {
    pyodide.globals.delete("__BROWSER_INPUT__");
    pyodide.globals.delete("__BROWSER_CACHE__");
    pyodide.globals.delete("__BROWSER_MIN_SCORE__");
  }
}

async function runCliCompatMap(payload) {
  await ensureRuntimeBundle(payload?.runtimeBundleBytes, payload?.runtimeBundleName);
  resetRunLogs();

  const taskType = parseTaskType(payload);
  const inputText = String(payload?.inputText ?? "");
  const minScore = Number(payload?.minScore ?? 0.82);
  const topK = Math.max(1, Number(payload?.topK ?? 1));
  const entityType = parseEntityType(payload);
  const measurementContext = parseMeasurementContext(payload);
  const matrixPriority = String(payload?.matrixPriority ?? "plasma,blood,serum").trim() || "plasma,blood,serum";

  pyodide.FS.mkdirTree(RUNTIME_WORK_ROOT);
  pyodide.FS.writeFile(`${RUNTIME_WORK_ROOT}/input.tsv`, inputText, { encoding: "utf8" });

  pyodide.globals.set("__CLI_RUNTIME_ROOT__", RUNTIME_REPO_ROOT);
  pyodide.globals.set("__CLI_WORK_ROOT__", RUNTIME_WORK_ROOT);
  pyodide.globals.set("__CLI_MIN_SCORE__", minScore);
  pyodide.globals.set("__CLI_TASK_TYPE__", taskType);
  pyodide.globals.set("__CLI_TOP_K__", topK);
  pyodide.globals.set("__CLI_ENTITY_TYPE__", entityType);
  pyodide.globals.set("__CLI_MEASUREMENT_CONTEXT__", measurementContext);
  pyodide.globals.set("__CLI_MATRIX_PRIORITY__", matrixPriority);

  try {
    postStatus(`Running CLI-compatible ${taskType === "analyte" ? "map" : "trait-map"} in worker...`);
    const resultJson = await pyodide.runPythonAsync(`
import json
import pathlib
import runpy
import sys
import types

# Pyodide does not provide _multiprocessing, so importing
# concurrent.futures.process fails. Inject a thread-backed shim so
# "from concurrent.futures import ProcessPoolExecutor" still works.
if "concurrent.futures.process" not in sys.modules:
    from concurrent.futures.thread import ThreadPoolExecutor as _ThreadPoolExecutor

    _process_mod = types.ModuleType("concurrent.futures.process")

    class ProcessPoolExecutor(_ThreadPoolExecutor):
        def __init__(self, max_workers=None, mp_context=None, initializer=None, initargs=()):
            # Ignore process-specific kwargs in browser mode.
            super().__init__(max_workers=max_workers)

    class BrokenProcessPool(RuntimeError):
        pass

    _process_mod.ProcessPoolExecutor = ProcessPoolExecutor
    _process_mod.BrokenProcessPool = BrokenProcessPool
    _process_mod.__all__ = ["ProcessPoolExecutor", "BrokenProcessPool"]
    sys.modules["concurrent.futures.process"] = _process_mod

runtime_root = pathlib.Path(__CLI_RUNTIME_ROOT__)
work_root = pathlib.Path(__CLI_WORK_ROOT__)
task_type = str(__CLI_TASK_TYPE__)
script_path = runtime_root / "skills" / "pqtl-measurement-mapper" / "scripts" / "map_measurement_efo.py"
if not script_path.exists():
    raise FileNotFoundError(f"Mapper script not found in runtime bundle: {script_path}")

input_path = work_root / "input.tsv"
output_path = work_root / "mapped.tsv"
review_path = work_root / "review.tsv"
qc_risk_path = work_root / "qc_risk.tsv"
qc_summary_path = work_root / "qc_summary.json"
unmapped_path = work_root / "unmapped.tsv"

trait_cache = runtime_root / "skills" / "pqtl-measurement-mapper" / "references" / "trait_mapping_cache.tsv"
efo_obo = runtime_root / "skills" / "pqtl-measurement-mapper" / "references" / "efo.obo"
ukb_fieldsum = runtime_root / "references" / "ukb" / "fieldsum.txt"
icd10_supp = runtime_root / "references" / "icd10" / "icd10_trait_supplement_cache.tsv"
icd10_labels = runtime_root / "references" / "icd10" / "icd10_label_index.tsv"
measurement_index = runtime_root / "skills" / "pqtl-measurement-mapper" / "references" / "measurement_index.json"

if task_type == "analyte":
    if not measurement_index.exists():
        raise FileNotFoundError(
            "measurement_index.json missing in runtime bundle. Rebuild runtime bundle with "
            "--include-analyte-assets."
        )
    argv = [
        str(script_path),
        "map",
        "--input", str(input_path),
        "--output", str(output_path),
        "--index", str(measurement_index),
        "--unmapped-output", str(unmapped_path),
        "--review-output", str(review_path),
        "--top-k", str(max(1, int(__CLI_TOP_K__))),
        "--min-score", str(float(__CLI_MIN_SCORE__)),
        "--workers", "1",
        "--parallel-mode", "thread",
        "--entity-type", str(__CLI_ENTITY_TYPE__),
        "--measurement-context", str(__CLI_MEASUREMENT_CONTEXT__),
        "--matrix-priority", str(__CLI_MATRIX_PRIORITY__),
        "--name-mode", "strict",
        "--progress",
    ]
else:
    argv = [
        str(script_path),
        "trait-map",
        "--input", str(input_path),
        "--output", str(output_path),
        "--review-output", str(review_path),
        "--qc-risk-output", str(qc_risk_path),
        "--qc-summary-output", str(qc_summary_path),
        "--trait-cache", str(trait_cache),
        "--efo-obo", str(efo_obo),
        "--ukb-field-catalog", str(ukb_fieldsum),
        "--icd10-supplement-cache", str(icd10_supp),
        "--icd10-label-cache", str(icd10_labels),
        "--min-score", str(float(__CLI_MIN_SCORE__)),
        "--stream-output",
        "--memoize-queries",
        "--progress",
    ]

old_argv = sys.argv[:]
try:
    sys.argv = argv
    try:
        runpy.run_path(str(script_path), run_name="__main__")
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 0
        if code not in (0, None):
            raise RuntimeError(f"mapper exited with code {code}")
finally:
    sys.argv = old_argv

if not output_path.exists():
    raise RuntimeError("CLI-compatible mapper did not produce output TSV")

json.dumps(
    {
        "output_text": output_path.read_text(encoding="utf-8"),
        "task_type": task_type,
    }
)
`);
    const parsed = JSON.parse(resultJson);
    self.postMessage({
      type: "mapped",
      payload: {
        outputText: String(parsed.output_text || ""),
        logs: getRunLogsText(),
        engine: "cli_compat",
        taskType: String(parsed.task_type || taskType),
      },
    });
  } finally {
    pyodide.globals.delete("__CLI_RUNTIME_ROOT__");
    pyodide.globals.delete("__CLI_WORK_ROOT__");
    pyodide.globals.delete("__CLI_MIN_SCORE__");
    pyodide.globals.delete("__CLI_TASK_TYPE__");
    pyodide.globals.delete("__CLI_TOP_K__");
    pyodide.globals.delete("__CLI_ENTITY_TYPE__");
    pyodide.globals.delete("__CLI_MEASUREMENT_CONTEXT__");
    pyodide.globals.delete("__CLI_MATRIX_PRIORITY__");
  }
}

async function runMap(payload) {
  const engine = parseEngine(payload);
  if (engine === "cli_compat") {
    await runCliCompatMap(payload);
    return;
  }
  await runFastMap(payload);
}

self.onmessage = async (event) => {
  const type = event?.data?.type;
  const payload = event?.data?.payload ?? {};
  try {
    if (type === "init") {
      await loadPyodideRuntime();
      self.postMessage({ type: "ready", payload: { ok: true } });
      return;
    }
    if (type === "map") {
      await runMap(payload);
      return;
    }
  } catch (error) {
    self.postMessage({
      type: "error",
      payload: { message: String(error?.message || error) },
    });
  }
};
