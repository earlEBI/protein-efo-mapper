/* global loadPyodide */

const PYODIDE_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.3/full/";
const CORE_PATH = "./mapper_core.py";
const DEFAULT_CACHE_PATH = "./cache_sample.json";
const DEFAULT_RUNTIME_BUNDLE_PATH = "./trait_runtime_bundle.zip";
const RUNTIME_ZIP_FS_PATH = "/runtime/trait_runtime_bundle.zip";
const RUNTIME_REPO_ROOT = "/runtime/repo";
const RUNTIME_WORK_ROOT = "/runtime/work";
const SUPPORTED_ENGINES = new Set(["fast", "cli_compat"]);

let pyodideReadyPromise = null;
let pyodide = null;
let defaultCacheText = "";
let runtimeBundleReady = false;
let runtimeBundleKey = "";

function postStatus(message) {
  self.postMessage({ type: "status", payload: { message } });
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
        logs: "",
        engine: "fast",
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

  const inputText = String(payload?.inputText ?? "");
  const minScore = Number(payload?.minScore ?? 0.82);
  pyodide.FS.mkdirTree(RUNTIME_WORK_ROOT);
  pyodide.FS.writeFile(`${RUNTIME_WORK_ROOT}/input.tsv`, inputText, { encoding: "utf8" });

  pyodide.globals.set("__CLI_RUNTIME_ROOT__", RUNTIME_REPO_ROOT);
  pyodide.globals.set("__CLI_WORK_ROOT__", RUNTIME_WORK_ROOT);
  pyodide.globals.set("__CLI_MIN_SCORE__", minScore);

  try {
    postStatus("Running CLI-compatible trait-map in worker...");
    const resultJson = await pyodide.runPythonAsync(`
import contextlib
import io
import json
import pathlib
import runpy
import sys

runtime_root = pathlib.Path(__CLI_RUNTIME_ROOT__)
work_root = pathlib.Path(__CLI_WORK_ROOT__)
script_path = runtime_root / "skills" / "pqtl-measurement-mapper" / "scripts" / "map_measurement_efo.py"
if not script_path.exists():
    raise FileNotFoundError(f"Mapper script not found in runtime bundle: {script_path}")

input_path = work_root / "input.tsv"
output_path = work_root / "mapped.tsv"
review_path = work_root / "review.tsv"
qc_risk_path = work_root / "qc_risk.tsv"
qc_summary_path = work_root / "qc_summary.json"

trait_cache = runtime_root / "skills" / "pqtl-measurement-mapper" / "references" / "trait_mapping_cache.tsv"
efo_obo = runtime_root / "skills" / "pqtl-measurement-mapper" / "references" / "efo.obo"
ukb_fieldsum = runtime_root / "references" / "ukb" / "fieldsum.txt"
icd10_supp = runtime_root / "references" / "icd10" / "icd10_trait_supplement_cache.tsv"
icd10_labels = runtime_root / "references" / "icd10" / "icd10_label_index.tsv"

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

log_buffer = io.StringIO()
old_argv = sys.argv[:]
try:
    sys.argv = argv
    with contextlib.redirect_stdout(log_buffer), contextlib.redirect_stderr(log_buffer):
        try:
            runpy.run_path(str(script_path), run_name="__main__")
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 0
            if code not in (0, None):
                raise RuntimeError(f"trait-map exited with code {code}")
finally:
    sys.argv = old_argv

if not output_path.exists():
    raise RuntimeError("CLI-compatible trait-map did not produce output TSV")

json.dumps(
    {
        "output_text": output_path.read_text(encoding="utf-8"),
        "logs": log_buffer.getvalue()[-50000:],
    }
)
`);
    const parsed = JSON.parse(resultJson);
    self.postMessage({
      type: "mapped",
      payload: {
        outputText: String(parsed.output_text || ""),
        logs: String(parsed.logs || ""),
        engine: "cli_compat",
      },
    });
  } finally {
    pyodide.globals.delete("__CLI_RUNTIME_ROOT__");
    pyodide.globals.delete("__CLI_WORK_ROOT__");
    pyodide.globals.delete("__CLI_MIN_SCORE__");
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
