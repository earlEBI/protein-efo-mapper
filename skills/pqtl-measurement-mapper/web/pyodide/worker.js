/* global loadPyodide */

const PYODIDE_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.3/full/";
const CORE_PATH = "./mapper_core.py";
const DEFAULT_CACHE_PATH = "./cache_sample.json";

let pyodideReadyPromise = null;
let pyodide = null;
let defaultCacheText = "";

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

    postStatus("Loading mapper core...");
    const [coreSource, cacheResponse] = await Promise.all([
      fetch(CORE_PATH).then((r) => {
        if (!r.ok) {
          throw new Error(`Failed to load ${CORE_PATH}`);
        }
        return r.text();
      }),
      fetch(DEFAULT_CACHE_PATH),
    ]);

    if (!cacheResponse.ok) {
      throw new Error(`Failed to load ${DEFAULT_CACHE_PATH}`);
    }
    defaultCacheText = await cacheResponse.text();

    await pyodide.runPythonAsync(coreSource);
    postStatus("Pyodide mapper ready.");
  })();

  return pyodideReadyPromise;
}

async function runMap(payload) {
  await loadPyodideRuntime();

  const inputText = String(payload?.inputText ?? "");
  const minScore = Number(payload?.minScore ?? 0.82);
  const cacheText = payload?.cacheText ? String(payload.cacheText) : defaultCacheText;

  pyodide.globals.set("__BROWSER_INPUT__", inputText);
  pyodide.globals.set("__BROWSER_CACHE__", cacheText);
  pyodide.globals.set("__BROWSER_MIN_SCORE__", minScore);

  try {
    postStatus("Running trait mapping in worker...");
    const resultText = await pyodide.runPythonAsync(
      "map_tsv(__BROWSER_INPUT__, __BROWSER_CACHE__, float(__BROWSER_MIN_SCORE__))"
    );
    self.postMessage({
      type: "mapped",
      payload: {
        outputText: resultText,
      },
    });
  } finally {
    pyodide.globals.delete("__BROWSER_INPUT__");
    pyodide.globals.delete("__BROWSER_CACHE__");
    pyodide.globals.delete("__BROWSER_MIN_SCORE__");
  }
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
