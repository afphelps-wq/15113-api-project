/* Front end: collects inputs, calls the Flask API, renders results.
   No API keys are ever handled here - all external calls happen on the server. */
"use strict";

const $ = (id) => document.getElementById(id);
const state = { session: null, columns: [] };

/* ---------- small helpers ------------------------------------------------ */

function showMessage(text, kind = "error") {
  const box = document.createElement("div");
  box.className = `message ${kind}`;
  box.innerHTML = `<span></span><button aria-label="Dismiss">&times;</button>`;
  box.firstChild.textContent = text;
  box.querySelector("button").onclick = () => box.remove();
  $("messages").append(box);
}

function clearMessages() { $("messages").replaceChildren(); }

function busy(on, text = "Working…") {
  $("busy-text").textContent = text;
  $("busy").classList.toggle("hidden", !on);
}

async function api(url, options) {
  const response = await fetch(url, options);
  let body;
  try { body = await response.json(); }
  catch { throw new Error(`The server returned an unexpected response (HTTP ${response.status}).`); }
  if (!response.ok) throw new Error(body.error || `Request failed (HTTP ${response.status}).`);
  return body;
}

function fillSelect(select, values, { placeholder = null } = {}) {
  select.replaceChildren();
  if (placeholder !== null) select.append(new Option(placeholder, ""));
  for (const [value, label] of values) select.append(new Option(label, value));
}

/* ---------- step 1: upload ----------------------------------------------- */

$("upload-btn").onclick = async () => {
  const counts = $("counts-file").files[0];
  const metadata = $("metadata-file").files[0];
  clearMessages();
  if (!counts || !metadata) {
    showMessage("Please choose both a counts file and a metadata file.");
    return;
  }

  const form = new FormData();
  form.append("counts", counts);
  form.append("metadata", metadata);

  busy(true, "Reading and checking your files…");
  try {
    const data = await api("/api/upload", { method: "POST", body: form });
    state.session = data.session;
    state.columns = data.columns;

    $("dataset-summary").innerHTML =
      `<strong>${data.n_samples.toLocaleString()}</strong> samples and ` +
      `<strong>${data.n_genes.toLocaleString()}</strong> genes loaded.`;
    data.warnings.forEach((w) => showMessage(w, "warn"));

    fillSelect($("column-select"), data.columns.map((c) => [c.name, c.name]));
    fillSelect($("batch-select"), data.columns.map((c) => [c.name, c.name]), { placeholder: "None" });
    onColumnChange();
    $("step-setup").classList.remove("hidden");
    $("step-setup").scrollIntoView({ behavior: "smooth", block: "start" });
  } catch (error) {
    showMessage(error.message);
  } finally {
    busy(false);
  }
};

/* ---------- step 2: pick groups ------------------------------------------ */

function onColumnChange() {
  const column = state.columns.find((c) => c.name === $("column-select").value);
  if (!column) return;
  const levels = Object.entries(column.levels).map(([name, n]) => [name, `${name} (${n})`]);
  fillSelect($("group-a"), levels);
  fillSelect($("group-b"), levels);
  $("group-b").selectedIndex = Math.min(1, levels.length - 1);
}

$("column-select").onchange = onColumnChange;

$("analyze-btn").onclick = async () => {
  clearMessages();
  const payload = {
    session: state.session,
    column: $("column-select").value,
    group_a: $("group-a").value,
    group_b: $("group-b").value,
    batch: $("batch-select").value,
    padj: parseFloat($("padj").value),
    lfc: parseFloat($("lfc").value),
  };
  if (payload.group_a === payload.group_b) {
    showMessage("Choose two different groups to compare.");
    return;
  }

  busy(true, "Fitting the model — this usually takes 10–30 seconds…");
  try {
    renderResults(await api("/api/analyze", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    }));
  } catch (error) {
    showMessage(error.message);
  } finally {
    busy(false);
  }
};

/* ---------- step 3: results ---------------------------------------------- */

function renderResults(data) {
  data.warnings.forEach((w) => showMessage(w, "warn"));

  const batchNote = data.batch ? `, adjusted for <strong>${data.batch}</strong>` : "";
  $("result-summary").innerHTML =
    `<strong>${data.group_a}</strong> (${data.n_a} samples) vs ` +
    `<strong>${data.group_b}</strong> (${data.n_b} samples)${batchNote}. ` +
    `Tested ${data.n_genes_tested.toLocaleString()} of ${data.n_genes_input.toLocaleString()} genes ` +
    `(the rest had too few reads). ` +
    `<strong>${data.n_up.toLocaleString()}</strong> higher and ` +
    `<strong>${data.n_down.toLocaleString()}</strong> lower in ${data.group_a}.`;

  const body = $("results-table").querySelector("tbody");
  body.replaceChildren();
  for (const row of data.top) {
    const tr = document.createElement("tr");
    const tag = row.regulation === "ns" ? "" :
      `<span class="tag ${row.regulation}">${row.regulation}</span>`;
    tr.innerHTML =
      `<td></td><td class="num">${row.log2FoldChange ?? "—"}</td>` +
      `<td class="num">${row.padj ?? "—"}</td>` +
      `<td class="num">${row.baseMean.toLocaleString()}</td><td>${tag}</td>`;
    tr.firstChild.textContent = row.gene;   // gene names come from the user's file: never as HTML
    body.append(tr);
  }

  const cacheBuster = Date.now();
  $("volcano").src = `/api/figure/${state.session}/volcano.png?t=${cacheBuster}`;
  $("download-png").href = `/api/figure/${state.session}/volcano.png?download=1&t=${cacheBuster}`;
  $("download-csv").href = `/api/results/${state.session}.csv?t=${cacheBuster}`;

  $("step-results").classList.remove("hidden");
  $("step-results").scrollIntoView({ behavior: "smooth", block: "start" });
}
