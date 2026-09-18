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
  $("busy").classList.toggle("is-hidden", !on);
}

async function api(url, options) {
  const response = await fetch(url, options);
  let body;
  try { body = await response.json(); }
  catch { throw new Error(`The server returned an unexpected response (HTTP ${response.status}).`); }
  if (!response.ok) {
    const error = new Error(body.error || `Request failed (HTTP ${response.status}).`);
    Object.assign(error, body);        // keep any partial results the server returned
    throw error;
  }
  return body;
}

function fillSelect(select, values, { placeholder = null } = {}) {
  select.replaceChildren();
  if (placeholder !== null) select.append(new Option(placeholder, ""));
  for (const [value, label] of values) select.append(new Option(label, value));
}

/* ---------- panels and sidebar navigation -------------------------------- */

const navItem = (panelId) => document.querySelector(`.nav-item[data-target="${panelId}"]`);

/** Make a panel reachable in the sidebar without showing it yet. */
function unlock(panelId) {
  $(panelId).classList.remove("is-hidden");
  navItem(panelId)?.classList.remove("is-locked");
}

/** Reveal a panel, mark it current, and scroll to it. */
function reveal(panelId) {
  unlock(panelId);
  setCurrent(panelId);
  $(panelId).scrollIntoView({ behavior: "smooth", block: "start" });
}

function setCurrent(panelId) {
  for (const item of document.querySelectorAll(".nav-item")) {
    item.classList.toggle("is-active", item.dataset.target === panelId);
  }
}

for (const item of document.querySelectorAll(".nav-item")) {
  item.onclick = (event) => {
    event.preventDefault();
    const panel = $(item.dataset.target);
    if (!panel || panel.classList.contains("is-hidden")) return;
    setCurrent(item.dataset.target);
    panel.scrollIntoView({ behavior: "smooth", block: "start" });
  };
}

/* Keep the sidebar in step with whichever panel the reader is looking at. */
const spy = new IntersectionObserver((entries) => {
  const visible = entries.filter((e) => e.isIntersecting)
    .sort((a, b) => b.intersectionRatio - a.intersectionRatio)[0];
  if (visible) setCurrent(visible.target.id);
}, { rootMargin: "-25% 0px -60% 0px", threshold: [0.1, 0.5] });

for (const panel of document.querySelectorAll(".panel")) spy.observe(panel);

/* ---------- right rail ---------------------------------------------------- */

function setDataset(nSamples, nGenes) {
  $("stat-samples").textContent = nSamples.toLocaleString();
  $("stat-genes").textContent = nGenes.toLocaleString();
  $("rail-title").textContent = "Dataset loaded";
  $("rail-sub").textContent = `${nSamples.toLocaleString()} samples ready to compare`;
}

function setOutcome(data) {
  $("stat-up").textContent = data.n_up.toLocaleString();
  $("stat-down").textContent = data.n_down.toLocaleString();
  $("outcome-empty").classList.add("is-hidden");
  $("outcome-figures").classList.remove("is-hidden");
  const note = $("outcome-note");
  note.textContent = `${data.group_a} vs ${data.group_b} · ` +
    `${data.n_genes_tested.toLocaleString()} genes tested`;
  note.classList.remove("is-hidden");
  $("rail-title").textContent = `${data.group_a} vs ${data.group_b}`;
  $("rail-sub").textContent = `${data.n_a} vs ${data.n_b} samples`;
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
    setDataset(data.n_samples, data.n_genes);
    data.warnings.forEach((w) => showMessage(w, "warn"));

    fillSelect($("column-select"), data.columns.map((c) => [c.name, c.name]));
    fillSelect($("batch-select"), data.columns.map((c) => [c.name, c.name]), { placeholder: "None" });
    onColumnChange();
    reveal("panel-setup");
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

  state.stamp = Date.now();
  state.comparisonColumn = data.column;
  $("download-csv").href = `/api/results/${state.session}.csv?t=${state.stamp}`;
  fillPcaColors();
  for (const name of Object.keys(FIGURES)) delete $(name).dataset.loaded;
  showFigure(state.figure || "volcano");

  setOutcome(data);
  reveal("panel-results");
  unlock("panel-pathways");
  unlock("panel-interpret");
}

/* Four views of the same analysis share one tab strip and one download button.
   Each figure is fetched the first time its tab is opened, not all at once. */
const FIGURES = {
  volcano: {
    url: () => "volcano.png?",
    caption: "Fold change against significance. Points far left or right and high up are the " +
      "strongest changes.",
  },
  ma: {
    url: () => "ma.png?",
    caption: "Fold change against mean expression. The black running median should hug zero; " +
      "a curve away from zero at low or high expression suggests a normalization bias.",
  },
  heatmap: {
    url: () => `heatmap.png?genes=${$("heatmap-genes").value}&`,
    control: "heatmap-genes-field",
    caption: "Each gene z-scored across samples and clustered. Shows whether the groups separate " +
      "consistently, sample by sample.",
  },
  pca: {
    url: () => `pca.png?color_by=${encodeURIComponent($("pca-color").value)}&`,
    control: "pca-color-field",
    caption: "Overall similarity between samples, from the 500 most variable genes. Colour by " +
      "another column to see what else drives the variation.",
  },
};

function figureUrl(name) {
  return `/api/figure/${state.session}/${FIGURES[name].url()}t=${state.stamp}`;
}

function showFigure(name) {
  state.figure = name;
  document.querySelectorAll(".seg[data-figure]")
    .forEach((t) => t.classList.toggle("is-active", t.dataset.figure === name));
  for (const [key, figure] of Object.entries(FIGURES)) {
    $(key).hidden = key !== name;
    if (figure.control) $(figure.control).hidden = key !== name;
  }
  const img = $(name);
  if (!img.dataset.loaded) {
    img.src = figureUrl(name);
    img.dataset.loaded = "1";
  }
  $("figure-caption").textContent = FIGURES[name].caption;
  $("download-figure").href = `${figureUrl(name)}&download=1`;
}

function reloadFigure(name) {
  delete $(name).dataset.loaded;
  if (state.figure === name) showFigure(name);
}

/* PCA can be coloured by the comparison or by any other grouping column in the metadata. */
function fillPcaColors() {
  const others = state.columns.map((c) => c.name).filter((n) => n !== state.comparisonColumn);
  fillSelect($("pca-color"), [[state.comparisonColumn, `${state.comparisonColumn} (comparison)`],
                              ...others.map((n) => [n, n])]);
}

for (const tab of document.querySelectorAll(".seg[data-figure]")) {
  tab.onclick = () => showFigure(tab.dataset.figure);
}

$("heatmap-genes").onchange = () => reloadFigure("heatmap");
$("pca-color").onchange = () => reloadFigure("pca");

for (const name of Object.keys(FIGURES)) {
  $(name).onerror = () => {
    if (!$(name).hidden) showMessage("That figure could not be drawn for this comparison.");
  };
}

/* ---------- step 4: pathway enrichment ----------------------------------- */

$("enrich-btn").onclick = async () => {
  clearMessages();
  busy(true, "Asking g:Profiler which pathways are over-represented…");
  try {
    const data = await api("/api/enrich", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ session: state.session, organism: $("organism").value }),
    });
    data.notes.forEach((n) => showMessage(n, "warn"));

    const body = $("pathways-table").querySelector("tbody");
    body.replaceChildren();
    for (const term of data.terms) {
      const tr = document.createElement("tr");
      const link = document.createElement("a");
      link.href = term.url;
      link.target = "_blank";
      link.rel = "noopener";
      link.textContent = term.name;

      const nameCell = document.createElement("td");
      nameCell.append(link);
      tr.append(nameCell);
      tr.insertAdjacentHTML("beforeend",
        `<td>${term.source}</td>` +
        `<td class="num">${term.intersection_size}/${term.term_size}</td>` +
        `<td class="num">${term.p_value.toExponential(1)}</td>` +
        `<td><span class="tag ${term.direction}">${term.direction}</span></td>`);
      body.append(tr);
    }
    $("pathways").classList.toggle("is-hidden", data.terms.length === 0);
    if (data.terms.length) refreshPathwayFigure();
  } catch (error) {
    showMessage(error.message);
  } finally {
    busy(false);
  }
};

/* The enrichment results can be drawn three ways, or read as a table. */
function refreshPathwayFigure() {
  if (state.pathwayView === "table") return;
  const terms = $("pathway-terms").value;
  const url = `/api/figure/${state.session}/pathways.png` +
    `?kind=${state.pathwayView || "dot"}&terms=${terms}&t=${state.stamp}`;
  $("pathway-figure").src = url;
  $("download-pathway-figure").href = `${url}&download=1`;
}

state.pathwayView = "dot";

for (const tab of document.querySelectorAll(".seg[data-pathfig]")) {
  tab.onclick = () => {
    document.querySelectorAll(".seg[data-pathfig]")
      .forEach((t) => t.classList.toggle("is-active", t === tab));
    state.pathwayView = tab.dataset.pathfig;
    const showTable = state.pathwayView === "table";
    $("pathway-table-wrap").classList.toggle("is-hidden", !showTable);
    $("pathway-figure-wrap").classList.toggle("is-hidden", showTable);
    $("download-pathway-figure").classList.toggle("is-hidden", showTable);
    refreshPathwayFigure();
  };
}

$("pathway-terms").onchange = refreshPathwayFigure;

$("pathway-figure").onerror = () => {
  if (state.pathwayView !== "table") showMessage("That pathway figure could not be drawn.");
};

/* ---------- step 5: literature + interpretation -------------------------- */

/* Minimal Markdown for the model's reply: headings, bullets, bold, paragraphs.
   Everything is escaped first, so model output can never inject HTML. */
function renderMarkdown(text) {
  const escape = (s) => s.replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const inline = (s) => escape(s)
    .replace(/\*\*(.+?)\*\*/g, "<strong>$1</strong>")
    .replace(/\bPMID:?\s*(\d{6,9})/g,
      '<a href="https://pubmed.ncbi.nlm.nih.gov/$1/" target="_blank" rel="noopener">PMID $1</a>');

  const html = [];
  let inList = false;
  for (const line of text.split("\n")) {
    const trimmed = line.trim();
    const bullet = trimmed.match(/^[-*]\s+(.*)$/);
    const heading = trimmed.match(/^(#{2,4})\s+(.*)$/);
    if (bullet) {
      if (!inList) { html.push("<ul>"); inList = true; }
      html.push(`<li>${inline(bullet[1])}</li>`);
      continue;
    }
    if (inList) { html.push("</ul>"); inList = false; }
    if (heading) html.push(`<h3>${inline(heading[2])}</h3>`);
    else if (trimmed) html.push(`<p>${inline(trimmed)}</p>`);
  }
  if (inList) html.push("</ul>");
  return html.join("");
}

function renderEvidence(evidence) {
  const list = $("evidence-list");
  list.replaceChildren();
  for (const item of evidence) {
    const block = document.createElement("div");
    block.className = "evidence-gene";
    const heading = document.createElement("h4");
    heading.textContent = `${item.gene} (${item.direction}, log2FC ${item.log2fc})`;
    block.append(heading);

    if (item.articles.length) {
      const ul = document.createElement("ul");
      for (const article of item.articles) {
        const li = document.createElement("li");
        const link = document.createElement("a");
        link.href = `https://pubmed.ncbi.nlm.nih.gov/${article.pmid}/`;
        link.target = "_blank";
        link.rel = "noopener";
        link.textContent = `PMID ${article.pmid}`;
        li.append(link, ` — ${article.title} (${article.journal} ${article.year})`);
        ul.append(li);
      }
      block.append(ul);
    } else {
      const note = document.createElement("p");
      note.className = "hint";
      note.textContent = item.note || "No abstracts retrieved.";
      block.append(note);
    }
    list.append(block);
  }
  $("evidence").classList.remove("is-hidden");
}

$("interpret-btn").onclick = async () => {
  clearMessages();
  busy(true, "Searching PubMed and interpreting the results…");
  try {
    const data = await api("/api/interpret", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        session: state.session,
        disease: $("disease").value,
        n_genes: parseInt($("n-genes").value, 10),
        api_key: $("api-key").value,
      }),
    });
    const box = $("summary");
    box.innerHTML = renderMarkdown(data.summary);
    const footnote = document.createElement("p");
    footnote.className = "hint";
    footnote.textContent = `Generated by ${data.model}. These are hypotheses from the retrieved ` +
      `abstracts, not conclusions — check the sources before relying on them.`;
    box.append(footnote);
    box.classList.remove("is-hidden");
    renderEvidence(data.evidence);
  } catch (error) {
    // The server still returns the PubMed results when only the model call failed
    showMessage(error.message);
    if (error.evidence) renderEvidence(error.evidence);
  } finally {
    busy(false);
  }
};

/* ---------- about ---------------------------------------------------------- */

$("help-btn").onclick = () => {
  clearMessages();
  showMessage("Upload a count matrix and sample metadata, pick two groups, and the app runs " +
    "DESeq2, finds enriched pathways via g:Profiler, and summarises the literature from PubMed. " +
    "Your counts never leave this machine. Demo files are in the repository's demo_data folder.",
    "warn");
};
