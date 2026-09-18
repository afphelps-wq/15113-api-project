# Build Plan: Bulk RNA-seq Explorer (Flask)

**Deadline:** Sunday, Sep 20, 2026, 8:00 pm. Today is Friday, Sep 18, so about 2 days.
**Goal:** researchers upload a count matrix and sample metadata, choose two groups, and get differential
expression results, a volcano plot, PubMed evidence for the top genes (NCBI E-utilities), and a grounded
AI summary (OpenAI API). Decisions behind this plan are in `project_context.txt`.

## 1. What I verified before planning

| Check | Result | Effect on plan |
|---|---|---|
| Python | 3.14.2, nothing installed yet | Pydeseq2 0.5.4, Flask 3.1, pandas 3, scipy, openai all install cleanly on 3.14. Use a project venv and pin versions. |
| DE runtime | Full data (289 samples, 22.5k filtered genes): **14 s**. 60-sample subset: about 11 s. | A background job queue is not needed. Run synchronously with a loading spinner. (Changes an earlier answer; revisit only if uploads are much larger.) |
| Gene IDs | `ID` is versioned Ensembl (`ENSG…14`). `hgnc_symbol` is blank for 21,976 of 60,554 rows. 16 duplicate symbols. | Use symbol, else `gene_name`, else Ensembl ID. Collapse duplicates by keeping the highest-expressed row. |
| Counts | Non-integer estimates | Round to int before PyDESeq2 (already done in the benchmark). |
| Sanity check | Primary vs. Met: GCG, ABCC8, CPB1, CTRC lower in Met; HP, HPX, SERPINC1, PLG, F2 higher in Met | The pipeline finds the expected biology. This is mostly tissue composition (normal pancreas vs. liver tissue), so the UI and AI prompt must warn about it. |

## 2. Scope: what ships and what gets cut

The deadline is short, so scope is tiered. Each tier is only started when the one above it is done and tested.

- **Tier 1: Must ship (Fri to Sat).** CSV upload, validation, PyDESeq2, results table, volcano plot, CSV and PNG download, PubMed abstracts through NCBI, grounded OpenAI summary, README, prompt log.
- **Tier 2: Submission tasks (Sun).** Portfolio entry, demo video, Google form, code walkthrough, final secret scan.
- **Tier 3: Stretch, in this order:** (a) g:Profiler enrichment (keyless, about 1 h); (b) QC/PCA plot (about 1.5 h); (c) HTML report download; (d) GEO accession fetch; (e) mouse support.
- **Realistic expectation:** (a) and maybe (b) fit. (c) to (e) probably don't. Those cuts are fine; the assignment grades API integration, README, prompt log, portfolio and video.

## 3. Architecture

```
15113-api-project/
  app.py                  Flask app: routes only, thin
  rnaseq/
    io_utils.py           parse + validate counts/metadata, gene ID cleanup, filtering
    analysis.py           PyDESeq2 wrapper -> results DataFrame
    plots.py              volcano plot (matplotlib -> PNG bytes)
    ncbi.py               PubMed esearch/efetch client (requests)
    llm.py                prompt builder + OpenAI call
  templates/index.html    single-page UI
  static/app.js, style.css
  demo_data/              ~60-sample subset + metadata (committed)
  scripts/make_demo_data.py   documents how the subset came from GSE205154
  tests/                  pytest
  requirements.txt  .env.example  README.md  prompt_log.md  WALKTHROUGH.md
```

**Routes (keys never reach browser JS):**
- `POST /api/upload`: parse and validate, return metadata columns and their levels.
- `POST /api/analyze`: run DE for chosen column, groups, optional batch column, thresholds.
- `GET /api/results/<id>.csv` and `GET /api/figure/<id>/volcano.png`: downloads.
- `POST /api/interpret`: fetch PubMed abstracts for the top genes, then call OpenAI, return the summary.

State is held server-side in memory keyed by a run ID (single local user). Nothing is written to disk except temp files.

**Input format:** counts is genes as rows and samples as columns, first column is the gene ID, and extra non-numeric annotation columns are auto-detected. Metadata's first column matches sample names, and other columns are group variables.

## 4. Milestones (in order, each with an acceptance test)

| # | Milestone | Est. | Done when |
|---|---|---|---|
| M0 | venv, `requirements.txt`, folder skeleton, `make_demo_data.py`, demo subset committed | 1 h | `pip install -r requirements.txt` works; demo files under about 6 MB |
| M1 | `io_utils` + `analysis` + unit tests | 3 h | Demo run reproduces the sanity-check genes; bad inputs give clear error messages |
| M2 | Flask routes + UI: upload, pick column/groups, results table, volcano, downloads | 3 h | Full flow works in the browser on demo data |
| M3 | `ncbi.py`: PubMed search and fetch for up and down genes | 2 h | Returns titles, abstracts and PMIDs for top 10 to 15 genes; handles empty results and errors |
| M4 | `llm.py`: grounded summary with PMID citations | 2 h | Summary cites only supplied PMIDs; works with the key in `.env` and with a pasted key |
| M5 | UI polish, error states, privacy notice in UI, README, prompt log | 2 h | Fresh-clone instructions work end to end |
| Stretch | g:Profiler, then QC/PCA | 1 to 2.5 h | Only if M0 to M5 are done by Saturday night |
| S1 | Portfolio entry, screen-recorded video, Google form, secret scan | 3 h | Video link tested in incognito; form submitted |
| S2 | `WALKTHROUGH.md` (code tour) | 1 h | Each module explained for discussion |

**Suggested schedule:** Friday: M0 to M2. Saturday: M3 to M5 (and stretch if ahead). Sunday: S1, S2 and buffer, finished by about 5 pm so there are three hours of slack before 8 pm.

## 5. API integration details

**NCBI E-utilities (PubMed).** For each top gene, `esearch` (`db=pubmed`, term = gene symbol in title/abstract AND the user's disease term, `retmax` about 3, sorted by date), then one batched `efetch` (`rettype=abstract`, XML) for all PMIDs. Uses `requests`. Limit is 3 requests/second without a key and 10/second with a free `NCBI_API_KEY`. So 12 genes takes about 13 calls. Returns XML, parsed into a list of `{pmid, title, abstract, year}`. Known issue: short symbols like HP, F2 and PLG return noisy hits, so the model is told to ignore irrelevant abstracts.

**OpenAI.** One call per analysis, capped output tokens, cheap model. I will check the current model names and API shape in OpenAI's docs at build time instead of assuming. The key comes from `OPENAI_API_KEY` in `.env`, or from a password field in the UI, held in memory only and never logged.

**Prompt guardrails.**
- Input to the model is limited to: gene symbols, log2FC, adjusted p-values, the disease term, and the fetched abstracts.
- The model is told to cite PMIDs, say "insufficient evidence" where abstracts don't support a claim, and label everything as hypotheses.
- The prompt reminds it that bulk tissue differences can reflect tissue composition.

## 6. Privacy

Raw counts and sample IDs are never sent to any external API. A unit test checks that the LLM payload builder can't include them, and the UI states this plainly.

## 7. Testing

- **Unit tests (pytest):** count rounding, ID/sample mismatch errors, duplicate-symbol collapse, PubMed XML parsing on a saved sample, and the payload privacy check.
- **Sanity test:** on the demo subset, GCG and CTRC are down and HP and HPX are up for Met vs. Primary.
- **Manual:** run the full 289-sample file once, and one deliberately broken upload.

## 8. Submission checklist (Sunday, before 8 pm)

- [ ] Repo public, `.env` never committed; run a secret scan of the full git history
- [ ] README: 3 to 5 sentences on how the APIs are called, key instructions (no key shown), install and run steps
- [ ] `prompt_log.md` lists tools and key prompts
- [ ] Project added to afphelps-wq.github.io
- [ ] Video link works in an incognito window
- [ ] Google form submitted

## 9. Risks

| Risk | Mitigation |
|---|---|
| Flask front end plus all features is a lot for 2 days | Tier order above; UI kept to one page and plain JS |
| No OpenAI API key yet | You create it and set a spend limit before Saturday. Until then M4 is built against a mock response |
| Tissue-composition confounding could mislead the AI summary | Prompt guardrail and a visible UI note |
| Noisy PubMed matches for short symbols | Disease term in query; model ignores irrelevant abstracts |
| Python 3.14 is new | Verified installs work; pin exact versions |
