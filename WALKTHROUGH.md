# Code Walkthrough

A tour of how the app works, following one analysis from uploaded file to AI summary. It explains
the parts that are not obvious from reading the code, and why several choices were made the way they
were.

## The shape of the project

The Flask routes in `app.py` are deliberately thin. Each one validates its input, calls into the
`rnaseq` package, and turns the result into JSON. All of the real work lives in six modules that
know nothing about the web:

| Module | Responsibility |
|---|---|
| `rnaseq/io_utils.py` | Read and validate user files; produce clean counts and metadata |
| `rnaseq/analysis.py` | Filter genes and run the PyDESeq2 comparison |
| `rnaseq/plots.py` | Render the volcano plot and heatmap as PNG bytes |
| `rnaseq/enrichment.py` | Ask g:Profiler which pathways are over-represented |
| `rnaseq/ncbi.py` | Search PubMed and parse the abstracts |
| `rnaseq/llm.py` | Build the prompt and call OpenAI |

The benefit is testability: 62 tests exercise these modules directly, with no HTTP server and no
network. The only network calls in the whole project are in `enrichment.py`, `ncbi.py` and `llm.py`.

## 1. Reading the files (`io_utils.py`)

A count matrix is genes in rows and samples in columns, but real files vary, so `parse_counts()`
works out the structure rather than assuming it:

```python
id_col = df.columns[0]
sample_cols = [c for c in rest if pd.api.types.is_numeric_dtype(df[c])]
annotation_cols = [c for c in rest if c not in sample_cols]
```

The first column is the gene ID, any numeric column is a sample, and anything else is a gene
annotation. That single rule handles the GEO file (which has `ensembl_gene_id`, `hgnc_symbol`,
`gene_biotype` and `gene_name` sitting between the ID and the samples) without needing a
configuration option.

**Counts get rounded.** The GEO file contains values like `613.70206`, because tools such as RSEM
and Salmon estimate counts probabilistically rather than counting whole reads. DESeq2's negative
binomial model requires integers, so the values are rounded and the user is told it happened.

**Gene labels fall back through three options.** About 22,000 of the 60,554 genes in GSE205154 have
no HGNC symbol. The label becomes the symbol, or the descriptive gene name, or finally the Ensembl
ID, so every row keeps a usable identity.

**Duplicate labels are collapsed, order preserved.** Several Ensembl IDs can share one symbol, and
DESeq2 needs unique row names. The most highly expressed row wins:

```python
table = (table.sort_values("_total", ascending=False, kind="stable")
         .drop_duplicates("_label").sort_values("_row", kind="stable"))
```

The final `sort_values("_row")` matters more than it looks. Without it the whole matrix stays sorted
by expression level, so every gene comes back in a different order than it went in. A test caught
this; nothing about the results would have looked obviously wrong.

Everything a user can get wrong raises `ValidationError`, which carries a message written for them
rather than a stack trace. `app.py` registers one handler that turns any such error into a clean
HTTP 400, so no route needs its own error handling.

## 2. Running the statistics (`analysis.py`)

Most genes in a typical matrix are not usably expressed, and testing them costs time while making
the multiple-testing correction harsher. `filter_low_counts()` keeps a gene only if it has at least
10 reads in at least 10% of samples (minimum 3). On the demo data this drops nothing, because the
demo file was already filtered; on the full 60,554-gene series it removes about two thirds.

`run_deseq()` builds the design matrix. The comparison is always two groups from one metadata
column, optionally with a batch covariate, which becomes the formula `~ batch + condition`.

The interesting part is what it refuses to do:

```python
elif ((crosstab > 0).sum(axis=1) < 2).all():
    raise ValidationError(f"'{batch}' is completely confounded with '{column}' ...")
```

If every batch contains only one of the two groups, the batch effect and the biological effect are
mathematically indistinguishable. DESeq2 would either fail cryptically or return nonsense, so the
app explains the problem instead.

Thresholds are deliberately **not** applied inside `run_deseq()`. A separate `classify()` function
labels each gene `up`, `down` or `ns`, which means changing a cutoff does not require refitting the
model.

## 3. Drawing the figures (`plots.py`)

Matplotlib is set to the `Agg` backend at import time, before `pyplot` is imported. Without this,
matplotlib tries to open a GUI window and crashes inside a web server process. Both figures are
returned as PNG bytes and never written to disk.

**One colour language across both figures.** Red always means higher expression and blue always
means lower, in the volcano plot and the heatmap alike. The heatmap's group annotation deliberately
uses two *different* hues (orange and aqua), so "which group is this sample in" can never be
misread as "is this gene up or down". The heatmap's z-score scale is a diverging ramp: two opposite
hues with a neutral grey midpoint, because zero means "average for this gene" and should read as
nothing at all.

### The volcano plot

Two details came directly from looking at the rendered image:

**Infinite values.** Very small adjusted p-values can underflow to exactly zero, and `-log10(0)` is
infinity, which breaks the y-axis. The fix clips them to the smallest non-zero value present.

**Label collisions.** The strongest genes cluster in the corners, so their labels overlap into an
unreadable smear. `_label_points()` walks the ranked genes and skips any whose label would land too
close to one already placed, measuring distance in axes-relative coordinates so the rule holds at any
data scale.

### The heatmap

The volcano plot shows *which* genes changed; the heatmap shows *how consistently*, one sample at a
time. That answers a question the volcano cannot: does this comparison actually separate the two
groups, or is it driven by a handful of samples?

**Each gene is z-scored across samples**, after a log2 transform. Without that, the plot would be
dominated by a few abundant genes and every other row would look flat. The z-score asks a more
useful question: is this sample high or low *for this gene*?

**Genes are clustered, samples are not.** Gene clustering uses correlation distance, which groups by
the shape of the pattern rather than the expression level. Samples stay in group order instead, so
the two conditions form contiguous blocks — which is what makes a clean split visible, or reveals
that there isn't one. On the demo data the result is two sharp blocks: liver and drug-metabolism
genes high in metastases, pancreatic acinar and islet genes high in primaries.

**Constant genes are removed before clustering.** A gene with the same value in every sample has a
standard deviation of zero, so z-scoring divides by zero, and its correlation with anything else is
undefined. That crashed `scipy.cluster.hierarchy.linkage` with "The condensed distance matrix must
contain only finite values." A test now covers it. Such genes carry no information anyway.

**The colour scale is clipped to the 98th percentile** of absolute z-scores rather than the maximum.
One extreme outlier would otherwise compress every other cell toward the neutral midpoint and wash
the figure out.

## 4. Finding the pathways (`enrichment.py`)

A list of 2,900 significant genes is hard to read. Enrichment asks a different question: which
biological processes appear more often in that list than chance would predict?

**Up and down are queried separately.** A pathway enriched among genes that went up means something
quite different from the same pathway going down, so merging them would destroy the signal.

**The background matters more than it looks.** g:Profiler compares your gene list against a
background set, and by default that is every annotated gene in the genome. But an experiment can
only detect genes it actually measured, so the app sends the tested genes as a custom background:

```python
"background": background,
"domain_scope": "custom_annotated",
```

On the demo data this is the difference between 69 "enriched" terms and 31. The larger number is
partly an artifact of comparing against genes the experiment never had a chance to detect.

**Broad terms are filtered out.** The first working version returned this as its top hit:

> multicellular organismal process — 533 of 6,617 genes, p = 1.6e-34

Statistically unarguable and completely useless. With thousands of significant genes, the largest GO
categories are always "enriched". Keeping only terms with at most 500 genes, and at least 3 matches,
replaced that with specific, interpretable biology — drug metabolism, biological oxidations and
complement/coagulation cascades up; pancreatic secretion and protein digestion down. That reads as a
liver-versus-pancreas tissue signature, which is exactly what these samples are.

The enriched terms are also passed into the OpenAI prompt, so the summary can describe processes
rather than reciting gene names.

## 5. Finding the literature (`ncbi.py`)

Two E-utilities calls per analysis cycle: one `esearch` per gene, then a **single** batched `efetch`
for every PubMed ID found. Fetching abstracts one at a time would work, but it multiplies the request
count for no benefit — NCBI explicitly prefers batched requests.

**Search terms are narrowed on purpose:**

```python
term = f'"{gene}"[Title/Abstract]'
if disease_term:
    term += f" AND ({disease_term})"
```

Restricting to title and abstract, and requiring the disease term, keeps results relevant. It is not
perfect: `HP` still matches "hereditary pancreatitis" and `TF` matches "transcription factor". That
residual noise is handled in the prompt, described below.

**Not every label is searchable.** GEO sometimes supplies a description instead of a symbol, such as
"SAGA complex associated factor 29 pseudogene". `is_searchable_symbol()` requires a short,
space-free token, so descriptions and Ensembl IDs are skipped rather than returning junk.

**Rate limiting is built in.** `RateLimiter` spaces requests to just under NCBI's published limits —
3 per second anonymously, 10 with a free key. The app sets 2.5 and 9 to leave headroom.

**Abstract parsing handles structured records.** Many PubMed abstracts are split into labelled
sections, so the parser rejoins them with their labels:

```python
for part in citation.iter("AbstractText"):
    label = part.get("Label")
    body = "".join(part.itertext()).strip()
```

`part.itertext()` rather than `part.text` matters, because titles and abstracts contain nested markup
like `<i>` for species names; reading `.text` alone truncates at the first tag.

**Network failures degrade rather than crash.** If a search fails, that gene gets a note explaining
why and the analysis continues. The literature is a bonus feature; losing it should not lose the
differential expression results.

## 6. Asking the model (`llm.py`)

`build_payload()` is the single place the OpenAI request is assembled. That is a deliberate
constraint: because every piece of outgoing text is built in one function, a test can assert that
sample identifiers and raw counts never appear in it. Funnelling all external communication through
one auditable function is easier to trust than scattering string formatting across the codebase.

The system prompt does most of the work. Three instructions matter most:

1. **Use only the supplied material, and cite PMIDs.** This is what separates a grounded summary from
   a plausible-sounding fabrication.
2. **Ignore irrelevant abstracts.** The prompt names the failure mode explicitly, warning that short
   symbols like `HP`, `TF` and `F2` often match unrelated topics.
3. **Flag tissue composition.** Bulk RNA-seq measures whatever tissue was in the sample, so a
   comparison across organs partly measures the organ. The model is told to say so when the gene list
   looks like a tissue signature.

**The token budget bug.** The first working version returned an empty string every time. The API
response explained it: `status=incomplete, reason=max_output_tokens`, with all 1,400 output tokens
spent and zero characters written. GPT-5 models reason internally before producing visible text, and
that reasoning is billed and capped as output. The fix raises the cap to 4,000 and sets reasoning
effort to `low`, since summarizing supplied text does not need deep reasoning:

```python
MAX_OUTPUT_TOKENS = 4000
REASONING_EFFORT = "low"
```

A typical run uses about 11,600 input and 2,000 output tokens, costing roughly $0.007.

**Errors are translated.** `_readable()` turns SDK exceptions into something actionable. A 429 is
ambiguous — it means either "out of credits" or "too many requests" — so the code inspects the error
body and gives the right advice, without ever echoing the API key.

## 7. Serving it (`app.py`)

Analyses are kept in a module-level dictionary keyed by an unguessable token from
`secrets.token_urlsafe()`, with a lock around mutation and a cap of 8 sessions. This is the right
weight for a single-user local tool: uploaded data never touches disk, and the process holds nothing
after it exits. A multi-user deployment would need a real store, and the README says so.

One route behaviour is worth noting. `/api/interpret` returns the PubMed evidence **even when the
OpenAI call fails**, using HTTP 502 with both an error and the evidence:

```python
except RuntimeError as exc:
    return jsonify({"error": str(exc), "evidence": _evidence_json(genes)}), 502
```

The front end keeps that partial data (`Object.assign(error, body)`) and renders the source list
anyway. A missing API key then costs the user the summary, not the literature search.

## 8. The front end (`static/app.js`)

Plain JavaScript, no framework — the interactivity is four buttons and some `fetch` calls, which does
not justify a build step in a project meant to be read.

**No API key ever reaches the browser's code.** The optional key field posts to the server and is used
for that request only. This is why the app is Flask rather than a static page: a keyed API called from
front-end JavaScript exposes the key to every visitor who opens developer tools.

**Model output is escaped before rendering.** `renderMarkdown()` escapes HTML first, then applies a
small set of Markdown rules and turns PMIDs into PubMed links. Gene names from the user's file go
through `textContent`, never `innerHTML`. Both are untrusted input in the security sense, and a gene
named `<script>` should never execute.

## Where to look first

- **To understand the analysis:** `analysis.py`, then the sanity-check tests in `tests/test_pipeline.py`
- **To understand the API integration:** `ncbi.py` and `llm.py`, then `tests/test_literature.py`
- **To understand the flow:** `app.py` top to bottom; it is under 200 lines
