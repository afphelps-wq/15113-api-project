"""Flask web app: upload counts, run differential expression, and interpret the top genes.

Routes are deliberately thin; the analysis lives in the rnaseq package. Uploaded data stays in
this process's memory and is never written to disk or sent to an external service. Only gene
symbols and summary statistics leave the server (see rnaseq/llm.py).
"""
import io
import secrets
import threading

from dotenv import load_dotenv
from flask import Flask, jsonify, request, send_file, render_template

from rnaseq import enrichment, llm, ncbi, plots
from rnaseq.analysis import classify, run_deseq
from rnaseq.io_utils import (ValidationError, align_samples, groupable_columns, parse_counts,
                             parse_metadata, read_table)

load_dotenv()                                   # reads .env if present; keys are never logged

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024        # 200 MB upload ceiling

# Analyses held in memory, keyed by an unguessable id. Single-user local app, so a dict is enough.
SESSIONS = {}
SESSIONS_LOCK = threading.Lock()
MAX_SESSIONS = 8


def store(payload):
    token = secrets.token_urlsafe(12)
    with SESSIONS_LOCK:
        SESSIONS[token] = payload
        for old in list(SESSIONS)[:-MAX_SESSIONS]:          # keep memory bounded
            del SESSIONS[old]
    return token


def load(token):
    with SESSIONS_LOCK:
        session = SESSIONS.get(token)
    if session is None:
        raise ValidationError("That analysis has expired. Please upload your files again.")
    return session


@app.errorhandler(ValidationError)
def handle_validation_error(exc):
    return jsonify({"error": str(exc)}), 400


@app.errorhandler(413)
def handle_too_large(_):
    return jsonify({"error": "That file is larger than the 200 MB limit."}), 413


@app.route("/")
def index():
    return render_template("index.html")


@app.post("/api/upload")
def upload():
    """Validate both files and report the grouping columns the user can compare."""
    if "counts" not in request.files or "metadata" not in request.files:
        raise ValidationError("Please choose both a counts file and a metadata file.")

    parsed = parse_counts(read_table(request.files["counts"]))
    meta = parse_metadata(read_table(request.files["metadata"]))
    counts, meta, align_warnings = align_samples(parsed.counts, meta)

    columns = groupable_columns(meta)
    if not columns:
        raise ValidationError("No metadata column has between 2 and 20 distinct values, so there is "
                              "nothing to compare. Check that the metadata describes sample groups.")

    token = store({"counts": counts, "meta": meta, "gene_ids": parsed.gene_ids})
    return jsonify({
        "session": token,
        "n_samples": int(counts.shape[1]),
        "n_genes": int(counts.shape[0]),
        "columns": columns,
        "warnings": parsed.warnings + align_warnings,
    })


@app.post("/api/analyze")
def analyze():
    """Run the two-group comparison. Takes about 15 seconds on 289 samples."""
    body = request.get_json(silent=True) or {}
    session = load(body.get("session"))
    padj_cutoff = float(body.get("padj", 0.05))
    lfc_cutoff = abs(float(body.get("lfc", 1.0)))

    result = run_deseq(session["counts"], session["meta"],
                       column=body.get("column"), group_a=body.get("group_a"),
                       group_b=body.get("group_b"), batch=body.get("batch") or None)
    table = classify(result.table, padj_cutoff, lfc_cutoff)
    session.update({"result": result, "table": table,
                    "cutoffs": {"padj": padj_cutoff, "lfc": lfc_cutoff}})

    counts = table["regulation"].value_counts()
    significant = table[table["regulation"] != "ns"]
    display = significant.head(50) if len(significant) else table.head(50)
    return jsonify({
        "comparison": f"{result.group_a} vs {result.group_b}",
        "column": result.column,
        "group_a": result.group_a, "group_b": result.group_b,
        "n_a": result.n_a, "n_b": result.n_b,
        "batch": result.batch,
        "n_genes_input": result.n_genes_input, "n_genes_tested": result.n_genes_tested,
        "n_up": int(counts.get("up", 0)), "n_down": int(counts.get("down", 0)),
        "warnings": result.warnings,
        "top": [
            {"gene": row.gene,
             "log2FoldChange": None if row.log2FoldChange != row.log2FoldChange else round(row.log2FoldChange, 3),
             "padj": None if row.padj != row.padj else float(f"{row.padj:.3g}"),
             "baseMean": round(row.baseMean, 1),
             "regulation": row.regulation}
            for row in display.itertuples()
        ],
    })


@app.post("/api/enrich")
def enrich():
    """Pathway and GO enrichment for the significant genes, via g:Profiler (no key needed)."""
    body = request.get_json(silent=True) or {}
    session = load(body.get("session"))
    if "table" not in session:
        raise ValidationError("Run the analysis before looking for pathways.")

    try:
        terms, notes = enrichment.enrich(session["table"], organism=body.get("organism", "human"))
    except enrichment.EnrichmentError as exc:
        return jsonify({"error": str(exc)}), 502

    session["pathways"] = terms
    return jsonify({
        "notes": notes,
        "terms": [{"source": t.source, "id": t.term_id, "name": t.name, "p_value": t.p_value,
                   "intersection_size": t.intersection_size, "term_size": t.term_size,
                   "direction": t.direction, "url": t.url}
                  for t in terms],
    })


@app.post("/api/interpret")
def interpret():
    """Look up PubMed evidence for the top genes, then ask the model to interpret it.

    Only gene symbols, their statistics and the retrieved abstracts are sent to OpenAI.
    """
    body = request.get_json(silent=True) or {}
    session = load(body.get("session"))
    if "table" not in session:
        raise ValidationError("Run the analysis before asking for an interpretation.")
    result = session["result"]
    disease = (body.get("disease") or "").strip()
    per_direction = max(1, min(int(body.get("n_genes", 7)), 15))

    genes = ncbi.select_top_genes(session["table"], n_each=per_direction)
    if not genes:
        raise ValidationError("No significant genes with recognisable symbols were found, so there "
                              "is nothing to look up. Try relaxing the thresholds.")
    genes = ncbi.gather_evidence(genes, disease, per_gene=int(body.get("per_gene", 3)))

    counts = session["table"]["regulation"].value_counts()
    try:
        summary = llm.summarize(
            genes, f"{result.group_a} vs {result.group_b}", disease,
            n_up=int(counts.get("up", 0)), n_down=int(counts.get("down", 0)),
            n_tested=result.n_genes_tested,
            pathway_lines=enrichment.as_prompt_lines(session.get("pathways", [])),
            api_key=(body.get("api_key") or "").strip() or None)
    except RuntimeError as exc:
        # The literature is still useful even when the model call fails, so return it either way
        return jsonify({"error": str(exc), "evidence": _evidence_json(genes)}), 502

    return jsonify({**summary, "evidence": _evidence_json(genes)})


def _evidence_json(genes):
    return [{"gene": g.gene, "direction": g.direction, "log2fc": g.log2fc, "padj": g.padj,
             "note": g.note,
             "articles": [{"pmid": a.pmid, "title": a.title, "journal": a.journal, "year": a.year}
                          for a in g.articles]}
            for g in genes]


@app.get("/api/results/<token>.csv")
def download_results(token):
    session = load(token)
    if "table" not in session:
        raise ValidationError("Run the analysis before downloading results.")
    result = session["result"]
    csv_bytes = session["table"].to_csv(index=False).encode()
    name = f"DE_{result.group_a}_vs_{result.group_b}.csv".replace(" ", "_")
    return send_file(_as_stream(csv_bytes), mimetype="text/csv",
                     as_attachment=True, download_name=name)


@app.get("/api/figure/<token>/volcano.png")
def download_volcano(token):
    session = load(token)
    if "table" not in session:
        raise ValidationError("Run the analysis before downloading the figure.")
    result, cutoffs = session["result"], session["cutoffs"]
    png = plots.volcano(session["table"], result.group_a, result.group_b,
                        cutoffs["padj"], cutoffs["lfc"])
    name = f"volcano_{result.group_a}_vs_{result.group_b}.png".replace(" ", "_")
    return send_file(_as_stream(png), mimetype="image/png",
                     as_attachment=request.args.get("download") == "1", download_name=name)


@app.get("/api/figure/<token>/ma.png")
def download_ma(token):
    session = load(token)
    if "table" not in session:
        raise ValidationError("Run the analysis before downloading the figure.")
    result, cutoffs = session["result"], session["cutoffs"]
    png = plots.ma_plot(session["table"], result.group_a, result.group_b, cutoffs["lfc"])
    name = f"MA_{result.group_a}_vs_{result.group_b}.png".replace(" ", "_")
    return send_file(_as_stream(png), mimetype="image/png",
                     as_attachment=request.args.get("download") == "1", download_name=name)


@app.get("/api/figure/<token>/pca.png")
def download_pca(token):
    """PCA of the compared samples, coloured by the comparison or by any metadata column."""
    session = load(token)
    if "table" not in session:
        raise ValidationError("Run the analysis before downloading the figure.")
    result = session["result"]
    color_by = request.args.get("color_by") or result.column
    if color_by == result.column:
        labels, heading = result.conditions, f"{result.group_a} vs {result.group_b}"
    elif color_by in session["meta"].columns:
        labels, heading = session["meta"][color_by].astype("string"), f"Coloured by {color_by}"
    else:
        raise ValidationError(f"'{color_by}' is not a column in the metadata.")
    png = plots.pca_plot(result.pca, result.pca_variance, labels, heading, result.pca_method)
    name = f"PCA_{result.group_a}_vs_{result.group_b}_by_{color_by}.png".replace(" ", "_")
    return send_file(_as_stream(png), mimetype="image/png",
                     as_attachment=request.args.get("download") == "1", download_name=name)


PATHWAY_FIGURES = {"dot": plots.enrichment_dot, "bar": plots.enrichment_bar,
                   "network": plots.enrichment_network}


@app.get("/api/figure/<token>/pathways.png")
def download_pathway_figure(token):
    """Dot plot, bar plot or enrichment network for the terms found in /api/enrich."""
    session = load(token)
    terms = session.get("pathways")
    if not terms:
        raise ValidationError("Find the enriched pathways before drawing this figure.")
    kind = request.args.get("kind", "dot")
    if kind not in PATHWAY_FIGURES:
        raise ValidationError(f"Unknown figure type '{kind}'.")
    result = session["result"]
    top_n = max(3, min(int(request.args.get("terms", 10)), 20))
    png = PATHWAY_FIGURES[kind](terms, result.group_a, result.group_b, top_n=top_n)
    name = f"{kind}_{result.group_a}_vs_{result.group_b}.png".replace(" ", "_")
    return send_file(_as_stream(png), mimetype="image/png",
                     as_attachment=request.args.get("download") == "1", download_name=name)


@app.get("/api/figure/<token>/heatmap.png")
def download_heatmap(token):
    session = load(token)
    if "table" not in session:
        raise ValidationError("Run the analysis before downloading the figure.")
    result = session["result"]
    top_n = max(5, min(int(request.args.get("genes", 30)), 100))
    significant = session["table"][session["table"]["regulation"] != "ns"]
    try:
        png = plots.heatmap(result.normalized, significant if len(significant) >= 5 else session["table"],
                            result.conditions, result.group_a, result.group_b, top_n=top_n)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc
    name = f"heatmap_{result.group_a}_vs_{result.group_b}.png".replace(" ", "_")
    return send_file(_as_stream(png), mimetype="image/png",
                     as_attachment=request.args.get("download") == "1", download_name=name)


def _as_stream(data):
    return io.BytesIO(data)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
