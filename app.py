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

from rnaseq import plots
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


def _as_stream(data):
    return io.BytesIO(data)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
