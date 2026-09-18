"""Grounded interpretation of the top genes using the OpenAI API.

The model is given only three things: gene symbols with their statistics, the user's disease
term, and the PubMed abstracts already retrieved in rnaseq.ncbi. Counts, sample names and any
other part of the uploaded files are never included - build_payload() is the single place where
the request is assembled, and tests assert it cannot leak them.

The output is explicitly framed as hypotheses: the model must cite the PMIDs it relied on and
say when the retrieved literature does not support a claim.
"""
import os

from openai import OpenAI, OpenAIError

DEFAULT_MODEL = "gpt-5-mini"        # $0.25/M input, $2/M output: a run costs well under a cent

# gpt-5 models spend output tokens on internal reasoning before writing anything, so the cap has
# to cover both. At 1400 the reasoning used the whole budget and the reply came back empty.
MAX_OUTPUT_TOKENS = 4000
REASONING_EFFORT = "low"            # this is summarisation, not a problem that needs deep reasoning

SYSTEM_PROMPT = """\
You are helping a wet-lab researcher interpret a bulk RNA-seq differential expression result.

Ground rules:
1. Use ONLY the gene statistics and the PubMed abstracts provided in this message. Do not add facts
   from memory, and do not invent citations.
2. Cite evidence as (PMID: 12345678). Only cite PMIDs that appear in the provided material.
3. The abstracts came from an automated PubMed search on the gene symbol, so some are irrelevant
   (short symbols such as HP, TF or F2 often match unrelated topics). Silently ignore any abstract
   that is not about the gene in question; never cite it.
4. If the provided abstracts do not support an interpretation for a gene, say so plainly rather
   than speculating.
5. These are correlations from bulk tissue. Bulk RNA-seq measures whatever tissue was in the
   sample, so differences between groups from different organs or tissue sites often reflect
   tissue composition (normal cells captured alongside tumour) rather than tumour biology. Where
   the gene list looks like a tissue signature, say that explicitly.
6. Frame everything as hypotheses to test, never as conclusions.

Write for a researcher who knows the biology but not the statistics. Use GitHub-flavoured Markdown
with these sections and nothing else:

## What stands out
Three to five bullets on the clearest patterns in the gene list and, where pathway enrichment is
provided, the biological processes those genes belong to.

## Possible interpretations
Two to four short paragraphs, each citing PMIDs where the abstracts support the point. Refer to
enriched pathways by name where they help explain the gene list.

## Caveats
Bullets covering tissue composition, sample size, and anything the literature did not cover.
"""


def build_payload(evidence, comparison, disease_term, n_up, n_down, n_tested, pathway_lines=""):
    """Assemble the user message. This is the only text sent to OpenAI."""
    lines = [
        f"Comparison: {comparison} (positive log2 fold change = higher in the first group).",
        f"Disease/context term used for the literature search: {disease_term or 'none given'}.",
        f"{n_tested:,} genes tested; {n_up:,} significantly higher and {n_down:,} significantly "
        f"lower in the first group.",
    ]
    if pathway_lines:
        lines += ["", "PATHWAY ENRICHMENT", pathway_lines]
    lines += ["", "TOP GENES AND THEIR RETRIEVED LITERATURE"]
    for item in evidence:
        lines.append(
            f"\n### {item.gene} ({'higher' if item.direction == 'up' else 'lower'}, "
            f"log2FC {item.log2fc:+.2f}, adjusted p {item.padj:.2e})")
        if not item.articles:
            lines.append(f"  No usable abstracts retrieved. {item.note}".rstrip())
            continue
        for article in item.articles:
            lines.append(f"  - PMID {article.pmid} ({article.year}, {article.journal}): "
                         f"{article.title}\n    {article.truncated()}")
    return "\n".join(lines)


def summarize(evidence, comparison, disease_term, n_up, n_down, n_tested,
              pathway_lines="", api_key=None, model=DEFAULT_MODEL):
    """Return {'summary', 'model', 'payload_chars', 'usage'}; raises RuntimeError with a
    user-facing message when the API cannot be reached."""
    key = (api_key or os.getenv("OPENAI_API_KEY") or "").strip()
    if not key:
        raise RuntimeError("No OpenAI API key was provided. Add OPENAI_API_KEY to your .env file, "
                           "or paste a key into the box above.")
    payload = build_payload(evidence, comparison, disease_term, n_up, n_down, n_tested, pathway_lines)
    try:
        response = OpenAI(api_key=key, timeout=120).responses.create(
            model=model,
            instructions=SYSTEM_PROMPT,
            input=payload,
            max_output_tokens=MAX_OUTPUT_TOKENS,
            reasoning={"effort": REASONING_EFFORT},
        )
    except OpenAIError as exc:
        raise RuntimeError(f"The OpenAI request failed: {_readable(exc)}") from exc

    text = (response.output_text or "").strip()
    if not text:
        reason = getattr(getattr(response, "incomplete_details", None), "reason", None)
        if reason == "max_output_tokens":
            raise RuntimeError("The model ran out of output tokens before writing a summary. "
                               "Try fewer genes per direction.")
        raise RuntimeError("The model returned an empty response. Please try again.")
    usage = getattr(response, "usage", None)
    return {
        "summary": text,
        "model": getattr(response, "model", model),
        "payload_chars": len(payload),
        "usage": {"input_tokens": getattr(usage, "input_tokens", None),
                  "output_tokens": getattr(usage, "output_tokens", None)} if usage else None,
    }


def _readable(exc):
    """Turn SDK errors into something a user can act on, without echoing the key."""
    text = str(exc)
    # A 429 means either "out of credits" or "too many requests", which need different fixes
    if "insufficient_quota" in text or "credit_balance_exhausted" in text:
        return ("the OpenAI account has no credits left. API usage is billed separately from a "
                "ChatGPT subscription - add credits at platform.openai.com under Billing.")
    status = getattr(exc, "status_code", None)
    known = {401: "the API key was rejected. Check it was copied in full.",
             404: "that model is not available to this account.",
             429: "too many requests in a short time. Wait a moment and try again.",
             500: "OpenAI had a server error. Try again in a moment.",
             503: "OpenAI is temporarily unavailable. Try again in a moment."}
    return known.get(status, text.split("\n")[0][:200])
