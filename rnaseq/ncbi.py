"""PubMed literature lookup through NCBI's E-utilities.

Two calls per request cycle:
  esearch  - one per gene, returning PubMed IDs for "<gene> AND <disease>"
  efetch   - a single batched call returning the abstracts for every ID found

Only gene symbols and the user's disease term are sent. Counts and sample names never leave
the machine. NCBI asks that clients identify themselves with `tool` and `email`, and limits
anonymous use to 3 requests/second (10/second with a free API key).
"""
import os
import re
import time
from dataclasses import dataclass, field
from xml.etree import ElementTree

import requests

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "bulk-rnaseq-explorer"
CONTACT_EMAIL = os.getenv("NCBI_EMAIL", "")
TIMEOUT = 20

# A gene symbol is short and has no spaces. GEO sometimes supplies a description instead
# (e.g. "SAGA complex associated factor 29 pseudogene"), which would return junk from PubMed.
SYMBOL_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9._-]{0,14}$")


@dataclass
class Article:
    pmid: str
    title: str
    abstract: str
    journal: str
    year: str

    def truncated(self, limit=1200):
        text = self.abstract
        return text if len(text) <= limit else text[:limit].rsplit(" ", 1)[0] + "…"


@dataclass
class GeneEvidence:
    gene: str
    direction: str                      # "up" or "down"
    log2fc: float
    padj: float
    articles: list = field(default_factory=list)
    note: str = ""                      # set when a gene was skipped or returned nothing


class RateLimiter:
    """Spaces out requests so we stay under NCBI's published limit."""

    def __init__(self, per_second):
        self._interval = 1.0 / per_second
        self._last = 0.0

    def wait(self):
        gap = time.monotonic() - self._last
        if gap < self._interval:
            time.sleep(self._interval - gap)
        self._last = time.monotonic()


def is_searchable_symbol(label):
    return bool(SYMBOL_PATTERN.match(str(label).strip()))


def _params(api_key, **extra):
    params = {"tool": TOOL_NAME, **extra}
    if CONTACT_EMAIL:
        params["email"] = CONTACT_EMAIL
    if api_key:
        params["api_key"] = api_key
    return params


def _text(node, path):
    """Joined text of a node, including nested markup such as <i> inside a title."""
    found = node.find(path)
    return "".join(found.itertext()).strip() if found is not None else ""


def search_gene(gene, disease_term, session, limiter, api_key, retmax=3):
    """Return PubMed IDs for a gene, most recent first.

    The gene is restricted to title/abstract and combined with the disease term, because short
    symbols (HP, F2, TF) match unrelated text elsewhere in a record.
    """
    term = f'"{gene}"[Title/Abstract]'
    if disease_term:
        term += f" AND ({disease_term})"
    limiter.wait()
    response = session.get(f"{BASE}/esearch.fcgi", timeout=TIMEOUT, params=_params(
        api_key, db="pubmed", term=term, retmax=retmax, sort="pub_date", retmode="json"))
    response.raise_for_status()
    return response.json().get("esearchresult", {}).get("idlist", [])


def fetch_articles(pmids, session, limiter, api_key):
    """Fetch many abstracts in one call. Returns {pmid: Article}."""
    if not pmids:
        return {}
    limiter.wait()
    response = session.post(f"{BASE}/efetch.fcgi", timeout=TIMEOUT, data=_params(
        api_key, db="pubmed", id=",".join(pmids), rettype="abstract", retmode="xml"))
    response.raise_for_status()
    return parse_pubmed_xml(response.content)


def parse_pubmed_xml(xml_bytes):
    """Parse an efetch PubMed response into {pmid: Article}."""
    articles = {}
    root = ElementTree.fromstring(xml_bytes)
    for citation in root.iter("MedlineCitation"):
        pmid = _text(citation, "PMID")
        if not pmid:
            continue
        # An abstract may be split into labelled sections (BACKGROUND, METHODS, ...)
        pieces = []
        for part in citation.iter("AbstractText"):
            label = part.get("Label")
            body = "".join(part.itertext()).strip()
            if body:
                pieces.append(f"{label}: {body}" if label else body)
        articles[pmid] = Article(
            pmid=pmid,
            title=_text(citation, "Article/ArticleTitle"),
            abstract=" ".join(pieces),
            journal=_text(citation, "Article/Journal/ISOAbbreviation"),
            year=_text(citation, "Article/Journal/JournalIssue/PubDate/Year")
                 or _text(citation, "Article/Journal/JournalIssue/PubDate/MedlineDate")[:4],
        )
    return articles


def select_top_genes(table, n_each=7):
    """Pick the most significant up- and down-regulated genes with usable symbols."""
    significant = table[table["regulation"] != "ns"].dropna(subset=["padj"])
    chosen = []
    for direction in ("up", "down"):
        subset = significant[significant["regulation"] == direction]
        subset = subset[subset["gene"].map(is_searchable_symbol)]
        for row in subset.head(n_each).itertuples():
            chosen.append(GeneEvidence(gene=row.gene, direction=direction,
                                       log2fc=round(row.log2FoldChange, 3), padj=float(row.padj)))
    return chosen


def gather_evidence(genes, disease_term, api_key=None, per_gene=3):
    """Look up literature for each gene. Network failures degrade to an empty article list."""
    api_key = api_key or os.getenv("NCBI_API_KEY") or None
    limiter = RateLimiter(9 if api_key else 2.5)          # just under NCBI's 10/s and 3/s limits
    session = requests.Session()
    session.headers["User-Agent"] = TOOL_NAME

    wanted = {}
    for evidence in genes:
        try:
            pmids = search_gene(evidence.gene, disease_term, session, limiter, api_key, per_gene)
        except requests.RequestException as exc:
            evidence.note = f"PubMed search failed ({type(exc).__name__})."
            continue
        if not pmids:
            evidence.note = "No PubMed articles matched this gene with the chosen disease term."
        for pmid in pmids:
            wanted.setdefault(pmid, []).append(evidence)

    if wanted:
        try:
            articles = fetch_articles(list(wanted), session, limiter, api_key)
        except (requests.RequestException, ElementTree.ParseError) as exc:
            articles = {}
            for evidence in genes:
                if not evidence.note:
                    evidence.note = f"Could not download abstracts ({type(exc).__name__})."
        for pmid, holders in wanted.items():
            article = articles.get(pmid)
            if article and article.abstract:
                for evidence in holders:
                    evidence.articles.append(article)
    return genes
