"""Pathway and GO enrichment through the g:Profiler API (g:GOSt).

No API key is needed. One POST per direction (genes higher and genes lower in the first group),
because a pathway enriched among up-regulated genes means something different from the same pathway
among down-regulated ones.

Only gene symbols are sent - the same privacy rule as the PubMed and OpenAI calls.
"""
from dataclasses import dataclass

import requests

ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
TOOL_NAME = "bulk-rnaseq-explorer"
TIMEOUT = 90

# GO biological process, plus two curated pathway databases. Deliberately excludes noisier
# sources (transcription factor motifs, miRNA targets) that crowd out interpretable results.
DEFAULT_SOURCES = ["GO:BP", "KEGG", "REAC"]

ORGANISMS = {"human": "hsapiens", "mouse": "mmusculus"}
MAX_QUERY_GENES = 1000          # g:Profiler slows down on very long lists
MAX_TERMS_PER_DIRECTION = 15

# With a few thousand significant genes, the most significant terms are enormous parent categories
# ("multicellular organismal process", 6,617 genes) that describe nothing useful. Keeping terms
# below a size ceiling surfaces specific biology instead.
MAX_TERM_SIZE = 500
MIN_INTERSECTION = 3


@dataclass
class Term:
    source: str                 # "GO:BP", "KEGG" or "REAC"
    term_id: str                # e.g. "GO:0042730"
    name: str
    p_value: float
    intersection_size: int      # how many of our genes are in this term
    term_size: int              # how many genes the term contains in total
    direction: str              # "up" or "down"

    @property
    def url(self):
        if self.term_id.startswith("GO:"):
            return f"https://amigo.geneontology.org/amigo/term/{self.term_id}"
        if self.term_id.startswith("REAC:"):
            return f"https://reactome.org/content/detail/{self.term_id.removeprefix('REAC:')}"
        if self.term_id.startswith("KEGG:"):
            return f"https://www.kegg.jp/pathway/map{self.term_id.removeprefix('KEGG:')}"
        return f"https://biit.cs.ut.ee/gprofiler/gost?query={self.term_id}"


class EnrichmentError(RuntimeError):
    """g:Profiler could not be reached or rejected the request; message is user-facing."""


def _post(payload, session):
    try:
        response = session.post(ENDPOINT, json=payload, timeout=TIMEOUT,
                                headers={"User-Agent": TOOL_NAME})
        response.raise_for_status()
        return response.json()
    except requests.Timeout as exc:
        raise EnrichmentError("g:Profiler took too long to respond. Try again in a moment.") from exc
    except requests.RequestException as exc:
        raise EnrichmentError(f"Could not reach g:Profiler ({type(exc).__name__}).") from exc
    except ValueError as exc:
        raise EnrichmentError("g:Profiler returned a response that could not be read.") from exc


def enrich(table, organism="human", sources=None, max_terms=MAX_TERMS_PER_DIRECTION):
    """Run enrichment separately for up- and down-regulated genes.

    `table` is the classified DE table. Every tested gene is sent as the statistical background,
    so terms are judged against what the experiment could actually have detected rather than
    against the whole genome.
    """
    species = ORGANISMS.get(organism, organism)
    background = [str(g) for g in table["gene"].tolist()]
    session = requests.Session()
    results, notes = [], []

    for direction in ("up", "down"):
        genes = table.loc[table["regulation"] == direction, "gene"].astype(str).tolist()
        if len(genes) < 3:
            notes.append(f"Too few {direction}-regulated genes ({len(genes)}) to test for enrichment.")
            continue
        payload = {
            "organism": species,
            "query": genes[:MAX_QUERY_GENES],
            "background": background,
            "domain_scope": "custom_annotated",     # background limited to genes we actually tested
            "sources": sources or DEFAULT_SOURCES,
            "user_threshold": 0.05,
            "significance_threshold_method": "g_SCS",   # g:Profiler's own multiple-testing correction
            "no_evidences": True,                        # smaller response; we only need the terms
        }
        found = _post(payload, session).get("result", [])
        if not found:
            notes.append(f"No enriched pathways found among the {direction}-regulated genes.")
        specific = [t for t in found
                    if t["term_size"] <= MAX_TERM_SIZE and t["intersection_size"] >= MIN_INTERSECTION]
        if found and not specific:
            notes.append(f"Only very broad categories were enriched among the {direction}-regulated "
                         "genes, which usually means the gene list is too long to be specific.")
        for item in sorted(specific, key=lambda x: x["p_value"])[:max_terms]:
            results.append(Term(source=item["source"], term_id=item["native"], name=item["name"],
                                p_value=float(item["p_value"]),
                                intersection_size=int(item["intersection_size"]),
                                term_size=int(item["term_size"]), direction=direction))
    return results, notes


def as_prompt_lines(terms, limit_per_direction=8):
    """Compact text for the OpenAI prompt: the strongest terms in each direction."""
    if not terms:
        return ""
    lines = []
    for direction in ("up", "down"):
        subset = [t for t in terms if t.direction == direction][:limit_per_direction]
        if not subset:
            continue
        lines.append(f"\nEnriched pathways among genes {'higher' if direction == 'up' else 'lower'} "
                     f"in the first group (g:Profiler, corrected p-value):")
        lines.extend(f"  - {t.name} ({t.source} {t.term_id}), p={t.p_value:.1e}, "
                     f"{t.intersection_size} of {t.term_size} genes" for t in subset)
    return "\n".join(lines)
