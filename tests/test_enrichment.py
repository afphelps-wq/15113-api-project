"""Tests for the g:Profiler client. No network: _post is replaced with a canned response."""
import pandas as pd
import pytest

from rnaseq import enrichment
from rnaseq.enrichment import Term, as_prompt_lines, enrich

GPROFILER_RESULT = {"result": [
    {"source": "GO:BP", "native": "GO:0042730", "name": "fibrinolysis", "p_value": 6.2e-9,
     "intersection_size": 5, "term_size": 27},
    {"source": "REAC", "native": "REAC:R-HSA-114608", "name": "Platelet degranulation",
     "p_value": 2.2e-8, "intersection_size": 7, "term_size": 126},
    {"source": "KEGG", "native": "KEGG:04610", "name": "Complement and coagulation cascades",
     "p_value": 1.0e-5, "intersection_size": 4, "term_size": 85},
]}


def table_with(n_up=5, n_down=5, n_ns=10):
    rows = ([("U%d" % i, "up") for i in range(n_up)] +
            [("D%d" % i, "down") for i in range(n_down)] +
            [("N%d" % i, "ns") for i in range(n_ns)])
    return pd.DataFrame(rows, columns=["gene", "regulation"])


@pytest.fixture
def captured(monkeypatch):
    """Record every payload sent to g:Profiler instead of making a request."""
    sent = []

    def fake_post(payload, session):
        sent.append(payload)
        return GPROFILER_RESULT

    monkeypatch.setattr(enrichment, "_post", fake_post)
    return sent


class TestEnrich:
    def test_queries_each_direction_separately(self, captured):
        terms, notes = enrich(table_with())
        assert len(captured) == 2
        assert captured[0]["query"] == ["U0", "U1", "U2", "U3", "U4"]
        assert captured[1]["query"] == ["D0", "D1", "D2", "D3", "D4"]
        assert {t.direction for t in terms} == {"up", "down"}

    def test_sends_all_tested_genes_as_background(self, captured):
        enrich(table_with(n_up=5, n_down=5, n_ns=10))
        assert len(captured[0]["background"]) == 20          # every row, not just the significant ones
        assert captured[0]["domain_scope"] == "custom_annotated"

    def test_skips_a_direction_with_too_few_genes(self, captured):
        terms, notes = enrich(table_with(n_up=2, n_down=5))
        assert len(captured) == 1                             # only the down-regulated query ran
        assert any("Too few up-regulated genes (2)" in n for n in notes)
        assert {t.direction for t in terms} == {"down"}

    def test_organism_name_is_translated(self, captured):
        enrich(table_with(), organism="mouse")
        assert captured[0]["organism"] == "mmusculus"

    def test_unknown_organism_passes_straight_through(self, captured):
        enrich(table_with(), organism="rnorvegicus")
        assert captured[0]["organism"] == "rnorvegicus"

    def test_results_are_sorted_and_capped(self, captured):
        terms, _ = enrich(table_with(), max_terms=2)
        up = [t for t in terms if t.direction == "up"]
        assert len(up) == 2
        assert up[0].p_value < up[1].p_value

    def test_broad_parent_categories_are_dropped(self, monkeypatch):
        """A term covering thousands of genes describes nothing useful, however significant it is."""
        monkeypatch.setattr(enrichment, "_post", lambda payload, session: {"result": [
            {"source": "GO:BP", "native": "GO:0032501", "name": "multicellular organismal process",
             "p_value": 1e-34, "intersection_size": 533, "term_size": 6617},
            {"source": "GO:BP", "native": "GO:0042730", "name": "fibrinolysis",
             "p_value": 6.2e-9, "intersection_size": 5, "term_size": 27},
        ]})
        terms, _ = enrich(table_with())
        assert [t.name for t in terms if t.direction == "up"] == ["fibrinolysis"]

    def test_terms_with_too_few_matching_genes_are_dropped(self, monkeypatch):
        monkeypatch.setattr(enrichment, "_post", lambda payload, session: {"result": [
            {"source": "KEGG", "native": "KEGG:1", "name": "barely there",
             "p_value": 1e-3, "intersection_size": 2, "term_size": 40},
        ]})
        terms, notes = enrich(table_with())
        assert terms == []
        assert any("Only very broad categories" in n for n in notes)

    def test_empty_result_is_reported(self, monkeypatch):
        monkeypatch.setattr(enrichment, "_post", lambda payload, session: {"result": []})
        terms, notes = enrich(table_with())
        assert terms == []
        assert len(notes) == 2 and "No enriched pathways" in notes[0]


class TestTermLinks:
    @pytest.mark.parametrize("term_id,expected", [
        ("GO:0042730", "amigo.geneontology.org"),
        ("REAC:R-HSA-114608", "reactome.org/content/detail/R-HSA-114608"),
        ("KEGG:04610", "kegg.jp/pathway/map04610"),
        ("WP:WP123", "biit.cs.ut.ee"),                        # unknown source falls back to g:Profiler
    ])
    def test_each_source_links_to_its_own_database(self, term_id, expected):
        term = Term("X", term_id, "name", 0.01, 3, 50, "up")
        assert expected in term.url


class TestPromptLines:
    def test_empty_terms_produce_no_text(self):
        assert as_prompt_lines([]) == ""

    def test_lines_label_direction_and_include_statistics(self):
        terms = [Term("GO:BP", "GO:0042730", "fibrinolysis", 6.2e-9, 5, 27, "up"),
                 Term("KEGG", "KEGG:04610", "coagulation", 1e-5, 4, 85, "down")]
        text = as_prompt_lines(terms)
        assert "higher in the first group" in text and "lower in the first group" in text
        assert "fibrinolysis (GO:BP GO:0042730), p=6.2e-09, 5 of 27 genes" in text
