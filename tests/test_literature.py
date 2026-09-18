"""Tests for the PubMed client and the OpenAI prompt builder.

No network calls: the PubMed parser runs against a saved response, and the prompt builder is
pure text assembly. The privacy test is the important one - it asserts that nothing from the
uploaded data can reach OpenAI.
"""
import pandas as pd
import pytest

from rnaseq.llm import build_payload
from rnaseq.ncbi import Article, GeneEvidence, is_searchable_symbol, parse_pubmed_xml, select_top_genes

# Trimmed from a real efetch response: one structured abstract, one plain, one with no abstract.
PUBMED_XML = b"""<?xml version="1.0"?>
<PubmedArticleSet>
 <PubmedArticle><MedlineCitation>
  <PMID Version="1">42468715</PMID>
  <Article>
   <Journal><ISOAbbrev/><ISOAbbreviation>Cell Signal</ISOAbbreviation>
    <JournalIssue><PubDate><Year>2026</Year></PubDate></JournalIssue></Journal>
   <ArticleTitle>Tissue factor drives <i>PDAC</i> invasion</ArticleTitle>
   <Abstract>
    <AbstractText Label="BACKGROUND">Thrombosis is common.</AbstractText>
    <AbstractText Label="RESULTS">TF promotes invasion.</AbstractText>
   </Abstract>
  </Article>
 </MedlineCitation></PubmedArticle>
 <PubmedArticle><MedlineCitation>
  <PMID Version="1">12345678</PMID>
  <Article>
   <Journal><ISOAbbreviation>Pancreas</ISOAbbreviation>
    <JournalIssue><PubDate><MedlineDate>2019 Spring</MedlineDate></PubDate></JournalIssue></Journal>
   <ArticleTitle>A plain abstract</ArticleTitle>
   <Abstract><AbstractText>One single block of text.</AbstractText></Abstract>
  </Article>
 </MedlineCitation></PubmedArticle>
 <PubmedArticle><MedlineCitation>
  <PMID Version="1">99999999</PMID>
  <Article>
   <Journal><ISOAbbreviation>Nature</ISOAbbreviation>
    <JournalIssue><PubDate><Year>2020</Year></PubDate></JournalIssue></Journal>
   <ArticleTitle>Editorial with no abstract</ArticleTitle>
  </Article>
 </MedlineCitation></PubmedArticle>
</PubmedArticleSet>"""


class TestPubMedParsing:
    @pytest.fixture(scope="class")
    def articles(self):
        return parse_pubmed_xml(PUBMED_XML)

    def test_parses_every_record(self, articles):
        assert set(articles) == {"42468715", "12345678", "99999999"}

    def test_joins_labelled_abstract_sections(self, articles):
        assert articles["42468715"].abstract == "BACKGROUND: Thrombosis is common. RESULTS: TF promotes invasion."

    def test_keeps_text_inside_markup_in_titles(self, articles):
        assert articles["42468715"].title == "Tissue factor drives PDAC invasion"

    def test_falls_back_to_medline_date_for_year(self, articles):
        assert articles["12345678"].year == "2019"

    def test_record_without_abstract_is_kept_but_empty(self, articles):
        assert articles["99999999"].abstract == ""

    def test_truncation_keeps_whole_words(self):
        article = Article("1", "t", "alpha beta gamma delta", "J", "2020")
        assert article.truncated(limit=12).endswith("…")
        assert "delta" not in article.truncated(limit=12)


class TestSymbolFiltering:
    @pytest.mark.parametrize("label,expected", [
        ("HP", True), ("SERPINC1", True), ("GCG", True),
        ("ENSG00000257017.8", False),                       # Ensembl IDs are too long
        ("SAGA complex associated factor 29 pseudogene", False),   # a description, not a symbol
        ("", False),
    ])
    def test_only_real_symbols_are_searchable(self, label, expected):
        assert is_searchable_symbol(label) is expected

    def test_select_top_genes_balances_directions_and_skips_bad_symbols(self):
        table = pd.DataFrame({
            "gene": ["A1", "A2", "A3", "B1", "B2", "long description of a pseudogene", "NS1"],
            "log2FoldChange": [5, 4, 3, -5, -4, -3, 0.1],
            "padj": [1e-9, 1e-8, 1e-7, 1e-9, 1e-8, 1e-7, 0.9],
            "regulation": ["up"] * 3 + ["down"] * 3 + ["ns"],
        })
        chosen = select_top_genes(table, n_each=2)
        assert [(g.gene, g.direction) for g in chosen] == [
            ("A1", "up"), ("A2", "up"), ("B1", "down"), ("B2", "down")]


class TestPromptPrivacy:
    """The payload is the only text sent to OpenAI, so it must not contain patient data."""

    @pytest.fixture
    def payload(self):
        evidence = [
            GeneEvidence("GCG", "down", -9.6, 1e-50,
                         articles=[Article("111", "Glucagon in pancreas", "Abstract text.", "J", "2024")]),
            GeneEvidence("HP", "up", 8.7, 1e-56, note="No PubMed articles matched this gene."),
        ]
        return build_payload(evidence, "Met vs Primary", "pancreatic cancer",
                             n_up=1563, n_down=1340, n_tested=22587)

    def test_contains_genes_stats_and_abstracts(self, payload):
        assert "GCG" in payload and "-9.60" in payload and "Abstract text." in payload
        assert "Met vs Primary" in payload and "pancreatic cancer" in payload

    def test_excludes_sample_identifiers_and_raw_counts(self, payload):
        for forbidden in ["ST-00018525-T", "ST-000", "baseMean", "lfcSE", "counts"]:
            assert forbidden not in payload

    def test_reports_genes_with_no_literature(self, payload):
        assert "No PubMed articles matched this gene." in payload

    def test_payload_size_stays_reasonable(self, payload):
        assert len(payload) < 4000            # ~1k tokens for two genes; 14 genes stays well under the cap
