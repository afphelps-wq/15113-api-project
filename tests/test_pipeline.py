"""Tests for input handling and differential expression.

The sanity check at the bottom runs the real demo dataset, so it is slower than the rest.
"""
from pathlib import Path

import pandas as pd
import pytest

from rnaseq.analysis import classify, filter_low_counts, run_deseq
from rnaseq.io_utils import (ValidationError, align_samples, groupable_columns, parse_counts,
                             parse_metadata, read_table)

DEMO = Path(__file__).resolve().parent.parent / "demo_data"


def counts_frame(**columns):
    """Build a raw counts table: gene_id, gene_symbol, then sample columns."""
    n = len(next(iter(columns.values())))
    return pd.DataFrame({"gene_id": [f"ENSG{i:05d}" for i in range(n)],
                         "gene_symbol": [f"GENE{i}" for i in range(n)], **columns})


class TestParseCounts:
    def test_rounds_fractional_counts(self):
        parsed = parse_counts(counts_frame(S1=[613.7, 10.2], S2=[368.4, 9.8]))
        assert parsed.counts.loc["GENE0", "S1"] == 614
        assert parsed.counts.dtypes.eq("int64").all()
        assert any("rounded" in w for w in parsed.warnings)

    def test_keeps_gene_ids_and_labels(self):
        parsed = parse_counts(counts_frame(S1=[5, 6], S2=[7, 8]))
        assert list(parsed.counts.index) == ["GENE0", "GENE1"]
        assert parsed.gene_ids["GENE0"] == "ENSG00000"

    def test_blank_symbol_falls_back_to_id(self):
        raw = counts_frame(S1=[5, 6], S2=[7, 8])
        raw.loc[0, "gene_symbol"] = None
        parsed = parse_counts(raw)
        assert "ENSG00000" in parsed.counts.index

    def test_duplicate_labels_keep_highest_expressed(self):
        raw = counts_frame(S1=[10, 900], S2=[10, 900])
        raw["gene_symbol"] = "SAME"
        parsed = parse_counts(raw)
        assert len(parsed.counts) == 1
        assert parsed.counts.loc["SAME", "S1"] == 900
        assert any("shared a gene label" in w for w in parsed.warnings)

    def test_negative_counts_rejected(self):
        with pytest.raises(ValidationError, match="negative"):
            parse_counts(counts_frame(S1=[-5, 6], S2=[7, 8]))

    def test_missing_values_rejected(self):
        with pytest.raises(ValidationError, match="missing values"):
            parse_counts(counts_frame(S1=[None, 6], S2=[7, 8]))

    def test_too_few_sample_columns_rejected(self):
        with pytest.raises(ValidationError, match="at least two sample columns"):
            parse_counts(pd.DataFrame({"gene_id": ["A"], "gene_symbol": ["B"]}))

    def test_normalized_data_warns(self):
        parsed = parse_counts(counts_frame(S1=[1.2, 0.4], S2=[2.5, 0.9]))
        assert any("normalized" in w for w in parsed.warnings)


class TestParseMetadata:
    def test_finds_sample_column_by_name(self):
        meta = parse_metadata(pd.DataFrame({"tumor_type": ["A", "B"], "sample": ["S1", "S2"]}))
        assert list(meta.index) == ["S1", "S2"]
        assert "tumor_type" in meta.columns

    def test_duplicate_samples_rejected(self):
        with pytest.raises(ValidationError, match="unique"):
            parse_metadata(pd.DataFrame({"sample": ["S1", "S1"], "group": ["a", "b"]}))

    def test_groupable_columns_reports_level_counts(self):
        meta = parse_metadata(pd.DataFrame({"sample": ["S1", "S2", "S3"],
                                            "group": ["a", "a", "b"],
                                            "note": ["x", "y", "z"]}))
        groups = {c["name"]: c["levels"] for c in groupable_columns(meta)}
        assert groups["group"] == {"a": 2, "b": 1}


class TestAlignSamples:
    def test_drops_samples_missing_from_either_table(self):
        counts = parse_counts(counts_frame(S1=[5, 6], S2=[7, 8], S3=[9, 10])).counts
        meta = parse_metadata(pd.DataFrame({"sample": ["S1", "S2", "S9"], "g": ["a", "b", "a"]}))
        aligned_counts, aligned_meta, warnings = align_samples(counts, meta)
        assert list(aligned_counts.columns) == ["S1", "S2"] == list(aligned_meta.index)
        assert len(warnings) == 2

    def test_no_overlap_is_an_error(self):
        counts = parse_counts(counts_frame(S1=[5, 6], S2=[7, 8])).counts
        meta = parse_metadata(pd.DataFrame({"sample": ["X1", "X2"], "g": ["a", "b"]}))
        with pytest.raises(ValidationError, match="No sample names match"):
            align_samples(counts, meta)


class TestFilteringAndClassification:
    def test_filter_requires_reads_in_several_samples(self):
        counts = pd.DataFrame({f"S{i}": [100, 0, 0] for i in range(6)},
                              index=["expressed", "silent", "blip"])
        counts.loc["blip", "S0"] = 500                     # only one sample has reads
        assert list(filter_low_counts(counts).index) == ["expressed"]

    def test_classify_labels_up_down_and_ns(self):
        table = pd.DataFrame({"gene": ["u", "d", "small", "weak"],
                              "log2FoldChange": [3.0, -3.0, 0.2, 5.0],
                              "padj": [0.001, 0.001, 0.001, 0.5]})
        labels = dict(zip(table["gene"], classify(table)["regulation"]))
        assert labels == {"u": "up", "d": "down", "small": "ns", "weak": "ns"}


class TestRunDeseqValidation:
    @pytest.fixture
    def tiny(self):
        counts = parse_counts(counts_frame(**{f"S{i}": [100 + i, 50 + i, 30 + i] for i in range(6)})).counts
        meta = parse_metadata(pd.DataFrame({"sample": [f"S{i}" for i in range(6)],
                                            "g": ["a", "a", "a", "b", "b", "b"],
                                            "batch": ["1", "1", "1", "1", "1", "1"]}))
        return counts, meta

    def test_same_group_twice_rejected(self, tiny):
        with pytest.raises(ValidationError, match="two different groups"):
            run_deseq(*tiny, column="g", group_a="a", group_b="a")

    def test_unknown_column_rejected(self, tiny):
        with pytest.raises(ValidationError, match="not a column"):
            run_deseq(*tiny, column="nope", group_a="a", group_b="b")

    def test_group_with_one_sample_rejected(self):
        counts = parse_counts(counts_frame(**{f"S{i}": [100 + i, 50 + i, 30 + i] for i in range(4)})).counts
        meta = parse_metadata(pd.DataFrame({"sample": [f"S{i}" for i in range(4)],
                                            "g": ["a", "a", "a", "b"]}))
        with pytest.raises(ValidationError, match="at least 2 samples"):
            run_deseq(counts, meta, column="g", group_a="a", group_b="b")


@pytest.mark.slow
class TestDemoDataSanityCheck:
    """The known biology of GSE205154: metastases carry liver/plasma genes, primaries carry pancreas genes."""

    @pytest.fixture(scope="class")
    def result(self):
        counts = parse_counts(read_table(DEMO / "demo_counts.csv")).counts
        meta = parse_metadata(read_table(DEMO / "demo_metadata.csv"))
        counts, meta, _ = align_samples(counts, meta)
        return run_deseq(counts, meta, column="tumor_type", group_a="Met", group_b="Primary")

    def test_groups_and_gene_counts(self, result):
        assert (result.n_a, result.n_b) == (30, 30)
        assert result.n_genes_tested > 10_000

    def test_pancreas_genes_are_lower_in_metastases(self, result):
        table = result.table.set_index("gene")
        for gene in ["GCG", "CTRC", "CPB1"]:                # pancreas-specific
            assert table.loc[gene, "log2FoldChange"] < -1, gene
            assert table.loc[gene, "padj"] < 0.05, gene

    def test_liver_plasma_genes_are_higher_in_metastases(self, result):
        table = result.table.set_index("gene")
        for gene in ["HP", "HPX", "SERPINC1"]:              # liver-secreted plasma proteins
            assert table.loc[gene, "log2FoldChange"] > 1, gene
            assert table.loc[gene, "padj"] < 0.05, gene

    def test_reversing_groups_flips_the_sign(self, result):
        counts = parse_counts(read_table(DEMO / "demo_counts.csv")).counts
        meta = parse_metadata(read_table(DEMO / "demo_metadata.csv"))
        counts, meta, _ = align_samples(counts, meta)
        flipped = run_deseq(counts, meta, column="tumor_type", group_a="Primary", group_b="Met")
        forward = result.table.set_index("gene").loc["GCG", "log2FoldChange"]
        backward = flipped.table.set_index("gene").loc["GCG", "log2FoldChange"]
        assert forward == pytest.approx(-backward, rel=1e-6)
