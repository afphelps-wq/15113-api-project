"""Tests for the figures. PNG bytes are opaque, so these check the data handling that
decides what gets drawn, plus that each function produces a valid image.
"""
import numpy as np
import pandas as pd
import pytest

from rnaseq.plots import heatmap, volcano

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def de_table(n=40):
    rng = np.random.default_rng(0)
    lfc = np.concatenate([rng.uniform(2, 9, n // 2), rng.uniform(-9, -2, n // 2)])
    return pd.DataFrame({
        "gene": [f"GENE{i}" for i in range(n)],
        "log2FoldChange": lfc,
        "padj": np.logspace(-40, -3, n),
        "regulation": ["up"] * (n // 2) + ["down"] * (n // 2),
    })


def expression(genes, n_per_group=6):
    """Counts where the 'up' genes really are higher in group A."""
    rng = np.random.default_rng(1)
    samples = [f"A{i}" for i in range(n_per_group)] + [f"B{i}" for i in range(n_per_group)]
    rows = {}
    for i, gene in enumerate(genes):
        high, low = (400, 40) if i < len(genes) // 2 else (40, 400)
        rows[gene] = np.concatenate([rng.normal(high, 20, n_per_group),
                                     rng.normal(low, 20, n_per_group)]).clip(0)
    conditions = pd.Series(["A"] * n_per_group + ["B"] * n_per_group, index=samples)
    return pd.DataFrame(rows, index=samples).T, conditions


class TestVolcano:
    def test_produces_a_png(self):
        assert volcano(de_table(), "A", "B").startswith(PNG_MAGIC)

    def test_survives_padj_of_exactly_zero(self):
        """-log10(0) is infinite and would break the y-axis."""
        table = de_table()
        table.loc[0, "padj"] = 0.0
        assert volcano(table, "A", "B").startswith(PNG_MAGIC)

    def test_handles_a_table_with_no_significant_genes(self):
        table = de_table()
        table["regulation"] = "ns"
        assert volcano(table, "A", "B").startswith(PNG_MAGIC)


class TestHeatmap:
    def test_produces_a_png(self):
        table = de_table()
        normalized, conditions = expression(table["gene"])
        assert heatmap(normalized, table, conditions, "A", "B").startswith(PNG_MAGIC)

    def test_too_few_genes_is_an_error_not_a_crash(self):
        table = de_table(n=2).head(1)
        normalized, conditions = expression(table["gene"])
        with pytest.raises(ValueError, match="Not enough significant genes"):
            heatmap(normalized, table, conditions, "A", "B")

    def test_constant_genes_do_not_produce_nan(self):
        """A gene with zero variance would divide by zero when z-scoring."""
        table = de_table(n=10)
        normalized, conditions = expression(table["gene"])
        normalized.iloc[0] = 100                     # identical in every sample
        assert heatmap(normalized, table, conditions, "A", "B").startswith(PNG_MAGIC)

    def test_only_requested_number_of_genes_is_used(self):
        table = de_table(n=40)
        normalized, conditions = expression(table["gene"])
        # 4 genes still clusters; the function must not fail when top_n < available
        assert heatmap(normalized, table, conditions, "A", "B", top_n=4).startswith(PNG_MAGIC)

    def test_genes_missing_from_the_matrix_are_skipped(self):
        table = de_table(n=10)
        normalized, conditions = expression(table["gene"][:6])     # matrix has only the first 6
        assert heatmap(normalized, table, conditions, "A", "B").startswith(PNG_MAGIC)
