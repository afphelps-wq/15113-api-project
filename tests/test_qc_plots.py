"""Tests for the MA plot and the sample PCA."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rnaseq.analysis import run_deseq, sample_pca
from rnaseq.io_utils import align_samples, parse_counts, parse_metadata, read_table
from rnaseq.plots import ma_plot, pca_plot

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
DEMO = Path(__file__).resolve().parent.parent / "demo_data"


def de_table(n=400):
    rng = np.random.default_rng(0)
    lfc = rng.normal(0, 1.5, n)
    padj = np.where(np.abs(lfc) > 2, 1e-6, 0.5)
    return pd.DataFrame({
        "gene": [f"G{i}" for i in range(n)],
        "baseMean": rng.lognormal(5, 2, n),
        "log2FoldChange": lfc,
        "padj": padj,
        "regulation": np.where(padj < 0.05, np.where(lfc > 0, "up", "down"), "ns"),
    })


class TestMaPlot:
    def test_renders(self):
        assert ma_plot(de_table(), "A", "B").startswith(PNG_MAGIC)

    def test_genes_with_zero_mean_are_skipped_not_logged(self):
        """log10(0) is -inf; such genes must be dropped rather than break the x-axis."""
        table = de_table()
        table.loc[:5, "baseMean"] = 0.0
        assert ma_plot(table, "A", "B").startswith(PNG_MAGIC)

    def test_small_tables_skip_the_running_median(self):
        assert ma_plot(de_table(n=30), "A", "B").startswith(PNG_MAGIC)

    def test_no_significant_genes(self):
        table = de_table()
        table["regulation"] = "ns"
        assert ma_plot(table, "A", "B").startswith(PNG_MAGIC)


class TestPcaPlot:
    @pytest.fixture
    def pca(self):
        rng = np.random.default_rng(1)
        samples = [f"S{i}" for i in range(12)]
        return pd.DataFrame(rng.normal(0, 5, (12, 2)), index=samples, columns=["PC1", "PC2"])

    def test_renders_with_two_groups(self, pca):
        labels = pd.Series(["a"] * 6 + ["b"] * 6, index=pca.index)
        assert pca_plot(pca, (0.4, 0.2), labels, "a vs b").startswith(PNG_MAGIC)

    def test_many_categories_fold_into_other(self, pca):
        """Scatter plots cap at three colours; the rest become 'Other' instead of new hues."""
        labels = pd.Series([f"tissue{i % 6}" for i in range(12)], index=pca.index)
        assert pca_plot(pca, (0.4, 0.2), labels, "by tissue").startswith(PNG_MAGIC)

    def test_samples_missing_a_label_are_still_drawn(self, pca):
        labels = pd.Series(["a"] * 6, index=pca.index[:6])        # half the samples unlabelled
        assert pca_plot(pca, (0.4, 0.2), labels, "partial").startswith(PNG_MAGIC)


class _FailingVst:
    """Stands in for a DeseqDataSet whose variance-stabilizing fit fails."""
    layers = {}

    def vst(self, use_design=False):
        raise RuntimeError("fit failed")


class TestSamplePca:
    @pytest.fixture
    def counts(self):
        rng = np.random.default_rng(2)
        samples = [f"S{i}" for i in range(8)]
        return pd.DataFrame(rng.poisson(100, (8, 50)), index=samples,
                            columns=[f"G{i}" for i in range(50)])

    def test_falls_back_to_log_counts_when_vst_fails(self, counts):
        pca, variance, method = sample_pca(_FailingVst(), counts, counts.astype(float))
        assert method.startswith("log2 normalized")
        assert pca.shape == (8, 2)

    def test_variance_fractions_are_ordered_and_bounded(self, counts):
        _, variance, _ = sample_pca(_FailingVst(), counts, counts.astype(float))
        assert 0 <= variance[1] <= variance[0] <= 1

    def test_uses_at_most_the_requested_number_of_genes(self, counts):
        _, _, method = sample_pca(_FailingVst(), counts, counts.astype(float), top_n=20)
        assert "top 20 most variable" in method


@pytest.mark.slow
class TestDemoPca:
    """On GSE205154, liver metastases should separate from everything else on PC1."""

    @pytest.fixture(scope="class")
    def result_and_meta(self):
        counts = parse_counts(read_table(DEMO / "demo_counts.csv")).counts
        meta = parse_metadata(read_table(DEMO / "demo_metadata.csv"))
        counts, meta, _ = align_samples(counts, meta)
        return run_deseq(counts, meta, "tumor_type", "Met", "Primary"), meta

    def test_every_sample_gets_coordinates(self, result_and_meta):
        result, _ = result_and_meta
        assert result.pca.shape == (60, 2)
        assert np.isfinite(result.pca.to_numpy()).all()
        assert result.pca_method.startswith("variance-stabilized")

    def test_liver_samples_separate_from_pancreas_on_pc1(self, result_and_meta):
        result, meta = result_and_meta
        tissue = meta["tissue"].astype(str).reindex(result.pca.index)
        liver = result.pca.loc[tissue == "Liver", "PC1"].median()
        pancreas = result.pca.loc[tissue == "Pancreas", "PC1"].median()
        spread = result.pca["PC1"].std()
        assert abs(liver - pancreas) > spread               # far apart relative to the spread

    def test_non_liver_metastases_sit_with_the_primaries(self, result_and_meta):
        """The finding the PCA exists to show: the Met/Primary split is mostly liver tissue."""
        result, meta = result_and_meta
        tissue = meta["tissue"].astype(str).reindex(result.pca.index)
        pc1 = result.pca["PC1"]
        liver = pc1[tissue == "Liver"].median()
        pancreas = pc1[tissue == "Pancreas"].median()
        other_mets = pc1[~tissue.isin(["Liver", "Pancreas"])].median()
        assert abs(other_mets - pancreas) < abs(other_mets - liver)
