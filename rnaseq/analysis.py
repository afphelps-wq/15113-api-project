"""Differential expression with PyDESeq2 (a Python port of R's DESeq2)."""
import warnings as pywarnings
from dataclasses import dataclass, field

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from .io_utils import ValidationError


HEATMAP_GENE_POOL = 300            # normalized counts kept for plotting, most significant first


@dataclass
class DEResult:
    table: pd.DataFrame            # one row per tested gene; log2FoldChange is group_a vs group_b
    normalized: pd.DataFrame       # genes x samples, size-factor normalized, top genes only
    conditions: pd.Series          # sample -> group label, for annotating plots
    column: str
    group_a: str                   # numerator: positive log2FC means higher in this group
    group_b: str                   # reference group
    batch: str | None
    n_a: int
    n_b: int
    n_genes_input: int
    n_genes_tested: int
    warnings: list = field(default_factory=list)


def filter_low_counts(counts, min_count=10):
    """Drop genes that are too sparse to test (genes x samples in, genes x samples out).

    A gene must have at least `min_count` reads in at least 10% of samples (minimum 3).
    """
    min_samples = max(3, int(0.1 * counts.shape[1]))
    return counts[(counts >= min_count).sum(axis=1) >= min_samples]


def run_deseq(counts, meta, column, group_a, group_b, batch=None, min_count=10, n_cpus=4):
    """Compare group_a with group_b within `meta[column]`, optionally adjusting for a batch column.

    counts is genes x samples (as returned by io_utils.parse_counts); meta is indexed by sample.
    Thresholds are not applied here: use classify() so they can change without refitting.
    """
    notes = []
    if column not in meta.columns:
        raise ValidationError(f"'{column}' is not a column in the metadata.")
    if str(group_a) == str(group_b):
        raise ValidationError("Choose two different groups to compare.")
    if batch == column:
        raise ValidationError("The batch column must be different from the comparison column.")

    condition = meta[column].astype("string")
    selected = condition[condition.isin([str(group_a), str(group_b)])].index
    selected = [s for s in counts.columns if s in set(selected)]
    design_meta = pd.DataFrame({"condition": condition.loc[selected].astype(str)}, index=selected)
    n_a = int((design_meta["condition"] == str(group_a)).sum())
    n_b = int((design_meta["condition"] == str(group_b)).sum())
    if min(n_a, n_b) < 2:
        raise ValidationError(f"Each group needs at least 2 samples ('{group_a}': {n_a}, '{group_b}': {n_b}).")
    if min(n_a, n_b) < 3:
        notes.append("A group has only 2 samples, so statistical power is very low.")

    formula = "~ condition"
    if batch:
        if batch not in meta.columns:
            raise ValidationError(f"'{batch}' is not a column in the metadata.")
        batch_values = meta[batch].loc[selected]
        if batch_values.isna().any():
            raise ValidationError(f"The batch column '{batch}' has missing values for "
                                  f"{int(batch_values.isna().sum())} of the compared samples.")
        design_meta["batch"] = batch_values.astype(str)
        crosstab = pd.crosstab(design_meta["batch"], design_meta["condition"])
        if design_meta["batch"].nunique() < 2:
            notes.append(f"Batch column '{batch}' has a single value, so it was not used.")
            design_meta = design_meta.drop(columns="batch")
        elif ((crosstab > 0).sum(axis=1) < 2).all():
            raise ValidationError(f"'{batch}' is completely confounded with '{column}' (each batch "
                                  "contains only one group), so the effect can't be separated.")
        else:
            formula = "~ batch + condition"

    sub_counts = counts[selected].T                       # PyDESeq2 wants samples x genes
    n_input = sub_counts.shape[1]
    sub_counts = filter_low_counts(sub_counts.T, min_count).T
    if sub_counts.shape[1] < 10:
        raise ValidationError("Fewer than 10 genes have enough reads to test. "
                              "Check that the file contains raw counts.")

    try:
        with pywarnings.catch_warnings():
            pywarnings.simplefilter("ignore")
            dds = DeseqDataSet(counts=sub_counts, metadata=design_meta, design=formula,
                               quiet=True, n_cpus=n_cpus)
            dds.deseq2()
            stats = DeseqStats(dds, contrast=["condition", str(group_a), str(group_b)],
                               quiet=True, n_cpus=n_cpus)
            stats.summary()
    except Exception as exc:
        raise ValidationError(f"The statistical model could not be fitted: {exc}") from exc

    table = (stats.results_df.rename_axis("gene").reset_index()
             .sort_values("padj", na_position="last", kind="stable").reset_index(drop=True))

    # Size-factor normalized counts make samples comparable to each other; keep only the most
    # significant genes, since that is all the heatmap can show and the full matrix is large.
    normalized = pd.DataFrame(dds.layers["normed_counts"],
                              index=sub_counts.index, columns=sub_counts.columns)
    keep = [g for g in table["gene"].head(HEATMAP_GENE_POOL) if g in normalized.columns]
    normalized = normalized[keep].T                       # back to genes x samples for plotting

    return DEResult(table=table, normalized=normalized,
                    conditions=design_meta["condition"],
                    column=column, group_a=str(group_a), group_b=str(group_b),
                    batch=batch if "batch" in design_meta else None, n_a=n_a, n_b=n_b,
                    n_genes_input=n_input, n_genes_tested=sub_counts.shape[1], warnings=notes)


def classify(table, padj_cutoff=0.05, lfc_cutoff=1.0):
    """Add a 'regulation' column: 'up', 'down' or 'ns' (not significant)."""
    out = table.copy()
    significant = out["padj"] < padj_cutoff
    out["regulation"] = "ns"
    out.loc[significant & (out["log2FoldChange"] >= lfc_cutoff), "regulation"] = "up"
    out.loc[significant & (out["log2FoldChange"] <= -lfc_cutoff), "regulation"] = "down"
    return out
