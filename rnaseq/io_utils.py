"""Read and validate user-supplied count matrices and sample metadata.

Every problem a user could plausibly cause raises ValidationError with a plain-language message,
so the web app can show it directly instead of a stack trace.
"""
import io
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


class ValidationError(ValueError):
    """The uploaded data can't be analysed; the message is safe to show to the user."""


# Columns we prefer as the displayed gene label, in priority order (case-insensitive)
LABEL_COLUMNS = ["gene_symbol", "symbol", "hgnc_symbol", "gene_name", "gene"]
SAMPLE_ID_COLUMNS = ["sample", "sample_id", "sampleid", "sample_name", "id"]
MAX_GROUP_LEVELS = 20      # a column with more distinct values than this isn't a useful grouping


@dataclass
class ParsedCounts:
    counts: pd.DataFrame          # genes x samples, integer counts, indexed by gene label
    gene_ids: pd.Series           # gene label -> original ID from the first column
    warnings: list = field(default_factory=list)


def read_table(source):
    """Read a CSV or TSV from a path, bytes, or file-like object (e.g. a Flask upload)."""
    if isinstance(source, (str, Path)):
        raw = Path(source).read_bytes()
    elif hasattr(source, "read"):
        raw = source.read()
    else:
        raw = bytes(source)
    if not raw.strip():
        raise ValidationError("The file is empty.")
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    first_line = text.split("\n", 1)[0]
    sep = "\t" if first_line.count("\t") > first_line.count(",") else ","
    try:
        return pd.read_csv(io.StringIO(text), sep=sep, low_memory=False)
    except Exception as exc:                       # pandas raises many parser-specific errors
        raise ValidationError(f"Could not read the file as a CSV/TSV table: {exc}") from exc


def _pick_column(columns, candidates):
    lowered = {str(c).strip().lower(): c for c in columns}
    for name in candidates:
        if name in lowered:
            return lowered[name]
    return None


def parse_counts(df):
    """Turn a raw table (genes as rows, samples as columns) into clean integer counts.

    The first column is the gene ID. Other non-numeric columns are treated as gene annotations.
    """
    warnings = []
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 3:
        raise ValidationError("The counts file needs a gene ID column plus at least two sample columns.")

    id_col = df.columns[0]
    rest = df.columns[1:]
    sample_cols = [c for c in rest if pd.api.types.is_numeric_dtype(df[c])]
    annotation_cols = [c for c in rest if c not in sample_cols]
    if len(sample_cols) < 2:
        raise ValidationError(
            "Fewer than two numeric sample columns were found. Genes should be rows and samples "
            "should be columns, with the gene ID in the first column.")
    if annotation_cols:
        warnings.append(
            f"Ignored non-numeric columns: {', '.join(annotation_cols[:5])}"
            f"{'...' if len(annotation_cols) > 5 else ''}. If any of these are samples, check that "
            "the column has only numbers (no 'NA' text).")

    values = df[sample_cols]
    if values.isna().any().any():
        bad = values.columns[values.isna().any()][0]
        raise ValidationError(f"The counts contain missing values (first seen in column '{bad}').")
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise ValidationError("The counts contain infinite values.")
    if (values < 0).any().any():
        raise ValidationError("The counts contain negative numbers; raw counts must be zero or higher.")

    rounded = values.round().astype("int64")
    fractional = bool((values != rounded).any().any())
    if fractional:
        warnings.append("Counts were not whole numbers (e.g. RSEM or Salmon estimates), "
                        "so they were rounded to integers.")
        if float(values.max().max()) < 50:
            warnings.append("Values are small and fractional, which suggests normalized or "
                            "log-transformed data. DESeq2 needs raw counts, so results may be wrong.")

    # Gene label: preferred annotation column, falling back to the ID for blanks
    label_col = _pick_column(annotation_cols, LABEL_COLUMNS) or id_col
    ids = df[id_col].astype(str)
    labels = df[label_col].astype("string").str.strip()
    labels = labels.mask(labels.isna() | (labels == ""), ids).astype(str)

    table = rounded.copy()
    table["_total"] = table.sum(axis=1)
    table["_label"] = labels.to_numpy()
    table["_id"] = ids.to_numpy()
    table["_row"] = range(len(table))
    n_before = len(table)
    # Keep the most highly expressed row per label, then restore the file's original gene order
    table = (table.sort_values("_total", ascending=False, kind="stable")
             .drop_duplicates("_label").sort_values("_row", kind="stable"))
    if len(table) < n_before:
        warnings.append(f"{n_before - len(table)} rows shared a gene label with another row; "
                        "kept the most highly expressed row for each.")
    gene_ids = pd.Series(table["_id"].to_numpy(), index=table["_label"].to_numpy(), name="gene_id")
    counts = table.drop(columns=["_total", "_label", "_id", "_row"])
    counts.index = table["_label"].to_numpy()
    counts.index.name = "gene"
    return ParsedCounts(counts=counts, gene_ids=gene_ids, warnings=warnings)


def parse_metadata(df):
    """Return sample metadata indexed by sample name. The sample column is found by name or is first."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 2:
        raise ValidationError("The metadata file needs a sample column plus at least one grouping column.")
    sample_col = _pick_column(df.columns, SAMPLE_ID_COLUMNS) or df.columns[0]
    for col in df.columns:                          # tidy text: strip whitespace, blanks become missing
        if not pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].astype("string").str.strip().replace("", pd.NA)
    df[sample_col] = df[sample_col].astype(str).str.strip()
    dupes = df[sample_col][df[sample_col].duplicated()].unique()
    if len(dupes):
        raise ValidationError(f"Sample names must be unique; repeated in metadata: {', '.join(dupes[:5])}")
    return df.set_index(sample_col).rename_axis("sample")


def groupable_columns(meta):
    """Columns usable for a two-group comparison, with the number of samples at each level."""
    result = []
    for col in meta.columns:
        sizes = meta[col].dropna().astype(str).value_counts()
        if 2 <= len(sizes) <= MAX_GROUP_LEVELS:
            result.append({"name": col, "levels": {str(k): int(v) for k, v in sizes.items()}})
    return result


def align_samples(counts, meta):
    """Keep only samples present in both tables. Returns (counts, meta, warnings)."""
    warnings = []
    in_meta = set(meta.index)
    common = [s for s in counts.columns if s in in_meta]
    if not common:
        raise ValidationError(
            "No sample names match between the counts and metadata files. "
            f"Counts columns look like: {', '.join(map(str, counts.columns[:3]))}. "
            f"Metadata samples look like: {', '.join(map(str, list(meta.index)[:3]))}.")
    only_counts = [s for s in counts.columns if s not in in_meta]
    only_meta = [s for s in meta.index if s not in set(counts.columns)]
    if only_counts:
        warnings.append(f"{len(only_counts)} samples in the counts file have no metadata and were "
                        f"dropped (e.g. {', '.join(only_counts[:3])}).")
    if only_meta:
        warnings.append(f"{len(only_meta)} samples in the metadata have no counts and were "
                        f"dropped (e.g. {', '.join(only_meta[:3])}).")
    return counts[common], meta.loc[common], warnings
