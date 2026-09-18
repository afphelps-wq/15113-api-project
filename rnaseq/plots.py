"""Figures rendered server-side as PNG bytes, so the browser only ever receives an image."""
import io

import matplotlib
matplotlib.use("Agg")                      # no GUI: required when running inside a web server
import matplotlib.pyplot as plt
import numpy as np

COLORS = {"up": "#c0392b", "down": "#2471a3", "ns": "#b8bfc7"}
MAX_LABELS = 12


def _to_png(fig):
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buffer.getvalue()


def _label_points(ax, candidates, limit):
    """Label the strongest hits, skipping any that would overlap an existing label.

    Overlap is judged in axes-relative coordinates, so it holds at any data scale.
    """
    if candidates.empty:
        return
    x_range = max(candidates["log2FoldChange"].abs().max() * 2, 1e-9)
    y_range = max(candidates["_y"].max(), 1e-9)
    placed = []
    for _, row in candidates.iterrows():
        if len(placed) >= limit:
            break
        x, y = row["log2FoldChange"] / x_range, row["_y"] / y_range
        if any(abs(x - px) < 0.06 and abs(y - py) < 0.035 for px, py in placed):
            continue
        placed.append((x, y))
        ax.annotate(row["gene"], (row["log2FoldChange"], row["_y"]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points", color="#2c3e50")


def volcano(table, group_a, group_b, padj_cutoff=0.05, lfc_cutoff=1.0, label_top=MAX_LABELS):
    """Volcano plot: log2 fold change against -log10 adjusted p-value.

    `table` must already have a 'regulation' column (see analysis.classify).
    """
    data = table.dropna(subset=["padj", "log2FoldChange"]).copy()
    # padj of exactly 0 (underflow in very large studies) would be infinite on a log scale
    smallest = data.loc[data["padj"] > 0, "padj"].min() if (data["padj"] > 0).any() else 1e-300
    data["_y"] = -np.log10(data["padj"].clip(lower=smallest))

    fig, ax = plt.subplots(figsize=(7.5, 6))
    for regulation in ("ns", "down", "up"):                 # grey first so colours sit on top
        subset = data[data["regulation"] == regulation]
        ax.scatter(subset["log2FoldChange"], subset["_y"], s=9, alpha=0.65,
                   c=COLORS[regulation], edgecolors="none",
                   label=f"{regulation} ({len(subset):,})" if regulation != "ns" else f"not significant ({len(subset):,})")

    ax.axhline(-np.log10(padj_cutoff), color="#7f8c8d", lw=0.8, ls="--")
    for x in (-lfc_cutoff, lfc_cutoff):
        ax.axvline(x, color="#7f8c8d", lw=0.8, ls="--")

    _label_points(ax, data[data["regulation"] != "ns"], label_top)

    ax.set_xlabel(f"log$_2$ fold change  ({group_a} vs {group_b})")
    ax.set_ylabel("$-$log$_{10}$ adjusted p-value")
    ax.set_title(f"{group_a} vs {group_b}", fontsize=12, pad=26)
    # Legend above the axes so it can never sit on top of the data
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=3, fontsize=8,
              frameon=False, markerscale=1.8, borderaxespad=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _to_png(fig)
