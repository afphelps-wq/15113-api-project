"""Figures rendered server-side as PNG bytes, so the browser only ever receives an image."""
import io

import matplotlib
matplotlib.use("Agg")                      # no GUI: required when running inside a web server
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage

# One colour system across both figures. Red and blue always mean expression direction: red is
# higher, blue is lower. The group annotation on the heatmap deliberately uses two different hues
# so that "which group" is never confused with "which direction".
COLORS = {"up": "#e34948", "down": "#2a78d6", "ns": "#c3c2b7"}
GROUP_COLORS = ["#eb6834", "#1baf7a"]      # orange, aqua
NEUTRAL = "#f0efec"
TEXT = "#52514e"

# Diverging ramp for z-scores: two opposite hues with a neutral grey midpoint, equal steps per arm.
EXPRESSION_CMAP = LinearSegmentedColormap.from_list("expression", [
    "#104281", "#2a78d6", "#9ec5f4", NEUTRAL, "#f3a3a2", "#e34948", "#8c2322"])

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


def heatmap(normalized, table, conditions, group_a, group_b, top_n=30):
    """Clustered heatmap of the most significant genes.

    Each gene is z-scored across samples, so the colour shows whether a sample is high or low
    *for that gene* rather than how abundant the gene is. Genes are clustered so similar
    patterns sit together; samples stay grouped by condition, which is what makes a clean
    split visible (or reveals that there isn't one).
    """
    ranked = [g for g in table["gene"] if g in normalized.index][:top_n]
    if len(ranked) < 2:
        raise ValueError("Not enough significant genes to draw a heatmap.")

    # Samples ordered by group, so the two conditions form contiguous blocks
    order = [s for s in conditions.index if conditions[s] == group_a] + \
            [s for s in conditions.index if conditions[s] == group_b]
    order = [s for s in order if s in normalized.columns]
    data = np.log2(normalized.loc[ranked, order] + 1)

    # A gene with the same value in every sample has nothing to show, and its correlation with
    # anything else is undefined, which would break the clustering below.
    spread = data.std(axis=1)
    data = data[spread > 0]
    if len(data) < 2:
        raise ValueError("Not enough significant genes vary between samples to draw a heatmap.")
    z = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

    # Cluster genes on their z-scored profiles; correlation distance groups by shape, not level
    if len(z) > 2:
        z = z.iloc[leaves_list(linkage(z.to_numpy(), method="average", metric="correlation"))]

    limit = float(np.nanpercentile(np.abs(z.to_numpy()), 98)) or 1.0   # symmetric, outlier-robust
    height = max(4.0, 0.22 * len(z) + 2.2)
    fig, (bar_ax, ax) = plt.subplots(
        2, 1, figsize=(max(7.0, 0.11 * len(order) + 3.4), height),
        gridspec_kw={"height_ratios": [1, max(12, len(z))], "hspace": 0.02})

    counts = {group_a: order.count(group_a), group_b: 0}
    group_index = [0 if conditions[s] == group_a else 1 for s in order]
    bar_ax.imshow([group_index], aspect="auto", interpolation="nearest",
                  cmap=LinearSegmentedColormap.from_list("groups", GROUP_COLORS, N=2))
    bar_ax.set_axis_off()
    bar_ax.legend(handles=[Patch(facecolor=GROUP_COLORS[0], label=group_a),
                           Patch(facecolor=GROUP_COLORS[1], label=group_b)],
                  loc="lower left", bbox_to_anchor=(0, 1.4), ncol=2, fontsize=8,
                  frameon=False, borderaxespad=0)

    image = ax.imshow(z.to_numpy(), aspect="auto", interpolation="nearest",
                      cmap=EXPRESSION_CMAP, vmin=-limit, vmax=limit)
    ax.set_yticks(range(len(z)), z.index, fontsize=7, color=TEXT)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(order)} samples, grouped by condition", fontsize=8, color=TEXT)
    for spine in ax.spines.values():
        spine.set_visible(False)

    bar = fig.colorbar(image, ax=[bar_ax, ax], fraction=0.035, pad=0.02)
    bar.set_label("expression relative to the gene's mean (z-score)", fontsize=8, color=TEXT)
    bar.ax.tick_params(labelsize=7, colors=TEXT)
    bar.outline.set_visible(False)
    return _to_png(fig)
