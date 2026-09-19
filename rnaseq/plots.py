"""Figures rendered server-side as PNG bytes, so the browser only ever receives an image.

Every figure can be drawn for a light or a dark page. Dark mode is not an inversion of light mode:
each palette has its own colour steps, checked for colour-blind separation and contrast against the
surface it is drawn on (the white card in light mode, #24242b in dark mode).
"""
import contextvars
import functools
import io
import textwrap
import threading
from dataclasses import dataclass

import matplotlib
matplotlib.use("Agg")                      # no GUI: required when running inside a web server
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage


@dataclass(frozen=True)
class Palette:
    """Every colour a figure uses, for one theme.

    Red and blue always mean expression direction: red is higher, blue is lower. Sample groups and
    PCA categories use other hues (orange, aqua, violet) so that "which group" is never confused
    with "which direction".
    """
    name: str
    surface: str            # figure background; matches the card the figure sits on
    strong: str             # titles
    text: str               # axis labels, ticks and secondary text
    label: str              # gene names written on scatter plots
    guide: str              # threshold lines
    neutral: str            # gridlines, zero lines
    edge: str               # network edges
    colors: dict            # "up", "down", "ns"
    groups: tuple           # the two sample groups on the heatmap
    categories: tuple       # PCA categories; three at most, validated all-pairs
    other: str              # categories folded into "Other", size-legend markers
    diverging: object       # heatmap z-scores: two hues with a neutral midpoint
    sequential: dict        # per-direction significance ramps for the dot plot

    @property
    def rc(self):
        """matplotlib settings that give every figure this palette's background and text."""
        return {"figure.facecolor": self.surface, "axes.facecolor": self.surface,
                "savefig.facecolor": self.surface, "text.color": self.strong,
                "axes.labelcolor": self.strong, "axes.edgecolor": self.text,
                "xtick.color": self.text, "ytick.color": self.text,
                "xtick.labelcolor": self.strong, "ytick.labelcolor": self.strong}


def _ramp(name, colours):
    return LinearSegmentedColormap.from_list(name, colours)


LIGHT = Palette(
    name="light", surface="#ffffff", strong="#15151c", text="#52514e", label="#2c3e50",
    guide="#7f8c8d", neutral="#f0efec", edge="#c9c9d2",
    colors={"up": "#e34948", "down": "#2a78d6", "ns": "#c3c2b7"},
    groups=("#eb6834", "#1baf7a"),
    categories=("#eb6834", "#1baf7a", "#4a3aa7"),
    other="#b9b9c4",
    # Light arms darken toward the extremes: more contrast against white means more extreme
    diverging=_ramp("expression-light", ["#104281", "#2a78d6", "#9ec5f4", "#f0efec",
                                         "#f3a3a2", "#e34948", "#8c2322"]),
    sequential={"up": _ramp("up-light", ["#fbd5d4", "#f3a3a2", "#e34948", "#a82a2a"]),
                "down": _ramp("down-light", ["#cde2fb", "#86b6ef", "#2a78d6", "#154a8a"])},
)

DARK = Palette(
    name="dark", surface="#24242b", strong="#ececf1", text="#a9a9b6", label="#d6d6de",
    guide="#6b6b78", neutral="#34343d", edge="#4a4a56",
    colors={"up": "#e66767", "down": "#3987e5", "ns": "#4b4b55"},
    groups=("#d95926", "#199e70"),
    categories=("#d95926", "#199e70", "#9085e9"),
    other="#5d5d68",
    # Dark arms brighten toward the extremes, for the same reason in reverse; the midpoint is a
    # dark neutral so an average value reads as nothing against the dark card
    diverging=_ramp("expression-dark", ["#a6cbf6", "#3987e5", "#1f4a80", "#383840",
                                        "#7a2a2e", "#e66767", "#f7b3b2"]),
    sequential={"up": _ramp("up-dark", ["#4a2226", "#9c3a3c", "#e66767", "#f7b3b2"]),
                "down": _ramp("down-dark", ["#1a2c45", "#25589a", "#3987e5", "#a6cbf6"])},
)

PALETTES = {"light": LIGHT, "dark": DARK}

# The palette for the figure being drawn right now. A ContextVar rather than a plain global, so
# two requests drawing at the same time can never see each other's theme.
_ACTIVE = contextvars.ContextVar("palette", default=LIGHT)

# pyplot and matplotlib's rcParams are process-wide and not thread-safe, and Flask serves requests
# on several threads, so figures are drawn one at a time.
_RENDER_LOCK = threading.Lock()


def _pal():
    return _ACTIVE.get()


def themed(draw):
    """Let a figure function be called with theme="light" or theme="dark"."""
    @functools.wraps(draw)
    def wrapper(*args, theme="light", **kwargs):
        palette = PALETTES.get(theme, LIGHT)
        token = _ACTIVE.set(palette)
        try:
            with _RENDER_LOCK, plt.rc_context(palette.rc):
                return draw(*args, **kwargs)
        finally:
            _ACTIVE.reset(token)
    return wrapper


DIRECTION_LABELS = {"up": "Higher in {a}", "down": "Lower in {a}"}

# Enrichment figures are drawn at the width they are displayed at (roughly the 700-800px card in
# the Pathways panel). Drawing them wider and letting the browser shrink them made the pathway
# names unreadably small.
ENRICH_WIDTH = 7.4                          # inches; at 150 dpi about 1,100px, shown at ~0.7x

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
                    xytext=(4, 3), textcoords="offset points", color=_pal().label)


@themed
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
                   c=_pal().colors[regulation], edgecolors="none",
                   label=f"{regulation} ({len(subset):,})" if regulation != "ns" else f"not significant ({len(subset):,})")

    ax.axhline(-np.log10(padj_cutoff), color=_pal().guide, lw=0.8, ls="--")
    for x in (-lfc_cutoff, lfc_cutoff):
        ax.axvline(x, color=_pal().guide, lw=0.8, ls="--")

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


@themed
def ma_plot(table, group_a, group_b, lfc_cutoff=1.0, label_top=MAX_LABELS):
    """MA plot: mean expression (log scale) against log2 fold change.

    A healthy result is a cloud centred on zero at every expression level. A cloud that bends
    away from zero at low or high expression points to an intensity-dependent bias, such as a
    normalization problem, which the volcano plot cannot show.
    """
    data = table.dropna(subset=["log2FoldChange", "baseMean"])
    data = data[data["baseMean"] > 0].copy()
    data["_x"] = np.log10(data["baseMean"])

    fig, ax = plt.subplots(figsize=(7.5, 6))
    for regulation in ("ns", "down", "up"):
        subset = data[data["regulation"] == regulation]
        ax.scatter(subset["_x"], subset["log2FoldChange"], s=9, alpha=0.6,
                   c=_pal().colors[regulation], edgecolors="none",
                   label=f"{regulation} ({len(subset):,})" if regulation != "ns"
                   else f"not significant ({len(subset):,})")

    ax.axhline(0, color=_pal().text, lw=0.9)
    for y in (-lfc_cutoff, lfc_cutoff):
        ax.axhline(y, color=_pal().guide, lw=0.8, ls="--")

    # Running median of the fold change: should hug zero if there is no intensity bias
    if len(data) > 200:
        ordered = data.sort_values("_x")
        window = max(51, len(ordered) // 40)
        trend = ordered["log2FoldChange"].rolling(window, center=True, min_periods=window // 3).median()
        ax.plot(ordered["_x"], trend, color=_pal().strong, lw=1.4, label="running median")

    # Label the strongest hits, reusing the volcano's collision rule
    labelled = data[data["regulation"] != "ns"].rename(columns={"_x": "_plotx"})
    _label_points_xy(ax, labelled, "_plotx", "log2FoldChange", label_top)

    ax.set_xlabel("mean expression  (log$_{10}$ normalized counts)")
    ax.set_ylabel(f"log$_2$ fold change  ({group_a} vs {group_b})")
    ax.set_title(f"{group_a} vs {group_b}", fontsize=12, pad=26)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=4, fontsize=8,
              frameon=False, markerscale=1.8, borderaxespad=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _to_png(fig)


def _label_points_xy(ax, candidates, x, y, limit):
    """Label points in ranked order, skipping any that would crowd an existing label."""
    if candidates.empty:
        return
    x_range = max(float(candidates[x].max() - candidates[x].min()), 1e-9)
    y_range = max(float(candidates[y].abs().max()) * 2, 1e-9)
    placed = []
    for _, row in candidates.iterrows():
        if len(placed) >= limit:
            break
        px, py = row[x] / x_range, row[y] / y_range
        if any(abs(px - qx) < 0.06 and abs(py - qy) < 0.035 for qx, qy in placed):
            continue
        placed.append((px, py))
        ax.annotate(row["gene"], (row[x], row[y]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points", color=_pal().label)


# Marker shapes give PCA groups a second cue besides colour.
PCA_MARKERS = ["o", "s", "^", "D"]       # shape as well as colour, so identity never rests on colour


@themed
def pca_plot(pca, variance, labels, title, method=""):
    """Samples on the first two principal components, coloured by a metadata column.

    `labels` maps each sample to the category it should be coloured by.
    """
    labels = labels.reindex(pca.index).fillna("missing").astype(str)
    counts = labels.value_counts()
    shown = list(counts.index[:len(_pal().categories)])
    folded = [c for c in counts.index if c not in shown]
    groups = [(name, [name]) for name in shown]
    if folded:
        groups.append((f"Other ({', '.join(folded[:3])}{'…' if len(folded) > 3 else ''})", folded))

    fig, ax = plt.subplots(figsize=(7.5, 6))
    for i, (name, members) in enumerate(groups):
        mask = labels.isin(members)
        colour = _pal().categories[i] if i < len(_pal().categories) else _pal().other
        ax.scatter(pca.loc[mask, "PC1"], pca.loc[mask, "PC2"], s=58, c=colour,
                   marker=PCA_MARKERS[i % len(PCA_MARKERS)],
                   edgecolors=_pal().surface, linewidths=1.2, alpha=0.92, zorder=3,
                   label=f"{name}  (n={int(mask.sum())})")

    ax.axhline(0, color=_pal().neutral, lw=0.9, zorder=0)
    ax.axvline(0, color=_pal().neutral, lw=0.9, zorder=0)
    ax.set_xlabel(f"PC1  ({variance[0] * 100:.1f}% of variance)")
    ax.set_ylabel(f"PC2  ({variance[1] * 100:.1f}% of variance)")
    ax.set_title(title, fontsize=12, pad=26 + 12 * ((len(groups) - 1) // 2))
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2, fontsize=8,
              frameon=False, borderaxespad=0, handletextpad=0.4)
    if method:
        ax.text(1, -0.12, method, transform=ax.transAxes, ha="right", fontsize=7, color=_pal().text)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _to_png(fig)


@themed
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
                  cmap=LinearSegmentedColormap.from_list("groups", _pal().groups, N=2))
    bar_ax.set_axis_off()
    bar_ax.legend(handles=[Patch(facecolor=_pal().groups[0], label=group_a),
                           Patch(facecolor=_pal().groups[1], label=group_b)],
                  loc="lower left", bbox_to_anchor=(0, 1.4), ncol=2, fontsize=8,
                  frameon=False, borderaxespad=0)

    image = ax.imshow(z.to_numpy(), aspect="auto", interpolation="nearest",
                      cmap=_pal().diverging, vmin=-limit, vmax=limit)
    ax.set_yticks(range(len(z)), z.index, fontsize=7, color=_pal().text)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(order)} samples, grouped by condition", fontsize=8, color=_pal().text)
    for spine in ax.spines.values():
        spine.set_visible(False)

    bar = fig.colorbar(image, ax=[bar_ax, ax], fraction=0.035, pad=0.02)
    bar.set_label("expression relative to the gene's mean (z-score)", fontsize=8, color=_pal().text)
    bar.ax.tick_params(labelsize=7, colors=_pal().text)
    bar.outline.set_visible(False)
    return _to_png(fig)


# --------------------------------------------------------------------------------------
# Pathway enrichment figures
#
# These follow the conventions of clusterProfiler's enrichplot, the reference tool for this
# kind of figure in R, so the output is familiar to anyone who has read a paper using it:
#   dot plot      gene ratio on x, dot size = genes in the term, colour = adjusted p-value
#   bar plot      bars ranked by significance
#   network map   terms as nodes, edges where two terms share genes (emapplot)
# --------------------------------------------------------------------------------------

def _wrap(name, width=42, max_lines=2):
    lines = textwrap.wrap(name, width=width)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1][:width - 1] + "…"
    return "\n".join(lines)


def _split_directions(terms, top_n):
    """Most significant terms per direction, keeping only directions that have any.

    GO and Reactome often describe the same process under the same words (for example
    "Extracellular matrix organization" and "extracellular matrix organization"), which would
    otherwise fill the figure with repeats, so only the most significant of each name is kept.
    """
    groups = []
    for direction in ("up", "down"):
        subset = sorted((t for t in terms if t.direction == direction), key=lambda t: t.p_value)
        seen, unique = set(), []
        for term in subset:
            key = term.name.strip().lower()
            if key not in seen:
                seen.add(key)
                unique.append(term)
        if unique:
            groups.append((direction, unique[:top_n]))
    return groups


def _no_data(message):
    fig, ax = plt.subplots(figsize=(7, 2))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=10, color=_pal().text)
    ax.set_axis_off()
    return _to_png(fig)


@themed
def enrichment_dot(terms, group_a, group_b, top_n=10):
    """clusterProfiler-style dot plot, one panel per direction."""
    groups = _split_directions(terms, top_n)
    if not groups:
        return _no_data("No enriched pathways to plot.")

    heights = [len(items) for _, items in groups]
    fig, axes = plt.subplots(len(groups), 1, squeeze=False,
                             figsize=(ENRICH_WIDTH, 1.6 + 0.46 * sum(heights) + 1.1 * len(groups)),
                             gridspec_kw={"height_ratios": heights})
    axes = axes.ravel()

    for ax, (direction, items) in zip(axes, groups):
        items = list(reversed(items))                      # most significant at the top
        y = range(len(items))
        ratios = [t.gene_ratio for t in items]
        sizes = [t.intersection_size for t in items]
        scores = [-np.log10(max(t.p_value, 1e-300)) for t in items]

        # Dot area in points^2, scaled so the smallest term is still visible
        span = max(sizes) or 1
        areas = [40 + 260 * (s / span) for s in sizes]
        dots = ax.scatter(ratios, list(y), s=areas, c=scores, cmap=_pal().sequential[direction],
                          edgecolors=_pal().surface, linewidths=0.8, zorder=3)

        ax.set_yticks(list(y), [_wrap(t.name, width=32) for t in items], fontsize=9, color=_pal().text)
        ax.tick_params(axis="x", labelsize=9, colors=_pal().text)
        ax.grid(axis="x", color=_pal().neutral, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(DIRECTION_LABELS[direction].format(a=group_a), fontsize=11, loc="left", pad=8)
        ax.margins(x=0.16, y=0.12)
        for spine in ax.spines.values():
            spine.set_visible(False)

        bar = fig.colorbar(dots, ax=ax, fraction=0.03, pad=0.015)
        bar.set_label("$-$log$_{10}$ adjusted p", fontsize=8, color=_pal().text)
        bar.ax.tick_params(labelsize=7, colors=_pal().text)
        bar.outline.set_visible(False)

        # Legend for dot size: smallest and largest term in this panel
        for count in sorted({min(sizes), max(sizes)}):
            ax.scatter([], [], s=40 + 260 * (count / span), c=_pal().other,
                       edgecolors=_pal().surface, linewidths=0.8, label=f"{count} genes")
        ax.legend(loc="lower right", fontsize=7, frameon=False, labelspacing=1.1,
                  borderpad=0.6, handletextpad=0.9)

    axes[-1].set_xlabel("gene ratio  (term genes found / genes submitted)", fontsize=9, color=_pal().text)
    fig.tight_layout()
    return _to_png(fig)


@themed
def enrichment_bar(terms, group_a, group_b, top_n=10):
    """Bar plot ranked by significance, one panel per direction."""
    groups = _split_directions(terms, top_n)
    if not groups:
        return _no_data("No enriched pathways to plot.")

    heights = [len(items) for _, items in groups]
    fig, axes = plt.subplots(len(groups), 1, squeeze=False,
                             figsize=(ENRICH_WIDTH, 1.4 + 0.44 * sum(heights) + 1.0 * len(groups)),
                             gridspec_kw={"height_ratios": heights})
    axes = axes.ravel()

    for ax, (direction, items) in zip(axes, groups):
        items = list(reversed(items))
        scores = [-np.log10(max(t.p_value, 1e-300)) for t in items]
        bars = ax.barh(range(len(items)), scores, height=0.62,
                       color=_pal().colors[direction], zorder=3)

        for rect, term in zip(bars, items):
            ax.text(rect.get_width() + max(scores) * 0.015, rect.get_y() + rect.get_height() / 2,
                    f"{term.intersection_size}/{term.term_size}", va="center",
                    fontsize=7, color=_pal().text)

        ax.set_yticks(range(len(items)), [_wrap(t.name, width=32) for t in items], fontsize=9, color=_pal().text)
        ax.tick_params(axis="x", labelsize=9, colors=_pal().text)
        ax.grid(axis="x", color=_pal().neutral, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(DIRECTION_LABELS[direction].format(a=group_a), fontsize=11, loc="left", pad=8)
        ax.margins(x=0.14)
        for spine in ax.spines.values():
            spine.set_visible(False)

    axes[-1].set_xlabel("$-$log$_{10}$ adjusted p-value   (labels show genes found / term size)",
                        fontsize=8, color=_pal().text)
    fig.tight_layout()
    return _to_png(fig)


def _spring_layout(weights, iterations=320, seed=0, gravity=0.0):
    """Force-directed layout (Fruchterman-Reingold).

    Every pair of nodes pushes apart with k^2/d, connected nodes pull together with w*d^2/k,
    and the step size cools over time. Written out rather than pulling in a graph library for
    one figure.
    """
    n = len(weights)
    if n == 1:
        return np.zeros((1, 2))
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, 1, (n, 2))
    k = np.sqrt(1.0 / n)                       # preferred distance between nodes
    temperature = 0.12
    self_pairs = np.eye(n)
    for _ in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]
        distance = np.linalg.norm(delta, axis=-1)
        # Keep the diagonal finite: an infinite self-distance meets a zero self-weight in the
        # attraction term, and 0 * inf is NaN, which would wipe out every position.
        np.fill_diagonal(distance, 1.0)
        distance = np.maximum(distance, 1e-4)
        unit = delta / distance[..., None]                 # zero on the diagonal already
        repulsion = ((1.0 - self_pairs) * k ** 2 / distance)[..., None] * unit
        attraction = (weights * distance ** 2 / k)[..., None] * unit
        step = (repulsion - attraction).sum(axis=1)
        # Gravity pulls every node toward the centre. Without it, pathways that share no genes
        # drift to the far corners and, once the layout is scaled to fit, squeeze every real
        # cluster down to a point.
        step -= gravity * (pos - pos.mean(axis=0))
        length = np.linalg.norm(step, axis=1, keepdims=True)
        pos += step / np.maximum(length, 1e-9) * np.minimum(length, temperature)
        temperature *= 0.99
    return _normalize(pos)


def _normalize(pos):
    span = pos.max(axis=0) - pos.min(axis=0)
    return (pos - pos.min(axis=0)) / np.where(span < 1e-9, 1, span)


# Network labels are laid out in points - the units text is drawn in - so a layout that is free of
# overlaps on paper stays free of them on the page.
LABEL_FONT = 8
CHAR_WIDTH = 0.6          # average glyph width as a fraction of the font size
LINE_HEIGHT = 1.3
LABEL_GAP = 3             # points between a dot and its label
NETWORK_WIDTH = 500       # points; the figure is drawn about this wide so text stays legible


def _node_boxes(labels, areas):
    """Each node's footprint in points: the dot, plus its label hanging underneath.

    Returns half-widths, half-heights, and how far each box's centre sits from its dot's centre
    (negative, because the label is below the dot).
    """
    radius = np.sqrt(areas / np.pi)
    lines = [label.split("\n") for label in labels]
    text_w = np.array([max(len(line) for line in ls) * LABEL_FONT * CHAR_WIDTH for ls in lines])
    text_h = np.array([len(ls) * LABEL_FONT * LINE_HEIGHT for ls in lines])
    above = radius                                   # the dot reaches this far up
    below = radius + LABEL_GAP + text_h              # the dot and its label reach this far down
    half_w = np.maximum(radius, text_w / 2) + 4
    half_h = (above + below) / 2 + 3
    return half_w, half_h, (above - below) / 2


def _separate(pos, half_w, half_h, offset=None, iterations=500):
    """Push apart any two nodes whose footprints overlap.

    A force-directed layout places clusters but stacks their members, so labels collide. Footprints
    are boxes rather than circles because a label is much wider than it is tall. Each overlapping
    pair moves along whichever axis needs the smaller shift, which in practice is usually vertical,
    so the network grows downward instead of getting too wide to read.
    """
    pos = np.asarray(pos, dtype=float).copy()
    n = len(pos)
    offset = np.zeros(n) if offset is None else np.asarray(offset, dtype=float)
    for _ in range(iterations):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                dx = pos[i, 0] - pos[j, 0]
                dy = (pos[i, 1] + offset[i]) - (pos[j, 1] + offset[j])
                overlap_x = half_w[i] + half_w[j] - abs(dx)
                overlap_y = half_h[i] + half_h[j] - abs(dy)
                if overlap_x <= 0 or overlap_y <= 0:
                    continue
                moved = True
                if overlap_x < overlap_y:
                    shift = overlap_x / 2 + 0.5
                    pos[i, 0] += shift if dx >= 0 else -shift
                    pos[j, 0] -= shift if dx >= 0 else -shift
                else:
                    shift = overlap_y / 2 + 0.5
                    pos[i, 1] += shift if dy >= 0 else -shift
                    pos[j, 1] -= shift if dy >= 0 else -shift
        if not moved:
            break
    return pos


@themed
def enrichment_network(terms, group_a, group_b, top_n=12, min_overlap=0.2):
    """Enrichment map (emapplot): pathways that share genes are drawn connected.

    Edge weight is the Jaccard index of the two terms' gene sets, so clusters of edges mark
    groups of pathways describing the same underlying biology.
    """
    chosen = []
    for _, items in _split_directions(terms, top_n):
        chosen.extend(items)
    chosen = [t for t in chosen if t.genes]
    if len(chosen) < 2:
        return _no_data("Not enough pathways with shared genes to draw a network.")

    gene_sets = [set(t.genes) for t in chosen]
    n = len(chosen)
    weights = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            union = gene_sets[i] | gene_sets[j]
            if union:
                weights[i, j] = weights[j, i] = len(gene_sets[i] & gene_sets[j]) / len(union)
    weights[weights < min_overlap] = 0

    counts = np.array([t.intersection_size for t in chosen], dtype=float)
    areas = 110 + 700 * (counts / counts.max())               # dot area in points^2
    labels = [_wrap(t.name, width=20, max_lines=3) for t in chosen]
    half_w, half_h, offset = _node_boxes(labels, areas)

    # Spread the force-directed layout over a canvas in points, then remove every overlap
    # Gravity 3 and 20pt of height per node were chosen by measuring the result on the demo data:
    # together they give a roughly square canvas that the overlap pass barely has to adjust.
    start = _normalize(_spring_layout(weights, gravity=3.0)) * [NETWORK_WIDTH - 110, max(300, 20 * n)]
    pos = _separate(start, half_w, half_h, offset)

    # The canvas is exactly the bounding box of every footprint, so nothing is clipped and one
    # data unit is one point: dot areas and label sizes come out as they were measured.
    pad = 10
    left, right = (pos[:, 0] - half_w).min() - pad, (pos[:, 0] + half_w).max() + pad
    bottom = (pos[:, 1] + offset - half_h).min() - pad
    top = (pos[:, 1] + offset + half_h).max() + pad
    width, height, header = right - left, top - bottom, 46
    fig = plt.figure(figsize=(width / 72, (height + header) / 72))
    ax = fig.add_axes([0, 0, 1, height / (height + header)])
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_axis_off()

    strongest = weights.max() or 1
    for i in range(n):
        for j in range(i + 1, n):
            if weights[i, j]:
                ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
                        color=_pal().edge, linewidth=0.6 + 3.0 * weights[i, j] / strongest,
                        zorder=1, alpha=0.85, solid_capstyle="round")

    handles = []
    for direction in ("up", "down"):
        index = [i for i, t in enumerate(chosen) if t.direction == direction]
        if index:
            handles.append(ax.scatter(pos[index, 0], pos[index, 1], s=areas[index],
                                      c=_pal().colors[direction], edgecolors=_pal().surface, linewidths=1.4,
                                      zorder=2, alpha=0.92,
                                      label=DIRECTION_LABELS[direction].format(a=group_a)))

    for i, label in enumerate(labels):
        # A soft white backing keeps edges that pass behind a label from striking through it
        ax.annotate(label, (pos[i, 0], pos[i, 1]), fontsize=LABEL_FONT, ha="center", va="top",
                    xytext=(0, -(np.sqrt(areas[i] / np.pi) + LABEL_GAP)), textcoords="offset points",
                    color=_pal().text, zorder=3, linespacing=1.1,
                    bbox={"boxstyle": "round,pad=0.12", "fc": _pal().surface, "ec": "none", "alpha": 0.8})

    fig.text(0.012, 1 - 12 / (height + header), "Pathways connected where they share genes",
             fontsize=10, color=_pal().strong, va="top")
    fig.text(0.012, 1 - 28 / (height + header),
             "line width = gene overlap   ·   dot size = genes found", fontsize=8, color=_pal().text, va="top")
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 1 - 6 / (height + header)),
               fontsize=8, frameon=False, markerscale=0.5, ncol=len(handles))
    return _to_png(fig)
