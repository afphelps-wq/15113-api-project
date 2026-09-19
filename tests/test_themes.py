"""Tests for light and dark figure themes."""
import io
import threading

import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import pytest

import app as webapp
from rnaseq import plots
from rnaseq.enrichment import Term

HEX = {"light": plots.LIGHT.surface, "dark": plots.DARK.surface}


def corner_colour(png):
    """The top-left pixel, as a hex string: always background, never data."""
    pixel = mpimg.imread(io.BytesIO(png))[1, 1, :3]
    return "#" + "".join(f"{round(float(c) * 255):02x}" for c in pixel)


def de_table(n=300):
    rng = np.random.default_rng(0)
    lfc = rng.normal(0, 1.5, n)
    padj = np.where(np.abs(lfc) > 2, 1e-6, 0.5)
    return pd.DataFrame({"gene": [f"G{i}" for i in range(n)], "baseMean": rng.lognormal(5, 2, n),
                         "log2FoldChange": lfc, "padj": padj,
                         "regulation": np.where(padj < 0.05, np.where(lfc > 0, "up", "down"), "ns")})


def terms():
    return [Term("GO:BP", f"GO:{i}", f"pathway {i}", 10 ** -(20 - i), 10 + i, 100, 500,
                 "up" if i % 2 else "down", tuple(f"G{j}" for j in range(i, i + 12)))
            for i in range(8)]


def expression():
    samples = [f"S{i}" for i in range(8)]
    rng = np.random.default_rng(1)
    genes = de_table()["gene"][:20]
    frame = pd.DataFrame(rng.lognormal(5, 1, (20, 8)), index=genes, columns=samples)
    return frame, pd.Series(["a"] * 4 + ["b"] * 4, index=samples)


def every_figure(theme):
    table = de_table()
    normalized, conditions = expression()
    pca = pd.DataFrame(np.random.default_rng(2).normal(0, 3, (8, 2)),
                       index=conditions.index, columns=["PC1", "PC2"])
    return {
        "volcano": plots.volcano(table, "a", "b", theme=theme),
        "ma": plots.ma_plot(table, "a", "b", theme=theme),
        "heatmap": plots.heatmap(normalized, table, conditions, "a", "b", theme=theme),
        "pca": plots.pca_plot(pca, (0.4, 0.2), conditions, "a vs b", theme=theme),
        "dot": plots.enrichment_dot(terms(), "a", "b", theme=theme),
        "bar": plots.enrichment_bar(terms(), "a", "b", theme=theme),
        "network": plots.enrichment_network(terms(), "a", "b", theme=theme),
    }


class TestThemedFigures:
    @pytest.mark.parametrize("theme", ["light", "dark"])
    def test_every_figure_paints_its_theme_background(self, theme):
        for name, png in every_figure(theme).items():
            assert corner_colour(png) == HEX[theme], f"{name} in {theme} mode"

    def test_unknown_theme_falls_back_to_light(self):
        assert corner_colour(plots.volcano(de_table(), "a", "b", theme="neon")) == HEX["light"]

    def test_the_active_palette_is_restored_after_drawing(self):
        plots.volcano(de_table(), "a", "b", theme="dark")
        assert plots._pal() is plots.LIGHT

    def test_the_palette_is_restored_even_if_drawing_fails(self):
        with pytest.raises(ValueError):
            plots.heatmap(pd.DataFrame(), de_table().head(1), pd.Series(dtype=str), "a", "b",
                          theme="dark")
        assert plots._pal() is plots.LIGHT

    def test_simultaneous_light_and_dark_requests_do_not_mix(self):
        """Flask serves requests on several threads; each must get the palette it asked for."""
        results, errors = {}, []

        def draw(key, theme):
            try:
                results[key] = corner_colour(plots.volcano(de_table(), "a", "b", theme=theme))
            except Exception as exc:                          # pragma: no cover - surfaced below
                errors.append(exc)

        threads = [threading.Thread(target=draw, args=(i, "dark" if i % 2 else "light"))
                   for i in range(8)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        assert not errors
        for key, colour in results.items():
            assert colour == HEX["dark" if key % 2 else "light"]


class TestPalettes:
    def test_direction_colours_stay_red_and_blue_in_both_themes(self):
        for palette in (plots.LIGHT, plots.DARK):
            up, down = palette.colors["up"], palette.colors["down"]
            assert int(up[1:3], 16) > int(up[5:7], 16)      # more red than blue
            assert int(down[5:7], 16) > int(down[1:3], 16)  # more blue than red

    def test_scatter_categories_never_reuse_the_direction_colours(self):
        for palette in (plots.LIGHT, plots.DARK):
            assert not set(palette.categories) & {palette.colors["up"], palette.colors["down"]}

    def test_at_most_three_scatter_categories(self):
        """Only three categorical colours pass the all-pairs colour-blind check."""
        assert len(plots.LIGHT.categories) == len(plots.DARK.categories) == 3


class TestThemeParameter:
    @pytest.mark.parametrize("query,expected", [
        ("theme=dark", "dark"),
        ("theme=light", "light"),
        ("", "light"),
        ("theme=neon", "light"),                           # unknown values are not trusted
        ("theme=dark&download=1", "light"),                # downloads are always light
    ])
    def test_route_theme_selection(self, query, expected):
        with webapp.app.test_request_context(f"/api/figure/x/volcano.png?{query}"):
            assert webapp._theme() == expected
