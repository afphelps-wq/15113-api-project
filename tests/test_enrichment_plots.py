"""Tests for the pathway enrichment figures and the layout maths behind the network map."""
import numpy as np
import pytest

from rnaseq.enrichment import Term, _genes_in_term
from rnaseq.plots import (_separate, _split_directions, _spring_layout, enrichment_bar,
                          enrichment_dot, enrichment_network)

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def term(name, direction="up", p=1e-10, hits=10, size=100, genes=None):
    return Term(source="GO:BP", term_id=f"GO:{abs(hash(name)) % 10**7:07d}", name=name,
                p_value=p, intersection_size=hits, term_size=size, query_size=500,
                direction=direction, genes=tuple(genes or [f"G{i}" for i in range(hits)]))


def a_few_terms():
    return [
        term("biological oxidations", "up", 1e-23, 58, 207, [f"U{i}" for i in range(20)]),
        term("drug ADME", "up", 1e-20, 40, 103, [f"U{i}" for i in range(5, 25)]),
        term("complement cascade", "up", 1e-18, 38, 86, [f"C{i}" for i in range(10)]),
        term("extracellular matrix organization", "down", 1e-13, 55, 296, [f"D{i}" for i in range(20)]),
        term("pancreatic secretion", "down", 1e-8, 27, 99, [f"D{i}" for i in range(15, 35)]),
    ]


class TestGeneRatio:
    def test_ratio_is_hits_over_query_size(self):
        assert term("x", hits=25).gene_ratio == pytest.approx(25 / 500)

    def test_zero_query_size_does_not_divide_by_zero(self):
        empty = Term("GO:BP", "GO:1", "x", 0.01, 0, 10, 0, "up")
        assert empty.gene_ratio == 0.0


class TestIntersectionParsing:
    def test_genes_with_evidence_are_kept_in_query_order(self):
        item = {"intersections": [[], ["REAC"], [], ["GO:BP", "IEA"]]}
        assert _genes_in_term(item, ["A", "B", "C", "D"]) == ("B", "D")

    def test_missing_intersections_is_empty_not_an_error(self):
        assert _genes_in_term({}, ["A", "B"]) == ()

    def test_extra_query_genes_are_ignored(self):
        """zip stops at the shorter sequence if g:Profiler returns fewer entries."""
        assert _genes_in_term({"intersections": [["REAC"]]}, ["A", "B", "C"]) == ("A",)


class TestSplitDirections:
    def test_terms_are_ordered_by_significance(self):
        (_, up), _ = _split_directions(a_few_terms(), top_n=5)
        assert [t.name for t in up] == ["biological oxidations", "drug ADME", "complement cascade"]

    def test_names_differing_only_in_case_are_collapsed(self):
        terms = [term("Extracellular matrix organization", "up", 1e-20),
                 term("extracellular matrix organization", "up", 1e-13),
                 term("something else", "up", 1e-9)]
        (_, kept), = _split_directions(terms, top_n=5)
        assert [t.name for t in kept] == ["Extracellular matrix organization", "something else"]

    def test_a_direction_with_no_terms_is_dropped(self):
        groups = _split_directions([term("only up", "up")], top_n=5)
        assert [direction for direction, _ in groups] == ["up"]


class TestSpringLayout:
    def test_positions_are_finite(self):
        """A zero self-weight against an infinite self-distance once produced NaN everywhere."""
        weights = np.array([[0, 0.5, 0], [0.5, 0, 0], [0, 0, 0]], dtype=float)
        assert np.isfinite(_spring_layout(weights)).all()

    def test_layout_is_deterministic(self):
        weights = np.array([[0, 0.4], [0.4, 0]])
        assert np.allclose(_spring_layout(weights), _spring_layout(weights))

    def test_connected_nodes_end_up_closer_than_unconnected_ones(self):
        weights = np.array([[0, 0.9, 0], [0.9, 0, 0], [0, 0, 0]], dtype=float)
        pos = _spring_layout(weights)
        assert np.hypot(*(pos[0] - pos[1])) < np.hypot(*(pos[0] - pos[2]))

    def test_single_node_does_not_crash(self):
        assert _spring_layout(np.zeros((1, 1))).shape == (1, 2)


class TestSeparate:
    def test_overlapping_nodes_are_pushed_apart(self):
        pos = np.array([[0.5, 0.5], [0.51, 0.5], [0.9, 0.9]])
        spread = _separate(pos, np.array([0.1, 0.1, 0.1]))
        assert np.hypot(*(spread[0] - spread[1])) > np.hypot(*(pos[0] - pos[1]))

    def test_identical_positions_do_not_divide_by_zero(self):
        spread = _separate(np.array([[0.5, 0.5], [0.5, 0.5]]), np.array([0.1, 0.1]))
        assert np.isfinite(spread).all()
        assert not np.allclose(spread[0], spread[1])


class TestFigures:
    @pytest.mark.parametrize("draw", [enrichment_dot, enrichment_bar, enrichment_network])
    def test_each_figure_renders(self, draw):
        assert draw(a_few_terms(), "Met", "Primary").startswith(PNG_MAGIC)

    @pytest.mark.parametrize("draw", [enrichment_dot, enrichment_bar, enrichment_network])
    def test_no_terms_gives_a_placeholder_not_a_crash(self, draw):
        assert draw([], "Met", "Primary").startswith(PNG_MAGIC)

    @pytest.mark.parametrize("draw", [enrichment_dot, enrichment_bar])
    def test_one_direction_only(self, draw):
        only_up = [t for t in a_few_terms() if t.direction == "up"]
        assert draw(only_up, "Met", "Primary").startswith(PNG_MAGIC)

    def test_network_needs_at_least_two_terms_with_genes(self):
        """Falls back to a message rather than failing when there is nothing to connect."""
        assert enrichment_network([term("lonely", genes=[])], "Met", "Primary").startswith(PNG_MAGIC)

    def test_network_handles_terms_with_no_overlap(self):
        apart = [term("one", "up", genes=["A", "B"]), term("two", "up", genes=["C", "D"])]
        assert enrichment_network(apart, "Met", "Primary").startswith(PNG_MAGIC)

    def test_p_value_of_zero_does_not_break_the_log_scale(self):
        terms = [term("underflowed", "up", p=0.0), term("normal", "up", p=1e-5)]
        assert enrichment_bar(terms, "Met", "Primary").startswith(PNG_MAGIC)
        assert enrichment_dot(terms, "Met", "Primary").startswith(PNG_MAGIC)
