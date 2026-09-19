"""Tests for the pathway enrichment figures and the layout maths behind the network map."""
import numpy as np
import pytest

from rnaseq.enrichment import Term, _genes_in_term
from rnaseq.plots import (_node_boxes, _separate, _split_directions, _spring_layout, enrichment_bar,
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


def _boxes_overlap(pos, half_w, half_h, offset):
    n = len(pos)
    for i in range(n):
        for j in range(i + 1, n):
            dx = abs(pos[i, 0] - pos[j, 0])
            dy = abs((pos[i, 1] + offset[i]) - (pos[j, 1] + offset[j]))
            if dx < half_w[i] + half_w[j] - 1e-6 and dy < half_h[i] + half_h[j] - 1e-6:
                return True
    return False


class TestSeparate:
    def test_no_two_footprints_overlap_afterwards(self):
        """The guarantee that keeps network labels readable."""
        rng = np.random.default_rng(3)
        pos = rng.uniform(0, 60, (15, 2))                  # deliberately crowded
        half_w, half_h, offset = np.full(15, 45.0), np.full(15, 18.0), np.full(15, -10.0)
        assert _boxes_overlap(pos, half_w, half_h, offset)
        separated = _separate(pos, half_w, half_h, offset)
        assert not _boxes_overlap(separated, half_w, half_h, offset)

    def test_identical_positions_do_not_divide_by_zero(self):
        pos = np.array([[50.0, 50.0], [50.0, 50.0]])
        separated = _separate(pos, np.array([20.0, 20.0]), np.array([10.0, 10.0]))
        assert np.isfinite(separated).all()
        assert not np.allclose(separated[0], separated[1])

    def test_nodes_already_apart_are_left_alone(self):
        pos = np.array([[0.0, 0.0], [500.0, 500.0]])
        assert np.allclose(_separate(pos, np.array([10.0, 10.0]), np.array([10.0, 10.0])), pos)


class TestNodeBoxes:
    def test_long_labels_get_wider_footprints(self):
        half_w, _, _ = _node_boxes(["short", "a much longer pathway name"], np.array([200.0, 200.0]))
        assert half_w[1] > half_w[0]

    def test_extra_lines_make_footprints_taller_and_centre_lower(self):
        _, half_h, offset = _node_boxes(["one line", "two\nlines"], np.array([200.0, 200.0]))
        assert half_h[1] > half_h[0]
        assert offset[1] < offset[0] < 0                    # the label hangs below the dot


class TestGravity:
    def test_isolated_nodes_stay_close_to_the_clusters(self):
        """Without gravity an unconnected node drifts away and squeezes the real clusters."""
        weights = np.zeros((6, 6))
        for i, j in [(0, 1), (1, 2), (0, 2), (3, 4)]:       # two clusters and one loner (5)
            weights[i, j] = weights[j, i] = 0.6

        def loner_distance(gravity):
            pos = _spring_layout(weights, gravity=gravity)
            return np.hypot(*(pos[5] - pos[:5].mean(axis=0)))

        assert loner_distance(3.0) < loner_distance(0.0)


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
