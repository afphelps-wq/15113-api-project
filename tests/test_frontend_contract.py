"""The contract between the page markup and static/app.js.

The JavaScript finds everything by id, class or data attribute. A redesign that renames or drops one
of them breaks a feature silently: the page still loads, but a button stops working. These tests
read both files and check that everything the script relies on is present in the template.
"""
import re
from html.parser import HTMLParser
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
JS = (ROOT / "static" / "app.js").read_text()
TEMPLATE = (ROOT / "templates" / "index.html").read_text()


class _Elements(HTMLParser):
    """Collects every element's tag and attributes."""

    def __init__(self):
        super().__init__()
        self.elements = []

    def handle_starttag(self, tag, attrs):
        self.elements.append((tag, dict(attrs)))


@pytest.fixture(scope="module")
def elements():
    parser = _Elements()
    parser.feed(TEMPLATE)
    return parser.elements


@pytest.fixture(scope="module")
def ids(elements):
    return {attrs["id"] for _, attrs in elements if "id" in attrs}


def with_class(elements, name):
    return [attrs for _, attrs in elements if name in attrs.get("class", "").split()]


def test_every_id_the_script_looks_up_exists(ids):
    wanted = set(re.findall(r'\$\("([^"]+)"\)', JS))
    assert wanted, "the script should look elements up by id"
    assert sorted(wanted - ids) == []


def test_ids_are_unique(elements):
    all_ids = [attrs["id"] for _, attrs in elements if "id" in attrs]
    assert len(all_ids) == len(set(all_ids))


def test_every_figure_tab_has_an_image_and_its_controls(elements, ids):
    """FIGURES in app.js maps each results tab to an <img> of the same id, plus optional controls."""
    figures = re.search(r"const FIGURES = \{(.*?)\n\};", JS, re.S).group(1)
    names = re.findall(r"^  (\w+): \{", figures, re.M)
    controls = re.findall(r'control: "([^"]+)"', figures)
    tabs = [attrs["data-figure"] for attrs in with_class(elements, "seg") if "data-figure" in attrs]
    assert sorted(names) == sorted(tabs) == ["heatmap", "ma", "pca", "volcano"]
    assert set(names) | set(controls) <= ids


def test_pathway_tabs_match_the_views_the_script_handles(elements):
    tabs = [attrs["data-pathfig"] for attrs in with_class(elements, "seg") if "data-pathfig" in attrs]
    assert sorted(tabs) == ["bar", "dot", "network", "table"]


def test_each_nav_item_targets_a_step_panel(elements, ids):
    """The sidebar and the scroll-spy both assume nav targets and .panel sections match one-to-one."""
    targets = [attrs.get("data-target") for attrs in with_class(elements, "nav-item")]
    panels = [attrs.get("id") for attrs in with_class(elements, "panel")]
    assert targets == panels
    assert len(panels) == 5


def test_only_the_first_step_starts_unlocked(elements):
    items = with_class(elements, "nav-item")
    locked = ["is-locked" in attrs.get("class", "") for attrs in items]
    assert locked == [False, True, True, True, True]
    panels = with_class(elements, "panel")
    hidden = ["is-hidden" in attrs.get("class", "") for attrs in panels]
    assert hidden == [False, True, True, True, True]


@pytest.mark.parametrize("table", ["results-table", "pathways-table"])
def test_tables_have_a_body_for_rows(table):
    block = re.search(rf'<table id="{table}">(.*?)</table>', TEMPLATE, re.S)
    assert block and "<tbody>" in block.group(1)


def test_classes_the_script_toggles_are_styled():
    css = (ROOT / "static" / "style.css").read_text()
    for name in set(re.findall(r'classList\.(?:toggle|add|remove)\("([^"]+)"', JS)):
        assert f".{name}" in css, name


def test_theme_is_applied_before_the_stylesheet_paints():
    """The saved theme must be set in <head>, or the page flashes the wrong colours on load."""
    head = TEMPLATE.split("</head>")[0]
    assert 'localStorage.getItem("theme")' in head
