"""Click through every feature of the running app in a real browser and report what works.

The unit tests cover the Python; this covers the page: uploads, every tab and figure, the pathway
views, navigation, the theme toggle, and (optionally) the AI summary. It found the missing favicon
that no unit test could see.

Optional - it needs Playwright, which is not in requirements.txt:
    pip install playwright
    python app.py                          # in another terminal
    python scripts/browser_check.py        # add --ai to include the OpenAI step (about $0.005)

It uses your installed Google Chrome, so no browser download is needed.
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

APP = "http://127.0.0.1:5000/"
REPO = Path(__file__).resolve().parent.parent
RUN_AI = "--ai" in sys.argv
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))


def image_loaded(page, selector):
    return page.eval_on_selector(selector, "el => el.complete && el.naturalWidth > 0")


def wait_image(page, selector, timeout=60000):
    page.wait_for_function(
        "s => { const el = document.querySelector(s); return el && el.complete && el.naturalWidth > 0 }",
        arg=selector, timeout=timeout)


with sync_playwright() as p:
    browser = p.chromium.launch(channel="chrome", headless=True)
    page = browser.new_page(viewport={"width": 1440, "height": 1000})
    errors = []
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.on("console", lambda m: m.type == "error" and errors.append(m.text))
    failed_requests = []
    page.on("response", lambda r: r.status >= 400 and failed_requests.append(f"{r.status} {r.request.method} {r.url}"))

    page.goto(APP)
    page.evaluate("localStorage.clear()")
    page.reload()

    # 1. Initial state -----------------------------------------------------------------------
    locked = page.eval_on_selector_all(".nav-item", "els => els.map(e => e.classList.contains('is-locked'))")
    check("only step 1 is unlocked at start", locked == [False, True, True, True, True], str(locked))
    check("later panels start hidden", page.is_hidden("#panel-setup") and page.is_hidden("#panel-results"))

    # 2. Help button and validation message ---------------------------------------------------
    page.click("#help-btn")
    check("help button shows a message", page.locator("#messages .message").count() == 1)
    page.click("#upload-btn")
    check("upload with no files shows an error",
          "both a counts file" in page.inner_text("#messages"))

    # 3. Upload --------------------------------------------------------------------------------
    page.set_input_files("#counts-file", str(REPO / "demo_data/demo_counts.csv"))
    page.set_input_files("#metadata-file", str(REPO / "demo_data/demo_metadata.csv"))
    page.click("#upload-btn")
    page.wait_for_selector("#panel-setup:not(.is-hidden)", timeout=30000)
    check("upload reveals the comparison step", page.is_visible("#panel-setup"))
    check("samples stat card filled", page.inner_text("#stat-samples") == "60", page.inner_text("#stat-samples"))
    check("genes stat card filled", page.inner_text("#stat-genes") == "22,587", page.inner_text("#stat-genes"))
    check("comparison card updated", page.inner_text("#rail-title") == "Dataset loaded", page.inner_text("#rail-title"))
    check("comparison nav unlocked",
          not page.eval_on_selector('.nav-item[data-target="panel-setup"]', "e => e.classList.contains('is-locked')"))
    check("group selects populated", page.locator("#group-a option").count() >= 2)
    check("upload warning shown", "shared a gene label" in page.inner_text("#messages"))

    # 4. Comparison ----------------------------------------------------------------------------
    page.select_option("#column-select", "tumor_type")
    page.select_option("#group-a", "Met")
    page.select_option("#group-b", "Primary")
    page.click("#analyze-btn")
    page.wait_for_selector("#panel-results:not(.is-hidden)", timeout=120000)
    check("analysis reveals the results step", page.is_visible("#panel-results"))
    check("significant-genes card shows counts",
          page.inner_text("#stat-up") == "1,563" and page.inner_text("#stat-down") == "1,340",
          f"{page.inner_text('#stat-up')} / {page.inner_text('#stat-down')}")
    check("significant-genes placeholder hidden", page.is_hidden("#outcome-empty"))
    check("comparison card shows the groups", page.inner_text("#rail-title") == "Met vs Primary")
    check("results table filled", page.locator("#results-table tbody tr").count() == 50)
    check("pathways and interpretation unlocked", page.is_visible("#panel-pathways") and page.is_visible("#panel-interpret"))

    # 5. The four figures ----------------------------------------------------------------------
    wait_image(page, "#volcano")
    check("volcano loads", image_loaded(page, "#volcano"))
    for name in ["ma", "heatmap", "pca"]:
        page.click(f'.seg[data-figure="{name}"]')
        wait_image(page, f"#{name}")
        others = [n for n in ["volcano", "ma", "heatmap", "pca"] if n != name]
        check(f"{name} tab shows only its figure",
              page.is_visible(f"#{name}") and all(page.is_hidden(f"#{o}") for o in others))
        check(f"download button follows the {name} tab",
              f"/{name}.png" in page.get_attribute("#download-figure", "href"))
    check("PCA colour-by control visible on PCA tab", page.is_visible("#pca-color-field"))
    before = page.get_attribute("#pca", "src")
    page.select_option("#pca-color", "tissue")
    page.wait_for_function("b => document.querySelector('#pca').src !== b", arg=before)
    wait_image(page, "#pca")
    check("PCA recolours by tissue", "color_by=tissue" in page.get_attribute("#pca", "src"))
    page.click('.seg[data-figure="heatmap"]')
    page.select_option("#heatmap-genes", "50")
    page.wait_for_function("() => document.querySelector('#heatmap').src.includes('genes=50')")
    wait_image(page, "#heatmap")
    check("heatmap gene count changes the figure", image_loaded(page, "#heatmap"))
    check("CSV download link set", "/api/results/" in page.get_attribute("#download-csv", "href"))

    # 6. Pathways ------------------------------------------------------------------------------
    page.click("#enrich-btn")
    page.wait_for_selector("#pathways:not(.is-hidden)", timeout=120000)
    wait_image(page, "#pathway-figure")
    check("enrichment shows the dot plot", "kind=dot" in page.get_attribute("#pathway-figure", "src"))
    for view in ["bar", "network"]:
        page.click(f'.seg[data-pathfig="{view}"]')
        page.wait_for_function("v => document.querySelector('#pathway-figure').src.includes('kind=' + v)", arg=view)
        wait_image(page, "#pathway-figure")
        check(f"pathway {view} view loads", image_loaded(page, "#pathway-figure"))
    page.click('.seg[data-pathfig="table"]')
    check("pathway table view shows the table and hides the figure",
          page.is_visible("#pathway-table-wrap") and page.is_hidden("#pathway-figure-wrap"))
    check("pathway table has rows", page.locator("#pathways-table tbody tr").count() > 0)
    results_tab = page.eval_on_selector('.seg[data-figure="heatmap"]', "e => e.classList.contains('is-active')")
    check("pathway tabs don't disturb the results tabs", results_tab)

    # 7. Navigation ----------------------------------------------------------------------------
    page.click('.nav-item[data-target="panel-results"]')
    page.wait_for_timeout(900)
    active = page.eval_on_selector_all(".nav-item.is-active", "els => els.map(e => e.dataset.target)")
    check("rail navigation marks the chosen step", active == ["panel-results"], str(active))

    # 8. Theme ---------------------------------------------------------------------------------
    page.click('.seg[data-figure="volcano"]')
    wait_image(page, "#volcano")
    start = page.evaluate("document.documentElement.dataset.theme || 'system'")
    page.click("#theme-btn")
    theme = page.evaluate("document.documentElement.dataset.theme")
    page.wait_for_function("t => document.querySelector('#volcano').src.includes('theme=' + t)", arg=theme)
    check("theme toggle switches the page", theme in ("light", "dark") and theme != start, f"{start} -> {theme}")
    check("figures redraw in the new theme", f"theme={theme}" in page.get_attribute("#volcano", "src"))
    check("downloads stay light", "download=1" in page.get_attribute("#download-figure", "href"))
    page.reload()
    check("theme choice survives a reload", page.evaluate("document.documentElement.dataset.theme") == theme)

    # 9. Interpretation (costs about $0.005, so opt in) -----------------------------------------
    if RUN_AI:
        page.set_input_files("#counts-file", str(REPO / "demo_data/demo_counts.csv"))
        page.set_input_files("#metadata-file", str(REPO / "demo_data/demo_metadata.csv"))
        page.click("#upload-btn")
        page.wait_for_selector("#panel-setup:not(.is-hidden)")
        page.select_option("#group-a", "Met"); page.select_option("#group-b", "Primary")
        page.click("#analyze-btn")
        page.wait_for_selector("#panel-interpret:not(.is-hidden)", timeout=120000)
        page.fill("#n-genes", "3")
        page.click("#interpret-btn")
        page.wait_for_selector("#summary:not(.is-hidden)", timeout=180000)
        check("AI summary rendered", "What stands out" in page.inner_text("#summary"))
        check("PMID links rendered", page.locator("#summary a[href*='pubmed']").count() > 0)
        check("sources list rendered", page.locator("#evidence-list .evidence-gene").count() > 0)

    check("no failed requests", not failed_requests, "; ".join(failed_requests[:3]))
    check("no JavaScript errors", not errors, "; ".join(errors[:3]))
    browser.close()

failed = [r for r in results if not r[1]]
print(f"\n{len(results) - len(failed)}/{len(results)} checks passed")
sys.exit(1 if failed else 0)
