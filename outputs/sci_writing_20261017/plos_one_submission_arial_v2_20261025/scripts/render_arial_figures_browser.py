"""Browser-only rendering using a task-scoped Fontconfig configuration."""
from pathlib import Path
from collections import Counter
import json
import re
import sys
from playwright.sync_api import sync_playwright


def main():
    root = Path(sys.argv[1]).resolve()
    figures = root / "figures"
    validation = root / "validation"
    manifest = json.loads((validation / "main_figure_transform_audit.json").read_text())
    results = []
    with sync_playwright() as pw:
        browser = pw.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1000, "height": 1100}, device_scale_factor=2)
        for item in manifest["figures"]:
            number = item["number"]
            target = figures / f"Fig{number}_editable.html"
            page.goto(target.as_uri(), wait_until="load")
            page.evaluate("document.fonts.ready")
            dom = page.evaluate('''() => {
              const rect = b => ({x:b.x,y:b.y,width:b.width,height:b.height});
              const svg = document.querySelector('svg');
              return {canvas: svg.getBoundingClientRect().toJSON(),
                texts:[...svg.querySelectorAll('text')].map(t=>({text:t.textContent,
                  bbox:rect(t.getBBox()), font:getComputedStyle(t).fontFamily,
                  size:getComputedStyle(t).fontSize,weight:getComputedStyle(t).fontWeight})),
                workflowNodeOverflows:[...svg.querySelectorAll('[data-node]')].flatMap(n=>{
                  const r=n.querySelector('rect').getBBox();
                  return [...n.querySelectorAll('text')].flatMap(t=>{
                    const b=t.getBBox();return b.x<r.x+8 || b.x+b.width>r.x+r.width-8 || b.y<r.y+5 || b.y+b.height>r.y+r.height-5 ?
                      [{node:n.id,text:t.textContent,box:rect(b),nodebox:rect(r)}]:[];});})};
            }''')
            session = page.context.new_cdp_session(page)
            session.send("DOM.enable")
            session.send("CSS.enable")
            document = session.send("DOM.getDocument")
            ids = session.send("DOM.querySelectorAll", {"nodeId":document["root"]["nodeId"],"selector":"svg text"})["nodeIds"]
            platform = []
            for index, node_id in enumerate(ids):
                fonts = session.send("CSS.getPlatformFontsForNode", {"nodeId":node_id})["fonts"]
                platform.append({"text_index": index, "text": dom["texts"][index]["text"], "fonts":fonts})
            session.detach()
            page.pdf(path=str(figures / f"Fig{number}.pdf"), prefer_css_page_size=True,
                     print_background=True, margin={"top":"0","bottom":"0","left":"0","right":"0"})
            page.locator("svg").screenshot(path=str(validation / f"main_figure_{number}_browser.png"))
            # Hiding text cannot change manually positioned SVG geometry. Identical
            # screenshots establish renderer-level equality of all non-text marks.
            page.evaluate("document.querySelectorAll('svg text').forEach(t=>t.style.visibility='hidden')")
            revised_no_text = validation / f"main_figure_{number}_revised_no_text.png"
            page.locator("svg").screenshot(path=str(revised_no_text))
            baseline = Path(item["source"]).read_text()
            baseline = re.sub(r"<\?xml[^>]*\?>\s*", "", baseline)
            baseline = re.sub(r"<!DOCTYPE[\s\S]*?>\s*", "", baseline)
            page.evaluate("""markup => {
              document.body.innerHTML = markup;
              document.querySelectorAll('svg text').forEach(t=>t.style.visibility='hidden');
            }""", baseline)
            source_no_text = validation / f"main_figure_{number}_source_no_text.png"
            page.locator("svg").screenshot(path=str(source_no_text))
            results.append({"number":number,"dom":dom,"platform_fonts":platform,
                            "no_text_images":[str(source_no_text.relative_to(root)),str(revised_no_text.relative_to(root))]})
        browser.close()
    (validation / "main_figure_browser_font_audit.json").write_text(json.dumps(results,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
