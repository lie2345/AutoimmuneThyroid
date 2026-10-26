"""Render the diagram with the established local Playwright runtime."""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'figures'
with sync_playwright() as pw:
    browser=pw.chromium.launch(headless=True)
    page=browser.new_page(viewport={'width':900,'height':1000},device_scale_factor=2)
    page.goto((FIG/'Figure_1_editable.html').as_uri(),wait_until='load')
    page.evaluate('document.fonts.ready')
    geom=page.evaluate('''() => {
      const nodes=[...document.querySelectorAll('[data-node]')];
      const rect=b=>({x:b.x,y:b.y,width:b.width,height:b.height});
      const violations=[];
      for(const n of nodes){const r=n.querySelector('rect').getBBox();for(const t of n.querySelectorAll('text')){
        const b=t.getBBox();if(b.x<r.x+8||b.x+b.width>r.x+r.width-8||b.y<r.y+5||b.y+b.height>r.y+r.height-5)violations.push({node:n.id,text:t.textContent,box:rect(b),nodebox:rect(r)});}}
      return {violations,texts:[...document.querySelectorAll('svg text')].map(t=>({text:t.textContent,bbox:rect(t.getBBox()),font:getComputedStyle(t).fontFamily})),svg:document.querySelector('svg').getBoundingClientRect().toJSON()};
    }''')
    (ROOT/'validation/figure1_browser_geometry.json').write_text(json.dumps(geom,indent=2))
    assert not geom['violations'],geom['violations']
    page.locator('svg').screenshot(path=str(FIG/'Figure_1_preview.png'))
    page.pdf(path=str(FIG/'Figure_1.pdf'),prefer_css_page_size=True,print_background=True,margin={'top':'0','bottom':'0','left':'0','right':'0'})
    browser.close()
print('Rendered Figure 1 PDF and preview; no text overflow.')
