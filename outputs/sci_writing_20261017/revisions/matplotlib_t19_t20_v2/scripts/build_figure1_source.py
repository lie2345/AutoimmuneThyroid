"""Wrap the final Matplotlib SVG in accessible, self-contained HTML."""
from pathlib import Path
import re

BASE = Path(__file__).resolve().parents[1]
svg = (BASE / 'figures/figure1.svg').read_text()
svg = svg[svg.index('<svg'):]
svg = re.sub(r'<title>.*?</title>', '', svg, count=1, flags=re.S)
end = svg.index('>')
svg = (svg[:end] + ' role="img" aria-labelledby="figure1-title figure1-desc"' +
       svg[end:end + 1] +
       '<title id="figure1-title">Non-HLA signal architecture and analysis scope</title>' +
       '<desc id="figure1-desc">Graves disease and autoimmune hypothyroidism '
       'summary statistics contribute to 16 overlapping non-HLA region pairs. '
       'The radial display classifies all 42 component-pair comparisons using '
       'both marker shape and color, with an equal sector for each region pair.</desc>' +
       svg[end + 1:])
html = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Figure 1 — editable SVG source</title>
<style>body{margin:0;background:white}svg{display:block;width:180mm;max-width:100%;height:auto}</style>
</head><body>''' + svg + '</body></html>\n'
(BASE / 'figures/figure1_source.html').write_text(html)
