"""Build a source-grounded study workflow; no statistical analysis is run."""
from pathlib import Path
from html import escape
import json

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'figures'
INK = '#202b3a'
MUTED = '#53657a'
ACCENT = '#315e83'
PALE = '#f1f6fa'

CAPTION_TITLE = 'Figure 1. Study workflow for comparing detected association components in Graves’ disease and autoimmune hypothyroidism.'
CAPTION = (
    '(A) FinnGen release 13 (R13) GWAS summary statistics and released SuSiE fine-mapping outputs were processed in separate tracks. '
    'The component audit retained the released definitions, including chromosome X; general summary-statistic quality control was restricted to autosomes. '
    'All component pairs within 16 overlapping non-HLA region pairs were compared, yielding 42 comparisons classified using posterior support, component quality, variant coverage, and credible-set overlap. '
    'H4 denotes a shared-variant hypothesis and H3 a distinct-variant hypothesis for the compared components. '
    'The full catalogue of 312 released components, together with eight strong-H4 pairings, defined 304 organizational analysis units: paired components were combined and all other components were retained separately. '
    'These units were not assumed to be independent loci. Marginal effect directions and equivalence criteria identified ten selected non-HLA representatives. '
    '(B) Complementary analyses used the relevant primary outputs and the indicated molecular, external, or phenotype resources. '
    'BioBank Japan (BBJ) provided selected-variant association context; external component-level correspondence was not established. '
    'The regulatory analysis did not establish a complete chromatin-to-expression-to-disease chain. '
    'Ophthalmopathy queries compared cases with population noncases, rather than estimating modifiers within Graves’ disease. '
    'Unresolved, inconclusive, and untestable evidence remained distinct. '
    'AIHT, autoimmune hypothyroidism; ALT, alternative allele; CS, credible set; caQTL, chromatin-accessibility quantitative trait locus; '
    'eQTL, expression quantitative trait locus; CC-GWAS, case–case GWAS using summary statistics. '
    'Detailed definitions and thresholds are given in Methods and Supplementary Methods.'
)

parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="216mm" viewBox="0 0 720 864" role="img" aria-labelledby="workflow-title workflow-desc">
<title id="workflow-title">Study workflow for Graves’ disease and autoimmune hypothyroidism</title>
<desc id="workflow-desc">Separate released-component and summary-statistic tracks feed component comparisons and marginal association classification. The complete component catalogue also feeds analysis-unit assembly directly. Complementary molecular, external, sensitivity, and phenotype analyses are displayed as separate modules.</desc>
<defs>
<marker id="arrow" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto" markerUnits="userSpaceOnUse"><polygon points="0 0, 8 3, 0 6" fill="{MUTED}"/></marker>
<marker id="arrow-accent" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="{ACCENT}"/></marker>
<marker id="arrow-link" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="{MUTED}"/></marker>
</defs>
<rect x="0" y="0" width="720" height="864" fill="#ffffff"/>
<g font-family="Arial, Helvetica, Liberation Sans, sans-serif" fill="{INK}">
''']
nodes = []
text_items = []

def text(x, y, value, size=12, bold=False, anchor='start', owner=None):
    parts.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{700 if bold else 400}" text-anchor="{anchor}">{escape(value)}</text>')
    text_items.append(dict(x=x,y=y,text=value,size=size,owner=owner))

def node(key, x, y, w, h, title, lines, focal=False):
    nodes.append(dict(id=key,x=x,y=y,w=w,h=h,title=title,lines=lines))
    parts.append(f'<g id="{key}" data-node="true"><rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="{PALE if focal else "#ffffff"}" stroke="{ACCENT if focal else MUTED}" stroke-width="1"/>')
    text(x+16,y+24,title,16,True,owner=key)
    for i,line in enumerate(lines): text(x+16,y+48+i*16,line,owner=key)
    parts.append('</g>')

# Six explicit input relationships. The direct component-catalogue route
# prevents the 304 units from appearing to originate from just 42 comparisons.
paths = [
    ('source','catalogue','M 220 128 V 160'),
    ('source','summary','M 544 128 V 160'),
    ('catalogue','coloc','M 220 256 V 296'),
    ('catalogue','units','M 80 208 H 48 Q 40 208 40 216 V 504 Q 40 512 48 512 H 80'),
    ('coloc','units','M 220 408 V 456'),
    ('summary','units','M 544 256 V 456'),
]
for source,target,path in paths:
    parts.append(f'<path data-from="{source}" data-to="{target}" d="{path}" fill="none" stroke="{MUTED}" stroke-width="1.2" marker-end="url(#arrow)"/>')

text(40,28,'A',20,True)
text(64,28,'Primary comparative framework',16,True)
node('source',40,48,640,80,'FinnGen R13: public aggregate data',[
    'Graves’ disease (E4_GRAVES_STRICT) and AIHT (E4_HYTHY_AI_STRICT)',
    'GWAS summary statistics + released SuSiE outputs · GRCh38',
],True)
node('catalogue',80,160,280,96,'Released component catalogue',[
    '21 Graves + 291 AIHT components',
    'CS audit: purity and Bayes factors',
    'Released definitions include chromosome X',
])
node('summary',408,160,272,96,'Marginal GWAS associations',[
    'Autosomal summary-statistic QC',
    'Exact variant and allele harmonization',
    'ALT-oriented marginal effect estimates',
])
node('coloc',80,296,280,112,'Non-HLA component comparison',[
    '16 overlapping region pairs',
    '42 component-pair colocalization tests',
    'Strong H4 / strong H3 / suggestive H4',
    'Unresolved comparisons retained',
])
node('units',80,456,600,108,'Analysis units and marginal association patterns',[
    '312 released components − 8 strong-H4 pairings = 304 analysis units',
    'Combine paired components; retain all other components separately',
    'Marginal direction + equivalence criteria → 10 selected non-HLA representatives',
],True)

text(40,604,'B',20,True)
text(64,604,'Complementary analyses and evidence context',16,True)
node('regulatory',40,624,312,100,'Regulatory context',[
    'FinnGen immune multiome caQTL / eQTL',
    '26 preselected cell–feature units',
    'GWAS–QTL comparisons; candidate chains',
])
node('bbj',376,624,304,100,'External associations',[
    'BioBank Japan Graves / Hashimoto’s',
    '10 variants × 2 diseases: 20 checks',
    'Direction, source P values and precision',
])
node('sensitivity',40,748,312,100,'Additional summary analyses',[
    'CC-GWAS control-overlap sensitivity',
    'Available genetic-correlation results',
    'Available cell-type annotation results',
])
node('phenotype',376,748,304,100,'HLA and ophthalmopathy',[
    'HLA regional associations: chr6, 25–34 Mb',
    '12 variants × 2 ophthalmopathy endpoints',
    'Population case–noncase comparisons',
])
parts.append('</g></svg>')
svg='\n'.join(parts)
(FIG/'Figure_1.svg').write_text(svg,encoding='utf-8')
html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure 1 Study workflow</title><style>
@page { size:180mm 216mm; margin:0; }
html,body {margin:0;padding:0;background:white;} main{width:180mm;height:216mm;margin:0;} svg{display:block;width:180mm;height:216mm;} @media screen and (max-width:680px){main{width:100%;height:auto;}svg{width:100%;height:auto;}}
</style></head><body><main>'''+svg+'</main></body></html>'
(FIG/'Figure_1_editable.html').write_text(html,encoding='utf-8')
(ROOT/'Figure_1_caption.txt').write_text(CAPTION_TITLE+'\n\n'+CAPTION+'\n',encoding='utf-8')
(ROOT/'source/workflow_spec.json').write_text(json.dumps(dict(
    size_mm=[180,216],viewBox=[0,0,720,864],font='Arial with Liberation Sans fallback',
    caption_title=CAPTION_TITLE,caption=CAPTION,nodes=nodes,edges=[dict(source=a,target=b,path=c) for a,b,c in paths],texts=text_items,
    source_paragraphs_zero_based={'source':[23,24],'catalogue':[27,28,47],'summary':[27,34,35],'coloc':[30,31,48],
    'units':[32,34,35,47,54],'regulatory':[37,58,67],'bbj':[38,59,68],'sensitivity':[25,39,60],'phenotype':[41,42,61,69]},
    limitations=['The 304 units are organizational units, not independent loci.','BBJ lookups do not establish external component correspondence.',
    'No complete regulatory chain was established.','Ophthalmopathy endpoints use population noncases, not within-Graves controls.'],
),ensure_ascii=False,indent=2),encoding='utf-8')
print('Built editable SVG/HTML, caption and evidence specification.')
