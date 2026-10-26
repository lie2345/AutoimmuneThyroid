"""Display-only Figure 2 at final journal width; never fit statistical models."""
from pathlib import Path
import csv, hashlib, json, math, re, shutil, subprocess
from collections import Counter
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath
from PIL import Image

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[2]
OLD = BASE.parent / 'bmc_submission_metadata_v2_20261021'
SRC = BASE.parent / 'revisions/matplotlib_clean_v3_20261018'
for name in ['figures', 'source_data', 'validation', 'evidence', 'reviews', 'rendered']:
    (BASE/name).mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))
def dump(p, data): p.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
protected = [p for p in OLD.rglob('*') if p.is_file()]
protected += [p for p in SRC.rglob('*') if p.is_file()]
protected += [ROOT/'outcome_data/t06'/n for n in ['gwas_gwas_susie_coloc.tsv','susie_signal_catalog.tsv','gwas_gwas_region_architecture.tsv']]
before = {str(p.relative_to(ROOT)):sha(p) for p in protected}
dump(BASE/'evidence/protected_before.json',before)
for n in ['main_figure2_nodes.tsv','main_figure2_edges.tsv']:
    shutil.copy2(SRC/'source_data'/n, BASE/'source_data'/n)
nodes = read(BASE/'source_data/main_figure2_nodes.tsv')
edges = read(BASE/'source_data/main_figure2_edges.tsv')
regions = {r['pair_id']:r for r in read(ROOT/'outcome_data/t06/gwas_gwas_region_architecture.tsv')}
raw_n = {(r['trait_id'],r['region'],r['cs']):r for r in read(ROOT/'outcome_data/t06/susie_signal_catalog.tsv')}
raw_e = {(r['pair_id'],r['graves_cs'],r['aiht_cs']):r for r in read(ROOT/'outcome_data/t06/gwas_gwas_susie_coloc.tsv')}
checks=[]
def check(name, value, detail=None):
    checks.append(dict(name=name,passed=bool(value),detail=detail))
    assert value, (name,detail)
for n in nodes:
    raw=raw_n[(n['trait_id'],n['region'],n['cs'])]
    check('Node source fields '+n['pair_id']+'/'+n['node_code'],all(n[k]==v for k,v in raw.items() if k in n))
for e in edges:
    raw=raw_e[(e['pair_id'],e['graves_cs'],e['aiht_cs'])]
    check('Edge source fields '+e['pair_id']+'/'+e['graves_cs']+'/'+e['aiht_cs'],all(e[k]==v for k,v in raw.items() if k in e))
classes=Counter(e['classification'] for e in edges)
check('27 unique nodes',len(nodes)==len({(n['pair_id'],n['node_key']) for n in nodes})==27)
check('24 unique comparisons',len(edges)==len(raw_e.keys() & {(e['pair_id'],e['graves_cs'],e['aiht_cs']) for e in edges})==24)
STYLES={
 'shared_strong_H4':('#0072B2','o','#0072B2','solid','Strong H4'),
 'distinct_signals_H3':('#D55E00','D','#D55E00','solid','Strong H3'),
 'unresolved_low_purity':('#777777','x','#777777',(0,(4,2)),'Low purity'),
 'unresolved_insufficient_overlap':('#777777','^','white',(0,(1,2)),'Low overlap'),
}
# Exact class identifiers are required; no outcome is recategorized.
check('Frozen class totals',classes==Counter({'shared_strong_H4':6,'distinct_signals_H3':9,'unresolved_low_purity':7,'unresolved_insufficient_overlap':2}),dict(classes))
check('All strong comparisons use local Graves CS1',all(e['graves_cs']=='1' for e in edges if e['classification'] in ['shared_strong_H4','distinct_signals_H3']))
MM=72/25.4; W=170*MM; H=169*MM
INK='#18212B'; MUTED='#505965'
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'path.simplify':False,'savefig.facecolor':'white','axes.unicode_minus':True})
PIDS=['RP01','RP02','RP06','RP07','RP10','RP12']
NAMES=['PTPN22','RTKN2–ARID5B','CTLA4–ICOS','AC007563.2–DIRC3','ZAP70–TMEM131–COX5B','LPP']
ROWS=[80,74,84,74,60,80]
geometry=[]; descriptors=[]
def cubic(points,t):
    weights=((1-t)**3,3*(1-t)**2*t,3*(1-t)*t*t,t**3)
    return tuple(sum(w*p[i] for w,p in zip(weights,points)) for i in (0,1))
def draw(height=H):
    fig=plt.figure(figsize=(170/25.4,height/72),dpi=100)
    ax=fig.add_axes([0,0,1,1],xlim=(0,W),ylim=(height,0)); ax.set_axis_off()
    texts=[]
    def text(x,y,s,size=8,bold=False,color=INK,ha='left'):
        a=ax.text(x,y,s,fontsize=size,fontweight='bold' if bold else 'normal',color=color,ha=ha,va='baseline',zorder=6,clip_on=False)
        texts.append(a); return a
    def symbol(x,y,marker,color,fill,size,z=4):
        return ax.plot([x],[y],linestyle='none',marker=marker,markersize=size,markerfacecolor=fill,markeredgecolor=color,markeredgewidth=.9,zorder=z)[0]
    top=5
    for i,(pid,name,row_h) in enumerate(zip(PIDS,NAMES,ROWS)):
        cx,cy,rad=48,top+row_h/2,27
        ns=sorted([n for n in nodes if n['pair_id']==pid],key=lambda n:(n['side']!='graves',int(n['cs'])))
        gs=[n for n in ns if n['side']=='graves']; aa=[n for n in ns if n['side']=='aiht']
        angles=([180] if len(gs)==1 else [150,210])+{4:[72,24,-24,-72],3:[60,0,-60],2:[40,-40]}[len(aa)]
        pos={n['node_code']:(cx+rad*math.cos(math.radians(a)),cy-rad*math.sin(math.radians(a))) for n,a in zip(ns,angles)}
        ax.add_patch(Circle((cx,cy),rad,fill=False,edgecolor='#DCE1E7',linewidth=.6))
        es=sorted([e for e in edges if e['pair_id']==pid],key=lambda e:(not e['classification'].startswith('unresolved'),int(e['graves_cs']),int(e['aiht_cs'])))
        for e in es:
            g,a='G'+e['graves_cs'],'A'+e['aiht_cs']; p0,p3=pos[g],pos[a]
            pts=[p0,(cx-7,cy+(p0[1]-cy)*.75),(cx+7,cy+(p3[1]-cy)*.75),p3]
            color,marker,fill,dash,_=STYLES[e['classification']]
            patch=PathPatch(MplPath(pts,[MplPath.MOVETO]+[MplPath.CURVE4]*3),fill=False,edgecolor=color,linewidth=1,linestyle=dash,zorder=2)
            patch.set_gid(f'edge_{pid}_{g}_{a}'); ax.add_patch(patch)
            t=.73 if e['graves_cs']=='2' else .51; mx,my=cubic(pts,t)
            symbol(mx,my,marker,color,fill,4.2)
            if height==H: geometry.append(dict(pair_id=pid,start=g,end=a,classification=e['classification'],control_points=pts,marker=[mx,my],line_width_pt=1))
        for n,angle in zip(ns,angles):
            color='#777F87' if n['low_purity']=='TRUE' else '#26313B'
            x,y=pos[n['node_code']]
            a=symbol(x,y,'o' if n['side']=='graves' else 's',color,color if n['side']=='graves' else 'white',7,5)
            a.set_gid('node_'+pid+'_'+n['node_code'])
            theta=math.radians(angle)
            text(cx+37*math.cos(theta),cy-37*math.sin(theta)+2.5,n['node_code'],bold=True,color=color,ha='center')
        text(102,top+10,f'{chr(65+i)}  {pid} · {name}',9,bold=True)
        r=regions[pid]
        text(102,top+22,'Graves '+r['graves_region']+'   AIHT '+r['aiht_region'],color=MUTED)
        for j,n in enumerate(ns):
            label=f'{n["node_code"]}  {n["lead_variant"]}  CS size={n["cs_size"]}  min r²={float(n["cs_min_r2"]):.3f} '+('low' if n['low_purity']=='TRUE' else 'pass')
            text(102,top+34+j*9.6,label,color=MUTED if n['low_purity']=='TRUE' else INK)
            if height==H: descriptors.append(dict(pair_id=pid,node_code=n['node_code'],label=label))
        top+=row_h
    ax.add_line(Line2D([6,W-6],[top+2,top+2],color='#D9DEE4',linewidth=.6))
    for x,style in zip([9,131,253,365],STYLES.values()):
        color,marker,fill,dash,label=style
        ax.add_line(Line2D([x,x+19],[top+13,top+13],color=color,linewidth=1,linestyle=dash))
        symbol(x+9.5,top+13,marker,color,fill,4.2)
        text(x+25,top+16,label)
    return fig,ax,text,texts

fig,ax,text,texts=draw()
fig.canvas.draw(); renderer=fig.canvas.get_renderer()
outside=[t.get_text() for t in texts if not fig.bbox.contains(*t.get_window_extent(renderer).get_points()[0]) or not fig.bbox.contains(*t.get_window_extent(renderer).get_points()[1])]
overlap=[(a.get_text(),b.get_text()) for i,a in enumerate(texts) for b in texts[i+1:] if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer))]
check('All figure text inside canvas',not outside,outside)
check('No text-text overlaps',not overlap,overlap)
check('Final figure minimum 8 pt',min(t.get_fontsize() for t in texts)>=8)
metadata={'Title':'Mixed component correspondence occurs within six non-HLA region pairs','Creator':'Matplotlib; display-only layout of frozen T06 results'}
for ext in ['pdf','svg']:
    fig.savefig(BASE/'figures'/('Figure_2.'+ext),metadata=metadata)
for dpi,suffix in [(600,''),(180,'_preview')]:
    p=BASE/'figures'/f'Figure_2{suffix}.png';fig.savefig(p,dpi=dpi)
    with Image.open(p) as im:
        check(f'{dpi} dpi PNG',all(abs(v-dpi)<.02 for v in im.info['dpi']),dict(pixels=im.size,dpi=im.info['dpi']))
plt.close(fig)
caption_file=BASE.parent/'bmc_submission_20261021/editing/figure_legends_bmc.md'
full='Figure 2. '+caption_file.read_text().split('## Figure 2. ')[1].split('\n\n## Figure 3.')[0].strip()
title,caption=full.split('\n',1)
(BASE/'source_data/Figure_2_caption.txt').write_text(full+'\n')
# Wrap using actual renderer glyph metrics. Caption is unchanged and external to upload artwork.
fig,ax,text,texts=draw(225*MM)
fig.canvas.draw(); renderer=fig.canvas.get_renderer()
from matplotlib.font_manager import FontProperties
font=FontProperties(family='DejaVu Sans',size=8)
words=caption.split(); lines=[]; line=''
for word in words:
    candidate=(line+' '+word).strip()
    width=renderer.get_text_width_height_descent(candidate,font,False)[0]*72/fig.dpi
    if width>W-12 and line: lines.append(line);line=word
    else: line=candidate
lines.append(line)
y=H+11
text(6,y,title,8.5,bold=True)
for s in lines:
    y+=9.5;text(6,y,s,8)
check('Complete legend in proof', ' '.join(lines)== ' '.join(caption.split()))
proof_h=y+5
check('Figure plus full legend fits 225 mm',proof_h<=225*MM,dict(height_mm=proof_h/MM,caption_lines=len(lines)))
fig.set_size_inches(170/25.4,proof_h/72);ax.set_ylim(proof_h,0)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
check('Proof text inside canvas',all(fig.bbox.contains(*t.get_window_extent(renderer).get_points()[0]) and fig.bbox.contains(*t.get_window_extent(renderer).get_points()[1]) for t in texts))
fig.savefig(BASE/'figures/Figure_2_with_caption_proof.pdf',metadata=metadata)
fig.savefig(BASE/'figures/Figure_2_with_caption_proof.png',dpi=180)
plt.close(fig)
svg=(BASE/'figures/Figure_2.svg').read_text();svg=svg[svg.index('<svg'):]
svg=re.sub(r'<title>.*?</title>', '', svg, count=1)
svg=svg.replace('<svg ','<svg role="img" aria-labelledby="f2-title f2-desc" ',1)
pos=svg.index('>')+1
svg=svg[:pos]+'<title id="f2-title">Figure 2 component correspondence</title><desc id="f2-desc">Six region panels with 27 nodes and 24 tested comparisons. All scientific source values preserved.</desc>'+svg[pos:]
(BASE/'figures/Figure_2_editable.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 2 editable vector</title><style>body{margin:0;background:white}svg{display:block;max-width:100%;height:auto}</style>'+svg+'</html>')
root=ET.parse(BASE/'figures/Figure_2.svg').getroot()
ids={e.get('id','') for e in root.iter()}
check('SVG 27 rendered nodes',len([s for s in ids if s.startswith('node_')])==27)
check('SVG 24 rendered edges',len([s for s in ids if s.startswith('edge_')])==24)
check('Editable vector SVG',len(root.findall('.//{http://www.w3.org/2000/svg}text'))==len(descriptors)+27+12+4 and not root.findall('.//{http://www.w3.org/2000/svg}image'))
for p in [BASE/'figures/Figure_2.pdf',BASE/'figures/Figure_2_with_caption_proof.pdf']:
    fonts=subprocess.check_output(['pdffonts',str(p)],text=True)
    check(p.name+' embedded fonts',all(re.search(r'\byes\s+yes\s+yes\b',s) for s in fonts.splitlines()[2:]),fonts)
    check(p.name+' vector only',b'/Subtype /Image' not in p.read_bytes())
check('Protected files unchanged',all(sha(ROOT/p)==h for p,h in before.items()),len(before))
dump(BASE/'evidence/display_geometry.json',geometry)
dump(BASE/'evidence/display_descriptors.json',descriptors)
dump(BASE/'validation/figure2_build.json',dict(status='PASS',figure_mm=[170,169],proof_mm=[170,proof_h/MM],minimum_figure_font_pt=8,caption_font_pt=8,minimum_line_width_pt=.6,source_counts=dict(nodes=27,edges=24,classes=dict(classes)),caption_unchanged=True,protected_files=len(before),checks=checks))
print(json.dumps(dict(figure_mm=[170,169],proof_mm=[170,proof_h/MM],caption_lines=len(lines),checks=len(checks))))
