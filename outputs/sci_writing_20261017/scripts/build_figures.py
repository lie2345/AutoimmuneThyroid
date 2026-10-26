#!/usr/bin/env python3
"""Render the approved main figures from frozen tabular outputs.

No statistical models or hypothesis tests are run. Confidence intervals in
Figure 3A are display transforms of stored beta and SE; Figure 3B uses stored
90% confidence intervals and existing threshold classifications.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import platform
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET

from PIL import Image
import reportlab
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs/sci_writing_20261017"
FIG = OUT / "figures"
DATA = OUT / "source_data"
REV = OUT / "reviews"
for folder in [FIG, DATA, REV]:
    folder.mkdir(parents=True, exist_ok=True)

INPUTS = [
    "outcome_data/t06/gwas_gwas_region_architecture.tsv",
    "outcome_data/t06/gwas_gwas_susie_coloc.tsv",
    "outcome_data/t06/gwas_gwas_susie_coloc_annotated.tsv",
    "outcome_data/t06/susie_signal_catalog.tsv",
    "outcome_data/t07/genetic_signal_five_classification.tsv",
    "outcome_data/t07/classification_equivalence_sensitivity.tsv",
    "config/traits.tsv",
    "documentation/analysis/PHENOTYPE_AND_CONTROL_PROTOCOL.md",
    "outputs/locus_novelty_20261017/manuscript_blueprint.md",
    "outputs/locus_novelty_20261017/candidate_snapshot.tsv",
]
hash_before = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in INPUTS}

def read_tsv(path):
    with (ROOT / path).open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return [dict(r, source_path=path, source_line=i+2) for i, r in enumerate(rows)]

def write_tsv(name, rows):
    assert rows
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (DATA / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)

CHECKS = []
def check(name, truth, details=None):
    CHECKS.append(dict(check=name, passed=bool(truth), details=details))
    if not truth:
        raise AssertionError(f"{name}: {details}")

arch = read_tsv(INPUTS[0]); pairs = read_tsv(INPUTS[1]); annotated = read_tsv(INPUTS[2])
catalog = read_tsv(INPUTS[3]); units = read_tsv(INPUTS[4]); sens = read_tsv(INPUTS[5])
traits = read_tsv(INPUTS[6]); snapshot = read_tsv(INPUTS[9])
akey = lambda r: (r["pair_id"], r["graves_cs"], r["aiht_cs"])
catkey = lambda r: (r["trait_id"], r["region"], r["cs"])
cidx = {catkey(r): r for r in catalog}
check("catalog keys unique", len(cidx)==len(catalog))
check("pair keys unique", len({akey(r) for r in pairs})==42==len(pairs))
check("annotated matches raw pair identities", {akey(r) for r in pairs}=={akey(r) for r in annotated})
aidx = {akey(r):r for r in annotated}
for r in pairs:
    q = aidx[akey(r)]
    check("pair annotation classification agrees " + ":".join(akey(r)), r["classification"]==q["classification"])
    for field in ["PP.H3", "PP.H4", "sensitivity_PP.H3", "sensitivity_PP.H4"]:
        check("pair annotation posterior agrees "+":".join(akey(r))+":"+field, r[field]==q[field])

ORDER = ["shared_strong_H4", "distinct_signals_H3", "shared_suggestive_H4", "unresolved_low_purity", "unresolved_insufficient_overlap"]
CLASSES = {
    "shared_strong_H4": dict(label="Strong H4", code="H4", color="#0072B2", shape="circle"),
    "distinct_signals_H3": dict(label="Strong H3", code="H3", color="#D55E00", shape="diamond"),
    "shared_suggestive_H4": dict(label="Suggestive H4", code="S", color="#0072B2", shape="open_circle"),
    "unresolved_low_purity": dict(label="Low purity", code="L", color="#777777", shape="cross"),
    "unresolved_insufficient_overlap": dict(label="Insufficient overlap", code="O", color="#777777", shape="triangle"),
}
G_COLOR, A_COLOR = "#0072B2", "#D55E00"
INK, MUTED, RULE = "#18212B", "#505965", "#D9DEE4"
FONT_FILE = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
BOLD_FILE = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")
pdfmetrics.registerFont(TTFont("FigureSans", str(FONT_FILE)))
pdfmetrics.registerFont(TTFont("FigureSansBold", str(BOLD_FILE)))
MM = 72 / 25.4

class Drawing:
    """A small shared vector renderer in physical-point coordinates."""
    def __init__(self, slug, height_mm, title, desc):
        self.slug, self.w, self.h = slug, 180*MM, height_mm*MM
        self.svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{height_mm}mm" viewBox="0 0 {self.w:.5f} {self.h:.5f}" role="img" aria-labelledby="{slug}-title {slug}-desc">', f'<title id="{slug}-title">{html.escape(title)}</title>', f'<desc id="{slug}-desc">{html.escape(desc)}</desc>']
        self.pdf = canvas.Canvas(str(FIG / f"{slug}.pdf"), pagesize=(self.w,self.h), initialFontName="FigureSans")
        self.pdf.setTitle(title); self.pdf.setAuthor("AutoimmuneThyroid research project")
        self.text_bounds = []
        self.rect(0,0,self.w,self.h,fill="#FFFFFF",stroke=None)
    @staticmethod
    def rgb(color):
        return tuple(int(color[i:i+2],16)/255 for i in (1,3,5))
    def rect(self,x,y,w,h,fill="none",stroke=RULE,lw=.5,rx=0):
        self.svg.append(f'<rect x="{x:.4f}" y="{y:.4f}" width="{w:.4f}" height="{h:.4f}" rx="{rx}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{lw}"/>')
        self.pdf.setLineWidth(lw)
        if fill!="none":self.pdf.setFillColorRGB(*self.rgb(fill))
        if stroke:self.pdf.setStrokeColorRGB(*self.rgb(stroke))
        self.pdf.roundRect(x,self.h-y-h,w,h,rx,stroke=int(bool(stroke)),fill=int(fill!="none"))
    def line(self,x1,y1,x2,y2,color=RULE,lw=.5,dash=None):
        d = f' stroke-dasharray="{dash[0]},{dash[1]}"' if dash else ""
        self.svg.append(f'<line x1="{x1:.4f}" y1="{y1:.4f}" x2="{x2:.4f}" y2="{y2:.4f}" stroke="{color}" stroke-width="{lw}"{d}/>')
        self.pdf.setStrokeColorRGB(*self.rgb(color)); self.pdf.setLineWidth(lw)
        self.pdf.setDash(*(dash or ()))
        self.pdf.line(x1,self.h-y1,x2,self.h-y2); self.pdf.setDash()
    def text(self,x,y,text,size=8,color=INK,bold=False,anchor="start"):
        text=str(text); font="FigureSansBold" if bold else "FigureSans"
        width=pdfmetrics.stringWidth(text,font,size)
        left=x if anchor=="start" else x-width/2 if anchor=="middle" else x-width
        self.text_bounds.append((left,y-size,left+width,y+size*.25,text))
        self.svg.append(f'<text x="{x:.4f}" y="{y:.4f}" font-family="DejaVu Sans,Arial,Helvetica,sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}" text-anchor="{anchor}">{html.escape(text)}</text>')
        self.pdf.setFillColorRGB(*self.rgb(color));self.pdf.setFont(font,size)
        self.pdf.drawString(left,self.h-y,text)
    def polygon(self,points,fill,stroke=None,lw=.8):
        pts=" ".join(f"{x:.4f},{y:.4f}" for x,y in points)
        self.svg.append(f'<polygon points="{pts}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{lw}"/>')
        path=self.pdf.beginPath();path.moveTo(points[0][0],self.h-points[0][1])
        for x,y in points[1:]:path.lineTo(x,self.h-y)
        path.close()
        self.pdf.setLineWidth(lw)
        if stroke:self.pdf.setStrokeColorRGB(*self.rgb(stroke))
        if fill!="none":self.pdf.setFillColorRGB(*self.rgb(fill))
        self.pdf.drawPath(path,stroke=int(bool(stroke)),fill=int(fill!="none"))
    def symbol(self,x,y,kind,color,r=3):
        if kind in ("circle","open_circle"):
            fill=color if kind=="circle" else "#FFFFFF"
            self.svg.append(f'<circle cx="{x:.4f}" cy="{y:.4f}" r="{r}" fill="{fill}" stroke="{color}" stroke-width=".9"/>')
            self.pdf.setStrokeColorRGB(*self.rgb(color));self.pdf.setFillColorRGB(*self.rgb(fill));self.pdf.setLineWidth(.9);self.pdf.circle(x,self.h-y,r,stroke=1,fill=1)
        elif kind=="diamond":self.polygon([(x,y-r*1.25),(x+r*1.25,y),(x,y+r*1.25),(x-r*1.25,y)],color)
        elif kind=="square":self.rect(x-r,y-r,2*r,2*r,fill=color,stroke=color,lw=.8)
        elif kind=="triangle":self.polygon([(x,y-r*1.2),(x+r*1.15,y+r),(x-r*1.15,y+r)],"#FFFFFF",color,.9)
        elif kind=="cross":
            self.line(x-r,y-r,x+r,y+r,color,.9);self.line(x-r,y+r,x+r,y-r,color,.9)
    def arrow(self,x,y1,y2):
        self.line(x,y1,x,y2-4,MUTED,.9)
        self.polygon([(x-2.8,y2-4),(x+2.8,y2-4),(x,y2)],MUTED)
    def finish(self,html_source=False):
        self.svg.append("</svg>")
        svg="\n".join(self.svg)
        (FIG/f"{self.slug}.svg").write_text('<?xml version="1.0" encoding="UTF-8"?>\n'+svg,encoding="utf-8")
        self.pdf.showPage();self.pdf.save()
        ET.fromstring(svg)
        bad=[b for b in self.text_bounds if b[0]<-.2 or b[1]<-.2 or b[2]>self.w+.2 or b[3]>self.h+.2]
        check(self.slug+" text within physical canvas",not bad,bad)
        if html_source:
            (FIG/"figure1_source.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure 1 editable vector source</title><style>html,body{margin:0;background:white}svg{display:block;max-width:100%;height:auto}</style><body>'+svg+'</body></html>',encoding="utf-8")
        return dict(figure=self.slug,width_mm=180,height_mm=self.h/MM,minimum_font_pt=min(float(e.attrib["font-size"]) for e in ET.fromstring(svg).iter() if e.tag.endswith("text")),text_elements=len(self.text_bounds))

def panel(d, letter, title, y):
    d.text(16,y,letter,11,bold=True);d.text(32,y,title,9.5,bold=True)

def legend(d,y,omit_suggestive=False):
    d.line(16,y-12,d.w-16,y-12)
    starts=[22,116,210,315,402]
    labels=["Strong H4","Strong H3","Suggestive H4","Low purity","Low overlap"]
    order=ORDER
    if omit_suggestive:
        starts=[22,142,264,386]
        order=[ORDER[0],ORDER[1],ORDER[3],ORDER[4]]
        labels=["Strong H4","Strong H3","Low purity","Low overlap"]
    for x,key,label in zip(starts,order,labels):
        c=CLASSES[key];d.symbol(x,y-2,c["shape"],c["color"],2.6);d.text(x+7,y+1,label,8)

METADATA=[]
counts=Counter(r["classification"] for r in pairs)
check("all five 42-comparison class totals", [counts[c] for c in ORDER]==[8,11,1,17,5],dict(counts))
check("16 region pairs and unique regions",len(arch)==16 and len({r["graves_region"] for r in arch})==15 and len({r["aiht_region"] for r in arch})==16)
traitidx={r["trait_id"]:r for r in traits}
check("endpoint sample sizes",[traitidx[k][v] for k in ["E4_GRAVES_STRICT","E4_HYTHY_AI_STRICT"] for v in ["n_cases","n_controls"]]==["4247","495939","56098","358456"])
f1_rows=[]
for r in arch:
    group=[p for p in pairs if p["pair_id"]==r["pair_id"]]
    cc=Counter(p["classification"] for p in group)
    check("region aggregate matches pairs "+r["pair_id"],len(group)==int(r["n_signal_pair_tests"]) and cc[ORDER[0]]==int(r["n_shared_strong_H4"]) and cc[ORDER[1]]==int(r["n_distinct_signals_H3"]) and cc[ORDER[2]]==int(r["n_shared_suggestive_H4"]) and cc[ORDER[3]]+cc[ORDER[4]]==int(r["n_unresolved"]))
    f1_rows.append(dict(r,**{key:cc[key] for key in ORDER},pair_source_path=INPUTS[1],pair_source_lines=";".join(str(p["source_line"]) for p in group)))
write_tsv("figure1_region_counts.tsv",f1_rows)
write_tsv("figure1_all42_comparisons.tsv",pairs)
write_tsv("figure1_endpoints.tsv",[traitidx[k] for k in ["E4_GRAVES_STRICT","E4_HYTHY_AI_STRICT"]])
d=Drawing("figure1",215,"Scope and outcomes of signal-pair comparisons","FinnGen R13 endpoints and all 42 comparisons across 16 overlapping fine-mapping region pairs; quality failures remain unresolved.")
panel(d,"A","Scope of the existing FinnGen R13 comparison",23)
for x,k,title in [(16,"E4_GRAVES_STRICT","Graves’ disease"),(274,"E4_HYTHY_AI_STRICT","Autoimmune hypothyroidism")]:
    r=traitidx[k];d.rect(x,40,220,72,fill="#F8FAFC",stroke=RULE,lw=.6,rx=3)
    d.text(x+10,57,title,9,bold=True);d.text(x+10,73,k,8)
    d.text(x+10,91,f'{int(r["n_cases"]):,} cases / {int(r["n_controls"]):,} controls',8.5)
    d.arrow(x+110,112,136)
d.rect(16,136,478,61,fill="#F3F7FA",stroke="#A9BBC9",rx=3)
d.text(d.w/2,153,"16 overlapping fine-mapping region pairs",9,bold=True,anchor="middle")
d.text(d.w/2,169,"15 unique Graves regions and 16 unique AIHT regions",8.5,anchor="middle")
d.text(d.w/2,185,"42 pairwise comparisons of detected association components",8.5,anchor="middle")
d.text(16,214,"Finnish-enriched European ancestry · GRCh38 · summary statistics",8,color=MUTED)
panel(d,"B","Outcome counts in every overlapping region pair",239)
cols=[221,280,339,399,459]
d.text(16,259,"Region pair",8,bold=True);d.text(158,259,"Total",8,bold=True,anchor="middle")
for x,key in zip(cols,ORDER):d.symbol(x,255,CLASSES[key]["shape"],CLASSES[key]["color"],2.7)
d.line(16,266,494,266)
for i,r in enumerate(f1_rows):
    y=280+i*16
    if i%2==0:d.rect(16,y-11,478,16,fill="#F7F8FA",stroke=None)
    ch=r["chromosome"];extra="*" if r["pair_id"] in ("RP08","RP09") else ""
    d.text(23,y,r["pair_id"]+extra+" · chr"+ch,8.3)
    d.text(158,y,r["n_signal_pair_tests"],8.3,anchor="middle")
    for x,key in zip(cols,ORDER):d.text(x,y,str(r[key]) if r[key] else "–",8.3,color=CLASSES[key]["color"] if r[key] else "#89929B",anchor="middle",bold=bool(r[key]))
y=280+16*16;d.line(16,y-10,494,y-10);d.text(23,y+3,"All comparisons",8.5,bold=True);d.text(158,y+3,"42",8.5,bold=True,anchor="middle")
for x,key in zip(cols,ORDER):d.text(x,y+3,counts[key],8.5,bold=True,anchor="middle")
legend(d,579)
d.text(16,601,"* RP08 and RP09 share one Graves region. Low overlap = insufficient variant overlap.",8,color=MUTED)
METADATA.append(d.finish(html_source=True))

# Figure 2: exact trait/region/CS catalog joins, including every unresolved edge.
mixed=[r for r in arch if r["architecture"]=="mixed_shared_and_distinct"]
mixedids=[r["pair_id"] for r in mixed]
check("six approved mixed region pairs",mixedids==["RP01","RP02","RP06","RP07","RP10","RP12"])
edges=[r for r in pairs if r["pair_id"] in mixedids]
check("all24 mixed comparisons retained",len(edges)==24 and Counter(r["classification"] for r in edges)==Counter({ORDER[0]:6,ORDER[1]:9,ORDER[3]:7,ORDER[4]:2}))
nodes=[]
for r in edges:
    for side,trait in [("graves","E4_GRAVES_STRICT"),("aiht","E4_HYTHY_AI_STRICT")]:
        key=(trait,r[side+"_region"],r[side+"_cs"])
        check("exact catalog match "+r["pair_id"]+":"+side+":"+r[side+"_cs"],key in cidx)
        q=cidx[key]
        check("pair/catalog CS size and purity agree "+r["pair_id"]+":"+side+":"+r[side+"_cs"],q["cs_size"]==r[side+"_cs_size"] and q["low_purity"]==r[side+"_low_purity"])
        node=dict(q,pair_id=r["pair_id"],side=side,node_key="|".join(key))
        if not any(n["pair_id"]==node["pair_id"] and n["node_key"]==node["node_key"] for n in nodes):nodes.append(node)
write_tsv("figure2_nodes.tsv",nodes);write_tsv("figure2_all24_comparisons.tsv",edges)
check("all15 strong edges reuse local GravesCS1",all(r["graves_cs"]=="1" for r in edges if r["classification"] in ORDER[:2]))
check("RP06 all five unresolved retained",sum(r["pair_id"]=="RP06" and r["classification"].startswith("unresolved") for r in edges)==5)
labels={"RP01":"PTPN22","RP02":"RTKN2–ARID5B","RP06":"CTLA4–ICOS","RP07":"AC007563.2–DIRC3","RP10":"ZAP70–TMEM131–COX5B","RP12":"LPP"}
d=Drawing("figure2",245,"Signal correspondence within six known risk regions","Matrices retain all 24 comparisons in six region pairs. Full-region lead variants, credible-set size and minimum r squared purity are shown for every component. Symbols denote existing pair-level classifications.")
yp=23
gcols={"1":396,"2":458}
for letter,ar in zip("ABCDEF",mixed):
    pid=ar["pair_id"];panel(d,letter,pid+" · "+labels[pid]+" region",yp)
    d.text(16,yp+14,"Regions: Graves "+ar["graves_region"]+"; AIHT "+ar["aiht_region"],8,color=MUTED)
    gn=sorted([n for n in nodes if n["pair_id"]==pid and n["side"]=="graves"],key=lambda n:int(n["cs"]))
    an=sorted([n for n in nodes if n["pair_id"]==pid and n["side"]=="aiht"],key=lambda n:int(n["cs"]))
    yy=yp+28
    for n in gn:
        status="low" if n["low_purity"]=="TRUE" else "pass"
        s=f'Graves G{n["cs"]}: {n["lead_variant"]} · n={n["cs_size"]} · min r²={float(n["cs_min_r2"]):.3f} ({status})'
        d.text(16,yy,s,8,color=MUTED if status=="low" else INK);yy+=11
    d.text(16,yy+1,"AIHT CS / full-region lead variant",8,bold=True)
    d.text(278,yy+1,"n",8,bold=True,anchor="middle");d.text(327,yy+1,"Min r²",8,bold=True,anchor="middle")
    for n in gn:d.text(gcols[n["cs"]],yy+1,"G"+n["cs"],8.5,bold=True,anchor="middle",color=MUTED if n["low_purity"]=="TRUE" else INK)
    yy+=6;d.line(16,yy,494,yy)
    for j,n in enumerate(an):
        yy+=13
        if j%2==0:d.rect(16,yy-9,478,13,fill="#F7F8FA",stroke=None)
        status=n["low_purity"]=="TRUE"
        d.text(20,yy,"A"+n["cs"]+"  "+n["lead_variant"],8,color=MUTED if status else INK)
        check("node label clears CS-size column "+pid+":A"+n["cs"],20+pdfmetrics.stringWidth("A"+n["cs"]+"  "+n["lead_variant"],"FigureSans",8)<263)
        d.text(278,yy,n["cs_size"],8,anchor="middle")
        d.text(327,yy,f'{float(n["cs_min_r2"]):.3f}'+(" L" if status else ""),8,color=MUTED if status else INK,anchor="middle")
        for g in gn:
            matches=[r for r in edges if r["pair_id"]==pid and r["graves_cs"]==g["cs"] and r["aiht_cs"]==n["cs"]]
            check("one plotted matrix cell "+pid+":"+g["cs"]+":"+n["cs"],len(matches)==1)
            c=CLASSES[matches[0]["classification"]];d.symbol(gcols[g["cs"]],yy-2,c["shape"],c["color"],3.2)
    yp=yy+14
legend(d,665,omit_suggestive=True)
d.text(16,685,"n = credible-set size; L = low purity. Gene names locate regions; coordinates are GRCh38.",8,color=MUTED)
METADATA.append(d.finish())

# Figure 3: no candidate omission; same stored ALT allele for both marginal estimates.
priority=[r for r in units if r["five_class_primary"]!="unresolved"]
priority.sort(key=lambda r:int(r["signal_unit_id"][3:]))
check("ten prioritized candidates agree with frozen snapshot",len(priority)==10 and {r["signal_unit_id"] for r in priority}=={r["signal_unit_id"] for r in snapshot})
sidx={r["signal_unit_id"]:r for r in snapshot}
for r in priority:
    q=sidx[r["signal_unit_id"]]
    check("candidate identity/estimates agree "+r["signal_unit_id"],all(r[f]==q[f] for f in ["representative_variant","representative_reported_rsids","graves_beta","graves_se","aiht_beta","aiht_se"]))
    allele=r["representative_variant"].split("_")
    check("candidate effect allele matches ALT "+r["signal_unit_id"],r["ref"]==allele[-2] and r["alt"]==allele[-1])
f3=[]
for r in priority:
    for trait in ["graves","aiht"]:
        b,se=float(r[trait+"_beta"]),float(r[trait+"_se"])
        f3.append(dict(signal_unit_id=r["signal_unit_id"],classification=r["five_class_primary"],representative_variant=r["representative_variant"],rsid=r["representative_reported_rsids"],ref=r["ref"],effect_allele_ALT=r["alt"],trait=trait,beta=b,se=se,ci95_lower=b-1.96*se,ci95_upper=b+1.96*se,display_formula="beta +/- 1.96*SE",source_path=r["source_path"],source_line=r["source_line"],snapshot_source_path=sidx[r["signal_unit_id"]]["source_path"],snapshot_source_line=sidx[r["signal_unit_id"]]["source_line"]))
write_tsv("figure3_marginal_estimates.tsv",f3)
f3b=[]
for r in priority:
    if r["signal_unit_id"] not in ["GSU138","GSU231"]:continue
    for bound in [1.05,1.10,1.15]:
        ss=[s for s in sens if s["signal_unit_id"]==r["signal_unit_id"] and float(s["equivalence_or_upper"])==bound]
        check("one existing equivalence class "+r["signal_unit_id"]+":"+str(bound),len(ss)==1)
        lo,hi=float(r["other_ci90_lower"]),float(r["other_ci90_upper"])
        flag="other_equivalent_or1_"+{1.05:"05",1.1:"10",1.15:"15"}[bound]
        inside=lo > -math.log(bound) and hi < math.log(bound)
        check("plotted bounds reproduce stored equivalence flag "+r["signal_unit_id"]+":"+str(bound),inside==(r[flag]=="TRUE"))
        f3b.append(dict(signal_unit_id=r["signal_unit_id"],other_trait=r["other_trait"],beta=r["other_beta"],stored_ci90_lower=lo,stored_ci90_upper=hi,equivalence_or_upper=bound,log_bound=math.log(bound),stored_equivalent=r[flag],existing_five_class=ss[0]["five_class"],source_path=r["source_path"],source_line=r["source_line"],sensitivity_source_path=ss[0]["source_path"],sensitivity_source_line=ss[0]["source_line"]))
write_tsv("figure3_equivalence_sensitivity.tsv",f3b)
d=Drawing("figure3",205,"Marginal allelic associations and equivalence-bound sensitivity","Both traits are aligned to the stored ALT allele at all ten prioritized discovery variants. Figure A shows marginal beta and display-derived 95% confidence intervals. Figure B uses stored AIHT 90% confidence intervals and existing equivalence-bound classifications for NRG1 and TSHR.")
panel(d,"A","Marginal allelic associations at all ten prioritized variants",23)
for x,c,shape,label in [(252,G_COLOR,"circle","Graves"),(353,A_COLOR,"square","AIHT")]:d.symbol(x,40,shape,c,2.8);d.text(x+8,43,label,8)
left,right=232,490;lo,hi=-.75,.60
xf=lambda b:left+(b-lo)/(hi-lo)*(right-left)
top,bottom=61,338
for tick in [-.6,-.3,0,.3,.6]:
    x=xf(tick);d.line(x,top-8,x,bottom+6,"#ADB6BF" if tick==0 else "#E6E9ED",.7 if tick==0 else .4,dash=(2,2) if tick==0 else None)
    d.text(x,bottom+20,f"{tick:.1f}".replace("-","−"),8,anchor="middle")
for i,r in enumerate(priority):
    y=top+i*29
    if i%2==0:d.rect(16,y-10,208,28,fill="#F7F8FA",stroke=None)
    gene={"GSU052":"CTLA4–ICOS","GSU042":"TMEM131–ZAP70","GSU056":"AC007563.2/DIRC3","GSU254":"LINC01229/MAFTRR"}.get(r["signal_unit_id"],r["representative_nearest_genes"])
    check("candidate label clears forest axis "+r["signal_unit_id"],19+pdfmetrics.stringWidth(r["signal_unit_id"]+" · "+r["representative_reported_rsids"]+" · "+gene,"FigureSans",8)<left-4)
    d.text(19,y-1,r["signal_unit_id"]+" · "+r["representative_reported_rsids"]+" · "+gene,8)
    d.text(19,y+11,r["representative_variant"]+" · ALT="+r["alt"],8,color=MUTED)
    for trait,offset,color,shape in [("graves",-3,G_COLOR,"circle"),("aiht",5,A_COLOR,"square")]:
        row=next(q for q in f3 if q["signal_unit_id"]==r["signal_unit_id"] and q["trait"]==trait)
        check("95CI inside display "+r["signal_unit_id"]+":"+trait,lo<row["ci95_lower"]<=row["beta"]<=row["ci95_upper"]<hi)
        x1,x2=xf(row["ci95_lower"]),xf(row["ci95_upper"]);yy=y+offset
        radius=2.2 if trait=="graves" else 1.2
        check("CI caps extend beyond marker "+r["signal_unit_id"]+":"+trait,x2-x1>2*radius)
        d.line(x1,yy,x2,yy,color,1);d.line(x1,yy-2,x1,yy+2,color,.7);d.line(x2,yy-2,x2,yy+2,color,.7);d.symbol(xf(row["beta"]),yy,shape,color,radius)
d.text((left+right)/2,376,"Marginal beta (log odds ratio) per ALT allele; 95% CI",8,anchor="middle")
panel(d,"B","AIHT equivalence at NRG1 and TSHR: stored 90% CIs",406)
bl,br=104,289;bx=lambda b:bl+(b+.16)/.32*(br-bl)
for bound,fill in [(1.15,"#F0F2F5"),(1.10,"#E1E6EB"),(1.05,"#CAD3DC")]:
    d.rect(bx(-math.log(bound)),430,bx(math.log(bound))-bx(-math.log(bound)),72,fill=fill,stroke=None)
    for v in [-math.log(bound),math.log(bound)]:d.line(bx(v),430,bx(v),502,"#9DA9B5",.5,dash=(2,2))
d.line(bx(0),430,bx(0),502,MUTED,.7)
d.text(402,426,"Existing class at OR upper bound",8,bold=True,anchor="middle")
for x,txt in zip([330,399,468],["1.05","1.10*","1.15"]):d.text(x,443,txt,8,bold=True,anchor="middle")
for y,sid,gene in [(461,"GSU138","NRG1"),(492,"GSU231","TSHR")]:
    r=next(q for q in f3b if q["signal_unit_id"]==sid)
    d.text(16,y-2,gene,8.5,bold=True);d.text(16,y+9,sid,8,color=MUTED)
    x1,x2=bx(r["stored_ci90_lower"]),bx(r["stored_ci90_upper"])
    d.line(x1,y-1,x2,y-1,A_COLOR,1);d.line(x1,y-3,x1,y+1,A_COLOR,.7);d.line(x2,y-3,x2,y+1,A_COLOR,.7);d.symbol(bx(float(r["beta"])),y-1,"square",A_COLOR,1.5)
    for x,bound in zip([330,399,468],[1.05,1.10,1.15]):
        cls=next(q for q in f3b if q["signal_unit_id"]==sid and q["equivalence_or_upper"]==bound)["existing_five_class"]
        d.text(x,y+2,"Unresolved" if cls=="unresolved" else "Graves-pred.",8,anchor="middle",color=MUTED if cls=="unresolved" else INK)
for tick in [-.15,0,.15]:d.text(bx(tick),519,f"{tick:.2f}".replace("-","−"),8,anchor="middle")
d.text((bl+br)/2,536,"AIHT marginal beta (log odds ratio)",8,anchor="middle")
d.text(16,558,"Shaded bounds: ±ln(1.05), ±ln(1.10), ±ln(1.15). * Primary bound. Graves-pred. = Graves-predominant.",8,color=MUTED)
d.text(16,574,"Class labels are from the existing T07 analysis; the plots display marginal associations.",8,color=MUTED)
METADATA.append(d.finish())

# Render fixed-size PDF exports with the existing Poppler CLI, retain PNG pHYs.
for meta in METADATA:
    slug=meta["figure"]
    for dpi,suffix in [(600,""),(180,"_preview")]:
        dest=FIG/(slug+suffix)
        result=subprocess.run(["pdftoppm","-singlefile","-png","-r",str(dpi),str(FIG/(slug+".pdf")),str(dest)],capture_output=True,text=True)
        check(slug+f" PNG render {dpi}dpi",result.returncode==0,result.stderr)
        path=dest.with_suffix(".png")
        with Image.open(path) as im:
            im.load();px=im.size
            im.save(path,dpi=(dpi,dpi))
        with Image.open(path) as im:
            actual=im.info.get("dpi")
        expected=(math.ceil(180/25.4*dpi),math.ceil(meta["height_mm"]/25.4*dpi))
        check(slug+f" PNG dimensions {dpi}dpi",px==expected,dict(actual=px,expected=expected))
        check(slug+f" PNG resolution metadata {dpi}dpi",actual is not None and all(abs(v-dpi)<.02 for v in actual),actual)
        meta["png"+suffix]=dict(path=str(path.relative_to(ROOT)),dimensions=px,dpi=actual)

hash_after={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
check("all authoritative inputs unchanged",hash_before==hash_after)
record=dict(status="PASS_DISPLAY_AND_SOURCE_CHECKS_PENDING_INDEPENDENT_FIGURE_REVIEW",generated_utc=datetime.now(timezone.utc).isoformat(),scope="Three figures only; no models or hypothesis tests rerun.",command="uv run --no-project --no-sync --python python3 python outputs/sci_writing_20261017/scripts/build_figures.py",python=sys.version,executable=sys.executable,reportlab_version=reportlab.Version,pillow_version=Image.__version__,platform=platform.platform(),font_files=[str(FONT_FILE),str(BOLD_FILE)],font_family="DejaVu Sans",inputs_sha256=hash_after,figures=METADATA,color_symbol_metadata=CLASSES,checks=CHECKS)
(REV/"figure_build_validation.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
(REV/"figure_build_notes.md").write_text("# Main figure build record\n\nFigures were rendered from existing output tables. No GWAS, SuSiE, coloc, classification or equivalence tests were rerun. Figure 3A 95% CIs use beta ± 1.96 SE solely for display. Figure 3B uses stored 90% CIs and stored classifications.\n\nAll figures have a 180 mm final width; Figure 1 is 215 mm high, Figure 2 245 mm high, and Figure 3 205 mm high. Editable text is retained in SVG; PDF embeds DejaVu Sans. High-resolution PNGs are 600 dpi, previews 180 dpi. The SVG uses installed DejaVu Sans, with Arial/Helvetica fallbacks.\n\nSource extracts contain exact original file and physical line numbers including the header. Figure 2 nodes are joined on trait + full region + CS, not on the truncated common-set hit. Every source input is hashed before and after production.\n\nScientific interpretation remains at pair-level association correspondence; H3 is not disease specificity, marginal allelic direction is not conditional causal direction, and genes are positional context. AI-assisted source/display checks do not replace independent scientific or human review.\n",encoding="utf-8")
print(json.dumps(dict(figures=METADATA,checks=len(CHECKS),passed=sum(c["passed"] for c in CHECKS)),indent=2))
