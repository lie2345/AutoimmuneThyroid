"""Apply a font-family-only Arial change and produce verified PLOS TIFFs.

Uses frozen display SVGs. Never reads statistical data or invokes analysis
builders. Fontconfig and its caches are scoped to this task folder.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import os
import re
import subprocess
import xml.etree.ElementTree as ET

from PIL import Image, ImageChops
import pdfplumber
import pypdfium2 as pdfium
from pypdf import PdfReader
from normalize_pdf_text_sizes import normalize_all

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent / "plos_one_submission_20261025"
FIG = ROOT / "figures"
VALID = ROOT / "validation"
OUT = ROOT / "upload" / "main_figures"
FONTS = ROOT / "assets" / "fonts"
PW = os.environ.get("PLAYWRIGHT_PYTHON", os.sys.executable)
DPI = 600
SVGNS = "{http://www.w3.org/2000/svg}"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(element):
    return (element.tag, tuple(sorted(element.attrib.items())), element.text, element.tail,
            tuple(canonical(child) for child in element))


def strip_font_family(element):
    for node in element.iter():
        if "font-family" in node.attrib:
            node.attrib["font-family"] = "FONT_FAMILY"
        if "style" in node.attrib:
            node.attrib["style"] = re.sub(r"(?<=px )'[^']+'", "'FONT_FAMILY'", node.attrib["style"])
    return element


def css_size(svg_size):
    # Retain the original physical SVG canvas in the HTML print page.
    assert re.fullmatch(r"[0-9.]+(?:mm|pt)", svg_size), svg_size
    return svg_size


def embedded_font_details(reader):
    fonts = []
    for page in reader.pages:
        for key, reference in page["/Resources"]["/Font"].items():
            font = reference.get_object()
            descendant = font.get("/DescendantFonts", [font])[0].get_object()
            descriptor = descendant.get("/FontDescriptor")
            descriptor = descriptor.get_object() if descriptor else {}
            streams = [name for name in ("/FontFile", "/FontFile2", "/FontFile3") if name in descriptor]
            fonts.append({"resource":key,"base_font":str(font.get("/BaseFont")),
                          "descendant_font":str(descendant.get("/BaseFont")),
                          "embedded_streams":streams,
                          "embedded_stream_bytes":sum(len(descriptor[name].get_object().get_data()) for name in streams)})
    return fonts


def main():
    for folder in (FIG, VALID, OUT):
        folder.mkdir(parents=True,exist_ok=True)
    sources = [BASE / "reference_artwork" / "adjusted_export_sources" / "Figure_1_panel_letters_12pt.svg"]
    sources += [BASE / "reference_artwork" / "original_main_figures" / f"Figure_{n}.svg" for n in range(2,5)]
    figures = []
    all_chars = set()
    for number, source in enumerate(sources,1):
        original = source.read_text()
        if number == 1:
            old = 'font-family="Arial, Helvetica, Liberation Sans, sans-serif"'
            assert original.count(old) == 1
            revised = original.replace(old,'font-family="Arial"')
        else:
            assert "'DejaVu Sans'" in original
            revised = original.replace("'DejaVu Sans'", "'Arial'")
        old_tree,new_tree = ET.fromstring(original),ET.fromstring(revised)
        assert canonical(strip_font_family(ET.fromstring(original))) == canonical(strip_font_family(ET.fromstring(revised)))
        old_text = ["".join(t.itertext()) for t in old_tree.iter(SVGNS+"text")]
        assert old_text == ["".join(t.itertext()) for t in new_tree.iter(SVGNS+"text")]
        all_chars.update("".join(old_text))
        target = FIG / f"Fig{number}.svg"
        target.write_text(revised)
        width,height = css_size(new_tree.attrib["width"]),css_size(new_tree.attrib["height"])
        inline = re.sub(r"<\?xml[^>]*\?>\s*", "", revised)
        inline = re.sub(r"<!DOCTYPE[\s\S]*?>\s*", "", inline)
        html = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width,initial-scale=1">'
                f'<style>@page{{size:{width} {height};margin:0}}html,body{{margin:0;'
                f'width:{width};height:{height}}}svg{{display:block}}</style></head><body>' + inline + '</body></html>')
        (FIG / f"Fig{number}_editable.html").write_text(html)
        figures.append({"number":number,"source":str(source),"source_sha256":sha(source),
                        "revised_svg":str(target.relative_to(ROOT)),"revised_sha256":sha(target),
                        "text_nodes":len(old_text),"text_strings":old_text,
                        "font_only_XML_invariance":True,"canvas_attributes":{k:new_tree.attrib[k] for k in ("width","height","viewBox")}})
    transform_audit={"figures":figures,"transformation":"Font family only; sizes, weights, text strings, coordinates, transformations and non-text geometry unchanged.",
                     "source_package":str(BASE)}
    (VALID / "main_figure_transform_audit.json").write_text(json.dumps(transform_audit,ensure_ascii=False,indent=2))
    font_files=sorted(FONTS.rglob("*.ttf"))+sorted(FONTS.rglob("*.TTF"))
    assert font_files,"Arial font resources must be supplied before rendering"
    # fontTools is available in the existing system Python; no installation.
    coverage_script='''from fontTools.ttLib import TTFont
from pathlib import Path
import json,sys
chars=set(json.loads(sys.argv[1]));out=[]
for name in sys.argv[2:]:
 f=TTFont(name);cmap=f.getBestCmap();names=f['name'];out.append({'file':name,'family':names.getDebugName(1),'subfamily':names.getDebugName(2),'postscript':names.getDebugName(6),'missing_chars':[{'char':c,'codepoint':f'U+{ord(c):04X}'} for c in sorted(chars) if ord(c) not in cmap]})
print(json.dumps(out,ensure_ascii=False))
'''
    coverage=json.loads(subprocess.check_output(['/usr/bin/python3','-c',coverage_script,json.dumps(sorted(all_chars)),*map(str,font_files)],text=True))
    assert all(not record['missing_chars'] for record in coverage if record['family']=='Arial')
    assert {record['subfamily'] for record in coverage if record['family']=='Arial'} >= {'Regular','Bold'}
    config=ROOT / 'scripts' / 'figure_fontconfig.xml'
    cache=VALID / 'main_figure_font_cache'
    cache.mkdir(exist_ok=True)
    config.write_text('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd"><fontconfig>'
                      f'<dir>{FONTS}</dir><cachedir>{cache}</cachedir><config><rescan><int>0</int></rescan></config></fontconfig>')
    env=os.environ.copy();env['FONTCONFIG_FILE']=str(config);env['FONTCONFIG_PATH']=str(config.parent)
    font_list=subprocess.check_output(['fc-list',':','family','style','file'],env=env,text=True)
    assert 'Arial' in font_list
    font_audit={'font_files':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in font_files],
                'glyph_coverage':coverage,'task_fontconfig':str(config.relative_to(ROOT)),
                'fontconfig_inventory':font_list,'global_font_configuration_changed':False}
    (VALID / 'main_figure_font_resources.json').write_text(json.dumps(font_audit,ensure_ascii=False,indent=2))
    subprocess.run([PW,str(ROOT/'scripts'/'render_arial_figures_browser.py'),str(ROOT)],env=env,check=True)
    normalize_all(ROOT)
    browser=json.loads((VALID/'main_figure_browser_font_audit.json').read_text())
    exports=[]
    for number,source in enumerate(sources,1):
        record=browser[number-1]
        font_uses=[font for unit in record['platform_fonts'] for font in unit['fonts']]
        assert font_uses and all(font['familyName']=='Arial' for font in font_uses),(number,font_uses)
        assert not record['dom']['workflowNodeOverflows'],record['dom']['workflowNodeOverflows']
        baseline,revised=(ROOT/path for path in record['no_text_images'])
        with Image.open(baseline) as before,Image.open(revised) as after:
            assert before.size==after.size
            difference=ImageChops.difference(before.convert('RGB'),after.convert('RGB'))
            assert difference.getbbox() is None,(number,'non-text geometry changed')
        pdf=FIG/f'Fig{number}.pdf'
        embedded=embedded_font_details(PdfReader(pdf))
        assert embedded and all(f['embedded_streams'] and 'Arial' in f['base_font'] for f in embedded)
        with pdfplumber.open(pdf) as document:
            assert len(document.pages)==1
            page=document.pages[0]
            chars=[c for c in page.chars if c['text'].strip()]
            sizes=Counter(round(c['size'],4) for c in chars)
            families=Counter(c['fontname'] for c in chars)
            assert min(sizes)>=8.0 and max(sizes)<=12.0,(number,sizes)
            assert all('Arial' in name for name in families),(number,families)
            bounds=[min(c['x0'] for c in chars),min(c['top'] for c in chars),max(c['x1'] for c in chars),max(c['bottom'] for c in chars)]
            assert bounds[0]>=0 and bounds[1]>=0 and bounds[2]<=page.width and bounds[3]<=page.height,(number,bounds)
            pdf_mm=[page.width*25.4/72,page.height*25.4/72]
        doc=pdfium.PdfDocument(pdf)
        rendered=doc[0].render(scale=DPI/72).to_pil().convert('RGB')
        tif=OUT/f'Fig{number}.tif'
        rendered.save(tif,format='TIFF',compression='tiff_lzw',dpi=(DPI,DPI))
        rendered.save(FIG/f'Fig{number}.png',dpi=(DPI,DPI))
        preview=rendered.copy();preview.thumbnail((1100,1300),Image.Resampling.LANCZOS)
        preview.save(VALID/f'main_figure_{number}_preview.png')
        with Image.open(tif) as checked:
            checks={'single_frame':checked.n_frames==1,'RGB8':checked.mode=='RGB' and tuple(checked.tag_v2[258])==(8,8,8),
                    'no_alpha':338 not in checked.tag_v2,'LZW':checked.tag_v2[259]==5,
                    'dpi_600':tuple(checked.info['dpi'])==(600.0,600.0),'below_10MB':tif.stat().st_size<10_000_000,
                    'dimensions':66.8<=checked.width/DPI*25.4<=190.5 and checked.height/DPI*25.4<=222.3,
                    'lossless_readback':checked.tobytes()==rendered.tobytes(),'true_Arial_only':True,'embedded_fonts':True,
                    'font_size_8_to_12pt':True,'font_size_8_to_12pt_strict':True,'non_text_pixel_equality':True,'SVG_font_only_invariance':True}
            assert all(checks.values()),(number,checks)
            exports.append({'number':number,'file':str(tif.relative_to(ROOT)),'sha256':sha(tif),'bytes':tif.stat().st_size,
                            'pixels':list(checked.size),'dpi':[float(v) for v in checked.info['dpi']],
                            'physical_mm':[checked.width/DPI*25.4,checked.height/DPI*25.4],'pdf_mm':pdf_mm,
                            'font_sizes_pt':dict(sizes),'embedded_fonts':embedded,'actual_pdf_fonts':dict(families),
                            'text_bounds_pt':bounds,'checks':checks})
        assert sha(source)==figures[number-1]['source_sha256']
        print(f'Fig{number}: Arial only; embedded; {rendered.size}; {tif.stat().st_size} bytes; geometry unchanged')
    audit={'figures':exports,'font_only_source_invariance':True,'all_input_svg_files_unchanged':True,
           'renderer':'Chromium, task-only Fontconfig; PDFium rasterization at 600 dpi',
           'bounded_renderer_difference':'SVG canvas remains exact. Chromium paper-size quantization is retained; nominal 8 pt text initially at 7.995 pt is corrected to effective 8.00005 pt using only Tf operands. See main_figure_pdf_text_normalization.json.',
           'statistical_builders_or_data_used':False,'dependencies_installed':False,'external_upload':False}
    (VALID/'main_figure_export_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
    selfcheck=subprocess.run(['/usr/bin/python3','tools/diagram-design/self_check.py',str(FIG/'Fig1_editable.html')],text=True,capture_output=True)
    (VALID/'main_figure_workflow_self_check.txt').write_text(selfcheck.stdout+selfcheck.stderr)
    assert selfcheck.returncode==0,selfcheck.stdout+selfcheck.stderr


if __name__=='__main__':
    main()
