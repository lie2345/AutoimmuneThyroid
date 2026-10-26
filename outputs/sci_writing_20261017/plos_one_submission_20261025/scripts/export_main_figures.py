"""Export the frozen main figures to PLOS ONE-compatible raster containers.

Run using the existing Codex bundled Python through uv --no-sync. No analytical
builder or source data is used. The only display edit reduces Figure 1 panel
letters from 14.17 pt to 11.9 pt in a separate SVG copy. Actual installed font
families are recorded rather than misreported as Arial.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import pdfplumber
import pypdfium2 as pdfium
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "plos_one_study_flow_20261025" / "figures"
ARCHIVE = ROOT / "reference_artwork" / "original_main_figures"
EDITABLE = ROOT / "reference_artwork" / "adjusted_export_sources"
OUT = ROOT / "upload" / "main_figures"
VALID = ROOT / "validation"
PLAYWRIGHT_PYTHON = __import__("os").environ.get("PLAYWRIGHT_PYTHON", sys.executable)
DPI = 600


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def font_match(family):
    return subprocess.check_output(["fc-match", family], text=True).strip()


def main():
    for folder in (ARCHIVE, EDITABLE, OUT, VALID):
        folder.mkdir(parents=True, exist_ok=True)
    preserved = []
    for number in range(1, 5):
        for ext in ("svg", "pdf"):
            source = SOURCE / f"Figure_{number}.{ext}"
            target = ARCHIVE / source.name
            shutil.copy2(source, target)
            preserved.append({"source": str(source), "copy": str(target),
                              "sha256": sha(source), "copy_identical": sha(source) == sha(target)})
    original_svg = (SOURCE / "Figure_1.svg").read_text()
    edited_svg = original_svg.replace('font-size="20"', 'font-size="16.8"')
    assert original_svg.count('font-size="20"') == 2
    assert edited_svg.replace('font-size="16.8"', 'font-size="20"') == original_svg
    adjusted = EDITABLE / "Figure_1_panel_letters_12pt.svg"
    adjusted.write_text(edited_svg)
    markup = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
              '<style>@page{size:180mm 216mm;margin:0}html,body{margin:0;'
              'width:180mm;height:216mm}svg{display:block}</style></head><body>'
              + edited_svg + '</body></html>')
    html = EDITABLE / "Figure_1_panel_letters_12pt.html"
    html.write_text(markup)
    pdf1 = EDITABLE / "Figure_1_panel_letters_12pt.pdf"
    browser_code = '''from pathlib import Path
import sys
from playwright.sync_api import sync_playwright
with sync_playwright() as pw:
    browser = pw.chromium.launch(headless=True)
    page = browser.new_page()
    page.goto(Path(sys.argv[1]).as_uri(), wait_until="load")
    page.evaluate("document.fonts.ready")
    page.pdf(path=sys.argv[2], prefer_css_page_size=True, print_background=True,
             margin={"top":"0","bottom":"0","left":"0","right":"0"})
    browser.close()
'''
    subprocess.run([PLAYWRIGHT_PYTHON, "-c", browser_code, str(html), str(pdf1)], check=True)
    # All text, paths, numeric labels and arrow topology remain unchanged.
    svg_ns = "{http://www.w3.org/2000/svg}"
    before, after = ET.fromstring(original_svg), ET.fromstring(edited_svg)
    assert [t.text for t in before.iter(svg_ns + "text")] == [t.text for t in after.iter(svg_ns + "text")]
    audit = {
        "policy_url": "https://journals.plos.org/plosone/s/figures",
        "policy_checked_utc_date": "2026-10-26",
        "limits": {"max_width_mm": 190.5, "max_height_mm": 222.3,
                   "dpi": [300, 600], "max_bytes": 10_000_000,
                   "font_families": ["Arial", "Times", "Symbol"], "text_pt": [8, 12]},
        "font_probe": {name: font_match(name) for name in ("Arial", "Times New Roman", "Symbol")},
        "source_copies": preserved,
        "display_change": "Figure 1 A/B panel letters reduced from 20 to 16.8 SVG units (approximately 14.17 to 11.90 pt); all other SVG bytes unchanged.",
        "font_exception": "Arial/Times are not installed. Existing Liberation Sans in Figure 1 and DejaVu Sans in Figures 2-4 were preserved. TIFFs satisfy raster/container requirements, but journal named-family typography is not fully met.",
        "figures": [],
    }
    for number in range(1, 5):
        source_pdf = pdf1 if number == 1 else ARCHIVE / f"Figure_{number}.pdf"
        with pdfplumber.open(source_pdf) as document:
            page = document.pages[0]
            width_mm, height_mm = page.width * 25.4 / 72, page.height * 25.4 / 72
            sizes = Counter(round(c["size"], 4) for c in page.chars if c["text"].strip())
            fonts = Counter(c["fontname"] for c in page.chars if c["text"].strip())
            text_bounds = [min(c["x0"] for c in page.chars), min(c["top"] for c in page.chars),
                           max(c["x1"] for c in page.chars), max(c["bottom"] for c in page.chars)]
        document = pdfium.PdfDocument(source_pdf)
        page_image = document[0].render(scale=DPI / 72).to_pil().convert("RGB")
        target = OUT / f"Fig{number}.tif"
        page_image.save(target, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
        preview = VALID / f"main_figure_{number}_preview.png"
        preview_image = page_image.copy()
        preview_image.thumbnail((1100, 1300), Image.Resampling.LANCZOS)
        preview_image.save(preview)
        with Image.open(target) as check:
            checks = {
                "single_frame": check.n_frames == 1,
                "RGB_8_bit_channels": check.mode == "RGB" and tuple(check.tag_v2[258]) == (8, 8, 8),
                "no_alpha": 338 not in check.tag_v2,
                "LZW": check.tag_v2[259] == 5,
                "dpi_600": tuple(check.info["dpi"]) == (600.0, 600.0),
                "file_below_10MB": target.stat().st_size < 10_000_000,
                "dimensions_within_policy": 66.8 <= check.width / DPI * 25.4 <= 190.5 and check.height / DPI * 25.4 <= 222.3,
                "font_size_8_to_12pt": min(sizes) >= 8.0 and max(sizes) <= 12.0,
                "text_inside_page": text_bounds[0] >= 0 and text_bounds[1] >= 0 and text_bounds[2] <= page.width and text_bounds[3] <= page.height,
            }
            # Exact readback verifies compression does not alter rendered pixels.
            check.load()
            checks["lossless_pixel_readback"] = check.tobytes() == page_image.tobytes()
            assert all(checks.values()), (number, checks)
            record = {"number": number, "raster_source_pdf": str(source_pdf.relative_to(ROOT)),
                      "upload_file": str(target.relative_to(ROOT)), "sha256": sha(target),
                      "bytes": target.stat().st_size, "pixels": list(check.size),
                      "dpi": [float(value) for value in check.info["dpi"]], "physical_mm": [check.width / DPI * 25.4, check.height / DPI * 25.4],
                      "pdf_mm": [width_mm, height_mm], "actual_embedded_fonts": dict(fonts),
                      "text_sizes_pt": dict(sizes), "text_bounds_pt": text_bounds,
                      "checks": checks, "named_font_family_compliant": False,
                      "visual_preview": str(preview.relative_to(ROOT))}
        audit["figures"].append(record)
        print(f"Fig{number}: {record['pixels']}, {record['bytes']} bytes, all container/size checks passed")
    assert all(record["copy_identical"] for record in preserved)
    (VALID / "main_figure_export_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
