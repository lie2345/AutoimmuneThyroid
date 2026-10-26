#!/usr/bin/env python3
"""Prepare native-width PNG/PDF visual-check artifacts without changing plots."""
import hashlib
import json
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from PIL import Image

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[3]
FIG, REV = BASE / "figures", BASE / "reviews"
removed = json.loads((REV / "main_removed_text.json").read_text())["removed_strings"]
build = json.loads((REV / "main_build_validation.json").read_text())
checks, renders = [], []


def check(name, passed, detail=None):
    checks.append({"check": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(name)


for record in build["figures"]:
    slug = record["figure"]
    svg_root = ET.parse(FIG / (slug + ".svg")).getroot()
    texts = ["".join(e.itertext()) for e in svg_root.iter() if e.tag.endswith("}text")]
    check(slug + " all ten explanatory strings absent from actual SVG", all(r["text"] not in texts for r in removed))
    text = subprocess.run(["pdftotext", "-layout", str(FIG / (slug + ".pdf")), "-"], capture_output=True, text=True, check=True).stdout
    normalized = " ".join(text.split())
    check(slug + " all ten explanatory strings absent from actual PDF", all(" ".join(r["text"].split()) not in normalized for r in removed))
    if slug.startswith("figure1"):
        check(slug + " caption markers retained in SVG", all(any(s in t for t in texts) for s in ("RP08*", "RP09*")))
    if slug == "figure3":
        check(slug + " primary margin marker retained in SVG", "1.10*" in texts)
    prefix = REV / ("main_" + slug + "_pdf_native96")
    command = ["pdftoppm", "-r", "96", "-png", "-singlefile", str(FIG / (slug + ".pdf")), str(prefix)]
    subprocess.run(command, check=True, capture_output=True)
    pdf_png = prefix.with_suffix(".png")
    with Image.open(pdf_png) as image:
        pdf_pixels = image.size
    # 180 mm at 96 screen pixels/inch; PDF rasterization may add one edge pixel.
    width = round(180 / 25.4 * 96)
    height = round(record["height_mm"] / 25.4 * 96)
    png_native = REV / ("main_" + slug + "_png_native96.png")
    with Image.open(FIG / (slug + ".png")) as image:
        image.convert("RGB").resize((width, height), Image.Resampling.LANCZOS).save(png_native)
    check(slug + " native PDF page dimensions", abs(pdf_pixels[0] - width) <= 1 and abs(pdf_pixels[1] - height) <= 1, pdf_pixels)
    renders.append({"figure": slug, "width_mm": 180, "height_mm": record["height_mm"],
                    "native_screen_dpi": 96, "png_path": str(png_native.relative_to(ROOT)),
                    "pdf_render_path": str(pdf_png.relative_to(ROOT)), "png_pixels": [width, height],
                    "pdf_pixels": list(pdf_pixels), "pdf_render_command": command})
for name, expected in build["source_hashes_sha256"].items():
    check("protected input unchanged " + name, hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected)
(REV / "main_native_validation.json").write_text(json.dumps({"checks": checks, "renders": renders,
    "visual_status": "PENDING_IMPLEMENTATION_AGENT_VIEWING", "note": "96 dpi represents 180 mm physical width on a 96 dpi screen; screen scaling varies by display."}, indent=2) + "\n")
print(json.dumps({"checks": len(checks), "passed": sum(r["passed"] for r in checks), "renders": renders}, indent=2))
