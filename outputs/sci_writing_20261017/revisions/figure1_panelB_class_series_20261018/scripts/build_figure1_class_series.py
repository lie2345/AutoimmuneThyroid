#!/usr/bin/env python3
"""Display-only Figure 1B five-class count-series design alternatives.

The existing Matplotlib Figure 1A drawing block is copied verbatim below. All
scientific values and classifications are read from frozen tables and checked
against the original T06 outputs. No estimates, tests, or classifications run.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle, Wedge
from PIL import Image, ImageChops

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[3]
PREVIOUS = BASE.parent / "matplotlib_t19_t20_v2"
FROZEN = BASE.parent / "circular_v1" / "source_data"
FIG = BASE / "figures"
DATA = BASE / "source_data"
REV = BASE / "reviews"
for folder in (FIG, DATA, REV):
    folder.mkdir(parents=True, exist_ok=True)

FONT = "DejaVu Sans"
INK, MUTED, RULE = "#18212B", "#505965", "#D9DEE4"
BLUE, ORANGE, GRAY = "#0072B2", "#D55E00", "#777777"
MM = 72 / 25.4
WIDTH = 180 * MM
FULL_HEIGHT_MM = 245
B_OFFSET_PT = 220
B_HEIGHT_MM = FULL_HEIGHT_MM - B_OFFSET_PT / MM
plt.rcParams.update({
    "font.family": FONT, "font.size": 8, "text.color": INK,
    "axes.labelcolor": INK, "svg.fonttype": "none", "pdf.fonttype": 42,
    "ps.fonttype": 42, "savefig.facecolor": "white", "axes.unicode_minus": True,
    "path.simplify": False,
})
CLASS_ORDER = ["shared_strong_H4", "distinct_signals_H3", "shared_suggestive_H4",
               "unresolved_low_purity", "unresolved_insufficient_overlap"]
STYLES = {
    CLASS_ORDER[0]: dict(label="Strong H4", color=BLUE, marker="o", fill=BLUE),
    CLASS_ORDER[1]: dict(label="Strong H3", color=ORANGE, marker="D", fill=ORANGE),
    CLASS_ORDER[2]: dict(label="Suggestive H4", color=BLUE, marker="o", fill="white"),
    CLASS_ORDER[3]: dict(label="Low purity", color=GRAY, marker="x", fill=GRAY),
    CLASS_ORDER[4]: dict(label="Low overlap", color=GRAY, marker="^", fill="white"),
}
VARIANTS = {
    "A_five_class_lines": "Five separate class-count line series across 16 region pairs",
    "B_five_class_grouped_bars": "Five grouped radial class-count bars in each of 16 region pairs",
}
CHECKS, FIGURES = [], []


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def check(name, passed, detail=None):
    CHECKS.append({"check": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(name, rows):
    with (DATA / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def hashes(paths):
    return {rel(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def key(row):
    return row["pair_id"], row["graves_cs"], row["aiht_cs"]


source_paths = [
    PREVIOUS / "scripts/build_main_figures.py",
    PREVIOUS / "source_data/main_figure1_marks.tsv",
    PREVIOUS / "figures/figure1.png", PREVIOUS / "figures/figure1.svg", PREVIOUS / "figures/figure1.pdf",
    PREVIOUS / "figure_legends.md",
    FROZEN / "figure1_all42_comparisons.tsv", FROZEN / "figure1_region_counts.tsv",
    FROZEN / "figure1_endpoints.tsv",
    ROOT / "outcome_data/t06/gwas_gwas_susie_coloc.tsv",
    ROOT / "outcome_data/t06/gwas_gwas_region_architecture.tsv",
]
before = hashes(source_paths)
f1 = read_tsv(PREVIOUS / "source_data/main_figure1_marks.tsv")
frozen_pairs = read_tsv(FROZEN / "figure1_all42_comparisons.tsv")
regions = read_tsv(FROZEN / "figure1_region_counts.tsv")
endpoints = read_tsv(FROZEN / "figure1_endpoints.tsv")
raw_pairs_list = read_tsv(ROOT / "outcome_data/t06/gwas_gwas_susie_coloc.tsv")
raw_pairs = {key(row): row for row in raw_pairs_list}
raw_regions_list = read_tsv(ROOT / "outcome_data/t06/gwas_gwas_region_architecture.tsv")
raw_regions = {r["pair_id"]: r for r in raw_regions_list}
check("42 unique frozen comparison identities", len(f1) == len({key(r) for r in f1}) == 42)
check("frozen comparison order matches previous source", [key(r) for r in f1] == [key(r) for r in frozen_pairs])
check("16 regions in original RP order", [r["pair_id"] for r in regions] == [f"RP{i:02d}" for i in range(1, 17)])
check("five source classes have counts 8,11,1,17,5",
      [Counter(r["classification"] for r in f1)[c] for c in CLASS_ORDER] == [8, 11, 1, 17, 5])
check("non-HLA comparisons retain chrX", all(r["is_hla"] == "FALSE" for r in f1)
      and any(r["pair_id"] == "RP16" and r["chromosome"] == "X" for r in f1))
check("RP08 and RP09 share one Graves region",
      raw_regions["RP08"]["graves_region"] == raw_regions["RP09"]["graves_region"])
for row in f1:
    original = raw_pairs[key(row)]
    check("all original fields unchanged " + "/".join(key(row)),
          all(row[field] == value for field, value in original.items()))
    check("original source line retained " + "/".join(key(row)),
          raw_pairs_list[int(row["source_line"]) - 2] == original)
for row in regions:
    pid, original = row["pair_id"], raw_regions[row["pair_id"]]
    counts = Counter(r["classification"] for r in f1 if r["pair_id"] == pid)
    check("region source fields unchanged " + pid,
          all(row[field] == value for field, value in original.items()))
    check("region comparison count equals stored total " + pid,
          sum(counts.values()) == int(original["n_signal_pair_tests"]))
    check("all five region class counts equal frozen totals " + pid,
          all(counts[c] == int(row[c]) for c in CLASS_ORDER))
check("region totals are 4,3,1,4,2,8,3,1,2,2,1,4,2,2,2,1",
      [int(r["n_signal_pair_tests"]) for r in regions] == [4, 3, 1, 4, 2, 8, 3, 1, 2, 2, 1, 4, 2, 2, 2, 1])


class Drawing:
    """Physical-point canvas matching the preceding Figure 1 renderer."""
    def __init__(self, slug, height_mm, title, offset=0):
        self.slug, self.height_mm, self.height = slug, height_mm, height_mm * MM
        self.fig = plt.figure(figsize=(180 / 25.4, height_mm / 25.4), dpi=100)
        self.ax = self.fig.add_axes([0, 0, 1, 1], xlim=(0, WIDTH), ylim=(self.height + offset, offset))
        self.ax.set_axis_off()
        self.texts, self.title = [], title

    def text(self, x, y, label, size=8, color=INK, bold=False, ha="left"):
        item = self.ax.text(x, y, str(label), fontsize=size, color=color,
                            fontweight="bold" if bold else "normal", ha=ha, va="baseline",
                            clip_on=False, zorder=6)
        self.texts.append(item)
        return item

    def line(self, x1, y1, x2, y2, color=RULE, width=.6, dash="solid", zorder=1):
        self.ax.add_line(Line2D([x1, x2], [y1, y2], color=color, linewidth=width,
                                linestyle=dash, solid_capstyle="butt", zorder=zorder))

    def rect(self, x, y, width, height, fill="white", stroke=RULE, radius=0):
        if radius:
            patch = FancyBboxPatch((x, y), width, height,
                                   boxstyle=f"round,pad=0,rounding_size={radius}",
                                   facecolor=fill, edgecolor=stroke, linewidth=.6, zorder=0)
        else:
            patch = Rectangle((x, y), width, height, facecolor=fill, edgecolor=stroke,
                              linewidth=.6, zorder=0)
        self.ax.add_patch(patch)

    def circle(self, x, y, radius, color=RULE, width=.6, zorder=1):
        self.ax.add_patch(Circle((x, y), radius, fill=False, edgecolor=color, linewidth=width, zorder=zorder))

    def symbol(self, x, y, style, radius=2.6, zorder=4):
        return self.ax.plot([x], [y], marker=style["marker"], markersize=2 * radius,
                            markerfacecolor=style["fill"], markeredgecolor=style["color"],
                            markeredgewidth=.8, color=style["color"], linestyle="none", zorder=zorder)[0]

    def panel(self, letter, title, y):
        self.text(16, y, letter, 11, bold=True)
        self.text(32, y, title, 9.5, bold=True)

    def finish(self):
        self.fig.canvas.draw()
        renderer, extent = self.fig.canvas.get_renderer(), self.fig.bbox
        outliers = []
        for item in self.texts:
            bbox = item.get_window_extent(renderer)
            if not (bbox.x0 >= extent.x0 and bbox.x1 <= extent.x1 and bbox.y0 >= extent.y0 and bbox.y1 <= extent.y1):
                outliers.append(item.get_text())
        check(self.slug + " all text within physical canvas", not outliers, outliers)
        check(self.slug + " all text at least 8 pt", min(t.get_fontsize() for t in self.texts) >= 8)
        check(self.slug + " no ambiguous n= labels", not any("n=" in t.get_text() for t in self.texts))
        metadata = {"Title": self.title, "Creator": "Matplotlib display-only alternatives; frozen T06 outputs"}
        self.fig.savefig(FIG / (self.slug + ".svg"), format="svg", metadata=metadata)
        self.fig.savefig(FIG / (self.slug + ".pdf"), format="pdf", metadata=metadata)
        png_info = {}
        for dpi, suffix in ((600, ""), (180, "_preview")):
            path = FIG / f"{self.slug}{suffix}.png"
            self.fig.savefig(path, format="png", dpi=dpi, metadata={"Title": self.title, "Software": "Matplotlib"})
            with Image.open(path) as image:
                png_info[str(dpi)] = {"path": rel(path), "pixels": list(image.size), "dpi": list(image.info["dpi"])}
                check(self.slug + f" PNG {dpi} dimensions",
                      image.size == (int(180 / 25.4 * dpi), int(self.height_mm / 25.4 * dpi)), list(image.size))
                check(self.slug + f" PNG {dpi} dpi metadata", all(abs(v - dpi) < .02 for v in image.info["dpi"]))
        svg_root = ET.parse(FIG / (self.slug + ".svg")).getroot()
        tags = [e.tag.rsplit("}", 1)[-1] for e in svg_root.iter()]
        check(self.slug + " editable SVG text retained", tags.count("text") == len(self.texts))
        check(self.slug + " no raster images in SVG", "image" not in tags)
        check(self.slug + " SVG physical dimensions", abs(float(svg_root.attrib["width"][:-2]) - WIDTH) < 1e-5
              and abs(float(svg_root.attrib["height"][:-2]) - self.height) < 1e-5,
              {"width": svg_root.attrib["width"], "height": svg_root.attrib["height"]})
        fonts = subprocess.run(["pdffonts", str(FIG / (self.slug + ".pdf"))], text=True, capture_output=True, check=True).stdout
        check(self.slug + " PDF fonts embedded", len(fonts.splitlines()) > 2 and all("yes" in line for line in fonts.splitlines()[2:]), fonts)
        check(self.slug + " no raster images in PDF", b"/Subtype /Image" not in (FIG / (self.slug + ".pdf")).read_bytes())
        FIGURES.append({"figure": self.slug, "width_mm": 180, "height_mm": self.height_mm,
                        "minimum_font_pt": min(t.get_fontsize() for t in self.texts),
                        "svg_text_count": tags.count("text"), "png": png_info, "pdf_fonts": fonts})
        plt.close(self.fig)


def draw_panel_a(d):
    # The drawing statements below are verbatim from build_main_figures.py.
    d.panel("A", "Scope of the existing FinnGen R13 comparison", 23)
    for x, r, title in zip((16, 274), endpoints, ("Graves’ disease", "Autoimmune hypothyroidism")):
        d.rect(x, 40, 220, 72, "#F8FAFC", RULE, 3)
        d.text(x + 10, 57, title, 9, bold=True)
        d.text(x + 10, 73, r["trait_id"])
        d.text(x + 10, 91, f'{int(r["n_cases"]):,} cases / {int(r["n_controls"]):,} controls', 8.5)
        d.ax.add_patch(FancyArrowPatch((x + 110, 113), (x + 110, 136), arrowstyle="-|>",
                                       mutation_scale=8, linewidth=.9, color=MUTED))
    d.rect(16, 136, 478, 61, "#F3F7FA", "#A9BBC9", 3)
    d.text(WIDTH / 2, 153, "16 overlapping non-HLA fine-mapping region pairs", 9, bold=True, ha="center")
    d.text(WIDTH / 2, 169, "15 unique Graves regions and 16 unique AIHT regions", 8.5, ha="center")
    d.text(WIDTH / 2, 185, "42 pairwise comparisons of detected association components", 8.5, ha="center")
    d.text(16, 214, "Finnish-enriched European ancestry · GRCh38 · summary statistics", color=MUTED)


def xy(cx, cy, radius, angle):
    return cx + radius * math.cos(angle), cy - radius * math.sin(angle)


# The zero-count circle has a nonzero geometric radius so that all 47 zero
# values remain visible after a fixed class-specific angular layout offset.
ZERO_RADIUS = 85
COUNT_MAX = 4
OUTER_RADIUS = 143
COUNT_STEP = (OUTER_RADIUS - ZERO_RADIUS) / COUNT_MAX
CLASS_OFFSETS_DEG = [8, 4, 0, -4, -8]
BAR_WIDTH_DEG = 2
DASHES = ["solid", "solid", (0, (1, 1.5)), (0, (4, 2)), (0, (1, 2))]
BAR_FILLS = [BLUE, ORANGE, "white", GRAY, "white"]
BAR_HATCHES = [None, None, "///", None, "..."]
for index, classification in enumerate(CLASS_ORDER):
    STYLES[classification].update(dash=DASHES[index], bar_fill=BAR_FILLS[index],
                                  hatch=BAR_HATCHES[index], offset_degrees=CLASS_OFFSETS_DEG[index])


def comparison_id(row):
    return f'{row["pair_id"]}:G{row["graves_cs"]}:A{row["aiht_cs"]}'


cell_rows, comparison_rows, region_rows = [], [], []
for i, region in enumerate(regions):
    pid = region["pair_id"]
    subset = [r for r in f1 if r["pair_id"] == pid]
    for j, classification in enumerate(CLASS_ORDER):
        members = [r for r in subset if r["classification"] == classification]
        raw_members = [r for r in raw_pairs_list if r["pair_id"] == pid and r["classification"] == classification]
        cell_id = pid + ":" + classification
        check("80-cell exact raw comparison set " + cell_id,
              {key(r) for r in members} == {key(r) for r in raw_members})
        check("80-cell count equals frozen architecture " + cell_id,
              len(members) == int(region[classification]))
        cell_rows.append({
            "cell_id": cell_id, "pair_id": pid, "chromosome": region["chromosome"],
            "region_order": i + 1, "class_order": j + 1, "classification": classification,
            "class_label": STYLES[classification]["label"], "comparison_count": len(members),
            "region_total_comparisons": len(subset),
            "comparison_ids": ";".join(comparison_id(r) for r in members),
            "comparison_source_path": "outcome_data/t06/gwas_gwas_susie_coloc.tsv",
            "comparison_source_lines": ";".join(r["source_line"] for r in members),
            "architecture_source_path": region["source_path"], "architecture_source_line": region["source_line"],
            "frozen_architecture_source_path": rel(FROZEN / "figure1_region_counts.tsv"),
            "frozen_architecture_source_line": i + 2,
            "is_chr_x": region["is_chr_x"], "shared_graves_region_star": pid in ("RP08", "RP09"),
        })
    region_rows.append(dict(region, comparison_label=f'{len(subset)} comparison' + ("" if len(subset) == 1 else "s")))
for row in f1:
    comparison_rows.append(dict(row, comparison_id=comparison_id(row),
                               count_cell_id=row["pair_id"] + ":" + row["classification"]))
check("all 80 class-by-region cells including zeros", len(cell_rows) == len({r["cell_id"] for r in cell_rows}) == 80)
check("33 nonzero and 47 zero cells retained", sum(r["comparison_count"] > 0 for r in cell_rows) == 33
      and sum(r["comparison_count"] == 0 for r in cell_rows) == 47)
check("actual maximum single-class count is four", max(r["comparison_count"] for r in cell_rows) == COUNT_MAX == 4)
check("five column totals equal 8,11,1,17,5", [sum(r["comparison_count"] for r in cell_rows if r["classification"] == c)
      for c in CLASS_ORDER] == [8, 11, 1, 17, 5])
check("all 16 row totals equal stored region totals", all(sum(r["comparison_count"] for r in cell_rows if r["pair_id"] == q["pair_id"])
      == int(q["n_signal_pair_tests"]) for q in regions))
check("all 42 identities map exactly once to count cells", len(comparison_rows) == len({r["comparison_id"] for r in comparison_rows}) == 42
      and sum(r["comparison_count"] for r in cell_rows) == 42)
check("suggestive H4 is displayed only as RP16 count one", [(r["pair_id"], r["comparison_count"]) for r in cell_rows
      if r["classification"] == CLASS_ORDER[2] and r["comparison_count"]] == [("RP16", 1)])
write_tsv("class_counts_80.tsv", cell_rows)
write_tsv("comparison_sources_42.tsv", comparison_rows)
write_tsv("region_totals_16.tsv", region_rows)


def draw_panel_b(d, variant):
    is_line = variant == "A_five_class_lines"
    title = "Five class-count series across the 16 non-HLA region pairs" if is_line else "Five grouped class-count bars in each non-HLA region pair"
    d.panel("B", title, 239)
    d.text(16, 255, "Five classes per region pair; common radial scale: 0–4 comparisons.", color=MUTED)
    d.text(16, 267, "Counts describe tested CS pairs, not people or genome-wide sharing rates.", color=MUTED)
    cx, cy, label_radius = WIDTH / 2, 450, 164
    tick_texts, data_artists, bottom_labels, region_texts, geometry, segments = [], [], [], [], [], []
    for tick in range(COUNT_MAX + 1):
        radius = ZERO_RADIUS + COUNT_STEP * tick
        d.circle(cx, cy, radius, "#A7B5BF" if tick == 0 else "#D9E2E7", .7 if tick == 0 else .5, 0)
    for i, region in enumerate(regions):
        pid = region["pair_id"]
        angle = math.pi / 2 - (i + .5) * 2 * math.pi / 16
        boundary = math.pi / 2 - i * 2 * math.pi / 16
        d.line(*xy(cx, cy, ZERO_RADIUS, boundary), *xy(cx, cy, OUTER_RADIUS, boundary), "#E2E7EB", .45)
        region_label_radius = 176 if pid == "RP06" else label_radius
        lx, ly = xy(cx, cy, region_label_radius, angle)
        star = "*" if pid in ("RP08", "RP09") else ""
        region_texts.append(d.text(lx, ly, f'{pid}{star} · chr{region["chromosome"]}', bold=True, ha="center"))
        total = int(region["n_signal_pair_tests"])
        label = f'{total} comparison' + ("" if total == 1 else "s")
        total_text = d.text(lx, ly + 11, label, color=MUTED, ha="center")
        region_texts.append(total_text)
        if pid in ("RP08", "RP09"):
            bottom_labels.append(total_text)
        for cell in [r for r in cell_rows if r["pair_id"] == pid]:
            style = STYLES[cell["classification"]]
            theta = angle + math.radians(style["offset_degrees"])
            radius = ZERO_RADIUS + COUNT_STEP * cell["comparison_count"]
            x, y = xy(cx, cy, radius, theta)
            geometry.append(dict(cell, variant=variant, region_center_degrees=math.degrees(angle),
                                 region_label_radius_pt=region_label_radius,
                                 class_angular_offset_degrees=style["offset_degrees"],
                                 displayed_angle_degrees=math.degrees(theta), zero_radius_pt=ZERO_RADIUS,
                                 points_per_comparison=COUNT_STEP, radius_pt=radius, x_pt=x, y_pt=y,
                                 marker=style["marker"], marker_fill=style["fill"], color=style["color"],
                                 line_dash=str(style["dash"]),
                                 bar_angular_width_degrees=0 if is_line else BAR_WIDTH_DEG,
                                 bar_fill=style["bar_fill"], bar_hatch=style["hatch"] or "none"))
    if is_line:
        for classification in (CLASS_ORDER[3], CLASS_ORDER[4], CLASS_ORDER[1], CLASS_ORDER[0], CLASS_ORDER[2]):
            rows = [r for r in geometry if r["classification"] == classification]
            style = STYLES[classification]
            d.ax.add_line(Line2D([r["x_pt"] for r in rows], [r["y_pt"] for r in rows],
                                color=style["color"], linewidth=1.05, linestyle=style["dash"],
                                solid_capstyle="butt", zorder=2))
            for first, second in zip(rows[:-1], rows[1:]):
                segments.append({"classification": classification, "from_cell_id": first["cell_id"],
                                 "to_cell_id": second["cell_id"], "from_count": first["comparison_count"],
                                 "to_count": second["comparison_count"], "from_x_pt": first["x_pt"],
                                 "from_y_pt": first["y_pt"], "to_x_pt": second["x_pt"], "to_y_pt": second["y_pt"]})
        for row in geometry:
            data_artists.append(d.symbol(row["x_pt"], row["y_pt"], STYLES[row["classification"]], 1.6, 4))
        check(variant + " five open discrete polylines with 75 adjacent segments", len(segments) == 75
              and all(s["from_cell_id"].split(":")[0] != "RP16" for s in segments))
    else:
        for row in geometry:
            style = STYLES[row["classification"]]
            radius, angle = row["radius_pt"], row["displayed_angle_degrees"]
            if row["comparison_count"]:
                bar = Wedge((cx, cy), radius, -angle - BAR_WIDTH_DEG / 2, -angle + BAR_WIDTH_DEG / 2,
                            width=radius - ZERO_RADIUS, facecolor=style["bar_fill"],
                            edgecolor=style["color"], hatch=style["hatch"], linewidth=.65, zorder=2)
                d.ax.add_patch(bar)
                data_artists.append(bar)
            else:
                # Zero-height bars have an explicit class-specific baseline symbol.
                data_artists.append(d.symbol(row["x_pt"], row["y_pt"], style, 1.6, 4))
        check(variant + " five parallel bar positions per region", all(
            len({r["displayed_angle_degrees"] for r in geometry if r["pair_id"] == q["pair_id"]}) == 5 for q in regions))
        check(variant + " all bars share zero ring with exact nonnegative lengths", all(
            r["radius_pt"] - ZERO_RADIUS == COUNT_STEP * r["comparison_count"] for r in geometry))
    d.line(cx, cy - ZERO_RADIUS, cx, cy - OUTER_RADIUS, "#7C929F", .5, zorder=3)
    for tick in range(COUNT_MAX + 1):
        yy = cy - (ZERO_RADIUS + COUNT_STEP * tick)
        d.line(cx - 2, yy, cx + 2, yy, "#7C929F", .6, zorder=3)
        text = d.text(cx, yy + 2.6, str(tick), color=INK, ha="center")
        text.set_bbox(dict(facecolor="white", edgecolor="none", pad=.1))
        tick_texts.append(text)
    d.text(cx, 441, "42 comparisons", 10, bold=True, ha="center")
    d.text(cx, 457, "16 region pairs", 8.5, ha="center")
    d.text(cx, 472, "5 classes", 8, ha="center")
    totals_text = d.text(16, 633, "42 comparisons: 8 strong H4 · 11 strong H3 · 1 suggestive H4 · 17 low purity · 5 low overlap.")
    d.line(16, 639, 494, 639)
    for classification, x in zip(CLASS_ORDER, (18, 112, 206, 311, 398)):
        style = STYLES[classification]
        if is_line:
            d.line(x, 650, x + 15, 650, style["color"], 1.05, style["dash"])
            d.symbol(x + 7.5, 650, style, 1.8)
        else:
            d.ax.add_patch(Rectangle((x + 1, 646), 11, 7, facecolor=style["bar_fill"],
                                    edgecolor=style["color"], hatch=style["hatch"], linewidth=.65))
            d.symbol(x + 18, 650, style, 1.6)
        d.text(x + (20 if is_line else 25), 653, style["label"])
    if is_line:
        d.text(16, 667, "Class-specific angular offsets separate coincident values; all radial counts remain exact.", color=MUTED)
        d.text(16, 679, "Connections guide the eye only; seam open. Angles do not encode genomic distance or causality.", color=MUTED)
    else:
        d.text(16, 667, "Five positions per sector encode class in legend order; radial length = comparisons.", color=MUTED)
        d.text(16, 679, "Zero counts have a baseline symbol. Angles do not encode genomic distance or causality.", color=MUTED)
    d.text(16, 691, "RP order is organizational. * RP08/RP09 share a Graves region. The chrX pair is retained.", color=MUTED)
    d.fig.canvas.draw()
    renderer = d.fig.canvas.get_renderer()
    totals_box = totals_text.get_window_extent(renderer)
    bottom_gaps = [(t.get_window_extent(renderer).y0 - totals_box.y1) * 72 / d.fig.dpi for t in bottom_labels]
    check(variant + " bottom labels clear of totals text", min(bottom_gaps) >= 2,
          {"minimum_vertical_gap_pt": min(bottom_gaps)})
    tick_collisions = [t.get_text() for t in tick_texts if any(
        t.get_window_extent(renderer).padded(.1 * d.fig.dpi / 72).overlaps(artist.get_window_extent(renderer))
        for artist in data_artists)]
    check(variant + " seam count labels clear of data markers or bars", not tick_collisions, tick_collisions)
    region_collisions = [{"label": t.get_text(), "cell": cell["cell_id"]} for t in region_texts
                         for artist, cell in zip(data_artists, geometry)
                         if t.get_window_extent(renderer).padded(.3 * d.fig.dpi / 72).overlaps(
                             artist.get_window_extent(renderer))]
    check(variant + " region titles and totals clear of all data markers or bars", not region_collisions, region_collisions)
    check(variant + " all 80 displayed cell artists retained", len(data_artists) == 80)
    if is_line:
        marker_overlaps = [(geometry[i]["cell_id"], geometry[j]["cell_id"])
                           for i in range(80) for j in range(i + 1, 80)
                           if data_artists[i].get_window_extent(renderer).overlaps(data_artists[j].get_window_extent(renderer))]
        check(variant + " all 80 class markers separated including zeros", not marker_overlaps, marker_overlaps)
    check(variant + " 80 source cells mapped to 80 plotted values", len(geometry) == 80
          and {r["cell_id"] for r in geometry} == {r["cell_id"] for r in cell_rows})
    check(variant + " radial count values are exact for all 80 cells", all(
        math.isclose(math.hypot(r["x_pt"] - cx, r["y_pt"] - cy), ZERO_RADIUS + COUNT_STEP * r["comparison_count"], abs_tol=1e-10)
        for r in geometry))
    check(variant + " rare RP16 suggestive H4 has visible nonzero radius", any(
        r["pair_id"] == "RP16" and r["classification"] == CLASS_ORDER[2]
        and r["comparison_count"] == 1 and r["radius_pt"] == ZERO_RADIUS + COUNT_STEP for r in geometry))
    return geometry, segments


for variant, title in VARIANTS.items():
    for view in ("full", "panelB"):
        slug = f"figure1_{variant}" + ("_panelB" if view == "panelB" else "")
        d = Drawing(slug, B_HEIGHT_MM if view == "panelB" else FULL_HEIGHT_MM, title,
                    offset=B_OFFSET_PT if view == "panelB" else 0)
        if view == "full":
            draw_panel_a(d)
        geometry, segments = draw_panel_b(d, variant)
        if view == "full":
            write_tsv(slug + "_geometry_80.tsv", geometry)
            if segments:
                write_tsv(slug + "_segments_75.tsv", segments)
        d.finish()

with Image.open(PREVIOUS / "figures/figure1.png") as original:
    panel_a_box = (0, 0, original.width, int(B_OFFSET_PT / 72 * 600))
    reference = original.crop(panel_a_box).convert("RGB")
    for variant in VARIANTS:
        with Image.open(FIG / f"figure1_{variant}.png") as current:
            difference = ImageChops.difference(reference, current.crop(panel_a_box).convert("RGB"))
            check(variant + " Panel A pixel-identical to preceding 600dpi Figure 1", difference.getbbox() is None,
                  {"pixel_crop": list(panel_a_box)})

after = hashes(source_paths)
check("all authoritative inputs and preceding Figure 1 files unchanged", before == after)
manifest = {
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "scope": "Display-only Figure 1B five-class count series; frozen values and classifications; no statistical rerun.",
    "command": "uv run --no-sync --python /usr/bin/python3 python outputs/sci_writing_20261017/revisions/figure1_panelB_class_series_20261018/scripts/build_figure1_class_series.py",
    "runtime": {"python": sys.version, "executable": sys.executable, "matplotlib": matplotlib.__version__, "platform": platform.platform()},
    "variants": VARIANTS,
    "count_table": {"cells": 80, "nonzero": 33, "zero": 47, "comparison_identities": 42,
                    "class_order": CLASS_ORDER, "class_totals": [8, 11, 1, 17, 5]},
    "radial_count_scale": {"minimum": 0, "maximum": COUNT_MAX, "ticks": list(range(COUNT_MAX + 1)),
                           "zero_radius_pt": ZERO_RADIUS, "points_per_comparison": COUNT_STEP},
    "class_angular_offsets_degrees": dict(zip(CLASS_ORDER, CLASS_OFFSETS_DEG)),
    "grouped_bar_angular_width_degrees": BAR_WIDTH_DEG,
    "line_connections": {"series": 5, "segments_per_series": 15, "RP16_to_RP01_connected": False,
                         "smoothing": False, "meaning": "guide the eye only"},
    "protected_sources_before_sha256": before, "protected_sources_after_sha256": after,
    "figures": FIGURES, "checks": CHECKS,
    "summary": {"passed": sum(c["passed"] for c in CHECKS), "total": len(CHECKS)},
    "independent_review": "Not performed by this builder; primary agent is arranging independent figure review.",
}
(REV / "build_validation.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps(manifest["summary"]))


