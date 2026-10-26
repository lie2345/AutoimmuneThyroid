#!/usr/bin/env python3
"""Display-only Matplotlib redraw of the approved main Figures 1–3.

Frozen source tables and original T06/T07 outputs are read, never modified.
Only display coordinates and the already defined beta +/- 1.96*SE intervals
are calculated. No statistical models, classifications, or tests are rerun.
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
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, PathPatch, Rectangle, Wedge
from matplotlib.path import Path as MplPath
from matplotlib.text import Text
from PIL import Image

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[3]
OLD = BASE.parent / "circular_v1"
PREVIOUS = BASE.parent / "matplotlib_t19_t20_v2"
APPROVED = BASE.parent / "figure1_panelB_class_series_20261018"
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
plt.rcParams.update({
    "font.family": FONT, "font.size": 8, "text.color": INK,
    "axes.labelcolor": INK, "svg.fonttype": "none", "pdf.fonttype": 42,
    "ps.fonttype": 42, "savefig.facecolor": "white", "axes.unicode_minus": True,
    "path.simplify": False,
})
CLASS_ORDER = ["shared_strong_H4", "distinct_signals_H3", "shared_suggestive_H4",
               "unresolved_low_purity", "unresolved_insufficient_overlap"]
STYLES = {
    CLASS_ORDER[0]: dict(label="Strong H4", color=BLUE, marker="o", fill=BLUE, dash="solid"),
    CLASS_ORDER[1]: dict(label="Strong H3", color=ORANGE, marker="D", fill=ORANGE, dash="solid"),
    CLASS_ORDER[2]: dict(label="Suggestive H4", color=BLUE, marker="o", fill="white", dash="solid"),
    CLASS_ORDER[3]: dict(label="Low purity", color=GRAY, marker="x", fill=GRAY, dash=(0, (4, 2))),
    CLASS_ORDER[4]: dict(label="Low overlap", color=GRAY, marker="^", fill="white", dash=(0, (1, 2))),
}
CHECKS: list[dict] = []
REMOVED_TEXT = []
FIGURES: list[dict] = []


def check(name, passed, detail=None):
    CHECKS.append({"check": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(name, rows):
    path = DATA / name
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def hashes(paths):
    return {rel(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def pair_key(row):
    return row["pair_id"], row["graves_cs"], row["aiht_cs"]


def same_source_fields(label, frozen, original, extra=()):
    # Frozen rows extend original source rows only with provenance/display fields.
    fields = [key for key in original if key in frozen]
    check(label, all(frozen[key] == original[key] for key in fields),
          {key: [frozen[key], original[key]] for key in fields if frozen[key] != original[key]})
    for field in extra:
        check(label + " " + field, field in frozen and field in original)


SOURCE_NAMES = ["figure1_all42_comparisons.tsv", "figure1_region_counts.tsv", "figure1_endpoints.tsv",
                "figure2_all24_comparisons.tsv", "figure2_nodes.tsv", "figure3_marginal_estimates.tsv",
                "figure3_equivalence_sensitivity.tsv"]
source_paths = [OLD / "source_data" / name for name in SOURCE_NAMES]
source_paths += [ROOT / name for name in (
    "outcome_data/t06/gwas_gwas_susie_coloc.tsv", "outcome_data/t06/susie_signal_catalog.tsv",
    "outcome_data/t06/gwas_gwas_region_architecture.tsv", "config/traits.tsv",
    "outcome_data/t07/genetic_signal_five_classification.tsv",
    "outcome_data/t07/classification_equivalence_sensitivity.tsv",
    "outputs/locus_novelty_20261017/candidate_snapshot.tsv",
    "outputs/manuscript_increment_t19_t20_20261018/assessment.md",
    "outputs/sci_writing_20261017/scripts/build_figures.py")]
source_paths += [OLD / "figure_legends.md", OLD / "scripts/build_circular_figures.R"]
source_paths += [OLD / "figures" / f"figure{i}.{ext}" for i in (1, 2, 3) for ext in ("png", "svg")]
source_paths += [PREVIOUS / "scripts/build_main_figures.py", APPROVED / "scripts/build_figure1_class_series.py",
                 APPROVED / "source_data/class_counts_80.tsv",
                 APPROVED / "source_data/figure1_B_five_class_grouped_bars_geometry_80.tsv"]
source_paths += [PREVIOUS / "figures" / f"figure{i}.{ext}" for i in (1, 2, 3) for ext in ("png", "svg", "pdf")]
source_paths += [APPROVED / "figures" / f"figure1_B_five_class_grouped_bars.{ext}" for ext in ("png", "svg", "pdf")]
before = hashes(source_paths)
f1, regions, endpoints, f2, nodes, f3, sensitivity = [read_tsv(OLD / "source_data" / name) for name in SOURCE_NAMES]
raw_pairs = {pair_key(row): row for row in read_tsv(ROOT / "outcome_data/t06/gwas_gwas_susie_coloc.tsv")}
raw_nodes = {(r["trait_id"], r["region"], r["cs"]): r
             for r in read_tsv(ROOT / "outcome_data/t06/susie_signal_catalog.tsv")}
raw_regions = {r["pair_id"]: r for r in read_tsv(ROOT / "outcome_data/t06/gwas_gwas_region_architecture.tsv")}
raw_traits = {r["trait_id"]: r for r in read_tsv(ROOT / "config/traits.tsv")}
units = read_tsv(ROOT / "outcome_data/t07/genetic_signal_five_classification.tsv")
unit_index = {r["signal_unit_id"]: r for r in units}
raw_sensitivity = {(r["signal_unit_id"], float(r["equivalence_or_upper"])): r
                   for r in read_tsv(ROOT / "outcome_data/t07/classification_equivalence_sensitivity.tsv")}
snapshot = {r["signal_unit_id"]: r for r in read_tsv(ROOT / "outputs/locus_novelty_20261017/candidate_snapshot.tsv")}

check("Figure 1 all 42 unique comparisons", len(f1) == len({pair_key(r) for r in f1}) == 42)
check("Figure 1 all 16 region pairs", [r["pair_id"] for r in regions] == [f"RP{i:02d}" for i in range(1, 17)])
check("Figure 1 class totals", [Counter(r["classification"] for r in f1)[c] for c in CLASS_ORDER] == [8, 11, 1, 17, 5])
check("15 unique Graves and 16 unique AIHT regions", len({r["graves_region"] for r in f1}) == 15 and len({r["aiht_region"] for r in f1}) == 16)
check("Figure 1 non-HLA scope including chrX", all(r["is_hla"] == "FALSE" for r in f1) and any(r["pair_id"] == "RP16" and r["chromosome"] == "X" for r in f1))
for r in f1:
    same_source_fields("F1 source equality " + "/".join(pair_key(r)), r, raw_pairs[pair_key(r)])
for r in regions:
    same_source_fields("F1 region equality " + r["pair_id"], r, raw_regions[r["pair_id"]])
    subset = [p for p in f1 if p["pair_id"] == r["pair_id"]]
    check("F1 per-region class count " + r["pair_id"], int(r["n_signal_pair_tests"]) == len(subset) and
          all(int(r[c]) == Counter(p["classification"] for p in subset)[c] for c in CLASS_ORDER))
for r in endpoints:
    same_source_fields("F1 endpoint equality " + r["trait_id"], r, raw_traits[r["trait_id"]])
check("Figure 2 all 24 unique edges", len(f2) == len({pair_key(r) for r in f2}) == 24)
check("Figure 2 all 27 unique nodes", len(nodes) == len({(r["pair_id"], r["node_key"]) for r in nodes}) == 27)
check("Figure 2 exact five-class totals", [Counter(r["classification"] for r in f2)[c] for c in CLASS_ORDER] == [6, 9, 0, 7, 2])
check("Figure 2 no HLA or ophthalmopathy nodes", all(r["is_hla"] == "FALSE" and r["trait_id"] in raw_traits and r["trait_id"] in {"E4_GRAVES_STRICT", "E4_HYTHY_AI_STRICT"} for r in nodes))
for r in f2:
    same_source_fields("F2 source equality " + "/".join(pair_key(r)), r, raw_pairs[pair_key(r)])
for r in nodes:
    same_source_fields("F2 full-region node equality " + r["pair_id"] + "/" + r["node_key"], r,
                       raw_nodes[(r["trait_id"], r["region"], r["cs"])])
check("RP06 five unresolved comparisons retained", {(r["graves_cs"], r["aiht_cs"]) for r in f2 if r["pair_id"] == "RP06" and r["classification"].startswith("unresolved")} == {("1", "2"), ("2", "1"), ("2", "2"), ("2", "3"), ("2", "4")})
check("All 15 strong F2 edges use local Graves CS1", all(r["graves_cs"] == "1" for r in f2 if r["classification"] in CLASS_ORDER[:2]))
priority = [r for r in units if r["five_class_primary"] != "unresolved"]
priority.sort(key=lambda r: int(r["signal_unit_id"][3:]))
check("304 constructed units and original classification counts", len(units) == 304 and Counter(r["five_class_primary"] for r in units) == {"shared-concordant": 6, "shared-opposite": 2, "Graves-predominant": 2, "unresolved": 294})
check("Figure 3 exact ten candidate set", len(priority) == len(snapshot) == 10 and {r["signal_unit_id"] for r in priority} == set(snapshot) == {r["signal_unit_id"] for r in f3})
check("Figure 3 twenty marginal estimates", len(f3) == len({(r["signal_unit_id"], r["trait"]) for r in f3}) == 20)
for r in f3:
    sid, trait = r["signal_unit_id"], r["trait"]
    q, s = unit_index[sid], snapshot[sid]
    b, se = float(r["beta"]), float(r["se"])
    check("F3 identity and stored estimate " + sid + "/" + trait,
          r["representative_variant"] == q["representative_variant"] == s["representative_variant"] and
          r["rsid"] == q["representative_reported_rsids"] == s["representative_reported_rsids"] and
          r["classification"] == q["five_class_primary"] and
          b == float(q[trait + "_beta"]) == float(s[trait + "_beta"]) and
          se == float(q[trait + "_se"]) == float(s[trait + "_se"]) and
          r["ref"] == q["ref"] and r["effect_allele_ALT"] == q["alt"] == r["representative_variant"].split("_")[-1])
    check("F3 display-only 95% CI " + sid + "/" + trait,
          math.isclose(float(r["ci95_lower"]), b - 1.96 * se, abs_tol=1e-14) and
          math.isclose(float(r["ci95_upper"]), b + 1.96 * se, abs_tol=1e-14))
check("Figure 3 six stored equivalence entries", len(sensitivity) == 6)
for r in sensitivity:
    sid, bound = r["signal_unit_id"], float(r["equivalence_or_upper"])
    q = unit_index[sid]
    flag = "other_equivalent_or1_" + {1.05: "05", 1.10: "10", 1.15: "15"}[bound]
    check("F3 stored 90% CI and class " + sid + "/" + str(bound),
          float(r["beta"]) == float(q["other_beta"]) and
          float(r["stored_ci90_lower"]) == float(q["other_ci90_lower"]) and
          float(r["stored_ci90_upper"]) == float(q["other_ci90_upper"]) and
          r["stored_equivalent"] == q[flag] and
          r["existing_five_class"] == raw_sensitivity[(sid, bound)]["five_class"] and
          math.isclose(float(r["log_bound"]), math.log(bound), abs_tol=1e-15))
    check("F3 displayed bound reproduces stored flag " + sid + "/" + str(bound),
          (float(r["stored_ci90_lower"]) > -math.log(bound) and float(r["stored_ci90_upper"]) < math.log(bound)) == (r["stored_equivalent"] == "TRUE"))


class Drawing:
    """Use physical points for a stable publication-size layout."""
    def __init__(self, slug, height_mm, title, offset=0):
        self.slug, self.height_mm, self.height = slug, height_mm, height_mm * MM
        self.fig = plt.figure(figsize=(180 / 25.4, height_mm / 25.4), dpi=100)
        self.ax = self.fig.add_axes([0, 0, 1, 1], xlim=(0, WIDTH), ylim=(self.height + offset, offset))
        self.ax.set_axis_off()
        self.texts = []
        self.title = title

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
                     markeredgewidth=.9, color=style["color"], linestyle="none", zorder=zorder)[0]

    def panel(self, letter, title, y):
        self.text(16, y, letter, 11, bold=True)
        self.text(32, y, title, 9.5, bold=True)

    def finish(self):
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        extent = self.fig.bbox
        outliers = []
        for item in self.texts:
            bbox = item.get_window_extent(renderer)
            if not (bbox.x0 >= extent.x0 and bbox.x1 <= extent.x1 and bbox.y0 >= extent.y0 and bbox.y1 <= extent.y1):
                outliers.append(item.get_text())
        check(self.slug + " text within physical canvas", not outliers, outliers)
        check(self.slug + " all text at least 8 pt", min(t.get_fontsize() for t in self.texts) >= 8)
        visible_text = [t.get_text() for t in self.texts]
        check(self.slug + " removed explanatory text absent", not any(r["text"] in visible_text for r in REMOVED_TEXT))
        metadata = {"Title": self.title, "Creator": "Matplotlib display-only redraw; frozen T06/T07 outputs"}
        self.fig.savefig(FIG / (self.slug + ".svg"), format="svg", metadata=metadata)
        self.fig.savefig(FIG / (self.slug + ".pdf"), format="pdf", metadata=metadata)
        png_info = {}
        for dpi, suffix in ((600, ""), (180, "_preview")):
            path = FIG / f"{self.slug}{suffix}.png"
            self.fig.savefig(path, format="png", dpi=dpi, metadata={"Title": self.title, "Software": "Matplotlib"})
            with Image.open(path) as image:
                png_info[str(dpi)] = {"path": rel(path), "pixels": list(image.size), "dpi": list(image.info["dpi"])}
                check(self.slug + f" PNG {dpi} dimensions", image.size == (int(180 / 25.4 * dpi), int(self.height_mm / 25.4 * dpi)), image.size)
                check(self.slug + f" PNG {dpi} dpi metadata", all(abs(v - dpi) < .02 for v in image.info["dpi"]))
        svg_root = ET.parse(FIG / (self.slug + ".svg")).getroot()
        tags = [e.tag.rsplit("}", 1)[-1] for e in svg_root.iter()]
        check(self.slug + " editable SVG text retained", tags.count("text") == len(self.texts), tags.count("text"))
        check(self.slug + " no raster images in SVG", "image" not in tags)
        if self.slug in ("figure1", "figure1_panelB"):
            svg_cells = {e.attrib["id"]: e for e in svg_root.iter() if e.attrib.get("id", "").startswith("f1_cell_")}
            check(self.slug + " actual SVG retains all 80 cell groups", len(svg_cells) == 80)
            for row in geometry:
                group = svg_cells["f1_cell_" + row["cell_id"].replace(":", "_")]
                descendants = [e.tag.rsplit("}", 1)[-1] for e in group.iter()]
                check(self.slug + " SVG cell glyph " + row["cell_id"],
                      ("path" in descendants and "use" not in descendants) if row["comparison_count"] else "use" in descendants)
        fonts = subprocess.run(["pdffonts", str(FIG / (self.slug + ".pdf"))], text=True, capture_output=True, check=True).stdout
        check(self.slug + " PDF fonts embedded", len(fonts.splitlines()) > 2 and all("yes" in line for line in fonts.splitlines()[2:]), fonts)
        check(self.slug + " no raster images in PDF", b"/Subtype /Image" not in (FIG / (self.slug + ".pdf")).read_bytes())
        FIGURES.append({"figure": self.slug, "width_mm": 180, "height_mm": self.height_mm,
                        "minimum_font_pt": min(t.get_fontsize() for t in self.texts), "svg_text_count": tags.count("text"),
                        "png": png_info, "pdf_fonts": fonts})
        plt.close(self.fig)


def legend(d, y, classes, xs, chord=False):
    for classification, x in zip(classes, xs):
        style = STYLES[classification]
        if chord:
            d.line(x, y, x + 17, y, style["color"], 1, style["dash"])
            d.symbol(x + 8.5, y, style, 2)
            tx = x + 23
        else:
            d.symbol(x, y, style, 2.6)
            tx = x + 7
        d.text(tx, y + 3, style["label"])


REMOVED_TEXT.extend([
    {"figure": "Figure 1A", "text": "Finnish-enriched European ancestry · GRCh38 · summary statistics",
     "caption_location": "Figure 1A caption: source population, genome build and summary-statistic scope"},
    {"figure": "Figure 1B", "text": "Five classes per region pair; common radial scale: 0–4 comparisons.",
     "caption_location": "Figure 1B caption: five grouped class positions and shared integer radial scale"},
    {"figure": "Figure 1B", "text": "Counts describe tested CS pairs, not people or genome-wide sharing rates.",
     "caption_location": "Figure 1B caption: counting unit and interpretation boundary"},
    {"figure": "Figure 1B", "text": "Five positions per sector encode class in legend order; radial length = comparisons.",
     "caption_location": "Figure 1B caption: clockwise class order; radial length rather than area encodes counts"},
    {"figure": "Figure 1B", "text": "Zero counts have a baseline symbol. Angles do not encode genomic distance or causality.",
     "caption_location": "Figure 1B caption: 47 zero-height baseline symbols and angular-layout interpretation"},
    {"figure": "Figure 1B", "text": "RP order is organizational. * RP08/RP09 share a Graves region. The chrX pair is retained.",
     "caption_location": "Figure 1B caption: organizational order, retained shared-region asterisks and chrX"},
    {"figure": "Figure 2", "text": "Neutral circles = Graves; neutral squares = AIHT. n = CS size; min r² and purity are from the source CS.",
     "caption_location": "Figure 2 caption: node symbols, CS size, min r² and source low-purity flag"},
    {"figure": "Figure 2", "text": "Each equal-width chord is a tested pair; length/crossing has no quantitative or causal meaning. GRCh38.",
     "caption_location": "Figure 2 caption: tested-pair chord semantics, equal widths, noncausal layout and genome build"},
    {"figure": "Figure 3", "text": "Shaded bounds: ±ln(1.05), ±ln(1.10), ±ln(1.15). * Primary bound. Graves-pred. = Graves-predominant.",
     "caption_location": "Figure 3B caption: equivalence bounds, retained primary-margin asterisk and class abbreviation"},
    {"figure": "Figure 3", "text": "Class labels are from the existing T07 analysis; the plots display marginal associations.",
     "caption_location": "Figure 3 caption: original classifications and marginal-association scope"},
])


def xy(cx, cy, radius, angle):
    return cx + radius * math.cos(angle), cy - radius * math.sin(angle)


# The approved geometry is preserved exactly except for a global 40-pt upward
# translation. The independent panel export uses the same display coordinates.
ZERO_RADIUS, OUTER_RADIUS, COUNT_MAX = 85, 143, 4
COUNT_STEP, BAR_WIDTH_DEG, Y_TRANSLATION = 14.5, 2, -40
for classification, hatch, fill, offset in zip(CLASS_ORDER, [None, None, "///", None, "..."],
                                                [BLUE, ORANGE, "white", GRAY, "white"], [8, 4, 0, -4, -8]):
    STYLES[classification].update(hatch=hatch, bar_fill=fill, offset_degrees=offset)
approved_cells = read_tsv(APPROVED / "source_data/class_counts_80.tsv")
approved_geometry = read_tsv(APPROVED / "source_data/figure1_B_five_class_grouped_bars_geometry_80.tsv")
cell_rows, geometry = [], []
for cell, original in zip(approved_cells, approved_geometry):
    count = int(cell["comparison_count"])
    members = [r for r in f1 if r["pair_id"] == cell["pair_id"] and r["classification"] == cell["classification"]]
    check("F1 approved cell source equality " + cell["cell_id"],
          count == len(members) == int(next(r for r in regions if r["pair_id"] == cell["pair_id"])[cell["classification"]])
          and cell["comparison_ids"] == ";".join(f'{r["pair_id"]}:G{r["graves_cs"]}:A{r["aiht_cs"]}' for r in members))
    row = dict(cell, comparison_count=count)
    cell_rows.append(row)
    displayed = dict(original)
    displayed["comparison_count"] = count
    displayed["y_pt"] = float(original["y_pt"]) + Y_TRANSLATION
    displayed["x_pt"] = float(original["x_pt"])
    displayed["global_y_translation_pt"] = Y_TRANSLATION
    geometry.append(displayed)
check("F1 all 80 cells retained", len(cell_rows) == len({r["cell_id"] for r in cell_rows}) == 80)
check("F1 exact 33 bars and 47 zeros", Counter(r["comparison_count"] > 0 for r in cell_rows) == {True: 33, False: 47})
check("F1 five class sums", [sum(r["comparison_count"] for r in cell_rows if r["classification"] == c) for c in CLASS_ORDER] == [8, 11, 1, 17, 5])
check("F1 all comparison identities map once", sum(r["comparison_count"] for r in cell_rows) == 42 and
      len({v for r in cell_rows for v in r["comparison_ids"].split(";") if v}) == 42)
check("F1 0–4 count scale retained", max(r["comparison_count"] for r in cell_rows) == COUNT_MAX == 4)
check("F1 RP16 suggestive H4 remains one", [(r["pair_id"], r["comparison_count"]) for r in cell_rows if r["classification"] == CLASS_ORDER[2] and r["comparison_count"]] == [("RP16", 1)])
write_tsv("main_figure1_class_counts_80.tsv", cell_rows)
write_tsv("main_figure1_grouped_bars_geometry_80.tsv", geometry)
write_tsv("main_figure1_comparison_sources_42.tsv", f1)


def draw_panel_a(d):
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


def draw_grouped_bars(d):
    d.panel("B", "Five grouped class-count bars in each non-HLA region pair", 217)
    cx, cy = WIDTH / 2, 450 + Y_TRANSLATION
    artists, tick_texts, region_texts = [], [], []
    for tick in range(COUNT_MAX + 1):
        d.circle(cx, cy, ZERO_RADIUS + COUNT_STEP * tick,
                 "#A7B5BF" if tick == 0 else "#D9E2E7", .7 if tick == 0 else .5, 0)
    for i, region in enumerate(regions):
        pid = region["pair_id"]
        angle = math.pi / 2 - (i + .5) * 2 * math.pi / 16
        boundary = math.pi / 2 - i * 2 * math.pi / 16
        d.line(*xy(cx, cy, ZERO_RADIUS, boundary), *xy(cx, cy, OUTER_RADIUS, boundary), "#E2E7EB", .45)
        label_radius = 176 if pid == "RP06" else 164
        lx, ly = xy(cx, cy, label_radius, angle)
        star = "*" if pid in ("RP08", "RP09") else ""
        region_texts.append(d.text(lx, ly, f'{pid}{star} · chr{region["chromosome"]}', bold=True, ha="center"))
        total = int(region["n_signal_pair_tests"])
        region_texts.append(d.text(lx, ly + 11, str(total) + (" comparison" if total == 1 else " comparisons"), color=MUTED, ha="center"))
    for row, original in zip(geometry, approved_geometry):
        style = STYLES[row["classification"]]
        count, radius, angle = row["comparison_count"], float(row["radius_pt"]), float(row["displayed_angle_degrees"])
        if count:
            artist = Wedge((cx, cy), radius, -angle - BAR_WIDTH_DEG / 2, -angle + BAR_WIDTH_DEG / 2,
                           width=radius - ZERO_RADIUS, facecolor=style["bar_fill"], edgecolor=style["color"],
                           hatch=style["hatch"], linewidth=.65, zorder=2)
            d.ax.add_patch(artist)
            check(d.slug + " actual bar geometry " + row["cell_id"],
                  artist.center == (cx, cy) and artist.r == float(original["radius_pt"]) and
                  artist.width == COUNT_STEP * count and
                  artist.theta1 == -float(original["displayed_angle_degrees"]) - float(original["bar_angular_width_degrees"]) / 2 and
                  artist.theta2 == -float(original["displayed_angle_degrees"]) + float(original["bar_angular_width_degrees"]) / 2 and
                  artist.get_hatch() == style["hatch"])
        else:
            artist = d.symbol(row["x_pt"], row["y_pt"], style, 1.6, 4)
            artist.set_markeredgewidth(.8)
            check(d.slug + " actual zero glyph " + row["cell_id"],
                  artist.get_xdata()[0] == row["x_pt"] and artist.get_ydata()[0] == row["y_pt"] and
                  artist.get_marker() == original["marker"] and artist.get_markersize() == 3.2 and
                  artist.get_markerfacecolor() == original["marker_fill"])
        artist.set_gid("f1_cell_" + row["cell_id"].replace(":", "_"))
        artists.append(artist)
        check(d.slug + " translated radius " + row["cell_id"],
              math.isclose(math.hypot(row["x_pt"] - cx, row["y_pt"] - cy), ZERO_RADIUS + COUNT_STEP * count, abs_tol=1e-10))
    d.line(cx, cy - ZERO_RADIUS, cx, cy - OUTER_RADIUS, "#7C929F", .5, zorder=3)
    for tick in range(COUNT_MAX + 1):
        yy = cy - (ZERO_RADIUS + COUNT_STEP * tick)
        d.line(cx - 2, yy, cx + 2, yy, "#7C929F", .6, zorder=3)
        t = d.text(cx, yy + 2.6, str(tick), ha="center")
        t.set_bbox(dict(facecolor="white", edgecolor="none", pad=.1))
        tick_texts.append(t)
    d.text(cx, cy - 9, "42 comparisons", 10, bold=True, ha="center")
    d.text(cx, cy + 7, "16 region pairs", 8.5, ha="center")
    d.text(cx, cy + 22, "5 classes", ha="center")
    total_text = d.text(16, 633 + Y_TRANSLATION, "42 comparisons: 8 strong H4 · 11 strong H3 · 1 suggestive H4 · 17 low purity · 5 low overlap.")
    d.line(16, 639 + Y_TRANSLATION, 494, 639 + Y_TRANSLATION)
    for classification, x in zip(CLASS_ORDER, (18, 112, 206, 311, 398)):
        style = STYLES[classification]
        d.ax.add_patch(Rectangle((x + 1, 646 + Y_TRANSLATION), 11, 7, facecolor=style["bar_fill"],
                                edgecolor=style["color"], hatch=style["hatch"], linewidth=.65))
        d.symbol(x + 18, 650 + Y_TRANSLATION, style, 1.6).set_markeredgewidth(.8)
        d.text(x + 25, 653 + Y_TRANSLATION, style["label"])
    d.fig.canvas.draw()
    renderer = d.fig.canvas.get_renderer()
    tick_collisions = [t.get_text() for t in tick_texts if any(t.get_window_extent(renderer).padded(.1 * d.fig.dpi / 72).overlaps(a.get_window_extent(renderer)) for a in artists)]
    check(d.slug + " scale tick labels clear of all actual glyphs", not tick_collisions, tick_collisions)
    region_collisions = [(t.get_text(), row["cell_id"]) for t in region_texts for a, row in zip(artists, geometry)
                         if t.get_window_extent(renderer).padded(.3 * d.fig.dpi / 72).overlaps(a.get_window_extent(renderer))]
    check(d.slug + " region labels clear of all actual glyphs", not region_collisions, region_collisions)
    check(d.slug + " totals clear of region labels", not any(total_text.get_window_extent(renderer).overlaps(t.get_window_extent(renderer)) for t in region_texts))
    check(d.slug + " 80 actual data glyphs include 33 Wedges and 47 zero markers", len(artists) == 80 and
          sum(isinstance(a, Wedge) for a in artists) == 33 and sum(isinstance(a, Line2D) for a in artists) == 47)
    check(d.slug + " RP08/RP09 caption markers retained", all(any(t.get_text().startswith(pid + "*") for t in region_texts) for pid in ("RP08", "RP09")))


d = Drawing("figure1", 220, "Scope and grouped counts of five classes of non-HLA comparisons")
draw_panel_a(d)
draw_grouped_bars(d)
d.finish()
d = Drawing("figure1_panelB", 150, "Grouped counts of five classes of non-HLA comparisons", offset=200)
draw_grouped_bars(d)
d.finish()

def cubic_point(points, t):
    weights = ((1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3)
    return tuple(sum(w * p[c] for w, p in zip(weights, points)) for c in (0, 1))


# Figure 2: six complete chord panels with full source-region node descriptors.
d = Drawing("figure2", 235, "Complete component correspondence at six non-HLA regions")
pids = ["RP01", "RP02", "RP06", "RP07", "RP10", "RP12"]
names = ["PTPN22", "RTKN2–ARID5B", "CTLA4–ICOS", "AC007563.2–DIRC3", "ZAP70–TMEM131–COX5B", "LPP"]
node_map, edge_map = [], []
for i, (pid, name) in enumerate(zip(pids, names)):
    top, gcx, radius = 14 + i * 104, 63, 33
    gcy = top + 51
    panel_nodes = sorted([r for r in nodes if r["pair_id"] == pid], key=lambda r: (r["side"] != "graves", int(r["cs"])))
    gs = [r for r in panel_nodes if r["side"] == "graves"]
    ais = [r for r in panel_nodes if r["side"] == "aiht"]
    ga = [180] if len(gs) == 1 else [150, 210]
    aa = {4: [72, 24, -24, -72], 3: [60, 0, -60], 2: [40, -40]}[len(ais)]
    pn = []
    for row, degrees in zip(panel_nodes, ga + aa):
        theta = math.radians(degrees)
        n = dict(row, node_code=("G" if row["side"] == "graves" else "A") + row["cs"],
                 angle_radian=theta, x_pt=gcx + radius * math.cos(theta), y_pt=gcy - radius * math.sin(theta),
                 marker="o" if row["side"] == "graves" else "s")
        pn.append(n)
        node_map.append(n)
    d.circle(gcx, gcy, radius, "#DCE1E7")
    panel_edges = sorted([r for r in f2 if r["pair_id"] == pid],
                         key=lambda r: (not r["classification"].startswith("unresolved"), int(r["graves_cs"]), int(r["aiht_cs"])))
    for e in panel_edges:
        g = [n for n in pn if n["side"] == "graves" and n["cs"] == e["graves_cs"]]
        a = [n for n in pn if n["side"] == "aiht" and n["cs"] == e["aiht_cs"]]
        check("F2 one source node per edge " + "/".join(pair_key(e)), len(g) == len(a) == 1)
        g, a = g[0], a[0]
        points = [(g["x_pt"], g["y_pt"]), (gcx - 8, gcy + (g["y_pt"] - gcy) * .75),
                  (gcx + 8, gcy + (a["y_pt"] - gcy) * .75), (a["x_pt"], a["y_pt"])]
        style = STYLES[e["classification"]]
        patch = PathPatch(MplPath(points, [MplPath.MOVETO] + [MplPath.CURVE4] * 3),
                          fill=False, edgecolor=style["color"], linewidth=1, linestyle=style["dash"], zorder=2)
        d.ax.add_patch(patch)
        t = .73 if e["graves_cs"] == "2" else .51
        gx, gy = cubic_point(points, t)
        d.symbol(gx, gy, style, 2.1)
        edge_map.append(dict(e, start_node=g["node_code"], end_node=a["node_code"],
                             x0=points[0][0], y0=points[0][1], x1=points[1][0], y1=points[1][1],
                             x2=points[2][0], y2=points[2][1], x3=points[3][0], y3=points[3][1],
                             glyph_x=gx, glyph_y=gy, glyph_t=t, edge_width_pt=1,
                             edge_dash=str(style["dash"]), marker=style["marker"], color=style["color"]))
    for n in pn:
        color = "#9AA2AA" if n["low_purity"] == "TRUE" else "#26313B"
        d.symbol(n["x_pt"], n["y_pt"], dict(marker=n["marker"], color=color,
                                               fill=color if n["side"] == "graves" else "white"), 3.5, 5)
        theta = n["angle_radian"]
        lx, ly = gcx + (radius + 11) * math.cos(theta), gcy - (radius + 11) * math.sin(theta) + 2.5
        d.text(lx, ly, n["node_code"], color=color, bold=True, ha="center")
    d.text(128, top + 7, f"{chr(65 + i)}  {pid} · {name}", 9.5, bold=True)
    ar = raw_regions[pid]
    d.text(128, top + 20, "Graves " + ar["graves_region"], color=MUTED)
    d.text(128, top + 30, "AIHT   " + ar["aiht_region"], color=MUTED)
    for j, n in enumerate(pn):
        label = (f'{n["node_code"]}  {n["lead_variant"]}  n={n["cs_size"]}  '
                 f'min r²={float(n["cs_min_r2"]):.3f} ' + ("low" if n["low_purity"] == "TRUE" else "pass"))
        d.text(128, top + 44 + j * 10, label, color=MUTED if n["low_purity"] == "TRUE" else INK)
d.line(16, 640, 494, 640)
legend(d, 653, [CLASS_ORDER[i] for i in (0, 1, 3, 4)], (18, 139, 262, 382), chord=True)
check("F2 all 24 display edges match original set", len(edge_map) == 24 and {pair_key(r) for r in edge_map} == {pair_key(r) for r in f2})
check("F2 all 27 display nodes match original set", len(node_map) == 27 and {(r["pair_id"], r["node_key"]) for r in node_map} == {(r["pair_id"], r["node_key"]) for r in nodes})
check("F2 edge widths constant", {r["edge_width_pt"] for r in edge_map} == {1})
write_tsv("main_figure2_edges.tsv", edge_map)
write_tsv("main_figure2_nodes.tsv", node_map)
d.finish()


# Figure 3: paired marginal estimates and unchanged equivalence sensitivity.
d = Drawing("figure3", 195, "Marginal allelic associations and stored equivalence-bound sensitivity")
d.panel("A", "Marginal allelic associations at all ten prioritized variants", 23)
for x, color, marker, label in ((252, BLUE, "o", "Graves"), (353, ORANGE, "s", "AIHT")):
    d.symbol(x, 40, dict(marker=marker, color=color, fill=color), 2.8)
    d.text(x + 8, 43, label)
left, right, lo, hi, top, bottom = 232, 490, -.75, .60, 61, 338
xf = lambda beta: left + (beta - lo) / (hi - lo) * (right - left)
for tick in (-.6, -.3, 0, .3, .6):
    x = xf(tick)
    d.line(x, top - 8, x, bottom + 6, "#ADB6BF" if tick == 0 else "#E6E9ED",
           .7 if tick == 0 else .4, (0, (2, 2)) if tick == 0 else "solid")
    d.text(x, bottom + 20, f"{tick:.1f}".replace("-", "−"), ha="center")
f3_map = []
gene_names = {"GSU052": "CTLA4–ICOS", "GSU042": "TMEM131–ZAP70", "GSU056": "AC007563.2/DIRC3", "GSU254": "LINC01229/MAFTRR"}
for i, r in enumerate(priority):
    y, sid = top + i * 29, r["signal_unit_id"]
    if i % 2 == 0:
        d.rect(16, y - 10, 208, 28, "#F7F8FA", "none")
    gene = gene_names.get(sid, r["representative_nearest_genes"])
    d.text(19, y - 1, sid + " · " + r["representative_reported_rsids"] + " · " + gene)
    d.text(19, y + 11, r["representative_variant"] + " · ALT=" + r["alt"], color=MUTED)
    for trait, offset, color, marker in (("graves", -3, BLUE, "o"), ("aiht", 5, ORANGE, "s")):
        row = next(q for q in f3 if q["signal_unit_id"] == sid and q["trait"] == trait)
        b, lower, upper = float(row["beta"]), float(row["ci95_lower"]), float(row["ci95_upper"])
        check("F3 untruncated 95% CI " + sid + "/" + trait, lo < lower <= b <= upper < hi)
        x1, x2, yy = xf(lower), xf(upper), y + offset
        radius = 2.2 if trait == "graves" else 1.2
        check("F3 CI wider than marker " + sid + "/" + trait, x2 - x1 > radius * 2)
        d.line(x1, yy, x2, yy, color, 1, zorder=3)
        d.line(x1, yy - 2, x1, yy + 2, color, .7)
        d.line(x2, yy - 2, x2, yy + 2, color, .7)
        d.symbol(xf(b), yy, dict(marker=marker, color=color, fill=color), radius)
        f3_map.append(dict(row, displayed_gene_label=gene, x_pt=xf(b), y_pt=yy, ci95_lower_x_pt=x1,
                           ci95_upper_x_pt=x2, marker=marker, color=color))
d.text((left + right) / 2, 376, "Marginal beta (log odds ratio) per ALT allele; 95% CI", ha="center")
d.panel("B", "AIHT equivalence at NRG1 and TSHR: stored 90% CIs", 406)
bl, br = 104, 289
bx = lambda beta: bl + (beta + .16) / .32 * (br - bl)
for bound, fill in ((1.15, "#F0F2F5"), (1.10, "#E1E6EB"), (1.05, "#CAD3DC")):
    value = math.log(bound)
    d.rect(bx(-value), 430, bx(value) - bx(-value), 72, fill, "none")
    for v in (-value, value):
        d.line(bx(v), 430, bx(v), 502, "#9DA9B5", .5, (0, (2, 2)))
d.line(bx(0), 430, bx(0), 502, MUTED, .7)
d.text(402, 426, "Existing class at OR upper bound", bold=True, ha="center")
for x, text in zip((330, 399, 468), ("1.05", "1.10*", "1.15")):
    d.text(x, 443, text, bold=True, ha="center")
sensitivity_map = []
for y, sid, gene in ((461, "GSU138", "NRG1"), (492, "GSU231", "TSHR")):
    row = next(r for r in sensitivity if r["signal_unit_id"] == sid)
    d.text(16, y - 2, gene, 8.5, bold=True)
    d.text(16, y + 9, sid, color=MUTED)
    x1, x2, yy = bx(float(row["stored_ci90_lower"])), bx(float(row["stored_ci90_upper"])), y - 1
    d.line(x1, yy, x2, yy, ORANGE, 1)
    d.line(x1, yy - 2, x1, yy + 2, ORANGE, .7)
    d.line(x2, yy - 2, x2, yy + 2, ORANGE, .7)
    d.symbol(bx(float(row["beta"])), yy, dict(marker="s", color=ORANGE, fill=ORANGE), 1.5)
    for x, bound in zip((330, 399, 468), (1.05, 1.10, 1.15)):
        r = next(q for q in sensitivity if q["signal_unit_id"] == sid and float(q["equivalence_or_upper"]) == bound)
        cls = r["existing_five_class"]
        label = "Unresolved" if cls == "unresolved" else "Graves-pred."
        d.text(x, y + 2, label, color=MUTED if cls == "unresolved" else INK, ha="center")
        sensitivity_map.append(dict(r, displayed_gene_label=gene, x_pt=bx(float(r["beta"])), y_pt=yy,
                                    ci90_lower_x_pt=x1, ci90_upper_x_pt=x2, displayed_class_label=label,
                                    bound_lower_x_pt=bx(-math.log(bound)), bound_upper_x_pt=bx(math.log(bound))))
for tick in (-.15, 0, .15):
    d.text(bx(tick), 519, f"{tick:.2f}".replace("-", "−"), ha="center")
d.text((bl + br) / 2, 536, "AIHT marginal beta (log odds ratio)", ha="center")
write_tsv("main_figure3_marginal_estimates.tsv", f3_map)
write_tsv("main_figure3_equivalence_sensitivity.tsv", sensitivity_map)
d.finish()

after = hashes(source_paths)
check("All authoritative tables, old figures and prior builders unchanged", before == after)
record = {
    "status": "PASS_DISPLAY_AND_EXPORT_CHECKS_PENDING_INDEPENDENT_FIGURE_REVIEW",
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "scope": "DISPLAY ONLY: approved grouped radial bars and prose removal; no analysis or model rerun",
    "command": "uv run --no-sync --python /usr/bin/python3 python " + rel(Path(__file__)),
    "runtime": {"python": sys.executable, "python_version": platform.python_version(), "matplotlib": matplotlib.__version__,
                "environment_changes": "None; existing system runtime used as previously approved"},
    "font": FONT, "figures": FIGURES, "checks": CHECKS, "source_hashes_sha256": after,
    "removed_text": REMOVED_TEXT,
    "size_changes_mm": {"figure1": [245, 220], "figure2": [245, 235], "figure3": [205, 195], "figure1_panelB": [167.38888888888889, 150]},
    "geometry": {"zero_radius_pt": ZERO_RADIUS, "outer_radius_pt": OUTER_RADIUS, "points_per_count": COUNT_STEP,
                 "bar_width_degrees": BAR_WIDTH_DEG, "class_offsets_degrees": [8, 4, 0, -4, -8],
                 "global_y_translation_pt": Y_TRANSLATION, "cells": 80, "positive_bars": 33, "zero_markers": 47},
    "panel_mapping": {"figure1": {"A": "Two source endpoints and non-HLA scope", "B": "Five grouped radial class counts per region pair; 80 cells"},
                      "figure2": dict(zip(list("ABCDEF"), [p + " " + name for p, name in zip(pids, names)])),
                      "figure3": {"A": "All ten candidates and twenty marginal estimates", "B": "NRG1/TSHR stored 90% CIs and six classifications"}},
    "limitations": ["Implementation/source and visual self-checks are not independent scientific verification or author acceptance.",
                    "Formal caption content is owned by the primary agent; removed explanatory strings are mapped to caption locations."]}
(REV / "main_build_validation.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
(REV / "main_removed_text.json").write_text(json.dumps({"removed_strings": REMOVED_TEXT,
    "retained_caption_markers": ["RP08*", "RP09*", "1.10*"]}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
(REV / "main_build_notes.md").write_text(
    "# Clean main-figure build\n\n"
    "The user-selected five-class grouped radial-bar design B was redrawn with the prior approved geometry. "
    "All 80 class-by-region cells remain: 33 positive bars and 47 class-specific zero markers. "
    "The original zero radius (85 pt), 14.5 pt per comparison, 2° bar widths and class angles are unchanged. "
    "The geometry is translated upward by 40 pt to remove unused space.\n\n"
    "Ten explanatory strings were removed from the artwork and mapped to formal-caption locations in main_removed_text.json. "
    "Titles, legends, axes, units, ticks, data/node labels and RP08*/RP09*/1.10* caption markers remain. "
    "All original scientific sources and prior exports/builders retain their SHA-256 hashes.\n\n"
    "Figures 1/2/3 are 180×220/235/195 mm, respectively. Standalone Figure 1B is 180×150 mm. "
    "All fonts remain at least 8 pt. Main PNGs are 600 dpi, previews 180 dpi; SVG text is editable and PDF is vector. "
    "No model, classification or statistical test was rerun; no dependencies were installed or changed.\n", encoding="utf-8")
print(json.dumps({"figures": FIGURES, "checks": len(CHECKS), "passed": sum(c["passed"] for c in CHECKS)}, indent=2))
