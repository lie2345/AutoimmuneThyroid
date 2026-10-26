#!/usr/bin/env python3
"""Display-only redraw of stored supplementary results; never fit a model.

Run from the project root with the existing Matplotlib installation:
uv run --no-sync --python /usr/bin/python3 python \
  outputs/sci_writing_20261017/revisions/matplotlib_t19_t20_v2/scripts/build_supplement_figures.py
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist
from xml.etree import ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parents[1]
FIG = OUT / "figures"
DATA = OUT / "source_data"
REVIEW = OUT / "reviews"
for directory in (FIG, DATA, REVIEW):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .6, "xtick.major.width": .6,
    "ytick.major.width": .6, "svg.fonttype": "none",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})
BLUE = "#286895"
INK = "#263746"
GREY = "#7C8791"
PALE = "#F3F6F8"
ORANGE = "#A85828"
Z95 = NormalDist().inv_cdf(.975)
CHECKS = []
INPUTS = {}
EXPORTS = []


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(name, condition, detail):
    CHECKS.append({"check": name, "pass": bool(condition), "detail": detail})
    if not condition:
        raise AssertionError(f"{name}: {detail}")


def read(relative):
    path = ROOT / relative
    INPUTS[relative] = {"sha256": sha(path), "bytes": path.stat().st_size}
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def register(relative):
    path = ROOT / relative
    INPUTS[relative] = {"sha256": sha(path), "bytes": path.stat().st_size}


def num(value):
    return float(value) if value not in ("NA", "", None) else np.nan


def write_data(df, name):
    path = DATA / name
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    return path


def export(fig, stem):
    for extension in ("png", "svg", "pdf"):
        path = FIG / f"{stem}.{extension}"
        fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=.12)
        EXPORTS.append(path)
    plt.close(fig)


def plain_axis(ax):
    ax.set_facecolor("white")
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(axis="x", color="#E5E9ED", lw=.6, zorder=0)
    ax.set_axisbelow(True)


def p_string(value):
    return "NA" if value in ("NA", "") else f"{float(value):.3g}"


def build_s1():
    features = read("outcome_data/t11b/p1a_feature_results.tsv")
    screens = read("outcome_data/t11b/p1a_screen_reclassification.tsv.gz")
    t13 = read("outcome_data/t13/t13_candidate_summary.tsv")
    register("documentation/analysis/T11B_P1A_WAVE2.md")
    register("documentation/analysis/T13_LOCUS_REGULATORY_CASCADE.md")
    key = ["locus_id", "modality", "level", "cell_type", "feature_id"]
    check("S1_feature_counts", len(features) == 26 and Counter(features.feature_result) == {
        "confirmed_both_same_direction": 17, "downgraded_or_mixed": 7, "conflicting": 2},
        dict(Counter(features.feature_result)))
    check("S1_unique_keys", not features.duplicated(key).any() and not screens.duplicated(key + ["trait_id"]).any(),
          "26 feature keys and 52 trait-feature keys are unique")
    check("S1_row_status_counts", len(screens) == 52 and Counter(screens.screen_reclassification) == {
        "confirmed": 34, "conflicting": 4, "downgraded": 11, "untestable": 3},
        dict(Counter(screens.screen_reclassification)))
    left = features.sort_values(key).reset_index(drop=True)
    right = t13.sort_values(key).reset_index(drop=True)
    check("S1_T11B_T13_identity", left[key + ["feature_result"]].equals(right[key + ["feature_result"]]),
          "Frozen feature identities and historical classifications agree exactly")
    check("S1_confirmed_locus_counts", Counter(features.loc[features.feature_result.eq("confirmed_both_same_direction"), "locus_id"].str[:4]) == {"L033": 2, "L080": 7, "L103": 8},
          "Historical confirmed features: L033=2, L080=7, L103=8; distinct from later T21 quality gates")
    features = features.merge(t13[key + ["candidate_id", "candidate_evidence_class", "gene_symbol", "n_chain_rows", "n_complete_cascade"]], on=key, validate="one_to_one")
    features.insert(0, "display_row", range(1, 27))
    screens = screens.merge(features[key + ["display_row", "candidate_id"]], on=key, validate="many_to_one")
    write_data(features, "supp_S1_features.tsv")
    write_data(screens, "supp_S1_trait_feature_quality.tsv")

    fig = plt.figure(figsize=(14.8, 11.3))
    ax = fig.add_axes([.025, .13, .95, .76])
    ax.set_xlim(0, 14.2)
    ax.set_ylim(26.9, -.65)
    ax.axis("off")
    xedge = [0, 5.2, 7.0, 8.35, 10.15, 11.5, 12.85, 14.2]
    colors = {"confirmed": "#D7EAF4", "downgraded": "#E9ECEF", "conflicting": "#F6E2D4", "untestable": "#F5F5F5"}
    outcome = {"shared_strong_H4": "Confirmed H4", "distinct_signals_H3": "Conflicting H3", "unresolved_insufficient_overlap": "Low overlap", "unresolved_low_purity": "Low purity", "NA": "No eligible pair"}
    summary = {"confirmed_both_same_direction": "Confirmed\nboth / same dir.", "downgraded_or_mixed": "Downgraded\nor mixed", "conflicting": "Conflicting"}
    contexts = {"caQTL-only": "caQTL only", "eQTL-only": "eQTL only", "evidence_insufficient": "Insufficient\nchain evidence"}
    for row_index, row in features.iterrows():
        y = row_index
        if row_index % 2 == 0:
            ax.add_patch(patches.Rectangle((0, y-.48), 14.2, .96, facecolor=PALE, edgecolor="none"))
        locus = row.locus_id.split("_")[0]
        feature = row.feature_id
        if row.gene_symbol != "NA":
            feature += f" ({row.gene_symbol})"
        ax.text(.06, y, f"{locus} | {row.modality} | {row.cell_type.replace('_', ' ')} ({row.level}) | {feature}", va="center", fontsize=7.7, color=INK)
        subset = screens.loc[screens.display_row.eq(row.display_row)]
        for j, trait in enumerate(("E4_GRAVES_STRICT", "E4_HYTHY_AI_STRICT")):
            r = subset.loc[subset.trait_id.eq(trait)].iloc[0]
            x = xedge[1 + j*2]
            ax.add_patch(patches.Rectangle((x+.035, y-.41), 1.73, .82, facecolor=colors[r.screen_reclassification], edgecolor="white", lw=.6))
            ax.text(x+.9, y, outcome[r.direct_pair_classification], ha="center", va="center", fontsize=7.5, color=INK)
            a, b = num(r["direct_PP.H4"]), num(r["sensitivity_PP.H4"])
            text_value = "NA / NA" if not np.isfinite(a) else f"{a:.3f} / {b:.3f}"
            ax.text(xedge[2 + j*2]+.675, y, text_value, ha="center", va="center", fontsize=7.5)
        feature_color = colors[{"confirmed_both_same_direction": "confirmed", "downgraded_or_mixed": "downgraded", "conflicting": "conflicting"}[row.feature_result]]
        ax.add_patch(patches.Rectangle((xedge[5]+.035, y-.41), 1.28, .82, facecolor=feature_color, edgecolor="white"))
        ax.text(xedge[5]+.675, y, summary[row.feature_result], ha="center", va="center", fontsize=7.0)
        ax.text(xedge[6]+.675, y, contexts[row.candidate_evidence_class], ha="center", va="center", fontsize=7.0, color=GREY)
    headers = ["Preselected cell–feature (26)", "Graves: direct outcome", "H4 main / sens.", "AIHT: direct outcome", "H4 main / sens.", "Feature summary", "T13 context"]
    for i, title in enumerate(headers):
        ax.text(xedge[i]+(.06 if i == 0 else (xedge[i+1]-xedge[i])/2), -.9, title,
                ha="left" if i == 0 else "center", va="bottom", fontsize=8, fontweight="bold", color=INK)
    for x in xedge[1:-1]:
        ax.axvline(x, color="white", lw=1.1)
    fig.text(.025, .97, "Supplementary Figure S1 | Quality-aware matrix of the frozen QTL features", fontsize=13, weight="bold", color=INK)
    fig.text(.025, .94, "Historical T11B summary: 17 confirmed in both traits with the same direction; 7 downgraded/mixed; 2 conflicting", fontsize=9.5, color=INK)
    fig.text(.025, .915, "52 trait–feature rows: 34 confirmed, 11 downgraded, 4 conflicting and 3 untestable. All rows are retained.", fontsize=9, color=GREY)
    fig.text(.025, .092, "H4: stored primary / sensitivity-prior posterior (rounded to 3 decimals). Quality failures remain unresolved even when H4 is high.", fontsize=8.5, color=INK)
    fig.text(.025, .068, "Confirmed denotes the historical T11B rule; it is not later T21 admission, causal validation, or a complete regulatory cascade.", fontsize=8.5, color=INK)
    fig.text(.025, .044, "T13: 20 caQTL-only, 2 eQTL-only and 4 insufficient-chain features; 3 candidate chains, 0 complete cascades. Healthy-donor PBMC context.", fontsize=8.5, color=INK)
    export(fig, "figureS1")


def build_s2():
    forest = read("outcome_data/t18/limited_synthesis_final_20260904/forest_data.tsv")
    targets = read("outcome_data/t18/limited_synthesis_final_20260904/priority_variant_validation.tsv")
    original = read("outcome_data/external_validation/run_20260904/priority_variant_validation.tsv")
    register("codes/64_plot_t18_forest.R")
    register("outcome_data/t18/limited_synthesis_final_20260904/T18_LIMITED_SYNTHESIS.md")
    check("S2_20_planned_tests", len(targets) == 20 and Counter(targets.classification) == {
        "same_direction_association_support": 3, "inconclusive_low_precision_OR1_10_scenario": 7, "untestable_missing_aligned_target": 10}, dict(Counter(targets.classification)))
    check("S2_final_vs_external_source", targets.equals(original), "Final synthesis preserved the original 20-test table exactly")
    check("S2_plot_counts", len(forest) == 40 and forest.available.eq("True").sum() == 30 and forest.loc[forest.cohort.eq("BBJ"), "available"].eq("False").sum() == 10,
          "40 cohort rows: 20 FinnGen, 10 observed BBJ, 10 explicit missing BBJ")
    maximum_error = 0.
    for _, row in forest.iterrows():
        t = targets.loc[targets.signal_unit_id.eq(row.signal_unit_id) & targets.trait.eq(row.trait)].iloc[0]
        prefix = "discovery" if row.cohort == "FinnGen" else "external"
        check(f"S2_source_identity_{row.signal_unit_id}_{row.trait}_{row.cohort}", row.beta == t[f"{prefix}_beta"] and row.se == t[f"{prefix}_se"] and row.source_p == t[f"{prefix}_{'p' if prefix == 'discovery' else 'source_p'}"], "beta, SE and source-P strings are preserved")
        if row.available == "True":
            b, s = num(row.beta), num(row.se)
            expected = [math.exp(b), math.exp(b-Z95*s), math.exp(b+Z95*s)]
            actual = [num(row.odds_ratio), num(row.lower95), num(row.upper95)]
            maximum_error = max(maximum_error, max(abs(a-b) for a,b in zip(expected, actual)))
            check(f"S2_interval_order_{row.signal_unit_id}_{row.trait}_{row.cohort}", actual[1] <= actual[0] <= actual[2], "Observed OR lies inside its stored interval")
        else:
            check(f"S2_missing_{row.signal_unit_id}_{row.trait}", all(row[k] == "NA" for k in ("beta", "se", "odds_ratio", "lower95", "upper95", "source_p")), "Missing BBJ is not assigned a zero/null estimate")
    check("S2_CI_conversion", maximum_error < 1e-12, {"max_absolute_error": maximum_error, "z": Z95})
    write_data(forest, "supp_S2_forest_data.tsv")
    write_data(targets, "supp_S2_20_variant_checks.tsv")
    tags = {"same_direction_association_support": "Same-direction support", "inconclusive_low_precision_OR1_10_scenario": "Inconclusive; low precision*", "untestable_missing_aligned_target": "Untestable (BBJ missing)"}
    fig = plt.figure(figsize=(15.2, 13.5))
    fig.text(.025, .98, "Supplementary Figure S2 | Existing FinnGen–BBJ representative-variant checks", fontsize=13, weight="bold", color=INK)
    fig.text(.025, .958, "20 planned tests: 3 same-direction association supports, 7 inconclusive/low precision and 10 untestable", fontsize=10, color=INK)
    intervals = forest.loc[forest.available.eq("True")]
    lo = intervals.lower95.astype(float).min() * .93
    hi = intervals.upper95.astype(float).max() * 1.07
    for panel, trait in enumerate(("Graves", "AIHT")):
        bottom = .535 if panel == 0 else .115
        height = .335
        lab = fig.add_axes([.025, bottom, .32, height])
        ax = fig.add_axes([.365, bottom, .28, height])
        txt = fig.add_axes([.67, bottom, .305, height])
        for a in (lab, ax, txt):
            a.set_ylim(9.65, -.65)
        for a in (lab, txt):
            a.set_xlim(0, 1); a.axis("off")
        plain_axis(ax)
        ax.set_xscale("log"); ax.set_xlim(lo, hi)
        ticks = [.5, .6, .8, 1, 1.2, 1.5]
        ax.set_xticks(ticks, labels=[str(x) for x in ticks]); ax.minorticks_off()
        ax.axvline(1, ls="--", color=GREY, lw=.8)
        ax.set_yticks([])
        ax.set_xlabel("Odds ratio (log scale)")
        fig.text(.025, bottom+height+.042, f"{'AB'[panel]}  {trait}", fontsize=11, weight="bold", color=INK)
        lab.text(0, -.93, "Representative / nearest gene; effect / other allele", fontsize=8, weight="bold")
        txt.text(0, -.93, "BBJ source P / variant-level result", fontsize=8, weight="bold")
        subset = forest.loc[forest.trait.eq(trait)]
        for i, gsu in enumerate(sorted(subset.signal_unit_id.unique())):
            rows = subset.loc[subset.signal_unit_id.eq(gsu)]
            fg = rows.loc[rows.cohort.eq("FinnGen")].iloc[0]
            bbj = rows.loc[rows.cohort.eq("BBJ")].iloc[0]
            if i % 2 == 0:
                for a in (lab, ax, txt):
                    a.axhspan(i-.45, i+.45, color=PALE, zorder=0)
            gene = fg.nearest_genes if fg.nearest_genes != "NA" else "unannotated"
            lab.text(.005, i-.12, f"{gsu} | {gene}", fontsize=8.4, va="center", color=INK)
            lab.text(.005, i+.20, f"{fg.variant_id} | EA {fg.effect_allele} / OA {fg.other_allele}", fontsize=7.2, va="center", color=GREY)
            for _, r in rows.iterrows():
                if r.available == "True":
                    b, lower, upper = num(r.odds_ratio), num(r.lower95), num(r.upper95)
                    y = i + (-.13 if r.cohort == "FinnGen" else .13)
                    ax.errorbar(b, y, xerr=[[b-lower], [upper-b]], fmt="s" if r.cohort == "FinnGen" else "o", color=BLUE if r.cohort == "FinnGen" else INK, markersize=4.2, capsize=2, elinewidth=1.1, zorder=3)
            txt.text(.005, i-.12, f"P = {p_string(bbj.source_p)}", fontsize=8.3, va="center", color=INK)
            txt.text(.005, i+.20, tags[bbj.external_classification], fontsize=7.6, va="center", color=GREY)
        ax.plot([], [], "s", ms=4, color=BLUE, label="FinnGen")
        ax.plot([], [], "o", ms=4, color=INK, label="BBJ")
        ax.legend(frameon=False, ncol=2, fontsize=8, loc="lower right", bbox_to_anchor=(1, 1.06), borderaxespad=0)
    fig.text(.025, .066, "Intervals are stored ordinary normal-approximation 95% CIs; they are not SPA-calibrated. Association support uses source P ≤ 0.0025 (20 planned tests).", fontsize=8.2)
    fig.text(.025, .043, "*All 10 evaluable BBJ tests have <80% power in the fixed OR=1.10 precision scenario, including the 3 supported variants. Missing estimates are omitted explicitly.", fontsize=8.2)
    fig.text(.025, .020, "These are variant-level checks. Discovery signal classes are not externally revalidated; no pooled estimate or cross-disease difference test is shown.", fontsize=8.2)
    check("S2_all_intervals_fit", intervals.lower95.astype(float).min() > lo and intervals.upper95.astype(float).max() < hi, {"common_OR_axis": [lo, hi]})
    export(fig, "figureS2")


def build_s3():
    leads = read("outcome_data/ccgwas/summary/provisional_lead_loci.tsv")
    candidates = read("outcome_data/ccgwas/summary/scenario_candidate_variants.tsv.gz")
    summary = read("outcome_data/ccgwas/summary/scenario_summary.tsv")
    stability = read("outcome_data/ccgwas/summary/signal_stability.tsv")
    analysis = read("outcome_data/ccgwas/summary/analysis_summary.tsv")
    register("codes/15_summarize_ccgwas_path_b.R")
    scenarios = ["overlap_000pct", "overlap_050pct", "overlap_100pct"]
    display = candidates.merge(leads, left_on="SNP", right_on="lead_SNP", validate="many_to_one", suffixes=("", "_lead"))
    display["display_lower95"] = display.OLS_beta.astype(float) - 1.96*display.OLS_se.astype(float)
    display["display_upper95"] = display.OLS_beta.astype(float) + 1.96*display.OLS_se.astype(float)
    check("S3_lead_scenario_counts", len(leads) == 16 and len(display) == 48 and not display.duplicated(["locus_id", "scenario"]).any(), "16 provisional physical-region leads × 3 stored scenarios")
    check("S3_all_scenarios_present", all(set(x.scenario) == set(scenarios) for _,x in display.groupby("locus_id")), "All three scenarios occur for every lead")
    check("S3_stored_summary", summary.n_formal_ccgwas.tolist() == ["324", "328", "334"] and len(stability) == 324 and stability.formal_all_overlap_scenarios.eq("TRUE").all(), "Formal SNP totals 324 / 328 / 334; all 324 original SNPs remain formal")
    check("S3_robustness_counts", Counter(stability.final_sensitivity_status) == {"robust_candidate": 321, "unresolved_threshold_sensitive": 3}, dict(Counter(stability.final_sensitivity_status)))
    write_data(display, "supp_S3_48_lead_scenario_effects.tsv")
    write_data(summary, "supp_S3_scenario_summary.tsv")
    write_data(analysis, "supp_S3_analysis_summary.tsv")
    fig, axes = plt.subplots(4, 4, figsize=(13.2, 11.1))
    fig.subplots_adjust(left=.075, right=.98, top=.86, bottom=.13, wspace=.34, hspace=.58)
    for ax, (_, lead) in zip(axes.flat, leads.iterrows()):
        values = display.loc[display.locus_id.eq(lead.locus_id)].set_index("scenario").loc[scenarios]
        b = values.OLS_beta.astype(float).to_numpy()
        low = values.display_lower95.to_numpy(); high = values.display_upper95.to_numpy()
        ax.axhline(0, color=GREY, ls="--", lw=.7)
        ax.errorbar(range(3), b, yerr=[b-low, high-b], fmt="o-", color=BLUE, markersize=4.2, capsize=3, lw=1.0)
        bottom, top = min(0., low.min()), max(0., high.max())
        margin = (top-bottom)*.15
        ax.set_ylim(bottom-margin, top+margin)
        ax.set_xlim(-.35, 2.35)
        ax.set_xticks(range(3), ["0%", "50%", "100%"], fontsize=8)
        ax.tick_params(axis="y", labelsize=7.7)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        gene = lead.lead_nearest_genes_graves
        ax.set_title(f"{lead.locus_id.split('_')[0]} | {gene}\n{lead.lead_SNP}  (EA {lead.lead_EA})", fontsize=8.6, loc="left", color=INK, pad=7)
        ax.grid(axis="y", color="#E5E9ED", lw=.5)
    fig.text(.025, .97, "Supplementary Figure S3 | CC-GWAS lead stability under shared-control assumptions", fontsize=13, weight="bold", color=INK)
    fig.text(.025, .941, "16 provisional non-HLA physical regions; original stored OLS estimates and approximate 95% CIs (β ± 1.96 × SE)", fontsize=9.5)
    fig.text(.025, .915, "Formal SNP totals: 324 at 0% assumed overlap, 328 at 50%, 334 at 100%. All 324 primary SNPs remain formal in every scenario.", fontsize=9, color=GREY)
    fig.text(.015, .52, "CC-GWAS OLS beta (positive: Graves relative to AIHT)", rotation=90, va="center", fontsize=10, color=INK)
    fig.text(.52, .083, "Assumed shared controls as a fraction of the smaller (AIHT) control set", ha="center", fontsize=10)
    fig.text(.025, .050, "0%, 50%, 100% correspond to 0, 179,228 and 358,456 assumed shared controls. Per-panel y-axis ranges vary and all include zero.", fontsize=8.4)
    fig.text(.025, .027, "Summary-statistic sensitivity analysis; these are not direct case–case GWAS estimates, independent signals, or external replications.", fontsize=8.4)
    export(fig, "figureS3")


def build_s4():
    data = read("outputs/t19_t24_20261017/t20/lookup_results.tsv")
    endpoints = read("outputs/t19_t24_20261017/t20/endpoint_metadata.tsv")
    frozen = read("outputs/t19_t24_20261017/t20/frozen_panel.tsv")
    check("S4_exact_24_queries", len(data) == 24 and data.variant_id.nunique() == 12 and not data.duplicated(["endpoint", "variant_id"]).any(), "12 frozen variants × 2 endpoints; 24 unique queries")
    check("S4_match_QC_and_missing", data.match_status.eq("EXACT_GRCH38_REF_ALT").all() and data.qc_status.eq("PASS_GENERAL_QC").all() and data[["beta", "se", "source_p_string"]].ne("NA").all().all(), "24 exact matches; 24 general-QC passes; 0 missing/unresolved queries")
    check("S4_population_noncase_controls", endpoints.control_definition_verbatim.eq("Controls for this endpoint are individuals that are not cases.").all(), "Both endpoints use population noncase controls")
    check("S4_frozen_panel", set(data.variant_id) == set(frozen.variant_id) and set(data.panel_index) == set(frozen.panel_index), "All 12 prelisted variants retained")
    maximum_error = 0.
    for _, row in data.iterrows():
        b, s = num(row.beta), num(row.se)
        actual = [num(row.OR_ALT), num(row.OR_lower95_normal_approx), num(row.OR_upper95_normal_approx)]
        expected = [math.exp(b), math.exp(b-Z95*s), math.exp(b+Z95*s)]
        maximum_error = max(maximum_error, max(abs(a-b) for a,b in zip(actual, expected)))
        check(f"S4_source_strings_{row.panel_index}_{row.endpoint}", row.beta == row.beta_source_string and row.se == row.se_source_string and float(row.source_p_string) == float(row.p_source_numeric), "Stored beta, SE and source P are retained; no P is recalculated")
    check("S4_CI_conversion", maximum_error < 1e-12, {"max_absolute_error": maximum_error, "z": Z95})
    data["display_group"] = data.panel_index.astype(int).map(lambda x: "Non-HLA discovery-selected representative" if x <= 10 else "HLA regional lead")
    data["display_row"] = data.panel_index.astype(int).map(lambda x: x-1 if x <= 10 else x-.15)
    check("S4_endpoint_counts", Counter(data.endpoint) == {"GRAVES_OPHT": 12, "E4_GRAVES_OPHT_STRICT": 12}, dict(Counter(data.endpoint)))
    check("S4_source_P_threshold_counts", {ep: [sum(x.p_source_numeric.astype(float) < .05), sum(x.p_source_numeric.astype(float) < .05/24)] for ep,x in data.groupby("endpoint")} == {"GRAVES_OPHT": [10,8], "E4_GRAVES_OPHT_STRICT": [11,7]}, "Source P counts only; no threshold changes or new test")
    write_data(data, "supp_S4_24_ophthalmopathy_queries.tsv")
    write_data(endpoints, "supp_S4_endpoint_metadata.tsv")
    fig = plt.figure(figsize=(14.8, 9.6))
    fig.text(.025, .97, "Supplementary Figure S4 | Exploratory ophthalmopathy associations at selected variants", fontsize=13, weight="bold", color=INK)
    fig.text(.025, .94, "12 discovery-selected variants × 2 FinnGen R13 endpoints; all 24 exact-match queries are shown", fontsize=10)
    label = fig.add_axes([.025, .18, .30, .63])
    label.set_xlim(0, 1); label.set_ylim(12.65, -1.3); label.axis("off")
    base = data.loc[data.endpoint.eq("GRAVES_OPHT")].sort_values("panel_index", key=lambda x: x.astype(int))
    label.text(.0, -1.05, "Representative variant; effect / other allele", fontsize=8.3, weight="bold")
    label.text(.0, -.53, "Non-HLA representatives (n=10)", fontsize=8.2, weight="bold", color=BLUE)
    for _, r in base.iterrows():
        y = r.display_row
        group = "HLA regional lead" if int(r.panel_index) > 10 else (r.nearest_genes_discovery if r.nearest_genes_discovery != "NA" else "unannotated")
        label.text(0, y-.10, f"{r.rsids_source} | {group}", fontsize=8.4, va="center", color=INK)
        label.text(0, y+.22, f"{r.variant_id} | EA {r.effect_allele} / OA {r.other_allele}", fontsize=7.2, va="center", color=GREY)
    label.text(0, 10.20, "HLA regional leads (n=2)", fontsize=8.2, weight="bold", color=BLUE)
    lo = data.OR_lower95_normal_approx.astype(float).min()*.93
    hi = data.OR_upper95_normal_approx.astype(float).max()*1.07
    for panel, ep in enumerate(("GRAVES_OPHT", "E4_GRAVES_OPHT_STRICT")):
        left = .335 + panel*.33
        ax = fig.add_axes([left, .18, .23, .63])
        p_ax = fig.add_axes([left+.238, .18, .085, .63])
        ax.set_ylim(12.65, -1.3); p_ax.set_ylim(12.65, -1.3)
        p_ax.set_xlim(0,1); p_ax.axis("off")
        plain_axis(ax)
        ax.set_xscale("log"); ax.set_xlim(lo,hi)
        ticks=[.6, .8, 1, 1.2, 1.5, 2, 2.5]
        ticks=[x for x in ticks if lo < x < hi]
        ax.set_xticks(ticks, labels=[str(x) for x in ticks]); ax.minorticks_off()
        ax.axvline(1, color=GREY, ls="--", lw=.8)
        ax.axhline(9.90, color="#CBD2D8", lw=.7)
        ax.set_yticks([])
        ax.set_xlabel("Odds ratio (log scale)")
        meta = endpoints.loc[endpoints.endpoint.eq(ep)].iloc[0]
        fig.text(left, .885, f"{'AB'[panel]}  {ep}", fontsize=9.0, weight="bold", color=INK)
        fig.text(left, .86, f"{int(meta.n_cases):,} cases / {int(meta.n_controls):,} population noncases", fontsize=8.5)
        p_ax.text(0, -1.05, "Source P", fontsize=8.2, weight="bold")
        values = data.loc[data.endpoint.eq(ep)]
        for _, r in values.iterrows():
            y = r.display_row
            b, lower, upper = num(r.OR_ALT), num(r.OR_lower95_normal_approx), num(r.OR_upper95_normal_approx)
            ax.errorbar(b, y, xerr=[[b-lower],[upper-b]], fmt="o", color=BLUE, markersize=4.7, capsize=2.7, lw=1.15, zorder=3)
            p_ax.text(0, y, p_string(r.source_p_string), va="center", fontsize=8.1, color=INK)
            check(f"S4_interval_fits_{ep}_{r.panel_index}", lo < lower <= b <= upper < hi, "Every stored interval lies within the common axis")
    fig.text(.025, .108, "ALT effects and approximate normal 95% CIs are from stored beta/SE; the displayed source P values retain the GWAS calibration and are not inferred from the CIs.", fontsize=8.3)
    fig.text(.025, .080, "Both endpoints compare cases with population noncases. This display does not estimate ophthalmopathy risk conditional on Graves’ disease.", fontsize=8.3)
    fig.text(.025, .052, "Variants were selected after the discovery analyses; the endpoints share FinnGen R13 with discovery. These are exploratory associations, not independent replication.", fontsize=8.3)
    fig.text(.025, .024, "Non-HLA representatives and HLA regional leads follow different selection rules. Nearest-gene labels are positional annotations, not causal-gene assignments.", fontsize=8.3)
    export(fig, "figureS4")


def audit_exports():
    for path in EXPORTS:
        check(f"export_nonempty_{path.name}", path.stat().st_size > 1000, {"bytes": path.stat().st_size})
        if path.suffix == ".svg":
            tree = ET.parse(path)
            text_count = len(list(tree.iter("{http://www.w3.org/2000/svg}text")))
            image_count = len(list(tree.iter("{http://www.w3.org/2000/svg}image")))
            check(f"editable_text_vector_{path.name}", text_count > 10 and image_count == 0, {"text_elements": text_count, "raster_image_elements": image_count})
        if path.suffix == ".png":
            with Image.open(path) as im:
                dpi = im.info.get("dpi", (0,0))
                check(f"PNG_300dpi_{path.name}", all(abs(x-300) < .1 for x in dpi), {"dimensions": im.size, "dpi": dpi})
    for relative, before in INPUTS.items():
        check(f"input_preserved_{relative}", sha(ROOT/relative) == before["sha256"], "SHA-256 unchanged during build")
    images = []
    for index in range(1,5):
        with Image.open(FIG/f"figureS{index}.png") as im:
            copy = im.convert("RGB")
            copy.thumbnail((1000, 960))
            images.append(copy)
    canvas = Image.new("RGB", (2000, 1920), "#DDE3E8")
    for k, im in enumerate(images):
        canvas.paste(im, ((k%2)*1000+(1000-im.width)//2, (k//2)*960+(960-im.height)//2))
    canvas.save(REVIEW/"supp_build_contact_sheet.png")


def main():
    build_s1(); build_s2(); build_s3(); build_s4()
    audit_exports()
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Display-only Matplotlib redraw; no models, association tests, package changes, or new scientific classification.",
        "inputs": INPUTS,
        "script": str(Path(__file__).relative_to(ROOT)),
        "script_sha256": sha(Path(__file__)),
        "exports": [{"path": str(p.relative_to(ROOT)), "sha256": sha(p), "bytes": p.stat().st_size} for p in EXPORTS],
    }
    (DATA/"supp_source_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    result = {
        "created_utc": manifest["created_utc"], "status": "PASS_DISPLAY_CHECKS",
        "checks_passed": sum(x["pass"] for x in CHECKS), "checks_total": len(CHECKS),
        "checks": CHECKS,
        "runtime": {"executable": sys.executable, "python": platform.python_version(), "matplotlib": matplotlib.__version__, "pandas": pd.__version__, "numpy": np.__version__},
        "run_command": "uv run --no-sync --python /usr/bin/python3 python outputs/sci_writing_20261017/revisions/matplotlib_t19_t20_v2/scripts/build_supplement_figures.py",
        "runtime_note": "Bundled Python was checked and lacked Matplotlib; existing /usr/bin/python3 used with parent authorization. uv reports no project and no synchronization was performed.",
        "verification_boundary": "Source/display checks and AI visual inspection only; no model rerun, raw-GWAS re-audit, scientific approval, or human independent verification.",
    }
    (REVIEW/"supp_build_validation.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "checks": result["checks_total"], "exports": len(EXPORTS)}))


if __name__ == "__main__":
    main()
