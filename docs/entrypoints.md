# Ordered entrypoints and source coverage

The order below describes the source lineage. It is not a turnkey execution script. Restore the separately acquired inputs and stage-specific frozen contracts first. All paths are relative to the repository/analysis root.

| Stage | Source entrypoints | Role and execution boundary |
|---|---|---|
| Data acquisition and structural QC | `codes/download_finngen_batch1.py`, `download_batch2.py`, `00_validate_inputs.R` | Author acquisition code and checks. Manifests and acquired input bodies are absent. Downloaders are not dry-run wrappers. |
| General GWAS processing | `01_harmonize_gwas.R`, `02_define_loci.R` | General autosomal QC and locus definitions; released fine-mapping component auditing uses its own scope. |
| LDSC/LAVA preparation and historical core | `03_prepare_ld_methods.R`, `04_run_ldsc.sh`, `04_build_lava_overlap.R`, `05_run_lava.R` | Requires separate LD resources, LDSC environment and reviewed sampling-correlation input. |
| Regional molecular screening | `06_filter_official_r12_coloc.R`, `07_coloc_multiome_r13.R`, `08_integrate_regulatory_evidence.R` | Historical regional screening, not direct multi-signal proof. |
| Path B power and CC-GWAS | `12_power_equivalence_path_b.R`, `13_prepare_ccgwas_inputs.R`, `14_run_ccgwas_path_b.R`, `15_summarize_ccgwas_path_b.R` | Frozen prevalence/overlap sensitivity; official external CC-GWAS source required. |
| T05 refined LD analyses | `16_run_t05_ldsc.sh`, `17_parse_t05_ldsc.R`, `18_run_t05_lava.R`, `19_summarize_t05.R` | Manuscript-era non-HLA LDSC/LAVA and HLA separation. |
| T06 component-pair comparisons | `20_run_t06_gwas_gwas_susie_coloc.R`, `21_summarize_t06.R` | Official FinnGen SuSiE components; normalized full-region LBF coverage over shared coordinate span. |
| T07 classification | `22_run_t07_five_classification.R`, `23_summarize_t07.R` | Organizational analysis units, marginal effect direction and stored equivalence criteria; preserve unresolved status. |
| Targeted molecular acquisition | `09_plan_targeted_susie.R`, `10_acquire_targeted_susie.sh`, `11_verify_targeted_susie.R`, `24`-`31` | Gated plans, Tabix slices, ordinary/HLA transport and verification. Acquisition script defaults to dry run; review all individual CLIs. |
| T10/T11 direct GWAS-QTL pairs | `32_run_t10_qtl_susie_coloc.R`, `33_summarize_t10.R`, `34`-`39` | T11B wave 2 uses script 32 with profile `t11b_wave2`; exact matched-ID alpha mass differs from T06 coverage. |
| T12 partial stratified LDSC | `40`-`46` | Graves/AIHT global-cell, global-lineage and within-lineage families ran; missing AITD-union/direct case-case resources retain full-scope HOLD. |
| T13 regulatory chains | `47_run_t13_regulatory_cascade.R` | Cross-cell eQTL candidate matching is permitted by the implemented method; no complete chain established. |
| T14 contextual evidence | `48_run_t14_orthogonal_validation.py` | Stage-specific public thyroid datasets and aggregated donor/sample/library outputs; no raw matrices distributed. |
| T15-T17 availability/QC/harmonization | `49`-`55`, `58_validate_t16_harmonization.py`, `60_revalidate_t17_binding.py` | Source availability, full QC and roundtrip liftOver/allele alignment; source snapshots, databases and chains absent. |
| Frozen external variant checks | `56_freeze_external_validation.py`, `61_bind_external_execution_plan.py`, `57_run_external_validation.py`, `59_verify_external_statistics.R`, `62_report_external_validation.py` | Execution-time 20-test freeze and exact selected-variant statistics. Follow the bound plan, not old interrupted directories. |
| T18 limited synthesis | `63_synthesize_t18.py`, `64_plot_t18_forest.R`, `65_validate_t18.py` | Uses verified variant-level result boundary; no pooled estimate or external signal-class validation. |
| T19 descriptive HLA | `outputs/t19_t24_20261017/t19/code/run_t19_hla_descriptive.R`, `validate_t19_hla_descriptive.R`, `hash_t19_outputs.R` | Actual extra-directory entrypoints. Runner takes project root and its exact output directory; restore required subfolders/inputs before use. |
| T20 ophthalmopathy | `outputs/t19_t24_20261017/t20/t20_lookup.R` (`freeze`, then `lookup`), `validate_t20.R`, `finalize_t20.R`, `select_panel.awk` | Frozen 12-variant by 2-endpoint lookup; finalization binds unchanged HLA selections. Population noncase controls. |

## Final main and supplementary figures

1. Current Figure 1: `outputs/sci_writing_20261017/plos_one_study_flow_20261025/scripts/build_workflow.py`; browser PDF helper is `render_workflow.py`. This is a diagram builder with embedded labels/counts, not an analysis.
2. Current Figures 2 and 4: `revisions/matplotlib_clean_v3_20261018/scripts/build_main_figures.py` within the writing directory. The same builder also draws the older Figure 2 (now Figure 3), which is superseded by the next entry. The included historical source dependencies are base `scripts/build_figures.py`, `revisions/circular_v1/scripts/build_circular_figures.R`, `revisions/matplotlib_t19_t20_v2/scripts/build_main_figures.py`, and `revisions/figure1_panelB_class_series_20261018/scripts/build_figure1_class_series.py`. Old frozen TSV extracts, baseline PNG/PDF files and captions required by preservation audits are absent.
3. Current Figure 3: `bmc_submission_figure2_v3_20261021/scripts/build_figure2.py`. Its original output name is `Figure_2`; the final current figure is renamed Figure 3. It reads frozen T06 nodes/edges and catalogs. Caption/protected-package dependencies are absent.
4. S1-S4 Figs: clean-v3 `scripts/build_supplement_figures.py`, with its previous v2 builder included for provenance. It reads stored molecular, BBJ, CC-GWAS and T20 outputs; those result tables are absent.
5. Final Arial export: `plos_one_submission_arial_v2_20261025/scripts/export_arial_main_figures.py`, importing adjacent `normalize_pdf_text_sizes.py` and invoking adjacent `render_arial_figures_browser.py`. The four exact input SVGs are retained in `plos_one_submission_20261025/reference_artwork/`; the four final SVGs are in the Arial package `figures/`. Supply genuine Arial Regular/Bold fonts in that package's `assets/fonts/` and an existing Playwright-capable interpreter through `PLAYWRIGHT_PYTHON`. Fonts, Chromium and fontconfig tools are not distributed. The raw PDF Tf-size correction changes only validated numeric tokens and retains reverse-substitution checks.

The earlier `plos_one_submission_20261025/scripts/export_main_figures.py` is included to document the Figure 1 panel-letter adjustment and historical raster-export lineage. It requires undistributed historical PDFs and is not required when starting from the included exact pre-Arial SVGs.

## Tables

Base `outputs/sci_writing_20261017/scripts/build_tables.py` assembles stored result tables; `revisions/matplotlib_t19_t20_v2/scripts/build_tables.py` extends the retained circular-package table payloads with T19/T20. The numerical result inputs, publisher snapshots and historical table JSON/workbooks are absent. Manuscript/Word assembly, submission packaging and author-contact tools are not part of this code release.
