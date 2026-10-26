# Shared and Divergent AITD Regulatory Architecture: analysis pipeline

## Current scope

This pipeline is matched to the first downloaded batch in this project:

- eight FinnGen R13 GWAS endpoints and their SuSiE/FINEMAP outputs;
- 33-cell-type FinnGen multiome eQTL and caQTL ACAT results;
- tabix-indexed nominal eQTL and caQTL summary statistics;
- feature-level multiome SuSiE credible-set summaries;
- peak-to-gene links, CASCADE annotations, consensus ATAC peaks, and official R12 colocalization results.

The original FinnGen GWAS and QTL batch uses GRCh38 and ALT-oriented effects. Subsequent external BBJ Graves/Hashimoto sources are GRCh37; both now have independently verified GRCh38 SNP derivatives. See `T16_T17_FEASIBLE_VARIANT_VALIDATION.md` for the current harmonization and feasible variant-validation gate; the earlier preparation report is historical.

The current batch does **not** contain the 2025 large hypothyroidism GWAS, cis-pQTL, drug databases, or FinnGen target PheWAS data. GTEx thyroid v10/v11 tables are present but lack sufficient candidate-region statistics for multi-signal colocalization. Four HRA001684 author-filtered thyroid libraries are available for descriptive summaries but lack a thyroid control group. Local GSE182416, GSE218743, GSE248205, and GSE285196 archives are also available for the partial T14 orthogonal-validation gate, with sample/donor aggregation and conservative phenotype-label handling documented in `documentation/analysis/T14_ORTHOGONAL_VALIDATION.md`.

## Scripts

| Script | Purpose | Main output |
|---|---|---|
| `00_validate_inputs.R` | Check expected files, headers, indexes, and batch composition without hashing 120 GB again | `outputs/qc/input_validation.tsv` |
| `01_harmonize_gwas.R` | Standardize GRCh38 variant IDs, ALT-oriented effects, QC fields, sample sizes, and HLA flags | `outcome_data/harmonized_gwas/*.tsv.gz` |
| `02_define_loci.R` | Construct a distance-based union of primary-trait genome-wide significant loci and preliminary directional classes | `outcome_data/loci/primary_shared_divergent_loci.tsv` |
| `03_prepare_ld_methods.R` | Prepare HapMap3-compatible LDSC inputs and LAVA summary-statistic inputs | `outcome_data/{ldsc,lava}/input/` |
| `04_run_ldsc.sh` | Run cross-trait LDSC using separately configured regression-LD and weight-LD prefixes | `outcome_data/ldsc/results/` |
| `04_build_lava_overlap.R` | Parse reviewed LDSC self/cross intercepts into the LAVA sampling-correlation matrix | `outcome_data/lava/input/sample.overlap.txt` |
| `05_run_lava.R` | Run local univariate and bivariate LAVA after an rsID-compatible reference and overlap matrix are supplied | `outcome_data/lava/results/` |
| `06_filter_official_r12_coloc.R` | Retain the flagship paper's strong R12 colocalizations as a benchmark | `outcome_data/coloc_official_r12/` |
| `07_coloc_multiome_r13.R` | Regionally query nominal QTL statistics and screen R13 GWAS–QTL colocalization by cell type | `outcome_data/coloc_r13/{eQTL,caQTL}/` |
| `08_integrate_regulatory_evidence.R` | Pair Graves/AIHT colocalizations, map caPeaks to genes, and attach CASCADE/LAVA evidence | `outcome_data/integrated/` |
| `09_plan_targeted_susie.R` | Convert Step 07 hits into priority-stratified per-variant Multiome SuSiE files and Tabix regions | `outcome_data/targeted_susie/` |
| `10_acquire_targeted_susie.sh` | Dry-run by default; remotely extract planned loci or download complete chromosome files only with `--run` | targeted slices or `raw_data/finngen_multiome/` |
| `11_verify_targeted_susie.R` | Audit targeted slices/full files, indexes, region counts, sizes, and optional MD5 | `outcome_data/targeted_susie/verification_*.tsv` |
| `47_run_t13_regulatory_cascade.R` | Reconstruct the frozen P1A GWAS–caQTL–peak–gene–eQTL chains with bilateral Tier 1 and direction gates | `outcome_data/t13/` |
| `48_run_t14_orthogonal_validation.py` | Audit available thyroid single-cell/spatial archives for T13 genes at sample/donor level while preserving GTEx/HRA001684 HOLDs | `outcome_data/t14/` |
| `49_run_t15_external_cohort_audit.py` | Compile the external cohort registry and overlap matrix with explicit FinnGen/UKB exclusion and unresolved metadata | `outcome_data/t15/` |
| `50_acquire_t17_bbj.sh` | Resume the official BBJ Graves archive download, then invoke offline preflight | `raw_data/external_gwas/bbj_2019_graves/` |
| `50_preflight_t17_bbj.py` | Verify archive CRC and scan combined-sex files; retain variant-preparation/build/EAS-LD HOLDs, no replication | `outcome_data/t17/` |
| `51_audit_t16_availability.py` | Rebuild AIHT availability evidence including BBJ-only GCST90018635 metadata/prefix and T15 overlap addendum; no full GWAS download | `outcome_data/t16/` |
| `20_run_t06_gwas_gwas_susie_coloc.R` | Audit official R13 GWAS SuSiE/FINEMAP signals and run direct multi-signal Graves–AIHT colocalization | `outcome_data/t06/` |
| `21_summarize_t06.R` | Annotate T06 signal pairs and generate high-confidence tables and QC figures | `outcome_data/t06/`, `outcome_figures/t06/` |
| `22_run_t07_five_classification.R` | Merge strong-H4 signal pairs and apply the preregistered five-class model using H3/H4, equivalence, and case-difference gates | `outcome_data/t07/` |
| `23_summarize_t07.R` | Summarize five-class results, threshold sensitivity, T04B resolution, and publication figures | `outcome_data/t07/`, `outcome_figures/t07/` |
| `24_review_t08_targeted_susie_plan.R` | Reconcile the core-lineage P0 plan with T07, promote the four divergent loci, separate HLA, and run manifest/build checks | `outcome_data/t08/` |
| `25_preflight_t08_remote_tabix.sh` | Perform retry-bounded read-only remote Tabix point queries for every frozen source-locus target | `outcome_data/t08/remote_tabix_preflight*.tsv` |
| `26_finalize_t08_preflight.R` | Merge primary and targeted-retry preflights and enforce the 32/32 T08 completion gate | `outcome_data/t08/t08_{validation_summary,decision}.tsv` |
| `27_prepare_t09_batches.R` | Load-balance a frozen non-HLA P0 plan into mutually exclusive, restartable extraction batches | `outcome_data/targeted_susie_t08/non_hla/*plans/` |
| `28_summarize_t09.R` | Join both T09 verification reports to the frozen plans, calculate local MD5 values, test BGZF integrity, inventory partials, and enforce the 30-file gate | `outcome_data/t09/` |
| `29_prepare_t09_hla_chunks.R` | Split each frozen HLA interval into an exact, non-overlapping transport-only partition and write one-row restartable plans | `outcome_data/targeted_susie_t08/hla/chunk_transport/` |
| `30_run_t09_hla_chunks.sh` | Run HLA transport chunks with a bounded dynamic worker queue and preserve per-chunk attempt logs | HLA chunk data and logs |
| `31_merge_t09_hla_chunks.sh` | Merge ordered HLA transport chunks, retain one header, rebuild indexes, and check original-region record counts | `outcome_data/targeted_susie_t08/hla/slices/` |
| `run_core_pipeline.sh` | Print the core execution plan by default; execute only with `--run` | no output in dry-run mode |

## Required software

Additional preparation entrypoints (no replication statistics in these scripts): `52_acquire_t16_bbj.sh` / `52_qc_t16_bbj.py` acquire and fully scan BBJ Hashimoto; `53_clean_harmonize_t17_bbj.py` cleans and roundtrip-lifts BBJ Graves SNPs before FinnGen allele matching; `54_validate_t17_harmonization.py` performs independent full-stream output verification. These use Python standard-library SQLite, curl, and the existing `tools/liftOver` with both UCSC chain files. Successful outputs are `outcome_data/t16/full_qc/run_20260904T074623Z/` and `outcome_data/t17/harmonization_resume_20260904/`. Interrupted outputs are preserved and not reused as completed results; rerun into a new directory.

The subsequent feasible-validation workflow is `55` → `58` for T16 harmonization/full audit, `60` for fresh T17 audit/hash binding, and `56` for the execution-time statistical freeze. Script `61` records an input-routing/hash-binding amendment without altering frozen targets/statistics; unchanged script `57` executes this bound plan, `59` independently verifies numbers using base R, and `62` writes the report. T16 has 4,684,134 aligned SNPs (24/24 checks); the frozen 20 tests yield 10 evaluable rows, three Graves same-direction variant associations, seven inconclusive/low-precision classifications and ten untestable rows. Outputs are `outcome_data/t16/harmonization_verified_20260904/` and `outcome_data/external_validation/run_20260904/`; R verification passed 29/29. See `codes/README.md` for exact entrypoints. Independent EUR sources, EAS LD, signal-level coloc, external T07 class validation and full T18 synthesis remain HOLD. Source-P association evidence is distinct from normal-approximation Q/TOST/precision calculations; no cross-disease BBJ contrast or ancestry-pooled effect is calculated.

T18 limited synthesis now uses scripts `63`–`65` at the verified variant-result boundary, without repeating raw-data preparation. It preserves the 20-test family, groups results by the original discovery class, audits shared-opposite candidates conservatively, and generates paired same-trait FinnGen–BBJ forest plots and a claim-disposition ledger. The formal directory is `outcome_data/t18/limited_synthesis_final_20260904/`; independent base-R checks passed 29/29 and separate synthesis/plot-data validation passed 67/67. Figures use ordinary Wald 95% CIs and source-P decisions, not SPA-calibrated intervals or pooled effects. Independent review closure is separate from the execution-time decision. See `T18_VARIANT_SYNTHESIS.md`; full signal-level T18 remains HOLD and T19 has not started.

T18 Japanese/EAS LD eligibility is audited offline by script `66`, using small public metadata snapshots, not genotype bodies. The resource audit passed 18/18 checks and three boundary tests; Phase3 JPT/EAS subsets are technical-pilot candidates only. The 30x manifest requires relatedness-aware subsetting, BBJ processed references are controlled, ToMMo genetic maps are not signed LD, and jMorp co-occurrence requires format/coverage verification. The GSU052 proposal is result-aware technical QC, not an association or signal-level execution plan. Formal output: `outcome_data/t18/eas_ld_audit_20260904/results/`; details: `T18_EAS_LD_RESOURCE_AUDIT.md`.

The scripts do not install packages. Prepare a controlled R environment containing:

- CRAN: `data.table`, `yaml`, `coloc`;
- Bioconductor: `Rsamtools`, `GenomicRanges`, `IRanges`;
- GitHub: `josefin-werme/LAVA` for local genetic correlation.

The targeted SuSiE acquisition step additionally uses HTSlib `tabix` and `bgzip` for remote indexed extraction, plus `curl` and `md5sum` for complete-file downloads. The system HTSlib 1.19 build lacks HTTPS support. T08 therefore added a project-local HTSlib 1.19 build linked to `libcurl` under `tools/htslib-1.19-curl/`; the acquisition and preflight scripts prefer these binaries and fall back to `PATH` only when they are absent.

Primary software documentation: [LAVA](https://github.com/josefin-werme/LAVA), [LDSC](https://github.com/bulik/ldsc), and [coloc data structures](https://chr1swallace.github.io/coloc/articles/a02_data.html).

LDSC is external and should be run in its own documented environment. Its directory, European LD-score directory, and HapMap3 allele list are passed to `04_run_ldsc.sh` as environment variables.

## Execution order

The following command only prints the steps and is safe for reviewing the plan:

```bash
bash codes/run_core_pipeline.sh --dry-run
```

When ready, run the initial data-derived modules:

```bash
Rscript codes/00_validate_inputs.R
Rscript codes/01_harmonize_gwas.R
Rscript codes/02_define_loci.R
Rscript codes/03_prepare_ld_methods.R
Rscript codes/06_filter_official_r12_coloc.R
```

The R13 colocalization screen is restartable by cell type. A small pilot should precede the full 33-cell-type analysis:

```bash
Rscript codes/07_coloc_multiome_r13.R \
  --workers 4 \
  --hla_workers 1 \
  --modalities eQTL,caQTL \
  --cell_types l2.B_naive,l2.B_memory,l2.Plasmablast,l2.Treg,l2.CD4_TCM,l2.CD14_Mono
```

Run all cells only after checking pilot row counts, error logs, locus sizes, and posterior behavior:

```bash
Rscript codes/07_coloc_multiome_r13.R --workers 4 --hla_workers 1
```

The parallel unit is one modality-cell-type pair. Non-HLA tasks use dynamic fork scheduling with four workers by default; HLA tasks are run in a separate one-worker phase because the 9 Mb region can produce much larger decompressed Tabix slices. Every worker is restricted to one `data.table`/BLAS thread. Worker results are written to independent temporary gzip files, combined by the parent process, and atomically renamed to the final cell-type output. A task-level failure therefore leaves no apparently completed final file.

Existing cell-type outputs are skipped. Use `--force` only when replacement is intentional. `--workers` overrides the configured non-HLA worker count, `--hla_workers` overrides the HLA limit, and `--loci` accepts comma-separated locus IDs for targeted reruns. Do not launch overlapping manual runs for the same modality-cell-type output.

## LDSC and LAVA gates

For LDSC:

```bash
PROJECT_ROOT="/path/to/AutoimmuneThyroid" \
LDSC_DIR="/path/to/ldsc" \
LDSC_PYTHON="/path/to/ldsc/environment/python" \
bash codes/04_run_ldsc.sh
```

The following inputs are already expanded and configured in `config/analysis.yml`:

- regression LD prefix: `raw_data/reference_panels/ldsc/LDscore/LDscore.`;
- regression-weight prefix: `raw_data/reference_panels/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.`;
- HapMap3 merge-allele list: `raw_data/reference_panels/ldsc/w_hm3.snplist`.

The two LD prefixes must remain separate. `04_run_ldsc.sh` reads them from YAML and validates all autosomes 1–22 before running. Environment variables `LDSC_REF_PREFIX`, `LDSC_WEIGHT_PREFIX`, and `HAPMAP3_SNPLIST` can override the YAML values when needed. `ldsc.ldsc_dir` and `ldsc.python` point to the user-level LDSC checkout and its isolated Python environment; `LDSC_DIR` and `LDSC_PYTHON` can override those YAML values at runtime.

Before LAVA:

1. confirm the configured 1000 Genomes EUR PLINK prefix and GRCh37 build in `lava.reference_prefix` and `lava.reference_build`;
2. record the reference build; the generated locus file contains semicolon-delimited rsID lists, so LAVA ignores coordinates and can use a different-build reference after SNP matching. If the `SNPS` column is removed, the reference must be GRCh38;
3. review the LDSC self/self and Graves/AIHT logs, then run `Rscript codes/04_build_lava_overlap.R`;
4. inspect `outcome_data/lava/input/sample.overlap.txt` and its provenance table before LAVA.

Do not set the Graves–AIHT overlap to zero. Both endpoints come from FinnGen and share many controls. An unreviewed overlap matrix can bias local covariance and correlation.

## Colocalization interpretation

The official multiome R12 results used `coloc.susie.direct`. The first local batch downloaded feature-level multiome credible-set summaries but not the multiome per-variant `.SUSIE.snp` files. Therefore `07_coloc_multiome_r13.R` uses regional nominal QTL P values, MAF, and an assumed default donor count of 1,108 with `coloc.abf`.

This is a screening analysis, not an exact reproduction of the flagship SuSiE colocalization. The output records this explicitly in `method` and `qtl_n_assumption`. Strong hits are configured as `PP.H4 >= 0.80`; feature-level credible-set log Bayes factor and purity fields are retained as supporting QC. If per-variant multiome SuSiE files become available, replace this screen with a multi-signal SuSiE/direct analysis.

The HLA region (chr6:25–34 Mb) is retained but flagged and ranked separately because long-range LD and multiple causal signals violate simple single-signal interpretations.

After Step 07, use `09_plan_targeted_susie.R` and the stage-gated acquisition workflow for per-variant Multiome fine-mapping data. Strong, HLA, observed multi-CS, and cross-disease opposite-direction candidates are mapped to exact manifest entries without downloading the 441 GB full collection. Remote interval extraction is the default recommendation because the current P1 full-chromosome plan would still require approximately 366.45 GiB. See `TARGETED_SUSIE.md` (historical source reference; not distributed) for selection rules, current counts, commands, overwrite guards, and verification.

T08 supersedes the pre-T07 core-lineage P0 queue with two frozen plans: 22 non-HLA files covering 24 source-locus intervals, followed by 8 HLA files covering 8 intervals in an isolated directory. All 32 targets passed real remote Tabix point queries. See `T08_TARGETED_SUSIE_PLAN_REVIEW.md` (historical source reference; not distributed); do not execute the legacy core P0 plan unchanged.

T09 acquired and independently verified all 30 frozen P0 targeted files. The final gate requires readable indexes, non-zero planned-region counts, full BGZF integrity, local data/index MD5 values, and no residual partial files. HLA was transferred in exact non-overlapping coordinate chunks because of remote throughput, then merged and re-indexed before validation against the original full HLA region. See `T09_TARGETED_SUSIE_ACQUISITION.md` (historical source reference; not distributed).

T13 is a separate downstream gate implemented by `codes/47_run_t13_regulatory_cascade.R`. It audits the frozen T11 P1A features using direct signal-level GWAS–QTL pairs, exact significant peak–gene links, CASCADE supporting annotations, and Ensembl release-116 symbols. Both Graves and AIHT must independently pass the Tier 1 PIP/credible-set gate; all three direction relationships are retained. T13 does not clear the T12 HOLD, and unavailable pQTL/coding inputs remain explicit HOLD classes.

T14 uses `codes/48_run_t14_orthogonal_validation.py`. It is currently a partial gate: GSE182416 normal-thyroid localization, GSE218743 sample-level supporting summaries, GSE248205 spatial archive/candidate-spot inventory, GSE285196 targeted-panel audit, and four HRA001684 thyroid-library summaries (56,199 cells) are reproducible locally. GTEx multi-signal colocalization and HRA001684 donor-level disease-state comparisons remain explicit HOLDs because candidate-region statistics and a thyroid control group are insufficient. No cell-level pseudoreplication, disease-state effect, spatial-neighborhood mechanism, or causal claim is promoted. The 2026-09-04 UTC rerun passed 11/11 formal checks, reproduced all ten TSV outputs byte-for-byte, and passed the supplemental numeric/barcode checks recorded in `outcome_data/t14/rerun_20260904/verification.json`.

## Disease-classification safeguards

The 2026-09-04 GSU052 external-LD technical pilot uses `codes/67`–`72` and is frozen in `outcome_data/t18/gsu052_ld_pilot_20260904/`. Bounded Phase 3 genotype acquisition and per-SNP hg19/hg38 REF/roundtrip checks yielded 5,008 jointly aligned variants; JPT/EAS signed matrices share 4,152 QC-passing variants. Independent preparation (17/17) and full LD/ridge checks (39/39) passed. Fixed-ridge diagnostic flags are panel/regularization-sensitive and uncalibrated for BBJ binary GWAS: they must not select a preferred panel/lambda or justify post-hoc SNP removal. JPT is nested within EAS. See GSU052 pilot report: `T18_GSU052_LD_TECHNICAL_PILOT.md` (historical source reference; not distributed).

The subsequent BBJ input/applicability prerequisite audit uses `codes/73`–`74` and is frozen in `outcome_data/t18/bbj_input_applicability_audit_20260904/`. Both full source files passed field/numerical auditing, and corrected independent validation passed 16/16. Graves retains SAIGE/SPA fields and per-SNP N; AIHT lacks per-SNP N and raw SAIGE diagnostics. Published Graves beta/SE tracks source SPA P, not the unadjusted score P. Technical LD success does not establish the joint RSS model or BBJ cohort LD. Status is `PASS_PREREQUISITE_AUDIT_NO_GO_SIGNAL_ANALYSIS`; full T18 and prior gates remain HOLD. Next is public availability auditing for the original hum0197/PheWeb AIHT SAIGE fields and release-matched BBJ signed LD, not a signal run. See BBJ input audit: `T18_BBJ_INPUT_APPLICABILITY_AUDIT.md` (historical source reference; not distributed).

`02_define_loci.R` deliberately labels one-trait significance as `*_threshold_only_provisional`. A non-significant P value in the other disease does not establish disease specificity, especially when case counts differ markedly and controls overlap.

The defensible hierarchy is:

1. global genetic correlation from LDSC;
2. local genetic covariance/correlation from LAVA with sample overlap accounted for;
3. aligned R13 fine-mapping and cell-type-specific colocalization;
4. peak-to-gene/CASCADE support;
5. external GWAS validation only after the T15 cohort-overlap and metadata gates are satisfied, plus thyroid-tissue validation; large hypothyroidism meta-analyses containing FinnGen or UKB are comparison-only under the frozen project rule.

Only after these layers agree should loci be described as shared, opposing, or disease-dominant. “Coloc in one disease only” remains evidence asymmetry, not proof of trait specificity.

## Later modules

After the next data batches are available, add separate modules for:

- external hypothyroidism GWAS replication and cross-study harmonization;
- GTEx thyroid eQTL colocalization;
- Graves/Hashimoto thyroid single-cell expression validation;
- cis-eQTL/cis-pQTL Mendelian randomization with colocalization gating;
- Open Targets/ChEMBL target tractability;
- FinnGen PheWAS safety and other-autoimmune-disease checks.

MR should not be run merely because a cis-QTL is significant. Require allele harmonization, adequate instrument strength, locus colocalization, and sensitivity to horizontal pleiotropy.
