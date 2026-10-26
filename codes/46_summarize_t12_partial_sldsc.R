#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages("data.table")

t12_dir <- file.path(project_root, "outcome_data", "t12")
result_dir <- file.path(t12_dir, "results")
traits <- c("Graves", "AIHT")
models <- c("global_cell", "global_lineage", "within_lineage")
expected_rows <- c(global_cell = 28L, global_lineage = 6L, within_lineage = 28L)

results <- list()
for (trait in traits) {
  for (model in models) {
    path <- file.path(result_dir, paste0(trait, ".", model, ".cell_type_results.txt"))
    if (!file.exists(path)) stop("Missing h2-cts result: ", path, call. = FALSE)
    x <- data.table::fread(path)
    if (nrow(x) != expected_rows[[model]]) stop("Unexpected result rows for ", trait, "/", model, call. = FALSE)
    x[, `:=`(trait = trait, model = model)]
    results[[paste(trait, model)]] <- x
  }
}
all_results <- data.table::rbindlist(results, use.names = TRUE)
data.table::setnames(
  all_results,
  c("Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"),
  c("annotation", "coefficient", "coefficient_se", "p_one_sided")
)
all_results[, z := coefficient / coefficient_se]
all_results[, q_bh_model_family := p.adjust(p_one_sided, method = "BH"), by = .(trait, model)]
all_results[, q_bh_trait_all_models := p.adjust(p_one_sided, method = "BH"), by = trait]
all_results[, positive_fdr05_model_family := coefficient > 0 & q_bh_model_family < 0.05]
data.table::setcolorder(all_results, c(
  "trait", "model", "annotation", "coefficient", "coefficient_se", "z",
  "p_one_sided", "q_bh_model_family", "q_bh_trait_all_models", "positive_fdr05_model_family"
))
data.table::setorder(all_results, trait, model, p_one_sided)
write_tsv(all_results, file.path(t12_dir, "t12_h2cts_all_results.tsv"))

lineage_by_cell <- c(
  B_intermediate = "B", B_memory = "B", B_naive = "B", Plasmablast = "B",
  CD4_CTL = "T", CD4_Naive = "T", CD4_Proliferating = "T", CD4_TCM = "T",
  CD4_TEM = "T", CD8_Naive = "T", CD8_TCM = "T", CD8_TEM = "T",
  MAIT = "T", Treg = "T", dnT = "T", gdT = "T",
  CD14_Mono = "Mono", CD16_Mono = "Mono", cDC1 = "DC", cDC2 = "DC", pDC = "DC",
  NK = "NK", NK_CD56bright = "NK", NK_Proliferating = "NK",
  HSPC = "Other", ILC = "Other", Platelet = "Other", Eryth = "Other"
)
pick <- function(model, prefix) {
  model_value <- model
  all_results[model == model_value, .(
    trait, annotation,
    coefficient = coefficient,
    p = p_one_sided,
    q = q_bh_model_family,
    positive_fdr05 = positive_fdr05_model_family
  )][, setNames(.SD, c("trait", "annotation", paste0(prefix, c("_coefficient", "_p", "_q", "_positive_fdr05"))))]
}
global_cell <- pick("global_cell", "global")
within <- pick("within_lineage", "within")
lineage <- pick("global_lineage", "lineage")
interpretation <- merge(global_cell, within, by = c("trait", "annotation"), all = TRUE)
interpretation[, lineage := unname(lineage_by_cell[annotation])]
interpretation <- merge(
  interpretation,
  lineage,
  by.x = c("trait", "lineage"), by.y = c("trait", "annotation"), all.x = TRUE
)
interpretation[, evidence_class := data.table::fifelse(
  within_positive_fdr05,
  "independent_source_cell_enrichment_within_lineage",
  data.table::fifelse(
    lineage_positive_fdr05,
    "lineage_shared_enrichment_without_independent_source_cell_evidence",
    data.table::fifelse(
      global_positive_fdr05,
      "global_source_cell_signal_not_robust_to_lineage_conditioning",
      "no_positive_FDR05_evidence"
    )
  )
)]
interpretation[, evidence_boundary := "source-labelled consensus peak association; not differential accessibility and not causal cell assignment"]
data.table::setorder(interpretation, trait, evidence_class, within_q, global_q)
write_tsv(interpretation, file.path(t12_dir, "t12_cell_type_evidence_classes.tsv"))

summary_counts <- interpretation[, .N, by = .(trait, evidence_class)]
write_tsv(summary_counts, file.path(t12_dir, "t12_evidence_class_counts.tsv"))

resource_gate <- data.table::fread(file.path(t12_dir, "reference_resource_gate.tsv"))
annotation_validation <- data.table::fread(file.path(t12_dir, "annotations", "annotation_build_validation.tsv"))
split_validation <- data.table::fread(file.path(t12_dir, "ldscore_split_validation.tsv"))
phenotypes <- data.table::fread(file.path(t12_dir, "phenotype_scope.tsv"))
chr6_bim <- data.table::fread(
  file.path(project_root, "raw_data", "reference_panels", "ldsc", "1000G_EUR_Phase3_plink", "1000G.EUR.QC.6.bim"),
  header = FALSE, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2")
)
chr6_ann <- data.table::fread(file.path(t12_dir, "annotations", "t12_source_labels.6.annot.gz"))
mhc_rows <- which(chr6_bim$BP >= 25000000L & chr6_bim$BP <= 34000000L)

validation <- data.table::data.table(
  check = c(
    "reference_resources_all_ready", "annotations_22_chromosomes", "annotation_columns_35",
    "MHC_annotation_rows_all_zero", "custom_ldscores_22_chromosomes",
    "split_ldscore_files_complete", "eligible_h2cts_runs_6", "all_coefficients_finite",
    "AITD_union_input_ready", "direct_case_case_input_ready"
  ),
  passed = c(
    all(resource_gate$ready), nrow(annotation_validation) == 22L,
    all(annotation_validation$annotation_columns == 35L),
    all(rowSums(as.matrix(chr6_ann[mhc_rows])) == 0),
    sum(file.exists(file.path(t12_dir, "ldscores", "master", paste0("t12_source_labels.", 1:22, ".l2.ldscore.gz")))) == 22L,
    all(split_validation$split_files_observed == split_validation$split_files_expected),
    length(results) == 6L, all(is.finite(all_results$coefficient)) && all(is.finite(all_results$coefficient_se)),
    phenotypes[requested_analysis == "AITD_combined_case_control", usable_input],
    phenotypes[requested_analysis == "Graves_vs_AIHT_direct_case_case", usable_input]
  )
)
write_tsv(validation, file.path(t12_dir, "t12_final_validation.tsv"))
decision <- data.table::data.table(
  gate = "T12_global_cell_type_enrichment",
  status = "HOLD",
  decision = "PARTIAL_GRAVES_AIHT_EXECUTED_HOLD_MISSING_AITD_UNION_AND_DIRECT_CASE_CASE",
  eligible_runs_passed = "6/6",
  completed_method_families = "global_cell;global_lineage;within_lineage",
  unavailable_required_analyses = "AITD_combined_case_control;Graves_vs_AIHT_direct_case_case",
  rationale = paste0(
    "All reference, annotation, MHC, LD-score, and eligible-trait checks passed. ",
    "The frozen T12 completion scope remains incomplete because no AITD-union GWAS or direct Graves-vs-AIHT GWAS is available; CC-GWAS was not substituted."
  ),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_tsv(decision, file.path(t12_dir, "t12_final_decision.tsv"))
log_message("T12 partial S-LDSC summary complete: Graves/AIHT executed; full T12 remains HOLD.")
