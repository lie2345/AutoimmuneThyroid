#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)

out_root <- file.path(project_root, cfg$project$outcome_dir)
t11_dir <- file.path(out_root, "t11")
t11a_dir <- file.path(out_root, "t11a")
direct_dir <- file.path(t11a_dir, "direct_coloc")
dir.create(t11a_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  screen = file.path(direct_dir, "regional_screening_rows.tsv.gz"),
  pairs = file.path(direct_dir, "gwas_qtl_susie_signal_pairs.tsv.gz"),
  feature_qc = file.path(direct_dir, "qtl_feature_qc.tsv.gz"),
  wave1_plan = file.path(t11_dir, "p1a_wave1_acquisition_plan.tsv"),
  full_plan = file.path(t11_dir, "p1a_acquisition_plan.tsv"),
  full_features = file.path(t11_dir, "p1a_candidate_features.tsv"),
  full_rows = file.path(t11_dir, "p1a_candidate_rows.tsv.gz"),
  verification = file.path(t11a_dir, "p1a_wave1_verification.tsv"),
  preflight = file.path(t11a_dir, "remote_tabix_preflight.tsv"),
  acquisition_log = file.path(out_root, "targeted_susie_t11a", "p1a_wave1", "acquisition_log.tsv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing T11A inputs: ", paste(missing, collapse = ", "), call. = FALSE)

screen <- data.table::fread(paths$screen, showProgress = FALSE)
pairs <- data.table::fread(paths$pairs, showProgress = FALSE)
feature_qc <- data.table::fread(paths$feature_qc, showProgress = FALSE)
wave1_plan <- data.table::fread(paths$wave1_plan, showProgress = FALSE)
full_plan <- data.table::fread(paths$full_plan, showProgress = FALSE)
full_features <- data.table::fread(paths$full_features, showProgress = FALSE)
full_rows <- data.table::fread(paths$full_rows, showProgress = FALSE)
verification <- data.table::fread(paths$verification, showProgress = FALSE)
preflight <- data.table::fread(paths$preflight, showProgress = FALSE)
acquisition_log <- data.table::fread(paths$acquisition_log, showProgress = FALSE)

best_pair <- function(x) {
  if (!nrow(x)) return(data.table::data.table())
  if (any(x$classification == "shared_strong_H4")) {
    x[classification == "shared_strong_H4"][which.max(PP.H4)]
  } else if (any(x$classification == "distinct_signals_H3")) {
    x[classification == "distinct_signals_H3"][which.max(PP.H3)]
  } else if (any(x$classification == "shared_suggestive_H4")) {
    x[classification == "shared_suggestive_H4"][which.max(PP.H4)]
  } else {
    x[which.max(PP.H4)]
  }
}

# Apply exactly the same screen-row reclassification logic used in T10.
reclass_rows <- lapply(seq_len(nrow(screen)), function(i) {
  s <- screen[i]
  x <- pairs[screen_id == s$screen_id]
  best <- best_pair(x)
  has_strong <- nrow(x) && any(x$classification == "shared_strong_H4")
  has_h3 <- nrow(x) && any(x$classification == "distinct_signals_H3")
  direction_match <- if (!nrow(best) || !has_strong || is.na(best$effect_direction_if_shared) ||
                         is.na(s$effect_direction) || s$effect_direction == "undetermined") {
    NA
  } else identical(best$effect_direction_if_shared, s$effect_direction)
  status <- if (!nrow(x)) {
    "untestable"
  } else if (has_strong && !identical(direction_match, FALSE)) {
    "confirmed"
  } else if ((has_strong && identical(direction_match, FALSE)) || (!has_strong && has_h3)) {
    "conflicting"
  } else {
    "downgraded"
  }
  reason <- if (status == "untestable") {
    "no_eligible_direct_signal_pair"
  } else if (status == "confirmed") {
    if (has_h3) "strong_H4_confirmed_with_additional_distinct_signal" else "strong_prior_robust_H4_and_direction_consistent"
  } else if (status == "conflicting") {
    if (has_strong) "direct_shared_signal_direction_conflicts_with_regional_screen" else "direct_signal_pair_supports_H3_not_H4"
  } else if (any(x$classification == "shared_suggestive_H4")) {
    "only_suggestive_direct_H4"
  } else {
    "direct_testable_but_no_qualified_H4_or_H3"
  }
  data.table::data.table(
    screen_id = s$screen_id, trait_id = s$trait_id, locus_id = s$locus_id,
    modality = s$modality, level = s$level, cell_type = s$cell_type, feature_id = s$feature_id,
    is_hla = s$is_hla, regional_PP.H3 = s$PP.H3.abf, regional_PP.H4 = s$PP.H4.abf,
    regional_top_variant = s$top_h4_variant, regional_effect_direction = s$effect_direction,
    n_direct_pairs = nrow(x), n_strong_H4_pairs = if (nrow(x)) sum(x$classification == "shared_strong_H4") else 0L,
    n_suggestive_H4_pairs = if (nrow(x)) sum(x$classification == "shared_suggestive_H4") else 0L,
    n_distinct_H3_pairs = if (nrow(x)) sum(x$classification == "distinct_signals_H3") else 0L,
    best_gwas_region = if (nrow(best)) best$gwas_region else NA_character_,
    best_gwas_cs = if (nrow(best)) best$gwas_cs else NA_integer_,
    best_qtl_cs = if (nrow(best)) best$qtl_cs else NA_integer_,
    direct_PP.H3 = if (nrow(best)) best$PP.H3 else NA_real_,
    direct_PP.H4 = if (nrow(best)) best$PP.H4 else NA_real_,
    sensitivity_PP.H3 = if (nrow(best)) best$sensitivity_PP.H3 else NA_real_,
    sensitivity_PP.H4 = if (nrow(best)) best$sensitivity_PP.H4 else NA_real_,
    direct_top_variant = if (nrow(best)) best$top_h4_variant else NA_character_,
    direct_effect_direction = if (nrow(best)) best$direction_at_top_h4 else NA_character_,
    same_top_variant_as_regional = if (nrow(best)) identical(best$top_h4_variant, s$top_h4_variant) else NA,
    direction_matches_regional = direction_match,
    direct_pair_classification = if (nrow(best)) best$classification else NA_character_,
    screen_reclassification = status, reclassification_reason = reason,
    final_evidence_method = if (status == "untestable") "unresolved_no_direct_SuSiE_evidence" else "direct_SuSiE_BF_signal_pair_coloc"
  )
})
reclass <- data.table::rbindlist(reclass_rows, fill = TRUE)
data.table::setorder(reclass, locus_id, modality, level, cell_type, feature_id, trait_id)
write_tsv(reclass, file.path(t11a_dir, "wave1_screen_reclassification.tsv.gz"))

graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
feature_keys <- c("locus_id", "modality", "level", "cell_type", "feature_id", "is_hla")
feature_units <- unique(reclass[, ..feature_keys])
feature_results <- data.table::rbindlist(lapply(seq_len(nrow(feature_units)), function(i) {
  f <- feature_units[i]
  x <- reclass[
    locus_id == f$locus_id & modality == f$modality & level == f$level &
      cell_type == f$cell_type & feature_id == f$feature_id
  ]
  g <- x[trait_id == graves_id][1]
  a <- x[trait_id == aiht_id][1]
  both_confirmed <- nrow(g) == 1L && nrow(a) == 1L &&
    g$screen_reclassification == "confirmed" && a$screen_reclassification == "confirmed"
  same_direction <- both_confirmed && !is.na(g$direct_effect_direction) && !is.na(a$direct_effect_direction) &&
    g$direct_effect_direction == a$direct_effect_direction && g$direct_effect_direction != "undetermined"
  status <- if (both_confirmed && same_direction) {
    "confirmed_both_same_direction"
  } else if (nrow(g) != 1L || nrow(a) != 1L) {
    "incomplete_trait_pair"
  } else if (any(c(g$screen_reclassification, a$screen_reclassification) == "conflicting")) {
    "conflicting"
  } else if (both_confirmed) {
    "confirmed_both_direction_discordant"
  } else {
    "downgraded_or_mixed"
  }
  data.table::data.table(
    locus_id = f$locus_id, modality = f$modality, level = f$level,
    cell_type = f$cell_type, feature_id = f$feature_id, is_hla = f$is_hla,
    graves_status = if (nrow(g)) g$screen_reclassification else NA_character_,
    aiht_status = if (nrow(a)) a$screen_reclassification else NA_character_,
    graves_PP.H4 = if (nrow(g)) g$direct_PP.H4 else NA_real_,
    aiht_PP.H4 = if (nrow(a)) a$direct_PP.H4 else NA_real_,
    graves_PP.H3 = if (nrow(g)) g$direct_PP.H3 else NA_real_,
    aiht_PP.H3 = if (nrow(a)) a$direct_PP.H3 else NA_real_,
    graves_direction = if (nrow(g)) g$direct_effect_direction else NA_character_,
    aiht_direction = if (nrow(a)) a$direct_effect_direction else NA_character_,
    same_direction = same_direction, feature_result = status
  )
}), fill = TRUE)
data.table::setorder(feature_results, locus_id, modality, level, cell_type, feature_id)
write_tsv(feature_results, file.path(t11a_dir, "wave1_feature_results.tsv"))

n_confirmed_rows <- sum(reclass$screen_reclassification == "confirmed")
n_dual_confirmed <- sum(feature_results$feature_result == "confirmed_both_same_direction")
continue_wave2 <- n_dual_confirmed >= 1L || n_confirmed_rows >= 4L

# Freeze the incremental wave-2 scope without acquiring it.
wave2_paths <- setdiff(full_plan$source_path, wave1_plan$source_path)
wave2_plan <- full_plan[source_path %in% wave2_paths]
wave2_features <- full_features[source_path %in% wave2_paths]
wave2_feature_key <- paste(wave2_features$locus_id, wave2_features$modality, wave2_features$level,
                           wave2_features$cell_type, wave2_features$feature_id, sep = "||")
full_rows[, feature_key__ := paste(locus_id, modality, level, cell_type, feature_id, sep = "||")]
wave2_rows <- full_rows[feature_key__ %in% wave2_feature_key]
wave2_rows[, feature_key__ := NULL]
write_tsv(wave2_plan, file.path(t11a_dir, "p1a_wave2_acquisition_plan.tsv"))
write_tsv(wave2_features, file.path(t11a_dir, "p1a_wave2_candidate_features.tsv"))
write_tsv(wave2_rows, file.path(t11a_dir, "p1a_wave2_candidate_rows.tsv.gz"))

# File-level provenance and integrity audit for the seven acquired slices.
inventory <- merge(
  wave1_plan[, .(source_path, modality, level, cell_type, chrom, tabix_regions)],
  verification[, .(source_path, local_data_path, local_index_path, observed_size_bytes,
                   n_records_in_planned_regions, index_ok, verification_status = status)],
  by = "source_path", all.x = TRUE, sort = FALSE
)
inventory[, bgzip_ok := vapply(local_data_path, function(p) {
  if (is.na(p) || !file.exists(p)) return(FALSE)
  identical(system2("bgzip", c("-t", p), stdout = FALSE, stderr = FALSE), 0L)
}, logical(1))]
inventory[, data_md5 := vapply(local_data_path, function(p) {
  if (is.na(p) || !file.exists(p)) return(NA_character_)
  unname(tools::md5sum(p))
}, character(1))]
inventory[, index_md5 := vapply(local_index_path, function(p) {
  if (is.na(p) || !file.exists(p)) return(NA_character_)
  unname(tools::md5sum(p))
}, character(1))]
write_tsv(inventory, file.path(t11a_dir, "wave1_file_inventory.tsv"))

part_roots <- c(t11a_dir, file.path(out_root, "targeted_susie_t11a"))
part_files <- unique(unlist(lapply(part_roots, function(d) {
  if (!dir.exists(d)) return(character())
  list.files(d, pattern = "\\.part([.]|$)", recursive = TRUE, full.names = TRUE, all.files = TRUE)
})))
part_audit <- data.table::data.table(path = part_files)
write_tsv(part_audit, file.path(t11a_dir, "residual_part_files.tsv"))

posterior_ok <- !nrow(pairs) || all(abs(rowSums(pairs[, .(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]) - 1) < 1e-6)
pair_key_unique <- !anyDuplicated(pairs[, .(screen_id, gwas_region, gwas_cs, qtl_cs)])
validation <- data.table::data.table(
  check = c(
    "remote_preflight_7_of_7_ok", "acquisition_log_7_completed", "verification_7_of_7_ok",
    "bgzip_integrity_7_of_7_ok", "no_residual_part_files", "all_22_screen_rows_present",
    "all_11_feature_units_present", "all_22_rows_have_direct_pairs", "signal_pair_keys_unique",
    "posterior_probabilities_sum_to_one", "zero_allele_mismatches_excluded",
    "all_22_rows_reclassified", "stopping_rule_evaluated", "wave2_frozen_9_files",
    "wave2_frozen_15_features", "wave2_frozen_30_trait_feature_rows"
  ),
  passed = c(
    nrow(preflight) == 7L && all(preflight$status == "ok"),
    nrow(acquisition_log) == 7L && all(acquisition_log$status == "completed"),
    nrow(verification) == 7L && all(verification$status == "ok") && all(verification$index_ok),
    nrow(inventory) == 7L && all(inventory$bgzip_ok),
    length(part_files) == 0L,
    nrow(screen) == 22L,
    nrow(feature_results) == 11L,
    nrow(reclass[n_direct_pairs > 0]) == 22L,
    pair_key_unique,
    posterior_ok,
    sum(pairs$allele_mismatches_excluded, na.rm = TRUE) == 0,
    nrow(reclass) == 22L && !anyNA(reclass$screen_reclassification),
    is.logical(continue_wave2) && length(continue_wave2) == 1L,
    nrow(wave2_plan) == 9L,
    nrow(wave2_features) == 15L,
    nrow(wave2_rows) == 30L
  )
)
write_tsv(validation, file.path(t11a_dir, "t11a_validation.tsv"))

pair_classes <- pairs[, .N, by = classification][order(classification)]
summary <- data.table::rbindlist(list(
  data.table::data.table(metric = c(
    "preflight_files_ok", "acquired_files_completed", "verified_files_ok", "acquired_records",
    "acquired_bytes", "screen_rows", "direct_signal_pairs", "feature_units",
    "confirmed_rows", "conflicting_rows", "downgraded_rows", "untestable_rows",
    "dual_disease_same_direction_confirmed_features", "wave2_files", "wave2_features", "wave2_rows"
  ), value = c(
    sum(preflight$status == "ok"), sum(acquisition_log$status == "completed"), sum(verification$status == "ok"),
    sum(verification$n_records_in_planned_regions), sum(verification$observed_size_bytes), nrow(screen), nrow(pairs),
    nrow(feature_results), n_confirmed_rows, sum(reclass$screen_reclassification == "conflicting"),
    sum(reclass$screen_reclassification == "downgraded"), sum(reclass$screen_reclassification == "untestable"),
    n_dual_confirmed, nrow(wave2_plan), nrow(wave2_features), nrow(wave2_rows)
  )),
  pair_classes[, .(metric = paste0("pair_class_", classification), value = N)]
), fill = TRUE)
write_tsv(summary, file.path(t11a_dir, "t11a_summary.tsv"))

decision <- data.table::data.table(
  gate = "T11A_P1A_wave1_direct_multisignal_coloc",
  status = if (all(validation$passed)) "PASS" else "FAIL",
  decision = if (!all(validation$passed)) "HOLD_T11A" else if (continue_wave2) "GO_P1A_WAVE2" else "STOP_AFTER_P1A_WAVE1",
  stopping_rule = "at_least_1_feature_confirmed_in_both_diseases_same_direction_OR_at_least_4_of_22_trait_feature_rows_confirmed",
  confirmed_trait_feature_rows = n_confirmed_rows,
  total_trait_feature_rows = nrow(reclass),
  dual_disease_same_direction_confirmed_features = n_dual_confirmed,
  total_feature_units = nrow(feature_results),
  rationale = paste0(n_confirmed_rows, "/", nrow(reclass), " rows confirmed; ", n_dual_confirmed, "/",
                     nrow(feature_results), " features confirmed in both diseases with the same direction."),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_tsv(decision, file.path(t11a_dir, "t11a_decision.tsv"))

if (!all(validation$passed)) {
  failed <- validation[passed == FALSE, check]
  stop("T11A validation failed: ", paste(failed, collapse = ", "), call. = FALSE)
}
log_message("T11A summary complete: ", decision$decision,
            "; ", n_confirmed_rows, "/22 confirmed rows; ", n_dual_confirmed, "/11 dual-confirmed features")
