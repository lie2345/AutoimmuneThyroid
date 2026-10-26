#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)

out_root <- file.path(project_root, cfg$project$outcome_dir)
t11a_dir <- file.path(out_root, "t11a")
t11b_dir <- file.path(out_root, "t11b")
direct_dir <- file.path(t11b_dir, "direct_coloc")
dir.create(t11b_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  screen = file.path(direct_dir, "regional_screening_rows.tsv.gz"),
  pairs = file.path(direct_dir, "gwas_qtl_susie_signal_pairs.tsv.gz"),
  feature_qc = file.path(direct_dir, "qtl_feature_qc.tsv.gz"),
  wave2_plan = file.path(t11a_dir, "p1a_wave2_acquisition_plan.tsv"),
  wave2_features = file.path(t11a_dir, "p1a_wave2_candidate_features.tsv"),
  wave2_rows = file.path(t11a_dir, "p1a_wave2_candidate_rows.tsv.gz"),
  verification = file.path(t11b_dir, "p1a_wave2_verification.tsv"),
  preflight = file.path(t11b_dir, "remote_tabix_preflight.tsv"),
  acquisition_log = file.path(out_root, "targeted_susie_t11b", "p1a_wave2", "acquisition_log.tsv"),
  wave1_reclass = file.path(t11a_dir, "wave1_screen_reclassification.tsv.gz"),
  wave1_features = file.path(t11a_dir, "wave1_feature_results.tsv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing T11B inputs: ", paste(missing, collapse = ", "), call. = FALSE)

screen <- data.table::fread(paths$screen, showProgress = FALSE)
pairs <- data.table::fread(paths$pairs, showProgress = FALSE)
feature_qc <- data.table::fread(paths$feature_qc, showProgress = FALSE)
wave2_plan <- data.table::fread(paths$wave2_plan, showProgress = FALSE)
wave2_features_expected <- data.table::fread(paths$wave2_features, showProgress = FALSE)
wave2_rows_expected <- data.table::fread(paths$wave2_rows, showProgress = FALSE)
verification <- data.table::fread(paths$verification, showProgress = FALSE)
preflight <- data.table::fread(paths$preflight, showProgress = FALSE)
acquisition_log <- data.table::fread(paths$acquisition_log, showProgress = FALSE)
wave1_reclass <- data.table::fread(paths$wave1_reclass, showProgress = FALSE)
wave1_features <- data.table::fread(paths$wave1_features, showProgress = FALSE)

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

# Keep the row-level rule identical to T10 and T11A.
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
write_tsv(reclass, file.path(t11b_dir, "wave2_screen_reclassification.tsv.gz"))

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
write_tsv(feature_results, file.path(t11b_dir, "wave2_feature_results.tsv"))

inventory <- merge(
  wave2_plan[, .(source_path, modality, level, cell_type, chrom, tabix_regions)],
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
write_tsv(inventory, file.path(t11b_dir, "wave2_file_inventory.tsv"))

part_roots <- c(t11b_dir, file.path(out_root, "targeted_susie_t11b"))
part_files <- unique(unlist(lapply(part_roots, function(d) {
  if (!dir.exists(d)) return(character())
  list.files(d, pattern = "\\.part([.]|$)", recursive = TRUE, full.names = TRUE, all.files = TRUE)
})))
write_tsv(data.table::data.table(path = part_files), file.path(t11b_dir, "residual_part_files.tsv"))

wave1_reclass[, wave := "wave1"]
reclass[, wave := "wave2"]
p1a_reclass <- data.table::rbindlist(list(wave1_reclass, reclass), fill = TRUE)
data.table::setorder(p1a_reclass, locus_id, modality, level, cell_type, feature_id, trait_id)
write_tsv(p1a_reclass, file.path(t11b_dir, "p1a_screen_reclassification.tsv.gz"))

wave1_features[, wave := "wave1"]
feature_results[, wave := "wave2"]
p1a_features <- data.table::rbindlist(list(wave1_features, feature_results), fill = TRUE)
data.table::setorder(p1a_features, locus_id, modality, level, cell_type, feature_id)
write_tsv(p1a_features, file.path(t11b_dir, "p1a_feature_results.tsv"))

latest_log <- acquisition_log[order(timestamp), .SD[.N], by = source_path]
posterior_ok <- !nrow(pairs) || all(abs(rowSums(pairs[, .(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]) - 1) < 1e-6)
pair_key_unique <- !anyDuplicated(pairs[, .(screen_id, gwas_region, gwas_cs, qtl_cs)])
row_key <- c("trait_id", "locus_id", "modality", "level", "cell_type", "feature_id")
feature_key <- c("locus_id", "modality", "level", "cell_type", "feature_id")
validation <- data.table::data.table(
  check = c(
    "remote_preflight_9_of_9_ok", "acquisition_latest_9_completed_or_skipped", "verification_9_of_9_ok",
    "bgzip_integrity_9_of_9_ok", "local_checksums_present", "no_residual_part_files",
    "all_30_wave2_screen_rows_present", "all_15_wave2_features_present", "all_30_wave2_rows_reclassified",
    "wave2_row_keys_unique", "wave2_feature_keys_unique", "signal_pair_keys_unique",
    "posterior_probabilities_sum_to_one", "zero_allele_mismatches_excluded",
    "complete_p1a_52_rows", "complete_p1a_26_features", "complete_p1a_keys_unique"
  ),
  passed = c(
    nrow(preflight) == 9L && all(preflight$status == "ok"),
    nrow(latest_log) == 9L && all(latest_log$status %in% c("completed", "skipped_existing")),
    nrow(verification) == 9L && all(verification$status == "ok") && all(verification$index_ok),
    nrow(inventory) == 9L && all(inventory$bgzip_ok),
    nrow(inventory) == 9L && all(nzchar(inventory$data_md5)) && all(nzchar(inventory$index_md5)),
    length(part_files) == 0L,
    nrow(screen) == 30L,
    nrow(feature_results) == 15L && nrow(wave2_features_expected) == 15L,
    nrow(reclass) == 30L && !anyNA(reclass$screen_reclassification) && nrow(wave2_rows_expected) == 30L,
    !anyDuplicated(reclass[, ..row_key]),
    !anyDuplicated(feature_results[, ..feature_key]),
    pair_key_unique,
    posterior_ok,
    sum(pairs$allele_mismatches_excluded, na.rm = TRUE) == 0,
    nrow(p1a_reclass) == 52L,
    nrow(p1a_features) == 26L,
    !anyDuplicated(p1a_reclass[, ..row_key]) && !anyDuplicated(p1a_features[, ..feature_key])
  )
)
write_tsv(validation, file.path(t11b_dir, "t11b_validation.tsv"))

pair_classes <- pairs[, .N, by = classification][order(classification)]
summary <- data.table::rbindlist(list(
  data.table::data.table(metric = c(
    "preflight_files_ok", "acquired_files_completed_or_skipped", "verified_files_ok", "acquired_records",
    "acquired_bytes", "wave2_screen_rows", "wave2_direct_signal_pairs", "wave2_feature_units",
    "wave2_confirmed_rows", "wave2_conflicting_rows", "wave2_downgraded_rows", "wave2_untestable_rows",
    "wave2_dual_disease_same_direction_confirmed_features", "p1a_total_rows", "p1a_confirmed_rows",
    "p1a_conflicting_rows", "p1a_downgraded_rows", "p1a_untestable_rows", "p1a_total_features",
    "p1a_dual_disease_same_direction_confirmed_features"
  ), value = c(
    sum(preflight$status == "ok"), sum(latest_log$status %in% c("completed", "skipped_existing")),
    sum(verification$status == "ok"), sum(verification$n_records_in_planned_regions),
    sum(verification$observed_size_bytes), nrow(screen), nrow(pairs), nrow(feature_results),
    sum(reclass$screen_reclassification == "confirmed"), sum(reclass$screen_reclassification == "conflicting"),
    sum(reclass$screen_reclassification == "downgraded"), sum(reclass$screen_reclassification == "untestable"),
    sum(feature_results$feature_result == "confirmed_both_same_direction"), nrow(p1a_reclass),
    sum(p1a_reclass$screen_reclassification == "confirmed"), sum(p1a_reclass$screen_reclassification == "conflicting"),
    sum(p1a_reclass$screen_reclassification == "downgraded"), sum(p1a_reclass$screen_reclassification == "untestable"),
    nrow(p1a_features), sum(p1a_features$feature_result == "confirmed_both_same_direction")
  )),
  pair_classes[, .(metric = paste0("wave2_pair_class_", classification), value = N)]
), fill = TRUE)
write_tsv(summary, file.path(t11b_dir, "t11b_summary.tsv"))

decision <- data.table::data.table(
  gate = "T11B_P1A_wave2_direct_multisignal_coloc",
  status = if (all(validation$passed)) "PASS" else "FAIL",
  decision = if (all(validation$passed)) "COMPLETE_P1A_STOP_EXPANSION" else "HOLD_T11B",
  rationale = if (all(validation$passed)) {
    paste0(
      nrow(reclass), " wave-2 rows and ", nrow(feature_results), " wave-2 features audited; complete P1A has ",
      nrow(p1a_reclass), " rows and ", nrow(p1a_features),
      " features. Frozen full P1/P2 expansion remains NO-GO under T11."
    )
  } else {
    paste0(sum(!validation$passed), " validation checks failed; T11B remains on hold.")
  },
  next_task = if (all(validation$passed)) "T12_global_cell_type_enrichment" else "resolve_T11B_validation_failures",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_tsv(decision, file.path(t11b_dir, "t11b_decision.tsv"))

if (!all(validation$passed)) {
  failed <- validation[passed == FALSE, check]
  stop("T11B validation failed: ", paste(failed, collapse = ", "), call. = FALSE)
}
log_message("T11B complete: PASS / COMPLETE_P1A_STOP_EXPANSION")
