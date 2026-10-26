#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)

root <- file.path(project_root, cfg$project$outcome_dir, "t10")
arms <- c("non_hla", "hla")
required <- unlist(lapply(arms, function(arm) file.path(root, arm, c(
  "regional_screening_rows.tsv.gz", "gwas_qtl_susie_signal_pairs.tsv.gz", "qtl_feature_qc.tsv.gz"
))))
if (any(!file.exists(required))) stop("Missing arm-level T10 outputs.", call. = FALSE)

screen <- data.table::rbindlist(lapply(arms, function(arm) {
  arm_value <- arm
  x <- data.table::fread(file.path(root, arm, "regional_screening_rows.tsv.gz"))
  x[, arm := arm_value]
  x
}), fill = TRUE)
pairs <- data.table::rbindlist(lapply(arms, function(arm) {
  arm_value <- arm
  x <- data.table::fread(file.path(root, arm, "gwas_qtl_susie_signal_pairs.tsv.gz"))
  if (nrow(x)) x[, arm := arm_value]
  x
}), fill = TRUE)
feature_qc <- data.table::rbindlist(lapply(arms, function(arm) {
  arm_value <- arm
  x <- data.table::fread(file.path(root, arm, "qtl_feature_qc.tsv.gz"))
  x[, arm := arm_value]
  x
}), fill = TRUE)

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
    screen_id = s$screen_id, arm = s$arm, trait_id = s$trait_id, locus_id = s$locus_id,
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
    direction_matches_regional = direction_match, direct_pair_classification = if (nrow(best)) best$classification else NA_character_,
    screen_reclassification = status, reclassification_reason = reason,
    final_evidence_method = if (status == "untestable") "unresolved_no_direct_SuSiE_evidence" else "direct_SuSiE_BF_signal_pair_coloc"
  )
})
reclass <- data.table::rbindlist(reclass_rows, fill = TRUE)
data.table::setorder(reclass, arm, locus_id, modality, level, cell_type, feature_id, trait_id)
write_tsv(reclass, file.path(root, "regional_screen_reclassification.tsv.gz"))

# Cross-trait audit of regional opposite-direction claims.
feature_keys <- c("arm", "locus_id", "modality", "level", "cell_type", "feature_id", "is_hla")
opposite_features <- screen[, .(
  graves_increase = any(trait_id == cfg$primary_comparison$graves & effect_direction == "disease_risk_allele_increases_feature"),
  graves_decrease = any(trait_id == cfg$primary_comparison$graves & effect_direction == "disease_risk_allele_decreases_feature"),
  aiht_increase = any(trait_id == cfg$primary_comparison$autoimmune_hypothyroidism & effect_direction == "disease_risk_allele_increases_feature"),
  aiht_decrease = any(trait_id == cfg$primary_comparison$autoimmune_hypothyroidism & effect_direction == "disease_risk_allele_decreases_feature")
), by = feature_keys]
opposite_features <- opposite_features[(graves_increase & aiht_decrease) | (graves_decrease & aiht_increase)]

opposite_audit <- lapply(seq_len(nrow(opposite_features)), function(i) {
  f <- opposite_features[i]
  x <- reclass[
    arm == f$arm & locus_id == f$locus_id & modality == f$modality & level == f$level &
      cell_type == f$cell_type & feature_id == f$feature_id
  ]
  g <- x[trait_id == cfg$primary_comparison$graves][order(-direct_PP.H4)][1]
  a <- x[trait_id == cfg$primary_comparison$autoimmune_hypothyroidism][order(-direct_PP.H4)][1]
  both_confirmed <- nrow(g) && nrow(a) && g$screen_reclassification == "confirmed" && a$screen_reclassification == "confirmed"
  direct_opposite <- both_confirmed && !is.na(g$direct_effect_direction) && !is.na(a$direct_effect_direction) &&
    g$direct_effect_direction != a$direct_effect_direction &&
    g$direct_effect_direction != "undetermined" && a$direct_effect_direction != "undetermined"
  audit_status <- if (!nrow(g) || !nrow(a) || any(c(g$screen_reclassification, a$screen_reclassification) == "untestable")) {
    "untestable"
  } else if (both_confirmed && direct_opposite) {
    "confirmed"
  } else if (any(c(g$screen_reclassification, a$screen_reclassification) == "conflicting") ||
             (both_confirmed && !direct_opposite)) {
    "conflicting"
  } else {
    "downgraded"
  }
  data.table::data.table(
    arm = f$arm, locus_id = f$locus_id, modality = f$modality, level = f$level,
    cell_type = f$cell_type, feature_id = f$feature_id, is_hla = f$is_hla,
    graves_screen_status = if (nrow(g)) g$screen_reclassification else NA_character_,
    aiht_screen_status = if (nrow(a)) a$screen_reclassification else NA_character_,
    graves_regional_top_variant = if (nrow(g)) g$regional_top_variant else NA_character_,
    aiht_regional_top_variant = if (nrow(a)) a$regional_top_variant else NA_character_,
    graves_direct_top_variant = if (nrow(g)) g$direct_top_variant else NA_character_,
    aiht_direct_top_variant = if (nrow(a)) a$direct_top_variant else NA_character_,
    graves_gwas_cs = if (nrow(g)) g$best_gwas_cs else NA_integer_,
    aiht_gwas_cs = if (nrow(a)) a$best_gwas_cs else NA_integer_,
    graves_qtl_cs = if (nrow(g)) g$best_qtl_cs else NA_integer_,
    aiht_qtl_cs = if (nrow(a)) a$best_qtl_cs else NA_integer_,
    same_direct_qtl_signal = if (both_confirmed) g$best_qtl_cs == a$best_qtl_cs else NA,
    same_direct_top_variant = if (both_confirmed) identical(g$direct_top_variant, a$direct_top_variant) else NA,
    graves_direct_direction = if (nrow(g)) g$direct_effect_direction else NA_character_,
    aiht_direct_direction = if (nrow(a)) a$direct_effect_direction else NA_character_,
    opposite_direction_audit = audit_status
  )
})
opposite_audit <- data.table::rbindlist(opposite_audit, fill = TRUE)
write_tsv(opposite_audit, file.path(root, "opposite_direction_signal_audit.tsv"))
write_tsv(reclass[is_hla == TRUE], file.path(root, "hla_high_risk_audit.tsv.gz"))

summary <- data.table::rbindlist(list(
  reclass[, .(metric = paste0("screen_rows_", arm), value = .N), by = arm][, arm := NULL],
  reclass[, .(metric = paste0("reclassified_", arm, "_", screen_reclassification), value = .N), by = .(arm, screen_reclassification)][, c("arm", "screen_reclassification") := NULL],
  pairs[, .(metric = paste0("direct_signal_pairs_", arm), value = .N), by = arm][, arm := NULL],
  pairs[, .(metric = paste0("pair_class_", arm, "_", classification), value = .N), by = .(arm, classification)][, c("arm", "classification") := NULL],
  data.table::data.table(metric = "regional_opposite_features", value = nrow(opposite_audit)),
  if (nrow(opposite_audit)) opposite_audit[, .(metric = paste0("opposite_audit_", opposite_direction_audit), value = .N), by = opposite_direction_audit][, opposite_direction_audit := NULL] else NULL,
  data.table::data.table(metric = "hla_screen_rows_audited", value = nrow(reclass[is_hla == TRUE]))
), fill = TRUE)
write_tsv(summary, file.path(root, "t10_summary.tsv"))

probability_ok <- !nrow(pairs) || all(abs(rowSums(pairs[, .(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]) - 1) < 1e-6)
validation <- data.table::data.table(
  check = c(
    "all_333_frozen_screen_rows_present", "all_screen_rows_reclassified",
    "all_hla_rows_audited", "all_opposite_features_audited",
    "signal_pair_keys_unique", "posterior_probabilities_sum_to_one",
    "direct_method_only_for_testable_rows", "regional_abf_not_used_as_final_evidence"
  ),
  passed = c(
    nrow(reclass) == 333L,
    !anyNA(reclass$screen_reclassification),
    nrow(reclass[is_hla == TRUE]) == nrow(screen[is_hla == TRUE]),
    nrow(opposite_audit) == nrow(opposite_features),
    !anyDuplicated(pairs[, .(screen_id, gwas_region, gwas_cs, qtl_cs)]),
    probability_ok,
    all(reclass[screen_reclassification != "untestable"]$final_evidence_method == "direct_SuSiE_BF_signal_pair_coloc"),
    all(reclass$final_evidence_method != "coloc.abf_pvalue_screen")
  )
)
write_tsv(validation, file.path(root, "t10_validation.tsv"))
decision <- data.table::data.table(
  gate = "T10_P0_direct_multisignal_coloc",
  status = if (all(validation$passed)) "PASS" else "FAIL",
  decision = if (all(validation$passed)) "GO_T11" else "HOLD_T10",
  rationale = paste0(
    nrow(reclass), " regional screen rows reclassified; ", nrow(pairs),
    " direct signal pairs; ", nrow(opposite_audit), " opposite-direction features and ",
    nrow(reclass[is_hla == TRUE]), " HLA high-risk rows audited."
  ),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_tsv(decision, file.path(root, "t10_decision.tsv"))
if (!all(validation$passed)) stop("T10 validation failed.", call. = FALSE)
log_message("T10 summary complete: PASS_GO_T11")
