#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)

out_dir <- file.path(project_root, cfg$project$outcome_dir, "t11")
ensure_dir(out_dir)

paths <- list(
  reclass = file.path(project_root, cfg$project$outcome_dir, "t10", "regional_screen_reclassification.tsv.gz"),
  opposite = file.path(project_root, cfg$project$outcome_dir, "t10", "opposite_direction_signal_audit.tsv"),
  qtl_non_hla = file.path(project_root, cfg$project$outcome_dir, "t10", "non_hla", "qtl_feature_qc.tsv.gz"),
  qtl_hla = file.path(project_root, cfg$project$outcome_dir, "t10", "hla", "qtl_feature_qc.tsv.gz"),
  target_rows = file.path(project_root, cfg$project$outcome_dir, "targeted_susie_core_lineages", "target_rows.tsv.gz"),
  target_regions = file.path(project_root, cfg$project$outcome_dir, "targeted_susie_core_lineages", "target_regions.tsv"),
  target_files = file.path(project_root, cfg$project$outcome_dir, "targeted_susie_core_lineages", "target_files.tsv"),
  all_regions = file.path(project_root, cfg$project$outcome_dir, "targeted_susie", "target_regions.tsv"),
  all_files = file.path(project_root, cfg$project$outcome_dir, "targeted_susie", "target_files.tsv"),
  signal_catalog = file.path(project_root, cfg$project$outcome_dir, "t06", "susie_signal_catalog.tsv"),
  p0_plan = file.path(project_root, cfg$project$outcome_dir, "t08", "p0_core_non_hla_acquisition_plan.tsv"),
  t09_inventory = file.path(project_root, cfg$project$outcome_dir, "t09", "t09_file_inventory.tsv"),
  t09_history = file.path(project_root, cfg$project$outcome_dir, "t09", "t09_acquisition_history.tsv")
)
missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing T11 input(s): ", paste(missing, collapse = ", "), call. = FALSE)

reclass <- data.table::fread(paths$reclass, showProgress = FALSE)
opposite <- data.table::fread(paths$opposite, showProgress = FALSE)
qtl_qc <- data.table::rbindlist(list(
  data.table::fread(paths$qtl_non_hla, showProgress = FALSE),
  data.table::fread(paths$qtl_hla, showProgress = FALSE)
), fill = TRUE)

rate_row <- function(scope, x) {
  n_all <- nrow(x)
  n_testable <- sum(x$n_direct_pairs > 0L)
  data.table::data.table(
    scope = scope,
    n_screen_rows = n_all,
    n_testable = n_testable,
    testability_rate = n_testable / n_all,
    n_confirmed = sum(x$screen_reclassification == "confirmed"),
    confirmation_rate_all = sum(x$screen_reclassification == "confirmed") / n_all,
    confirmation_rate_testable = if (n_testable) sum(x$screen_reclassification == "confirmed") / n_testable else NA_real_,
    n_conflicting = sum(x$screen_reclassification == "conflicting"),
    conflict_rate_testable = if (n_testable) sum(x$screen_reclassification == "conflicting") / n_testable else NA_real_,
    n_downgraded = sum(x$screen_reclassification == "downgraded"),
    downgrade_rate_testable = if (n_testable) sum(x$screen_reclassification == "downgraded") / n_testable else NA_real_,
    n_untestable = sum(x$screen_reclassification == "untestable"),
    untestable_rate = sum(x$screen_reclassification == "untestable") / n_all,
    n_decisive = sum(x$screen_reclassification %in% c("confirmed", "conflicting")),
    decisive_rate_testable = if (n_testable) sum(x$screen_reclassification %in% c("confirmed", "conflicting")) / n_testable else NA_real_
  )
}

p0_rates <- data.table::rbindlist(list(
  rate_row("all_P0", reclass),
  rate_row("non_hla", reclass[arm == "non_hla"]),
  rate_row("hla", reclass[arm == "hla"]),
  rate_row("non_hla_Graves", reclass[arm == "non_hla" & trait_id == cfg$primary_comparison$graves]),
  rate_row("non_hla_AIHT", reclass[arm == "non_hla" & trait_id == cfg$primary_comparison$autoimmune_hypothyroidism])
), fill = TRUE)
write_tsv(p0_rates, file.path(out_dir, "p0_reclassification_rates.tsv"))

qtl_availability <- qtl_qc[, .(
  n_features = .N,
  n_qtl_ready = sum(feature_status == "ready"),
  qtl_feature_availability_rate = sum(feature_status == "ready") / .N
), by = .(scope = ifelse(is_hla, "hla", "non_hla"))]
qtl_availability <- data.table::rbindlist(list(
  data.table::data.table(
    scope = "all_P0", n_features = nrow(qtl_qc),
    n_qtl_ready = sum(qtl_qc$feature_status == "ready"),
    qtl_feature_availability_rate = mean(qtl_qc$feature_status == "ready")
  ), qtl_availability
))
write_tsv(qtl_availability, file.path(out_dir, "p0_qtl_availability.tsv"))

information_gain <- data.table::data.table(
  information_class = c(
    "direct_H4_confirmation", "direct_H3_conflict", "downgraded_regional_screen",
    "explicitly_untestable", "opposite_direction_confirmed", "opposite_direction_prevented_from_overclaim"
  ),
  n = c(
    sum(reclass$screen_reclassification == "confirmed"),
    sum(reclass$screen_reclassification == "conflicting"),
    sum(reclass$screen_reclassification == "downgraded"),
    sum(reclass$screen_reclassification == "untestable"),
    sum(opposite$opposite_direction_audit == "confirmed"),
    sum(opposite$opposite_direction_audit != "confirmed")
  ),
  denominator = c(nrow(reclass), nrow(reclass), nrow(reclass), nrow(reclass), nrow(opposite), nrow(opposite))
)
information_gain[, fraction := n / denominator]
write_tsv(information_gain, file.path(out_dir, "p0_information_gain.tsv"))

# Identify rows with an overlapping official GWAS SuSiE credible-set region.
target_rows <- data.table::fread(paths$target_rows, showProgress = FALSE)
signal_catalog <- data.table::fread(paths$signal_catalog, showProgress = FALSE)
target_rows[, chrom_label := sub("^chr", "", as.character(chrom))]
gwas_regions <- unique(signal_catalog[, .(
  trait_id, chrom_label, gwas_region_start = region_start, gwas_region_end = region_end
)])
coverage_map <- merge(target_rows, gwas_regions, by = c("trait_id", "chrom_label"), allow.cartesian = TRUE)
coverage_map <- coverage_map[gwas_region_start <= end & gwas_region_end >= start]
row_keys <- c("trait_id", "locus_id", "modality", "level", "cell_type", "feature_id", "priority")
covered <- unique(coverage_map[, ..row_keys])
target_rows[, row_key := do.call(paste, c(.SD, sep = "||")), .SDcols = row_keys]
covered[, row_key := do.call(paste, c(.SD, sep = "||")), .SDcols = row_keys]
target_rows[, gwas_susie_covered := row_key %in% covered$row_key]

feature_keys <- c("locus_id", "modality", "level", "cell_type", "feature_id")
p1 <- target_rows[priority == "P1_strong" & is_hla == FALSE]
p1_features <- p1[, .(
  n_traits = data.table::uniqueN(trait_id),
  n_traits_with_gwas_susie = data.table::uniqueN(trait_id[gwas_susie_covered]),
  both_regional_h4_strong = all(PP.H4.abf >= as.numeric(cfg$coloc$pp_h4_strong)),
  both_qtl_pure = all(as.logical(qtl_has_pure_cs)),
  same_regional_direction = data.table::uniqueN(effect_direction) == 1L,
  regional_direction = paste(sort(unique(effect_direction)), collapse = ";"),
  min_regional_PP.H4 = min(PP.H4.abf),
  max_regional_PP.H4 = max(PP.H4.abf),
  source_path = unique(source_path)[1],
  chrom = unique(chrom)[1],
  start = min(start), end = max(end)
), by = feature_keys]

# P1A is deliberately narrow: dual-trait, both GWAS-finemapped, pure QTL CS,
# concordant regional direction, core lineage, non-HLA. Regional direction is
# still only a screening property; direct T10-style coloc remains mandatory.
p1a_features <- p1_features[
  n_traits == 2L & n_traits_with_gwas_susie == 2L & both_regional_h4_strong &
    both_qtl_pure & same_regional_direction
]
data.table::setorder(p1a_features, -min_regional_PP.H4, locus_id, modality, level, cell_type, feature_id)
p1a_features[, p1a_rank := seq_len(.N)]
p1a_rows <- merge(p1, p1a_features[, c(feature_keys, "p1a_rank"), with = FALSE], by = feature_keys)
data.table::setorder(p1a_rows, p1a_rank, trait_id)
write_tsv(p1a_features, file.path(out_dir, "p1a_candidate_features.tsv"))
write_tsv(p1a_rows, file.path(out_dir, "p1a_candidate_rows.tsv.gz"))

p1a_regions <- unique(p1a_rows[, .(
  source_path, modality, level, cell_type, chrom = manifest_chrom, locus_id,
  start, end, tabix_region
)])
p1a_regions[, query_bp := as.numeric(end) - as.numeric(start) + 1]
p1a_region_features <- p1a_rows[, .(
  n_features = data.table::uniqueN(feature_id),
  features = paste(sort(unique(feature_id)), collapse = ";"),
  min_p1a_rank = min(p1a_rank)
), by = .(source_path, locus_id)]
p1a_regions <- merge(p1a_regions, p1a_region_features, by = c("source_path", "locus_id"), sort = FALSE)
data.table::setorder(p1a_regions, min_p1a_rank, modality, level, cell_type, chrom)
write_tsv(p1a_regions, file.path(out_dir, "p1a_acquisition_regions.tsv"))

p1a_plan <- p1a_rows[, .(
  modality = unique(modality)[1], level = unique(level)[1], cell_type = unique(cell_type)[1],
  chrom = unique(manifest_chrom)[1], data_url = unique(data_url)[1], index_url = unique(index_url)[1],
  size_bytes = unique(as.numeric(size_bytes))[1], md5 = unique(md5)[1],
  index_size_bytes = unique(as.numeric(index_size_bytes))[1], index_md5 = unique(index_md5)[1]
), by = source_path]
p1a_region_plan <- p1a_regions[, .(
  n_loci = data.table::uniqueN(locus_id),
  tabix_regions = paste(unique(tabix_region[order(start, end)]), collapse = ";")
), by = source_path]
p1a_plan <- merge(p1a_plan, p1a_region_plan, by = "source_path")
p1a_plan[, `:=`(
  stage = "through_P1_strong", max_priority_rank = 1L,
  selection_reasons = "P1A_dual_trait_both_GWAS_SuSiE_pure_QTL_CS_concordant_core_nonHLA",
  recommended_mode = "remote_extract"
)]
data.table::setcolorder(p1a_plan, c(
  "stage", "max_priority_rank", "modality", "level", "cell_type", "chrom",
  "source_path", "data_url", "index_url", "size_bytes", "md5",
  "index_size_bytes", "index_md5", "n_loci", "tabix_regions",
  "selection_reasons", "recommended_mode"
))
data.table::setorder(p1a_plan, modality, level, cell_type, chrom)
write_tsv(p1a_plan, file.path(out_dir, "p1a_acquisition_plan.tsv"))

# Wave 1 includes the source files containing the top ten ranked features.
# Because extraction occurs per source-locus interval, one additional feature
# sharing those files is obtained at no extra query cost.
p1a_wave1_sources <- unique(p1a_features[p1a_rank <= 10L]$source_path)
p1a_wave1_features <- p1a_features[source_path %in% p1a_wave1_sources]
p1a_wave1_rows <- p1a_rows[source_path %in% p1a_wave1_sources]
p1a_wave1_regions <- p1a_regions[source_path %in% p1a_wave1_sources]
p1a_wave1_plan <- p1a_plan[source_path %in% p1a_wave1_sources]
write_tsv(p1a_wave1_features, file.path(out_dir, "p1a_wave1_candidate_features.tsv"))
write_tsv(p1a_wave1_rows, file.path(out_dir, "p1a_wave1_candidate_rows.tsv.gz"))
write_tsv(p1a_wave1_regions, file.path(out_dir, "p1a_wave1_acquisition_regions.tsv"))
write_tsv(p1a_wave1_plan, file.path(out_dir, "p1a_wave1_acquisition_plan.tsv"))

parse_query_bp <- function(x) {
  regions <- strsplit(x, ";", fixed = TRUE)[[1]]
  sum(vapply(regions, function(region) {
    pieces <- regmatches(region, regexec(":([0-9]+)-([0-9]+)$", region))[[1]]
    as.numeric(pieces[[3]]) - as.numeric(pieces[[2]]) + 1
  }, numeric(1)))
}
p0_plan <- data.table::fread(paths$p0_plan, showProgress = FALSE)
p0_plan[, query_bp := vapply(tabix_regions, parse_query_bp, numeric(1))]
t09_inventory <- data.table::fread(paths$t09_inventory, showProgress = FALSE)[arm == "non_hla"]
t09_inventory[, observed_size_bytes := as.numeric(observed_size_bytes)]
empirical_bytes_per_source_bp <- sum(t09_inventory$observed_size_bytes) / sum(p0_plan$query_bp)
t09_history <- data.table::fread(paths$t09_history, showProgress = FALSE)
p0_completed <- t09_history[arm == "non_hla" & status == "completed"]
p0_elapsed_hours <- as.numeric(difftime(max(p0_completed$timestamp), min(p0_completed$timestamp), units = "hours"))

core_regions <- data.table::fread(paths$target_regions, showProgress = FALSE)
core_regions[, query_bp := as.numeric(end) - as.numeric(start) + 1]
core_files <- data.table::fread(paths$target_files, showProgress = FALSE)
all_regions <- data.table::fread(paths$all_regions, showProgress = FALSE)
all_regions[, query_bp := as.numeric(end) - as.numeric(start) + 1]
all_files <- data.table::fread(paths$all_files, showProgress = FALSE)

scenario_row <- function(scenario, regions, files, selected_rows, features, loci, decision) {
  query_bp <- sum(regions$query_bp)
  serial_hours <- p0_elapsed_hours * query_bp / sum(p0_plan$query_bp)
  data.table::data.table(
    scenario = scenario, n_source_files = data.table::uniqueN(files$source_path),
    n_source_locus_regions = nrow(regions), n_selected_rows = selected_rows,
    n_features = features, n_loci = loci, query_source_bp = query_bp,
    full_source_gib = sum(unique(files[, .(source_path, size_bytes)])$size_bytes) / 1024^3,
    estimated_targeted_gib = query_bp * empirical_bytes_per_source_bp / 1024^3,
    estimated_acquisition_hours_4way = serial_hours / 4,
    estimated_acquisition_hours_staged = serial_hours,
    estimated_local_analysis_minutes = if (is.na(selected_rows)) NA_real_ else max(1, 2 * selected_rows / nrow(reclass)),
    decision = decision
  )
}

p1a_files <- unique(p1a_rows[, .(source_path, size_bytes = as.numeric(size_bytes))])
storage_runtime <- data.table::rbindlist(list(
  scenario_row(
    "P1A_wave1", p1a_wave1_regions, unique(p1a_wave1_rows[, .(source_path, size_bytes = as.numeric(size_bytes))]),
    nrow(p1a_wave1_rows), nrow(p1a_wave1_features), data.table::uniqueN(p1a_wave1_features$locus_id), "GO_plan_only"
  ),
  scenario_row(
    "P1A_dual_trait_rescue", p1a_regions, p1a_files, nrow(p1a_rows),
    nrow(p1a_features), data.table::uniqueN(p1a_features$locus_id), "HOLD_after_wave1"
  ),
  scenario_row(
    "full_core_incremental_P1", core_regions[priority == "P1_strong"], core_files[priority == "P1_strong"],
    nrow(target_rows[priority == "P1_strong"]), data.table::uniqueN(target_rows[priority == "P1_strong"]$feature_id),
    data.table::uniqueN(target_rows[priority == "P1_strong"]$locus_id), "NO_GO"
  ),
  scenario_row(
    "full_core_incremental_P2", core_regions[priority == "P2_multisignal"], core_files[priority == "P2_multisignal"],
    nrow(target_rows[priority == "P2_multisignal"]), data.table::uniqueN(target_rows[priority == "P2_multisignal"]$feature_id),
    data.table::uniqueN(target_rows[priority == "P2_multisignal"]$locus_id), "NO_GO"
  ),
  scenario_row(
    "all_cell_incremental_P1", all_regions[priority == "P1_strong"], all_files[priority == "P1_strong"],
    NA_integer_, NA_integer_, data.table::uniqueN(all_regions[priority == "P1_strong"]$locus_id), "NO_GO"
  ),
  scenario_row(
    "all_cell_incremental_P2", all_regions[priority == "P2_multisignal"], all_files[priority == "P2_multisignal"],
    NA_integer_, NA_integer_, data.table::uniqueN(all_regions[priority == "P2_multisignal"]$locus_id), "NO_GO"
  )
), fill = TRUE)
write_tsv(storage_runtime, file.path(out_dir, "expansion_storage_runtime.tsv"))

expansion_coverage <- target_rows[priority %in% c("P1_strong", "P2_multisignal"), .(
  n_rows = .N, n_gwas_susie_covered = sum(gwas_susie_covered),
  gwas_susie_coverage_rate = mean(gwas_susie_covered),
  n_features = data.table::uniqueN(feature_id), n_loci = data.table::uniqueN(locus_id)
), by = .(priority, is_hla, trait_id)]
write_tsv(expansion_coverage, file.path(out_dir, "expansion_gwas_coverage.tsv"))

decision <- data.table::data.table(
  scope = c("P1A_wave1", "P1A_wave2", "full_core_P1", "core_P2", "all_cell_P1_P2", "HLA_repeat_with_same_inputs"),
  decision = c("GO_PLAN_ONLY", "HOLD", "NO_GO", "NO_GO", "NO_GO", "NO_GO"),
  rationale = c(
    paste0(nrow(p1a_wave1_features), " dual-trait concordant features across ", data.table::uniqueN(p1a_wave1_features$locus_id),
           " loci and ", nrow(p1a_wave1_plan), " source files form the highest-yield pilot; T11 freezes the plan but performs no acquisition."),
    paste0("The remaining ", nrow(p1a_features) - nrow(p1a_wave1_features),
           " P1A features proceed only if wave 1 yields at least one feature confirmed in both diseases or at least four confirmed trait-feature rows."),
    "Most rows are AIHT-only and the full expansion would dilute the shared-versus-divergent question while multiplying remote queries.",
    "P2 is lower-priority multisignal evidence and P0 showed low purity/overlap as the dominant failure modes.",
    "The all-cell plan exposes hundreds of source files and a 300-plus-GiB full-download footprint without a proportionate primary-hypothesis gain.",
    "No new Graves MHC SuSiE/HLA conditional evidence is available; repeating QTL extraction cannot repair the GWAS-side limitation."
  )
)
write_tsv(decision, file.path(out_dir, "t11_expansion_decision.tsv"))

stopping_rule <- data.table::data.table(
  action = c("CONTINUE_P1A_WAVE2", "STOP_AFTER_P1A_WAVE1"),
  rule = c(
    "at_least_1_feature_confirmed_in_both_diseases_same_direction_OR_at_least_4_of_22_trait_feature_rows_confirmed",
    "zero_dual_disease_confirmed_features_AND_fewer_than_4_of_22_trait_feature_rows_confirmed"
  ),
  interpretation = c(
    "Wave 1 demonstrates signal-level value for the shared-concordant axis.",
    "Additional QTL extraction is unlikely to repair the causal-signal bottleneck."
  )
)
write_tsv(stopping_rule, file.path(out_dir, "p1a_wave1_stopping_rule.tsv"))

expected_header <- c(
  "stage", "max_priority_rank", "modality", "level", "cell_type", "chrom", "source_path",
  "data_url", "index_url", "size_bytes", "md5", "index_size_bytes", "index_md5", "n_loci",
  "tabix_regions", "selection_reasons", "recommended_mode"
)
validation <- data.table::data.table(
  check = c(
    "p0_screen_total_333", "p0_rates_bounded", "p1a_has_both_traits_per_feature",
    "p1a_all_gwas_susie_covered", "p1a_all_non_hla_core_P1", "p1a_all_qtl_pure",
    "p1a_regional_directions_concordant", "p1a_plan_header_exact", "p1a_wave1_plan_header_exact", "p1a_remote_extract_only",
    "full_P1_P2_and_all_cell_are_no_go", "no_acquisition_performed_in_T11"
  ),
  passed = c(
    nrow(reclass) == 333L,
    all(unlist(p0_rates[, .(testability_rate, confirmation_rate_all, confirmation_rate_testable,
                            conflict_rate_testable, downgrade_rate_testable, untestable_rate,
                            decisive_rate_testable)])[is.finite(unlist(p0_rates[, .(
                              testability_rate, confirmation_rate_all, confirmation_rate_testable,
                              conflict_rate_testable, downgrade_rate_testable, untestable_rate,
                              decisive_rate_testable)]))] >= 0 &
          unlist(p0_rates[, .(testability_rate, confirmation_rate_all, confirmation_rate_testable,
                              conflict_rate_testable, downgrade_rate_testable, untestable_rate,
                              decisive_rate_testable)])[is.finite(unlist(p0_rates[, .(
                                testability_rate, confirmation_rate_all, confirmation_rate_testable,
                                conflict_rate_testable, downgrade_rate_testable, untestable_rate,
                                decisive_rate_testable)]))] <= 1),
    all(p1a_rows[, data.table::uniqueN(trait_id), by = feature_keys]$V1 == 2L),
    all(p1a_rows$gwas_susie_covered),
    all(p1a_rows$priority == "P1_strong" & p1a_rows$is_hla == FALSE &
          paste(p1a_rows$level, p1a_rows$cell_type, sep = ".") %in% c("l1.CD4_T", "l2.B_naive", "l2.Treg", "l2.CD14_Mono")),
    all(as.logical(p1a_rows$qtl_has_pure_cs)),
    all(p1a_features$same_regional_direction),
    identical(names(p1a_plan), expected_header),
    identical(names(p1a_wave1_plan), expected_header),
    all(p1a_plan$recommended_mode == "remote_extract"),
    all(decision[scope %in% c("full_core_P1", "core_P2", "all_cell_P1_P2")]$decision == "NO_GO"),
    !dir.exists(file.path(project_root, cfg$project$outcome_dir, "targeted_susie_t11"))
  )
)
write_tsv(validation, file.path(out_dir, "t11_validation.tsv"))

gate <- data.table::data.table(
  gate = "T11_P1_P2_expansion_decision",
  status = if (all(validation$passed)) "PASS" else "FAIL",
  decision = if (all(validation$passed)) "GO_P1A_WAVE1_PLAN_ONLY_HOLD_WAVE2_NO_GO_FULL_P1_P2" else "HOLD_T11",
  rationale = paste0(
    "P0 direct testability ", sprintf("%.1f%%", 100 * p0_rates[scope == "all_P0"]$testability_rate),
    ", confirmation among testable ", sprintf("%.1f%%", 100 * p0_rates[scope == "all_P0"]$confirmation_rate_testable),
    "; wave 1 retains ", nrow(p1a_wave1_features), " dual-trait pure-CS features in ", nrow(p1a_wave1_plan),
    " source files, wave 2 is conditional, and full P1/P2/all-cell plans are no-go."
  ),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_tsv(gate, file.path(out_dir, "t11_decision.tsv"))
if (!all(validation$passed)) stop("T11 validation failed.", call. = FALSE)
log_message("T11 complete: GO P1A wave 1 plan only; HOLD wave 2; NO-GO full P1/P2.")
