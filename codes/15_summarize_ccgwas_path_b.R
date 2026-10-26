#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "ggplot2"))
cfg <- read_analysis_config(project_root)
cc_cfg <- cfg$ccgwas
out_root <- project_path(project_root, cc_cfg$output_dir)
results_dir <- file.path(out_root, "results")
summary_dir <- file.path(out_root, "summary")
figure_dir <- file.path(project_root, cfg$project$figure_dir, "ccgwas")
ensure_dir(summary_dir)
ensure_dir(figure_dir)

parameters <- data.table::fread(file.path(results_dir, "scenario_parameters.tsv"))
if (!all(parameters$scenario %in% c("overlap_000pct", "overlap_050pct", "overlap_100pct"))) {
  stop("Expected overlap scenarios are missing", call. = FALSE)
}

scenario_summaries <- list()
candidate_rows <- list()
for (i in seq_len(nrow(parameters))) {
  scenario <- parameters$scenario[[i]]
  path <- file.path(results_dir, paste0(scenario, ".results.gz"))
  if (!file.exists(path)) stop("Missing CC-GWAS result: ", path, call. = FALSE)
  x <- data.table::fread(path)
  scenario_summaries[[i]] <- data.table::data.table(
    scenario = scenario,
    n_tested = nrow(x),
    n_ols_gws = sum(x$OLS_pval < 5e-8, na.rm = TRUE),
    n_formal_ccgwas = sum(x$CCGWAS_signif == 1L, na.rm = TRUE),
    min_ols_p = min(x$OLS_pval, na.rm = TRUE),
    min_exact_p = min(x$Exact_pval, na.rm = TRUE)
  )
  candidate_rows[[i]] <- x[OLS_pval < 5e-8 | CCGWAS_signif == 1L, .(
    scenario, SNP, CHR, BP, EA, NEA,
    OLS_beta, OLS_se, OLS_pval,
    Exact_beta, Exact_se, Exact_pval,
    CCGWAS_signif
  )]
  rm(x)
  invisible(gc())
}
scenario_summary <- merge(
  data.table::rbindlist(scenario_summaries),
  parameters[, .(
    scenario,
    overlap_fraction_of_smaller_control_set,
    N_overlap_A0B0,
    modelled_intercept,
    primary_analysis
  )],
  by = "scenario",
  all.x = TRUE
)

all_path <- file.path(results_dir, "overlap_000pct.results.ALL.gz")
all_primary <- data.table::fread(all_path)
valid_p <- all_primary$OLS_pval[is.finite(all_primary$OLS_pval) & all_primary$OLS_pval > 0 & all_primary$OLS_pval <= 1]
lambda_gc <- stats::median(stats::qchisq(1 - valid_p, df = 1), na.rm = TRUE) / stats::qchisq(0.5, df = 1)
primary_diagnostics <- data.table::data.table(
  metric = c(
    "n_tested", "lambda_gc_ols", "n_ols_gws", "n_exact_failed",
    "n_differential_tagging_failed", "n_prevalence_range_failed", "n_formal_ccgwas"
  ),
  value = c(
    nrow(all_primary),
    lambda_gc,
    sum(all_primary$OLS_pval < 5e-8, na.rm = TRUE),
    sum(all_primary$OLS_pval < 5e-8 & all_primary$Exact_pval >= 1e-4, na.rm = TRUE),
    sum(all_primary$OLS_pval < 5e-8 & all_primary$potential_tagging_stresstest == 1L, na.rm = TRUE),
    sum(
      all_primary$OLS_pval < 5e-8 &
        (all_primary$Exact_ll_pval >= 1e-4 |
           all_primary$Exact_lh_pval >= 1e-4 |
           all_primary$Exact_hl_pval >= 1e-4 |
           all_primary$Exact_hh_pval >= 1e-4),
      na.rm = TRUE
    ),
    sum(all_primary$CCGWAS_signif == 1L, na.rm = TRUE)
  )
)
primary_components <- all_primary[CCGWAS_signif == 1L, .(
  SNP, CHR, BP, EA, NEA,
  A1A0_beta, A1A0_se, A1A0_pval,
  B1B0_beta, B1B0_se, B1B0_pval,
  OLS_beta, OLS_se, OLS_pval,
  Exact_beta, Exact_se, Exact_pval,
  Exact_ll_pval, Exact_lh_pval, Exact_hl_pval, Exact_hh_pval,
  potential_tagging_stresstest, CCGWAS_signif
)]
rm(all_primary, valid_p)
invisible(gc())

candidate_long <- data.table::rbindlist(candidate_rows)
data.table::setorder(candidate_long, CHR, BP, scenario)
data.table::fwrite(
  candidate_long,
  file.path(summary_dir, "scenario_candidate_variants.tsv.gz"),
  sep = "\t", quote = FALSE, compress = "gzip"
)

primary <- candidate_long[scenario == "overlap_000pct" & CCGWAS_signif == 1L]
stability <- primary[, .(
  SNP, CHR, BP, EA, NEA,
  primary_OLS_beta = OLS_beta,
  primary_OLS_se = OLS_se,
  primary_OLS_pval = OLS_pval,
  primary_Exact_beta = Exact_beta,
  primary_Exact_pval = Exact_pval
)]

for (scenario in parameters$scenario) {
  suffix <- sub("^overlap_", "overlap_", scenario)
  requested_scenario <- scenario
  value <- candidate_long[scenario == requested_scenario, .(
    SNP,
    formal = CCGWAS_signif == 1L,
    OLS_beta,
    OLS_pval
  )]
  data.table::setnames(
    value,
    c("formal", "OLS_beta", "OLS_pval"),
    paste0(c("formal_", "OLS_beta_", "OLS_pval_"), suffix)
  )
  stability <- merge(stability, value, by = "SNP", all.x = TRUE)
}
formal_columns <- grep("^formal_overlap_", names(stability), value = TRUE)
beta_columns <- grep("^OLS_beta_overlap_", names(stability), value = TRUE)
stability[, formal_all_overlap_scenarios := rowSums(as.matrix(.SD) == 1, na.rm = TRUE) == length(formal_columns), .SDcols = formal_columns]
stability[, same_ols_direction_all_overlap_scenarios := apply(.SD, 1L, function(z) {
  z <- z[is.finite(z)]
  length(z) == length(beta_columns) && all(sign(z) == sign(z[[1L]]))
}), .SDcols = beta_columns]

prev_g <- as.numeric(cc_cfg$prevalence$graves$main)
prev_a <- as.numeric(cc_cfg$prevalence$aiht$main)
contamination_values <- as.numeric(unlist(cc_cfg$contamination_fraction_sensitivity))
contamination_grid <- data.table::CJ(
  graves_control_contamination = contamination_values,
  aiht_control_contamination = contamination_values
)
contamination_rows <- vector("list", nrow(contamination_grid))
for (i in seq_len(nrow(contamination_grid))) {
  grid <- contamination_grid[i]
  g_retention <- 1 - grid$graves_control_contamination
  a_retention <- 1 - grid$aiht_control_contamination
  temp <- data.table::copy(primary_components)
  temp[, `:=`(
    graves_control_contamination = grid$graves_control_contamination,
    aiht_control_contamination = grid$aiht_control_contamination,
    corrected_exact_beta = (1 - prev_g) * A1A0_beta / g_retention -
      (1 - prev_a) * B1B0_beta / a_retention,
    corrected_exact_se = sqrt(
      ((1 - prev_g) * A1A0_se / g_retention)^2 +
        ((1 - prev_a) * B1B0_se / a_retention)^2
    )
  )]
  temp[, corrected_exact_pval := 2 * stats::pnorm(-abs(corrected_exact_beta / corrected_exact_se))]
  temp[, same_direction_as_primary := sign(corrected_exact_beta) == sign(Exact_beta)]
  temp[, passes_exact_threshold := corrected_exact_pval < 1e-4]
  contamination_rows[[i]] <- temp[, .(
    SNP, CHR, BP, EA, NEA,
    graves_control_contamination,
    aiht_control_contamination,
    corrected_exact_beta,
    corrected_exact_se,
    corrected_exact_pval,
    same_direction_as_primary,
    passes_exact_threshold
  )]
}
contamination <- data.table::rbindlist(contamination_rows)
data.table::fwrite(
  contamination,
  file.path(summary_dir, "control_contamination_sensitivity.tsv.gz"),
  sep = "\t", quote = FALSE, compress = "gzip"
)
contamination_status <- contamination[, .(
  contamination_direction_stable = all(same_direction_as_primary),
  contamination_exact_threshold_stable = all(passes_exact_threshold),
  worst_contamination_pval = max(corrected_exact_pval),
  min_corrected_beta = min(corrected_exact_beta),
  max_corrected_beta = max(corrected_exact_beta)
), by = SNP]
stability <- merge(stability, contamination_status, by = "SNP", all.x = TRUE)
stability[, final_sensitivity_status := data.table::fcase(
  formal_all_overlap_scenarios &
    same_ols_direction_all_overlap_scenarios &
    contamination_direction_stable &
    contamination_exact_threshold_stable,
  "robust_candidate",
  same_ols_direction_all_overlap_scenarios & contamination_direction_stable,
  "unresolved_threshold_sensitive",
  default = "unresolved"
)]
data.table::setorder(stability, CHR, BP)
data.table::fwrite(stability, file.path(summary_dir, "signal_stability.tsv"), sep = "\t", quote = FALSE)

annotation <- data.table::fread(file.path(out_root, "input", "matched_variant_annotation.tsv.gz"))
stability <- merge(stability, annotation, by = c("SNP", "CHR", "BP", "EA", "NEA"), all.x = TRUE)
data.table::setorder(stability, CHR, BP)
stability[, new_locus := is.na(data.table::shift(BP)) |
  CHR != data.table::shift(CHR) |
  (BP - data.table::shift(BP) > 1000000L)]
stability[, provisional_locus := cumsum(new_locus)]
loci <- stability[, {
  lead_index <- which.min(primary_OLS_pval)
  .(
    CHR = CHR[[1L]],
    start = min(BP),
    end = max(BP),
    n_primary_significant_snps = .N,
    n_robust_snps = sum(final_sensitivity_status == "robust_candidate"),
    lead_SNP = SNP[[lead_index]],
    lead_EA = EA[[lead_index]],
    lead_NEA = NEA[[lead_index]],
    lead_OLS_beta = primary_OLS_beta[[lead_index]],
    lead_OLS_se = primary_OLS_se[[lead_index]],
    lead_OLS_pval = primary_OLS_pval[[lead_index]],
    lead_Exact_beta = primary_Exact_beta[[lead_index]],
    lead_Exact_pval = primary_Exact_pval[[lead_index]],
    lead_nearest_genes_graves = nearest_genes_graves[[lead_index]],
    lead_nearest_genes_aiht = nearest_genes_aiht[[lead_index]],
    lead_sensitivity_status = final_sensitivity_status[[lead_index]],
    locus_sensitivity_status = if (any(final_sensitivity_status == "robust_candidate")) "robust_candidate_locus" else "unresolved_locus"
  )
}, by = provisional_locus]
loci[, locus_id := sprintf("CC%02d_chr%d_%d_%d", provisional_locus, CHR, start, end)]
graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
graves_lead_effects <- data.table::fread(
  file.path(out_root, "input", paste0(graves_id, ".ccgwas.tsv.gz")),
  select = c("SNP", "OR", "SE", "P")
)[SNP %in% loci$lead_SNP]
aiht_lead_effects <- data.table::fread(
  file.path(out_root, "input", paste0(aiht_id, ".ccgwas.tsv.gz")),
  select = c("SNP", "OR", "SE", "P")
)[SNP %in% loci$lead_SNP]
data.table::setnames(
  graves_lead_effects,
  c("SNP", "OR", "SE", "P"),
  c("lead_SNP", "graves_case_control_OR", "graves_case_control_SE_logOR", "graves_case_control_P")
)
data.table::setnames(
  aiht_lead_effects,
  c("SNP", "OR", "SE", "P"),
  c("lead_SNP", "aiht_case_control_OR", "aiht_case_control_SE_logOR", "aiht_case_control_P")
)
loci <- merge(loci, graves_lead_effects, by = "lead_SNP", all.x = TRUE)
loci <- merge(loci, aiht_lead_effects, by = "lead_SNP", all.x = TRUE)
loci[, observed_case_control_sign_pattern := data.table::fcase(
  sign(log(graves_case_control_OR)) != sign(log(aiht_case_control_OR)),
  "opposite_observed_signs",
  sign(log(graves_case_control_OR)) == sign(log(aiht_case_control_OR)),
  "concordant_observed_signs",
  default = "unresolved_observed_signs"
)]
data.table::setorder(loci, provisional_locus)
data.table::setcolorder(loci, c("locus_id", setdiff(names(loci), "locus_id")))
data.table::fwrite(loci, file.path(summary_dir, "provisional_lead_loci.tsv"), sep = "\t", quote = FALSE)

analysis_summary <- data.table::data.table(
  metric = c(
    "primary_tested_snps",
    "primary_formal_significant_snps",
    "primary_provisional_physical_loci",
    "primary_snps_formal_in_all_overlap_scenarios",
    "primary_snps_direction_stable_under_contamination",
    "primary_robust_candidate_snps",
    "robust_candidate_loci",
    "lead_loci_opposite_observed_case_control_signs",
    "lead_loci_concordant_observed_case_control_signs",
    "lambda_gc_ols"
  ),
  value = c(
    scenario_summary[scenario == "overlap_000pct", n_tested],
    nrow(stability),
    nrow(loci),
    sum(stability$formal_all_overlap_scenarios),
    sum(stability$contamination_direction_stable),
    sum(stability$final_sensitivity_status == "robust_candidate"),
    sum(loci$locus_sensitivity_status == "robust_candidate_locus"),
    sum(loci$observed_case_control_sign_pattern == "opposite_observed_signs"),
    sum(loci$observed_case_control_sign_pattern == "concordant_observed_signs"),
    lambda_gc
  )
)
data.table::fwrite(scenario_summary, file.path(summary_dir, "scenario_summary.tsv"), sep = "\t", quote = FALSE)
data.table::fwrite(primary_diagnostics, file.path(summary_dir, "primary_diagnostics.tsv"), sep = "\t", quote = FALSE)
data.table::fwrite(analysis_summary, file.path(summary_dir, "analysis_summary.tsv"), sep = "\t", quote = FALSE)

lead_long <- merge(
  candidate_long,
  loci[, .(lead_SNP, locus_id)],
  by.x = "SNP", by.y = "lead_SNP",
  all = FALSE
)
lead_long[, scenario := factor(scenario, levels = parameters$scenario)]
plot <- ggplot2::ggplot(
  lead_long,
  ggplot2::aes(x = scenario, y = OLS_beta, group = locus_id)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.35) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = OLS_beta - 1.96 * OLS_se, ymax = OLS_beta + 1.96 * OLS_se),
    width = 0.12, colour = "grey35", linewidth = 0.35
  ) +
  ggplot2::geom_point(size = 1.6, colour = "#1F5A7A") +
  ggplot2::facet_wrap(~locus_id, scales = "free_y", ncol = 4) +
  ggplot2::labs(
    x = "Assumed overlap of the smaller control set",
    y = "CC-GWAS OLS beta (positive: Graves relative to AIHT)",
    title = "Lead-SNP stability across control-overlap assumptions"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 7),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  )
ggplot2::ggsave(
  file.path(figure_dir, "ccgwas_overlap_sensitivity_leads.png"),
  plot,
  width = 9,
  height = 8,
  dpi = 300,
  bg = "white"
)

log_message("Primary formal SNPs: ", nrow(stability))
log_message("Provisional physical loci: ", nrow(loci))
log_message("Robust candidate SNPs: ", sum(stability$final_sensitivity_status == "robust_candidate"))
log_message("OLS lambda GC: ", round(lambda_gc, 4))
