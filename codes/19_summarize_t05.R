#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "ggplot2"))
cfg <- read_analysis_config(project_root)
graves <- cfg$primary_comparison$graves
aiht <- cfg$primary_comparison$autoimmune_hypothyroidism
t05_root <- file.path(project_root, "outcome_data", "t05")
lava_root <- file.path(t05_root, "lava")
summary_dir <- file.path(t05_root, "summary")
figure_dir <- file.path(project_root, cfg$project$figure_dir, "t05")
ensure_dir(summary_dir)
ensure_dir(figure_dir)

ldsc <- data.table::fread(file.path(t05_root, "ldsc", "results", "t05_ldsc_summary.tsv"))
cross <- ldsc[analysis == "cross"]
if (nrow(cross) != 1L) stop("T05 cross-trait LDSC result is missing", call. = FALSE)

read_optional <- function(path) {
  if (file.exists(path)) data.table::fread(path) else data.table::data.table()
}

build_panel_status <- function(panel) {
  panel_dir <- file.path(lava_root, panel)
  provenance <- data.table::fread(file.path(panel_dir, "run_provenance.tsv"))
  if (panel == "discovery") {
    planned <- data.table::fread(file.path(lava_root, "input", "discovery_non_hla.lava"))[, .(
      locus_id = LOC, chrom = CHR, planned_start = START, planned_end = STOP
    )]
    planned[, `:=`(
      ccgwas_locus_id = NA_character_,
      lead_SNP = NA_character_,
      lead_nearest_gene = NA_character_,
      hla_adjacent_within_2mb = FALSE
    )]
  } else {
    planned <- data.table::fread(file.path(lava_root, "input", "targeted_ccgwas_locus_map.tsv"))[
      overlaps_hla == FALSE,
      .(
        locus_id = LOC, chrom = CHR, planned_start = START, planned_end = STOP,
        ccgwas_locus_id, lead_SNP, lead_nearest_gene, hla_adjacent_within_2mb
      )
    ]
  }
  univ <- data.table::fread(file.path(panel_dir, "univariate.tsv"))
  bivar <- read_optional(file.path(panel_dir, "bivariate.tsv"))
  failures <- read_optional(file.path(panel_dir, "failures.tsv"))
  threshold <- provenance$univariate_bonferroni_threshold[[1L]]

  g <- univ[phen == graves, .(
    locus_id,
    processed_start = start,
    processed_end = end,
    n_snps,
    n_pcs,
    graves_h2_obs = h2.obs,
    graves_h2_latent = h2.latent,
    graves_univ_p = p
  )]
  a <- univ[phen == aiht, .(
    locus_id,
    aiht_h2_obs = h2.obs,
    aiht_h2_latent = h2.latent,
    aiht_univ_p = p
  )]
  status <- merge(planned, g, by = "locus_id", all.x = TRUE)
  status <- merge(status, a, by = "locus_id", all.x = TRUE)
  status[, `:=`(
    panel = panel,
    univariate_threshold = threshold,
    graves_local_h2_significant = is.finite(graves_univ_p) & graves_univ_p < threshold,
    aiht_local_h2_significant = is.finite(aiht_univ_p) & aiht_univ_p < threshold
  )]

  if (nrow(bivar)) {
    bivar[, fdr_q := stats::p.adjust(p, method = "BH")]
    bivar_small <- bivar[, .(
      locus_id,
      local_rg = rho,
      local_rg_ci_lower = rho.lower,
      local_rg_ci_upper = rho.upper,
      local_r2 = r2,
      local_rg_p = p,
      local_rg_fdr_q = fdr_q
    )]
    status <- merge(status, bivar_small, by = "locus_id", all.x = TRUE)
  } else {
    status[, `:=`(
      local_rg = NA_real_, local_rg_ci_lower = NA_real_, local_rg_ci_upper = NA_real_,
      local_r2 = NA_real_, local_rg_p = NA_real_, local_rg_fdr_q = NA_real_
    )]
  }

  if (nrow(failures)) {
    failure_summary <- failures[, .(
      failure_stage = paste(unique(stage), collapse = ";"),
      failure_reason = paste(unique(reason), collapse = " | ")
    ), by = locus_id]
    status <- merge(status, failure_summary, by = "locus_id", all.x = TRUE)
  } else {
    status[, `:=`(failure_stage = NA_character_, failure_reason = NA_character_)]
  }

  status[, final_status := data.table::fcase(
    failure_stage == "preflight", "insufficient_data_build_mismatch_no_rsid",
    !is.na(failure_stage), "insufficient_data_processing_failure",
    !is.finite(graves_univ_p) | !is.finite(aiht_univ_p), "insufficient_local_h2_estimation",
    !graves_local_h2_significant | !aiht_local_h2_significant, "insufficient_joint_local_heritability",
    !is.finite(local_rg_p), "insufficient_bivariate_result",
    hla_adjacent_within_2mb, "hla_adjacent_tested_separately",
    local_rg_fdr_q < 0.05, "local_rg_evidence_fdr",
    default = "tested_no_local_rg_evidence"
  )]
  status[, interpretation := data.table::fcase(
    final_status == "tested_no_local_rg_evidence", "Both traits had local signal; local rg was tested but not significant after FDR",
    final_status == "local_rg_evidence_fdr", "Both traits had local signal and local rg passed panel-wise FDR",
    final_status == "hla_adjacent_tested_separately", "Statistically tested but interpreted separately because the window lies within 2 Mb of HLA",
    final_status == "insufficient_joint_local_heritability", "At least one trait did not pass the panel-wise local h2 gate; local rg was not tested",
    final_status == "insufficient_local_h2_estimation", "At least one trait lacked an estimable non-negative local variance; local rg was not tested",
    final_status == "insufficient_data_build_mismatch_no_rsid", "No eligible rsID list allowed safe GRCh38-to-GRCh37 reference use",
    default = "Local rg could not be interpreted"
  )]
  data.table::setorder(status, chrom, planned_start)

  summary <- status[, .N, by = final_status]
  summary[, panel := panel]
  summary[, planned_loci := nrow(status)]
  summary[, univariate_threshold := threshold]
  summary[, bivariate_tests := sum(is.finite(status$local_rg_p))]
  summary[, fdr_significant_all := sum(status$local_rg_fdr_q < 0.05, na.rm = TRUE)]
  summary[, fdr_significant_ordinary_non_hla := sum(
    status$local_rg_fdr_q < 0.05 & !status$hla_adjacent_within_2mb,
    na.rm = TRUE
  )]
  list(status = status, summary = summary, bivar = bivar)
}

discovery <- build_panel_status("discovery")
targeted <- build_panel_status("targeted")
all_status <- data.table::rbindlist(list(discovery$status, targeted$status), fill = TRUE)
panel_summary <- data.table::rbindlist(list(discovery$summary, targeted$summary), fill = TRUE)
write_tsv(all_status, file.path(summary_dir, "t05_locus_testability.tsv"))
write_tsv(panel_summary, file.path(summary_dir, "t05_panel_summary.tsv"))

targeted_status <- targeted$status[, .(
  ccgwas_locus_id,
  lead_SNP,
  lead_nearest_gene,
  hla_adjacent_within_2mb,
  graves_h2_obs,
  graves_univ_p,
  aiht_h2_obs,
  aiht_univ_p,
  local_rg,
  local_rg_ci_lower,
  local_rg_ci_upper,
  local_rg_p,
  local_rg_fdr_q,
  final_status,
  interpretation
)]
data.table::setorder(targeted_status, ccgwas_locus_id)
write_tsv(targeted_status, file.path(summary_dir, "t04b_loci_lava_followup.tsv"))

hla_excluded <- data.table::fread(file.path(lava_root, "hla_excluded_from_ordinary_lava.tsv"))
hla_summary <- data.table::data.table(
  excluded_locus_id = hla_excluded$LOC,
  chrom = hla_excluded$CHR,
  start = hla_excluded$START,
  end = hla_excluded$STOP,
  reason = "Overlaps prespecified HLA interval and was excluded from ordinary LAVA"
)
write_tsv(hla_summary, file.path(summary_dir, "t05_hla_separation.tsv"))

overall <- data.table::data.table(
  metric = c(
    "global_rg", "global_rg_se", "global_rg_ci_lower", "global_rg_ci_upper",
    "global_rg_p", "cross_trait_intercept", "lava_overlap_proxy",
    "discovery_planned_non_hla_loci", "discovery_preflight_excluded",
    "discovery_bivariate_tested", "discovery_local_rg_fdr_significant",
    "targeted_planned_loci", "targeted_bivariate_tested",
    "targeted_local_rg_fdr_significant_all", "targeted_local_rg_fdr_significant_ordinary_non_hla",
    "hla_regions_excluded_from_ordinary_lava"
  ),
  value = c(
    cross$rg, cross$se, cross$rg_ci_lower, cross$rg_ci_upper,
    cross$p, cross$gcov_int,
    data.table::fread(file.path(t05_root, "lava", "input", "sample.overlap.txt"))[[2L]][[2L]],
    nrow(discovery$status),
    sum(discovery$status$final_status == "insufficient_data_build_mismatch_no_rsid"),
    sum(is.finite(discovery$status$local_rg_p)),
    sum(discovery$status$local_rg_fdr_q < 0.05, na.rm = TRUE),
    nrow(targeted$status),
    sum(is.finite(targeted$status$local_rg_p)),
    sum(targeted$status$local_rg_fdr_q < 0.05, na.rm = TRUE),
    sum(targeted$status$local_rg_fdr_q < 0.05 & !targeted$status$hla_adjacent_within_2mb, na.rm = TRUE),
    nrow(hla_summary)
  )
)
write_tsv(overall, file.path(summary_dir, "t05_overall_summary.tsv"))

plot_data <- all_status[, .N, by = .(panel, final_status)]
plot_data[, final_status := factor(final_status, levels = c(
  "local_rg_evidence_fdr",
  "tested_no_local_rg_evidence",
  "hla_adjacent_tested_separately",
  "insufficient_joint_local_heritability",
  "insufficient_local_h2_estimation",
  "insufficient_data_build_mismatch_no_rsid",
  "insufficient_data_processing_failure",
  "insufficient_bivariate_result"
))]
p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = panel, y = N, fill = final_status)) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = N),
    position = ggplot2::position_stack(vjust = 0.5),
    colour = "white",
    size = 3
  ) +
  ggplot2::scale_fill_manual(values = c(
    local_rg_evidence_fdr = "#1B7837",
    tested_no_local_rg_evidence = "#4575B4",
    hla_adjacent_tested_separately = "#D73027",
    insufficient_joint_local_heritability = "#BDBDBD",
    insufficient_local_h2_estimation = "#737373",
    insufficient_data_build_mismatch_no_rsid = "#252525",
    insufficient_data_processing_failure = "#8C510A",
    insufficient_bivariate_result = "#DFC27D"
  ), labels = c(
    local_rg_evidence_fdr = "Local rg evidence after FDR",
    tested_no_local_rg_evidence = "Tested: no local rg evidence",
    hla_adjacent_tested_separately = "HLA-adjacent: separate interpretation",
    insufficient_joint_local_heritability = "Insufficient joint local h2",
    insufficient_local_h2_estimation = "Local h2 not estimable",
    insufficient_data_build_mismatch_no_rsid = "Build mismatch without rsID list",
    insufficient_data_processing_failure = "Processing failure",
    insufficient_bivariate_result = "Bivariate result unavailable"
  ), drop = TRUE) +
  ggplot2::labs(
    x = NULL,
    y = "Number of loci",
    fill = "Final testability status",
    title = "T05 local genetic-correlation testability"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "right")
ggplot2::ggsave(
  file.path(figure_dir, "t05_lava_testability.png"),
  p, width = 9, height = 5.5, dpi = 300, bg = "white"
)

tested_targeted <- targeted$status[is.finite(local_rg_p)]
if (nrow(tested_targeted)) {
  tested_targeted[, label := paste0(ccgwas_locus_id, " (", lead_nearest_gene, ")")]
  forest <- ggplot2::ggplot(
    tested_targeted,
    ggplot2::aes(x = local_rg, y = reorder(label, local_rg), colour = hla_adjacent_within_2mb)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = local_rg_ci_lower, xmax = local_rg_ci_upper),
      width = 0.12, linewidth = 0.6, orientation = "y"
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#1F5A7A", `TRUE` = "#B2182B"),
      labels = c(`FALSE` = "Ordinary non-HLA", `TRUE` = "Within 2 Mb of HLA")
    ) +
    ggplot2::coord_cartesian(xlim = c(-1, 1)) +
    ggplot2::labs(x = "Local genetic correlation (rho), 95% CI", y = NULL, colour = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
  ggplot2::ggsave(
    file.path(figure_dir, "t05_targeted_local_rg.png"),
    forest, width = 8, height = 3.8, dpi = 300, bg = "white"
  )
}

log_message("T05 summary complete: discovery bivariate tests = ", sum(is.finite(discovery$status$local_rg_p)),
            "; targeted bivariate tests = ", sum(is.finite(targeted$status$local_rg_p)))
