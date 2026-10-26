#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "ggplot2", "scales"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)

power_cfg <- cfg$path_b_power
if (is.null(power_cfg)) stop("Missing path_b_power in config/analysis.yml", call. = FALSE)

graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
primary <- traits[trait_id %in% c(graves_id, aiht_id)]
if (nrow(primary) != 2L) stop("Primary trait rows are missing from config/traits.tsv", call. = FALSE)

trait_value <- function(trait_id, column) {
  requested_trait <- trait_id
  value <- primary[trait_id == requested_trait, get(column)]
  if (length(value) != 1L || !is.finite(as.numeric(value))) {
    stop("Invalid ", column, " for ", trait_id, call. = FALSE)
  }
  as.numeric(value)
}

n_graves_case <- trait_value(graves_id, "n_cases")
n_graves_control <- trait_value(graves_id, "n_controls")
n_aiht_case <- trait_value(aiht_id, "n_cases")
n_aiht_control <- trait_value(aiht_id, "n_controls")

alpha_gws <- as.numeric(power_cfg$genome_wide_alpha)
power_targets <- as.numeric(unlist(power_cfg$power_targets))
maf_values <- as.numeric(unlist(power_cfg$maf_values))
primary_rho <- as.numeric(power_cfg$primary_error_correlation)
effect_retention <- as.numeric(unlist(power_cfg$effect_retention_sensitivity))
equivalence_alpha <- as.numeric(power_cfg$equivalence_alpha)
equivalence_or_primary <- as.numeric(power_cfg$equivalence_or_upper_primary)
equivalence_or_sensitivity <- as.numeric(unlist(power_cfg$equivalence_or_upper_sensitivity))

if (!is.finite(alpha_gws) || alpha_gws <= 0 || alpha_gws >= 1) {
  stop("genome_wide_alpha must be in (0, 1)", call. = FALSE)
}
if (any(!is.finite(power_targets)) || any(power_targets <= 0 | power_targets >= 1)) {
  stop("power_targets must all be in (0, 1)", call. = FALSE)
}
if (any(!is.finite(maf_values)) || any(maf_values <= 0 | maf_values > 0.5)) {
  stop("maf_values must all be in (0, 0.5]", call. = FALSE)
}
if (!is.finite(primary_rho) || abs(primary_rho) >= 1) {
  stop("primary_error_correlation must be in (-1, 1)", call. = FALSE)
}
if (any(!is.finite(effect_retention)) || any(effect_retention <= 0 | effect_retention > 1)) {
  stop("effect_retention_sensitivity must be in (0, 1]", call. = FALSE)
}
if (!is.finite(equivalence_alpha) || equivalence_alpha <= 0 || equivalence_alpha >= 0.5) {
  stop("equivalence_alpha must be in (0, 0.5)", call. = FALSE)
}
equivalence_or_values <- sort(unique(c(equivalence_or_primary, equivalence_or_sensitivity)))
if (any(!is.finite(equivalence_or_values)) || any(equivalence_or_values <= 1)) {
  stop("Equivalence OR upper bounds must exceed 1", call. = FALSE)
}

case_control_information <- function(n_case, n_control) {
  n_case * n_control / (n_case + n_control)
}

genotype_variance <- function(maf) {
  2 * maf * (1 - maf)
}

se_case_control <- function(maf, n_case, n_control) {
  1 / sqrt(genotype_variance(maf) * case_control_information(n_case, n_control))
}

se_direct_case_case <- function(maf, n_case_a, n_case_b) {
  1 / sqrt(genotype_variance(maf) * case_control_information(n_case_a, n_case_b))
}

se_summary_contrast <- function(maf, rho) {
  se_a <- se_case_control(maf, n_graves_case, n_graves_control)
  se_b <- se_case_control(maf, n_aiht_case, n_aiht_control)
  variance <- se_a^2 + se_b^2 - 2 * rho * se_a * se_b
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Non-positive summary-contrast variance", call. = FALSE)
  }
  sqrt(variance)
}

two_sided_power <- function(beta, se, alpha) {
  critical <- stats::qnorm(1 - alpha / 2)
  ncp <- abs(beta) / se
  stats::pnorm(-critical, mean = ncp) + stats::pnorm(critical, mean = ncp, lower.tail = FALSE)
}

minimum_detectable_beta <- function(se, target_power, alpha) {
  vapply(seq_along(se), function(i) {
    upper <- max(1, 20 * se[[i]])
    while (two_sided_power(upper, se[[i]], alpha) < target_power) upper <- upper * 2
    stats::uniroot(
      function(beta) two_sided_power(beta, se[[i]], alpha) - target_power,
      interval = c(0, upper),
      tol = 1e-12
    )$root
  }, numeric(1))
}

overlap_path <- project_path(project_root, cfg$lava$sample_overlap_file)
if (!file.exists(overlap_path)) stop("Missing LAVA overlap matrix: ", overlap_path, call. = FALSE)
overlap_dt <- data.table::fread(overlap_path)
if (ncol(overlap_dt) < 3L) stop("Malformed LAVA overlap matrix", call. = FALSE)
overlap_ids <- as.character(overlap_dt[[1L]])
overlap_matrix <- as.matrix(overlap_dt[, -1L, with = FALSE])
rownames(overlap_matrix) <- overlap_ids
if (!all(c(graves_id, aiht_id) %in% rownames(overlap_matrix)) ||
    !all(c(graves_id, aiht_id) %in% colnames(overlap_matrix))) {
  stop("Primary traits missing from LAVA overlap matrix", call. = FALSE)
}
ldsc_rho_proxy <- as.numeric(overlap_matrix[graves_id, aiht_id])
if (!is.finite(ldsc_rho_proxy) || abs(ldsc_rho_proxy) >= 1) {
  stop("Invalid LDSC overlap proxy", call. = FALSE)
}

scenario_table <- data.table::data.table(
  scenario = c(
    "ideal_direct_case_case",
    "summary_contrast_primary_unknown_overlap",
    "summary_contrast_ldsc_proxy_sensitivity"
  ),
  scenario_label = c(
    "Ideal direct case-case",
    "Summary contrast, primary (rho=0)",
    sprintf("Summary contrast, LDSC proxy (rho=%.6f)", ldsc_rho_proxy)
  ),
  error_correlation = c(NA_real_, primary_rho, ldsc_rho_proxy),
  primary_analysis = c(FALSE, TRUE, FALSE),
  interpretation = c(
    "Individual-level benchmark using Graves cases versus AIHT cases",
    "Conservative summary-statistic benchmark; exact overlapping controls are unknown",
    "Sensitivity benchmark only; LDSC-normalized proxy is not an overlapping-control count"
  )
)

power_rows <- vector("list", nrow(scenario_table) * length(power_targets))
index <- 0L
for (scenario_index in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table[scenario_index]
  if (scenario$scenario == "ideal_direct_case_case") {
    se_values <- se_direct_case_case(maf_values, n_graves_case, n_aiht_case)
  } else {
    se_values <- se_summary_contrast(maf_values, scenario$error_correlation)
  }
  for (target_power in power_targets) {
    index <- index + 1L
    mde <- minimum_detectable_beta(se_values, target_power, alpha_gws)
    power_rows[[index]] <- data.table::data.table(
      scenario = scenario$scenario,
      scenario_label = scenario$scenario_label,
      primary_analysis = scenario$primary_analysis,
      error_correlation = scenario$error_correlation,
      maf = maf_values,
      maf_class = ifelse(maf_values < 0.05, "low_frequency", "common"),
      target_power = target_power,
      alpha = alpha_gws,
      se_log_or_difference = se_values,
      mde_log_or = mde,
      mde_or_above_1 = exp(mde),
      mde_or_below_1 = exp(-mde),
      power_scope = "Necessary OLS-like two-sided threshold; not full CC-GWAS significance power"
    )
  }
}
power_table <- data.table::rbindlist(power_rows)

primary_power <- power_table[scenario == "summary_contrast_primary_unknown_overlap"]
distortion_table <- primary_power[rep(seq_len(.N), each = length(effect_retention))]
distortion_table[, effect_retention_fraction := rep(effect_retention, times = nrow(primary_power))]
distortion_table[, `:=`(
  true_mde_log_or_if_attenuated = mde_log_or / effect_retention_fraction,
  true_mde_or_above_1_if_attenuated = exp(mde_log_or / effect_retention_fraction),
  note = "Illustrative multiplicative attenuation only; control-definition bias can also be nonmultiplicative or directional"
)]

equivalence_power_null <- function(se, delta, alpha) {
  z <- stats::qnorm(1 - alpha)
  margin <- delta / se - z
  pmax(0, 2 * stats::pnorm(margin) - 1)
}

trait_table <- data.table::data.table(
  trait_id = c(graves_id, aiht_id),
  trait_label = c("Graves", "AIHT"),
  n_cases = c(n_graves_case, n_aiht_case),
  n_controls = c(n_graves_control, n_aiht_control)
)

equivalence_rows <- list()
index <- 0L
for (trait_index in seq_len(nrow(trait_table))) {
  trait <- trait_table[trait_index]
  se_values <- se_case_control(maf_values, trait$n_cases, trait$n_controls)
  for (or_upper in equivalence_or_values) {
    index <- index + 1L
    delta <- log(or_upper)
    equivalence_rows[[index]] <- data.table::data.table(
      trait_id = trait$trait_id,
      trait_label = trait$trait_label,
      n_cases = trait$n_cases,
      n_controls = trait$n_controls,
      maf = maf_values,
      maf_class = ifelse(maf_values < 0.05, "low_frequency", "common"),
      equivalence_alpha = equivalence_alpha,
      equivalence_or_lower = exp(-delta),
      equivalence_or_upper = exp(delta),
      primary_equivalence_bound = abs(or_upper - equivalence_or_primary) < 1e-12,
      se_log_or = se_values,
      equivalence_power_if_true_beta_zero = equivalence_power_null(se_values, delta, equivalence_alpha),
      equivalence_rule = "90% CI must lie wholly inside the symmetric log-OR bounds"
    )
  }
}
equivalence_table <- data.table::rbindlist(equivalence_rows)

minimum_maf_for_equivalence <- function(n_case, n_control, delta, target_power, alpha) {
  objective <- function(maf) {
    se <- se_case_control(maf, n_case, n_control)
    equivalence_power_null(se, delta, alpha) - target_power
  }
  if (objective(0.5) < 0) return(NA_real_)
  stats::uniroot(objective, interval = c(1e-8, 0.5), tol = 1e-12)$root
}

threshold_rows <- list()
index <- 0L
for (trait_index in seq_len(nrow(trait_table))) {
  trait <- trait_table[trait_index]
  for (or_upper in equivalence_or_values) {
    for (target_power in power_targets) {
      index <- index + 1L
      delta <- log(or_upper)
      max_se <- delta / (
        stats::qnorm(1 - equivalence_alpha) +
          stats::qnorm((target_power + 1) / 2)
      )
      threshold_rows[[index]] <- data.table::data.table(
        trait_id = trait$trait_id,
        trait_label = trait$trait_label,
        equivalence_or_lower = exp(-delta),
        equivalence_or_upper = exp(delta),
        primary_equivalence_bound = abs(or_upper - equivalence_or_primary) < 1e-12,
        target_power = target_power,
        equivalence_alpha = equivalence_alpha,
        maximum_se_for_target_power = max_se,
        minimum_maf_if_true_beta_zero = minimum_maf_for_equivalence(
          trait$n_cases, trait$n_controls, delta, target_power, equivalence_alpha
        )
      )
    }
  }
}
equivalence_thresholds <- data.table::rbindlist(threshold_rows)

neff <- function(n_case, n_control) 4 * case_control_information(n_case, n_control)
neff_graves <- neff(n_graves_case, n_graves_control)
neff_aiht <- neff(n_aiht_case, n_aiht_control)
neff_direct <- neff(n_graves_case, n_aiht_case)
neff_summary <- function(rho) {
  1 / (1 / neff_graves + 1 / neff_aiht - 2 * rho / sqrt(neff_graves * neff_aiht))
}

parameter_table <- data.table::data.table(
  parameter = c(
    "n_graves_cases", "n_graves_controls", "n_aiht_cases", "n_aiht_controls",
    "neff_graves_case_control", "neff_aiht_case_control",
    "neff_ideal_direct_case_case", "neff_summary_primary_rho0",
    "neff_summary_ldsc_proxy", "primary_error_correlation",
    "ldsc_normalized_error_correlation_proxy", "genome_wide_alpha",
    "equivalence_alpha", "primary_equivalence_or_upper",
    "primary_equivalence_or_lower"
  ),
  value = c(
    n_graves_case, n_graves_control, n_aiht_case, n_aiht_control,
    neff_graves, neff_aiht, neff_direct, neff_summary(primary_rho),
    neff_summary(ldsc_rho_proxy), primary_rho, ldsc_rho_proxy, alpha_gws,
    equivalence_alpha, equivalence_or_primary, 1 / equivalence_or_primary
  ),
  role = c(
    rep("input", 4), rep("derived_effective_sample_size", 5),
    "primary_conservative_assumption", "sensitivity_only", "testing_threshold",
    "testing_threshold", "equivalence_threshold", "equivalence_threshold"
  )
)

output_dir <- file.path(project_root, cfg$project$outcome_dir, "power")
figure_dir <- file.path(project_root, cfg$project$figure_dir, "power")
ensure_dir(output_dir)
ensure_dir(figure_dir)

write_tsv(power_table, file.path(output_dir, "case_case_power.tsv"))
write_tsv(distortion_table, file.path(output_dir, "control_definition_sensitivity.tsv"))
write_tsv(equivalence_table, file.path(output_dir, "equivalence_power.tsv"))
write_tsv(equivalence_thresholds, file.path(output_dir, "equivalence_maf_thresholds.tsv"))
write_tsv(parameter_table, file.path(output_dir, "analysis_parameters.tsv"))
write_tsv(scenario_table, file.path(output_dir, "scenario_definitions.tsv"))

curve_maf <- exp(seq(log(0.005), log(0.5), length.out = 300L))
curve_rows <- list()
index <- 0L
for (scenario_index in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table[scenario_index]
  if (scenario$scenario == "ideal_direct_case_case") {
    se_values <- se_direct_case_case(curve_maf, n_graves_case, n_aiht_case)
  } else {
    se_values <- se_summary_contrast(curve_maf, scenario$error_correlation)
  }
  for (target_power in power_targets) {
    index <- index + 1L
    mde <- minimum_detectable_beta(se_values, target_power, alpha_gws)
    curve_rows[[index]] <- data.table::data.table(
      scenario_label = scenario$scenario_label,
      maf = curve_maf,
      target_power = paste0(round(100 * target_power), "% power"),
      mde_or_above_1 = exp(mde)
    )
  }
}
power_curve <- data.table::rbindlist(curve_rows)
power_curve[, scenario_label := factor(scenario_label, levels = scenario_table$scenario_label)]

power_plot <- ggplot2::ggplot(
  power_curve,
  ggplot2::aes(x = maf, y = mde_or_above_1, color = scenario_label, linetype = target_power)
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::scale_x_log10(
    breaks = c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
    labels = c("0.5%", "1%", "2%", "5%", "10%", "20%", "50%")
  ) +
  ggplot2::scale_color_manual(values = c("#4D4D4D", "#0072B2", "#D55E00")) +
  ggplot2::labs(
    x = "Minor allele frequency",
    y = "Minimum detectable case-case OR (>1)",
    color = NULL,
    linetype = NULL,
    title = "Path B analytic detection benchmark",
    subtitle = "Two-sided alpha = 5e-8; CC-GWAS stress-test filters make this a necessary, not sufficient, benchmark"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title.position = "plot"
  )

ggplot2::ggsave(
  file.path(figure_dir, "case_case_power.png"),
  power_plot, width = 8.2, height = 5.6, units = "in", dpi = 300, bg = "white"
)
ggplot2::ggsave(
  file.path(figure_dir, "case_case_power.pdf"),
  power_plot, width = 8.2, height = 5.6, units = "in", device = grDevices::cairo_pdf
)

equivalence_curve_rows <- list()
index <- 0L
for (trait_index in seq_len(nrow(trait_table))) {
  trait <- trait_table[trait_index]
  se_values <- se_case_control(curve_maf, trait$n_cases, trait$n_controls)
  index <- index + 1L
  equivalence_curve_rows[[index]] <- data.table::data.table(
    trait_label = trait$trait_label,
    maf = curve_maf,
    equivalence_power = equivalence_power_null(
      se_values, log(equivalence_or_primary), equivalence_alpha
    )
  )
}
equivalence_curve <- data.table::rbindlist(equivalence_curve_rows)

equivalence_plot <- ggplot2::ggplot(
  equivalence_curve,
  ggplot2::aes(x = maf, y = equivalence_power, color = trait_label)
) +
  ggplot2::geom_hline(yintercept = c(0.8, 0.9), color = "#9E9E9E", linetype = c("dashed", "dotted")) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_x_log10(
    breaks = c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
    labels = c("0.5%", "1%", "2%", "5%", "10%", "20%", "50%")
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_color_manual(values = c("Graves" = "#0072B2", "AIHT" = "#D55E00")) +
  ggplot2::labs(
    x = "Minor allele frequency",
    y = "Power to conclude equivalence when true beta = 0",
    color = NULL,
    title = "Power for the primary near-null equivalence rule",
    subtitle = sprintf("TOST alpha = %.2f; symmetric log-OR interval %.3f to %.2f", equivalence_alpha, 1 / equivalence_or_primary, equivalence_or_primary)
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", plot.title.position = "plot")

ggplot2::ggsave(
  file.path(figure_dir, "equivalence_power.png"),
  equivalence_plot, width = 7.5, height = 5.2, units = "in", dpi = 300, bg = "white"
)
ggplot2::ggsave(
  file.path(figure_dir, "equivalence_power.pdf"),
  equivalence_plot, width = 7.5, height = 5.2, units = "in", device = grDevices::cairo_pdf
)

expected_power_rows <- length(maf_values) * length(power_targets) * nrow(scenario_table)
if (nrow(power_table) != expected_power_rows) stop("Unexpected power-table row count", call. = FALSE)
if (any(power_table$mde_or_above_1 <= 1) || any(power_table$mde_or_below_1 >= 1)) {
  stop("Invalid MDE OR values", call. = FALSE)
}
if (any(equivalence_table$equivalence_power_if_true_beta_zero < 0 |
        equivalence_table$equivalence_power_if_true_beta_zero > 1)) {
  stop("Invalid equivalence power", call. = FALSE)
}
if (abs(neff_direct - 4 * n_graves_case * n_aiht_case / (n_graves_case + n_aiht_case)) > 1e-8) {
  stop("Direct case-case effective sample size validation failed", call. = FALSE)
}

log_message("T03 outputs written to ", output_dir)
log_message("Ideal direct case-case effective N = ", round(neff_direct, 1))
log_message("Primary summary-contrast effective N = ", round(neff_summary(primary_rho), 1))
log_message("LDSC-proxy summary-contrast effective N = ", round(neff_summary(ldsc_rho_proxy), 1))
