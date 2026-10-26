#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "R.utils"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)
cc_cfg <- cfg$ccgwas
if (is.null(cc_cfg)) stop("Missing ccgwas configuration", call. = FALSE)

local_library <- project_path(project_root, cc_cfg$r_library)
.libPaths(c(local_library, .libPaths()))
if (!requireNamespace("CCGWAS", quietly = TRUE)) {
  stop("CCGWAS is not installed in ", local_library, call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(CCGWAS))

graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
meta <- traits[trait_id %in% c(graves_id, aiht_id)]
if (nrow(meta) != 2L) stop("Primary trait metadata are incomplete", call. = FALSE)
trait_value <- function(id, column) as.numeric(meta[trait_id == id, get(column)])

n_g_case <- trait_value(graves_id, "n_cases")
n_g_control <- trait_value(graves_id, "n_controls")
n_a_case <- trait_value(aiht_id, "n_cases")
n_a_control <- trait_value(aiht_id, "n_controls")

ldsc_path <- file.path(project_root, "outcome_data", "ldsc", "results", "ldsc_rg_summary.tsv")
ldsc <- data.table::fread(ldsc_path)
row_one <- function(label) {
  value <- ldsc[analysis == label]
  if (nrow(value) != 1L) stop("Missing LDSC row: ", label, call. = FALSE)
  value
}
g_ldsc <- row_one("graves_self")
a_ldsc <- row_one("aiht_self")
cross_ldsc <- row_one("cross")

prev_g <- cc_cfg$prevalence$graves
prev_a <- cc_cfg$prevalence$aiht
liability_h2 <- function(h2_observed, n_case, n_control, prevalence) {
  sample_prevalence <- n_case / (n_case + n_control)
  density <- stats::dnorm(stats::qnorm(1 - prevalence))
  h2_observed * prevalence^2 * (1 - prevalence)^2 /
    (sample_prevalence * (1 - sample_prevalence) * density^2)
}
h2l_g <- liability_h2(g_ldsc$h2_obs, n_g_case, n_g_control, as.numeric(prev_g$main))
h2l_a <- liability_h2(a_ldsc$h2_obs, n_a_case, n_a_control, as.numeric(prev_a$main))
if (any(!is.finite(c(h2l_g, h2l_a))) || any(c(h2l_g, h2l_a) <= 0 | c(h2l_g, h2l_a) >= 1)) {
  stop("Invalid liability-scale heritability", call. = FALSE)
}

out_root <- project_path(project_root, cc_cfg$output_dir)
input_dir <- file.path(out_root, "input")
results_dir <- file.path(out_root, "results")
ensure_dir(results_dir)
g_file <- file.path(input_dir, paste0(graves_id, ".ccgwas.tsv.gz"))
a_file <- file.path(input_dir, paste0(aiht_id, ".ccgwas.tsv.gz"))
if (!file.exists(g_file) || !file.exists(a_file)) stop("Run codes/13_prepare_ccgwas_inputs.R first", call. = FALSE)

source_dir <- project_path(project_root, cc_cfg$source_dir)
git_commit <- tryCatch(
  system2("git", c("-C", source_dir, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1L],
  error = function(e) NA_character_
)

overlap_fractions <- as.numeric(unlist(cc_cfg$overlap_fraction_sensitivity))
if (!0 %in% overlap_fractions) stop("Primary overlap fraction 0 is required", call. = FALSE)
scenario_table <- data.table::data.table(
  scenario = paste0("overlap_", sprintf("%03d", round(100 * overlap_fractions)), "pct"),
  overlap_fraction_of_smaller_control_set = overlap_fractions,
  N_overlap_A0B0 = as.integer(round(overlap_fractions * min(n_g_control, n_a_control))),
  primary_analysis = overlap_fractions == 0
)
scenario_table[, `:=`(
  K_graves_main = as.numeric(prev_g$main),
  K_graves_low = as.numeric(prev_g$low),
  K_graves_high = as.numeric(prev_g$high),
  K_aiht_main = as.numeric(prev_a$main),
  K_aiht_low = as.numeric(prev_a$low),
  K_aiht_high = as.numeric(prev_a$high),
  h2l_graves = h2l_g,
  h2l_aiht = h2l_a,
  rg = as.numeric(cross_ldsc$rg),
  empirical_cross_trait_intercept = as.numeric(cross_ldsc$gcov_int),
  m_effective_loci = as.numeric(cc_cfg$m_effective_loci),
  subtype_data = isTRUE(cc_cfg$subtype_data),
  ccgwas_git_commit = git_commit,
  input_scope = "HapMap3 MAF>0.01 autosomal non-HLA"
)]

analytical_intercept <- function(n_overlap) {
  (n_overlap / (4 * n_g_control * n_a_control)) /
    sqrt(
      (1 / (4 * n_g_case) + 1 / (4 * n_g_control)) *
        (1 / (4 * n_a_case) + 1 / (4 * n_a_control))
    )
}
scenario_table[, analytical_overlap_intercept := analytical_intercept(N_overlap_A0B0)]
scenario_table[, modelled_intercept := pmin(empirical_cross_trait_intercept, analytical_overlap_intercept)]
data.table::fwrite(
  scenario_table,
  file.path(results_dir, "scenario_parameters.tsv"),
  sep = "\t", quote = FALSE
)

status <- vector("list", nrow(scenario_table))
for (i in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table[i]
  prefix <- file.path(results_dir, scenario$scenario)
  expected <- paste0(prefix, ".results.gz")
  if (file.exists(expected)) stop("Refusing to overwrite existing result: ", expected, call. = FALSE)
  log_message("Running CC-GWAS scenario ", scenario$scenario, " (N overlap = ", scenario$N_overlap_A0B0, ")")
  started <- Sys.time()
  error_message <- NA_character_
  tryCatch(
    CCGWAS(
      outcome_file = prefix,
      A_name = "GD",
      B_name = "AIHT",
      sumstats_fileA1A0 = g_file,
      sumstats_fileB1B0 = a_file,
      K_A1A0 = as.numeric(prev_g$main),
      K_A1A0_high = as.numeric(prev_g$high),
      K_A1A0_low = as.numeric(prev_g$low),
      K_B1B0 = as.numeric(prev_a$main),
      K_B1B0_high = as.numeric(prev_a$high),
      K_B1B0_low = as.numeric(prev_a$low),
      h2l_A1A0 = h2l_g,
      h2l_B1B0 = h2l_a,
      rg_A1A0_B1B0 = as.numeric(cross_ldsc$rg),
      intercept_A1A0_B1B0 = as.numeric(cross_ldsc$gcov_int),
      m = as.numeric(cc_cfg$m_effective_loci),
      N_A1 = n_g_case,
      N_B1 = n_a_case,
      N_A0 = n_g_control,
      N_B0 = n_a_control,
      N_overlap_A0B0 = scenario$N_overlap_A0B0,
      subtype_data = isTRUE(cc_cfg$subtype_data),
      sumstats_fileA1B1 = NA,
      save.all = isTRUE(scenario$primary_analysis)
    ),
    error = function(e) error_message <<- conditionMessage(e)
  )
  status[[i]] <- data.table::data.table(
    scenario = scenario$scenario,
    success = is.na(error_message) && file.exists(expected),
    error = error_message,
    started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    ended = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    result_file = expected
  )
  data.table::fwrite(
    data.table::rbindlist(status[seq_len(i)], fill = TRUE),
    file.path(results_dir, "run_status.tsv"),
    sep = "\t", quote = FALSE, na = "NA"
  )
  if (!is.na(error_message)) stop("CC-GWAS failed for ", scenario$scenario, ": ", error_message, call. = FALSE)
}

log_message("All ", nrow(scenario_table), " CC-GWAS scenarios completed")
