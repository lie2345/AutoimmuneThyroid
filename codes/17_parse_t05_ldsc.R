#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
graves <- cfg$primary_comparison$graves
aiht <- cfg$primary_comparison$autoimmune_hypothyroidism
result_dir <- file.path(project_root, "outcome_data", "t05", "ldsc", "results")
lava_input_dir <- file.path(project_root, "outcome_data", "t05", "lava", "input")
ensure_dir(lava_input_dir)

extract_rg_row <- function(path) {
  if (!file.exists(path)) stop("Missing LDSC log: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  trimmed <- trimws(lines)
  header_index <- grep("^p1[[:space:]]+p2[[:space:]]+rg[[:space:]]+se", trimmed)
  if (!length(header_index)) stop("Could not find LDSC rg table in ", path, call. = FALSE)
  header_index <- tail(header_index, 1L)
  header <- strsplit(trimmed[header_index], "[[:space:]]+")[[1L]]
  candidates <- trimmed[seq.int(header_index + 1L, length(trimmed))]
  candidates <- candidates[nzchar(candidates) & grepl("sumstats\\.gz", candidates)]
  candidates <- candidates[lengths(strsplit(candidates, "[[:space:]]+")) == length(header)]
  if (!length(candidates)) stop("Could not parse LDSC result row in ", path, call. = FALSE)
  result <- data.table::fread(
    text = paste(
      c(paste(header, collapse = "\t"), gsub("[[:space:]]+", "\t", candidates[[1L]])),
      collapse = "\n"
    ),
    sep = "\t"
  )
  detailed_p <- suppressWarnings(as.numeric(sub(
    "^P:[[:space:]]*", "", grep("^P:[[:space:]]*", trimmed, value = TRUE)
  )))
  if (length(detailed_p) && is.finite(tail(detailed_p, 1L))) result[, p := tail(detailed_p, 1L)]
  result
}

paths <- c(
  graves_self = file.path(result_dir, paste0(graves, "__", graves, ".log")),
  aiht_self = file.path(result_dir, paste0(aiht, "__", aiht, ".log")),
  cross = file.path(result_dir, paste0(graves, "__", aiht, ".log"))
)
rows <- lapply(paths, extract_rg_row)
summary <- data.table::rbindlist(rows, idcol = "analysis")
summary[, `:=`(
  rg_ci_lower = rg - 1.96 * se,
  rg_ci_upper = rg + 1.96 * se,
  h2_z = h2_obs / h2_obs_se,
  h2_p = 2 * stats::pnorm(-abs(h2_obs / h2_obs_se))
)]
write_tsv(summary, file.path(result_dir, "t05_ldsc_summary.tsv"))

diagonal_value <- function(x) {
  value <- as.numeric(x$gcov_int[[1L]])
  if (!is.finite(value) || value <= 0) value <- as.numeric(x$h2_int[[1L]])
  value
}
raw_intercepts <- matrix(
  c(
    diagonal_value(rows$graves_self), as.numeric(rows$cross$gcov_int[[1L]]),
    as.numeric(rows$cross$gcov_int[[1L]]), diagonal_value(rows$aiht_self)
  ),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(c(graves, aiht), c(graves, aiht))
)
overlap <- round(stats::cov2cor(raw_intercepts), 6L)
if (any(!is.finite(overlap)) || min(eigen(overlap, symmetric = TRUE, only.values = TRUE)$values) < -1e-6) {
  stop("Invalid LDSC-derived overlap matrix", call. = FALSE)
}
overlap_path <- file.path(lava_input_dir, "sample.overlap.txt")
write.table(overlap, overlap_path, quote = FALSE, sep = "\t", col.names = NA)

provenance <- data.table::data.table(
  analysis_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  genome_build = cfg$project$genome_build,
  input_scope = "autosomal SNPs MAF>=0.01, palindromic excluded, HLA excluded before LDSC munging",
  hla_interval_excluded = paste0("chr", cfg$hla$chromosome, ":", cfg$hla$start, "-", cfg$hla$end),
  ld_reference = project_path(project_root, cfg$ldsc$reference_ld_prefix),
  regression_weights = project_path(project_root, cfg$ldsc$regression_weight_prefix),
  regression_weights_scope = "HapMap3 no MHC",
  cross_trait_intercept = as.numeric(rows$cross$gcov_int[[1L]]),
  standardized_overlap_proxy = overlap[graves, aiht],
  overlap_proxy_note = "LDSC-derived error-correlation proxy for LAVA; not an overlapping-person count"
)
write_tsv(provenance, file.path(result_dir, "t05_ldsc_provenance.tsv"))

old_path <- file.path(project_root, "outcome_data", "ldsc", "results", "ldsc_rg_summary.tsv")
if (file.exists(old_path)) {
  old <- data.table::fread(old_path)
  comparison <- merge(
    old[, .(analysis, old_rg = rg, old_rg_se = se, old_rg_p = p, old_cross_intercept = gcov_int)],
    summary[, .(analysis, t05_rg = rg, t05_rg_se = se, t05_rg_p = p, t05_cross_intercept = gcov_int)],
    by = "analysis", all = TRUE
  )
  comparison[, `:=`(
    delta_rg = t05_rg - old_rg,
    delta_rg_se = t05_rg_se - old_rg_se,
    delta_cross_intercept = t05_cross_intercept - old_cross_intercept
  )]
  write_tsv(comparison, file.path(result_dir, "comparison_with_previous_ldsc.tsv"))
}

log_message("T05 LDSC parsed; cross-trait rg = ", summary[analysis == "cross", rg],
            ", intercept = ", summary[analysis == "cross", gcov_int])
