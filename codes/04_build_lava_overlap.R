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
result_dir <- file.path(project_root, cfg$project$outcome_dir, "ldsc", "results")

extract_rg_row <- function(path) {
  if (!file.exists(path)) stop("Missing LDSC log: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  trimmed <- trimws(lines)
  header_index <- grep("^p1[[:space:]]+p2[[:space:]]+rg[[:space:]]+se", trimmed)
  if (!length(header_index)) stop("Could not find LDSC rg table in ", path, call. = FALSE)
  header_index <- tail(header_index, 1L)
  header <- strsplit(trimmed[header_index], "[[:space:]]+")[[1]]
  candidates <- trimmed[seq.int(header_index + 1L, length(trimmed))]
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[grepl("sumstats\\.gz", candidates)]
  token_counts <- lengths(strsplit(candidates, "[[:space:]]+"))
  candidates <- candidates[token_counts == length(header)]
  if (!length(candidates)) stop("Could not parse LDSC rg result row in ", path, call. = FALSE)
  result <- data.table::fread(
    text = paste(c(paste(header, collapse = "\t"), gsub("[[:space:]]+", "\t", candidates[1])), collapse = "\n"),
    sep = "\t"
  )
  # The fixed-width summary table rounds small P values to 0.0000. Preserve
  # the exact genetic-correlation P value printed immediately above it.
  detailed_p <- suppressWarnings(as.numeric(sub(
    "^P:[[:space:]]*", "", grep("^P:[[:space:]]*", trimmed, value = TRUE)
  )))
  if (length(detailed_p) && is.finite(tail(detailed_p, 1L))) {
    result[, p := tail(detailed_p, 1L)]
  }
  result
}

paths <- c(
  graves_self = file.path(result_dir, paste0(graves, "__", graves, ".log")),
  aiht_self = file.path(result_dir, paste0(aiht, "__", aiht, ".log")),
  cross = file.path(result_dir, paste0(graves, "__", aiht, ".log"))
)
rows <- lapply(paths, extract_rg_row)
required <- c("h2_int", "gcov_int")
if (any(!vapply(rows, function(x) all(required %in% names(x)), logical(1)))) {
  stop("LDSC result tables must contain h2_int and gcov_int.", call. = FALSE)
}
write_tsv(
  data.table::rbindlist(rows, idcol = "analysis"),
  file.path(result_dir, "ldsc_rg_summary.tsv")
)

diagonal_value <- function(x) {
  value <- as.numeric(x$gcov_int[1])
  if (!is.finite(value) || value <= 0) value <- as.numeric(x$h2_int[1])
  value
}
raw_overlap <- matrix(
  c(
    diagonal_value(rows$graves_self), as.numeric(rows$cross$gcov_int[1]),
    as.numeric(rows$cross$gcov_int[1]), diagonal_value(rows$aiht_self)
  ),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(c(graves, aiht), c(graves, aiht))
)
if (any(!is.finite(raw_overlap)) || any(diag(raw_overlap) <= 0)) {
  stop("Non-finite or non-positive LDSC intercept values; inspect the logs before LAVA.", call. = FALSE)
}
overlap <- round(stats::cov2cor(raw_overlap), 6L)
if (min(eigen(overlap, symmetric = TRUE, only.values = TRUE)$values) < -1e-6) {
  stop("The LDSC-derived overlap matrix is not positive semidefinite.", call. = FALSE)
}

out_path <- project_path(project_root, cfg$lava$sample_overlap_file)
ensure_dir(dirname(out_path))
write.table(overlap, out_path, quote = FALSE, sep = "\t", col.names = NA)
provenance <- data.table::data.table(
  phenotype_1 = c(graves, graves, aiht),
  phenotype_2 = c(graves, aiht, aiht),
  raw_intercept = c(raw_overlap[1, 1], raw_overlap[1, 2], raw_overlap[2, 2]),
  source_log = unname(paths[c("graves_self", "cross", "aiht_self")])
)
write_tsv(provenance, file.path(dirname(out_path), "sample_overlap_provenance.tsv"))
log_message("Wrote LAVA sample-overlap matrix to ", out_path, ". Review it before running LAVA.")
