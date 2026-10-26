#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("data.table", "Rsamtools", "GenomicRanges", "IRanges"))
args <- parse_cli()

option_value <- function(name, default) {
  value <- args[[name]]
  if (is.null(value)) default else value
}

mode <- option_value("mode", "extract")
max_priority <- option_value("max-priority", "P1_strong")
plan_path <- project_path(
  project_root,
  option_value("plan", "outcome_data/targeted_susie/acquisition_plan.tsv")
)
slice_root <- project_path(
  project_root,
  option_value("slice-root", "outcome_data/targeted_susie/slices")
)
download_root <- project_path(
  project_root,
  option_value("download-root", "raw_data/finngen_multiome")
)
output_path <- project_path(
  project_root,
  option_value(
    "output",
    paste0("outcome_data/targeted_susie/verification_", mode, "_", max_priority, ".tsv")
  )
)
checksum <- as_flag(args$checksum, default = FALSE)
allow_missing <- as_flag(args[["allow-missing"]], default = FALSE)

if (!mode %in% c("extract", "download")) stop("--mode must be extract or download.", call. = FALSE)
if (!max_priority %in% c("P0_critical", "P1_strong", "P2_multisignal")) {
  stop("Invalid --max-priority: ", max_priority, call. = FALSE)
}
if (!file.exists(plan_path)) stop("Missing acquisition plan: ", plan_path, call. = FALSE)

stage <- paste0("through_", max_priority)
plan <- data.table::fread(plan_path, sep = "\t", showProgress = FALSE)
plan <- plan[plan[["stage"]] == stage]
if (!nrow(plan)) stop("No acquisition-plan rows for: ", stage, call. = FALSE)

parse_regions <- function(x) {
  tokens <- strsplit(x, ";", fixed = TRUE)[[1]]
  parts <- regexec("^([^:]+):([0-9]+)-([0-9]+)$", tokens)
  matches <- regmatches(tokens, parts)
  valid <- lengths(matches) == 4L
  if (!all(valid)) stop("Malformed tabix region(s): ", paste(tokens[!valid], collapse = ", "), call. = FALSE)
  matches <- matches[valid]
  GenomicRanges::GRanges(
    seqnames = vapply(matches, `[[`, character(1), 2L),
    ranges = IRanges::IRanges(
      start = as.integer(vapply(matches, `[[`, character(1), 3L)),
      end = as.integer(vapply(matches, `[[`, character(1), 4L))
    )
  )
}

verify_one <- function(row) {
  source_path <- row$source_path[[1]]
  if (mode == "extract") {
    relative <- sub("\\.SUSIE\\.snp\\.tsv\\.gz$", ".SUSIE.snp.targeted.tsv.gz", source_path)
    data_path <- file.path(slice_root, stage, relative)
    expected_size <- NA_real_
    expected_md5 <- NA_character_
  } else {
    data_path <- file.path(download_root, source_path)
    expected_size <- as.numeric(row$size_bytes[[1]])
    expected_md5 <- row$md5[[1]]
  }
  index_path <- paste0(data_path, ".tbi")
  data_exists <- file.exists(data_path)
  index_exists <- file.exists(index_path)
  observed_size <- if (data_exists) file.info(data_path)$size else NA_real_
  size_matches <- if (mode == "download" && data_exists) {
    identical(as.numeric(observed_size), expected_size)
  } else if (data_exists) {
    observed_size > 0
  } else {
    FALSE
  }
  observed_md5 <- NA_character_
  md5_matches <- NA
  if (checksum && mode == "download" && data_exists) {
    observed_md5 <- unname(tools::md5sum(data_path))
    md5_matches <- identical(observed_md5, expected_md5)
  }

  index_ok <- FALSE
  indexed_seqnames <- NA_character_
  n_records_in_planned_regions <- NA_real_
  index_error <- NA_character_
  if (data_exists && index_exists) {
    tryCatch({
      tabix_file <- Rsamtools::TabixFile(data_path)
      seqs <- Rsamtools::seqnamesTabix(tabix_file)
      indexed_seqnames <- if (length(seqs)) paste(seqs, collapse = ";") else ""
      target_regions <- parse_regions(row$tabix_regions[[1]])
      counts <- Rsamtools::countTabix(tabix_file, param = target_regions)
      n_records_in_planned_regions <- sum(unlist(counts, use.names = FALSE))
      target_chromosomes <- unique(as.character(GenomicRanges::seqnames(target_regions)))
      index_ok <- all(target_chromosomes %in% seqs)
    }, error = function(e) {
      index_error <<- conditionMessage(e)
    })
  }

  status <- if (!data_exists) {
    "missing_data"
  } else if (!index_exists) {
    "missing_index"
  } else if (!size_matches) {
    "size_mismatch"
  } else if (isFALSE(md5_matches)) {
    "md5_mismatch"
  } else if (!index_ok) {
    "index_error"
  } else if (!is.na(n_records_in_planned_regions) && n_records_in_planned_regions == 0) {
    "no_records_in_planned_regions"
  } else {
    "ok"
  }

  data.table::data.table(
    stage = stage,
    mode = mode,
    source_path = source_path,
    local_data_path = data_path,
    local_index_path = index_path,
    data_exists = data_exists,
    index_exists = index_exists,
    expected_size_bytes = expected_size,
    observed_size_bytes = observed_size,
    size_matches = size_matches,
    checksum_requested = checksum,
    expected_md5 = expected_md5,
    observed_md5 = observed_md5,
    md5_matches = md5_matches,
    index_ok = index_ok,
    indexed_seqnames = indexed_seqnames,
    n_records_in_planned_regions = n_records_in_planned_regions,
    status = status,
    detail = index_error
  )
}

log_message("Verifying ", nrow(plan), " files for ", stage, " in ", mode, " mode")
verification <- data.table::rbindlist(lapply(seq_len(nrow(plan)), function(i) verify_one(plan[i])))
ensure_dir(dirname(output_path))
data.table::fwrite(verification, output_path, sep = "\t", quote = FALSE, na = "NA")
log_message("Wrote verification report: ", output_path)
print(verification[, .N, by = status][order(status)])

failed <- verification[status != "ok"]
if (nrow(failed) && !allow_missing) {
  stop(nrow(failed), " file(s) did not pass verification. Use --allow-missing only for inventory audits.", call. = FALSE)
}
