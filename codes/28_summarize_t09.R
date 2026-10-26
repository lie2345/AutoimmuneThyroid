#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "codes/28_summarize_t09.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- match(name, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}

out_dir <- get_arg("--outdir", file.path(project_root, "outcome_data", "t09"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bgzip_bin <- file.path(project_root, "tools", "htslib-1.19-curl", "bin", "bgzip")
if (!file.exists(bgzip_bin)) {
  bgzip_bin <- Sys.which("bgzip")
}
if (!nzchar(bgzip_bin) || !file.exists(bgzip_bin)) stop("bgzip is required for integrity tests")

arms <- data.frame(
  arm = c("non_hla", "hla"),
  plan_path = c(
    file.path(project_root, "outcome_data", "t08", "p0_core_non_hla_acquisition_plan.tsv"),
    file.path(project_root, "outcome_data", "t08", "p0_core_hla_acquisition_plan.tsv")
  ),
  verification_path = c(
    file.path(project_root, "outcome_data", "targeted_susie_t08", "non_hla", "verification.tsv"),
    file.path(project_root, "outcome_data", "targeted_susie_t08", "hla", "verification.tsv")
  ),
  slice_root = c(
    file.path(project_root, "outcome_data", "targeted_susie_t08", "non_hla", "slices"),
    file.path(project_root, "outcome_data", "targeted_susie_t08", "hla", "slices")
  )
)

required_plan <- c(
  "stage", "max_priority_rank", "modality", "level", "cell_type", "chrom",
  "source_path", "data_url", "index_url", "size_bytes", "md5",
  "index_size_bytes", "index_md5", "n_loci", "tabix_regions",
  "selection_reasons", "recommended_mode"
)
required_verification <- c(
  "stage", "mode", "source_path", "local_data_path", "local_index_path",
  "data_exists", "index_exists", "observed_size_bytes", "index_ok",
  "n_records_in_planned_regions", "status"
)

md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

inventory_parts <- vector("list", nrow(arms))
partial_parts <- vector("list", nrow(arms))

for (i in seq_len(nrow(arms))) {
  arm <- arms$arm[[i]]
  plan_path <- arms$plan_path[[i]]
  verification_path <- arms$verification_path[[i]]
  if (!file.exists(plan_path)) stop("Missing plan: ", plan_path)
  if (!file.exists(verification_path)) stop("Missing verification: ", verification_path)

  plan <- read.delim(plan_path, check.names = FALSE, na.strings = c("", "NA"))
  verification <- read.delim(verification_path, check.names = FALSE, na.strings = c("", "NA"))
  if (!identical(names(plan), required_plan)) stop("Unexpected plan header: ", plan_path)
  missing_verification <- setdiff(required_verification, names(verification))
  if (length(missing_verification)) {
    stop("Verification is missing columns: ", paste(missing_verification, collapse = ", "))
  }
  if (anyDuplicated(plan$source_path) || anyDuplicated(verification$source_path)) {
    stop("Duplicate source_path in ", arm, " plan or verification")
  }

  merged <- merge(plan, verification, by = c("stage", "source_path"), all = TRUE,
                  sort = FALSE, suffixes = c("_plan", "_verification"))
  if (nrow(merged) != nrow(plan) || anyNA(merged$status)) {
    stop("Plan/verification mismatch for ", arm)
  }
  merged$arm <- arm
  merged$local_index_size_bytes <- vapply(
    merged$local_index_path,
    function(path) if (is.na(path) || !file.exists(path)) NA_real_ else file.info(path)$size,
    numeric(1)
  )
  merged$local_data_md5 <- vapply(merged$local_data_path, md5_or_na, character(1))
  merged$local_index_md5 <- vapply(merged$local_index_path, md5_or_na, character(1))
  merged$bgzip_integrity_ok <- vapply(merged$local_data_path, function(path) {
    if (is.na(path) || !file.exists(path)) return(FALSE)
    identical(system2(bgzip_bin, c("-t", path), stdout = FALSE, stderr = FALSE), 0L)
  }, logical(1))
  merged$n_planned_regions <- vapply(
    strsplit(merged$tabix_regions, ";", fixed = TRUE), length, integer(1)
  )

  keep <- c(
    "arm", "stage", "modality", "level", "cell_type", "chrom", "source_path",
    "data_url", "index_url", "size_bytes", "md5", "index_size_bytes", "index_md5",
    "n_loci", "n_planned_regions", "tabix_regions", "selection_reasons",
    "local_data_path", "local_index_path", "observed_size_bytes", "local_index_size_bytes",
    "local_data_md5", "local_index_md5", "bgzip_integrity_ok", "n_records_in_planned_regions",
    "data_exists", "index_exists", "index_ok", "status", "detail"
  )
  inventory_parts[[i]] <- merged[, keep, drop = FALSE]

  partial_paths <- list.files(arms$slice_root[[i]], pattern = "\\.part\\.", recursive = TRUE,
                              full.names = TRUE, all.files = TRUE)
  partial_parts[[i]] <- data.frame(
    arm = rep(arm, length(partial_paths)),
    path = partial_paths,
    size_bytes = if (length(partial_paths)) file.info(partial_paths)$size else numeric(),
    stringsAsFactors = FALSE
  )
}

inventory <- do.call(rbind, inventory_parts)
inventory <- inventory[order(inventory$arm, inventory$modality, inventory$level,
                             inventory$cell_type, inventory$chrom), ]
write.table(inventory, file.path(out_dir, "t09_file_inventory.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

partials <- do.call(rbind, partial_parts)
write.table(partials, file.path(out_dir, "t09_partial_files.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

log_rows <- list()
log_id <- 0L
for (i in seq_len(nrow(arms))) {
  log_dir <- dirname(arms$verification_path[[i]])
  log_paths <- list.files(log_dir, pattern = "^acquisition.*\\.tsv$", full.names = TRUE)
  chunk_log_dir <- file.path(log_dir, "chunk_transport", "logs")
  if (dir.exists(chunk_log_dir)) {
    log_paths <- c(
      log_paths,
      list.files(chunk_log_dir, pattern = "^source[0-9]+_chunk[0-9]+\\.tsv$",
                 full.names = TRUE)
    )
  }
  for (path in log_paths) {
    event <- read.delim(path, check.names = FALSE, na.strings = c("", "NA"))
    if (!nrow(event)) next
    log_id <- log_id + 1L
    event$arm <- arms$arm[[i]]
    event$log_file <- basename(path)
    log_rows[[log_id]] <- event
  }
}
if (length(log_rows)) {
  acquisition_history <- do.call(rbind, log_rows)
} else {
  acquisition_history <- data.frame()
}
write.table(acquisition_history, file.path(out_dir, "t09_acquisition_history.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

arm_summary <- do.call(rbind, lapply(split(inventory, inventory$arm), function(x) {
  regions <- unique(unlist(strsplit(x$tabix_regions, ";", fixed = TRUE)))
  data.frame(
    arm = x$arm[[1]],
    planned_source_files = nrow(x),
    unique_physical_regions = length(regions),
    source_region_targets = sum(x$n_planned_regions),
    verified_ok_files = sum(x$status == "ok"),
    bgzip_integrity_ok_files = sum(x$bgzip_integrity_ok),
    failed_files = sum(x$status != "ok"),
    targeted_data_bytes = sum(x$observed_size_bytes, na.rm = TRUE),
    targeted_index_bytes = sum(x$local_index_size_bytes, na.rm = TRUE),
    records_in_planned_regions = sum(x$n_records_in_planned_regions, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
arm_summary <- arm_summary[match(c("non_hla", "hla"), arm_summary$arm), ]
total <- data.frame(
  arm = "total",
  planned_source_files = sum(arm_summary$planned_source_files),
  unique_physical_regions = sum(arm_summary$unique_physical_regions),
  source_region_targets = sum(arm_summary$source_region_targets),
  verified_ok_files = sum(arm_summary$verified_ok_files),
  bgzip_integrity_ok_files = sum(arm_summary$bgzip_integrity_ok_files),
  failed_files = sum(arm_summary$failed_files),
  targeted_data_bytes = sum(arm_summary$targeted_data_bytes),
  targeted_index_bytes = sum(arm_summary$targeted_index_bytes),
  records_in_planned_regions = sum(arm_summary$records_in_planned_regions),
  stringsAsFactors = FALSE
)
summary <- rbind(arm_summary, total)
write.table(summary, file.path(out_dir, "t09_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

all_ok <- nrow(inventory) == 30L && all(inventory$status == "ok") &&
  all(inventory$bgzip_integrity_ok)
no_partials <- nrow(partials) == 0L
decision <- data.frame(
  task = "T09",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  expected_files = 30L,
  observed_files = nrow(inventory),
  verified_ok_files = sum(inventory$status == "ok"),
  bgzip_integrity_ok_files = sum(inventory$bgzip_integrity_ok),
  failed_files = sum(inventory$status != "ok"),
  partial_files = nrow(partials),
  decision = if (all_ok && no_partials) "PASS_GO_T10" else "FAIL_REMAIN_T09",
  interpretation = if (all_ok && no_partials) {
    "All frozen P0 targeted SuSiE slices and indexes passed the independent verification gate; T10 may start."
  } else {
    "At least one T09 acquisition or cleanup gate failed; do not start T10."
  },
  stringsAsFactors = FALSE
)
write.table(decision, file.path(out_dir, "t09_decision.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

print(summary)
print(decision)
if (!(all_ok && no_partials)) stop("T09 completion gate failed; see t09_decision.tsv")
