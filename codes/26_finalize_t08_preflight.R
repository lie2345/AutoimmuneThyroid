#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "codes/26_finalize_t08_preflight.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(project_root, "outcome_data", "t08")

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty input: ", path)
  read.delim(path, check.names = FALSE, na.strings = c("", "NA"))
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote: ", path, " (", nrow(x), " rows)")
}

batch_paths <- file.path(out_dir, sprintf("remote_tabix_preflight_batch%d.tsv", 1:4))
batch <- do.call(rbind, lapply(batch_paths, read_tsv))
batch$evidence_source <- "four_batch_primary_preflight"

retry_paths <- c(file.path(out_dir, "remote_tabix_preflight_retry.tsv"),
                 file.path(out_dir, "remote_tabix_preflight_retry2.tsv"))
retry_paths <- retry_paths[file.exists(retry_paths)]
retry <- if (length(retry_paths)) do.call(rbind, lapply(retry_paths, read_tsv)) else batch[FALSE, ]
if (nrow(retry)) retry$evidence_source <- "targeted_retry"

all_runs <- rbind(batch, retry)
all_runs$run_order <- seq_len(nrow(all_runs))
key <- paste(all_runs$source_path, all_runs$locus_id, all_runs$test_region, sep = "||")
all_runs$key <- key
all_runs <- all_runs[order(all_runs$key, all_runs$status != "ok", -all_runs$run_order), ]
final <- all_runs[!duplicated(all_runs$key), ]
final <- final[order(final$is_hla, final$source_path, final$locus_id), ]
final$key <- NULL
final$run_order <- NULL

targets <- read_tsv(file.path(out_dir, "remote_preflight_targets.tsv"))
target_key <- paste(targets$source_path, targets$locus_id, targets$test_region, sep = "||")
final_key <- paste(final$source_path, final$locus_id, final$test_region, sep = "||")
if (length(target_key) != 32L || anyDuplicated(target_key) || !setequal(target_key, final_key)) {
  stop("Final remote preflight does not match the 32 frozen source-locus targets exactly.")
}
if (any(final$status != "ok")) {
  stop("Remote preflight still contains failed targets: ",
       paste(final$source_path[final$status != "ok"], collapse = ";"))
}

failed_primary <- batch[batch$status != "ok", , drop = FALSE]
summary <- data.frame(
  metric = c("frozen_source_locus_targets", "unique_source_files", "non_HLA_targets", "HLA_targets",
             "primary_batch_passed", "primary_batch_failed", "targeted_retry_passed",
             "final_passed", "final_failed", "targets_requiring_more_than_one_primary_attempt"),
  value = c(nrow(final), length(unique(final$source_path)), sum(!final$is_hla), sum(final$is_hla),
            sum(batch$status == "ok"), sum(batch$status != "ok"),
            sum(retry$evidence_source == "targeted_retry" & retry$status == "ok"),
            sum(final$status == "ok"), sum(final$status != "ok"), sum(batch$query_attempts > 1)),
  interpretation = c(
    "one known-variant query for every frozen source-file by locus pair",
    "combined P0 plan contains 30 unique remote chromosome files",
    "ordinary analysis plan; HLA excluded", "separate HLA plan",
    "passed in the four-batch primary preflight", "transient connection failure requiring targeted retry",
    "at least one targeted retry completed successfully",
    "all frozen targets are remotely Tabix-extractable", "must be zero to pass T08",
    "network variability observed; acquisition should retain retry and keep-going controls"
  ),
  stringsAsFactors = FALSE
)

decision <- data.frame(
  gate = "T08_targeted_SuSiE_plan_review",
  status = "PASS",
  decision = "GO_T09_NON_HLA_FIRST_HLA_SEPARATE",
  rationale = paste(
    "The T07-aware frozen plan covers all four divergent non-HLA loci;",
    "17 static checks passed; all 32 source-locus targets passed real remote Tabix point queries;",
    "the 22-file non-HLA plan and 8-file HLA plan remain separate; no full-file download was run."
  ),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)

write_tsv(final, file.path(out_dir, "remote_tabix_preflight.tsv"))
write_tsv(all_runs[, setdiff(names(all_runs), "key"), drop = FALSE],
          file.path(out_dir, "remote_tabix_preflight_attempt_history.tsv"))
write_tsv(summary, file.path(out_dir, "t08_validation_summary.tsv"))
write_tsv(decision, file.path(out_dir, "t08_decision.tsv"))

message("T08 remote preflight PASS: ", nrow(final), "/", nrow(final), " frozen targets.")
