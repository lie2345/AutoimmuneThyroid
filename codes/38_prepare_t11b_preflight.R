#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages("data.table")

rows_path <- file.path(project_root, "outcome_data", "t11a", "p1a_wave2_candidate_rows.tsv.gz")
plan_path <- file.path(project_root, "outcome_data", "t11a", "p1a_wave2_acquisition_plan.tsv")
out_dir <- file.path(project_root, "outcome_data", "t11b")
out_path <- file.path(out_dir, "remote_preflight_targets.tsv")
if (!file.exists(rows_path) || !file.exists(plan_path)) stop("Missing frozen T11B inputs.", call. = FALSE)

rows <- data.table::fread(rows_path, showProgress = FALSE)
plan <- data.table::fread(plan_path, showProgress = FALSE)
rows[, test_position := suppressWarnings(as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", top_h4_variant)))]
rows[, test_region := paste0(manifest_chrom, ":", test_position, "-", test_position)]
rows[, in_planned_region := test_position >= start & test_position <= end]
if (anyNA(rows$test_position) || any(!rows$in_planned_region)) {
  stop("At least one selected regional top variant is outside its planned interval.", call. = FALSE)
}

best <- rows[order(-PP.H4.abf), .SD[1], by = source_path]
targets <- merge(
  best[, .(
    source_path, modality, level, cell_type, locus_id, is_hla,
    chrom = manifest_chrom, top_h4_variant, PP.H4.abf, test_position, test_region
  )],
  plan[, .(source_path, data_url, index_url)],
  by = "source_path", all.x = TRUE, sort = FALSE
)
data.table::setcolorder(targets, c(
  "source_path", "modality", "level", "cell_type", "locus_id", "is_hla",
  "data_url", "index_url", "chrom", "top_h4_variant", "PP.H4.abf",
  "test_position", "test_region"
))
data.table::setorder(targets, modality, level, cell_type, chrom)
if (nrow(targets) != 9L || nrow(plan) != 9L || anyNA(targets$data_url) ||
    anyDuplicated(targets$source_path) || any(targets$is_hla)) {
  stop("T11B preflight targets do not map one-to-one to the frozen 9-file non-HLA plan.", call. = FALSE)
}
ensure_dir(out_dir)
write_tsv(targets, out_path)
log_message("Prepared 9 T11B wave-2 preflight targets: ", out_path)
