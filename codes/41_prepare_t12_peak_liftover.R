#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages("data.table")

source_path <- file.path(
  project_root, "raw_data", "finngen_multiome", "atac_peaks",
  "finngen_multiome_v1.consensus_peaks.bed.gz"
)
out_dir <- file.path(project_root, "outcome_data", "t12", "annotations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

peaks <- data.table::fread(
  source_path, header = FALSE,
  col.names = c("chrom", "start0", "end", "source_peak_names", "read_count", "strand"),
  showProgress = FALSE
)
peaks[, source_cell_type := vapply(strsplit(source_peak_names, ",", fixed = TRUE), function(tokens) {
  values <- sort(unique(sub("_peak_.*$", "", tokens)))
  if (length(values) != 1L) stop("Peak has multiple source cell-type labels.", call. = FALSE)
  sub("^predicted[.]celltype[.]l2[.]", "", values)
}, character(1))]
peaks[, peak_uid := sprintf("T12PK%06d", .I)]

metadata <- peaks[, .(
  peak_uid, source_cell_type, source_peak_names, grch38_chrom = chrom,
  grch38_start0 = start0, grch38_end = end, read_count
)]
write_tsv(metadata, file.path(out_dir, "consensus_peaks.grch38.metadata.tsv.gz"))

bed <- peaks[, .(chrom, start0, end, peak_uid, score = 0L, strand)]
data.table::fwrite(
  bed, file.path(out_dir, "consensus_peaks.grch38.bed"),
  sep = "\t", col.names = FALSE, quote = FALSE
)
log_message("Prepared ", nrow(bed), " GRCh38 peaks for UCSC liftOver.")
