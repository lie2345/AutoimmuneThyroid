#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "codes/27_prepare_t09_batches.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- match(name, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}

plan_path <- get_arg("--plan", file.path(project_root, "outcome_data", "t08",
                                         "p0_core_non_hla_acquisition_plan.tsv"))
out_dir <- get_arg("--outdir", file.path(project_root, "outcome_data", "targeted_susie_t08",
                                         "non_hla", "batch_plans"))
workers <- as.integer(get_arg("--workers", "4"))
if (is.na(workers) || workers < 1L || workers > 8L) stop("--workers must be between 1 and 8")
if (!file.exists(plan_path)) stop("Missing plan: ", plan_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan <- read.delim(plan_path, check.names = FALSE, na.strings = c("", "NA"))
expected_header <- c("stage", "max_priority_rank", "modality", "level", "cell_type", "chrom",
                     "source_path", "data_url", "index_url", "size_bytes", "md5",
                     "index_size_bytes", "index_md5", "n_loci", "tabix_regions",
                     "selection_reasons", "recommended_mode")
if (!identical(names(plan), expected_header)) stop("Unexpected acquisition-plan header")
if (anyDuplicated(plan$source_path)) stop("Source paths must be unique before batching")

chrom_lengths <- c(chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
                   chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
                   chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
                   chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
                   chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
                   chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415)

interval_width <- function(x) {
  tokens <- strsplit(x, ";", fixed = TRUE)[[1]]
  bounds <- regmatches(tokens, regexec("^[^:]+:([0-9]+)-([0-9]+)$", tokens))
  if (any(lengths(bounds) != 3L)) stop("Malformed interval: ", x)
  sum(vapply(bounds, function(z) as.numeric(z[[3]]) - as.numeric(z[[2]]) + 1, numeric(1)))
}

plan$interval_bp <- vapply(plan$tabix_regions, interval_width, numeric(1))
plan$estimated_slice_weight <- plan$size_bytes * plan$interval_bp /
  unname(chrom_lengths[plan$chrom])
if (anyNA(plan$estimated_slice_weight)) stop("Missing chromosome length")

order_idx <- order(plan$estimated_slice_weight, decreasing = TRUE)
loads <- rep(0, workers)
assignment <- integer(nrow(plan))
for (idx in order_idx) {
  worker <- which.min(loads)
  assignment[[idx]] <- worker
  loads[[worker]] <- loads[[worker]] + plan$estimated_slice_weight[[idx]]
}
plan$batch <- assignment

for (batch_id in seq_len(workers)) {
  batch <- plan[plan$batch == batch_id, expected_header, drop = FALSE]
  batch <- batch[order(batch$modality, batch$level, batch$cell_type, batch$chrom), , drop = FALSE]
  path <- file.path(out_dir, sprintf("acquisition_plan_batch%d.tsv", batch_id))
  write.table(batch, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote ", path, " (", nrow(batch), " files)")
}

assignment_table <- plan[, c("batch", "source_path", "modality", "level", "cell_type", "chrom",
                             "n_loci", "interval_bp", "size_bytes", "estimated_slice_weight")]
assignment_table <- assignment_table[order(assignment_table$batch,
                                           -assignment_table$estimated_slice_weight), ]
write.table(assignment_table, file.path(out_dir, "batch_assignment.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

summary <- data.frame(
  batch = seq_len(workers),
  n_files = vapply(seq_len(workers), function(i) sum(plan$batch == i), integer(1)),
  estimated_slice_weight = loads,
  stringsAsFactors = FALSE
)
write.table(summary, file.path(out_dir, "batch_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
message("Prepared ", workers, " mutually exclusive batches for ", nrow(plan), " files.")
