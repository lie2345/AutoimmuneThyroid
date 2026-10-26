#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "codes/29_prepare_t09_hla_chunks.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- match(name, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}

plan_path <- get_arg("--plan", file.path(project_root, "outcome_data", "t08",
                                          "p0_core_hla_acquisition_plan.tsv"))
out_dir <- get_arg("--outdir", file.path(project_root, "outcome_data",
                                          "targeted_susie_t08", "hla", "chunk_transport"))
n_chunks <- as.integer(get_arg("--chunks", "4"))
if (is.na(n_chunks) || n_chunks < 2L || n_chunks > 16L) stop("--chunks must be between 2 and 16")
if (!file.exists(plan_path)) stop("Missing plan: ", plan_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

plan <- read.delim(plan_path, check.names = FALSE, na.strings = c("", "NA"))
expected_header <- c("stage", "max_priority_rank", "modality", "level", "cell_type", "chrom",
                     "source_path", "data_url", "index_url", "size_bytes", "md5",
                     "index_size_bytes", "index_md5", "n_loci", "tabix_regions",
                     "selection_reasons", "recommended_mode")
if (!identical(names(plan), expected_header)) stop("Unexpected acquisition-plan header")
if (anyDuplicated(plan$source_path)) stop("HLA source paths must be unique")
if (any(plan$n_loci != 1L)) stop("Chunk transport currently requires one region per source")

parse_region <- function(x) {
  m <- regmatches(x, regexec("^([^:]+):([0-9]+)-([0-9]+)$", x))[[1]]
  if (length(m) != 4L) stop("Malformed region: ", x)
  list(chrom = m[[2]], start = as.integer(m[[3]]), end = as.integer(m[[4]]))
}

manifest_rows <- vector("list", nrow(plan) * n_chunks)
row_id <- 0L
for (source_id in seq_len(nrow(plan))) {
  region <- parse_region(plan$tabix_regions[[source_id]])
  width <- region$end - region$start + 1L
  bounds <- region$start + floor((0:n_chunks) * width / n_chunks)
  chunk_starts <- bounds[seq_len(n_chunks)]
  chunk_ends <- bounds[seq_len(n_chunks) + 1L] - 1L
  chunk_ends[[n_chunks]] <- region$end
  if (chunk_starts[[1]] != region$start || chunk_ends[[n_chunks]] != region$end ||
      any(chunk_starts[-1L] != chunk_ends[-n_chunks] + 1L)) {
    stop("Chunk union failed for source ", source_id)
  }

  for (chunk_id in seq_len(n_chunks)) {
    row_id <- row_id + 1L
    chunk_region <- sprintf("%s:%d-%d", region$chrom, chunk_starts[[chunk_id]],
                            chunk_ends[[chunk_id]])
    job_id <- sprintf("source%02d_chunk%02d", source_id, chunk_id)
    job_dir <- file.path(out_dir, "jobs", job_id)
    slice_root <- file.path(job_dir, "slices")
    log_path <- file.path(out_dir, "logs", paste0(job_id, ".tsv"))
    plan_out <- file.path(out_dir, "plans", paste0(job_id, ".tsv"))
    dir.create(dirname(plan_out), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    dir.create(slice_root, recursive = TRUE, showWarnings = FALSE)

    chunk_plan <- plan[source_id, , drop = FALSE]
    chunk_plan$tabix_regions <- chunk_region
    chunk_plan$n_loci <- 1L
    chunk_plan$selection_reasons <- paste0(
      chunk_plan$selection_reasons, ";transport_chunk_", chunk_id, "_of_", n_chunks
    )
    write.table(chunk_plan, plan_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    relative <- sub("\\.SUSIE\\.snp\\.tsv\\.gz$", ".SUSIE.snp.targeted.tsv.gz",
                    chunk_plan$source_path)
    chunk_data_path <- file.path(slice_root, chunk_plan$stage, relative)
    manifest_rows[[row_id]] <- data.frame(
      source_id = source_id,
      chunk_id = chunk_id,
      n_chunks = n_chunks,
      source_path = chunk_plan$source_path,
      original_region = plan$tabix_regions[[source_id]],
      chunk_region = chunk_region,
      plan_path = plan_out,
      slice_root = slice_root,
      chunk_data_path = chunk_data_path,
      log_path = log_path,
      stringsAsFactors = FALSE
    )
  }
}

manifest <- do.call(rbind, manifest_rows)
manifest <- manifest[order(manifest$chunk_id, manifest$source_id), ]
manifest_path <- file.path(out_dir, "chunk_manifest.tsv")
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Prepared", nrow(manifest), "mutually exclusive HLA transport chunks in", out_dir, "\n")
cat("Manifest:", manifest_path, "\n")
