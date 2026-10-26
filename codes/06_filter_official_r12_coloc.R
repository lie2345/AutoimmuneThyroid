#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)[include_locus_analysis == TRUE]
args <- parse_cli()
requested <- split_csv(args$traits %||% "")
if (length(requested)) traits <- traits[trait_id %in% requested]

out <- list()
for (modality in c("eQTL", "caQTL")) {
  path <- file.path(
    project_root, "raw_data", "finngen_multiome", "coloc",
    paste0("finngen_multiome_v1.coloc.FinnGen_R12.", modality, ".tsv.gz")
  )
  if (!file.exists(path)) stop("Missing official R12 coloc table: ", path, call. = FALSE)
  log_message("Reading official R12 ", modality, " colocalization table")
  x <- data.table::fread(path, showProgress = TRUE)
  x <- x[
    trait1 %in% traits$trait_id &
      PP.H4.abf > 0.8 &
      cs1_log10bf > 0.9 & cs2_log10bf > 0.9 &
      probmass_1 > 0.5 & probmass_2 > 0.5 & cs_overlap > 0
  ]
  if (!nrow(x)) next
  x[, `:=`(
    modality = modality,
    level = sub(".*predicted\\.celltype\\.(l[12])\\..*", "\\1", dataset2),
    cell_type = sub(".*predicted\\.celltype\\.l[12]\\.([^.]+)\\..*", "\\1", dataset2),
    source_release = "FinnGen_R12_official"
  )]
  out[[length(out) + 1L]] <- x
}

result <- if (length(out)) data.table::rbindlist(out, fill = TRUE) else data.table::data.table()
write_tsv(result, file.path(project_root, cfg$project$outcome_dir, "coloc_official_r12", "thyroid_strong_coloc.tsv.gz"))
log_message("Retained ", nrow(result), " official R12 strong colocalizations as a benchmark, not as a substitute for R13 analysis.")
