#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages("data.table")

t12_dir <- file.path(project_root, "outcome_data", "t12")
master_dir <- file.path(t12_dir, "ldscores", "master")
split_dir <- file.path(t12_dir, "ldscores", "by_annotation")
model_dir <- file.path(t12_dir, "models")
dict <- data.table::fread(file.path(t12_dir, "annotations", "annotation_dictionary.tsv"))
annotations <- dict$annotation
ld_columns <- paste0(annotations, "L2")
dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

n_total <- 0
sum_x <- setNames(rep(0, length(annotations)), annotations)
cross_x <- matrix(0, length(annotations), length(annotations), dimnames = list(annotations, annotations))
validation <- vector("list", 22)

for (chr in 1:22) {
  master_prefix <- file.path(master_dir, paste0("t12_source_labels.", chr))
  ld_path <- paste0(master_prefix, ".l2.ldscore.gz")
  m_path <- paste0(master_prefix, ".l2.M")
  m550_path <- paste0(master_prefix, ".l2.M_5_50")
  if (!all(file.exists(c(ld_path, m_path, m550_path)))) stop("Missing master LD scores for chr", chr, call. = FALSE)
  ld <- data.table::fread(ld_path, showProgress = FALSE)
  if (!all(c("CHR", "SNP", "BP", ld_columns) %in% names(ld))) stop("Unexpected LD-score columns for chr", chr, call. = FALSE)
  m <- scan(m_path, quiet = TRUE)
  m550 <- scan(m550_path, quiet = TRUE)
  if (length(m) != length(annotations) || length(m550) != length(annotations)) stop("Unexpected M vector length for chr", chr, call. = FALSE)

  x <- as.matrix(ld[, ..ld_columns])
  colnames(x) <- annotations
  n_total <- n_total + nrow(x)
  sum_x <- sum_x + colSums(x)
  cross_x <- cross_x + crossprod(x)

  for (j in seq_along(annotations)) {
    annotation <- annotations[[j]]
    annotation_dir <- file.path(split_dir, annotation)
    dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
    out_prefix <- file.path(annotation_dir, paste0("t12.", annotation, ".", chr))
    out_ld <- ld[, .(CHR, SNP, BP, value = get(ld_columns[[j]]))]
    data.table::setnames(out_ld, "value", ld_columns[[j]])
    data.table::fwrite(out_ld, paste0(out_prefix, ".l2.ldscore.gz"), sep = "\t", quote = FALSE)
    writeLines(format(m[[j]], scientific = FALSE, trim = TRUE), paste0(out_prefix, ".l2.M"))
    writeLines(format(m550[[j]], scientific = FALSE, trim = TRUE), paste0(out_prefix, ".l2.M_5_50"))
  }
  validation[[chr]] <- data.table::data.table(
    chromosome = chr, hm3_ldscore_snps = nrow(ld), annotation_ldscores = length(annotations),
    split_files_expected = length(annotations) * 3L,
    split_files_observed = sum(vapply(annotations, function(a) {
      prefix <- file.path(split_dir, a, paste0("t12.", a, ".", chr))
      all(file.exists(paste0(prefix, c(".l2.ldscore.gz", ".l2.M", ".l2.M_5_50"))))
    }, logical(1))) * 3L
  )
  log_message("Split and validated T12 LD scores for chr", chr, ".")
}

cov_x <- (cross_x - tcrossprod(sum_x) / n_total) / (n_total - 1)
sd_x <- sqrt(diag(cov_x))
cor_x <- cov_x / outer(sd_x, sd_x)
cor_x[!is.finite(cor_x)] <- NA_real_
cor_out <- data.table::as.data.table(cor_x, keep.rownames = "annotation")
write_tsv(cor_out, file.path(t12_dir, "ldscore_correlation_matrix.tsv"))
write_tsv(data.table::rbindlist(validation), file.path(t12_dir, "ldscore_split_validation.tsv"))

prefix_for <- function(annotation) {
  file.path(split_dir, annotation, paste0("t12.", annotation, "."))
}
all_control <- "control__AllConsensus"
cell_rows <- dict[class == "source_cell_type"]
lineage_rows <- dict[class == "lineage"]

global_cell <- cell_rows[, .(
  Name = source_label,
  prefixes = paste(prefix_for(annotation), prefix_for(all_control), sep = ",")
)]
global_lineage <- lineage_rows[, .(
  Name = source_label,
  prefixes = paste(prefix_for(annotation), prefix_for(all_control), sep = ",")
)]
lineage_by_cell <- c(
  B_intermediate = "B", B_memory = "B", B_naive = "B", Plasmablast = "B",
  CD4_CTL = "T", CD4_Naive = "T", CD4_Proliferating = "T", CD4_TCM = "T",
  CD4_TEM = "T", CD8_Naive = "T", CD8_TCM = "T", CD8_TEM = "T",
  MAIT = "T", Treg = "T", dnT = "T", gdT = "T",
  CD14_Mono = "Mono", CD16_Mono = "Mono", cDC1 = "DC", cDC2 = "DC", pDC = "DC",
  NK = "NK", NK_CD56bright = "NK", NK_Proliferating = "NK",
  HSPC = "Other", ILC = "Other", Platelet = "Other", Eryth = "Other"
)
within_lineage <- cell_rows[, .(
  Name = source_label,
  prefixes = paste(
    prefix_for(annotation),
    prefix_for(paste0("lineage__", unname(lineage_by_cell[source_label]))),
    prefix_for(all_control), sep = ","
  )
)]

data.table::fwrite(global_cell, file.path(model_dir, "global_cell.ldcts"), sep = "\t", col.names = FALSE, quote = FALSE)
data.table::fwrite(global_lineage, file.path(model_dir, "global_lineage.ldcts"), sep = "\t", col.names = FALSE, quote = FALSE)
data.table::fwrite(within_lineage, file.path(model_dir, "within_lineage.ldcts"), sep = "\t", col.names = FALSE, quote = FALSE)
design <- data.table::data.table(
  model = c("global_cell", "global_lineage", "within_lineage"),
  tested_annotations = c(nrow(global_cell), nrow(global_lineage), nrow(within_lineage)),
  custom_controls = c("AllConsensus", "AllConsensus", "matched_lineage+AllConsensus"),
  shared_covariates = "baselineLD_v2.2",
  primary_use = c(
    "cell-source association beyond general consensus-peak annotation",
    "lineage association beyond general consensus-peak annotation",
    "cell-source association conditional on its lineage and general consensus-peak annotation"
  )
)
write_tsv(design, file.path(model_dir, "model_design.tsv"))
log_message("Prepared split LD scores, correlation matrix, and three T12 conditional model families.")
