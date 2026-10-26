#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages("data.table")

args <- commandArgs(trailingOnly = TRUE)
chr_arg <- sub("^--chr=", "", grep("^--chr=", args, value = TRUE))
chromosomes <- if (length(chr_arg)) as.integer(chr_arg) else 1:22
if (anyNA(chromosomes) || any(!chromosomes %in% 1:22)) stop("--chr must be 1-22.", call. = FALSE)

plink_dir <- file.path(project_root, "raw_data", "reference_panels", "ldsc", "1000G_EUR_Phase3_plink")
ann_dir <- file.path(project_root, "outcome_data", "t12", "annotations")
mapped_bed <- file.path(ann_dir, "consensus_peaks.grch37.bed")
metadata_path <- file.path(ann_dir, "consensus_peaks.grch38.metadata.tsv.gz")
if (!file.exists(mapped_bed) || !file.exists(metadata_path)) stop("Run T12 liftOver preparation first.", call. = FALSE)

lineage_map <- c(
  B_intermediate = "B", B_memory = "B", B_naive = "B", Plasmablast = "B",
  CD4_CTL = "T", CD4_Naive = "T", CD4_Proliferating = "T", CD4_TCM = "T",
  CD4_TEM = "T", CD8_Naive = "T", CD8_TCM = "T", CD8_TEM = "T",
  MAIT = "T", Treg = "T", dnT = "T", gdT = "T",
  CD14_Mono = "Mono", CD16_Mono = "Mono",
  cDC1 = "DC", cDC2 = "DC", pDC = "DC",
  NK = "NK", NK_CD56bright = "NK", NK_Proliferating = "NK",
  HSPC = "Other", ILC = "Other", Platelet = "Other", Eryth = "Other"
)
cell_types <- sort(names(lineage_map))
lineages <- sort(unique(unname(lineage_map)))
cell_cols <- paste0("cell__", cell_types)
lineage_cols <- paste0("lineage__", lineages)
annotation_cols <- c(cell_cols, lineage_cols, "control__AllConsensus")

metadata <- data.table::fread(metadata_path, select = c("peak_uid", "source_cell_type"), showProgress = FALSE)
if (!setequal(unique(metadata$source_cell_type), cell_types)) stop("Unexpected source cell-type labels.", call. = FALSE)
mapped <- data.table::fread(
  mapped_bed, header = FALSE,
  col.names = c("chrom", "start0", "end", "peak_uid", "score", "strand"),
  showProgress = FALSE
)
mapped <- mapped[, .(chrom, start0, end, peak_uid)]
mapped <- metadata[mapped, on = "peak_uid", nomatch = 0L]
if (nrow(mapped) != 297603L || anyDuplicated(mapped$peak_uid)) stop("Unexpected mapped peak inventory.", call. = FALSE)
mapped[, chrom_num := suppressWarnings(as.integer(sub("^chr", "", chrom)))]
mapped <- mapped[chrom_num %in% 1:22]
mapped[, `:=`(start = start0 + 1L, lineage = unname(lineage_map[source_cell_type]))]
if (anyNA(mapped$lineage)) stop("Unmapped lineage label.", call. = FALSE)

dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
validation_rows <- vector("list", length(chromosomes))
for (idx in seq_along(chromosomes)) {
  chr <- chromosomes[[idx]]
  bim_path <- file.path(plink_dir, paste0("1000G.EUR.QC.", chr, ".bim"))
  if (!file.exists(bim_path)) stop("Missing official PLINK BIM: ", bim_path, call. = FALSE)
  bim <- data.table::fread(
    bim_path, header = FALSE,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"), showProgress = FALSE
  )
  if (any(bim$CHR != chr) || any(diff(bim$BP) < 0)) stop("Unexpected BIM chromosome/order for chr", chr, call. = FALSE)

  p <- mapped[chrom_num == chr, .(chrom_num, start, end, source_cell_type, lineage)]
  data.table::setkey(p, chrom_num, start, end)
  q <- bim[, .(snp_index = .I, chrom_num = CHR, start = BP, end = BP)]
  hits <- data.table::foverlaps(q, p, by.x = c("chrom_num", "start", "end"), nomatch = 0L)

  mat <- matrix(0L, nrow = nrow(bim), ncol = length(annotation_cols),
                dimnames = list(NULL, annotation_cols))
  if (nrow(hits)) {
    mat[cbind(hits$snp_index, match(paste0("cell__", hits$source_cell_type), annotation_cols))] <- 1L
    mat[cbind(hits$snp_index, match(paste0("lineage__", hits$lineage), annotation_cols))] <- 1L
    mat[hits$snp_index, "control__AllConsensus"] <- 1L
  }
  mhc_mask <- chr == 6L & bim$BP >= 25000000L & bim$BP <= 34000000L
  if (any(mhc_mask)) mat[mhc_mask, ] <- 0L

  out_path <- file.path(ann_dir, paste0("t12_source_labels.", chr, ".annot.gz"))
  data.table::fwrite(data.table::as.data.table(mat), out_path, sep = "\t", quote = FALSE)
  validation_rows[[idx]] <- data.table::data.table(
    chromosome = chr,
    plink_snps = nrow(bim),
    mapped_peaks_on_chromosome = nrow(p),
    snps_overlapping_any_peak_non_mhc = sum(mat[, "control__AllConsensus"] == 1L),
    mhc_snps_zeroed = sum(mhc_mask),
    annotation_columns = ncol(mat),
    annotation_file = file.path("outcome_data", "t12", "annotations", basename(out_path))
  )
  log_message("T12 annotation chr", chr, ": ", nrow(bim), " SNPs; ",
              validation_rows[[idx]]$snps_overlapping_any_peak_non_mhc, " overlap source-labelled peaks.")
}

new_validation <- data.table::rbindlist(validation_rows)
validation_path <- file.path(ann_dir, "annotation_build_validation.tsv")
old_validation <- if (file.exists(validation_path)) data.table::fread(validation_path) else NULL
old_keep <- if (is.null(old_validation)) NULL else old_validation[!chromosome %in% chromosomes]
combined <- data.table::rbindlist(list(old_keep, new_validation), fill = TRUE)
data.table::setorder(combined, chromosome)
write_tsv(combined, validation_path)

map_table <- data.table::data.table(
  annotation = annotation_cols,
  class = c(rep("source_cell_type", length(cell_cols)), rep("lineage", length(lineage_cols)), "open_chromatin_control"),
  source_label = c(cell_types, lineages, "AllConsensus"),
  interpretation = c(
    rep("source-labelled consensus peak set; not a differential-accessibility contrast", length(cell_cols)),
    rep("union of source-labelled consensus peak sets within lineage", length(lineage_cols)),
    "union of all source-labelled consensus peaks"
  )
)
write_tsv(map_table, file.path(ann_dir, "annotation_dictionary.tsv"))
