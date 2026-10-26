#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
args <- parse_cli()
cfg <- read_analysis_config(project_root)

option_value <- function(name, default) {
  value <- args[[name]]
  if (is.null(value)) default else value
}

collapse_unique <- function(x) {
  values <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(values)) paste(values, collapse = ";") else "NA"
}

collapse_reason_tokens <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("NA")
  tokens <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  collapse_unique(tokens)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

coloc_dir <- project_path(
  project_root,
  option_value("coloc_dir", file.path(cfg$project$outcome_dir, "coloc_r13"))
)
manifest_path <- project_path(
  project_root,
  option_value("manifest", cfg$targeted_susie$manifest)
)
output_dir <- project_path(
  project_root,
  option_value("output_dir", cfg$targeted_susie$output_dir)
)
strong_pp <- as.numeric(option_value("strong_pp", cfg$coloc$pp_h4_strong))
report_pp <- as.numeric(option_value("report_pp", cfg$coloc$pp_h4_report))
qtl_cs_min <- as.integer(option_value("qtl_cs_min", cfg$targeted_susie$qtl_credible_sets_min))
padding_bp <- as.integer(option_value("padding_bp", cfg$targeted_susie$region_padding_bp))
base_url <- sub("/+$", "", option_value("base_url", cfg$targeted_susie$public_base_url))

if (!is.finite(strong_pp) || strong_pp <= 0 || strong_pp > 1) {
  stop("--strong_pp must be in (0, 1].", call. = FALSE)
}
if (!is.finite(report_pp) || report_pp <= 0 || report_pp > strong_pp) {
  stop("--report_pp must be in (0, strong_pp].", call. = FALSE)
}
if (is.na(qtl_cs_min) || qtl_cs_min < 2L) {
  stop("--qtl_cs_min must be at least 2.", call. = FALSE)
}
if (is.na(padding_bp) || padding_bp < 0L) {
  stop("--padding_bp must be a non-negative integer.", call. = FALSE)
}
if (!file.exists(manifest_path)) stop("Missing Multiome manifest: ", manifest_path, call. = FALSE)

coloc_files <- list.files(coloc_dir, pattern = "\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
if (!length(coloc_files)) {
  stop("No Step 07 outputs found under: ", coloc_dir, call. = FALSE)
}

log_message("Reading ", length(coloc_files), " Step 07 result files")
coloc <- data.table::rbindlist(
  lapply(coloc_files, function(path) {
    x <- data.table::fread(path, showProgress = FALSE)
    x[, source_coloc_file := sub(paste0("^", coloc_dir, "/?"), "", path)]
    x
  }),
  fill = TRUE,
  use.names = TRUE
)

required_coloc <- c(
  "trait_id", "locus_id", "chrom", "start", "end", "is_hla", "modality",
  "level", "cell_type", "feature_id", "PP.H4.abf", "effect_direction",
  "gwas_susie_cs", "qtl_n_credible_sets"
)
missing_coloc <- setdiff(required_coloc, names(coloc))
if (length(missing_coloc)) {
  stop("Step 07 outputs lack columns: ", paste(missing_coloc, collapse = ", "), call. = FALSE)
}
if (!nrow(coloc)) stop("Step 07 outputs contain no rows.", call. = FALSE)

graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
primary_traits <- c(graves_id, aiht_id)
coloc <- coloc[trait_id %in% primary_traits & is.finite(PP.H4.abf) & PP.H4.abf >= report_pp]

if (length(split_csv(args$modalities))) coloc <- coloc[modality %in% split_csv(args$modalities)]
if (length(split_csv(args$cell_types))) {
  requested_cells <- split_csv(args$cell_types)
  coloc <- coloc[paste(level, cell_type, sep = ".") %in% requested_cells]
}
if (length(split_csv(args$loci))) coloc <- coloc[locus_id %in% split_csv(args$loci)]
if (!nrow(coloc)) stop("No Step 07 rows remain after filters.", call. = FALSE)

coloc[, `:=`(
  chrom = sub("^chr", "", as.character(chrom), ignore.case = TRUE),
  start = as.integer(start),
  end = as.integer(end),
  is_hla = as.logical(is_hla),
  PP.H4.abf = as.numeric(PP.H4.abf),
  gwas_susie_cs = suppressWarnings(as.integer(gwas_susie_cs)),
  qtl_n_credible_sets = suppressWarnings(as.integer(qtl_n_credible_sets))
)]

feature_keys <- c("locus_id", "modality", "level", "cell_type", "feature_id")
direction_summary <- coloc[, .(
  graves_increase = any(
    trait_id == graves_id & effect_direction == "disease_risk_allele_increases_feature",
    na.rm = TRUE
  ),
  graves_decrease = any(
    trait_id == graves_id & effect_direction == "disease_risk_allele_decreases_feature",
    na.rm = TRUE
  ),
  aiht_increase = any(
    trait_id == aiht_id & effect_direction == "disease_risk_allele_increases_feature",
    na.rm = TRUE
  ),
  aiht_decrease = any(
    trait_id == aiht_id & effect_direction == "disease_risk_allele_decreases_feature",
    na.rm = TRUE
  ),
  graves_max_pp_h4 = safe_max(PP.H4.abf[trait_id == graves_id]),
  aiht_max_pp_h4 = safe_max(PP.H4.abf[trait_id == aiht_id])
), by = feature_keys]
direction_summary[, opposite_cross_trait_direction :=
  (graves_increase & aiht_decrease) | (graves_decrease & aiht_increase)]
direction_summary[, strong_without_strong_counterpart :=
  xor(
    !is.na(graves_max_pp_h4) & graves_max_pp_h4 >= strong_pp,
    !is.na(aiht_max_pp_h4) & aiht_max_pp_h4 >= strong_pp
  )]

gwas_cs_summary <- coloc[
  !is.na(gwas_susie_cs) & gwas_susie_cs > 0L,
  .(gwas_credible_sets_observed = data.table::uniqueN(gwas_susie_cs)),
  by = .(trait_id, locus_id)
]

coloc <- merge(coloc, direction_summary, by = feature_keys, all.x = TRUE, sort = FALSE)
coloc <- merge(coloc, gwas_cs_summary, by = c("trait_id", "locus_id"), all.x = TRUE, sort = FALSE)
coloc[is.na(gwas_credible_sets_observed), gwas_credible_sets_observed := 0L]

coloc[, `:=`(
  reason_hla = is_hla == TRUE,
  reason_strong_coloc = PP.H4.abf >= strong_pp,
  reason_qtl_multisignal = !is.na(qtl_n_credible_sets) & qtl_n_credible_sets >= qtl_cs_min,
  reason_gwas_multisignal = gwas_credible_sets_observed >= 2L,
  reason_opposite_cross_trait_direction = opposite_cross_trait_direction == TRUE,
  note_strong_without_strong_counterpart = strong_without_strong_counterpart == TRUE
)]

reason_columns <- c(
  "reason_hla", "reason_strong_coloc", "reason_qtl_multisignal",
  "reason_gwas_multisignal", "reason_opposite_cross_trait_direction"
)
coloc[, selected := Reduce(`|`, .SD), .SDcols = reason_columns]
selected <- coloc[selected == TRUE]
if (!nrow(selected)) stop("No rows meet the targeted SuSiE selection rules.", call. = FALSE)

reason_labels <- c(
  reason_hla = "HLA_complex_LD",
  reason_strong_coloc = "strong_PP.H4",
  reason_qtl_multisignal = "QTL_multiple_credible_sets",
  reason_gwas_multisignal = "GWAS_multiple_CS_observed_in_reported_hits",
  reason_opposite_cross_trait_direction = "opposite_cross_trait_regulatory_direction"
)
reason_matrix <- as.matrix(selected[, ..reason_columns])
selected[, selection_reasons := vapply(seq_len(.N), function(i) {
  labels <- reason_labels[reason_matrix[i, ]]
  if (note_strong_without_strong_counterpart[[i]]) {
    labels <- c(labels, "strong_without_strong_counterpart_in_reported_results")
  }
  paste(unique(labels), collapse = ";")
}, character(1))]

selected[, priority_rank := data.table::fcase(
  reason_hla | reason_opposite_cross_trait_direction, 0L,
  reason_strong_coloc, 1L,
  default = 2L
)]
priority_labels <- c("P0_critical", "P1_strong", "P2_multisignal")
selected[, priority := priority_labels[priority_rank + 1L]]
selected[, manifest_chrom := paste0("chr", chrom)]

manifest <- data.table::fread(manifest_path, sep = "\t", showProgress = FALSE)
required_manifest <- c("path", "product", "modality", "level", "celltype", "chrom", "analysis", "size_bytes", "md5")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest)) {
  stop("Multiome manifest lacks columns: ", paste(missing_manifest, collapse = ", "), call. = FALSE)
}

manifest_data <- manifest[
  product == "finemapping" & analysis == "SUSIE.snp" & !endsWith(path, ".tbi")
]
manifest_index <- manifest[
  product == "finemapping" & analysis == "SUSIE.snp" & endsWith(path, ".tbi"),
  .(index_path = path, index_size_bytes = as.numeric(size_bytes), index_md5 = md5)
]
manifest_data[, `:=`(
  size_bytes = as.numeric(size_bytes),
  index_path = paste0(path, ".tbi")
)]
manifest_data <- merge(manifest_data, manifest_index, by = "index_path", all.x = TRUE, sort = FALSE)
manifest_data[, data_url := paste0(base_url, "/", path)]
manifest_data[, index_url := paste0(base_url, "/", index_path)]

mapped <- merge(
  selected,
  manifest_data[, .(
    modality, level, cell_type = celltype, manifest_chrom = chrom,
    source_path = path, data_url, index_url, size_bytes, md5,
    index_size_bytes, index_md5
  )],
  by = c("modality", "level", "cell_type", "manifest_chrom"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(mapped$source_path)) {
  missing_map <- unique(mapped[is.na(source_path), .(modality, level, cell_type, manifest_chrom)])
  stop(
    "Manifest mapping failed for: ",
    paste(apply(missing_map, 1, paste, collapse = "/"), collapse = ", "),
    call. = FALSE
  )
}
if (anyNA(mapped$index_md5)) {
  stop("One or more targeted SUSIE.snp files lack a matching .tbi manifest row.", call. = FALSE)
}

mapped[, `:=`(
  region_start = pmax(1L, start - padding_bp),
  region_end = end + padding_bp
)]
mapped[, tabix_region := paste0(manifest_chrom, ":", region_start, "-", region_end)]

region_plan <- mapped[, .(
  priority_rank = min(priority_rank),
  priority = priority_labels[min(priority_rank) + 1L],
  chrom = unique(manifest_chrom)[1],
  start = min(region_start),
  end = max(region_end),
  tabix_region = unique(tabix_region)[1],
  n_selected_rows = .N,
  n_features = data.table::uniqueN(feature_id),
  traits = collapse_unique(trait_id),
  features = collapse_unique(feature_id),
  selection_reasons = collapse_reason_tokens(selection_reasons)
), by = .(source_path, modality, level, cell_type, locus_id)]
data.table::setorder(region_plan, priority_rank, modality, level, cell_type, chrom, start)

file_plan <- mapped[, .(
  priority_rank = min(priority_rank),
  priority = priority_labels[min(priority_rank) + 1L],
  chrom = unique(manifest_chrom)[1],
  data_url = unique(data_url)[1],
  index_url = unique(index_url)[1],
  size_bytes = unique(size_bytes)[1],
  size_gib = unique(size_bytes)[1] / 1024^3,
  md5 = unique(md5)[1],
  index_size_bytes = unique(index_size_bytes)[1],
  index_md5 = unique(index_md5)[1],
  n_loci = data.table::uniqueN(locus_id),
  n_features = data.table::uniqueN(feature_id),
  n_selected_rows = .N,
  traits = collapse_unique(trait_id),
  loci = collapse_unique(locus_id),
  selection_reasons = collapse_reason_tokens(selection_reasons)
), by = .(source_path, modality, level, cell_type)]

region_strings <- region_plan[, .(
  tabix_regions = paste(unique(tabix_region[order(start, end)]), collapse = ";")
), by = source_path]
file_plan <- merge(file_plan, region_strings, by = "source_path", all.x = TRUE, sort = FALSE)
file_plan[, recommended_mode := "remote_extract"]
data.table::setcolorder(file_plan, c(
  "priority", "priority_rank", "modality", "level", "cell_type", "chrom",
  "source_path", "data_url", "index_url", "size_bytes", "size_gib", "md5",
  "index_size_bytes", "index_md5", "n_loci", "n_features", "n_selected_rows",
  "traits", "loci", "selection_reasons", "tabix_regions", "recommended_mode"
))
data.table::setorder(file_plan, priority_rank, modality, level, cell_type, chrom)

acquisition_plan <- data.table::rbindlist(lapply(0:2, function(rank) {
  stage_label <- paste0("through_", priority_labels[rank + 1L])
  stage_regions <- region_plan[priority_rank <= rank, .(
    n_loci = data.table::uniqueN(locus_id),
    tabix_regions = paste(unique(tabix_region[order(start, end)]), collapse = ";"),
    selection_reasons = collapse_reason_tokens(selection_reasons)
  ), by = source_path]
  stage_files <- file_plan[priority_rank <= rank]
  ans <- merge(
    stage_files[, .(
      modality, level, cell_type, chrom, source_path, data_url, index_url,
      size_bytes, md5, index_size_bytes, index_md5, recommended_mode
    )],
    stage_regions,
    by = "source_path",
    all = FALSE,
    sort = FALSE
  )
  ans[, `:=`(stage = stage_label, max_priority_rank = rank)]
  data.table::setcolorder(ans, c(
    "stage", "max_priority_rank", "modality", "level", "cell_type", "chrom",
    "source_path", "data_url", "index_url", "size_bytes", "md5",
    "index_size_bytes", "index_md5", "n_loci", "tabix_regions",
    "selection_reasons", "recommended_mode"
  ))
  ans
}))
data.table::setorder(acquisition_plan, max_priority_rank, modality, level, cell_type, chrom)

summary_rows <- data.table::rbindlist(lapply(0:2, function(rank) {
  stage_files <- file_plan[priority_rank <= rank]
  stage_rows <- mapped[priority_rank <= rank]
  data.table::data.table(
    stage = paste0("through_", priority_labels[rank + 1L]),
    max_priority_rank = rank,
    n_files = nrow(stage_files),
    n_loci = data.table::uniqueN(stage_rows$locus_id),
    n_features = data.table::uniqueN(stage_rows$feature_id),
    n_selected_rows = nrow(stage_rows),
    full_chromosome_download_bytes = sum(stage_files$size_bytes),
    full_chromosome_download_gib = sum(stage_files$size_bytes) / 1024^3
  )
}))

ensure_dir(output_dir)
rows_path <- file.path(output_dir, "target_rows.tsv.gz")
regions_path <- file.path(output_dir, "target_regions.tsv")
files_path <- file.path(output_dir, "target_files.tsv")
acquisition_path <- file.path(output_dir, "acquisition_plan.tsv")
summary_path <- file.path(output_dir, "target_summary.tsv")
data.table::fwrite(mapped, rows_path, sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(region_plan, regions_path, sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(file_plan, files_path, sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(acquisition_plan, acquisition_path, sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(summary_rows, summary_path, sep = "\t", quote = FALSE, na = "NA")

summary_md <- c(
  "# Targeted Multiome SuSiE acquisition plan",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Step 07 input files: ", length(coloc_files)),
  paste0("- Report threshold: PP.H4 >= ", format(report_pp)),
  paste0("- Strong threshold: PP.H4 >= ", format(strong_pp)),
  paste0("- QTL multi-signal threshold: at least ", qtl_cs_min, " credible sets"),
  paste0("- Region padding: ", padding_bp, " bp"),
  "",
  "| Stage | Files | Loci | Features | Selected rows | Full-chromosome GiB |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(summary_rows)), function(i) {
    paste0(
      "| ", summary_rows$stage[[i]], " | ", summary_rows$n_files[[i]], " | ",
      summary_rows$n_loci[[i]], " | ", summary_rows$n_features[[i]], " | ",
      summary_rows$n_selected_rows[[i]], " | ",
      sprintf("%.2f", summary_rows$full_chromosome_download_gib[[i]]), " |"
    )
  }, character(1)),
  "",
  "The GiB column is the worst-case size for full chromosome-file downloads. Remote Tabix extraction transfers and stores only planned locus slices."
)
writeLines(summary_md, file.path(output_dir, "TARGET_PLAN_SUMMARY.md"), useBytes = TRUE)

log_message("Wrote selected rows: ", rows_path, " (", nrow(mapped), " rows)")
log_message("Wrote target regions: ", regions_path, " (", nrow(region_plan), " rows)")
log_message("Wrote exact file plan: ", files_path, " (", nrow(file_plan), " files)")
log_message("Wrote stage-specific acquisition plan: ", acquisition_path)
log_message("Wrote cumulative summary: ", summary_path)
