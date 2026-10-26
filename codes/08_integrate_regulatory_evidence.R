#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
cell_types <- read_cell_types(project_root, cfg)
graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism

coloc_dir <- file.path(project_root, cfg$project$outcome_dir, "coloc_r13")
coloc_files <- list.files(coloc_dir, pattern = "\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
if (!length(coloc_files)) stop("No R13 coloc output found. Run codes/07_coloc_multiome_r13.R first.", call. = FALSE)
coloc <- data.table::rbindlist(lapply(coloc_files, data.table::fread, showProgress = FALSE), fill = TRUE)
if (!nrow(coloc)) stop("R13 coloc output files are empty.", call. = FALSE)
coloc <- coloc[PP.H4.abf >= as.numeric(cfg$coloc$pp_h4_strong)]
if (!nrow(coloc)) stop("No colocalizations pass the configured strong PP.H4 threshold.", call. = FALSE)

keys <- c("locus_id", "chrom", "start", "end", "is_hla", "modality", "level", "cell_type", "feature_id")
value_cols <- c(
  "PP.H4.abf", "effect_direction", "top_h4_variant", "top_snp_pp_h4",
  "gwas_susie_pip", "gwas_susie_cs", "gwas_susie_low_purity",
  "qtl_n_credible_sets", "qtl_max_cs_log10bf", "qtl_has_pure_cs", "n_shared_snps"
)
g <- coloc[trait_id == graves_id, c(keys, value_cols), with = FALSE]
h <- coloc[trait_id == aiht_id, c(keys, value_cols), with = FALSE]
paired <- merge(g, h, by = keys, all = TRUE, suffixes = c("_graves", "_aiht"))
paired[, `:=`(
  graves_coloc = !is.na(PP.H4.abf_graves),
  aiht_coloc = !is.na(PP.H4.abf_aiht)
)]
paired[, regulatory_class := data.table::fcase(
  graves_coloc & aiht_coloc & effect_direction_graves == effect_direction_aiht,
  "shared_coloc_same_regulatory_direction",
  graves_coloc & aiht_coloc & effect_direction_graves != effect_direction_aiht,
  "shared_coloc_opposite_regulatory_direction",
  graves_coloc & !aiht_coloc,
  "graves_coloc_only_not_trait_specificity_proof",
  !graves_coloc & aiht_coloc,
  "aiht_coloc_only_not_trait_specificity_proof",
  default = "unclassified"
)]
paired[, `:=`(
  qtl_n_credible_sets = data.table::fcoalesce(qtl_n_credible_sets_graves, qtl_n_credible_sets_aiht),
  qtl_max_cs_log10bf = data.table::fcoalesce(qtl_max_cs_log10bf_graves, qtl_max_cs_log10bf_aiht),
  qtl_has_pure_cs = data.table::fcoalesce(qtl_has_pure_cs_graves, qtl_has_pure_cs_aiht),
  lineage = cell_types$lineage[match(paste(level, cell_type), paste(cell_types$level, cell_types$cell_type))]
)]

link_rows <- list()
ca_pairs <- paired[modality == "caQTL"]
if (nrow(ca_pairs)) {
  for (cell_id in unique(paste(ca_pairs$level, ca_pairs$cell_type, sep = "."))) {
    parts <- strsplit(cell_id, "\\.")[[1]]
    level_value <- parts[1]
    cell_value <- paste(parts[-1], collapse = ".")
    peaks <- unique(ca_pairs[level == level_value & cell_type == cell_value, feature_id])
    path <- file.path(
      project_root, "raw_data", "finngen_multiome", "peak_gene_links",
      paste0("finngen_multiome_v1.peak_gene_links.", level_value, ".", cell_value, ".tsv.gz")
    )
    if (!file.exists(path)) next
    links <- data.table::fread(
      path,
      select = c("peak_id", "gene_id", "cell_type", "q_joint", "p_adj_zero", "p_adj_count", "mode", "sig"),
      showProgress = FALSE
    )
    links <- links[peak_id %in% peaks & as.logical(sig)]
    if (!nrow(links)) next
    links[, `:=`(level = level_value, short_cell_type = cell_value)]
    data.table::setnames(links, c("peak_id", "gene_id"), c("feature_id", "linked_gene_id"))
    link_rows[[length(link_rows) + 1L]] <- links
  }
}
links <- if (length(link_rows)) data.table::rbindlist(link_rows, fill = TRUE) else data.table::data.table(
  level = character(), short_cell_type = character(), feature_id = character(), linked_gene_id = character(),
  q_joint = numeric(), p_adj_zero = numeric(), p_adj_count = numeric(), mode = character(), sig = logical()
)

eqtl_evidence <- paired[modality == "eQTL"]
eqtl_evidence[, `:=`(
  candidate_gene_id = feature_id,
  peak_gene_link_mode = NA_character_,
  peak_gene_link_q = NA_real_
)]
ca_evidence <- merge(
  paired[modality == "caQTL"],
  links[, .(level, cell_type = short_cell_type, feature_id, linked_gene_id, peak_gene_link_mode = mode, peak_gene_link_q = q_joint)],
  by = c("level", "cell_type", "feature_id"),
  all.x = TRUE
)
ca_evidence[, candidate_gene_id := linked_gene_id]
evidence <- data.table::rbindlist(list(eqtl_evidence, ca_evidence), fill = TRUE, use.names = TRUE)

locus_meta_path <- file.path(project_root, cfg$project$outcome_dir, "loci", "primary_shared_divergent_loci.tsv")
if (file.exists(locus_meta_path)) {
  locus_meta <- data.table::fread(
    locus_meta_path,
    select = c("locus_id", "lead_variant", "provisional_class", "requires_local_rg_confirmation")
  )
  evidence <- merge(evidence, locus_meta, by = "locus_id", all.x = TRUE)
}

gene_map_path <- project_path(project_root, cfg$optional_annotations$gene_map_tsv)
if (nzchar(gene_map_path) && file.exists(gene_map_path)) {
  gene_map <- data.table::fread(gene_map_path)
  if (!all(c("gene_id", "gene_symbol") %in% names(gene_map))) {
    stop("optional gene map must contain gene_id and gene_symbol", call. = FALSE)
  }
  evidence <- merge(evidence, gene_map[, .(candidate_gene_id = gene_id, gene_symbol)], by = "candidate_gene_id", all.x = TRUE)
} else {
  evidence[, gene_symbol := NA_character_]
}

evidence[, evidence_tier := data.table::fcase(
  is_hla, "HLA_separate_complex_LD",
  graves_coloc & aiht_coloc & qtl_has_pure_cs %in% TRUE &
    pmax(gwas_susie_pip_graves, gwas_susie_pip_aiht, na.rm = TRUE) >= 0.10,
  "tier_1_shared_coloc_plus_finemap",
  graves_coloc & aiht_coloc,
  "tier_2_shared_coloc",
  default = "tier_3_single_trait_coloc"
)]

lava_path <- file.path(project_root, cfg$project$outcome_dir, "lava", "results", "primary.bivar.lava.tsv")
if (file.exists(lava_path)) {
  lava <- data.table::fread(lava_path)
  if ("locus" %in% names(lava) && !"locus_id" %in% names(lava)) data.table::setnames(lava, "locus", "locus_id")
  if (all(c("phen1", "phen2") %in% names(lava))) {
    lava <- lava[(phen1 == graves_id & phen2 == aiht_id) | (phen1 == aiht_id & phen2 == graves_id)]
  }
  if ("locus_id" %in% names(lava)) evidence <- merge(evidence, lava, by = "locus_id", all.x = TRUE, suffixes = c("", "_lava"))
}

selected_variants <- unique(na.omit(c(evidence$top_h4_variant_graves, evidence$top_h4_variant_aiht)))
cascade_path <- file.path(project_root, "raw_data", "finngen_multiome", "cascade", "finngen_multiome_v1.cascade.variant.tsv.gz")
if (length(selected_variants) && file.exists(cascade_path)) {
  cascade <- data.table::fread(cascade_path, showProgress = TRUE)
  cascade <- cascade[variant_id %in% selected_variants]
  write_tsv(cascade, file.path(project_root, cfg$project$outcome_dir, "integrated", "selected_variant_cascade_annotations.tsv.gz"))
}

out_dir <- file.path(project_root, cfg$project$outcome_dir, "integrated")
write_tsv(evidence, file.path(out_dir, "shared_divergent_regulatory_evidence.tsv.gz"))

gene_summary <- evidence[!is.na(candidate_gene_id), .(
  gene_symbol = gene_symbol[which.max(!is.na(gene_symbol))],
  strongest_evidence_tier = sort(unique(evidence_tier))[1],
  regulatory_classes = paste(sort(unique(regulatory_class)), collapse = ";"),
  modalities = paste(sort(unique(modality)), collapse = ";"),
  lineages = paste(sort(unique(lineage)), collapse = ";"),
  cell_types = paste(sort(unique(paste(level, cell_type, sep = "."))), collapse = ";"),
  loci = paste(sort(unique(locus_id)), collapse = ";"),
  max_pp_h4_graves = if (all(is.na(PP.H4.abf_graves))) NA_real_ else max(PP.H4.abf_graves, na.rm = TRUE),
  max_pp_h4_aiht = if (all(is.na(PP.H4.abf_aiht))) NA_real_ else max(PP.H4.abf_aiht, na.rm = TRUE),
  includes_hla = any(is_hla)
), by = candidate_gene_id]
write_tsv(gene_summary, file.path(out_dir, "candidate_gene_celltype_summary.tsv"))
log_message(
  "Integrated ", nrow(evidence), " regulatory evidence rows across ",
  data.table::uniqueN(stats::na.omit(evidence$candidate_gene_id)), " candidate genes."
)
