#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "codes/24_review_t08_targeted_susie_plan.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

core_dir <- file.path(project_root, "outcome_data", "targeted_susie_core_lineages")
out_dir <- file.path(project_root, "outcome_data", "t08")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty input: ", path)
  read.delim(path, check.names = FALSE, na.strings = c("", "NA"))
}

write_tsv <- function(x, name) {
  path <- file.path(out_dir, name)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote: ", path, " (", nrow(x), " rows)")
  invisible(path)
}

collapse_unique <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")

regions <- read_tsv(file.path(core_dir, "target_regions.tsv"))
files <- read_tsv(file.path(core_dir, "target_files.tsv"))
rows <- read.delim(gzfile(file.path(core_dir, "target_rows.tsv.gz")), check.names = FALSE,
                   na.strings = c("", "NA"))
t07 <- read_tsv(file.path(project_root, "outcome_data", "t07", "genetic_signal_five_classification.tsv"))
loci <- read_tsv(file.path(project_root, "outcome_data", "loci", "primary_shared_divergent_loci.tsv"))
manifest <- read_tsv(file.path(project_root, "raw_data", "manifests", "finngen_multiome_v1.manifest.tsv"))

required_region_cols <- c("source_path", "modality", "level", "cell_type", "locus_id", "priority_rank",
                          "priority", "chrom", "start", "end", "tabix_region", "selection_reasons")
required_file_cols <- c("source_path", "data_url", "index_url", "size_bytes", "md5",
                        "index_size_bytes", "index_md5", "recommended_mode")
required_row_cols <- c("source_path", "locus_id", "top_h4_variant", "PP.H4.abf", "genome_build")
stopifnot(all(required_region_cols %in% names(regions)),
          all(required_file_cols %in% names(files)),
          all(required_row_cols %in% names(rows)))

target_classes <- c("shared-opposite", "Graves-predominant", "AIHT-predominant")
critical <- t07[t07$five_class_primary %in% target_classes, , drop = FALSE]
critical$is_hla <- critical$chromosome == 6 & critical$position >= 24239950 & critical$position <= 37114944
critical <- critical[!critical$is_hla, , drop = FALSE]

map_locus <- function(chr, pos) {
  hit <- loci[loci$chrom == chr & loci$start <= pos & loci$end >= pos, , drop = FALSE]
  if (nrow(hit) != 1L) return(NA_character_)
  hit$locus_id[[1]]
}
critical$locus_id <- mapply(map_locus, critical$chromosome, critical$position)
if (nrow(critical) != 4L || anyNA(critical$locus_id)) {
  stop("Expected four non-HLA T07 divergent signals with unique physical-locus mappings; observed ",
       nrow(critical), " signals and ", sum(is.na(critical$locus_id)), " unmapped.")
}

legacy <- regions[regions$priority_rank == 0, , drop = FALSE]
promoted <- regions[regions$locus_id %in% critical$locus_id, , drop = FALSE]
freeze_keys <- unique(c(paste(legacy$source_path, legacy$locus_id, sep = "||"),
                        paste(promoted$source_path, promoted$locus_id, sep = "||")))
frozen <- regions[paste(regions$source_path, regions$locus_id, sep = "||") %in% freeze_keys, , drop = FALSE]

critical_class <- setNames(critical$five_class_primary, critical$locus_id)
frozen$original_priority <- frozen$priority
frozen$original_priority_rank <- frozen$priority_rank
frozen$t07_class <- unname(critical_class[frozen$locus_id])
frozen$freeze_source <- ifelse(frozen$locus_id %in% critical$locus_id,
                               paste0("T07_", frozen$t07_class, "_promotion"),
                               ifelse(frozen$locus_id == "L067_chr6_24239950_37114944",
                                      "legacy_P0_HLA_separate", "legacy_P0_screen"))
frozen$frozen_priority <- "P0_critical"
frozen$frozen_priority_rank <- 0L
frozen$is_hla <- frozen$locus_id == "L067_chr6_24239950_37114944"
frozen <- frozen[order(frozen$is_hla, frozen$chrom, frozen$source_path, frozen$start), ]

file_meta <- files[!duplicated(files$source_path), required_file_cols, drop = FALSE]

make_plan <- function(region_subset) {
  if (!nrow(region_subset)) {
    return(data.frame(stage = character(), max_priority_rank = integer(), modality = character(),
                      level = character(), cell_type = character(), chrom = character(),
                      source_path = character(), data_url = character(), index_url = character(),
                      size_bytes = numeric(), md5 = character(), index_size_bytes = numeric(),
                      index_md5 = character(), n_loci = integer(), tabix_regions = character(),
                      selection_reasons = character(), recommended_mode = character()))
  }
  groups <- split(region_subset, region_subset$source_path)
  z <- do.call(rbind, lapply(groups, function(g) {
    data.frame(
      stage = "through_P0_critical",
      max_priority_rank = 0L,
      modality = g$modality[[1]],
      level = g$level[[1]],
      cell_type = g$cell_type[[1]],
      chrom = g$chrom[[1]],
      source_path = g$source_path[[1]],
      n_loci = length(unique(g$locus_id)),
      tabix_regions = collapse_unique(g$tabix_region),
      selection_reasons = collapse_unique(c(g$selection_reasons, g$freeze_source)),
      stringsAsFactors = FALSE
    )
  }))
  z <- merge(z, file_meta, by = "source_path", all.x = TRUE, sort = FALSE)
  expected <- c("stage", "max_priority_rank", "modality", "level", "cell_type", "chrom",
                "source_path", "data_url", "index_url", "size_bytes", "md5", "index_size_bytes",
                "index_md5", "n_loci", "tabix_regions", "selection_reasons", "recommended_mode")
  z <- z[, expected, drop = FALSE]
  z[order(z$modality, z$level, z$cell_type, z$chrom, z$source_path), , drop = FALSE]
}

plan_combined <- make_plan(frozen)
plan_non_hla <- make_plan(frozen[!frozen$is_hla, , drop = FALSE])
plan_hla <- make_plan(frozen[frozen$is_hla, , drop = FALSE])

coverage <- critical[, c("signal_unit_id", "five_class_primary", "chromosome", "position",
                         "representative_variant", "representative_nearest_genes", "locus_id")]
coverage$legacy_p0_n_source_files <- vapply(coverage$locus_id, function(x) {
  length(unique(legacy$source_path[legacy$locus_id == x]))
}, integer(1))
coverage$frozen_p0_n_source_files <- vapply(coverage$locus_id, function(x) {
  length(unique(frozen$source_path[frozen$locus_id == x]))
}, integer(1))
coverage$legacy_p0_covered <- coverage$legacy_p0_n_source_files > 0
coverage$frozen_p0_covered <- coverage$frozen_p0_n_source_files > 0

audit <- data.frame(
  audit_item = c("legacy_core_P0_files", "legacy_core_P0_loci", "legacy_core_P0_full_file_GiB",
                 "T07_divergent_loci_total", "T07_divergent_loci_covered_by_legacy_P0",
                 "T07_divergent_loci_covered_by_frozen_P0", "legacy_P0_is_current_for_T07"),
  value = c(length(unique(legacy$source_path)), length(unique(legacy$locus_id)),
            sprintf("%.5f", sum(unique(files[files$source_path %in% legacy$source_path,
                                                 c("source_path", "size_bytes")])$size_bytes) / 2^30),
            nrow(coverage), sum(coverage$legacy_p0_covered), sum(coverage$frozen_p0_covered), "FALSE"),
  interpretation = c("regenerated core-lineage plan", "regenerated core-lineage plan",
                     "hypothetical whole-file footprint; remote slices are required",
                     "shared-opposite plus disease-predominant non-HLA signals from T07",
                     "old plan misses T07-priority loci", "all current T07-priority loci included",
                     "old P0 predates T07 and must not be executed unchanged")
)

manifest_data <- manifest[manifest$path %in% plan_combined$source_path, c("path", "size_bytes", "md5")]
manifest_index <- manifest[manifest$path %in% paste0(plan_combined$source_path, ".tbi"),
                           c("path", "size_bytes", "md5")]
data_match <- match(plan_combined$source_path, manifest_data$path)
index_match <- match(paste0(plan_combined$source_path, ".tbi"), manifest_index$path)

grch38_lengths <- c(chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
                    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
                    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
                    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
                    chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
                    chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415)

region_syntax_ok <- grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+-[0-9]+$", frozen$tabix_region)
coordinate_ok <- frozen$start >= 1 & frozen$end >= frozen$start &
  frozen$end <= unname(grch38_lengths[frozen$chrom])
official_base <- "https://storage.googleapis.com/finngen-public-data-multiome/"
expected_cells <- c("CD4_T", "B_naive", "Treg", "CD14_Mono")

static_checks <- data.frame(
  check_id = c("S01_exact_plan_header", "S02_unique_combined_source_paths", "S03_official_https_urls",
               "S04_manifest_data_rows", "S05_manifest_index_rows", "S06_manifest_data_metadata",
               "S07_manifest_index_metadata", "S08_chr_prefix", "S09_region_syntax",
               "S10_grch38_coordinate_bounds", "S11_input_build_grch38", "S12_four_core_lineages",
               "S13_both_modalities", "S14_all_T07_divergent_loci", "S15_HLA_excluded_non_HLA_plan",
               "S16_HLA_only_HLA_plan", "S17_remote_extract_mode"),
  passed = c(
    identical(names(plan_combined), c("stage", "max_priority_rank", "modality", "level", "cell_type", "chrom",
                                      "source_path", "data_url", "index_url", "size_bytes", "md5",
                                      "index_size_bytes", "index_md5", "n_loci", "tabix_regions",
                                      "selection_reasons", "recommended_mode")),
    !anyDuplicated(plan_combined$source_path),
    all(startsWith(plan_combined$data_url, official_base)) && all(startsWith(plan_combined$index_url, official_base)),
    !anyNA(data_match), !anyNA(index_match),
    !anyNA(data_match) && all(plan_combined$size_bytes == manifest_data$size_bytes[data_match]) &&
      all(plan_combined$md5 == manifest_data$md5[data_match]),
    !anyNA(index_match) && all(plan_combined$index_size_bytes == manifest_index$size_bytes[index_match]) &&
      all(plan_combined$index_md5 == manifest_index$md5[index_match]),
    all(grepl("^chr", frozen$chrom)), all(region_syntax_ok), all(coordinate_ok),
    all(rows$genome_build == "GRCh38"), setequal(unique(plan_combined$cell_type), expected_cells),
    setequal(unique(plan_combined$modality), c("caQTL", "eQTL")), all(coverage$frozen_p0_covered),
    !any(grepl("chr6:24239950-37114944", plan_non_hla$tabix_regions, fixed = TRUE)),
    nrow(plan_hla) > 0 && all(grepl("chr6:24239950-37114944", plan_hla$tabix_regions, fixed = TRUE)),
    all(plan_combined$recommended_mode == "remote_extract")
  ),
  detail = c("17-column header required by acquisition script", "one extraction command per source file",
             "official FinnGen Multiome GCS endpoint", "every data file exists in manifest",
             "every Tabix index exists in manifest", "data byte counts and MD5 match manifest",
             "index byte counts and MD5 match manifest", "Tabix chromosome naming is chr-prefixed",
             "all target intervals have valid Tabix syntax", "intervals fit GRCh38 chromosome bounds",
             "all selected Step07 regulatory rows declare GRCh38", paste(expected_cells, collapse = ";"),
             "caQTL and eQTL are both represented", "all four T07 divergent loci are included",
             "ordinary plan contains no HLA interval", "separate plan contains HLA interval only",
             "no whole-file acquisition is recommended"),
  stringsAsFactors = FALSE
)

candidate_rows <- rows[paste(rows$source_path, rows$locus_id, sep = "||") %in%
                         paste(frozen$source_path, frozen$locus_id, sep = "||") &
                         !is.na(rows$top_h4_variant), , drop = FALSE]
candidate_rows <- candidate_rows[order(candidate_rows$source_path, candidate_rows$locus_id,
                                       -candidate_rows$PP.H4.abf), , drop = FALSE]
candidate_rows <- candidate_rows[!duplicated(paste(candidate_rows$source_path, candidate_rows$locus_id,
                                                   sep = "||")), , drop = FALSE]
variant_pos <- suppressWarnings(as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1",
                                              candidate_rows$top_h4_variant)))
preflight <- candidate_rows[, c("source_path", "modality", "level", "cell_type", "locus_id",
                                "data_url", "index_url", "top_h4_variant", "PP.H4.abf")]
frozen_chrom <- setNames(frozen$chrom, paste(frozen$source_path, frozen$locus_id, sep = "||"))
preflight$chrom <- unname(frozen_chrom[paste(preflight$source_path, preflight$locus_id, sep = "||")])
preflight$test_position <- variant_pos
preflight$test_region <- paste0(preflight$chrom, ":", variant_pos, "-", variant_pos)
preflight$is_hla <- preflight$locus_id == "L067_chr6_24239950_37114944"
preflight <- preflight[, c("source_path", "modality", "level", "cell_type", "locus_id", "is_hla",
                           "data_url", "index_url", "chrom", "top_h4_variant", "PP.H4.abf",
                           "test_position", "test_region")]
if (nrow(preflight) != nrow(frozen) || anyNA(preflight$test_position)) {
  stop("Could not assign one exact known-variant preflight target to every frozen source-locus row.")
}

summary_table <- do.call(rbind, lapply(list(combined = plan_combined, non_hla = plan_non_hla,
                                            hla = plan_hla), function(p) {
  data.frame(n_source_files = nrow(p), n_loci = length(unique(unlist(strsplit(p$tabix_regions, ";", fixed = TRUE)))),
             full_file_bytes = sum(p$size_bytes), full_file_gib = sum(p$size_bytes) / 2^30,
             n_region_requests = sum(p$n_loci), stringsAsFactors = FALSE)
}))
summary_table$plan <- rownames(summary_table)
rownames(summary_table) <- NULL
summary_table <- summary_table[, c("plan", "n_source_files", "n_loci", "n_region_requests",
                                   "full_file_bytes", "full_file_gib")]

decision <- data.frame(
  gate = "T08_static_plan_review",
  status = if (all(static_checks$passed)) "PASS" else "FAIL",
  decision = if (all(static_checks$passed)) "PROCEED_TO_REMOTE_TABIX_PREFLIGHT" else "STOP_AND_REPAIR_PLAN",
  rationale = if (all(static_checks$passed))
    "Frozen core-lineage P0 plan is manifest-exact, GRCh38-consistent, T07-aware, and HLA-separated."
  else paste(static_checks$check_id[!static_checks$passed], collapse = ";"),
  stringsAsFactors = FALSE
)

write_tsv(audit, "legacy_core_p0_audit.tsv")
write_tsv(coverage, "t07_critical_locus_coverage.tsv")
write_tsv(frozen, "p0_core_frozen_regions.tsv")
write_tsv(plan_combined, "p0_core_combined_acquisition_plan.tsv")
write_tsv(plan_non_hla, "p0_core_non_hla_acquisition_plan.tsv")
write_tsv(plan_hla, "p0_core_hla_acquisition_plan.tsv")
write_tsv(preflight, "remote_preflight_targets.tsv")
for (batch_id in 1:4) {
  batch_rows <- preflight[((seq_len(nrow(preflight)) - 1L) %% 4L) + 1L == batch_id, , drop = FALSE]
  write_tsv(batch_rows, sprintf("remote_preflight_targets_batch%d.tsv", batch_id))
}
write_tsv(summary_table, "p0_core_summary.tsv")
write_tsv(static_checks, "t08_static_checks.tsv")
write_tsv(decision, "t08_decision.tsv")

if (!all(static_checks$passed)) {
  stop("T08 static checks failed: ", paste(static_checks$check_id[!static_checks$passed], collapse = ", "))
}
message("T08 static review PASS. Remote Tabix preflight is required before the T09 gate can open.")
