#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "LAVA"))
cfg <- read_analysis_config(project_root)
graves <- cfg$primary_comparison$graves
aiht <- cfg$primary_comparison$autoimmune_hypothyroidism
phenos <- c(graves, aiht)

ref_prefix <- project_path(project_root, cfg$lava$reference_prefix)
ref_build <- cfg$lava$reference_build
gwas_build <- cfg$project$genome_build
overlap_file <- file.path(project_root, "outcome_data", "t05", "lava", "input", "sample.overlap.txt")
source_info_file <- file.path(project_root, "outcome_data", "lava", "input", "input.info.tsv")
t05_root <- file.path(project_root, "outcome_data", "t05", "lava")
input_dir <- file.path(t05_root, "input")
ensure_dir(input_dir)

if (!file.exists(overlap_file)) stop("Run codes/16_run_t05_ldsc.sh first", call. = FALSE)
if (!file.exists(source_info_file)) stop("Missing LAVA input metadata", call. = FALSE)
if (!all(file.exists(paste0(ref_prefix, c(".bed", ".bim", ".fam"))))) {
  stop("Missing PLINK LAVA reference for prefix: ", ref_prefix, call. = FALSE)
}

info <- data.table::fread(source_info_file)[phenotype %in% phenos]
if (nrow(info) != 2L) stop("Primary phenotype rows are missing from LAVA input.info", call. = FALSE)
info_file <- file.path(input_dir, "input.info.tsv")
write_tsv(info, info_file)

source_loci_path <- file.path(project_root, "outcome_data", "loci", "primary_shared_divergent_loci.lava")
source_metadata_path <- file.path(project_root, "outcome_data", "loci", "primary_shared_divergent_loci.tsv")
source_loci <- data.table::fread(source_loci_path)
source_metadata <- data.table::fread(source_metadata_path, select = c("locus_id", "is_hla"))
discovery <- merge(source_loci, source_metadata, by.x = "LOC", by.y = "locus_id", all.x = TRUE)
if (any(is.na(discovery$is_hla))) stop("Locus HLA annotation is incomplete", call. = FALSE)
hla_excluded <- discovery[is_hla == TRUE]
discovery <- discovery[is_hla == FALSE, .(LOC, CHR, START, STOP, SNPS)]
discovery_path <- file.path(input_dir, "discovery_non_hla.lava")
write_tsv(discovery, discovery_path)
write_tsv(hla_excluded, file.path(t05_root, "hla_excluded_from_ordinary_lava.tsv"))

cc_loci <- data.table::fread(file.path(project_root, "outcome_data", "ccgwas", "summary", "provisional_lead_loci.tsv"))
annotation <- data.table::fread(file.path(project_root, "outcome_data", "ccgwas", "input", "matched_variant_annotation.tsv.gz"))
lead_positions <- annotation[SNP %in% cc_loci$lead_SNP, .(lead_SNP = SNP, lead_BP = BP)]
cc_loci <- merge(cc_loci, lead_positions, by = "lead_SNP", all.x = TRUE)
if (any(is.na(cc_loci$lead_BP))) stop("Missing lead positions for targeted CC-GWAS loci", call. = FALSE)
window <- as.integer(cfg$gwas_qc$locus_half_window_bp)
targeted_rows <- lapply(seq_len(nrow(cc_loci)), function(i) {
  chr <- as.integer(cc_loci$CHR[[i]])
  start <- max(1L, as.integer(cc_loci$lead_BP[[i]]) - window)
  stop <- as.integer(cc_loci$lead_BP[[i]]) + window
  snps <- annotation[CHR == chr & BP >= start & BP <= stop, unique(SNP)]
  data.table::data.table(
    LOC = paste0("T", sprintf("%02d", i), "_", cc_loci$locus_id[[i]]),
    CHR = chr,
    START = start,
    STOP = stop,
    SNPS = paste(snps, collapse = ";"),
    ccgwas_locus_id = cc_loci$locus_id[[i]],
    lead_SNP = cc_loci$lead_SNP[[i]],
    lead_nearest_gene = cc_loci$lead_nearest_genes_graves[[i]],
    n_candidate_snps = length(snps),
    overlaps_hla = chr == as.integer(cfg$hla$chromosome) &
      start <= as.integer(cfg$hla$end) & stop >= as.integer(cfg$hla$start),
    hla_adjacent_within_2mb = chr == as.integer(cfg$hla$chromosome) &
      start > as.integer(cfg$hla$end) & start <= as.integer(cfg$hla$end) + 2000000L
  )
})
targeted_map <- data.table::rbindlist(targeted_rows)
targeted_hla_excluded <- targeted_map[overlaps_hla == TRUE]
if (nrow(targeted_hla_excluded)) {
  write_tsv(targeted_hla_excluded, file.path(t05_root, "targeted_hla_excluded.tsv"))
}
targeted <- targeted_map[overlaps_hla == FALSE]
if (any(!nzchar(targeted$SNPS))) stop("A targeted locus has no eligible SNPs", call. = FALSE)
write_tsv(targeted_map, file.path(input_dir, "targeted_ccgwas_locus_map.tsv"))
targeted_path <- file.path(input_dir, "targeted_ccgwas_non_hla.lava")
write_tsv(targeted[, .(LOC, CHR, START, STOP, SNPS)], targeted_path)

input <- LAVA::process.input(
  input.info.file = info_file,
  sample.overlap.file = overlap_file,
  ref.prefix = ref_prefix,
  phenos = phenos
)

run_panel <- function(panel_name, locus_file) {
  out_dir <- file.path(t05_root, panel_name)
  ensure_dir(out_dir)
  planned <- data.table::fread(locus_file)
  planned_count <- nrow(planned)
  if (!planned_count) stop("No planned loci for panel ", panel_name, call. = FALSE)
  has_snps <- "SNPS" %in% names(planned)
  missing_snps <- if (has_snps) is.na(planned$SNPS) | !nzchar(planned$SNPS) else rep(TRUE, planned_count)
  preflight <- data.table::data.table()
  eligible_preview <- planned
  if (!identical(ref_build, gwas_build)) {
    if (!has_snps) stop("Build mismatch requires SNP lists", call. = FALSE)
    if (any(missing_snps)) {
      preflight <- planned[missing_snps, .(
        locus_id = LOC,
        stage = "preflight",
        reason = paste0(
          "No eligible rsID SNP list; ", gwas_build,
          " coordinates cannot be applied to the ", ref_build, " reference"
        )
      )]
      eligible_preview <- planned[!missing_snps]
    }
  }
  loci <- LAVA::read.loci(locus_file)
  loci <- loci[loci$LOC %in% eligible_preview$LOC, , drop = FALSE]
  threshold <- as.numeric(cfg$lava$univariate_familywise_alpha) / planned_count

  provenance <- data.table::data.table(
    run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    panel = panel_name,
    phenotypes = paste(phenos, collapse = ";"),
    planned_non_hla_loci = planned_count,
    loci_after_preflight = nrow(loci),
    preflight_excluded_loci = nrow(preflight),
    univariate_bonferroni_threshold = threshold,
    reference_prefix = ref_prefix,
    reference_build = ref_build,
    gwas_build = gwas_build,
    sample_overlap_file = overlap_file,
    lava_version = as.character(utils::packageVersion("LAVA")),
    hla_included = FALSE
  )
  write_tsv(provenance, file.path(out_dir, "run_provenance.tsv"))

  univ <- list()
  bivar <- list()
  failures <- list()
  if (nrow(preflight)) failures[[length(failures) + 1L]] <- preflight
  for (i in seq_len(nrow(loci))) {
    if (i == 1L || i %% 20L == 0L || i == nrow(loci)) {
      log_message("T05 LAVA ", panel_name, " ", i, "/", nrow(loci), ": ", loci$LOC[[i]])
    }
    loc <- tryCatch(LAVA::process.locus(loci[i, ], input), error = identity)
    if (inherits(loc, "error") || is.null(loc)) {
      failures[[length(failures) + 1L]] <- data.table::data.table(
        locus_id = loci$LOC[[i]],
        stage = "process.locus",
        reason = if (inherits(loc, "error")) conditionMessage(loc) else "process.locus returned NULL"
      )
      next
    }
    result <- tryCatch(
      LAVA::run.univ.bivar(loc, univ.thresh = threshold),
      error = identity
    )
    if (inherits(result, "error") || is.null(result)) {
      failures[[length(failures) + 1L]] <- data.table::data.table(
        locus_id = loci$LOC[[i]],
        stage = "run.univ.bivar",
        reason = if (inherits(result, "error")) conditionMessage(result) else "run.univ.bivar returned NULL"
      )
      next
    }
    locus_info <- data.frame(
      locus_id = loc$id,
      chrom = loc$chr,
      start = loc$start,
      end = loc$stop,
      n_snps = loc$n.snps,
      n_pcs = loc$K
    )
    univ[[length(univ) + 1L]] <- data.table::as.data.table(cbind(locus_info, result$univ))
    if (!is.null(result$bivar)) {
      bivar[[length(bivar) + 1L]] <- data.table::as.data.table(cbind(locus_info, result$bivar))
    }
  }
  univ_result <- if (length(univ)) data.table::rbindlist(univ, fill = TRUE) else data.table::data.table()
  bivar_result <- if (length(bivar)) data.table::rbindlist(bivar, fill = TRUE) else data.table::data.table()
  failure_result <- if (length(failures)) data.table::rbindlist(failures, fill = TRUE) else data.table::data.table()
  if (nrow(univ_result)) write_tsv(univ_result, file.path(out_dir, "univariate.tsv"))
  if (nrow(bivar_result)) write_tsv(bivar_result, file.path(out_dir, "bivariate.tsv"))
  if (nrow(failure_result)) write_tsv(failure_result, file.path(out_dir, "failures.tsv"))
  log_message("Completed T05 LAVA panel ", panel_name, ": ", nrow(bivar_result), " bivariate test(s)")
}

run_panel("discovery", discovery_path)
run_panel("targeted", targeted_path)
