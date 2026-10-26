#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "LAVA"))
cfg <- read_analysis_config(project_root)
args <- parse_cli()
phenos <- split_csv(args$traits %||% paste(
  cfg$primary_comparison$graves,
  cfg$primary_comparison$autoimmune_hypothyroidism,
  sep = ","
))

ref_prefix <- project_path(project_root, cfg$lava$reference_prefix)
ref_build <- cfg$lava$reference_build %||% ""
overlap_file <- project_path(project_root, cfg$lava$sample_overlap_file)
if (!nzchar(ref_prefix)) stop("Set lava.reference_prefix in config/analysis.yml", call. = FALSE)
if (!nzchar(ref_build)) stop("Set lava.reference_build for provenance in config/analysis.yml", call. = FALSE)
if (!nzchar(overlap_file) || !file.exists(overlap_file)) {
  stop("Configure an LDSC-derived LAVA sample-overlap matrix; these FinnGen traits are not independent.", call. = FALSE)
}

plink_ok <- all(file.exists(paste0(ref_prefix, c(".bed", ".bim", ".fam"))))
single_ld_ok <- all(file.exists(paste0(ref_prefix, c(".bcor", ".info"))))
multi_bcor <- Sys.glob(paste0(sub("_chr$", "", ref_prefix), "_chr*.bcor"))
multi_info <- Sys.glob(paste0(sub("_chr$", "", ref_prefix), "_chr*.info"))
multi_ld_ok <- length(multi_bcor) > 0L && setequal(sub("\\.bcor$", "", multi_bcor), sub("\\.info$", "", multi_info))
reference_ok <- plink_ok || single_ld_ok || multi_ld_ok
if (!reference_ok) stop("Could not find PLINK or LAVA reference files for prefix: ", ref_prefix, call. = FALSE)

info_file <- file.path(project_root, cfg$project$outcome_dir, "lava", "input", "input.info.tsv")
locus_file <- file.path(project_root, cfg$project$outcome_dir, "loci", "primary_shared_divergent_loci.lava")
if (!file.exists(info_file) || !file.exists(locus_file)) {
  stop("Run codes/02_define_loci.R and codes/03_prepare_ld_methods.R first.", call. = FALSE)
}

locus_preview <- data.table::fread(locus_file)
planned_locus_count <- nrow(locus_preview)
has_snp_column <- "SNPS" %in% names(locus_preview)
missing_snp_list <- if (has_snp_column) is.na(locus_preview$SNPS) | !nzchar(locus_preview$SNPS) else rep(TRUE, planned_locus_count)
preflight_failures <- data.table::data.table()
if (!identical(ref_build, cfg$project$genome_build)) {
  if (!has_snp_column) {
    stop("LAVA loci without SNPS require reference build ", cfg$project$genome_build, "; configured: ", ref_build, call. = FALSE)
  }
  if (any(missing_snp_list)) {
    preflight_failures <- locus_preview[missing_snp_list, .(
      locus_id = LOC,
      stage = "preflight",
      reason = paste0(
        "No eligible rsID SNP list; excluded because the ", ref_build,
        " LAVA reference cannot use ", cfg$project$genome_build, " locus coordinates."
      )
    )]
    locus_preview <- locus_preview[!missing_snp_list]
    log_message(
      "Excluding ", nrow(preflight_failures), " locus/loci without eligible rsID lists; ",
      nrow(locus_preview), " loci remain for SNP-list analysis."
    )
  }
  if (!nrow(locus_preview)) stop("No LAVA loci have eligible rsID lists.", call. = FALSE)
  log_message("LAVA reference build differs from GWAS coordinates; eligible loci will use SNP lists so coordinates are ignored.")
}

loci <- LAVA::read.loci(locus_file)
loci <- loci[loci$LOC %in% locus_preview$LOC, , drop = FALSE]
input <- LAVA::process.input(
  input.info.file = info_file,
  sample.overlap.file = overlap_file,
  ref.prefix = ref_prefix,
  phenos = phenos
)
univ_threshold <- as.numeric(cfg$lava$univariate_familywise_alpha) / planned_locus_count

out_dir <- file.path(project_root, cfg$project$outcome_dir, "lava", "results")
ensure_dir(out_dir)
write_tsv(data.table::data.table(
  run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  phenotypes = paste(phenos, collapse = ";"),
  planned_loci = planned_locus_count,
  analyzed_loci = nrow(loci),
  preflight_excluded_loci = nrow(preflight_failures),
  univariate_threshold = univ_threshold,
  reference_prefix = ref_prefix,
  reference_build = ref_build,
  overlap_file = overlap_file
), file.path(out_dir, "lava_run_provenance.tsv"))

univ <- list()
bivar <- list()
failures <- list()
if (nrow(preflight_failures)) failures[[length(failures) + 1L]] <- preflight_failures
for (i in seq_len(nrow(loci))) {
  log_message("LAVA locus ", i, "/", nrow(loci), ": ", loci$LOC[i])
  loc <- tryCatch(LAVA::process.locus(loci[i, ], input), error = identity)
  if (inherits(loc, "error") || is.null(loc)) {
    failures[[length(failures) + 1L]] <- data.table::data.table(
      locus_id = loci$LOC[i],
      stage = "process.locus",
      reason = if (inherits(loc, "error")) conditionMessage(loc) else "process.locus returned NULL"
    )
    next
  }
  result <- tryCatch(LAVA::run.univ.bivar(loc, univ.thresh = univ_threshold), error = identity)
  if (inherits(result, "error")) {
    failures[[length(failures) + 1L]] <- data.table::data.table(
      locus_id = loci$LOC[i], stage = "run.univ.bivar", reason = conditionMessage(result)
    )
    next
  }
  locus_info <- data.frame(
    locus_id = loc$id, chrom = loc$chr, start = loc$start, end = loc$stop,
    n_snps = loc$n.snps, n_pcs = loc$K
  )
  univ[[length(univ) + 1L]] <- data.table::as.data.table(cbind(locus_info, result$univ))
  if (!is.null(result$bivar)) bivar[[length(bivar) + 1L]] <- data.table::as.data.table(cbind(locus_info, result$bivar))
}

univ_result <- if (length(univ)) data.table::rbindlist(univ, fill = TRUE) else data.table::data.table()
bivar_result <- if (length(bivar)) data.table::rbindlist(bivar, fill = TRUE) else data.table::data.table()
failure_result <- if (length(failures)) data.table::rbindlist(failures, fill = TRUE) else data.table::data.table()
if (nrow(univ_result)) write_tsv(univ_result, file.path(out_dir, "primary.univ.lava.tsv"))
if (nrow(bivar_result)) write_tsv(bivar_result, file.path(out_dir, "primary.bivar.lava.tsv"))
if (nrow(failure_result)) write_tsv(failure_result, file.path(out_dir, "failed_loci.tsv"))

locus_phenotype_counts <- if (nrow(univ_result)) univ_result[, .N, by = locus_id] else data.table::data.table(N = integer())
result_qc <- data.table::data.table(
  metric = c(
    "planned_loci", "preflight_excluded_loci", "loci_with_univariate_results",
    "loci_with_both_phenotypes", "loci_with_one_phenotype", "univariate_rows",
    paste0("significant_univariate_", phenos[1]), paste0("significant_univariate_", phenos[2]),
    "bivariate_rows", "significant_bivariate_p_lt_0.05"
  ),
  value = c(
    planned_locus_count, nrow(preflight_failures), data.table::uniqueN(univ_result$locus_id),
    locus_phenotype_counts[N == length(phenos), .N], locus_phenotype_counts[N == 1L, .N], nrow(univ_result),
    univ_result[phen == phenos[1] & p < univ_threshold, .N],
    univ_result[phen == phenos[2] & p < univ_threshold, .N],
    nrow(bivar_result), bivar_result[p < 0.05, .N]
  )
)
write_tsv(result_qc, file.path(out_dir, "lava_result_qc.tsv"))

locus_status <- data.table::data.table(locus_id = data.table::fread(locus_file, select = "LOC")$LOC)
phenotype_status <- if (nrow(univ_result)) univ_result[, .(
  phenotypes_with_univariate_result = paste(sort(unique(phen)), collapse = ";"),
  n_phenotypes_with_univariate_result = data.table::uniqueN(phen)
), by = locus_id] else data.table::data.table(
  locus_id = character(), phenotypes_with_univariate_result = character(),
  n_phenotypes_with_univariate_result = integer()
)
locus_status <- merge(locus_status, phenotype_status, by = "locus_id", all.x = TRUE)
locus_status[, `:=`(
  status = data.table::fcase(
    locus_id %in% preflight_failures$locus_id, "preflight_excluded",
    n_phenotypes_with_univariate_result == length(phenos), "both_phenotypes_univariate",
    n_phenotypes_with_univariate_result == 1L, "one_phenotype_univariate",
    default = "no_univariate_result"
  ),
  bivariate_tested = locus_id %in% bivar_result$locus_id,
  bivariate_p = bivar_result$p[match(locus_id, bivar_result$locus_id)]
)]
write_tsv(locus_status, file.path(out_dir, "lava_locus_status.tsv"))
log_message("LAVA analysis finished.")
