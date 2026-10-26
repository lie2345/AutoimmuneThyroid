#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)
args <- parse_cli()

requested <- split_csv(args$traits %||% "")
if (length(requested)) {
  unknown <- setdiff(requested, traits$trait_id)
  if (length(unknown)) stop("Unknown trait IDs: ", paste(unknown, collapse = ", "), call. = FALSE)
  traits <- traits[trait_id %in% requested]
}

out_dir <- file.path(project_root, cfg$project$outcome_dir, "harmonized_gwas")
ensure_dir(out_dir)
qc_rows <- list()

for (i in seq_len(nrow(traits))) {
  meta <- traits[i]
  input <- gwas_path(project_root, meta$trait_id)
  output <- file.path(out_dir, paste0(meta$trait_id, ".harmonized.tsv.gz"))
  if (file.exists(output) && !as_flag(args$force)) {
    log_message("Skipping existing output: ", output)
    next
  }
  if (!file.exists(input)) stop("Missing GWAS file: ", input, call. = FALSE)
  log_message("Reading ", meta$trait_id)
  x <- data.table::fread(
    input,
    select = c("#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "beta", "sebeta", "af_alt"),
    showProgress = TRUE
  )
  data.table::setnames(x, "#chrom", "chrom")
  n_raw <- nrow(x)
  x[, chrom := sub("^chr", "", as.character(chrom), ignore.case = TRUE)]
  x[, `:=`(
    pos = as.integer(pos), ref = toupper(ref), alt = toupper(alt),
    rsid = first_rsid(rsids), p = as.numeric(pval), se = as.numeric(sebeta),
    eaf = as.numeric(af_alt), beta = as.numeric(beta)
  )]
  x <- x[
    !is.na(pos) & pos > 0L & nzchar(ref) & nzchar(alt) &
      is.finite(p) & p > 0 & p <= 1 & is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(eaf) & eaf > 0 & eaf < 1
  ]
  if (isTRUE(cfg$gwas_qc$autosomes_only)) x <- x[chrom %in% as.character(1:22)]
  x[, maf := pmin(eaf, 1 - eaf)]
  x <- x[maf >= as.numeric(cfg$gwas_qc$maf_min_general)]
  x[, `:=`(
    variant_id = canonical_variant_id(chrom, pos, ref, alt),
    effect_allele = alt,
    other_allele = ref,
    z = beta / se,
    n_cases = as.integer(meta$n_cases),
    n_controls = as.integer(meta$n_controls),
    n_total = as.integer(meta$n_cases + meta$n_controls),
    is_snp = nchar(ref) == 1L & nchar(alt) == 1L,
    palindromic = is_palindromic(ref, alt),
    is_hla = is_hla_position(chrom, pos, cfg)
  )]
  data.table::setorder(x, variant_id, p)
  x <- x[!duplicated(variant_id)]
  keep <- c(
    "variant_id", "chrom", "pos", "rsid", "ref", "alt", "effect_allele", "other_allele",
    "beta", "se", "z", "p", "eaf", "maf", "n_cases", "n_controls", "n_total",
    "nearest_genes", "is_snp", "palindromic", "is_hla"
  )
  data.table::setcolorder(x, keep)
  data.table::fwrite(x, output, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
  qc_rows[[length(qc_rows) + 1L]] <- data.table::data.table(
    trait_id = meta$trait_id,
    raw_rows = n_raw,
    retained_rows = nrow(x),
    retained_fraction = nrow(x) / n_raw,
    genome_wide_significant = sum(x$p < as.numeric(cfg$gwas_qc$genome_wide_p)),
    hla_rows = sum(x$is_hla),
    output = sub(paste0("^", project_root, "/"), "", output)
  )
  rm(x)
  invisible(gc())
}

if (length(qc_rows)) {
  write_tsv(data.table::rbindlist(qc_rows), file.path(project_root, cfg$project$output_dir, "qc", "gwas_harmonization.tsv"))
}
log_message("GWAS harmonization finished.")
