#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)[include_global_rg == TRUE]
args <- parse_cli()
requested <- split_csv(args$traits %||% "")
if (length(requested)) traits <- traits[trait_id %in% requested]
if (!nrow(traits)) stop("No traits selected for LDSC/LAVA preparation.", call. = FALSE)

harmonized_dir <- file.path(project_root, cfg$project$outcome_dir, "harmonized_gwas")
ldsc_dir <- file.path(project_root, cfg$project$outcome_dir, "ldsc", "input")
lava_dir <- file.path(project_root, cfg$project$outcome_dir, "lava", "input")
ensure_dir(ldsc_dir)
ensure_dir(lava_dir)
info_rows <- list()
qc_rows <- list()

for (i in seq_len(nrow(traits))) {
  meta <- traits[i]
  source_path <- file.path(harmonized_dir, paste0(meta$trait_id, ".harmonized.tsv.gz"))
  if (!file.exists(source_path)) stop("Missing harmonized GWAS: ", source_path, call. = FALSE)
  log_message("Preparing LD-method input for ", meta$trait_id)
  x <- data.table::fread(
    source_path,
    select = c("rsid", "ref", "alt", "beta", "z", "p", "eaf", "maf", "n_total", "is_snp", "palindromic", "is_hla"),
    showProgress = TRUE
  )
  x <- x[
    is_snp & !is.na(rsid) & grepl("^rs[0-9]+$", rsid) &
      maf >= as.numeric(cfg$gwas_qc$maf_min_ld_methods) & !is_hla
  ]
  if (isTRUE(cfg$gwas_qc$remove_palindromic_for_ld_methods)) x <- x[palindromic == FALSE]
  data.table::setorder(x, rsid, p)
  x <- x[!duplicated(rsid)]

  ldsc <- x[, .(
    SNP = rsid, A1 = alt, A2 = ref, P = p, BETA = beta, N = n_total,
    N_CAS = as.integer(meta$n_cases), N_CON = as.integer(meta$n_controls), FRQ = eaf
  )]
  ldsc_path <- file.path(ldsc_dir, paste0(meta$trait_id, ".ldsc_input.tsv.gz"))
  data.table::fwrite(ldsc, ldsc_path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

  lava <- x[, .(SNP = rsid, A1 = alt, A2 = ref, Z = z, N = n_total)]
  lava_path <- file.path(lava_dir, paste0(meta$trait_id, ".lava.stats"))
  data.table::fwrite(lava, lava_path, sep = "\t", quote = FALSE, na = "NA")
  info_rows[[length(info_rows) + 1L]] <- data.table::data.table(
    phenotype = meta$trait_id,
    cases = as.integer(meta$n_cases),
    controls = as.integer(meta$n_controls),
    filename = lava_path
  )
  qc_rows[[length(qc_rows) + 1L]] <- data.table::data.table(
    trait_id = meta$trait_id,
    variants_for_ld_methods = nrow(x),
    ldsc_input = ldsc_path,
    lava_input = lava_path
  )
  rm(x, ldsc, lava)
  invisible(gc())
}

write_tsv(data.table::rbindlist(info_rows), file.path(lava_dir, "input.info.tsv"))
write_tsv(data.table::rbindlist(qc_rows), file.path(project_root, cfg$project$output_dir, "qc", "ld_method_preparation.tsv"))
log_message("Prepared LDSC and LAVA inputs. LAVA remains blocked until a build-matched LD reference and an LDSC-derived sample-overlap matrix are configured.")
