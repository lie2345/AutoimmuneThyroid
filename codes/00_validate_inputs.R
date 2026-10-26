#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)
cell_types <- read_cell_types(project_root, cfg)

result_rows <- list()
add_check <- function(component, item, expected, observed, ok, note = "") {
  result_rows[[length(result_rows) + 1L]] <<- data.table::data.table(
    component = component,
    item = item,
    expected = as.character(expected),
    observed = as.character(observed),
    status = if (isTRUE(ok)) "PASS" else "FAIL",
    note = note
  )
}

expected_gwas_header <- c(
  "chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval",
  "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"
)

for (trait_id in traits$trait_id) {
  sumstats <- gwas_path(project_root, trait_id)
  add_check("GWAS", paste0(trait_id, " summary"), "file", file.exists(sumstats), file.exists(sumstats))
  add_check("GWAS", paste0(trait_id, " tabix index"), "file", file.exists(paste0(sumstats, ".tbi")), file.exists(paste0(sumstats, ".tbi")))
  if (file.exists(sumstats)) {
    observed <- read_gzip_header(sumstats)
    add_check(
      "GWAS", paste0(trait_id, " header"), paste(expected_gwas_header, collapse = ","),
      paste(observed, collapse = ","), identical(observed, expected_gwas_header),
      "FinnGen R13 beta is oriented to ALT."
    )
  }
  for (method in c("susie", "finemap")) {
    prefix <- if (method == "susie") "SUSIE" else "FINEMAP"
    fm_dir <- file.path(project_root, "raw_data", "finngen_r13", "finemap", method)
    expected_files <- if (method == "susie") {
      paste0("finngen_R13_", trait_id, ".", prefix, c(".cred.bgz", ".cred_99.bgz", ".snp.bgz", ".snp.bgz.tbi"))
    } else {
      paste0("finngen_R13_", trait_id, ".", prefix, c(".config.bgz", ".region.bgz", ".snp.bgz", ".snp.bgz.tbi"))
    }
    present <- file.exists(file.path(fm_dir, expected_files))
    add_check("GWAS fine-mapping", paste(trait_id, method), length(expected_files), sum(present), all(present))
  }
}

for (modality in c("eQTL", "caQTL")) {
  cis_dir <- file.path(project_root, "raw_data", "finngen_multiome", modality, "cis")
  nominal_dir <- file.path(project_root, "raw_data", "finngen_multiome", modality, "cis_nominal")
  fm_dir <- file.path(project_root, "raw_data", "finngen_multiome", modality, "finemapping")
  cis_files <- list.files(cis_dir, pattern = "\\.cis_acat\\.tsv\\.gz$", full.names = TRUE)
  nominal_files <- list.files(nominal_dir, pattern = "\\.cis_nominal\\.tsv\\.gz$", full.names = TRUE)
  nominal_indexes <- list.files(nominal_dir, pattern = "\\.cis_nominal\\.tsv\\.gz\\.tbi$", full.names = TRUE)
  cred_files <- list.files(fm_dir, pattern = "\\.SUSIE\\.cred\\.tsv\\.gz$", full.names = TRUE)
  cred99_files <- list.files(fm_dir, pattern = "\\.SUSIE\\.cred_99\\.tsv\\.gz$", full.names = TRUE)
  expected_cell_n <- sum(cell_types$qtl_available)
  expected_fm_n <- expected_cell_n * 23L
  add_check("Multiome", paste(modality, "ACAT files"), expected_cell_n, length(cis_files), length(cis_files) == expected_cell_n)
  add_check("Multiome", paste(modality, "nominal files"), expected_cell_n, length(nominal_files), length(nominal_files) == expected_cell_n)
  add_check("Multiome", paste(modality, "nominal indexes"), expected_cell_n, length(nominal_indexes), length(nominal_indexes) == expected_cell_n)
  add_check("Multiome", paste(modality, "95% credible-set files"), expected_fm_n, length(cred_files), length(cred_files) == expected_fm_n)
  add_check("Multiome", paste(modality, "99% credible-set files"), expected_fm_n, length(cred99_files), length(cred99_files) == expected_fm_n)
  add_check(
    "Multiome", paste(modality, "variant-level SuSiE files"), "not in batch 1", 0,
    !any(grepl("\\.SUSIE\\.snp\\.tsv\\.gz$", list.files(fm_dir))),
    "R13 colocalization will use regional nominal QTL statistics; credible-set tables are feature-level QC annotations."
  )
}

required_products <- c(
  "raw_data/finngen_multiome/cascade/finngen_multiome_v1.cascade.variant.tsv.gz",
  "raw_data/finngen_multiome/cascade/finngen_multiome_v1.cascade.gene.tsv.gz",
  "raw_data/finngen_multiome/cascade/finngen_multiome_v1.cascade.peak.tsv.gz",
  "raw_data/finngen_multiome/coloc/finngen_multiome_v1.coloc.FinnGen_R12.eQTL.tsv.gz",
  "raw_data/finngen_multiome/coloc/finngen_multiome_v1.coloc.FinnGen_R12.caQTL.tsv.gz",
  "raw_data/finngen_multiome/atac_peaks/finngen_multiome_v1.consensus_peaks.bed.gz"
)
for (relative in required_products) {
  add_check("Multiome product", relative, "file", file.exists(file.path(project_root, relative)), file.exists(file.path(project_root, relative)))
}

manifest <- file.path(project_root, "raw_data", "manifests", "local_data_manifest.tsv")
add_check("Manifest", "local_data_manifest.tsv", "file", file.exists(manifest), file.exists(manifest))

all_files <- list.files(file.path(project_root, "raw_data"), recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
inventory <- data.table::data.table(
  relative_path = sub(paste0("^", project_root, "/"), "", all_files),
  size_bytes = file.info(all_files)$size,
  modified_time = format(file.info(all_files)$mtime, "%Y-%m-%d %H:%M:%S")
)
inventory[, extension := data.table::fifelse(grepl("\\.tbi$", relative_path), "tbi",
  data.table::fifelse(grepl("\\.(bgz|gz)$", relative_path), "compressed_table", tools::file_ext(relative_path))
)]

checks <- data.table::rbindlist(result_rows, fill = TRUE)
qc_dir <- file.path(project_root, cfg$project$output_dir, "qc")
write_tsv(checks, file.path(qc_dir, "input_validation.tsv"))
write_tsv(inventory, file.path(qc_dir, "raw_file_inventory.tsv"))

if (any(checks$status == "FAIL")) {
  stop("Input validation found failures. Inspect outputs/qc/input_validation.tsv", call. = FALSE)
}
log_message("Input validation completed with ", nrow(checks), " checks.")
