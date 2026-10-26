#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)
cc_cfg <- cfg$ccgwas
if (is.null(cc_cfg)) stop("Missing ccgwas configuration", call. = FALSE)

graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
primary_ids <- c(graves_id, aiht_id)
primary <- traits[trait_id %in% primary_ids]
if (nrow(primary) != 2L) stop("Primary trait metadata are incomplete", call. = FALSE)

out_root <- project_path(project_root, cc_cfg$output_dir)
input_dir <- file.path(out_root, "input")
qc_dir <- file.path(out_root, "qc")
ensure_dir(input_dir)
ensure_dir(qc_dir)

hm3_path <- project_path(project_root, cc_cfg$hapmap3_snplist)
if (!file.exists(hm3_path)) stop("Missing HapMap3 SNP list: ", hm3_path, call. = FALSE)
hm3 <- data.table::fread(hm3_path, select = c("SNP", "A1", "A2"))
hm3[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
hm3 <- hm3[grepl("^rs[0-9]+$", SNP) & nchar(A1) == 1L & nchar(A2) == 1L]
hm3 <- hm3[!is_palindromic(A1, A2)]
hm3 <- hm3[!duplicated(SNP)]
data.table::setkey(hm3, SNP)

complement <- c(A = "T", T = "A", C = "G", G = "C")
qc_rows <- list()

record_qc <- function(trait_id, stage, n) {
  qc_rows[[length(qc_rows) + 1L]] <<- data.table::data.table(
    trait_id = trait_id,
    stage = stage,
    n_variants = as.integer(n)
  )
}

prepare_trait <- function(trait_id) {
  path <- gwas_path(project_root, trait_id)
  if (!file.exists(path)) stop("Missing GWAS file: ", path, call. = FALSE)
  log_message("Reading ", trait_id, " from ", path)
  x <- data.table::fread(
    path,
    select = c(
      "#chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
      "pval", "beta", "sebeta", "af_alt_controls"
    ),
    showProgress = TRUE
  )
  data.table::setnames(x, "#chrom", "chrom")
  record_qc(trait_id, "raw_rows", nrow(x))

  x[, `:=`(
    chrom = as.integer(chrom),
    pos = as.integer(pos),
    ref = toupper(ref),
    alt = toupper(alt),
    SNP = first_rsid(rsids)
  )]
  x <- x[
    chrom %in% 1:22 &
      is.finite(pos) & pos > 0 &
      !is_hla_position(chrom, pos, cfg) &
      nchar(ref) == 1L & nchar(alt) == 1L &
      ref %in% names(complement) & alt %in% names(complement)
  ]
  record_qc(trait_id, "autosomal_biallelic_non_hla", nrow(x))

  x <- x[!is.na(SNP) & SNP %in% hm3$SNP]
  record_qc(trait_id, "hapmap3_rsid", nrow(x))
  x <- hm3[x, on = "SNP", nomatch = 0L]

  x[, `:=`(
    direct = alt == A1 & ref == A2,
    swapped = ref == A1 & alt == A2,
    strand_direct = unname(complement[alt]) == A1 & unname(complement[ref]) == A2,
    strand_swapped = unname(complement[ref]) == A1 & unname(complement[alt]) == A2
  )]
  x <- x[direct | swapped | strand_direct | strand_swapped]
  record_qc(trait_id, "alleles_matched_to_hapmap3", nrow(x))

  x[, flip := swapped | strand_swapped]
  x[, `:=`(
    oriented_beta = ifelse(flip, -beta, beta),
    oriented_frequency = ifelse(flip, 1 - af_alt_controls, af_alt_controls)
  )]
  maf_min <- as.numeric(cc_cfg$input_maf_min)
  x <- x[
    is.finite(oriented_beta) &
      is.finite(sebeta) & sebeta > 0 &
      is.finite(pval) & pval >= 0 & pval <= 1 &
      is.finite(oriented_frequency) &
      oriented_frequency > maf_min & oriented_frequency < (1 - maf_min) &
      exp(oriented_beta) >= 0.5 & exp(oriented_beta) <= 2
  ]
  record_qc(trait_id, "ccgwas_numeric_maf_effect_qc", nrow(x))

  duplicate_ids <- x[duplicated(SNP) | duplicated(SNP, fromLast = TRUE), unique(SNP)]
  if (length(duplicate_ids)) x <- x[!SNP %in% duplicate_ids]
  record_qc(trait_id, "unique_rsid", nrow(x))

  requested_trait_id <- trait_id
  metadata <- primary[trait_id == requested_trait_id]
  if (nrow(metadata) != 1L) stop("Missing trait metadata for ", trait_id, call. = FALSE)
  neff <- 4 / (1 / as.numeric(metadata$n_cases) + 1 / as.numeric(metadata$n_controls))
  result <- x[, .(
    SNP,
    CHR = chrom,
    BP = pos,
    EA = A1,
    NEA = A2,
    FRQ = oriented_frequency,
    OR = exp(oriented_beta),
    SE = sebeta,
    P = pmax(pval, .Machine$double.xmin),
    Neff = neff,
    nearest_genes
  )]
  data.table::setkey(result, SNP)
  rm(x)
  invisible(gc())
  result
}

prepared <- lapply(primary_ids, prepare_trait)
names(prepared) <- primary_ids
graves <- prepared[[graves_id]]
aiht <- prepared[[aiht_id]]

matched <- merge(
  graves,
  aiht,
  by = "SNP",
  suffixes = c("_graves", "_aiht"),
  all = FALSE,
  sort = FALSE
)
record_qc("PAIR", "rsid_intersection", nrow(matched))
matched <- matched[
  CHR_graves == CHR_aiht & BP_graves == BP_aiht &
    EA_graves == EA_aiht & NEA_graves == NEA_aiht
]
record_qc("PAIR", "coordinate_and_allele_concordant", nrow(matched))
data.table::setorder(matched, CHR_graves, BP_graves)

make_input <- function(prefix) {
  matched[, .(
    SNP,
    CHR = get(paste0("CHR_", prefix)),
    BP = get(paste0("BP_", prefix)),
    EA = get(paste0("EA_", prefix)),
    NEA = get(paste0("NEA_", prefix)),
    FRQ = get(paste0("FRQ_", prefix)),
    OR = get(paste0("OR_", prefix)),
    SE = get(paste0("SE_", prefix)),
    P = get(paste0("P_", prefix)),
    Neff = get(paste0("Neff_", prefix))
  )]
}

graves_out <- file.path(input_dir, paste0(graves_id, ".ccgwas.tsv.gz"))
aiht_out <- file.path(input_dir, paste0(aiht_id, ".ccgwas.tsv.gz"))
data.table::fwrite(make_input("graves"), graves_out, sep = "\t", quote = FALSE, compress = "gzip")
data.table::fwrite(make_input("aiht"), aiht_out, sep = "\t", quote = FALSE, compress = "gzip")

annotation <- matched[, .(
  SNP,
  CHR = CHR_graves,
  BP = BP_graves,
  EA = EA_graves,
  NEA = NEA_graves,
  nearest_genes_graves,
  nearest_genes_aiht
)]
data.table::fwrite(
  annotation,
  file.path(input_dir, "matched_variant_annotation.tsv.gz"),
  sep = "\t", quote = FALSE, compress = "gzip"
)

qc <- data.table::rbindlist(qc_rows)
qc[, fraction_of_raw := n_variants / n_variants[stage == "raw_rows"][1L], by = trait_id]
qc[trait_id == "PAIR", fraction_of_raw := NA_real_]
data.table::fwrite(qc, file.path(qc_dir, "input_qc_counts.tsv"), sep = "\t", quote = FALSE)

manifest <- data.table::data.table(
  trait_id = primary_ids,
  n_cases = primary[match(primary_ids, trait_id), n_cases],
  n_controls = primary[match(primary_ids, trait_id), n_controls],
  n_variants = nrow(matched),
  neff = c(unique(make_input("graves")$Neff), unique(make_input("aiht")$Neff)),
  genome_build = cfg$project$genome_build,
  variant_scope = "HapMap3, MAF>0.01, non-palindromic, autosomal, non-HLA",
  effect_allele = "HapMap3 A1",
  file = c(graves_out, aiht_out)
)
data.table::fwrite(manifest, file.path(input_dir, "input_manifest.tsv"), sep = "\t", quote = FALSE)

log_message("Prepared ", format(nrow(matched), big.mark = ","), " matched SNPs for CC-GWAS")
log_message("Graves input: ", graves_out)
log_message("AIHT input: ", aiht_out)
