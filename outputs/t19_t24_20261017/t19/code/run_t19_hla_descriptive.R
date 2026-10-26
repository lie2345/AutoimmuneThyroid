#!/usr/bin/env Rscript
# T19 Path B: descriptive HLA-region association reporting only.
# The established general GWAS preprocessing in codes/01_harmonize_gwas.R is
# reproduced within the configured HLA interval and checked against frozen data.
options(stringsAsFactors = FALSE, scipen = 0)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript run_t19_hla_descriptive.R PROJECT_ROOT OUTPUT_DIR")
root <- normalizePath(args[1], mustWork = TRUE)
out <- normalizePath(args[2], mustWork = TRUE)
allowed <- file.path(root, "outputs/t19_t24_20261017/t19")
stopifnot(identical(out, allowed))
for (p in c("data.table", "yaml", "Rsamtools", "GenomicRanges", "IRanges", "digest")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Required installed package unavailable: ", p)
}
library(data.table)
setDTthreads(2L)
source(file.path(root, "codes/lib/common.R"))
options(scipen = 0)
path <- function(folder, name) file.path(out, folder, name)
write_tab <- function(x, name, folder = "tables") {
  fwrite(x, path(folder, name), sep = "\t", quote = FALSE, na = "NA")
}
log_con <- file(path("logs", "execution.log"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
cat("Started UTC: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
cat("Command: Rscript ", shQuote(file.path(out, "code/run_t19_hla_descriptive.R")),
    " ", shQuote(root), " ", shQuote(out), "\n", sep = "")
cat("Task progress: running; source verification: in progress; inference: descriptive associations\n")
cfg <- read_analysis_config(root)
traits <- read_traits(root, cfg)[trait_id %in% c("E4_GRAVES_STRICT", "E4_HYTHY_AI_STRICT")]
stopifnot(nrow(traits) == 2L, cfg$project$genome_build == "GRCh38",
          cfg$hla$chromosome == 6, cfg$hla$start == 25000000, cfg$hla$end == 34000000)
capture.output(sessionInfo(), file = path("logs", "sessionInfo.txt"))
pkgs <- c("data.table", "yaml", "Rsamtools", "GenomicRanges", "IRanges", "digest")
write_tab(data.table(package = pkgs, version = vapply(pkgs, function(p) as.character(packageVersion(p)), "")),
          "package_versions.tsv", "logs")
sha <- function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE)
input_files <- c("config/analysis.yml", "config/traits.tsv", "codes/01_harmonize_gwas.R", "codes/lib/common.R",
                 "raw_data/manifests/local_data_manifest.tsv", "raw_data/manifests/finngen_R13_manifest.tsv",
                 "outputs/qc/gwas_harmonization.tsv", "documentation/analysis/PIPELINE.md")
for (t in traits$trait_id) input_files <- c(input_files,
  sub(paste0("^", root, "/"), "", gwas_path(root, t)),
  paste0(sub(paste0("^", root, "/"), "", gwas_path(root, t)), ".tbi"),
  paste0("outcome_data/harmonized_gwas/", t, ".harmonized.tsv.gz"))
input_hashes <- data.table(path = input_files)
input_hashes[, `:=`(bytes = file.info(file.path(root, path))$size,
                    sha256 = vapply(file.path(root, path), sha, ""))]
write_tab(input_hashes, "input_sha256.tsv", "validation")
# Preserve every pre-existing file under outputs/ except this new multi-task run.
# Other agents own siblings in t19_t24_20261017 and may legitimately create them.
old_files <- list.files(file.path(root, "outputs"), recursive = TRUE, full.names = TRUE)
old_files <- old_files[!startsWith(old_files, file.path(root, "outputs/t19_t24_20261017/"))]
old_hashes <- data.table(path = sub(paste0("^", root, "/"), "", old_files),
                         sha256 = vapply(old_files, sha, ""))
write_tab(old_hashes, "preexisting_output_sha256_before.tsv", "validation")
manifest <- fread(file.path(root, "raw_data/manifests/local_data_manifest.tsv"))
checks <- list()
check <- function(id, pass, details) {
  checks[[length(checks) + 1L]] <<- data.table(check_id = id, passed = isTRUE(pass), details = details)
  if (!isTRUE(pass)) stop("Validation failed: ", id, ": ", details)
}
qc <- list(); summary_rows <- list(); all_regions <- list(); leads <- list(); provenance <- list()
for (i in seq_len(nrow(traits))) {
  meta <- traits[i]
  trait <- meta$trait_id
  raw_file <- gwas_path(root, trait)
  cat("Processing ", trait, "\n", sep = "")
  header <- read_gzip_header(raw_file)
  query <- GenomicRanges::GRanges(seqnames = "6", ranges = IRanges::IRanges(cfg$hla$start, cfg$hla$end))
  lines <- Rsamtools::scanTabix(Rsamtools::TabixFile(raw_file), param = query)[[1]]
  raw <- fread(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t", colClasses = "character")
  setnames(raw, header)
  raw[, variant_id := canonical_variant_id(chrom, as.integer(pos), ref, alt)]
  write_tab(raw, paste0(trait, ".raw_hla_indexed.tsv"), "inputs")
  check(paste0(trait, "_raw_boundary"), all(raw$chrom == "6" & as.integer(raw$pos) >= cfg$hla$start & as.integer(raw$pos) <= cfg$hla$end),
        paste0("Indexed raw rows=", nrow(raw), "; inclusive chr6:", cfg$hla$start, "-", cfg$hla$end))
  x <- copy(raw)
  x[, `:=`(chrom = sub("^chr", "", chrom, ignore.case = TRUE), pos = as.integer(pos),
           ref = toupper(ref), alt = toupper(alt), rsid = first_rsid(rsids),
           p = as.numeric(pval), se = as.numeric(sebeta), eaf = as.numeric(af_alt), beta = as.numeric(beta))]
  n_raw <- nrow(x)
  filter_counts <- data.table(trait_id = trait, step = "raw_hla", rows_before = n_raw, rows_after = n_raw, removed = 0L)
  apply_filter <- function(label, keep) {
    before <- nrow(x)
    x <<- x[which(keep)]
    filter_counts <<- rbind(filter_counts, data.table(trait_id = trait, step = label,
                                                   rows_before = before, rows_after = nrow(x), removed = before - nrow(x)))
  }
  apply_filter("position_and_nonempty_alleles", !is.na(x$pos) & x$pos > 0 & nzchar(x$ref) & nzchar(x$alt))
  apply_filter("source_p_finite_gt0_le1", is.finite(x$p) & x$p > 0 & x$p <= 1)
  apply_filter("finite_beta_and_se_gt0", is.finite(x$beta) & is.finite(x$se) & x$se > 0)
  apply_filter("finite_eaf_gt0_lt1", is.finite(x$eaf) & x$eaf > 0 & x$eaf < 1)
  if (isTRUE(cfg$gwas_qc$autosomes_only)) apply_filter("autosomes", x$chrom %in% as.character(1:22))
  x[, maf := pmin(eaf, 1 - eaf)]
  apply_filter("general_maf_min", x$maf >= cfg$gwas_qc$maf_min_general)
  x[, `:=`(variant_id = canonical_variant_id(chrom, pos, ref, alt), effect_allele = alt,
            other_allele = ref, n_cases = meta$n_cases, n_controls = meta$n_controls,
            is_snp = nchar(ref) == 1L & nchar(alt) == 1L, palindromic = is_palindromic(ref, alt))]
  setorder(x, variant_id, p)
  apply_filter("deduplicate_variant_id_min_source_p", !duplicated(x$variant_id))
  x[, `:=`(trait_id = trait, source_p = pval, source_mlogp = mlogp,
            source_p_status = fifelse(as.numeric(pval) == 0, "SOURCE_ZERO_UNDERFLOW_EXCLUDED", "SOURCE_POSITIVE"),
            odds_ratio = exp(beta), or_ci95_lower = exp(beta - qnorm(.975) * se),
            or_ci95_upper = exp(beta + qnorm(.975) * se))]
  keep <- c("trait_id", "variant_id", "chrom", "pos", "rsid", "ref", "alt", "effect_allele", "other_allele",
            "beta", "se", "odds_ratio", "or_ci95_lower", "or_ci95_upper", "source_p", "source_mlogp",
            "source_p_status", "eaf", "maf", "n_cases", "n_controls", "is_snp", "palindromic")
  write_tab(x[, ..keep], paste0(trait, ".hla_retained.tsv.gz"))
  # Extract the frozen harmonized HLA records without loading a whole GWAS in R.
  harm_file <- file.path(root, "outcome_data/harmonized_gwas", paste0(trait, ".harmonized.tsv.gz"))
  harm_extract <- path("inputs", paste0(trait, ".frozen_harmonized_hla.tsv"))
  awk <- sprintf('NR == 1 || ($2 == "%s" && $3 >= %d && $3 <= %d)', cfg$hla$chromosome, cfg$hla$start, cfg$hla$end)
  cmd <- paste("gzip -cd --", shQuote(harm_file), "| awk -F '\t'", shQuote(awk), ">", shQuote(harm_extract))
  cat("Extract command: ", cmd, "\n", sep = "")
  stopifnot(system(cmd) == 0)
  h <- fread(harm_extract)
  setorder(h, variant_id)
  setorder(x, variant_id)
  check(paste0(trait, "_frozen_variant_set"), identical(x$variant_id, h$variant_id), paste0(nrow(x), " retained HLA rows equal frozen harmonized HLA rows"))
  fields <- c("chrom", "pos", "ref", "alt", "effect_allele", "other_allele", "beta", "se", "p", "eaf", "maf", "n_cases", "n_controls", "is_snp", "palindromic")
  for (f in fields) check(paste0(trait, "_frozen_", f),
    if (is.numeric(x[[f]]) || is.numeric(h[[f]])) isTRUE(all.equal(as.numeric(x[[f]]), as.numeric(h[[f]]), tolerance = 2e-14)) else identical(as.character(x[[f]]), as.character(h[[f]])),
    paste0("All ", nrow(x), " rows checked against frozen harmonized ", f))
  check(paste0(trait, "_alt_alignment"), all(x$effect_allele == x$alt & x$other_allele == x$ref), "Effects remain source ALT-oriented")
  check(paste0(trait, "_source_sample_counts"), all(h$n_cases == meta$n_cases & h$n_controls == meta$n_controls), "No sample count changes")
  source_manifest <- manifest[filename == basename(raw_file)]
  stopifnot(nrow(source_manifest) == 1L)
  actual_md5 <- unname(tools::md5sum(raw_file))
  check(paste0(trait, "_manifest_md5"), actual_md5 == source_manifest$checksum, "Complete raw body MD5 equals frozen local manifest")
  provenance[[i]] <- source_manifest[, .(trait_id = trait_or_celltype, release, filename, source_url = url_or_gcs_path,
     genome_build, effect_allele, download_date, manifest_md5 = checksum, actual_md5 = actual_md5,
     raw_path = sub(paste0("^", root, "/"), "", raw_file),
     harmonized_path = sub(paste0("^", root, "/"), "", harm_file))]
  qc[[i]] <- filter_counts
  summary_rows[[i]] <- data.table(trait_id = trait, label = meta$label, genome_build = cfg$project$genome_build,
    hla_chr = cfg$hla$chromosome, hla_start_inclusive = cfg$hla$start, hla_end_inclusive = cfg$hla$end,
    n_cases = meta$n_cases, n_controls = meta$n_controls, raw_hla_rows = n_raw,
    raw_source_p_zero_rows = sum(as.numeric(raw$pval) == 0, na.rm = TRUE), retained_hla_rows = nrow(x),
    retained_snp_rows = sum(x$is_snp), retained_indel_or_other_rows = sum(!x$is_snp),
    retained_palindromic_rows = sum(x$palindromic), general_maf_min = cfg$gwas_qc$maf_min_general,
    source_p_lt_genomewide = sum(x$p < cfg$gwas_qc$genome_wide_p), genomewide_threshold = cfg$gwas_qc$genome_wide_p,
    source_p_lt_suggestive = sum(x$p < cfg$gwas_qc$suggestive_p), suggestive_threshold = cfg$gwas_qc$suggestive_p,
    min_source_p = min(x$p), min_source_p_text = x[which.min(p), source_p],
    retained_min_pos = min(x$pos), retained_max_pos = max(x$pos),
    rows_at_start_boundary = sum(x$pos == cfg$hla$start), rows_at_end_boundary = sum(x$pos == cfg$hla$end))
  lead <- x[p == min(p), ..keep]
  lead[, `:=`(lead_selection = "Minimum positive source P among established general-QC retained HLA variants; all exact ties retained",
               min_p_ties = .N)]
  leads[[i]] <- lead
  all_regions[[trait]] <- x
  cat("Raw rows: ", n_raw, "; retained: ", nrow(x), "; source-P lead: ", paste(lead$variant_id, collapse = ","), "\n", sep = "")
}
write_tab(rbindlist(qc), "hla_qc_filter_flow.tsv")
write_tab(rbindlist(summary_rows), "hla_within_trait_summary.tsv")
lead_table <- rbindlist(leads)
write_tab(lead_table, "hla_regional_lead_variants.tsv")
write_tab(rbindlist(provenance), "input_provenance.tsv", "validation")
cross <- list()
for (j in seq_len(nrow(lead_table))) {
  l <- lead_table[j]
  for (target in names(all_regions)) {
    a <- all_regions[[target]][variant_id == l$variant_id]
    stopifnot(nrow(a) <= 1L)
    if (nrow(a)) {
      cross[[length(cross) + 1L]] <- data.table(selected_as_lead_for = l$trait_id, evaluated_trait = target,
        variant_id = l$variant_id, rsid = a$rsid, effect_allele = a$effect_allele, other_allele = a$other_allele,
        beta = a$beta, se = a$se, odds_ratio = a$odds_ratio, or_ci95_lower = a$or_ci95_lower, or_ci95_upper = a$or_ci95_upper,
        source_p = a$source_p, source_mlogp = a$source_mlogp, eaf = a$eaf, maf = a$maf,
        n_cases = a$n_cases, n_controls = a$n_controls, matching_status = "EXACT_GRCh38_CHR_POS_REF_ALT",
        effect_direction_relative_to_selecting_trait = if (sign(a$beta) == sign(l$beta)) "same_sign" else "different_sign",
        interpretation = "Descriptive same-variant effects only; no cross-trait effect test or independent-signal claim")
    } else {
      cross[[length(cross) + 1L]] <- data.table(selected_as_lead_for = l$trait_id, evaluated_trait = target,
        variant_id = l$variant_id, matching_status = "NOT_IN_GENERAL_QC_RETAINED_HLA_SET")
    }
  }
}
write_tab(rbindlist(cross, fill = TRUE), "hla_cross_trait_at_regional_leads.tsv")
check("input_sha256_unchanged", identical(input_hashes$sha256, unname(vapply(file.path(root, input_hashes$path), sha, ""))),
      "All authoritative input bodies/configuration/code remain unchanged")
old_after <- vapply(file.path(root, old_hashes$path), sha, "")
check("preexisting_outputs_unchanged", identical(old_hashes$sha256, unname(old_after)),
      paste0(length(old_after), " pre-existing files under outputs/ outside the new multi-task run hashed before/after"))
write_tab(data.table(path = old_hashes$path, sha256 = unname(old_after)), "preexisting_output_sha256_after.tsv", "validation")
write_tab(rbindlist(checks), "implementation_checks.tsv", "validation")
cat("Implementation checks: ", length(checks), "/", length(checks), " PASS\n", sep = "")
cat("Completed UTC: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
cat("Task progress: completed implementation; source verification: local frozen inputs verified; inference: descriptive associations\n")
cat("Separate validation code must pass before handoff. No classical HLA, amino-acid, haplotype, conditional, case-case, gene-causal, or independent-signal inference performed. T18 remains HOLD.\n")
sink(type = "message"); sink(); close(log_con)
