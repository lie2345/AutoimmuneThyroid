#!/usr/bin/env Rscript
# Separate base-R checks against sequential gzip extraction of original sources.
# This is an alternative implementation check, not independent scientific review.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- normalizePath(args[2], mustWork = TRUE)
stopifnot(identical(out, file.path(root, "outputs/t19_t24_20261017/t19")))
options(stringsAsFactors = FALSE, scipen = 0)
log_file <- file(file.path(out, "logs/validation.log"), "wt")
sink(log_file, split = TRUE); sink(log_file, type = "message")
cat("Started UTC: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
cat("Command: Rscript ", shQuote(file.path(out, "code/validate_t19_hla_descriptive.R")), " ",
    shQuote(root), " ", shQuote(out), "\n", sep = "")
read_tab <- function(f, all_char = FALSE) read.delim(f, check.names = FALSE, quote = "", comment.char = "",
    colClasses = if (all_char) "character" else NA, na.strings = "NA")
results <- list()
check <- function(id, pass, details) {
  results[[length(results) + 1L]] <<- data.frame(check_id = id, passed = isTRUE(pass), details = details)
  if (!isTRUE(pass)) stop("Alternative validation failed: ", id, ": ", details)
}
near <- function(x, y, tol = 1e-12) all(is.finite(x) & is.finite(y) & abs(x - y) <= tol * pmax(1, abs(y)))
summary <- read_tab(file.path(out, "tables/hla_within_trait_summary.tsv"))
leads <- read_tab(file.path(out, "tables/hla_regional_lead_variants.tsv"), TRUE)
cross <- read_tab(file.path(out, "tables/hla_cross_trait_at_regional_leads.tsv"), TRUE)
flow <- read_tab(file.path(out, "tables/hla_qc_filter_flow.tsv"))
trait_meta <- read_tab(file.path(root, "config/traits.tsv"))
raws <- list(); filtered <- list()
for (trait in c("E4_GRAVES_STRICT", "E4_HYTHY_AI_STRICT")) {
  raw_file <- file.path(root, "raw_data/finngen_r13/summary_stats", paste0("finngen_R13_", trait, ".gz"))
  seq_file <- file.path(out, "inputs", paste0(trait, ".raw_hla_sequential.tsv"))
  cmd <- paste("gzip -cd --", shQuote(raw_file), "| awk -F '\t'",
               shQuote('NR == 1 || ($1 == "6" && $2 >= 25000000 && $2 <= 34000000)'), ">", shQuote(seq_file))
  cat("Full sequential source extraction command: ", cmd, "\n", sep = "")
  status <- system2("bash", c("-o", "pipefail", "-c", shQuote(cmd)))
  check(paste0(trait, "_sequential_gzip_status"), status == 0L, "Full source body decompressed successfully with pipefail")
  seq <- read_tab(seq_file, TRUE)
  names(seq)[names(seq) == "#chrom"] <- "chrom"
  indexed <- read_tab(file.path(out, "inputs", paste0(trait, ".raw_hla_indexed.tsv")), TRUE)
  check(paste0(trait, "_indexed_equals_sequential"), identical(seq, indexed[names(seq)]),
        paste0(nrow(seq), " source rows and all source fields exactly equal indexed extraction; timestamp warning resolved by body comparison"))
  check(paste0(trait, "_sequential_boundary"), all(seq$chrom == "6" & as.integer(seq$pos) >= 25000000 & as.integer(seq$pos) <= 34000000),
        "Independent sequential extraction uses inclusive region limits")
  rows <- nrow(seq)
  p <- as.numeric(seq$pval); b <- as.numeric(seq$beta); se <- as.numeric(seq$sebeta); af <- as.numeric(seq$af_alt)
  pos <- as.integer(seq$pos)
  mask <- !is.na(pos) & pos > 0 & nzchar(seq$ref) & nzchar(seq$alt) &
          is.finite(p) & p > 0 & p <= 1 & is.finite(b) & is.finite(se) & se > 0 &
          is.finite(af) & af > 0 & af < 1 & seq$chrom %in% as.character(1:22) & pmin(af, 1 - af) >= .005
  expected <- seq[which(mask), ]
  expected$variant_id <- paste0("chr", expected$chrom, "_", expected$pos, "_", toupper(expected$ref), "_", toupper(expected$alt))
  # Match the bytewise variant identity order used by data.table, independent
  # of the user's locale-specific collation of indel strings.
  expected <- expected[order(expected$variant_id, as.numeric(expected$pval), method = "radix"), ]
  expected <- expected[!duplicated(expected$variant_id), ]
  rownames(expected) <- NULL
  retained <- read_tab(file.path(out, "tables", paste0(trait, ".hla_retained.tsv.gz")), TRUE)
  check(paste0(trait, "_all_filter_variant_ids"), identical(expected$variant_id, retained$variant_id),
        paste0(rows, " raw rows -> ", nrow(retained), " established general-QC variants reproduced from sequential raw source"))
  s <- summary[summary$trait_id == trait, ]
  check(paste0(trait, "_count_summary"), s$raw_hla_rows == rows & s$retained_hla_rows == nrow(expected) &
        s$raw_source_p_zero_rows == sum(p == 0, na.rm = TRUE), "Raw, retained and raw source-P-zero counts independently reproduced")
  expected_snp <- nchar(expected$ref) == 1L & nchar(expected$alt) == 1L
  expected_pal <- paste0(toupper(expected$ref), toupper(expected$alt)) %in% c("AT", "TA", "CG", "GC")
  check(paste0(trait, "_variant_type_counts"), s$retained_snp_rows == sum(expected_snp) &
        s$retained_indel_or_other_rows == sum(!expected_snp) & s$retained_palindromic_rows == sum(expected_pal),
        "SNP, indel/other and palindromic counts reproduced; LD-specific filters were not introduced")
  expected_p <- as.numeric(expected$pval)
  check(paste0(trait, "_source_p_threshold_counts"), s$source_p_lt_genomewide == sum(expected_p < 5e-8) &
        s$source_p_lt_suggestive == sum(expected_p < 1e-5), "Threshold counts use original source P, never beta/SE-reconstructed significance")
  check(paste0(trait, "_raw_effect_fields"), identical(expected$alt, retained$effect_allele) &
        identical(expected$ref, retained$other_allele) & near(as.numeric(expected$beta), as.numeric(retained$beta)) &
        near(as.numeric(expected$sebeta), as.numeric(retained$se)), "All retained ALT/ref and beta/SE values equal original sequential raw source")
  check(paste0(trait, "_exact_source_p_and_mlogp"), identical(expected$pval, retained$source_p) & identical(expected$mlogp, retained$source_mlogp),
        "All original P and source mlogp text preserved byte-for-byte at field level")
  expected_or <- exp(as.numeric(expected$beta))
  expected_lo <- exp(log(expected_or) - qnorm(.025, lower.tail = FALSE) * as.numeric(expected$sebeta))
  expected_hi <- exp(log(expected_or) + qnorm(.025, lower.tail = FALSE) * as.numeric(expected$sebeta))
  check(paste0(trait, "_or_and_normal_ci"), near(expected_or, as.numeric(retained$odds_ratio)) &
        near(expected_lo, as.numeric(retained$or_ci95_lower)) & near(expected_hi, as.numeric(retained$or_ci95_upper)),
        "All retained OR=exp(source beta), normal-approximation 95% CI independently recomputed")
  check(paste0(trait, "_ci_order"), all(as.numeric(retained$or_ci95_lower) < as.numeric(retained$odds_ratio) &
        as.numeric(retained$odds_ratio) < as.numeric(retained$or_ci95_upper)), "Every CI brackets its OR")
  lead_expected <- expected[expected_p == min(expected_p), ]
  l <- leads[leads$trait_id == trait, ]
  check(paste0(trait, "_lead_selection"), identical(lead_expected$variant_id, l$variant_id) &
        identical(lead_expected$pval, l$source_p), "All exact minimum-source-P ties retained")
  meta <- trait_meta[trait_meta$trait_id == trait, ]
  check(paste0(trait, "_meta_counts"), s$n_cases == meta$n_cases & s$n_controls == meta$n_controls &
        all(as.integer(retained$n_cases) == meta$n_cases & as.integer(retained$n_controls) == meta$n_controls), "All sample counts agree with config/traits.tsv")
  f <- flow[flow$trait_id == trait, ]
  check(paste0(trait, "_flow_accounting"), all(f$rows_before - f$rows_after == f$removed) &
        all(tail(f$rows_before, -1) == head(f$rows_after, -1)) & tail(f$rows_after, 1) == nrow(expected),
        "Filter-flow rows reconcile without lost records")
  raws[[trait]] <- seq
  filtered[[trait]] <- retained
}
for (j in seq_len(nrow(cross))) {
  c <- cross[j, ]
  x <- filtered[[c$evaluated_trait]]
  x <- x[x$variant_id == c$variant_id, ]
  check(paste0("cross_trait_row_", j), nrow(x) == 1L &
    identical(x$effect_allele, c$effect_allele) & identical(x$other_allele, c$other_allele) &
    identical(x$source_p, c$source_p) & near(as.numeric(x$beta), as.numeric(c$beta)) & near(as.numeric(x$se), as.numeric(c$se)) &
    near(as.numeric(x$odds_ratio), as.numeric(c$odds_ratio)) & near(as.numeric(x$or_ci95_lower), as.numeric(c$or_ci95_lower)) &
    near(as.numeric(x$or_ci95_upper), as.numeric(c$or_ci95_upper)), "Exact same-variant ALT alignment, effects, source P and CI verified")
}
check("cross_trait_row_count", nrow(cross) == nrow(leads) * 2L, "Each selected lead is described in both traits")
result <- do.call(rbind, results)
write.table(result, file.path(out, "validation/alternative_checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat("Alternative implementation checks: ", nrow(result), "/", nrow(result), " PASS\n", sep = "")
cat("Completed UTC: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
cat("These same-agent alternative-code checks are not independent scientific verification. No cross-trait test performed.\n")
sink(type = "message"); sink(); close(log_file)
