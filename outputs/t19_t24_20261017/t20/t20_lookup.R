#!/usr/bin/env Rscript
# Frozen-panel descriptive lookup only; no GWAS, score, interaction, or case-case fitting.
options(stringsAsFactors = FALSE, digits = 17)
argv <- commandArgs(trailingOnly = TRUE)
mode <- if (length(argv)) argv[[1]] else stop("Specify freeze or lookup")
root <- normalizePath(getwd(), mustWork = TRUE)
out <- file.path(root, "outputs/t19_t24_20261017/t20")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
relative <- function(x) sub(paste0("^", root, "/"), "", x)
tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
read_tsv <- function(path) read.delim(path, colClasses = "character", check.names = FALSE, quote = "", na.strings = "NA")
sha256 <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "[[:space:]]+")[[1]][[1]]
log_command <- function(args) {
  cat(paste(c("tabix", vapply(args, shQuote, character(1))), collapse = " "), "\n", file = file.path(out, "tabix_commands.log"), append = TRUE)
  x <- system2("tabix", vapply(args, shQuote, character(1)), stdout = TRUE, stderr = file.path(out, "tabix_stderr.log"))
  if (!is.null(attr(x, "status")) && attr(x, "status") != 0L) stop("tabix failed")
  x
}
priority_path <- file.path(root, "outcome_data/external_validation/run_20260904/priority_variant_validation.tsv")
panel_path <- file.path(out, "frozen_panel.tsv")
freeze_path <- file.path(out, "freeze_record.tsv")
endpoints <- c("GRAVES_OPHT", "E4_GRAVES_OPHT_STRICT")

if (mode == "freeze") {
  if (file.exists(panel_path) || file.exists(freeze_path)) stop("Refusing to replace an existing frozen panel")
  priority <- read_tsv(priority_path)
  ids <- unique(priority$variant_id)
  stopifnot(length(ids) == 10L, nrow(priority) == 20L, all(table(priority$variant_id) == 2L))
  hla_ids <- c("chr6_32703430_T_G", "chr6_32659200_T_C")
  stopifnot(!any(hla_ids %in% ids))
  ids <- c(ids, hla_ids)
  parts <- strsplit(ids, "_", fixed = TRUE)
  stopifnot(all(lengths(parts) == 4L))
  panel <- data.frame(panel_index = seq_along(ids), variant_id = ids, genome_build = "GRCh38",
    chrom = vapply(parts, function(z) sub("^chr", "", z[[1]]), character(1)),
    pos = vapply(parts, `[[`, character(1), 2L), ref = vapply(parts, `[[`, character(1), 3L),
    alt = vapply(parts, `[[`, character(1), 4L),
    selection_source = c(rep(relative(priority_path), 10L), rep("Primary-agent assignment: two Graves/AIHT marginal HLA leads; T19 confirmation pending at freeze", 2L)),
    selection_class = c(rep("historical_discovery_selected_representative", 10L), "Graves_marginal_HLA_lead", "AIHT_marginal_HLA_lead"),
    nearest_genes_discovery = NA_character_, five_class_discovery = NA_character_,
    graves_discovery_beta_ALT = NA_character_, aiht_discovery_beta_ALT = NA_character_,
    effect_allele = vapply(parts, `[[`, character(1), 4L), other_allele = vapply(parts, `[[`, character(1), 3L))
  for (i in seq_len(10L)) {
    rows <- priority[priority$variant_id == ids[[i]], , drop = FALSE]
    stopifnot(setequal(rows$trait, c("AIHT", "Graves")), all(rows$effect_allele == panel$alt[[i]]), all(rows$other_allele == panel$ref[[i]]))
    panel$nearest_genes_discovery[[i]] <- rows$nearest_genes[[1]]
    panel$five_class_discovery[[i]] <- rows$five_class_discovery[[1]]
    panel$graves_discovery_beta_ALT[[i]] <- rows$discovery_beta[rows$trait == "Graves"]
    panel$aiht_discovery_beta_ALT[[i]] <- rows$discovery_beta[rows$trait == "AIHT"]
  }
  stopifnot(nrow(panel) == 12L, !anyDuplicated(panel$variant_id))
  tsv(panel, panel_path)
  tests <- merge(panel[c("panel_index", "variant_id")], data.frame(endpoint = endpoints), all = TRUE)
  tests <- tests[order(match(tests$endpoint, endpoints), tests$panel_index), ]
  tests$test_index <- seq_len(nrow(tests))
  tests$bonferroni_family_n <- 24L
  tsv(tests, file.path(out, "frozen_tests.tsv"))
  freeze <- data.frame(frozen_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    frozen_before_this_run_eye_extraction = TRUE,
    source_priority_path = relative(priority_path), source_priority_sha256 = sha256(priority_path),
    panel_sha256 = sha256(panel_path), tests_sha256 = sha256(file.path(out, "frozen_tests.tsv")),
    n_variants = 12L, n_endpoints = 2L, n_planned_lookups = 24L)
  tsv(freeze, freeze_path)
  print(freeze)
  quit(status = 0L)
}
if (mode != "lookup") stop("Unknown mode")
freeze <- read_tsv(freeze_path)
stopifnot(sha256(panel_path) == freeze$panel_sha256, sha256(priority_path) == freeze$source_priority_sha256,
  sha256(file.path(out, "frozen_tests.tsv")) == freeze$tests_sha256)
panel <- read_tsv(panel_path)
traits <- read_tsv(file.path(root, "config/traits.tsv"))
manifest <- read_tsv(file.path(root, "raw_data/manifests/finngen_R13_manifest.tsv"))
stopifnot(nrow(panel) == 12L, !anyDuplicated(panel$variant_id))
source_files <- c(priority_path, panel_path, freeze_path, file.path(out, "frozen_tests.tsv"),
  file.path(root, "config/traits.tsv"), file.path(root, "config/analysis.yml"),
  file.path(root, "codes/01_harmonize_gwas.R"), file.path(root, "raw_data/manifests/finngen_R13_manifest.tsv"),
  file.path(out, "t20_lookup.R"),
  file.path(root, "documentation/analysis/PIPELINE.md"),
  file.path(root, "outputs/t19_t24_20261017/sources/source_manifest.json"))
results <- list()
captures <- list()
for (endpoint in endpoints) {
  message("Looking up ", endpoint, " at ", format(Sys.time(), tz = "UTC"))
  input <- file.path(root, "raw_data/finngen_r13/summary_stats", paste0("finngen_R13_", endpoint, ".gz"))
  harmonized <- file.path(root, "outcome_data/harmonized_gwas", paste0(endpoint, ".harmonized.tsv.gz"))
  snapshot_paths <- file.path(root, "outputs/t19_t24_20261017/sources", paste0(endpoint, c("_source.txt", "_visible.txt")))
  stopifnot(all(file.exists(c(input, paste0(input, ".tbi"), harmonized, snapshot_paths))))
  source_files <- c(source_files, input, paste0(input, ".tbi"), harmonized, snapshot_paths)
  header <- log_command(c("-H", input))
  fields <- strsplit(tail(header, 1L), "\t")[[1]]
  stopifnot(all(c("#chrom", "pos", "ref", "alt", "pval", "beta", "sebeta", "af_alt") %in% fields))
  captured <- character()
  for (i in seq_len(nrow(panel))) {
    v <- panel[i, ]
    region <- paste0(v$chrom, ":", v$pos, "-", v$pos)
    lines <- log_command(c(input, region))
    captured <- c(captured, lines)
    queried <- if (length(lines)) read.delim(text = paste(lines, collapse = "\n"), header = FALSE, col.names = fields, colClasses = "character", check.names = FALSE, quote = "") else data.frame()
    exact <- if (nrow(queried)) queried[sub("^chr", "", queried[["#chrom"]]) == v$chrom & queried$pos == v$pos & queried$ref == v$ref & queried$alt == v$alt, , drop = FALSE] else data.frame()
    meta <- traits[traits$trait_id == endpoint, , drop = FALSE]
    stopifnot(nrow(meta) == 1L)
    res <- cbind(v, data.frame(endpoint = endpoint, release = "FinnGen_R13", n_cases = meta$n_cases, n_controls = meta$n_controls,
      raw_path = relative(input), source_query = region, source_position_row_count = nrow(queried), exact_ref_alt_row_count = nrow(exact),
      match_status = "MISSING_EXACT_BUILD_POSITION_REF_ALT", qc_status = "NOT_EVALUABLE_MISSING", qc_reason = "NO_EXACT_ROW",
      beta_source_string = NA_character_, se_source_string = NA_character_, source_p_string = NA_character_, eaf_source_string = NA_character_,
      beta = NA_real_, se = NA_real_, p_source_numeric = NA_real_, eaf_ALT = NA_real_, maf = NA_real_, OR_ALT = NA_real_,
      OR_lower95_normal_approx = NA_real_, OR_upper95_normal_approx = NA_real_, p_bonferroni24_optional_descriptive = NA_real_,
      source_p_below_0_05 = NA, source_p_below_0_05_over_24 = NA, direction_vs_graves_discovery = NA_character_,
      direction_vs_aiht_discovery = NA_character_, rsids_source = NA_character_, nearest_genes_source = NA_character_))
    if (nrow(exact)) {
      # Duplicate handling matches the existing harmonizer: smallest finite source P.
      ord <- order(suppressWarnings(as.numeric(exact$pval)), na.last = TRUE)
      x <- exact[ord[[1]], , drop = FALSE]
      b <- as.numeric(x$beta); s <- as.numeric(x$sebeta); p <- as.numeric(x$pval); f <- as.numeric(x$af_alt)
      problems <- c(if (!is.finite(b)) "BETA_NONFINITE", if (!is.finite(s) || s <= 0) "SE_INVALID",
        if (!is.finite(p) || p <= 0 || p > 1) "P_INVALID", if (!is.finite(f) || f <= 0 || f >= 1) "EAF_INVALID",
        if (!v$chrom %in% as.character(1:22)) "NONAUTOSOMAL",
        if (is.finite(f) && min(f, 1-f) < 0.005) "MAF_LT_0_005")
      res$match_status <- if (nrow(exact) == 1L) "EXACT_GRCH38_REF_ALT" else "EXACT_DUPLICATE_LOWEST_SOURCE_P"
      res$qc_status <- if (length(problems)) "EXCLUDED_GENERAL_QC" else "PASS_GENERAL_QC"
      res$qc_reason <- if (length(problems)) paste(problems, collapse = ";") else "PASS"
      res$beta_source_string <- x$beta; res$se_source_string <- x$sebeta; res$source_p_string <- x$pval; res$eaf_source_string <- x$af_alt
      res$beta <- b; res$se <- s; res$p_source_numeric <- p; res$eaf_ALT <- f; res$maf <- min(f, 1-f)
      res$rsids_source <- x$rsids; res$nearest_genes_source <- x$nearest_genes
      if (!length(problems)) {
        res$OR_ALT <- exp(b); res$OR_lower95_normal_approx <- exp(b-qnorm(0.975)*s); res$OR_upper95_normal_approx <- exp(b+qnorm(0.975)*s)
        res$p_bonferroni24_optional_descriptive <- min(1, 24*p)
        res$source_p_below_0_05 <- p < 0.05; res$source_p_below_0_05_over_24 <- p < 0.05/24
        direction <- function(d) if (is.na(d)) "NOT_ASSESSED_DISCOVERY_CONTEXT_NOT_IN_PANEL_SOURCE" else if (b == 0 || as.numeric(d) == 0) "ZERO_BETA" else if (sign(b) == sign(as.numeric(d))) "same" else "opposite"
        res$direction_vs_graves_discovery <- direction(v$graves_discovery_beta_ALT)
        res$direction_vs_aiht_discovery <- direction(v$aiht_discovery_beta_ALT)
      }
    }
    results[[length(results)+1L]] <- res
  }
  capture_path <- file.path(out, paste0(endpoint, ".raw_position_rows.tsv"))
  writeLines(c(tail(header, 1L), unique(captured)), capture_path)
  captures[[endpoint]] <- capture_path
}
result <- do.call(rbind, results)
stopifnot(nrow(result) == 24L, !anyDuplicated(paste(result$endpoint, result$variant_id)), all(result$effect_allele == result$alt))
tsv(result, file.path(out, "lookup_results.tsv"))
summary <- do.call(rbind, lapply(endpoints, function(ep) {
  z <- result[result$endpoint == ep, ]
  data.frame(endpoint = ep, n_prelisted = nrow(z), n_exact_match = sum(z$exact_ref_alt_row_count > 0),
    n_qc_pass = sum(z$qc_status == "PASS_GENERAL_QC"), n_missing = sum(z$exact_ref_alt_row_count == 0),
    n_excluded = sum(z$qc_status == "EXCLUDED_GENERAL_QC"), n_nominal_source_p_below_0_05 = sum(z$source_p_below_0_05, na.rm = TRUE),
    n_source_p_below_0_05_over_24_optional = sum(z$source_p_below_0_05_over_24, na.rm = TRUE))
}))
tsv(summary, file.path(out, "lookup_summary.tsv"))
source_files <- unique(source_files)
message("Hashing ", length(source_files), " input/code files")
hashes <- data.frame(path = relative(source_files), bytes = file.info(source_files)$size,
  sha256 = vapply(source_files, sha256, character(1)))
tsv(hashes, file.path(out, "source_checksum_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
print(summary)
cat("Task progress: lookup completed; source metadata verified; inference descriptive only.\n")
