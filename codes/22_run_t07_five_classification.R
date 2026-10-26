#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "Rsamtools", "GenomicRanges", "IRanges"))
cfg <- read_analysis_config(project_root)
if (is.null(cfg$t07_five_classification)) stop("Missing t07_five_classification config", call. = FALSE)
t07 <- cfg$t07_five_classification
graves <- cfg$primary_comparison$graves
aiht <- cfg$primary_comparison$autoimmune_hypothyroidism
out_dir <- file.path(project_root, cfg$project$outcome_dir, "t07")
ensure_dir(out_dir)

t06_dir <- file.path(project_root, cfg$project$outcome_dir, "t06")
signals <- data.table::fread(file.path(t06_dir, "susie_signal_catalog.tsv"))
members <- data.table::fread(file.path(t06_dir, "susie_credible_set_members.tsv.gz"))
coloc <- data.table::fread(file.path(t06_dir, "gwas_gwas_susie_coloc_annotated.tsv"))
cc_stability <- data.table::fread(file.path(
  project_root, cfg$project$outcome_dir, "ccgwas", "summary", "signal_stability.tsv"
))
cc_annotation <- data.table::fread(file.path(
  project_root, cfg$project$outcome_dir, "ccgwas", "input", "matched_variant_annotation.tsv.gz"
))
cc_loci <- data.table::fread(file.path(
  project_root, cfg$project$outcome_dir, "ccgwas", "summary", "provisional_lead_loci.tsv"
))

posterior_threshold <- as.numeric(t07$classification_posterior)
suggestive_h4_threshold <- as.numeric(t07$block_suggestive_h4)
equivalence_alpha <- as.numeric(t07$equivalence_alpha)
primary_or_upper <- as.numeric(t07$equivalence_or_upper_primary)
sensitivity_or_upper <- as.numeric(unlist(t07$equivalence_or_upper_sensitivity))
or_bounds <- sort(unique(c(sensitivity_or_upper, primary_or_upper)))
case_difference_p <- as.numeric(t07$case_difference_p)
z_equivalence <- stats::qnorm(1 - equivalence_alpha)

signal_id <- function(trait_id, region, cs) paste(trait_id, region, cs, sep = "|")
signals[, source_signal_id := signal_id(trait_id, region, cs)]
members[, source_signal_id := signal_id(trait_id, region, cs)]
signals[, source_trait_label := ifelse(trait_id == graves, "Graves", "AIHT")]

parse_variant <- function(x) {
  pieces <- regmatches(x, regexec("^chr([^_]+)_([0-9]+)_([^_]+)_([^_]+)$", x))
  valid <- lengths(pieces) == 5L
  data.table::data.table(
    variant_id = x,
    chromosome = ifelse(valid, vapply(pieces, function(z) z[[2L]] %||% NA_character_, character(1)), NA_character_),
    position = ifelse(valid, as.integer(vapply(pieces, function(z) z[[3L]] %||% NA_character_, character(1))), NA_integer_),
    ref = ifelse(valid, vapply(pieces, function(z) z[[4L]] %||% NA_character_, character(1)), NA_character_),
    alt = ifelse(valid, vapply(pieces, function(z) z[[5L]] %||% NA_character_, character(1)), NA_character_)
  )
}

annotate_variants <- function(variant_ids, trait_id) {
  variants <- unique(parse_variant(unique(variant_ids[!is.na(variant_ids)])))
  path <- gwas_path(project_root, trait_id)
  header <- read_gzip_header(path)
  tabix_file <- Rsamtools::TabixFile(path)
  available <- suppressWarnings(Rsamtools::seqnamesTabix(tabix_file))
  suppressWarnings(Rsamtools::open.TabixFile(tabix_file))
  on.exit(Rsamtools::close.TabixFile(tabix_file), add = TRUE)
  ans <- lapply(seq_len(nrow(variants)), function(i) {
    v <- variants[i]
    seqname <- c(paste0("chr", v$chromosome), v$chromosome)
    seqname <- seqname[seqname %in% available][1L]
    empty <- data.table::data.table(
      variant_id = v$variant_id, beta = NA_real_, se = NA_real_, p = NA_real_,
      af_alt = NA_real_, nearest_genes = NA_character_, reported_rsids = NA_character_,
      annotation_status = if (is.na(seqname)) "chromosome_absent_from_tabix_index" else "variant_not_found"
    )
    if (is.na(seqname)) return(empty)
    param <- GenomicRanges::GRanges(seqnames = seqname, ranges = IRanges::IRanges(v$position, v$position))
    lines <- suppressWarnings(Rsamtools::scanTabix(tabix_file, param = param)[[1L]])
    if (!length(lines)) return(empty)
    x <- data.table::fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE, showProgress = FALSE)
    data.table::setnames(x, header)
    x <- x[pos == v$position & toupper(ref) == toupper(v$ref) & toupper(alt) == toupper(v$alt)]
    if (!nrow(x)) return(empty)
    data.table::data.table(
      variant_id = v$variant_id,
      beta = as.numeric(x$beta[[1L]]), se = as.numeric(x$sebeta[[1L]]),
      p = as.numeric(x$pval[[1L]]), af_alt = as.numeric(x$af_alt[[1L]]),
      nearest_genes = as.character(x$nearest_genes[[1L]]),
      reported_rsids = as.character(x$rsids[[1L]]), annotation_status = "matched"
    )
  })
  data.table::rbindlist(ans, fill = TRUE)
}

# Strong H4 pairs are the only cross-trait signals collapsed into one unit.
shared_pairs <- coloc[classification == "shared_strong_H4"]
shared_pairs[, `:=`(
  graves_source_signal_id = signal_id(graves, graves_region, graves_cs),
  aiht_source_signal_id = signal_id(aiht, aiht_region, aiht_cs),
  conservative_posterior = pmin(PP.H4, sensitivity_PP.H4)
)]
if (shared_pairs[, anyDuplicated(graves_source_signal_id)] || shared_pairs[, anyDuplicated(aiht_source_signal_id)]) {
  stop("A source signal is assigned to more than one strong H4 pair", call. = FALSE)
}
matched_source_ids <- c(shared_pairs$graves_source_signal_id, shared_pairs$aiht_source_signal_id)

# Create a signal-pair view from the perspective of every source signal.
pair_long <- data.table::rbindlist(list(
  coloc[, .(
    source_signal_id = signal_id(graves, graves_region, graves_cs),
    counterpart_source_signal_id = signal_id(aiht, aiht_region, aiht_cs),
    pair_id, PP.H3, PP.H4, sensitivity_PP.H3, sensitivity_PP.H4,
    quality_h3, quality_suggestive, pair_classification = classification
  )],
  coloc[, .(
    source_signal_id = signal_id(aiht, aiht_region, aiht_cs),
    counterpart_source_signal_id = signal_id(graves, graves_region, graves_cs),
    pair_id, PP.H3, PP.H4, sensitivity_PP.H3, sensitivity_PP.H4,
    quality_h3, quality_suggestive, pair_classification = classification
  )]
))
pair_long[, `:=`(
  conservative_h3 = pmin(PP.H3, sensitivity_PP.H3),
  conservative_h4 = pmin(PP.H4, sensitivity_PP.H4)
)]
pair_evidence <- pair_long[, {
  h3_index <- which.max(conservative_h3)
  eligible_h4 <- conservative_h4[as.logical(quality_suggestive) & is.finite(conservative_h4)]
  .(
    n_direct_signal_pairs = .N,
    best_h3_pair_id = pair_id[h3_index],
    best_h3_counterpart = counterpart_source_signal_id[h3_index],
    best_h3_primary = PP.H3[h3_index],
    best_h3_sensitivity = sensitivity_PP.H3[h3_index],
    best_h3_conservative = conservative_h3[h3_index],
    best_h3_quality = as.logical(quality_h3[h3_index]),
    max_h4_conservative = if (length(eligible_h4)) max(eligible_h4) else NA_real_,
    suggestive_shared_model = any(
      as.logical(quality_suggestive) & conservative_h4 >= suggestive_h4_threshold,
      na.rm = TRUE
    )
  )
}, by = source_signal_id]
pair_evidence[!is.finite(max_h4_conservative), max_h4_conservative := NA_real_]

# Map robust formal CC-GWAS SNPs directly into each 95% credible set.
cc <- merge(
  cc_stability[final_sensitivity_status == "robust_candidate"],
  cc_annotation,
  by = c("SNP", "CHR", "BP", "EA", "NEA")
)
cc[, ccgwas_locus_id := NA_character_]
for (i in seq_len(nrow(cc_loci))) {
  loc <- cc_loci[i]
  cc[CHR == loc$CHR & BP >= loc$start & BP <= loc$end, ccgwas_locus_id := loc$locus_id]
}
cc_small <- cc[, .(
  chromosome = paste0("chr", CHR), position = BP, SNP, EA, NEA,
  ccgwas_locus_id, primary_OLS_beta, primary_OLS_pval,
  primary_Exact_beta, primary_Exact_pval
)]
member_cc <- merge(members, cc_small, by = c("chromosome", "position"), allow.cartesian = TRUE)
member_cc[, allele_match := (toupper(allele1) == toupper(EA) & toupper(allele2) == toupper(NEA)) |
  (toupper(allele1) == toupper(NEA) & toupper(allele2) == toupper(EA))]
member_cc <- member_cc[allele_match == TRUE]
member_cc[, ccgwas_beta_alt := data.table::fifelse(
  toupper(EA) == toupper(allele2), primary_OLS_beta,
  data.table::fifelse(toupper(NEA) == toupper(allele2), -primary_OLS_beta, NA_real_)
)]
member_cc[, ccgwas_direction_consistent := data.table::fifelse(
  trait_id == graves,
  ccgwas_beta_alt * marginal_beta > 0,
  ccgwas_beta_alt * marginal_beta < 0
)]
cc_support <- member_cc[, .(
  ccgwas_robust_member_n = data.table::uniqueN(SNP),
  ccgwas_robust_member_pip = sum(cs_specific_pip_recomputed, na.rm = TRUE),
  ccgwas_direction_consistent_pip = sum(
    cs_specific_pip_recomputed[ccgwas_direction_consistent %in% TRUE], na.rm = TRUE
  ),
  ccgwas_direction_consistency = {
    denom <- sum(cs_specific_pip_recomputed, na.rm = TRUE)
    if (denom > 0) sum(cs_specific_pip_recomputed[ccgwas_direction_consistent %in% TRUE], na.rm = TRUE) / denom else NA_real_
  },
  ccgwas_locus_ids = paste(sort(unique(stats::na.omit(ccgwas_locus_id))), collapse = ";")
), by = source_signal_id]

# Build 304 signal units: eight strong-H4 pairs are collapsed; all other source
# signals remain separate until an external dataset can resolve them.
shared_units <- shared_pairs[, .(
  unit_kind = "strong_H4_pair",
  graves_source_signal_id, aiht_source_signal_id,
  source_signal_id = NA_character_, source_trait_label = "Graves+AIHT",
  representative_variant = top_h4_variant,
  evidence_pair_id = pair_id,
  evidence_counterpart_signal_id = NA_character_,
  posterior_type = "H4",
  posterior_primary = PP.H4,
  posterior_sensitivity = sensitivity_PP.H4,
  classification_posterior = conservative_posterior,
  preliminary_shared_direction = effect_direction_if_shared
)]
unmatched <- signals[!source_signal_id %in% matched_source_ids]
unmatched_units <- unmatched[, .(
  unit_kind = "single_trait_signal",
  graves_source_signal_id = ifelse(trait_id == graves, source_signal_id, NA_character_),
  aiht_source_signal_id = ifelse(trait_id == aiht, source_signal_id, NA_character_),
  source_signal_id, source_trait_label,
  representative_variant = lead_variant,
  evidence_pair_id = NA_character_, evidence_counterpart_signal_id = NA_character_,
  posterior_type = "H3",
  posterior_primary = NA_real_, posterior_sensitivity = NA_real_,
  classification_posterior = NA_real_, preliminary_shared_direction = NA_character_
)]
units <- data.table::rbindlist(list(shared_units, unmatched_units), fill = TRUE)
parsed_units <- unique(parse_variant(unique(units$representative_variant)))
units <- merge(units, parsed_units, by.x = "representative_variant", by.y = "variant_id", all.x = TRUE)
units[, chromosome_order := match(chromosome, c(as.character(1:22), "X"))]
data.table::setorder(units, chromosome_order, position, unit_kind)
units[, chromosome_order := NULL]
units[, signal_unit_id := sprintf("GSU%03d", seq_len(.N))]

signal_unit_members <- data.table::rbindlist(list(
  units[unit_kind == "strong_H4_pair", .(
    signal_unit_id, source_signal_id = graves_source_signal_id, membership_role = "Graves member of strong H4 pair"
  )],
  units[unit_kind == "strong_H4_pair", .(
    signal_unit_id, source_signal_id = aiht_source_signal_id, membership_role = "AIHT member of strong H4 pair"
  )],
  units[unit_kind == "single_trait_signal", .(
    signal_unit_id, source_signal_id, membership_role = "unmatched single-trait signal"
  )]
))
signal_unit_members <- merge(
  signal_unit_members,
  signals[, .(source_signal_id, trait_id, region, cs, lead_variant, cs_log10bf, low_purity, cs_size)],
  by = "source_signal_id", all.x = TRUE
)
unit_member_cc <- merge(
  member_cc,
  signal_unit_members[, .(signal_unit_id, source_signal_id)],
  by = "source_signal_id", allow.cartesian = TRUE
)
unit_cc_support <- unit_member_cc[, .(
  unit_ccgwas_robust_member_n = data.table::uniqueN(SNP),
  unit_ccgwas_source_signal_n = data.table::uniqueN(source_signal_id),
  unit_ccgwas_robust_member_pip = sum(cs_specific_pip_recomputed, na.rm = TRUE),
  unit_ccgwas_direction_consistency = {
    denom <- sum(cs_specific_pip_recomputed, na.rm = TRUE)
    if (denom > 0) sum(cs_specific_pip_recomputed[ccgwas_direction_consistent %in% TRUE], na.rm = TRUE) / denom else NA_real_
  },
  unit_ccgwas_locus_ids = paste(sort(unique(stats::na.omit(ccgwas_locus_id))), collapse = ";")
), by = signal_unit_id]
units <- merge(units, unit_cc_support, by = "signal_unit_id", all.x = TRUE)

# Annotate every unit representative variant in both GWAS using exact REF/ALT.
all_variants <- unique(units$representative_variant)
g_ann <- annotate_variants(all_variants, graves)
a_ann <- annotate_variants(all_variants, aiht)
data.table::setnames(
  g_ann,
  setdiff(names(g_ann), "variant_id"),
  paste0("graves_", setdiff(names(g_ann), "variant_id"))
)
data.table::setnames(
  a_ann,
  setdiff(names(a_ann), "variant_id"),
  paste0("aiht_", setdiff(names(a_ann), "variant_id"))
)
units <- merge(units, g_ann, by.x = "representative_variant", by.y = "variant_id", all.x = TRUE)
units <- merge(units, a_ann, by.x = "representative_variant", by.y = "variant_id", all.x = TRUE)

# Attach source-signal and pair evidence to the unmatched units.
source_meta <- signals[, .(
  source_signal_id, source_trait_id = trait_id, source_region = region, source_cs = cs,
  source_cs_log10bf = cs_log10bf, source_low_purity = low_purity,
  source_cs_size = cs_size, source_lead_signal_pip = lead_signal_pip,
  source_lead_marginal_beta = lead_marginal_beta,
  source_lead_marginal_se = lead_marginal_se, source_lead_marginal_p = lead_marginal_p,
  source_maf = NA_real_, source_is_chr_x = is_chr_x
)]
lead_maf <- members[order(-cs_specific_pip_recomputed), .SD[1L], by = source_signal_id][, .(
  source_signal_id, source_maf_member = maf
)]
source_meta <- merge(source_meta, lead_maf, by = "source_signal_id", all.x = TRUE)
source_meta[, source_maf := source_maf_member][, source_maf_member := NULL]
units <- merge(units, source_meta, by = "source_signal_id", all.x = TRUE)
units <- merge(units, pair_evidence, by = "source_signal_id", all.x = TRUE)
units <- merge(units, cc_support, by = "source_signal_id", all.x = TRUE)

# Restore H4 evidence fields after the unmatched-only joins.
units[unit_kind == "single_trait_signal", `:=`(
  evidence_pair_id = best_h3_pair_id,
  evidence_counterpart_signal_id = best_h3_counterpart,
  posterior_primary = best_h3_primary,
  posterior_sensitivity = best_h3_sensitivity,
  classification_posterior = best_h3_conservative
)]
units[, source_signal_quality := data.table::fifelse(
  unit_kind == "strong_H4_pair", TRUE,
  !is.na(source_low_purity) & !source_low_purity & source_cs_log10bf > as.numeric(cfg$t06_gwas_gwas_coloc$min_cs_log10bf)
)]
units[, distinct_model_supported := unit_kind == "single_trait_signal" &
  best_h3_quality %in% TRUE & classification_posterior >= posterior_threshold]
units[, shared_model_not_excluded := unit_kind == "single_trait_signal" & suggestive_shared_model %in% TRUE]

# For focal single-trait units, the catalog values are authoritative for the
# focal GWAS. The exact Tabix lookup supplies the opposite-trait estimate.
units[unit_kind == "single_trait_signal" & source_trait_id == graves, `:=`(
  graves_beta = source_lead_marginal_beta,
  graves_se = source_lead_marginal_se,
  graves_p = source_lead_marginal_p,
  aiht_beta = aiht_beta,
  aiht_se = aiht_se,
  aiht_p = aiht_p
)]
units[unit_kind == "single_trait_signal" & source_trait_id == aiht, `:=`(
  graves_beta = graves_beta,
  graves_se = graves_se,
  graves_p = graves_p,
  aiht_beta = source_lead_marginal_beta,
  aiht_se = source_lead_marginal_se,
  aiht_p = source_lead_marginal_p
)]
units[unit_kind == "strong_H4_pair", `:=`(
  graves_beta = graves_beta, graves_se = graves_se, graves_p = graves_p,
  aiht_beta = aiht_beta, aiht_se = aiht_se, aiht_p = aiht_p
)]
units[, `:=`(
  representative_nearest_genes = data.table::fcoalesce(graves_nearest_genes, aiht_nearest_genes),
  representative_reported_rsids = data.table::fcoalesce(graves_reported_rsids, aiht_reported_rsids),
  difference_beta_rho0 = graves_beta - aiht_beta,
  difference_se_rho0 = sqrt(graves_se^2 + aiht_se^2)
)]
units[, difference_p_rho0 := 2 * stats::pnorm(-abs(difference_beta_rho0 / difference_se_rho0))]
units[, case_difference_gws := is.finite(difference_p_rho0) & difference_p_rho0 < case_difference_p]
units[, case_difference_direction_consistent := data.table::fifelse(
  unit_kind != "single_trait_signal", NA,
  data.table::fifelse(
    source_trait_id == graves,
    difference_beta_rho0 * source_lead_marginal_beta > 0,
    difference_beta_rho0 * source_lead_marginal_beta < 0
  )
)]
units[, `:=`(
  other_trait = data.table::fifelse(source_trait_id == graves, "AIHT", "Graves"),
  other_beta = data.table::fifelse(source_trait_id == graves, aiht_beta, graves_beta),
  other_se = data.table::fifelse(source_trait_id == graves, aiht_se, graves_se)
)]
units[, `:=`(
  other_ci90_lower = other_beta - z_equivalence * other_se,
  other_ci90_upper = other_beta + z_equivalence * other_se
)]
for (bound in or_bounds) {
  suffix <- gsub("\\.", "_", format(bound, nsmall = 2L))
  log_bound <- log(bound)
  units[, (paste0("other_equivalent_or", suffix)) :=
    is.finite(other_ci90_lower) & is.finite(other_ci90_upper) &
      other_ci90_lower > -log_bound & other_ci90_upper < log_bound]
}

classify_at_bound <- function(dt, bound) {
  suffix <- gsub("\\.", "_", format(bound, nsmall = 2L))
  equivalent <- dt[[paste0("other_equivalent_or", suffix)]] %in% TRUE
  predominant_ready <- dt$unit_kind == "single_trait_signal" &
    dt$source_signal_quality %in% TRUE &
    dt$distinct_model_supported %in% TRUE &
    !(dt$shared_model_not_excluded %in% TRUE) &
    equivalent &
    dt$case_difference_gws %in% TRUE &
    dt$case_difference_direction_consistent %in% TRUE
  out <- rep("unresolved", nrow(dt))
  out[dt$unit_kind == "strong_H4_pair" & dt$preliminary_shared_direction == "concordant"] <- "shared-concordant"
  out[dt$unit_kind == "strong_H4_pair" & dt$preliminary_shared_direction == "opposite"] <- "shared-opposite"
  out[predominant_ready & dt$source_trait_id == graves] <- "Graves-predominant"
  out[predominant_ready & dt$source_trait_id == aiht] <- "AIHT-predominant"
  out
}

sensitivity <- data.table::rbindlist(lapply(or_bounds, function(bound) {
  data.table::data.table(
    signal_unit_id = units$signal_unit_id,
    equivalence_or_upper = bound,
    primary_threshold = bound == primary_or_upper,
    five_class = classify_at_bound(units, bound)
  )
}))
primary_class <- sensitivity[primary_threshold == TRUE, .(signal_unit_id, five_class)]
data.table::setnames(primary_class, "five_class", "five_class_primary")
units <- merge(units, primary_class, by = "signal_unit_id", all.x = TRUE)
class_robust <- sensitivity[, .(
  classification_across_bounds = paste(five_class[order(equivalence_or_upper)], collapse = ";"),
  classification_robust_all_bounds = data.table::uniqueN(five_class) == 1L
), by = signal_unit_id]
units <- merge(units, class_robust, by = "signal_unit_id", all.x = TRUE)

unresolved_reason <- function(i) {
  if (units$five_class_primary[[i]] != "unresolved") return(NA_character_)
  reasons <- character()
  if (units$unit_kind[[i]] == "strong_H4_pair") return("shared_direction_undetermined")
  if (!(units$source_signal_quality[[i]] %in% TRUE)) reasons <- c(reasons, "weak_or_low_purity_source_signal")
  if (is.na(units$n_direct_signal_pairs[[i]]) || units$n_direct_signal_pairs[[i]] == 0L) {
    reasons <- c(reasons, "no_two_trait_finemap_comparison")
  }
  if (units$shared_model_not_excluded[[i]] %in% TRUE) reasons <- c(reasons, "suggestive_H4_not_excluded")
  if (!(units$distinct_model_supported[[i]] %in% TRUE)) reasons <- c(reasons, "H3_posterior_below_0.80_or_QC_failed")
  eq_col <- paste0("other_equivalent_or", gsub("\\.", "_", format(primary_or_upper, nsmall = 2L)))
  if (!(units[[eq_col]][[i]] %in% TRUE)) reasons <- c(reasons, "other_trait_not_equivalent_at_OR1.10")
  if (!(units$case_difference_gws[[i]] %in% TRUE)) reasons <- c(reasons, "case_difference_not_GWS")
  if (!(units$case_difference_direction_consistent[[i]] %in% TRUE)) reasons <- c(reasons, "case_difference_direction_inconsistent_or_missing")
  if (!is.finite(units$other_beta[[i]]) || !is.finite(units$other_se[[i]])) reasons <- c(reasons, "opposite_trait_effect_unavailable")
  paste(unique(reasons), collapse = ";")
}
units[, unresolved_reason := vapply(seq_len(.N), unresolved_reason, character(1))]

# Signal-unit support from T04B is reported separately from the mandatory gates.
units[, source_t04b_formal_cs_support := ccgwas_robust_member_n > 0 & ccgwas_direction_consistency >= posterior_threshold]
units[is.na(source_t04b_formal_cs_support), source_t04b_formal_cs_support := FALSE]
units[, t04b_formal_cs_support := unit_ccgwas_robust_member_n > 0 & unit_ccgwas_direction_consistency >= posterior_threshold]
units[is.na(t04b_formal_cs_support), t04b_formal_cs_support := FALSE]

# Map T04B loci to T07 signal units via credible-set members, rather than broad
# three-megabase region overlap.
member_units <- merge(
  members[, .(source_signal_id, chromosome, position)],
  signal_unit_members[, .(signal_unit_id, source_signal_id)],
  by = "source_signal_id", allow.cartesian = TRUE
)
locus_unit_rows <- list()
for (i in seq_len(nrow(cc_loci))) {
  loc <- cc_loci[i]
  hits <- unique(member_units[
    chromosome == paste0("chr", loc$CHR) & position >= loc$start & position <= loc$end,
    signal_unit_id
  ])
  if (length(hits)) {
    locus_unit_rows[[length(locus_unit_rows) + 1L]] <- data.table::data.table(
      locus_id = loc$locus_id, signal_unit_id = hits
    )
  }
}
locus_units <- data.table::rbindlist(locus_unit_rows, fill = TRUE)
if (nrow(locus_units)) {
  locus_units <- merge(
    locus_units,
    units[, .(signal_unit_id, five_class_primary, representative_variant, representative_nearest_genes,
      classification_posterior, classification_robust_all_bounds)],
    by = "signal_unit_id", all.x = TRUE
  )
}
t04b_resolution <- merge(
  cc_loci[, .(locus_id, lead_SNP, CHR, start, end, observed_case_control_sign_pattern)],
  if (nrow(locus_units)) locus_units[, {
    finite_posteriors <- classification_posterior[is.finite(classification_posterior)]
    .(
    t07_signal_unit_n = data.table::uniqueN(signal_unit_id),
    t07_signal_units = paste(sort(unique(signal_unit_id)), collapse = ";"),
    t07_classes = paste(sort(unique(five_class_primary)), collapse = ";"),
    t07_non_unresolved_n = sum(five_class_primary != "unresolved"),
    max_classification_posterior = if (length(finite_posteriors)) max(finite_posteriors) else NA_real_,
    any_threshold_robust = any(classification_robust_all_bounds & five_class_primary != "unresolved")
    )
  }, by = locus_id] else data.table::data.table(locus_id = character()),
  by = "locus_id", all.x = TRUE
)
t04b_resolution[!is.finite(max_classification_posterior), max_classification_posterior := NA_real_]
t04b_resolution[is.na(t07_signal_unit_n), `:=`(
  t07_signal_unit_n = 0L, t07_signal_units = "", t07_classes = "unresolved",
  t07_non_unresolved_n = 0L, any_threshold_robust = FALSE
)]

# Keep the main result readable while retaining every audit field.
preferred_order <- c(
  "signal_unit_id", "five_class_primary", "classification_robust_all_bounds",
  "classification_across_bounds", "unit_kind", "chromosome", "position",
  "representative_variant", "representative_reported_rsids", "representative_nearest_genes",
  "graves_source_signal_id", "aiht_source_signal_id", "source_signal_id", "source_trait_id",
  "evidence_pair_id", "evidence_counterpart_signal_id", "posterior_type",
  "posterior_primary", "posterior_sensitivity", "classification_posterior",
  "graves_beta", "graves_se", "graves_p", "aiht_beta", "aiht_se", "aiht_p",
  "difference_beta_rho0", "difference_se_rho0", "difference_p_rho0",
  "preliminary_shared_direction", "other_trait", "other_beta", "other_se",
  "other_ci90_lower", "other_ci90_upper",
  paste0("other_equivalent_or", gsub("\\.", "_", format(or_bounds, nsmall = 2L))),
  "source_signal_quality", "distinct_model_supported", "shared_model_not_excluded",
  "case_difference_gws", "case_difference_direction_consistent",
  "ccgwas_robust_member_n", "ccgwas_robust_member_pip", "ccgwas_direction_consistency",
  "ccgwas_locus_ids", "source_t04b_formal_cs_support",
  "unit_ccgwas_robust_member_n", "unit_ccgwas_robust_member_pip",
  "unit_ccgwas_direction_consistency", "unit_ccgwas_locus_ids",
  "t04b_formal_cs_support", "unresolved_reason"
)
data.table::setcolorder(units, c(intersect(preferred_order, names(units)), setdiff(names(units), preferred_order)))
data.table::setorder(units, signal_unit_id)

rules <- data.table::data.table(
  classification = c("shared-concordant", "shared-opposite", "Graves-predominant", "AIHT-predominant", "unresolved"),
  rule = c(
    "Strong T06 H4 under both p12 priors; both signals pass QC; effects at top H4 variant have the same ALT-oriented sign",
    "Strong T06 H4 under both p12 priors; both signals pass QC; effects at top H4 variant have opposite ALT-oriented signs",
    "Unmatched Graves SuSiE signal; source QC pass; no suggestive H4; conservative H3 posterior >=0.80; AIHT 90% CI equivalent to zero; rho=0 case difference P<5e-8 and direction consistent",
    "Unmatched AIHT SuSiE signal; source QC pass; no suggestive H4; conservative H3 posterior >=0.80; Graves 90% CI equivalent to zero; rho=0 case difference P<5e-8 and direction consistent",
    "Any mandatory criterion missing or failed; includes one-trait-only fine-mapping, weak/low-purity signals, prior-sensitive sharing, non-equivalence, or inadequate case-difference evidence"
  ),
  primary_equivalence_or_upper = c(NA, NA, primary_or_upper, primary_or_upper, NA),
  classification_posterior_threshold = posterior_threshold,
  therapeutic_interpretation = c(
    "Eligible for later common-target evaluation only after molecular direction and external replication",
    "Mechanistic-switch candidate; not a common therapeutic target",
    "Disease-predominant candidate; external replication required",
    "Disease-predominant candidate; external replication required",
    "No directional therapeutic interpretation"
  )
)

summary_table <- data.table::rbindlist(list(
  data.table::data.table(metric = "source_susie_signals", value = nrow(signals)),
  data.table::data.table(metric = "strong_H4_pairs_collapsed", value = nrow(shared_pairs)),
  data.table::data.table(metric = "independent_signal_units", value = nrow(units)),
  units[, .(value = .N), by = .(metric = paste0("primary_class_", five_class_primary))],
  data.table::data.table(metric = "classification_robust_all_equivalence_bounds", value = units[classification_robust_all_bounds == TRUE, .N]),
  data.table::data.table(metric = "predominant_primary", value = units[five_class_primary %in% c("Graves-predominant", "AIHT-predominant"), .N]),
  data.table::data.table(metric = "predominant_robust_all_bounds", value = units[
    five_class_primary %in% c("Graves-predominant", "AIHT-predominant") & classification_robust_all_bounds == TRUE, .N
  ]),
  data.table::data.table(metric = "t04b_loci_with_non_unresolved_t07_signal", value = t04b_resolution[t07_non_unresolved_n > 0, .N])
), fill = TRUE)

provenance <- data.table::data.table(
  field = c(
    "genome_build", "source_release", "shared_model", "classification_posterior",
    "primary_equivalence", "equivalence_sensitivity", "case_difference_model",
    "HLA_scope", "sample_overlap_limit", "external_validation"
  ),
  value = c(
    cfg$project$genome_build, "FinnGen R13", "T06 direct SuSiE H4/H3 with two p12 priors",
    paste0("minimum of primary and sensitivity H4/H3; threshold ", posterior_threshold),
    paste0("90% CI within +/-log(", primary_or_upper, ")"),
    paste(or_bounds, collapse = ","),
    paste0("beta_Graves-beta_AIHT; rho=0; P<", case_difference_p),
    "non-HLA only; no public official HLA SuSiE signal entered direct coloc",
    "same FinnGen release; correlated summary-statistic error is not explicitly modeled by direct coloc",
    "required before mechanistic or therapeutic claims"
  )
)

# Assertions encode the T07 completion gate.
stopifnot(
  nrow(signals) == 312L,
  nrow(shared_pairs) == 8L,
  nrow(units) == nrow(signals) - nrow(shared_pairs),
  nrow(signal_unit_members) == nrow(signals),
  !anyDuplicated(signal_unit_members$source_signal_id),
  all(units$five_class_primary %in% c(
    "shared-concordant", "shared-opposite", "Graves-predominant", "AIHT-predominant", "unresolved"
  )),
  units[five_class_primary == "shared-concordant", .N] == 6L,
  units[five_class_primary == "shared-opposite", .N] == 2L,
  all(units[five_class_primary %in% c("shared-concordant", "shared-opposite"), classification_posterior >= posterior_threshold]),
  all(units[five_class_primary %in% c("Graves-predominant", "AIHT-predominant"),
    source_signal_quality & distinct_model_supported & !shared_model_not_excluded &
      case_difference_gws & case_difference_direction_consistent]),
  all(units[five_class_primary == "Graves-predominant", source_trait_id == graves]),
  all(units[five_class_primary == "AIHT-predominant", source_trait_id == aiht])
)

write_tsv(units, file.path(out_dir, "genetic_signal_five_classification.tsv"))
write_tsv(signal_unit_members, file.path(out_dir, "signal_unit_members.tsv"))
write_tsv(sensitivity, file.path(out_dir, "classification_equivalence_sensitivity.tsv"))
write_tsv(rules, file.path(out_dir, "classification_rules.tsv"))
write_tsv(member_cc, file.path(out_dir, "credible_set_ccgwas_support.tsv.gz"))
write_tsv(locus_units, file.path(out_dir, "t04b_locus_signal_units.tsv"))
write_tsv(t04b_resolution, file.path(out_dir, "t04b_locus_t07_resolution.tsv"))
write_tsv(summary_table, file.path(out_dir, "t07_summary.tsv"))
write_tsv(provenance, file.path(out_dir, "t07_provenance.tsv"))

log_message(
  "T07 classification complete: ", nrow(units), " units; ",
  units[five_class_primary == "shared-concordant", .N], " shared-concordant; ",
  units[five_class_primary == "shared-opposite", .N], " shared-opposite; ",
  units[five_class_primary == "Graves-predominant", .N], " Graves-predominant; ",
  units[five_class_primary == "AIHT-predominant", .N], " AIHT-predominant; ",
  units[five_class_primary == "unresolved", .N], " unresolved"
)
