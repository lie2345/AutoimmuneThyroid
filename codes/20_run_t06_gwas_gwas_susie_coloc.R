#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c(
  "yaml", "data.table", "coloc", "Rsamtools", "GenomicRanges", "IRanges", "ggplot2"
))
cfg <- read_analysis_config(project_root)
if (is.null(cfg$t06_gwas_gwas_coloc)) stop("Missing t06_gwas_gwas_coloc config", call. = FALSE)
t06_cfg <- cfg$t06_gwas_gwas_coloc
graves <- cfg$primary_comparison$graves
aiht <- cfg$primary_comparison$autoimmune_hypothyroidism
out_dir <- file.path(project_root, cfg$project$outcome_dir, "t06")
figure_dir <- file.path(project_root, cfg$project$figure_dir, "t06")
ensure_dir(out_dir)
ensure_dir(figure_dir)

susie_dir <- file.path(project_root, "raw_data", "finngen_r13", "finemap", "susie")
finemap_dir <- file.path(project_root, "raw_data", "finngen_r13", "finemap", "finemap")

susie_file <- function(trait, suffix) {
  file.path(susie_dir, paste0("finngen_R13_", trait, ".SUSIE.", suffix, ".bgz"))
}
finemap_file <- function(trait, suffix) {
  file.path(finemap_dir, paste0("finngen_R13_", trait, ".FINEMAP.", suffix, ".bgz"))
}

required <- c(
  susie_file(graves, "snp"), susie_file(graves, "cred"), susie_file(graves, "cred_99"),
  susie_file(aiht, "snp"), susie_file(aiht, "cred"), susie_file(aiht, "cred_99"),
  finemap_file(graves, "region"), finemap_file(aiht, "region")
)
if (any(!file.exists(required))) {
  stop("Missing T06 official fine-mapping input(s): ", paste(required[!file.exists(required)], collapse = ", "), call. = FALSE)
}
for (path in c(susie_file(graves, "snp"), susie_file(aiht, "snp"))) {
  if (!file.exists(paste0(path, ".tbi"))) stop("Missing tabix index: ", path, ".tbi", call. = FALSE)
}

parse_regions <- function(region) {
  pieces <- regmatches(region, regexec("^chr([^:]+):([0-9]+)-([0-9]+)$", region))
  if (any(lengths(pieces) != 4L)) stop("Unparseable fine-mapping region", call. = FALSE)
  data.table::data.table(
    chrom_label = vapply(pieces, `[[`, character(1), 2L),
    region_start = as.integer(vapply(pieces, `[[`, character(1), 3L)),
    region_end = as.integer(vapply(pieces, `[[`, character(1), 4L))
  )
}

read_cred <- function(trait, suffix = "cred") {
  x <- data.table::fread(susie_file(trait, suffix), showProgress = FALSE)
  x <- cbind(x, parse_regions(x$region))
  x[, `:=`(
    trait_id = trait,
    credible_level = if (suffix == "cred") 0.95 else 0.99,
    is_hla = chrom_label == as.character(cfg$hla$chromosome) &
      region_start <= as.integer(cfg$hla$end) & region_end >= as.integer(cfg$hla$start),
    is_chr_x = chrom_label == "X"
  )]
  x
}

cred95_g <- read_cred(graves, "cred")
cred95_a <- read_cred(aiht, "cred")
cred99_g <- read_cred(graves, "cred_99")
cred99_a <- read_cred(aiht, "cred_99")
cred95 <- data.table::rbindlist(list(cred95_g, cred95_a), fill = TRUE)
cred99 <- data.table::rbindlist(list(cred99_g, cred99_a), fill = TRUE)

region_meta <- function(x, prefix) {
  y <- unique(x[, .(region, chrom_label, region_start, region_end, is_hla, is_chr_x)])
  data.table::setnames(
    y,
    c("region", "region_start", "region_end", "is_hla", "is_chr_x"),
    paste0(prefix, c("region", "start", "end", "is_hla", "is_chr_x"))
  )
  y
}
g_regions <- region_meta(cred95_g, "g_")
a_regions <- region_meta(cred95_a, "a_")
region_pairs <- merge(g_regions, a_regions, by = "chrom_label", allow.cartesian = TRUE)
region_pairs <- region_pairs[g_start <= a_end & g_end >= a_start]
region_pairs[, `:=`(
  pair_id = paste0("RP", sprintf("%02d", seq_len(.N))),
  overlap_start = pmax(g_start, a_start),
  overlap_end = pmin(g_end, a_end),
  is_hla = g_is_hla | a_is_hla,
  is_chr_x = g_is_chr_x | a_is_chr_x
)]
data.table::setcolorder(region_pairs, c(
  "pair_id", "chrom_label", "g_region", "a_region", "g_start", "g_end",
  "a_start", "a_end", "overlap_start", "overlap_end", "is_hla", "is_chr_x"
))
write_tsv(region_pairs, file.path(out_dir, "overlapping_finemap_region_pairs.tsv"))
write_tsv(region_pairs[is_hla == TRUE], file.path(out_dir, "hla_region_pairs_excluded.tsv"))
ordinary_pairs <- region_pairs[is_hla == FALSE]

read_finemap_region <- function(trait) {
  x <- data.table::fread(finemap_file(trait, "region"), showProgress = FALSE)
  x[, trait_id := trait]
  x
}
finemap_region <- data.table::rbindlist(list(
  read_finemap_region(graves), read_finemap_region(aiht)
), fill = TRUE)
write_tsv(finemap_region, file.path(out_dir, "finemap_region_support.tsv"))

lbf_to_alpha <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(rep(NA_real_, length(x)))
  anchor <- max(x[is.finite(x)])
  w <- exp(x - anchor)
  w[!is.finite(w)] <- 0
  total <- sum(w)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(x)))
  w / total
}

probmass_by_position <- function(x, lbf_col, start, end) {
  pip <- lbf_to_alpha(x[[lbf_col]])
  sum(pip[x$position >= start & x$position <= end], na.rm = TRUE)
}

minimal_snp_columns <- c(
  "trait", "region", "v", "rsid", "chromosome", "position", "allele1", "allele2",
  "maf", "beta", "se", "p", "mean", "sd", "prob", "cs", "cs_specific_prob",
  "low_purity", "lead_r2", paste0("alpha", 1:10), paste0("lbf_variable", 1:10)
)

process_trait <- function(trait, cred, cache_regions) {
  trait_value <- trait
  path <- susie_file(trait, "snp")
  header <- read_gzip_header(path)
  missing_columns <- setdiff(minimal_snp_columns, header)
  if (length(missing_columns)) stop("Missing SuSiE SNP columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  tabix_file <- Rsamtools::TabixFile(path)
  available <- suppressWarnings(Rsamtools::seqnamesTabix(tabix_file))
  suppressWarnings(Rsamtools::open.TabixFile(tabix_file))
  on.exit(Rsamtools::close.TabixFile(tabix_file), add = TRUE)

  catalog <- list()
  members <- list()
  region_qc <- list()
  cache <- list()
  regions <- unique(cred[, .(region, chrom_label, region_start, region_end, is_hla, is_chr_x)])
  for (i in seq_len(nrow(regions))) {
    meta <- regions[i]
    seqname <- paste0("chr", meta$chrom_label)
    if (!seqname %in% available) stop("Chromosome absent from SuSiE tabix index: ", seqname, call. = FALSE)
    param <- GenomicRanges::GRanges(
      seqnames = seqname,
      ranges = IRanges::IRanges(start = meta$region_start, end = meta$region_end)
    )
    lines <- suppressWarnings(Rsamtools::scanTabix(tabix_file, param = param)[[1]])
    if (!length(lines)) stop("No SuSiE SNP rows for ", trait, " ", meta$region, call. = FALSE)
    x <- data.table::fread(
      text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE,
      select = match(minimal_snp_columns, header), showProgress = FALSE
    )
    data.table::setnames(x, minimal_snp_columns)
    x <- x[region == meta$region & trait == trait_value]
    if (!nrow(x)) stop("Tabix rows did not match exact region: ", meta$region, call. = FALSE)
    if (anyDuplicated(x$rsid)) stop("Duplicate variant IDs within region: ", meta$region, call. = FALSE)

    cs_rows <- cred[region == meta$region]
    for (j in seq_len(nrow(cs_rows))) {
      cs_meta <- cs_rows[j]
      signal <- as.integer(cs_meta$cs)
      lbf_col <- paste0("lbf_variable", signal)
      alpha_col <- paste0("alpha", signal)
      if (!lbf_col %in% names(x)) stop("Missing signal LBF column ", lbf_col, call. = FALSE)
      signal_pip <- lbf_to_alpha(x[[lbf_col]])
      lead_index <- which.max(x[[lbf_col]])
      cs_members <- x[cs == signal]
      if (nrow(cs_members) != as.integer(cs_meta$cs_size)) {
        stop(
          "Credible-set size mismatch for ", trait, " ", meta$region, " CS", signal,
          ": expected ", cs_meta$cs_size, ", observed ", nrow(cs_members), call. = FALSE
        )
      }
      member_indices <- which(x$cs == signal)
      members[[length(members) + 1L]] <- data.table::data.table(
        trait_id = trait,
        region = meta$region,
        cs = signal,
        rsid = x$rsid[member_indices],
        v = x$v[member_indices],
        chromosome = x$chromosome[member_indices],
        position = x$position[member_indices],
        allele1 = x$allele1[member_indices],
        allele2 = x$allele2[member_indices],
        maf = x$maf[member_indices],
        marginal_beta = x$beta[member_indices],
        marginal_se = x$se[member_indices],
        marginal_p = x$p[member_indices],
        overall_pip = x$prob[member_indices],
        cs_specific_pip_official = x$cs_specific_prob[member_indices],
        cs_specific_pip_recomputed = signal_pip[member_indices],
        alpha_official = x[[alpha_col]][member_indices],
        lead_r2 = x$lead_r2[member_indices],
        cs_log10bf = cs_meta$cs_log10bf,
        cs_avg_r2 = cs_meta$cs_avg_r2,
        cs_min_r2 = cs_meta$cs_min_r2,
        low_purity = as.logical(cs_meta$low_purity),
        is_hla = meta$is_hla,
        is_chr_x = meta$is_chr_x
      )
      fm <- finemap_region[trait_id == trait & region == meta$region]
      catalog[[length(catalog) + 1L]] <- data.table::data.table(
        trait_id = trait,
        region = meta$region,
        chrom_label = meta$chrom_label,
        region_start = meta$region_start,
        region_end = meta$region_end,
        cs = signal,
        cs_log10bf = cs_meta$cs_log10bf,
        cs_avg_r2 = cs_meta$cs_avg_r2,
        cs_min_r2 = cs_meta$cs_min_r2,
        low_purity = as.logical(cs_meta$low_purity),
        cs_size = as.integer(cs_meta$cs_size),
        lead_variant = x$rsid[lead_index],
        lead_position = x$position[lead_index],
        lead_signal_pip = signal_pip[lead_index],
        lead_overall_pip = x$prob[lead_index],
        lead_marginal_beta = x$beta[lead_index],
        lead_marginal_se = x$se[lead_index],
        lead_marginal_p = x$p[lead_index],
        finemap_region_log10bf = if (nrow(fm)) fm$log10bf[[1L]] else NA_real_,
        finemap_expected_signals = if (nrow(fm)) fm$expectedvalue[[1L]] else NA_real_,
        convergence_status = "not_exposed_in_public_summary",
        ld_source = "FinnGen in-sample dosage LD (official pipeline)",
        genome_build = cfg$project$genome_build,
        is_hla = meta$is_hla,
        is_chr_x = meta$is_chr_x
      )
    }
    region_qc[[length(region_qc) + 1L]] <- data.table::data.table(
      trait_id = trait,
      region = meta$region,
      n_variants = nrow(x),
      n_credible_sets = nrow(cs_rows),
      n_pure_cs = sum(!as.logical(cs_rows$low_purity)),
      n_strong_pure_cs = sum(
        !as.logical(cs_rows$low_purity) & cs_rows$cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf)
      ),
      is_hla = meta$is_hla,
      is_chr_x = meta$is_chr_x
    )
    if (meta$region %in% cache_regions) cache[[meta$region]] <- x
    if (i == 1L || i %% 25L == 0L || i == nrow(regions)) {
      log_message("T06 audited ", trait, " region ", i, "/", nrow(regions))
    }
  }
  list(
    catalog = data.table::rbindlist(catalog, fill = TRUE),
    members = data.table::rbindlist(members, fill = TRUE),
    region_qc = data.table::rbindlist(region_qc, fill = TRUE),
    cache = cache
  )
}

g_processed <- process_trait(graves, cred95_g, unique(ordinary_pairs$g_region))
a_processed <- process_trait(aiht, cred95_a, unique(ordinary_pairs$a_region))
signal_catalog <- data.table::rbindlist(list(g_processed$catalog, a_processed$catalog), fill = TRUE)
credible_members <- data.table::rbindlist(list(g_processed$members, a_processed$members), fill = TRUE)
region_qc <- data.table::rbindlist(list(g_processed$region_qc, a_processed$region_qc), fill = TRUE)
write_tsv(signal_catalog, file.path(out_dir, "susie_signal_catalog.tsv"))
write_tsv(credible_members, file.path(out_dir, "susie_credible_set_members.tsv.gz"))
write_tsv(region_qc, file.path(out_dir, "susie_region_qc.tsv"))
write_tsv(cred99, file.path(out_dir, "susie_credible_sets_99_summary.tsv"))

coloc_bf_bf <- getFromNamespace("coloc.bf_bf", "coloc")
coloc_rows <- list()
h4_rows <- list()
pair_qc <- list()

for (pair_index in seq_len(nrow(ordinary_pairs))) {
  pair <- ordinary_pairs[pair_index]
  g <- data.table::copy(g_processed$cache[[pair$g_region]])
  a <- data.table::copy(a_processed$cache[[pair$a_region]])
  if (is.null(g) || is.null(a)) stop("Missing cached fine-mapping region", call. = FALSE)
  g_signals <- sort(cred95_g[region == pair$g_region, as.integer(cs)])
  a_signals <- sort(cred95_a[region == pair$a_region, as.integer(cs)])

  g_keep <- c(
    "rsid", "position", "allele1", "allele2", "beta", "se", "p", "maf", "cs",
    paste0("lbf_variable", g_signals)
  )
  a_keep <- c(
    "rsid", "position", "allele1", "allele2", "beta", "se", "p", "maf", "cs",
    paste0("lbf_variable", a_signals)
  )
  g_small <- g[, ..g_keep]
  a_small <- a[, ..a_keep]
  data.table::setnames(g_small, setdiff(names(g_small), "rsid"), paste0(setdiff(names(g_small), "rsid"), "_g"))
  data.table::setnames(a_small, setdiff(names(a_small), "rsid"), paste0(setdiff(names(a_small), "rsid"), "_a"))
  common <- merge(g_small, a_small, by = "rsid")
  if (!nrow(common)) stop("No shared variants for fine-mapping region pair ", pair$pair_id, call. = FALSE)
  allele_mismatch <- sum(
    common$allele1_g != common$allele1_a | common$allele2_g != common$allele2_a |
      common$position_g != common$position_a,
    na.rm = TRUE
  )
  if (allele_mismatch > 0L) stop("Allele mismatch in region pair ", pair$pair_id, call. = FALSE)

  g_lbf_cols <- paste0("lbf_variable", g_signals, "_g")
  a_lbf_cols <- paste0("lbf_variable", a_signals, "_a")
  bf_g <- t(as.matrix(common[, ..g_lbf_cols]))
  bf_a <- t(as.matrix(common[, ..a_lbf_cols]))
  colnames(bf_g) <- common$rsid
  colnames(bf_a) <- common$rsid
  primary <- coloc_bf_bf(
    bf_g, bf_a,
    p1 = as.numeric(cfg$coloc$p1), p2 = as.numeric(cfg$coloc$p2),
    p12 = as.numeric(t06_cfg$primary_p12), trim_by_posterior = FALSE
  )
  sensitivity <- coloc_bf_bf(
    bf_g, bf_a,
    p1 = as.numeric(cfg$coloc$p1), p2 = as.numeric(cfg$coloc$p2),
    p12 = as.numeric(t06_cfg$sensitivity_p12), trim_by_posterior = FALSE
  )
  primary_summary <- data.table::as.data.table(primary$summary)
  sensitivity_summary <- data.table::as.data.table(sensitivity$summary)[, .(
    idx1, idx2,
    sensitivity_PP.H0 = PP.H0.abf,
    sensitivity_PP.H1 = PP.H1.abf,
    sensitivity_PP.H2 = PP.H2.abf,
    sensitivity_PP.H3 = PP.H3.abf,
    sensitivity_PP.H4 = PP.H4.abf
  )]
  result <- merge(primary_summary, sensitivity_summary, by = c("idx1", "idx2"), sort = FALSE)
  result[, `:=`(g_cs = g_signals[idx1], a_cs = a_signals[idx2])]

  common_min <- min(common$position_g)
  common_max <- max(common$position_g)
  for (k in seq_len(nrow(result))) {
    g_cs <- result$g_cs[[k]]
    a_cs <- result$a_cs[[k]]
    g_meta <- cred95_g[region == pair$g_region & cs == g_cs]
    a_meta <- cred95_a[region == pair$a_region & cs == a_cs]
    g_lbf <- paste0("lbf_variable", g_cs)
    a_lbf <- paste0("lbf_variable", a_cs)
    alpha_g <- lbf_to_alpha(g[[g_lbf]])
    alpha_a <- lbf_to_alpha(a[[a_lbf]])
    cs_common <- common[cs_g == g_cs & cs_a == a_cs]
    if (nrow(cs_common)) {
      idx_g <- match(cs_common$rsid, g$rsid)
      idx_a <- match(cs_common$rsid, a$rsid)
      clpp <- sum(alpha_g[idx_g] * alpha_a[idx_a], na.rm = TRUE)
      clpa <- sum(pmin(alpha_g[idx_g], alpha_a[idx_a]), na.rm = TRUE)
    } else {
      clpp <- NA_real_
      clpa <- NA_real_
    }
    probmass_g <- probmass_by_position(g, g_lbf, common_min, common_max)
    probmass_a <- probmass_by_position(a, a_lbf, common_min, common_max)
    lead_g_index <- which.max(g[[g_lbf]])
    lead_a_index <- which.max(a[[a_lbf]])
    top_in_overlap <- paste0(
      as.integer(g$position[lead_g_index] >= common_min & g$position[lead_g_index] <= common_max), ",",
      as.integer(a$position[lead_a_index] >= common_min & a$position[lead_a_index] <= common_max)
    )
    h4_col <- if (nrow(result) == 1L) "SNP.PP.H4.abf" else paste0("SNP.PP.H4.row", k)
    h4_table <- data.table::as.data.table(primary$results)[, .(
      rsid = snp,
      SNP.PP.H4 = get(h4_col)
    )]
    top_h4 <- h4_table[which.max(SNP.PP.H4)]
    top_common <- common[rsid == top_h4$rsid][1]
    direction <- if (!nrow(top_common) || !is.finite(top_common$beta_g * top_common$beta_a)) {
      "undetermined"
    } else if (top_common$beta_g * top_common$beta_a > 0) {
      "concordant"
    } else if (top_common$beta_g * top_common$beta_a < 0) {
      "opposite"
    } else {
      "undetermined"
    }
    quality_core <- !as.logical(g_meta$low_purity) & !as.logical(a_meta$low_purity) &
      g_meta$cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf) &
      a_meta$cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf) &
      probmass_g > as.numeric(t06_cfg$min_probmass_core) &
      probmass_a > as.numeric(t06_cfg$min_probmass_core) & nrow(cs_common) > 0L
    quality_suggestive <- !as.logical(g_meta$low_purity) & !as.logical(a_meta$low_purity) &
      g_meta$cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf) &
      a_meta$cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf) &
      probmass_g > as.numeric(t06_cfg$min_probmass_suggestive) &
      probmass_a > as.numeric(t06_cfg$min_probmass_suggestive) & nrow(cs_common) > 0L
    quality_h3 <- !as.logical(g_meta$low_purity) & !as.logical(a_meta$low_purity) &
      g_meta$cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf) &
      a_meta$cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf) &
      probmass_g > as.numeric(t06_cfg$min_probmass_core) &
      probmass_a > as.numeric(t06_cfg$min_probmass_core)
    classification <- if (as.logical(g_meta$low_purity) || as.logical(a_meta$low_purity)) {
      "unresolved_low_purity"
    } else if (
      g_meta$cs_log10bf <= as.numeric(t06_cfg$min_cs_log10bf) ||
        a_meta$cs_log10bf <= as.numeric(t06_cfg$min_cs_log10bf)
    ) {
      "unresolved_weak_credible_set"
    } else if (
      quality_core & result$PP.H4.abf[[k]] >= as.numeric(t06_cfg$strong_h4) &
        result$sensitivity_PP.H4[[k]] >= as.numeric(t06_cfg$strong_h4)
    ) {
      "shared_strong_H4"
    } else if (
      quality_suggestive & result$PP.H4.abf[[k]] >= as.numeric(t06_cfg$suggestive_h4) &
        result$sensitivity_PP.H4[[k]] >= as.numeric(t06_cfg$suggestive_h4)
    ) {
      "shared_suggestive_H4"
    } else if (
      quality_h3 & result$PP.H3.abf[[k]] >= as.numeric(t06_cfg$distinct_h3) &
        result$sensitivity_PP.H3[[k]] >= as.numeric(t06_cfg$distinct_h3)
    ) {
      "distinct_signals_H3"
    } else if (
      probmass_g <= as.numeric(t06_cfg$min_probmass_suggestive) ||
        probmass_a <= as.numeric(t06_cfg$min_probmass_suggestive) || nrow(cs_common) == 0L
    ) {
      "unresolved_insufficient_overlap"
    } else {
      "unresolved_posterior"
    }
    coloc_rows[[length(coloc_rows) + 1L]] <- data.table::data.table(
      pair_id = pair$pair_id,
      chromosome = pair$chrom_label,
      graves_region = pair$g_region,
      aiht_region = pair$a_region,
      graves_cs = g_cs,
      aiht_cs = a_cs,
      n_snps_graves = nrow(g),
      n_snps_aiht = nrow(a),
      n_shared_snps = nrow(common),
      graves_hit = result$hit1[[k]],
      aiht_hit = result$hit2[[k]],
      PP.H0 = result$PP.H0.abf[[k]],
      PP.H1 = result$PP.H1.abf[[k]],
      PP.H2 = result$PP.H2.abf[[k]],
      PP.H3 = result$PP.H3.abf[[k]],
      PP.H4 = result$PP.H4.abf[[k]],
      sensitivity_PP.H3 = result$sensitivity_PP.H3[[k]],
      sensitivity_PP.H4 = result$sensitivity_PP.H4[[k]],
      graves_cs_log10bf = g_meta$cs_log10bf,
      aiht_cs_log10bf = a_meta$cs_log10bf,
      graves_low_purity = as.logical(g_meta$low_purity),
      aiht_low_purity = as.logical(a_meta$low_purity),
      graves_cs_size = as.integer(g_meta$cs_size),
      aiht_cs_size = as.integer(a_meta$cs_size),
      cs_overlap = nrow(cs_common),
      graves_probmass_in_overlap = probmass_g,
      aiht_probmass_in_overlap = probmass_a,
      top_in_overlap = top_in_overlap,
      clpp = clpp,
      clpa = clpa,
      top_h4_variant = top_h4$rsid,
      top_snp_pp_h4 = top_h4$SNP.PP.H4,
      graves_beta_at_top_h4 = top_common$beta_g,
      aiht_beta_at_top_h4 = top_common$beta_a,
      marginal_direction_at_top_h4 = direction,
      effect_direction_if_shared = if (classification %in% c("shared_strong_H4", "shared_suggestive_H4")) direction else NA_character_,
      quality_core = quality_core,
      quality_suggestive = quality_suggestive,
      quality_h3 = quality_h3,
      prior_sensitivity_abs_delta_h4 = abs(result$PP.H4.abf[[k]] - result$sensitivity_PP.H4[[k]]),
      classification = classification,
      is_hla = pair$is_hla,
      is_chr_x = pair$is_chr_x
    )
    h4_table[, `:=`(
      pair_id = pair$pair_id,
      graves_region = pair$g_region,
      aiht_region = pair$a_region,
      graves_cs = g_cs,
      aiht_cs = a_cs
    )]
    h4_rows[[length(h4_rows) + 1L]] <- h4_table
  }
  pair_qc[[length(pair_qc) + 1L]] <- data.table::data.table(
    pair_id = pair$pair_id,
    chromosome = pair$chrom_label,
    graves_region = pair$g_region,
    aiht_region = pair$a_region,
    graves_signals = length(g_signals),
    aiht_signals = length(a_signals),
    signal_pair_tests = length(g_signals) * length(a_signals),
    n_snps_graves = nrow(g),
    n_snps_aiht = nrow(a),
    n_shared_snps = nrow(common),
    shared_fraction_graves = nrow(common) / nrow(g),
    shared_fraction_aiht = nrow(common) / nrow(a),
    allele_mismatches = allele_mismatch,
    is_chr_x = pair$is_chr_x
  )
  log_message("T06 coloc ", pair_index, "/", nrow(ordinary_pairs), ": ", pair$pair_id)
}

coloc_result <- data.table::rbindlist(coloc_rows, fill = TRUE)
h4_variants <- data.table::rbindlist(h4_rows, fill = TRUE)
pair_qc <- data.table::rbindlist(pair_qc, fill = TRUE)
data.table::setorder(coloc_result, chromosome, graves_region, graves_cs, aiht_cs)
write_tsv(coloc_result, file.path(out_dir, "gwas_gwas_susie_coloc.tsv"))
write_tsv(h4_variants, file.path(out_dir, "gwas_gwas_h4_variant_posteriors.tsv.gz"))
write_tsv(pair_qc, file.path(out_dir, "gwas_gwas_region_pair_qc.tsv"))

region_architecture <- coloc_result[, .(
  n_signal_pair_tests = .N,
  n_shared_strong_H4 = sum(classification == "shared_strong_H4"),
  n_shared_suggestive_H4 = sum(classification == "shared_suggestive_H4"),
  n_distinct_signals_H3 = sum(classification == "distinct_signals_H3"),
  n_unresolved = sum(!classification %in% c("shared_strong_H4", "shared_suggestive_H4", "distinct_signals_H3")),
  max_PP.H4 = max(PP.H4, na.rm = TRUE),
  max_PP.H3 = max(PP.H3, na.rm = TRUE),
  shared_effect_directions = paste(sort(unique(effect_direction_if_shared[classification == "shared_strong_H4"])), collapse = ";")
), by = .(pair_id, chromosome, graves_region, aiht_region, is_chr_x)]
region_architecture[, architecture := data.table::fcase(
  n_shared_strong_H4 > 0 & n_distinct_signals_H3 > 0, "mixed_shared_and_distinct",
  n_shared_strong_H4 > 0, "shared_only",
  n_distinct_signals_H3 > 0, "distinct_only",
  n_shared_suggestive_H4 > 0, "suggestive_shared",
  default = "unresolved"
)]
data.table::setorder(region_architecture, chromosome, graves_region)
write_tsv(region_architecture, file.path(out_dir, "gwas_gwas_region_architecture.tsv"))

# Map the 16 T04B lead loci to the official single-trait fine-mapping and direct-coloc coverage.
cc_loci <- data.table::fread(file.path(project_root, "outcome_data", "ccgwas", "summary", "provisional_lead_loci.tsv"))
cc_annot <- data.table::fread(
  file.path(project_root, "outcome_data", "ccgwas", "input", "matched_variant_annotation.tsv.gz"),
  select = c("SNP", "CHR", "BP")
)
cc <- merge(
  cc_loci,
  unique(cc_annot[SNP %in% cc_loci$lead_SNP, .(lead_SNP = SNP, lead_chr = CHR, lead_position = BP)]),
  by = "lead_SNP", all.x = TRUE
)
coverage <- lapply(seq_len(nrow(cc)), function(i) {
  row <- cc[i]
  g_hit_regions <- unique(signal_catalog[
    trait_id == graves & chrom_label == as.character(row$lead_chr) &
      region_start <= row$lead_position & region_end >= row$lead_position,
    region
  ])
  a_hit_regions <- unique(signal_catalog[
    trait_id == aiht & chrom_label == as.character(row$lead_chr) &
      region_start <= row$lead_position & region_end >= row$lead_position,
    region
  ])
  coloc_hits <- coloc_result[
    chromosome == as.character(row$lead_chr) &
      graves_region %in% g_hit_regions & aiht_region %in% a_hit_regions
  ]
  architecture <- if (nrow(coloc_hits)) {
    has_shared <- any(coloc_hits$classification == "shared_strong_H4")
    has_distinct <- any(coloc_hits$classification == "distinct_signals_H3")
    if (has_shared && has_distinct) "mixed_shared_and_distinct" else if (has_shared) {
      "shared_only"
    } else if (has_distinct) {
      "distinct_only"
    } else if (any(coloc_hits$classification == "shared_suggestive_H4")) {
      "suggestive_shared"
    } else {
      "unresolved"
    }
  } else if (length(g_hit_regions) && length(a_hit_regions)) {
    "unresolved_no_signal_pair_result"
  } else if (length(g_hit_regions) || length(a_hit_regions)) {
    "unresolved_one_trait_finemapped"
  } else {
    "unresolved_no_official_finemap_region"
  }
  priority_hits <- if (nrow(coloc_hits)) {
    if (any(coloc_hits$classification == "shared_strong_H4")) {
      coloc_hits[classification == "shared_strong_H4"][which.max(PP.H4)]
    } else if (any(coloc_hits$classification == "distinct_signals_H3")) {
      coloc_hits[classification == "distinct_signals_H3"][which.max(PP.H3)]
    } else if (any(coloc_hits$classification == "shared_suggestive_H4")) {
      coloc_hits[classification == "shared_suggestive_H4"][which.max(PP.H4)]
    } else {
      coloc_hits[which.max(PP.H4)]
    }
  } else {
    data.table::data.table()
  }
  data.table::data.table(
    ccgwas_locus_id = row$locus_id,
    lead_SNP = row$lead_SNP,
    lead_nearest_gene = row$lead_nearest_genes_graves,
    chromosome = row$lead_chr,
    lead_position = row$lead_position,
    graves_finemap_regions = length(g_hit_regions),
    aiht_finemap_regions = length(a_hit_regions),
    graves_signal_count = signal_catalog[trait_id == graves & region %in% g_hit_regions, .N],
    aiht_signal_count = signal_catalog[trait_id == aiht & region %in% a_hit_regions, .N],
    direct_coloc_available = nrow(coloc_hits) > 0L,
    n_shared_strong_H4 = if (nrow(coloc_hits)) sum(coloc_hits$classification == "shared_strong_H4") else 0L,
    n_shared_suggestive_H4 = if (nrow(coloc_hits)) sum(coloc_hits$classification == "shared_suggestive_H4") else 0L,
    n_distinct_signals_H3 = if (nrow(coloc_hits)) sum(coloc_hits$classification == "distinct_signals_H3") else 0L,
    max_qualified_PP.H4 = if (nrow(coloc_hits[classification %in% c("shared_strong_H4", "shared_suggestive_H4")])) {
      max(coloc_hits[classification %in% c("shared_strong_H4", "shared_suggestive_H4")]$PP.H4)
    } else NA_real_,
    max_qualified_PP.H3 = if (nrow(coloc_hits[classification == "distinct_signals_H3"])) {
      max(coloc_hits[classification == "distinct_signals_H3"]$PP.H3)
    } else NA_real_,
    signal_architecture = architecture,
    priority_pair_PP.H3 = if (nrow(priority_hits)) priority_hits$PP.H3 else NA_real_,
    priority_pair_PP.H4 = if (nrow(priority_hits)) priority_hits$PP.H4 else NA_real_,
    priority_pair_classification = if (nrow(priority_hits)) priority_hits$classification else architecture,
    shared_effect_directions = if (nrow(coloc_hits[classification == "shared_strong_H4"])) {
      paste(sort(unique(coloc_hits[classification == "shared_strong_H4"]$effect_direction_if_shared)), collapse = ";")
    } else NA_character_
  )
})
coverage <- data.table::rbindlist(coverage, fill = TRUE)
data.table::setorder(coverage, ccgwas_locus_id)
write_tsv(coverage, file.path(out_dir, "t04b_locus_signal_coverage.tsv"))

summary <- data.table::rbindlist(list(
  data.table::data.table(metric = "graves_official_finemap_regions", value = data.table::uniqueN(cred95_g$region)),
  data.table::data.table(metric = "graves_susie_95pct_signals", value = nrow(cred95_g)),
  data.table::data.table(metric = "graves_strong_pure_signals", value = cred95_g[low_purity == FALSE & cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf), .N]),
  data.table::data.table(metric = "aiht_official_finemap_regions", value = data.table::uniqueN(cred95_a$region)),
  data.table::data.table(metric = "aiht_susie_95pct_signals", value = nrow(cred95_a)),
  data.table::data.table(metric = "aiht_strong_pure_signals", value = cred95_a[low_purity == FALSE & cs_log10bf > as.numeric(t06_cfg$min_cs_log10bf), .N]),
  data.table::data.table(metric = "overlapping_non_hla_region_pairs", value = nrow(ordinary_pairs)),
  data.table::data.table(metric = "gwas_gwas_signal_pair_tests", value = nrow(coloc_result)),
  data.table::data.table(metric = "shared_strong_H4", value = coloc_result[classification == "shared_strong_H4", .N]),
  data.table::data.table(metric = "shared_strong_H4_region_pairs", value = region_architecture[n_shared_strong_H4 > 0, .N]),
  data.table::data.table(metric = "shared_suggestive_H4", value = coloc_result[classification == "shared_suggestive_H4", .N]),
  data.table::data.table(metric = "distinct_signals_H3", value = coloc_result[classification == "distinct_signals_H3", .N]),
  data.table::data.table(metric = "distinct_signals_H3_region_pairs", value = region_architecture[n_distinct_signals_H3 > 0, .N]),
  data.table::data.table(metric = "mixed_shared_and_distinct_region_pairs", value = region_architecture[architecture == "mixed_shared_and_distinct", .N]),
  data.table::data.table(metric = "unresolved_signal_pairs", value = coloc_result[!classification %in% c("shared_strong_H4", "shared_suggestive_H4", "distinct_signals_H3"), .N]),
  data.table::data.table(metric = "hla_region_pairs_excluded", value = region_pairs[is_hla == TRUE, .N]),
  data.table::data.table(metric = "t04b_loci_both_traits_finemapped", value = coverage[graves_finemap_regions > 0 & aiht_finemap_regions > 0, .N]),
  data.table::data.table(metric = "t04b_loci_one_trait_finemapped", value = coverage[(graves_finemap_regions > 0) != (aiht_finemap_regions > 0), .N]),
  data.table::data.table(metric = "t04b_loci_neither_trait_finemapped", value = coverage[graves_finemap_regions == 0 & aiht_finemap_regions == 0, .N])
), fill = TRUE)
write_tsv(summary, file.path(out_dir, "t06_summary.tsv"))

classification_levels <- c(
  "shared_strong_H4", "shared_suggestive_H4", "distinct_signals_H3",
  "unresolved_posterior", "unresolved_insufficient_overlap",
  "unresolved_weak_credible_set", "unresolved_low_purity"
)
coloc_plot <- data.table::copy(coloc_result)
coloc_plot[, classification := factor(classification, levels = classification_levels)]
p1 <- ggplot2::ggplot(coloc_plot, ggplot2::aes(x = PP.H3, y = PP.H4, colour = classification)) +
  ggplot2::geom_point(size = 2.7, alpha = 0.85) +
  ggplot2::geom_vline(xintercept = as.numeric(t06_cfg$distinct_h3), linetype = 2, colour = "#777777") +
  ggplot2::geom_hline(yintercept = as.numeric(t06_cfg$strong_h4), linetype = 2, colour = "#777777") +
  ggplot2::scale_colour_manual(values = c(
    shared_strong_H4 = "#1B7837", shared_suggestive_H4 = "#5AAE61",
    distinct_signals_H3 = "#762A83", unresolved_posterior = "#4575B4",
    unresolved_insufficient_overlap = "#D8B365", unresolved_weak_credible_set = "#E08214",
    unresolved_low_purity = "#BDBDBD"
  ), drop = TRUE) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "T06 signal-level Graves–AIHT colocalization",
    subtitle = "FinnGen R13 official SuSiE signals; dashed lines mark PP = 0.80",
    x = "PP.H3 (distinct causal variants)", y = "PP.H4 (shared causal variant)", colour = "Classification"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(figure_dir, "t06_h3_h4_signal_pairs.png"), p1, width = 10.5, height = 6.2, dpi = 300)

qc_plot <- signal_catalog[, .N, by = .(
  trait_id,
  signal_qc = data.table::fcase(
    low_purity, "Low purity",
    cs_log10bf <= as.numeric(t06_cfg$min_cs_log10bf), "Weak BF",
    default = "Pure and BF > 0.9"
  )
)]
p2 <- ggplot2::ggplot(qc_plot, ggplot2::aes(x = trait_id, y = N, fill = signal_qc)) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = N), position = ggplot2::position_stack(vjust = 0.5),
    colour = "white", size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = c(
    "Pure and BF > 0.9" = "#2166AC", "Weak BF" = "#E08214", "Low purity" = "#BDBDBD"
  )) +
  ggplot2::labs(
    title = "T06 official SuSiE signal quality",
    x = NULL, y = "Number of 95% credible sets", fill = "Signal QC"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(figure_dir, "t06_susie_signal_qc.png"), p2, width = 9.2, height = 5.8, dpi = 300)

provenance <- data.table::data.table(
  run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  graves_trait = graves,
  aiht_trait = aiht,
  source_release = "FinnGen R13 official fine-mapping",
  fine_mapping_method = "SuSiE with L=10; FINEMAP region support retained",
  ld_source = "FinnGen in-sample dosage LD from official fine-mapping pipeline",
  genome_build = cfg$project$genome_build,
  primary_p1 = as.numeric(cfg$coloc$p1),
  primary_p2 = as.numeric(cfg$coloc$p2),
  primary_p12 = as.numeric(t06_cfg$primary_p12),
  sensitivity_p12 = as.numeric(t06_cfg$sensitivity_p12),
  strong_h4_threshold = as.numeric(t06_cfg$strong_h4),
  suggestive_h4_threshold = as.numeric(t06_cfg$suggestive_h4),
  distinct_h3_threshold = as.numeric(t06_cfg$distinct_h3),
  min_cs_log10bf = as.numeric(t06_cfg$min_cs_log10bf),
  min_probmass_core = as.numeric(t06_cfg$min_probmass_core),
  min_probmass_suggestive = as.numeric(t06_cfg$min_probmass_suggestive),
  coloc_version = as.character(utils::packageVersion("coloc")),
  convergence_note = "Iteration/ELBO convergence is not exposed in public summary files",
  sample_overlap_note = "Graves and AIHT are from the same FinnGen release; direct coloc does not explicitly model correlated summary-statistic error"
)
write_tsv(provenance, file.path(out_dir, "t06_provenance.tsv"))

log_message(
  "T06 complete: ", nrow(signal_catalog), " official SuSiE signals; ",
  nrow(coloc_result), " signal-pair tests; strong H4 = ",
  coloc_result[classification == "shared_strong_H4", .N], "; distinct H3 = ",
  coloc_result[classification == "distinct_signals_H3", .N]
)
