#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)
traits <- read_traits(project_root, cfg)
graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
harmonized_dir <- file.path(project_root, cfg$project$outcome_dir, "harmonized_gwas")

read_primary <- function(trait_id, suffix) {
  path <- file.path(harmonized_dir, paste0(trait_id, ".harmonized.tsv.gz"))
  if (!file.exists(path)) stop("Run codes/01_harmonize_gwas.R first; missing ", path, call. = FALSE)
  x <- data.table::fread(
    path,
    select = c("variant_id", "chrom", "pos", "rsid", "ref", "alt", "beta", "se", "p", "maf", "is_snp", "palindromic", "is_hla"),
    showProgress = TRUE
  )
  data.table::setnames(x, c("beta", "se", "p"), paste0(c("beta_", "se_", "p_"), suffix))
  x
}

g <- read_primary(graves_id, "graves")
h <- read_primary(aiht_id, "aiht")
p_threshold <- as.numeric(cfg$gwas_qc$genome_wide_p)
window <- as.integer(cfg$gwas_qc$locus_half_window_bp)

signals <- data.table::rbindlist(list(
  g[p_graves < p_threshold, .(variant_id, chrom, pos, p = p_graves, source_trait = graves_id, is_hla)],
  h[p_aiht < p_threshold, .(variant_id, chrom, pos, p = p_aiht, source_trait = aiht_id, is_hla)]
))
if (!nrow(signals)) stop("No genome-wide significant variants found at configured threshold.", call. = FALSE)
data.table::setorder(signals, p)
signals <- signals[!duplicated(variant_id)]

greedy_distance_leads <- function(x, half_window) {
  selected <- integer()
  for (i in seq_len(nrow(x))) {
    if (!length(selected) || !any(x$chrom[selected] == x$chrom[i] & abs(x$pos[selected] - x$pos[i]) <= half_window)) {
      selected <- c(selected, i)
    }
  }
  x[selected]
}

non_hla <- greedy_distance_leads(signals[is_hla == FALSE], window)
candidate_windows <- non_hla[, .(
  chrom = as.integer(chrom),
  start = pmax(1L, pos - window),
  end = pos + window,
  seed_variant = variant_id,
  seed_p = p
)]
if (any(signals$is_hla)) {
  candidate_windows <- data.table::rbindlist(list(
    candidate_windows,
    data.table::data.table(
      chrom = as.integer(cfg$hla$chromosome), start = as.integer(cfg$hla$start), end = as.integer(cfg$hla$end),
      seed_variant = signals[is_hla == TRUE][which.min(p), variant_id], seed_p = min(signals[is_hla == TRUE]$p)
    )
  ))
}
data.table::setorder(candidate_windows, chrom, start, end)

merge_windows <- function(w) {
  out <- list()
  for (chr_value in unique(w$chrom)) {
    z <- w[chrom == chr_value]
    current_start <- z$start[1]
    current_end <- z$end[1]
    seeds <- z$seed_variant[1]
    if (nrow(z) > 1L) {
      for (i in 2:nrow(z)) {
        if (z$start[i] <= current_end) {
          current_end <- max(current_end, z$end[i])
          seeds <- c(seeds, z$seed_variant[i])
        } else {
          out[[length(out) + 1L]] <- data.table::data.table(chrom = chr_value, start = current_start, end = current_end, seed_variants = paste(seeds, collapse = ";"))
          current_start <- z$start[i]
          current_end <- z$end[i]
          seeds <- z$seed_variant[i]
        }
      }
    }
    out[[length(out) + 1L]] <- data.table::data.table(chrom = chr_value, start = current_start, end = current_end, seed_variants = paste(seeds, collapse = ";"))
  }
  data.table::rbindlist(out)
}

loci <- merge_windows(candidate_windows)
data.table::setorder(loci, chrom, start)
loci[, locus_id := sprintf("L%03d_chr%s_%d_%d", .I, chrom, start, end)]
loci[, is_hla := chrom == as.integer(cfg$hla$chromosome) & start <= as.integer(cfg$hla$end) & end >= as.integer(cfg$hla$start)]

shared <- merge(
  g[, .(variant_id, chrom, pos, rsid, ref, alt, maf, is_snp, palindromic, beta_graves, se_graves, p_graves)],
  h[, .(variant_id, beta_aiht, se_aiht, p_aiht)],
  by = "variant_id",
  all = FALSE
)

describe_locus <- function(chr_value, start_value, end_value) {
  z <- shared[as.integer(chrom) == chr_value & pos >= start_value & pos <= end_value]
  if (!nrow(z)) {
    return(data.table::data.table(
      lead_variant = NA_character_, lead_pos = NA_integer_, ref = NA_character_, alt = NA_character_,
      beta_graves = NA_real_, se_graves = NA_real_, p_graves = NA_real_,
      beta_aiht = NA_real_, se_aiht = NA_real_, p_aiht = NA_real_,
      provisional_class = "no_shared_variant_after_qc",
      minimum_p_graves_in_locus = NA_real_, minimum_p_aiht_in_locus = NA_real_,
      lava_snps = NA_character_
    ))
  }
  z[, combined_p := pmin(p_graves, p_aiht, na.rm = TRUE)]
  lead <- z[which.min(combined_p)]
  g_sig <- lead$p_graves < p_threshold
  h_sig <- lead$p_aiht < p_threshold
  g_support <- lead$p_graves < as.numeric(cfg$gwas_qc$cross_trait_support_p)
  h_support <- lead$p_aiht < as.numeric(cfg$gwas_qc$cross_trait_support_p)
  same <- sign(lead$beta_graves) == sign(lead$beta_aiht)
  class <- if (g_sig && h_sig && same) {
    "shared_gws_same_direction"
  } else if (g_sig && h_sig && !same) {
    "shared_gws_opposite_direction"
  } else if ((g_sig || h_sig) && g_support && h_support && same) {
    "cross_trait_supported_same_direction"
  } else if ((g_sig || h_sig) && g_support && h_support && !same) {
    "cross_trait_supported_opposite_direction"
  } else if (g_sig && !h_sig) {
    "graves_threshold_only_provisional"
  } else if (h_sig && !g_sig) {
    "aiht_threshold_only_provisional"
  } else {
    "ambiguous"
  }
  data.table::data.table(
    lead_variant = lead$variant_id,
    lead_pos = lead$pos,
    ref = lead$ref,
    alt = lead$alt,
    beta_graves = lead$beta_graves,
    se_graves = lead$se_graves,
    p_graves = lead$p_graves,
    beta_aiht = lead$beta_aiht,
    se_aiht = lead$se_aiht,
    p_aiht = lead$p_aiht,
    provisional_class = class,
    minimum_p_graves_in_locus = min(z$p_graves, na.rm = TRUE),
    minimum_p_aiht_in_locus = min(z$p_aiht, na.rm = TRUE),
    lava_snps = paste(
      unique(na.omit(z[
        is_snp & !palindromic & maf >= as.numeric(cfg$gwas_qc$maf_min_ld_methods) & grepl("^rs[0-9]+$", rsid),
        rsid
      ])),
      collapse = ";"
    )
  )
}

descriptions <- data.table::rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  describe_locus(loci$chrom[i], loci$start[i], loci$end[i])
}), fill = TRUE)
loci <- cbind(loci, descriptions)
loci[, `:=`(
  genome_build = cfg$project$genome_build,
  locus_definition = "distance-based union of primary-trait GWS signals",
  requires_local_rg_confirmation = TRUE
)]

out_dir <- file.path(project_root, cfg$project$outcome_dir, "loci")
lava_loci <- loci[, .(LOC = locus_id, CHR = chrom, START = start, STOP = end, SNPS = lava_snps)]
write_tsv(lava_loci, file.path(out_dir, "primary_shared_divergent_loci.lava"))
loci[, lava_snps := NULL]
write_tsv(loci, file.path(out_dir, "primary_shared_divergent_loci.tsv"))
log_message("Defined ", nrow(loci), " loci. Threshold-only labels remain provisional until local-rg/heterogeneity evidence is added.")
