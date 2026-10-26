#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "coloc", "Rsamtools", "GenomicRanges", "IRanges"))
cfg <- read_analysis_config(project_root)
args <- parse_cli()
t10 <- cfg$t10_qtl_susie_coloc
if (is.null(t10)) stop("Missing t10_qtl_susie_coloc config.", call. = FALSE)

profile <- as.character(args$profile %||% "t10")
p1a_profiles <- c("t11a_wave1", "t11b_wave2")
if (!profile %in% c("t10", p1a_profiles)) {
  stop("--profile must be t10, t11a_wave1, or t11b_wave2.", call. = FALSE)
}
arm <- as.character(args$arm %||% if (profile %in% p1a_profiles) "non_hla" else "")
if (!arm %in% c("non_hla", "hla")) stop("--arm must be non_hla or hla.", call. = FALSE)
if (profile %in% p1a_profiles && arm != "non_hla") stop("P1A wave profiles are non-HLA only.", call. = FALSE)
arm_value <- arm
workers <- suppressWarnings(as.integer(args$workers %||% if (arm == "hla") 1L else 4L))
if (is.na(workers) || workers < 1L) stop("--workers must be a positive integer.", call. = FALSE)
if (arm == "hla" && workers != 1L) {
  log_message("HLA is forced to one worker.")
  workers <- 1L
}
if (.Platform$OS.type == "windows") workers <- 1L

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)
data.table::setDTthreads(1L)

out_dir <- project_path(
  project_root,
  args[["output-dir"]] %||% if (profile %in% p1a_profiles) {
    file.path(cfg$project$outcome_dir, if (profile == "t11a_wave1") "t11a" else "t11b", "direct_coloc")
  } else {
    file.path(cfg$project$outcome_dir, "t10", arm)
  }
)
ensure_dir(out_dir)
catalog_path <- file.path(project_root, cfg$project$outcome_dir, "t06", "susie_signal_catalog.tsv")
arm_is_hla <- identical(arm, "hla")
if (!file.exists(catalog_path)) stop("Missing GWAS SuSiE signal catalog.", call. = FALSE)

if (profile %in% p1a_profiles) {
  is_wave1 <- identical(profile, "t11a_wave1")
  stage_dir <- if (is_wave1) "t11a" else "t11b"
  rows_default <- if (is_wave1) {
    file.path(cfg$project$outcome_dir, "t11", "p1a_wave1_candidate_rows.tsv.gz")
  } else {
    file.path(cfg$project$outcome_dir, "t11a", "p1a_wave2_candidate_rows.tsv.gz")
  }
  verification_default <- file.path(
    cfg$project$outcome_dir, stage_dir,
    if (is_wave1) "p1a_wave1_verification.tsv" else "p1a_wave2_verification.tsv"
  )
  rows_path <- project_path(
    project_root, args[["target-rows"]] %||% rows_default
  )
  verification_path <- project_path(
    project_root, args$verification %||% verification_default
  )
  if (!file.exists(rows_path) || !file.exists(verification_path)) stop("Missing P1A wave rows or verification.", call. = FALSE)
  screen <- data.table::fread(rows_path, showProgress = FALSE)
  verification <- data.table::fread(verification_path, showProgress = FALSE)
  expected_files <- if (is_wave1) 7L else 9L
  if (nrow(verification) != expected_files || any(verification$status != "ok") || any(!verification$index_ok)) {
    stop(stage_dir, " verification gate is not ", expected_files, "/", expected_files, " ok.", call. = FALSE)
  }
  targets <- unique(screen[, .(
    source_path, modality, level, cell_type, locus_id,
    chrom = manifest_chrom, start, end, is_hla, feature_id
  )])
  targets[, selection_reasons := if (is_wave1) {
    "P1A_wave1_dual_trait_concordant_pure_CS"
  } else {
    "P1A_wave2_dual_trait_concordant_pure_CS"
  }]
  targets <- merge(
    targets, verification[, .(source_path, local_data_path, index_ok, status)],
    by = "source_path", all.x = TRUE, sort = FALSE
  )
  if (anyNA(targets$local_data_path) || any(targets$status != "ok") || any(!targets$index_ok)) {
    stop(stage_dir, " targets do not map cleanly to verified slices.", call. = FALSE)
  }
  screen_prefix <- if (is_wave1) "W1SCR" else "W2SCR"
  analysis_label <- if (is_wave1) "T11A wave 1" else "T11B wave 2"
} else {
  frozen_path <- file.path(project_root, cfg$project$outcome_dir, "t08", "p0_core_frozen_regions.tsv")
  inventory_path <- file.path(project_root, cfg$project$outcome_dir, "t09", "t09_file_inventory.tsv")
  for (path in c(frozen_path, inventory_path)) {
    if (!file.exists(path)) stop("Missing T10 input: ", path, call. = FALSE)
  }
  frozen <- data.table::fread(frozen_path, showProgress = FALSE)[is_hla == arm_is_hla]
  inventory <- data.table::fread(inventory_path, showProgress = FALSE)[arm == arm_value]
  if (!nrow(frozen) || !nrow(inventory)) stop("No frozen/inventory rows for arm: ", arm, call. = FALSE)
  if (any(inventory$status != "ok")) stop("T09 inventory contains non-ok inputs for arm: ", arm, call. = FALSE)
  targets <- frozen[, .(
    feature_id = unlist(strsplit(features, ";", fixed = TRUE))
  ), by = .(
    source_path, modality, level, cell_type, locus_id, chrom, start, end,
    is_hla, selection_reasons
  )]
  targets <- unique(targets[nzchar(feature_id)])
  targets <- merge(
    targets,
    inventory[, .(source_path, local_data_path, bgzip_integrity_ok, index_ok, status)],
    by = "source_path", all.x = TRUE, sort = FALSE
  )
  if (anyNA(targets$local_data_path) || any(!targets$bgzip_integrity_ok) || any(!targets$index_ok)) {
    stop("Frozen targets do not map cleanly to verified T09 files.", call. = FALSE)
  }
  screen_files <- list.files(
    file.path(project_root, cfg$project$outcome_dir, "coloc_r13"),
    pattern = "\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE
  )
  if (!length(screen_files)) stop("No R13 regional screen files found.", call. = FALSE)
  screen <- data.table::rbindlist(lapply(screen_files, data.table::fread, showProgress = FALSE), fill = TRUE)
  screen_keys <- unique(targets[, .(modality, level, cell_type, locus_id, feature_id)])
  screen <- merge(screen, screen_keys, by = c("modality", "level", "cell_type", "locus_id", "feature_id"))
  screen_prefix <- paste0(if (arm_is_hla) "H" else "N", "SCR")
  analysis_label <- paste0("T10 ", arm)
}
targets[, chrom_label := sub("^chr", "", chrom)]
screen <- screen[trait_id %in% c(cfg$primary_comparison$graves, cfg$primary_comparison$autoimmune_hypothyroidism)]
data.table::setorder(screen, trait_id, locus_id, modality, level, cell_type, feature_id, -PP.H4.abf)
screen[, screen_id := paste0(screen_prefix, sprintf("%04d", seq_len(.N)))]
if (!nrow(screen)) stop("No R13 screen rows matched frozen targets.", call. = FALSE)
write_tsv(screen, file.path(out_dir, "regional_screening_rows.tsv.gz"))

signal_catalog <- data.table::fread(catalog_path, showProgress = FALSE)
signal_catalog <- signal_catalog[trait_id %in% unique(screen$trait_id)]

# Map each regional screen row to official FinnGen GWAS fine-mapping regions.
screen_loci <- unique(screen[, .(
  screen_id, trait_id, locus_id, chrom_label = sub("^chr", "", as.character(chrom)),
  locus_start = as.integer(start), locus_end = as.integer(end)
)])
gwas_map <- merge(screen_loci, signal_catalog, by = c("trait_id", "chrom_label"), allow.cartesian = TRUE)
gwas_map <- gwas_map[region_start <= locus_end & region_end >= locus_start]
data.table::setorder(gwas_map, screen_id, region, cs)

required_snp <- c(
  "trait", "region", "rsid", "position", "allele1", "allele2", "beta", "se", "p",
  "prob", "cs", "low_purity", paste0("alpha", 1:10), paste0("lbf_variable", 1:10)
)

read_gwas_region <- function(trait_id, region_value, chrom_label, region_start, region_end) {
  path <- gwas_susie_path(project_root, trait_id)
  header <- read_gzip_header(path)
  missing <- setdiff(required_snp, header)
  if (length(missing)) stop("GWAS SuSiE columns missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x <- scan_tabix_dt(path, chrom_label, region_start, region_end)
  x <- x[trait == trait_id & region == region_value, ..required_snp]
  if (!nrow(x)) stop("No exact GWAS rows for ", trait_id, " ", region_value, call. = FALSE)
  if (anyDuplicated(x$rsid)) stop("Duplicate GWAS rsid in ", trait_id, " ", region_value, call. = FALSE)
  x
}

region_meta <- unique(gwas_map[, .(trait_id, region, chrom_label, region_start, region_end)])
gwas_cache <- list()
if (nrow(region_meta)) {
  for (i in seq_len(nrow(region_meta))) {
    m <- region_meta[i]
    key <- paste(m$trait_id, m$region, sep = "||")
    gwas_cache[[key]] <- read_gwas_region(m$trait_id, m$region, m$chrom_label, m$region_start, m$region_end)
  }
}
log_message("Loaded ", length(gwas_cache), " official GWAS SuSiE regions for ", arm, ".")

read_selected_qtl <- function(path, features) {
  feature_file <- tempfile("t10_features_", fileext = ".txt")
  on.exit(unlink(feature_file), add = TRUE)
  data.table::fwrite(data.table::data.table(feature_id = unique(features)), feature_file, col.names = FALSE)
  command <- paste(
    "gzip -dc", shQuote(path), "|",
    "awk -F '\\t' 'NR==FNR {keep[$1]=1; next} FNR==1 || ($2 in keep)'",
    shQuote(feature_file), "-"
  )
  x <- data.table::fread(cmd = command, select = required_snp, showProgress = FALSE)
  x
}

probmass <- function(x, signal, shared_ids) {
  alpha_col <- paste0("alpha", signal)
  if (!alpha_col %in% names(x)) return(NA_real_)
  sum(as.numeric(x[[alpha_col]])[x$rsid %in% shared_ids], na.rm = TRUE)
}

signal_low_purity <- function(x, signal) {
  values <- as.logical(x[cs == signal]$low_purity)
  values <- values[!is.na(values)]
  if (!length(values)) NA else any(values)
}

process_source <- function(source_path_value) {
  unit_targets <- targets[source_path == source_path_value]
  local_path <- unique(unit_targets$local_data_path)
  if (length(local_path) != 1L || !file.exists(local_path)) stop("Invalid local QTL path.", call. = FALSE)
  qtl <- read_selected_qtl(local_path, unit_targets$feature_id)
  if (anyDuplicated(qtl[, .(region, rsid)])) stop("Duplicate QTL region/rsid rows in ", local_path, call. = FALSE)
  coloc_bf_bf <- getFromNamespace("coloc.bf_bf", "coloc")
  pair_rows <- list()
  qc_rows <- list()

  for (i in seq_len(nrow(unit_targets))) {
    target <- unit_targets[i]
    q <- qtl[region == target$feature_id]
    screen_subset <- screen[
      modality == target$modality & level == target$level & cell_type == target$cell_type &
        locus_id == target$locus_id & feature_id == target$feature_id
    ]
    q_signals <- sort(unique(as.integer(q$cs[is.finite(q$cs) & q$cs > 0L])))
    qc_status <- if (!nrow(q)) "qtl_feature_absent" else if (!length(q_signals)) "qtl_no_credible_set_in_slice" else "ready"
    qc_rows[[length(qc_rows) + 1L]] <- data.table::data.table(
      source_path = source_path_value, modality = target$modality, level = target$level,
      cell_type = target$cell_type, locus_id = target$locus_id, feature_id = target$feature_id,
      is_hla = target$is_hla, n_qtl_rows = nrow(q), n_qtl_signals = length(q_signals),
      qtl_signals = if (length(q_signals)) paste(q_signals, collapse = ";") else NA_character_,
      n_screen_rows = nrow(screen_subset), feature_status = qc_status
    )
    if (qc_status != "ready" || !nrow(screen_subset)) next

    for (screen_index in seq_len(nrow(screen_subset))) {
      s <- screen_subset[screen_index]
      maps <- gwas_map[screen_id == s$screen_id]
      if (!nrow(maps)) next
      for (gwas_region in unique(maps$region)) {
        map_rows <- maps[region == gwas_region]
        g_signals <- sort(unique(as.integer(map_rows$cs)))
        g <- gwas_cache[[paste(s$trait_id, gwas_region, sep = "||")]]
        if (is.null(g) || !nrow(g) || !length(g_signals)) next

        g_keep <- c("rsid", "position", "allele1", "allele2", "beta", "se", "cs",
                    paste0("alpha", g_signals), paste0("lbf_variable", g_signals))
        q_keep <- c("rsid", "position", "allele1", "allele2", "beta", "se", "cs",
                    paste0("alpha", q_signals), paste0("lbf_variable", q_signals))
        gs <- g[, ..g_keep]
        qs <- q[, ..q_keep]
        data.table::setnames(gs, setdiff(names(gs), "rsid"), paste0(setdiff(names(gs), "rsid"), "_g"))
        data.table::setnames(qs, setdiff(names(qs), "rsid"), paste0(setdiff(names(qs), "rsid"), "_q"))
        common_all <- merge(gs, qs, by = "rsid")
        mismatch <- common_all[
          position_g != position_q | allele1_g != allele1_q | allele2_g != allele2_q
        ]
        common <- common_all[
          position_g == position_q & allele1_g == allele1_q & allele2_g == allele2_q
        ]
        if (nrow(common) < as.integer(t10$minimum_shared_snps)) next

        g_lbf_cols <- paste0("lbf_variable", g_signals, "_g")
        q_lbf_cols <- paste0("lbf_variable", q_signals, "_q")
        bf_g <- t(as.matrix(common[, ..g_lbf_cols]))
        bf_q <- t(as.matrix(common[, ..q_lbf_cols]))
        storage.mode(bf_g) <- "double"
        storage.mode(bf_q) <- "double"
        colnames(bf_g) <- common$rsid
        colnames(bf_q) <- common$rsid
        primary <- tryCatch(
          coloc_bf_bf(
            bf_g, bf_q, p1 = as.numeric(cfg$coloc$p1), p2 = as.numeric(cfg$coloc$p2),
            p12 = as.numeric(t10$primary_p12), trim_by_posterior = FALSE
          ), error = function(e) e
        )
        if (inherits(primary, "error")) next
        sensitivity <- coloc_bf_bf(
          bf_g, bf_q, p1 = as.numeric(cfg$coloc$p1), p2 = as.numeric(cfg$coloc$p2),
          p12 = as.numeric(t10$sensitivity_p12), trim_by_posterior = FALSE
        )
        result <- merge(
          data.table::as.data.table(primary$summary),
          data.table::as.data.table(sensitivity$summary)[, .(
            idx1, idx2, sensitivity_PP.H3 = PP.H3.abf, sensitivity_PP.H4 = PP.H4.abf
          )], by = c("idx1", "idx2"), sort = FALSE
        )
        result[, `:=`(gwas_cs = g_signals[idx1], qtl_cs = q_signals[idx2])]

        for (k in seq_len(nrow(result))) {
          g_cs <- result$gwas_cs[[k]]
          q_cs <- result$qtl_cs[[k]]
          g_meta <- map_rows[cs == g_cs][1]
          shared_ids <- common$rsid
          g_mass <- probmass(g, g_cs, shared_ids)
          q_mass <- probmass(q, q_cs, shared_ids)
          cs_overlap <- common[cs_g == g_cs & cs_q == q_cs, .N]
          h4_col <- if (nrow(result) == 1L) "SNP.PP.H4.abf" else paste0("SNP.PP.H4.row", k)
          h4 <- data.table::as.data.table(primary$results)[, .(rsid = snp, SNP.PP.H4 = get(h4_col))]
          top <- h4[which.max(SNP.PP.H4)]
          top_common <- common[rsid == top$rsid][1]
          product <- as.numeric(top_common$beta_g) * as.numeric(top_common$beta_q)
          direction <- if (!is.finite(product) || product == 0) "undetermined" else if (product > 0) {
            "disease_risk_allele_increases_feature"
          } else {
            "disease_risk_allele_decreases_feature"
          }
          g_low <- as.logical(g_meta$low_purity)
          q_low <- signal_low_purity(q, q_cs)
          quality_core <- !isTRUE(g_low) && !isTRUE(q_low) &&
            is.finite(g_meta$cs_log10bf) && g_meta$cs_log10bf > as.numeric(t10$min_gwas_cs_log10bf) &&
            is.finite(g_mass) && g_mass > as.numeric(t10$min_probmass_core) &&
            is.finite(q_mass) && q_mass > as.numeric(t10$min_probmass_core) && cs_overlap > 0L
          quality_suggestive <- !isTRUE(g_low) && !isTRUE(q_low) &&
            is.finite(g_meta$cs_log10bf) && g_meta$cs_log10bf > as.numeric(t10$min_gwas_cs_log10bf) &&
            is.finite(g_mass) && g_mass > as.numeric(t10$min_probmass_suggestive) &&
            is.finite(q_mass) && q_mass > as.numeric(t10$min_probmass_suggestive) && cs_overlap > 0L
          quality_h3 <- !isTRUE(g_low) && !isTRUE(q_low) &&
            is.finite(g_meta$cs_log10bf) && g_meta$cs_log10bf > as.numeric(t10$min_gwas_cs_log10bf) &&
            is.finite(g_mass) && g_mass > as.numeric(t10$min_probmass_core) &&
            is.finite(q_mass) && q_mass > as.numeric(t10$min_probmass_core)
          classification <- if (isTRUE(g_low) || isTRUE(q_low)) {
            "unresolved_low_purity"
          } else if (!is.finite(g_meta$cs_log10bf) || g_meta$cs_log10bf <= as.numeric(t10$min_gwas_cs_log10bf)) {
            "unresolved_weak_gwas_credible_set"
          } else if (quality_core && result$PP.H4.abf[[k]] >= as.numeric(t10$strong_h4) &&
                     result$sensitivity_PP.H4[[k]] >= as.numeric(t10$strong_h4)) {
            "shared_strong_H4"
          } else if (quality_suggestive && result$PP.H4.abf[[k]] >= as.numeric(t10$suggestive_h4) &&
                     result$sensitivity_PP.H4[[k]] >= as.numeric(t10$suggestive_h4)) {
            "shared_suggestive_H4"
          } else if (quality_h3 && result$PP.H3.abf[[k]] >= as.numeric(t10$distinct_h3) &&
                     result$sensitivity_PP.H3[[k]] >= as.numeric(t10$distinct_h3)) {
            "distinct_signals_H3"
          } else if (!quality_suggestive) {
            "unresolved_insufficient_overlap"
          } else {
            "unresolved_posterior"
          }
          pair_rows[[length(pair_rows) + 1L]] <- data.table::data.table(
            screen_id = s$screen_id, trait_id = s$trait_id, locus_id = target$locus_id,
            modality = target$modality, level = target$level, cell_type = target$cell_type,
            feature_id = target$feature_id, is_hla = target$is_hla,
            gwas_region = gwas_region, gwas_cs = g_cs, qtl_cs = q_cs,
            n_gwas_snps = nrow(g), n_qtl_snps_in_slice = nrow(q),
            n_shared_snps = nrow(common), allele_mismatches_excluded = nrow(mismatch),
            PP.H0 = result$PP.H0.abf[[k]], PP.H1 = result$PP.H1.abf[[k]],
            PP.H2 = result$PP.H2.abf[[k]], PP.H3 = result$PP.H3.abf[[k]],
            PP.H4 = result$PP.H4.abf[[k]], sensitivity_PP.H3 = result$sensitivity_PP.H3[[k]],
            sensitivity_PP.H4 = result$sensitivity_PP.H4[[k]],
            gwas_probmass_in_shared = g_mass, qtl_probmass_in_shared = q_mass,
            cs_overlap = cs_overlap, gwas_cs_log10bf = g_meta$cs_log10bf,
            gwas_low_purity = g_low, qtl_low_purity = q_low,
            top_h4_variant = top$rsid, top_snp_pp_h4 = top$SNP.PP.H4,
            gwas_beta_at_top_h4 = as.numeric(top_common$beta_g),
            qtl_beta_at_top_h4 = as.numeric(top_common$beta_q),
            effect_direction_if_shared = if (classification %in% c("shared_strong_H4", "shared_suggestive_H4")) direction else NA_character_,
            direction_at_top_h4 = direction, quality_core = quality_core,
            quality_suggestive = quality_suggestive, quality_h3 = quality_h3,
            classification = classification, method = "direct_SuSiE_BF_signal_pair_coloc",
            primary_p12 = as.numeric(t10$primary_p12), sensitivity_p12 = as.numeric(t10$sensitivity_p12),
            genome_build = cfg$project$genome_build
          )
        }
      }
    }
  }
  list(
    pairs = data.table::rbindlist(pair_rows, fill = TRUE),
    qc = data.table::rbindlist(qc_rows, fill = TRUE)
  )
}

source_units <- unique(targets$source_path)
log_message("Processing ", length(source_units), " ", arm, " QTL source files with ", workers, " worker(s).")
if (workers > 1L) {
  results <- parallel::mclapply(source_units, process_source, mc.cores = workers, mc.preschedule = FALSE)
} else {
  results <- lapply(source_units, process_source)
}
pairs <- data.table::rbindlist(lapply(results, `[[`, "pairs"), fill = TRUE)
feature_qc <- data.table::rbindlist(lapply(results, `[[`, "qc"), fill = TRUE)
if (nrow(pairs)) data.table::setorder(pairs, screen_id, gwas_region, gwas_cs, qtl_cs)
data.table::setorder(feature_qc, locus_id, modality, level, cell_type, feature_id)
write_tsv(pairs, file.path(out_dir, "gwas_qtl_susie_signal_pairs.tsv.gz"))
write_tsv(feature_qc, file.path(out_dir, "qtl_feature_qc.tsv.gz"))

screen_coverage <- screen[, .(screen_id, trait_id, locus_id, modality, level, cell_type, feature_id, is_hla)]
screen_coverage <- merge(
  screen_coverage,
  if (nrow(pairs)) pairs[, .(
    n_signal_pairs = .N, n_strong_h4 = sum(classification == "shared_strong_H4"),
    n_suggestive_h4 = sum(classification == "shared_suggestive_H4"),
    n_distinct_h3 = sum(classification == "distinct_signals_H3"),
    max_PP.H4 = max(PP.H4), max_PP.H3 = max(PP.H3)
  ), by = "screen_id"] else data.table::data.table(screen_id = character()),
  by = "screen_id", all.x = TRUE
)
screen_coverage[is.na(n_signal_pairs), `:=`(
  n_signal_pairs = 0L, n_strong_h4 = 0L, n_suggestive_h4 = 0L, n_distinct_h3 = 0L
)]
write_tsv(screen_coverage, file.path(out_dir, "screen_pair_coverage.tsv"))

software <- data.table::data.table(
  component = c("R", "data.table", "coloc", "Rsamtools", "arm", "workers"),
  version = c(
    R.version.string, as.character(packageVersion("data.table")), as.character(packageVersion("coloc")),
    as.character(packageVersion("Rsamtools")), arm, as.character(workers)
  )
)
write_tsv(software, file.path(out_dir, "software_provenance.tsv"))
log_message(analysis_label, " complete: ", nrow(screen), " screen rows; ", nrow(pairs), " direct signal pairs.")
