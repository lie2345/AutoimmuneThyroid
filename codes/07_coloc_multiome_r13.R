#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "coloc", "Rsamtools", "GenomicRanges", "IRanges"))
cfg <- read_analysis_config(project_root)
args <- parse_cli()

parse_positive_integer <- function(value, label) {
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop(label, " must be a positive integer.", call. = FALSE)
  }
  parsed
}

workers <- parse_positive_integer(args$workers %||% cfg$coloc$workers_default %||% 4L, "--workers")
hla_workers <- parse_positive_integer(args$hla_workers %||% cfg$coloc$hla_workers %||% 1L, "--hla_workers")
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (is.finite(detected_cores) && workers > detected_cores) {
  log_message("Reducing --workers from ", workers, " to the ", detected_cores, " detected logical CPUs.")
  workers <- as.integer(detected_cores)
}
hla_workers <- min(hla_workers, workers)
if (.Platform$OS.type == "windows" && (workers > 1L || hla_workers > 1L)) {
  log_message("Fork-based dynamic scheduling is unavailable on Windows; using one worker.")
  workers <- 1L
  hla_workers <- 1L
}

# Prevent nested data.table/BLAS parallelism inside each process. The expensive
# parallel unit is modality x cell type, not an individual coloc.abf call.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
data.table::setDTthreads(1L)

traits_all <- read_traits(project_root, cfg)
traits <- traits_all[trait_id %in% c(
  cfg$primary_comparison$graves,
  cfg$primary_comparison$autoimmune_hypothyroidism
)]
cell_types <- read_cell_types(project_root, cfg)[qtl_available == TRUE]

modalities <- split_csv(args$modalities %||% "eQTL,caQTL")
if (any(!modalities %in% c("eQTL", "caQTL"))) stop("--modalities must contain eQTL and/or caQTL", call. = FALSE)
requested_cells <- split_csv(args$cell_types %||% "")
if (length(requested_cells)) {
  cell_types[, short_id := paste(level, cell_type, sep = ".")]
  unknown <- setdiff(requested_cells, cell_types$short_id)
  if (length(unknown)) stop("Unknown cell types: ", paste(unknown, collapse = ", "), call. = FALSE)
  cell_types <- cell_types[short_id %in% requested_cells]
}

locus_path <- file.path(project_root, cfg$project$outcome_dir, "loci", "primary_shared_divergent_loci.tsv")
if (!file.exists(locus_path)) stop("Run codes/02_define_loci.R first.", call. = FALSE)
loci <- data.table::fread(locus_path)
requested_loci <- split_csv(args$loci %||% "")
if (length(requested_loci)) loci <- loci[locus_id %in% requested_loci]
if (!nrow(loci)) stop("No loci selected.", call. = FALSE)

standardize_gwas_region <- function(x) {
  if (!nrow(x)) return(x)
  x[, `:=`(
    variant_id = canonical_variant_id(chrom, pos, ref, alt),
    gwas_beta = as.numeric(beta),
    gwas_se = as.numeric(sebeta),
    gwas_p = as.numeric(pval),
    gwas_eaf = as.numeric(af_alt)
  )]
  x[, gwas_maf := pmin(gwas_eaf, 1 - gwas_eaf)]
  x[
    is.finite(gwas_p) & gwas_p > 0 & gwas_p <= 1 &
      is.finite(gwas_maf) & gwas_maf >= as.numeric(cfg$coloc$minimum_gwas_maf),
    .(variant_id, gwas_beta, gwas_se, gwas_p, gwas_eaf, gwas_maf)
  ]
}

standardize_qtl_region <- function(x, modality) {
  if (!nrow(x)) return(data.table::data.table())
  if (identical(modality, "eQTL")) {
    ans <- x[, .(
      feature_id = phenotype_id,
      variant_id = MarkerID,
      qtl_af = as.numeric(AF_Allele2),
      qtl_beta = as.numeric(BETA),
      qtl_se = as.numeric(SE),
      qtl_p = as.numeric(p.value)
    )]
  } else {
    ans <- x[, .(
      feature_id = phenotype_id,
      variant_id = variant_id,
      qtl_af = as.numeric(af),
      qtl_beta = as.numeric(slope),
      qtl_se = as.numeric(slope_se),
      qtl_p = as.numeric(pval_nominal)
    )]
  }
  ans[, qtl_maf := pmin(qtl_af, 1 - qtl_af)]
  ans[is.finite(qtl_p) & qtl_p > 0 & qtl_p <= 1 & is.finite(qtl_maf) & qtl_maf > 0]
}

read_significant_features <- function(modality, level, cell_type) {
  path <- qtl_file(project_root, modality, "cis", level, cell_type)
  if (!file.exists(path)) stop("Missing QTL ACAT file: ", path, call. = FALSE)
  if (identical(modality, "eQTL")) {
    x <- data.table::fread(path, select = c("phenotype_id", "ACAT_q"), showProgress = FALSE)
    data.table::setnames(x, "ACAT_q", "feature_q")
  } else {
    x <- data.table::fread(path, select = c("phenotype_id", "qval"), showProgress = FALSE)
    data.table::setnames(x, "qval", "feature_q")
  }
  if (isTRUE(cfg$multiome$use_significant_qtl_features_only)) x <- x[feature_q < as.numeric(cfg$multiome$qtl_fdr)]
  x[order(feature_q)]
}

read_qtl_cs_qc <- function(modality, level, cell_type, chromosome, features) {
  path <- qtl_file(project_root, modality, "finemapping", level, cell_type, chromosome)
  if (!file.exists(path) || !length(features)) return(data.table::data.table())
  x <- data.table::fread(path, showProgress = FALSE)
  x <- x[region %in% features]
  if (!nrow(x)) return(data.table::data.table())
  x[, .(
    qtl_n_credible_sets = .N,
    qtl_max_cs_log10bf = max(cs_log10bf, na.rm = TRUE),
    qtl_has_pure_cs = any(!as.logical(low_purity))
  ), by = .(feature_id = region)]
}

extract_gwas_finemap <- function(x, top_variant) {
  if (!nrow(x)) return(list(pip = NA_real_, cs = NA_integer_, low_purity = NA))
  z <- x[rsid == top_variant]
  if (!nrow(z)) return(list(pip = NA_real_, cs = NA_integer_, low_purity = NA))
  z <- z[which.max(as.numeric(prob))]
  list(pip = as.numeric(z$prob), cs = as.integer(z$cs), low_purity = as.logical(z$low_purity))
}

run_one_coloc <- function(shared, trait_meta, locus, modality, level, cell_type, feature_id, feature_q, qtl_qc, gwas_finemap) {
  minimum_snps <- as.integer(cfg$multiome$minimum_shared_snps)
  if (nrow(shared) < minimum_snps) return(NULL)
  d1 <- list(
    pvalues = shared$gwas_p,
    MAF = shared$gwas_maf,
    N = as.numeric(trait_meta$n_cases + trait_meta$n_controls),
    s = as.numeric(trait_meta$n_cases / (trait_meta$n_cases + trait_meta$n_controls)),
    type = "cc",
    snp = shared$variant_id
  )
  d2 <- list(
    pvalues = shared$qtl_p,
    MAF = shared$qtl_maf,
    N = as.numeric(cfg$multiome$donor_n_default),
    type = "quant",
    snp = shared$variant_id
  )
  invisible(utils::capture.output(
    coloc_result <- suppressMessages(coloc::coloc.abf(
      dataset1 = d1,
      dataset2 = d2,
      p1 = as.numeric(cfg$coloc$p1),
      p2 = as.numeric(cfg$coloc$p2),
      p12 = as.numeric(cfg$coloc$p12)
    ))
  ))
  summary <- coloc_result$summary
  snp_results <- data.table::as.data.table(coloc_result$results)
  top <- snp_results[which.max(SNP.PP.H4)]
  top_shared <- shared[variant_id == top$snp][1]
  fm <- extract_gwas_finemap(gwas_finemap, top$snp)
  direction_product <- top_shared$gwas_beta * top_shared$qtl_beta
  direction <- if (!is.finite(direction_product) || direction_product == 0) {
    "undetermined"
  } else if (direction_product > 0) {
    "disease_risk_allele_increases_feature"
  } else {
    "disease_risk_allele_decreases_feature"
  }
  data.table::data.table(
    trait_id = trait_meta$trait_id,
    locus_id = locus$locus_id,
    chrom = locus$chrom,
    start = locus$start,
    end = locus$end,
    is_hla = locus$is_hla,
    modality = modality,
    level = level,
    cell_type = cell_type,
    feature_id = feature_id,
    feature_q = feature_q,
    n_shared_snps = nrow(shared),
    PP.H0.abf = as.numeric(summary[["PP.H0.abf"]]),
    PP.H1.abf = as.numeric(summary[["PP.H1.abf"]]),
    PP.H2.abf = as.numeric(summary[["PP.H2.abf"]]),
    PP.H3.abf = as.numeric(summary[["PP.H3.abf"]]),
    PP.H4.abf = as.numeric(summary[["PP.H4.abf"]]),
    top_h4_variant = top$snp,
    top_snp_pp_h4 = as.numeric(top$SNP.PP.H4),
    gwas_beta_at_top = top_shared$gwas_beta,
    qtl_beta_at_top = top_shared$qtl_beta,
    effect_direction = direction,
    gwas_susie_pip = fm$pip,
    gwas_susie_cs = fm$cs,
    gwas_susie_low_purity = fm$low_purity,
    qtl_n_credible_sets = qtl_qc$qtl_n_credible_sets %||% NA_integer_,
    qtl_max_cs_log10bf = qtl_qc$qtl_max_cs_log10bf %||% NA_real_,
    qtl_has_pure_cs = qtl_qc$qtl_has_pure_cs %||% NA,
    method = "coloc.abf_pvalue_screen",
    qtl_n_assumption = as.integer(cfg$multiome$donor_n_default),
    genome_build = cfg$project$genome_build
  )
}

empty_result <- data.table::data.table(
  trait_id = character(), locus_id = character(), chrom = integer(), start = integer(), end = integer(),
  is_hla = logical(), modality = character(), level = character(), cell_type = character(), feature_id = character(),
  feature_q = numeric(), n_shared_snps = integer(), PP.H0.abf = numeric(), PP.H1.abf = numeric(),
  PP.H2.abf = numeric(), PP.H3.abf = numeric(), PP.H4.abf = numeric(), top_h4_variant = character(),
  top_snp_pp_h4 = numeric(), gwas_beta_at_top = numeric(), qtl_beta_at_top = numeric(), effect_direction = character(),
  gwas_susie_pip = numeric(), gwas_susie_cs = integer(), gwas_susie_low_purity = logical(),
  qtl_n_credible_sets = integer(), qtl_max_cs_log10bf = numeric(), qtl_has_pure_cs = logical(),
  method = character(), qtl_n_assumption = integer(), genome_build = character()
)

empty_error <- data.table::data.table(
  trait_id = character(), locus_id = character(), modality = character(), level = character(),
  cell_type = character(), feature_id = character(), phase = character(), error_stage = character(),
  worker_pid = integer(), error = character()
)

atomic_write_tsv <- function(x, path, gzip = grepl("\\.gz$", path, ignore.case = TRUE)) {
  ensure_dir(dirname(path))
  suffix <- if (gzip) ".tmp.tsv.gz" else ".tmp.tsv"
  temporary_path <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = suffix)
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  # data.table::fwrite() can leave a truncated gzip footer when writing a
  # zero-row table on some filesystems.  Write the header through R's gzip
  # connection in that case so even deliberately empty results pass gzip -t.
  if (gzip && !nrow(x)) {
    connection <- gzfile(temporary_path, open = "wt")
    on.exit(try(close(connection), silent = TRUE), add = TRUE)
    writeLines(paste(names(x), collapse = "\t"), connection)
    close(connection)
  } else {
    data.table::fwrite(
      x, temporary_path, sep = "\t", quote = FALSE, na = "NA",
      compress = if (gzip) "gzip" else "none"
    )
  }
  if (!file.rename(temporary_path, path)) {
    stop("Atomic rename failed for output: ", path, call. = FALSE)
  }
  invisible(path)
}

build_gwas_cache <- function(phase_loci, phase) {
  cache <- vector("list", nrow(phase_loci))
  names(cache) <- phase_loci$locus_id
  for (locus_i in seq_len(nrow(phase_loci))) {
    locus <- phase_loci[locus_i]
    log_message("Preparing shared GWAS cache for ", phase, " locus ", locus_i, "/", nrow(phase_loci), ": ", locus$locus_id)
    gwas_regions <- lapply(seq_len(nrow(traits)), function(i) {
      standardize_gwas_region(
        scan_tabix_dt(gwas_path(project_root, traits$trait_id[i]), locus$chrom, locus$start, locus$end)
      )
    })
    names(gwas_regions) <- traits$trait_id
    gwas_finemap_regions <- lapply(seq_len(nrow(traits)), function(i) {
      path <- gwas_susie_path(project_root, traits$trait_id[i])
      if (!file.exists(path) || !file.exists(paste0(path, ".tbi"))) return(data.table::data.table())
      available_seqnames <- Rsamtools::seqnamesTabix(Rsamtools::TabixFile(path))
      chromosome <- sub("^chr", "", as.character(locus$chrom), ignore.case = TRUE)
      if (!any(c(chromosome, paste0("chr", chromosome)) %in% available_seqnames)) {
        return(data.table::data.table())
      }
      scan_tabix_dt(path, locus$chrom, locus$start, locus$end)
    })
    names(gwas_finemap_regions) <- traits$trait_id
    cache[[as.character(locus$locus_id)]] <- list(
      gwas = gwas_regions,
      gwas_finemap = gwas_finemap_regions
    )
  }
  cache
}

process_task_phase <- function(task, phase_loci, gwas_cache, phase, task_tmp_root) {
  data.table::setDTthreads(1L)
  modality <- as.character(task$modality[[1]])
  level <- as.character(task$level[[1]])
  cell_type <- as.character(task$cell_type[[1]])
  task_id <- as.character(task$task_id[[1]])
  worker_pid <- Sys.getpid()
  result_part <- file.path(task_tmp_root, paste0(task_id, ".", phase, ".results.tsv.gz"))
  error_part <- file.path(task_tmp_root, paste0(task_id, ".", phase, ".errors.tsv.gz"))
  feature_set <- read_significant_features(modality, level, cell_type)
  task_results <- list()
  task_errors <- list()

  log_message("Worker ", worker_pid, " started ", phase, " task ", task_id, " with ", nrow(phase_loci), " loci.")
  for (locus_i in seq_len(nrow(phase_loci))) {
    locus <- phase_loci[locus_i]
    locus_id <- as.character(locus$locus_id[[1]])
    nominal_path <- qtl_file(project_root, modality, "cis_nominal", level, cell_type)
    qtl_raw <- scan_tabix_dt(nominal_path, locus$chrom, locus$start, locus$end)
    if (!nrow(qtl_raw)) next
    qtl <- standardize_qtl_region(qtl_raw, modality)
    qtl <- qtl[feature_id %in% feature_set$phenotype_id]
    if (!nrow(qtl)) next
    present_features <- unique(qtl$feature_id)
    feature_meta <- feature_set[phenotype_id %in% present_features][order(feature_q)]
    feature_meta <- head(feature_meta, as.integer(cfg$coloc$maximum_features_per_locus))
    qtl <- qtl[feature_id %in% feature_meta$phenotype_id]
    cs_qc <- read_qtl_cs_qc(modality, level, cell_type, locus$chrom, feature_meta$phenotype_id)
    qtl_groups <- split(qtl, by = "feature_id", keep.by = TRUE)
    locus_cache <- gwas_cache[[locus_id]]

    for (feature_id in names(qtl_groups)) {
      qtl_feature <- qtl_groups[[feature_id]]
      feature_q <- feature_meta[phenotype_id == feature_id, feature_q][1]
      feature_value <- feature_id
      qtl_qc <- if (nrow(cs_qc)) cs_qc[feature_id == feature_value] else list()
      if (is.data.frame(qtl_qc) && !nrow(qtl_qc)) qtl_qc <- list()
      for (trait_i in seq_len(nrow(traits))) {
        trait_meta <- traits[trait_i]
        trait_id <- as.character(trait_meta$trait_id[[1]])
        shared <- merge(locus_cache$gwas[[trait_id]], qtl_feature, by = "variant_id", all = FALSE)
        attempt <- tryCatch(
          list(
            value = run_one_coloc(
              shared, trait_meta, locus, modality, level, cell_type, feature_id, feature_q,
              qtl_qc, locus_cache$gwas_finemap[[trait_id]]
            ),
            error = NULL
          ),
          error = function(e) list(value = NULL, error = conditionMessage(e))
        )
        if (!is.null(attempt$error)) {
          task_errors[[length(task_errors) + 1L]] <- data.table::data.table(
            trait_id = trait_id, locus_id = locus_id, modality = modality, level = level,
            cell_type = cell_type, feature_id = feature_id, phase = phase,
            error_stage = "coloc.abf", worker_pid = worker_pid, error = attempt$error
          )
        }
        one <- attempt$value
        if (!is.null(one) && one$PP.H4.abf >= as.numeric(cfg$coloc$pp_h4_report)) {
          task_results[[length(task_results) + 1L]] <- one
        }
      }
    }
    rm(qtl_raw, qtl, qtl_groups)
  }

  value <- if (length(task_results)) data.table::rbindlist(task_results, fill = TRUE) else data.table::copy(empty_result)
  if (nrow(value)) data.table::setorder(value, chrom, start, trait_id, feature_id)
  errors <- if (length(task_errors)) data.table::rbindlist(task_errors, fill = TRUE) else data.table::copy(empty_error)
  atomic_write_tsv(value, result_part, gzip = TRUE)
  atomic_write_tsv(errors, error_part, gzip = TRUE)
  log_message("Worker ", worker_pid, " finished ", phase, " task ", task_id, ": ", nrow(value), " reported coloc rows; ", nrow(errors), " errors.")
  data.table::data.table(
    task_id = task_id, phase = phase, ok = TRUE, result_part = result_part, error_part = error_part,
    n_results = nrow(value), n_errors = nrow(errors), worker_pid = worker_pid, error_message = NA_character_
  )
}

run_phase <- function(task_table, phase_loci, phase, worker_limit, task_tmp_root) {
  if (!nrow(phase_loci) || !nrow(task_table)) return(data.table::data.table())
  active_workers <- min(as.integer(worker_limit), nrow(task_table))
  log_message("Starting ", phase, " phase with ", active_workers, " worker(s), ", nrow(task_table), " modality-cell tasks, and ", nrow(phase_loci), " loci.")
  gwas_cache <- build_gwas_cache(phase_loci, phase)
  run_safely <- function(i) {
    task <- task_table[i]
    tryCatch(
      process_task_phase(task, phase_loci, gwas_cache, phase, task_tmp_root),
      error = function(e) data.table::data.table(
        task_id = as.character(task$task_id[[1]]), phase = phase, ok = FALSE,
        result_part = NA_character_, error_part = NA_character_, n_results = NA_integer_,
        n_errors = 1L, worker_pid = Sys.getpid(), error_message = conditionMessage(e)
      )
    )
  }
  indexes <- seq_len(nrow(task_table))
  status_list <- if (active_workers == 1L) {
    lapply(indexes, run_safely)
  } else {
    parallel::mclapply(
      indexes, run_safely, mc.cores = active_workers,
      mc.preschedule = FALSE, mc.set.seed = FALSE
    )
  }
  status <- data.table::rbindlist(status_list, fill = TRUE)
  rm(gwas_cache)
  invisible(gc())
  status
}

coloc_root <- file.path(project_root, cfg$project$outcome_dir, "coloc_r13")
ensure_dir(coloc_root)
task_tables <- lapply(modalities, function(modality_value) {
  result_dir <- file.path(coloc_root, modality_value)
  ensure_dir(result_dir)
  ans <- data.table::copy(cell_types)[, .(level, cell_type)]
  ans[, `:=`(
    modality = modality_value,
    task_id = paste(modality_value, level, cell_type, sep = "."),
    output_path = file.path(result_dir, paste0(level, ".", cell_type, ".tsv.gz"))
  )]
  ans
})
tasks <- data.table::rbindlist(task_tables, use.names = TRUE)
tasks[, nominal_path := mapply(
  function(modality, level, cell_type) qtl_file(project_root, modality, "cis_nominal", level, cell_type),
  modality, level, cell_type, USE.NAMES = FALSE
)]
missing_nominal <- tasks[!file.exists(nominal_path), nominal_path]
if (length(missing_nominal)) stop("Missing nominal QTL files: ", paste(missing_nominal, collapse = ", "), call. = FALSE)
tasks[, nominal_bytes := as.numeric(file.info(nominal_path)$size)]
data.table::setorder(tasks, -nominal_bytes, modality, level, cell_type)
force_outputs <- as_flag(args$force)
if (!force_outputs) {
  skipped <- tasks[file.exists(output_path), .N]
  tasks <- tasks[!file.exists(output_path)]
  if (skipped) log_message("Skipping ", skipped, " completed modality-cell outputs; use --force to replace them.")
}

if (!nrow(tasks)) {
  log_message("All requested modality-cell outputs already exist; nothing to run.")
} else {
  run_token <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_pid", Sys.getpid())
  task_tmp_root <- file.path(coloc_root, ".task_tmp", run_token)
  ensure_dir(task_tmp_root)

  non_hla_loci <- loci[is.na(is_hla) | !is_hla]
  hla_loci <- loci[!is.na(is_hla) & is_hla]
  phase_names <- character()
  phase_statuses <- list()
  if (nrow(non_hla_loci)) {
    phase_names <- c(phase_names, "non_hla")
    phase_statuses[["non_hla"]] <- run_phase(tasks, non_hla_loci, "non_hla", workers, task_tmp_root)
  }
  if (nrow(hla_loci)) {
    phase_names <- c(phase_names, "hla")
    phase_statuses[["hla"]] <- run_phase(tasks, hla_loci, "hla", hla_workers, task_tmp_root)
  }
  status <- data.table::rbindlist(phase_statuses, fill = TRUE)

  worker_error_parts <- status[ok & !is.na(error_part), error_part]
  error_tables <- lapply(worker_error_parts[file.exists(worker_error_parts)], function(path) {
    data.table::fread(path, showProgress = FALSE)
  })
  task_failures <- status[ok == FALSE, .(
    trait_id = NA_character_, locus_id = NA_character_,
    modality = sub("\\..*$", "", task_id), level = NA_character_, cell_type = NA_character_,
    feature_id = NA_character_, phase = phase, error_stage = "task",
    worker_pid = worker_pid, error = error_message
  )]
  if (nrow(task_failures)) error_tables[[length(error_tables) + 1L]] <- task_failures
  all_errors <- if (length(error_tables)) data.table::rbindlist(error_tables, fill = TRUE) else data.table::copy(empty_error)
  error_output <- file.path(project_root, cfg$project$output_dir, "qc", "coloc_r13_errors.tsv")
  atomic_write_tsv(all_errors, error_output, gzip = FALSE)

  finalized <- character()
  for (task_i in seq_len(nrow(tasks))) {
    task <- tasks[task_i]
    task_status <- status[task_id == task$task_id]
    complete <- nrow(task_status) == length(phase_names) &&
      setequal(task_status$phase, phase_names) && all(task_status$ok)
    if (!complete) next
    result_parts <- task_status$result_part
    value <- data.table::rbindlist(lapply(result_parts, data.table::fread, showProgress = FALSE), fill = TRUE)
    if (!nrow(value)) value <- data.table::copy(empty_result)
    if (nrow(value)) data.table::setorder(value, chrom, start, trait_id, feature_id)
    atomic_write_tsv(value, task$output_path[[1]], gzip = TRUE)
    finalized <- c(finalized, as.character(task$task_id[[1]]))
  }

  normalized_tmp <- normalizePath(task_tmp_root, mustWork = TRUE)
  normalized_root <- normalizePath(file.path(coloc_root, ".task_tmp"), mustWork = TRUE)
  if (startsWith(normalized_tmp, paste0(normalized_root, .Platform$file.sep))) {
    unlink(task_tmp_root, recursive = TRUE, force = TRUE)
  }

  failed_tasks <- setdiff(tasks$task_id, finalized)
  if (length(failed_tasks)) {
    stop(
      "Step 07 did not finalize ", length(failed_tasks), " task(s): ", paste(failed_tasks, collapse = ", "),
      ". Review ", error_output, " and rerun; completed outputs will be skipped.",
      call. = FALSE
    )
  }
  log_message(
    "R13 multiome colocalization screening finished with ", workers,
    " non-HLA worker(s) and ", hla_workers, " HLA worker(s). Finalized ", length(finalized), " outputs."
  )
}
