options(stringsAsFactors = FALSE, scipen = 999)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". See documentation/analysis/PIPELINE.md.",
      call. = FALSE
    )
  }
}

script_project_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    return(normalizePath(file.path(dirname(script), ".."), mustWork = TRUE))
  }
  current <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "config", "analysis.yml"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Cannot locate project root containing config/analysis.yml", call. = FALSE)
}

read_analysis_config <- function(project_root) {
  assert_packages("yaml")
  yaml::read_yaml(file.path(project_root, "config", "analysis.yml"))
}

project_path <- function(project_root, relative_path) {
  if (is.null(relative_path) || !nzchar(relative_path)) return("")
  if (grepl("^/", relative_path)) return(relative_path)
  file.path(project_root, relative_path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create directory: ", path, call. = FALSE)
  }
  invisible(path)
}

read_traits <- function(project_root, cfg) {
  assert_packages("data.table")
  data.table::fread(project_path(project_root, cfg$project$traits_file), sep = "\t")
}

read_cell_types <- function(project_root, cfg) {
  assert_packages("data.table")
  data.table::fread(project_path(project_root, cfg$project$cell_types_file), sep = "\t")
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    key <- sub("^--", "", token)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

log_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

canonical_variant_id <- function(chrom, pos, ref, alt, chr_prefix = TRUE) {
  chrom <- sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  prefix <- if (chr_prefix) "chr" else ""
  paste0(prefix, chrom, "_", pos, "_", toupper(ref), "_", toupper(alt))
}

first_rsid <- function(x) {
  x <- sub("[,;].*$", "", as.character(x))
  x[!grepl("^rs[0-9]+$", x)] <- NA_character_
  x
}

is_palindromic <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

is_hla_position <- function(chrom, pos, cfg) {
  as.integer(sub("^chr", "", chrom)) == as.integer(cfg$hla$chromosome) &
    pos >= as.numeric(cfg$hla$start) & pos <= as.numeric(cfg$hla$end)
}

read_gzip_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line)) return(character())
  sub("^#", "", strsplit(line, "\t", fixed = TRUE)[[1]])
}

scan_tabix_dt <- function(path, chrom, start, end) {
  assert_packages(c("data.table", "Rsamtools", "GenomicRanges", "IRanges"))
  if (!file.exists(path)) stop("Missing tabix file: ", path, call. = FALSE)
  if (!file.exists(paste0(path, ".tbi"))) stop("Missing tabix index: ", path, ".tbi", call. = FALSE)
  tabix_file <- Rsamtools::TabixFile(path)
  available_seqnames <- Rsamtools::seqnamesTabix(tabix_file)
  chr_core <- sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  chr_candidates <- c(paste0("chr", chr_core), chr_core)
  chr <- chr_candidates[chr_candidates %in% available_seqnames][1]
  if (is.na(chr)) {
    stop(
      "Chromosome ", chrom, " is absent from tabix index ", path,
      ". Available examples: ", paste(head(available_seqnames), collapse = ", "),
      call. = FALSE
    )
  }
  region <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = max(1, as.integer(start)), end = as.integer(end))
  )
  lines <- tryCatch(
    Rsamtools::scanTabix(tabix_file, param = region)[[1]],
    error = function(e) stop("Tabix query failed for ", path, " at ", chr, ":", start, "-", end, ": ", conditionMessage(e), call. = FALSE)
  )
  header <- read_gzip_header(path)
  if (!length(lines)) return(data.table::as.data.table(setNames(replicate(length(header), logical(0), simplify = FALSE), header)))
  ans <- data.table::fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE, showProgress = FALSE)
  if (ncol(ans) != length(header)) {
    stop("Header/data column mismatch in ", path, call. = FALSE)
  }
  data.table::setnames(ans, header)
  ans
}

gwas_path <- function(project_root, trait_id) {
  file.path(project_root, "raw_data", "finngen_r13", "summary_stats", paste0("finngen_R13_", trait_id, ".gz"))
}

gwas_susie_path <- function(project_root, trait_id) {
  file.path(project_root, "raw_data", "finngen_r13", "finemap", "susie", paste0("finngen_R13_", trait_id, ".SUSIE.snp.bgz"))
}

qtl_file <- function(project_root, modality, product, level, cell_type, chromosome = NULL) {
  base <- file.path(project_root, "raw_data", "finngen_multiome", modality, product)
  stem <- paste0("finngen_multiome_v1.", modality, ".", level, ".", cell_type)
  if (identical(product, "cis")) return(file.path(base, paste0(stem, ".cis_acat.tsv.gz")))
  if (identical(product, "cis_nominal")) return(file.path(base, paste0(stem, ".cis_nominal.tsv.gz")))
  if (identical(product, "finemapping")) {
    if (is.null(chromosome)) stop("chromosome is required for fine-mapping files", call. = FALSE)
    chr <- sub("^chr", "", chromosome)
    return(file.path(base, paste0(stem, ".chr", chr, ".SUSIE.cred.tsv.gz")))
  }
  stop("Unsupported QTL product: ", product, call. = FALSE)
}

write_tsv <- function(x, path) {
  assert_packages("data.table")
  ensure_dir(dirname(path))
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}
