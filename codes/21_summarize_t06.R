#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "ggplot2", "Rsamtools", "GenomicRanges", "IRanges"))
cfg <- read_analysis_config(project_root)
graves <- cfg$primary_comparison$graves
aiht <- cfg$primary_comparison$autoimmune_hypothyroidism
out_dir <- file.path(project_root, cfg$project$outcome_dir, "t06")
figure_dir <- file.path(project_root, cfg$project$figure_dir, "t06")
ensure_dir(figure_dir)

coloc <- data.table::fread(file.path(out_dir, "gwas_gwas_susie_coloc.tsv"))
architecture <- data.table::fread(file.path(out_dir, "gwas_gwas_region_architecture.tsv"))
coverage <- data.table::fread(file.path(out_dir, "t04b_locus_signal_coverage.tsv"))

parse_variant <- function(x) {
  pieces <- regmatches(x, regexec("^chr([^_]+)_([0-9]+)_([^_]+)_([^_]+)$", x))
  valid <- lengths(pieces) == 5L
  data.table::data.table(
    variant_id = x,
    chrom = ifelse(valid, vapply(pieces, function(z) z[[2L]] %||% NA_character_, character(1)), NA_character_),
    pos = ifelse(valid, as.integer(vapply(pieces, function(z) z[[3L]] %||% NA_character_, character(1))), NA_integer_),
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
    seqname <- c(paste0("chr", v$chrom), v$chrom)
    seqname <- seqname[seqname %in% available][1L]
    if (is.na(seqname)) {
      return(data.table::data.table(
        variant_id = v$variant_id, nearest_genes = NA_character_, reported_rsids = NA_character_
      ))
    }
    param <- GenomicRanges::GRanges(
      seqnames = seqname,
      ranges = IRanges::IRanges(start = v$pos, end = v$pos)
    )
    lines <- suppressWarnings(Rsamtools::scanTabix(tabix_file, param = param)[[1L]])
    if (!length(lines)) {
      return(data.table::data.table(
        variant_id = v$variant_id, nearest_genes = NA_character_, reported_rsids = NA_character_
      ))
    }
    x <- data.table::fread(
      text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE, showProgress = FALSE
    )
    data.table::setnames(x, header)
    x <- x[pos == v$pos & toupper(ref) == toupper(v$ref) & toupper(alt) == toupper(v$alt)]
    if (!nrow(x)) {
      return(data.table::data.table(
        variant_id = v$variant_id, nearest_genes = NA_character_, reported_rsids = NA_character_
      ))
    }
    data.table::data.table(
      variant_id = v$variant_id,
      nearest_genes = x$nearest_genes[[1L]],
      reported_rsids = x$rsids[[1L]]
    )
  })
  data.table::rbindlist(ans, fill = TRUE)
}

g_ann <- annotate_variants(unique(c(coloc$graves_hit, coloc$top_h4_variant)), graves)
a_ann <- annotate_variants(unique(c(coloc$aiht_hit, coloc$top_h4_variant)), aiht)

annotated <- merge(
  coloc,
  g_ann[, .(
    graves_hit = variant_id,
    graves_hit_nearest_genes = nearest_genes,
    graves_hit_reported_rsids = reported_rsids
  )],
  by = "graves_hit", all.x = TRUE
)
annotated <- merge(
  annotated,
  a_ann[, .(
    aiht_hit = variant_id,
    aiht_hit_nearest_genes = nearest_genes,
    aiht_hit_reported_rsids = reported_rsids
  )],
  by = "aiht_hit", all.x = TRUE
)
annotated <- merge(
  annotated,
  g_ann[, .(
    top_h4_variant = variant_id,
    top_h4_nearest_genes = nearest_genes,
    top_h4_reported_rsids = reported_rsids
  )],
  by = "top_h4_variant", all.x = TRUE
)
data.table::setcolorder(annotated, c(
  "pair_id", "chromosome", "graves_region", "aiht_region", "graves_cs", "aiht_cs",
  "classification", "PP.H3", "PP.H4", "sensitivity_PP.H3", "sensitivity_PP.H4",
  "graves_hit", "graves_hit_nearest_genes", "aiht_hit", "aiht_hit_nearest_genes",
  "top_h4_variant", "top_h4_nearest_genes", "effect_direction_if_shared"
))
data.table::setorder(annotated, chromosome, graves_region, graves_cs, aiht_cs)
write_tsv(annotated, file.path(out_dir, "gwas_gwas_susie_coloc_annotated.tsv"))
write_tsv(
  annotated[classification %in% c("shared_strong_H4", "shared_suggestive_H4", "distinct_signals_H3")],
  file.path(out_dir, "gwas_gwas_high_confidence_signal_pairs.tsv")
)

class_labels <- c(
  shared_strong_H4 = "Strong shared signal (H4)",
  shared_suggestive_H4 = "Suggestive shared signal (H4)",
  distinct_signals_H3 = "Distinct signals (H3)",
  unresolved_posterior = "Unresolved posterior",
  unresolved_insufficient_overlap = "Insufficient CS overlap",
  unresolved_weak_credible_set = "Weak credible set",
  unresolved_low_purity = "Low-purity credible set"
)
class_colours <- c(
  shared_strong_H4 = "#1B7837", shared_suggestive_H4 = "#5AAE61",
  distinct_signals_H3 = "#762A83", unresolved_posterior = "#4575B4",
  unresolved_insufficient_overlap = "#D8B365", unresolved_weak_credible_set = "#E08214",
  unresolved_low_purity = "#BDBDBD"
)
plot_data <- data.table::copy(annotated)
plot_data[, classification := factor(classification, levels = names(class_labels))]
p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PP.H3, y = PP.H4, colour = classification)) +
  ggplot2::geom_point(size = 2.7, alpha = 0.85) +
  ggplot2::geom_vline(xintercept = 0.8, linetype = 2, colour = "#777777") +
  ggplot2::geom_hline(yintercept = 0.8, linetype = 2, colour = "#777777") +
  ggplot2::scale_colour_manual(values = class_colours, labels = class_labels, drop = TRUE) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "T06 signal-level Graves–AIHT colocalization",
    subtitle = "FinnGen R13 official SuSiE signals; dashed lines mark PP = 0.80",
    x = "PP.H3 (distinct causal variants)", y = "PP.H4 (shared causal variant)", colour = "Classification"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(figure_dir, "t06_h3_h4_signal_pairs.png"), p1, width = 10.5, height = 6.2, dpi = 300)

arch_labels <- c(
  mixed_shared_and_distinct = "Mixed shared\nand distinct",
  shared_only = "Shared only",
  distinct_only = "Distinct only",
  suggestive_shared = "Suggestive\nshared",
  unresolved = "Unresolved"
)
arch_plot <- architecture[, .N, by = architecture]
arch_plot[, architecture := factor(architecture, levels = names(arch_labels))]
p2 <- ggplot2::ggplot(arch_plot, ggplot2::aes(x = architecture, y = N, fill = architecture)) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.35, size = 4) +
  ggplot2::scale_x_discrete(labels = arch_labels) +
  ggplot2::scale_fill_manual(values = c(
    mixed_shared_and_distinct = "#2166AC", shared_only = "#1B7837",
    distinct_only = "#762A83", suggestive_shared = "#5AAE61", unresolved = "#BDBDBD"
  )) +
  ggplot2::labs(
    title = "T06 architecture of overlapping fine-mapped regions",
    x = NULL, y = "Number of region pairs"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    plot.margin = ggplot2::margin(8, 12, 8, 12)
  ) +
  ggplot2::expand_limits(y = max(arch_plot$N) + 1)
ggplot2::ggsave(file.path(figure_dir, "t06_region_architecture.png"), p2, width = 9.5, height = 5.8, dpi = 300)

coverage_plot <- coverage[, .N, by = signal_architecture]
coverage_plot[, signal_architecture := factor(signal_architecture, levels = c(
  "mixed_shared_and_distinct", "shared_only", "distinct_only", "suggestive_shared",
  "unresolved", "unresolved_one_trait_finemapped", "unresolved_no_official_finemap_region"
))]
p3 <- ggplot2::ggplot(coverage_plot, ggplot2::aes(x = signal_architecture, y = N, fill = signal_architecture)) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.35, size = 4) +
  ggplot2::scale_x_discrete(labels = c(
    mixed_shared_and_distinct = "Mixed",
    shared_only = "Shared only",
    distinct_only = "Distinct only",
    suggestive_shared = "Suggestive shared",
    unresolved = "Unresolved",
    unresolved_one_trait_finemapped = "Only one trait\nfine-mapped",
    unresolved_no_official_finemap_region = "Neither trait\nfine-mapped"
  )) +
  ggplot2::scale_fill_manual(values = c(
    mixed_shared_and_distinct = "#2166AC", shared_only = "#1B7837",
    distinct_only = "#762A83", suggestive_shared = "#5AAE61",
    unresolved = "#737373", unresolved_one_trait_finemapped = "#BDBDBD",
    unresolved_no_official_finemap_region = "#252525"
  ), drop = TRUE) +
  ggplot2::labs(
    title = "T06 signal-level coverage of 16 T04B candidate loci",
    x = NULL, y = "Number of T04B loci"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    plot.margin = ggplot2::margin(8, 12, 8, 12)
  ) +
  ggplot2::expand_limits(y = max(coverage_plot$N) + 1)
ggplot2::ggsave(file.path(figure_dir, "t06_t04b_signal_coverage.png"), p3, width = 10.2, height = 5.8, dpi = 300)

log_message(
  "T06 annotation complete: ",
  annotated[classification == "shared_strong_H4", .N], " strong H4; ",
  annotated[classification == "distinct_signals_H3", .N], " distinct H3 signal pairs"
)
