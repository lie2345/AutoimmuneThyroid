#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))

assert_packages(c("yaml", "data.table", "ggplot2", "scales"))
cfg <- read_analysis_config(project_root)
out_dir <- file.path(project_root, cfg$project$outcome_dir, "t07")
figure_dir <- file.path(project_root, cfg$project$figure_dir, "t07")
ensure_dir(figure_dir)

units <- data.table::fread(file.path(out_dir, "genetic_signal_five_classification.tsv"))
sensitivity <- data.table::fread(file.path(out_dir, "classification_equivalence_sensitivity.tsv"))
t04b <- data.table::fread(file.path(out_dir, "t04b_locus_t07_resolution.tsv"))

class_levels <- c(
  "shared-concordant", "shared-opposite", "Graves-predominant",
  "AIHT-predominant", "unresolved"
)
class_labels <- c(
  `shared-concordant` = "Shared–concordant",
  `shared-opposite` = "Shared–opposite",
  `Graves-predominant` = "Graves-predominant",
  `AIHT-predominant` = "AIHT-predominant",
  unresolved = "Unresolved"
)
class_colours <- c(
  `shared-concordant` = "#1B7837",
  `shared-opposite` = "#C51B7D",
  `Graves-predominant` = "#D95F02",
  `AIHT-predominant` = "#2166AC",
  unresolved = "#BDBDBD"
)

class_counts <- units[, .N, by = five_class_primary]
class_counts <- merge(
  data.table::data.table(five_class_primary = class_levels),
  class_counts,
  by = "five_class_primary", all.x = TRUE
)
class_counts[is.na(N), N := 0L]
class_counts[, five_class_primary := factor(five_class_primary, levels = class_levels)]

p1 <- ggplot2::ggplot(class_counts, ggplot2::aes(x = five_class_primary, y = N, fill = five_class_primary)) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.35, size = 4) +
  ggplot2::scale_x_discrete(labels = class_labels) +
  ggplot2::scale_fill_manual(values = class_colours, drop = FALSE) +
  ggplot2::scale_y_continuous(trans = "sqrt", breaks = c(0, 2, 6, 10, 50, 100, 200, 300)) +
  ggplot2::expand_limits(y = sqrt(max(class_counts$N)) + 1) +
  ggplot2::labs(
    title = "T07 primary five-class genetic signal model",
    subtitle = "304 independent signal units; y-axis uses a square-root transformation",
    x = NULL, y = "Number of signal units"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    plot.margin = ggplot2::margin(8, 12, 8, 12)
  )
ggplot2::ggsave(file.path(figure_dir, "t07_five_class_counts.png"), p1, width = 10.2, height = 5.8, dpi = 300)

h3_candidates <- units[unit_kind == "single_trait_signal" & distinct_model_supported == TRUE]
h3_candidates[, plot_label := paste0(
  signal_unit_id, " | ", data.table::fcoalesce(representative_nearest_genes, representative_reported_rsids),
  " | focal ", source_trait_label
)]
h3_candidates[, plot_label := factor(plot_label, levels = rev(plot_label[order(other_beta)]))]

main_bound <- log(as.numeric(cfg$t07_five_classification$equivalence_or_upper_primary))
strict_bound <- log(min(as.numeric(unlist(cfg$t07_five_classification$equivalence_or_upper_sensitivity))))
relaxed_bound <- log(max(as.numeric(unlist(cfg$t07_five_classification$equivalence_or_upper_sensitivity))))
p2 <- ggplot2::ggplot(h3_candidates, ggplot2::aes(
  y = plot_label, x = other_beta, xmin = other_ci90_lower, xmax = other_ci90_upper,
  colour = five_class_primary
)) +
  ggplot2::annotate("rect", xmin = -main_bound, xmax = main_bound, ymin = -Inf, ymax = Inf,
    fill = "#E5F5E0", alpha = 0.6) +
  ggplot2::geom_vline(xintercept = c(-relaxed_bound, relaxed_bound), linetype = "dotted", colour = "#636363") +
  ggplot2::geom_vline(xintercept = c(-main_bound, main_bound), linetype = "dashed", colour = "#238B45") +
  ggplot2::geom_vline(xintercept = c(-strict_bound, strict_bound), linetype = "dotdash", colour = "#006D2C") +
  ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  ggplot2::geom_errorbar(orientation = "y", width = 0.18, linewidth = 0.65) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_colour_manual(values = class_colours, labels = class_labels, drop = FALSE) +
  ggplot2::labs(
    title = "T07 equivalence gate among high-H3 signal units",
    subtitle = "Points and bars are the opposite-trait log-OR and 90% CI; green shading marks the primary OR 0.909–1.10 interval",
    x = "Opposite-trait log-OR at focal SuSiE lead variant", y = NULL, colour = "Primary class",
    caption = "Dot-dash: OR 1.05; dashed: OR 1.10; dotted: OR 1.15. Predominant labels also require a genome-wide case-difference test and direction consistency."
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0)
  )
ggplot2::ggsave(file.path(figure_dir, "t07_h3_equivalence_candidates.png"), p2, width = 11.2, height = 7.2, dpi = 300)

sensitivity_counts <- sensitivity[
  five_class %in% c("Graves-predominant", "AIHT-predominant"),
  .N,
  by = .(equivalence_or_upper, five_class)
]
sensitivity_counts <- merge(
  data.table::CJ(
    equivalence_or_upper = sort(unique(sensitivity$equivalence_or_upper)),
    five_class = c("Graves-predominant", "AIHT-predominant")
  ),
  sensitivity_counts,
  by = c("equivalence_or_upper", "five_class"), all.x = TRUE
)
sensitivity_counts[is.na(N), N := 0L]
p3 <- ggplot2::ggplot(sensitivity_counts, ggplot2::aes(
  x = factor(equivalence_or_upper), y = N, group = five_class, colour = five_class
)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.6, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = class_colours, labels = class_labels) +
  ggplot2::scale_y_continuous(breaks = 0:3, limits = c(0, 3)) +
  ggplot2::labs(
    title = "T07 disease-predominant classifications across equivalence bounds",
    x = "Symmetric equivalence upper OR", y = "Number of signal units", colour = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom")
ggplot2::ggsave(file.path(figure_dir, "t07_equivalence_sensitivity.png"), p3, width = 8.6, height = 5.4, dpi = 300)

t04b[, t07_resolution_class := data.table::fifelse(
  grepl("shared-opposite", t07_classes), "shared-opposite",
  data.table::fifelse(
    grepl("Graves-predominant", t07_classes), "Graves-predominant",
    data.table::fifelse(grepl("AIHT-predominant", t07_classes), "AIHT-predominant", "unresolved")
  )
)]
t04b_counts <- t04b[, .N, by = t07_resolution_class]
t04b_counts <- merge(
  data.table::data.table(t07_resolution_class = class_levels[-1L]),
  t04b_counts,
  by = "t07_resolution_class", all.x = TRUE
)
t04b_counts[is.na(N), N := 0L]
t04b_counts[, t07_resolution_class := factor(t07_resolution_class, levels = class_levels[-1L])]
p4 <- ggplot2::ggplot(t04b_counts, ggplot2::aes(x = t07_resolution_class, y = N, fill = t07_resolution_class)) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.35, size = 4) +
  ggplot2::scale_x_discrete(labels = class_labels) +
  ggplot2::scale_fill_manual(values = class_colours, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 16, 2)) +
  ggplot2::expand_limits(y = max(t04b_counts$N) + 1) +
  ggplot2::labs(
    title = "T07 resolution of the 16 T04B candidate loci",
    subtitle = "A locus is resolved only when at least one overlapping credible-set signal passes the complete five-class gate",
    x = NULL, y = "Number of T04B loci"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(figure_dir, "t07_t04b_locus_resolution.png"), p4, width = 9.6, height = 5.8, dpi = 300)

unmatched <- units[unit_kind == "single_trait_signal"]
gate_counts <- data.table::data.table(
  gate = c(
    "Unmatched source signals", "Source SuSiE QC passed", "Two-trait signal comparison available",
    "Conservative H3 posterior >=0.80", "No suggestive H4 competing model",
    "Opposite trait equivalent at OR 1.10", "Case difference P<5e-8",
    "Case-difference direction consistent", "Complete disease-predominant gate"
  ),
  n_signal_units = c(
    nrow(unmatched), sum(unmatched$source_signal_quality %in% TRUE),
    sum(unmatched$n_direct_signal_pairs > 0, na.rm = TRUE),
    sum(unmatched$distinct_model_supported %in% TRUE),
    sum(!(unmatched$shared_model_not_excluded %in% TRUE)),
    sum(unmatched$other_equivalent_or1_10 %in% TRUE),
    sum(unmatched$case_difference_gws %in% TRUE),
    sum(unmatched$case_difference_direction_consistent %in% TRUE),
    sum(unmatched$five_class_primary %in% c("Graves-predominant", "AIHT-predominant"))
  )
)

priority_units <- merge(
  units[
    five_class_primary != "unresolved" |
      signal_unit_id %in% sensitivity[five_class %in% c("Graves-predominant", "AIHT-predominant"), signal_unit_id]
  ],
  sensitivity[, .(
    sensitivity_classes = paste(paste0("OR", equivalence_or_upper, "=", five_class), collapse = ";")
  ), by = signal_unit_id],
  by = "signal_unit_id", all.x = TRUE
)

write_tsv(class_counts, file.path(out_dir, "classification_counts.tsv"))
write_tsv(sensitivity_counts, file.path(out_dir, "classification_sensitivity_counts.tsv"))
write_tsv(gate_counts, file.path(out_dir, "predominant_gate_counts.tsv"))
write_tsv(priority_units, file.path(out_dir, "priority_signal_units.tsv"))
write_tsv(t04b_counts, file.path(out_dir, "t04b_resolution_counts.tsv"))

stopifnot(
  nrow(units) == 304L,
  sum(class_counts$N) == 304L,
  units[five_class_primary == "shared-concordant", .N] == 6L,
  units[five_class_primary == "shared-opposite", .N] == 2L,
  units[five_class_primary == "Graves-predominant", .N] == 2L,
  units[five_class_primary == "AIHT-predominant", .N] == 0L,
  units[five_class_primary == "unresolved", .N] == 294L,
  t04b[t07_resolution_class == "shared-opposite", .N] == 2L,
  t04b[t07_resolution_class == "Graves-predominant", .N] == 2L,
  t04b[t07_resolution_class == "unresolved", .N] == 12L
)

log_message("T07 summary complete: figures and audit summaries written to ", out_dir)
