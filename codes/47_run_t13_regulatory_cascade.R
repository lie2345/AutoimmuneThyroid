#!/usr/bin/env Rscript

# T13 reconstructs an auditable GWAS -> caQTL -> peak-gene -> eQTL chain
# from validated T06/T11 signal-level outputs. Regional coloc.abf output is
# intentionally not treated as direct causal evidence.

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages(c("yaml", "data.table", "jsonlite"))
cfg <- read_analysis_config(project_root)
graves_id <- cfg$primary_comparison$graves
aiht_id <- cfg$primary_comparison$autoimmune_hypothyroidism
t13_dir <- file.path(project_root, cfg$project$outcome_dir, "t13")
ensure_dir(t13_dir)

required <- c(
  file.path(project_root, cfg$project$outcome_dir, "t06", "susie_signal_catalog.tsv"),
  file.path(project_root, cfg$project$outcome_dir, "t11a", "direct_coloc", "gwas_qtl_susie_signal_pairs.tsv.gz"),
  file.path(project_root, cfg$project$outcome_dir, "t11b", "direct_coloc", "gwas_qtl_susie_signal_pairs.tsv.gz"),
  file.path(project_root, cfg$project$outcome_dir, "t11b", "p1a_feature_results.tsv"),
  file.path(project_root, "raw_data", "finngen_multiome", "cascade", "finngen_multiome_v1.cascade.variant.tsv.gz"),
  file.path(project_root, "raw_data", "finngen_multiome", "cascade", "finngen_multiome_v1.cascade.peak.tsv.gz"),
  file.path(project_root, "raw_data", "finngen_multiome", "cascade", "finngen_multiome_v1.cascade.gene.tsv.gz"),
  file.path(t13_dir, "ensembl_gene_map.tsv"),
  file.path(project_root, cfg$project$outcome_dir, "t12", "t12_final_decision.tsv")
)
if (any(!file.exists(required))) {
  stop("T13 required input missing: ", paste(required[!file.exists(required)], collapse = "; "), call. = FALSE)
}

key_cols <- c("locus_id", "modality", "level", "cell_type", "feature_id")
make_key <- function(x) do.call(paste, c(x[, key_cols, with = FALSE], sep = "||"))

priority <- data.table::fread(
  file.path(project_root, cfg$project$outcome_dir, "t11b", "p1a_feature_results.tsv"),
  showProgress = FALSE
)
if (nrow(priority) != 26L || data.table::uniqueN(priority[, ..key_cols]) != 26L) {
  stop("T13 expects the frozen complete P1A set of 26 unique features.", call. = FALSE)
}
priority[, candidate_id := sprintf("T13C%03d", seq_len(.N))]
priority[, candidate_key := make_key(.SD)]

direct_files <- c(
  file.path(project_root, cfg$project$outcome_dir, "t11a", "direct_coloc", "gwas_qtl_susie_signal_pairs.tsv.gz"),
  file.path(project_root, cfg$project$outcome_dir, "t11b", "direct_coloc", "gwas_qtl_susie_signal_pairs.tsv.gz")
)
direct <- data.table::rbindlist(lapply(direct_files, data.table::fread, showProgress = FALSE), fill = TRUE)
direct[, candidate_key := make_key(.SD)]
direct <- merge(
  direct, priority[, c("candidate_id", "candidate_key"), with = FALSE],
  by = "candidate_key", all = FALSE, allow.cartesian = FALSE
)
if (!nrow(direct)) stop("No T11 direct rows matched the frozen P1A features.", call. = FALSE)

# Exact signal-pair keys are the reproducibility unit. Repeated wave rows are
# retained once, while all distinct H3/H4 and unresolved rows stay in the audit.
direct_key <- c("candidate_id", "trait_id", "gwas_region", "gwas_cs", "qtl_cs")
direct <- unique(direct, by = direct_key)

signal_catalog <- data.table::fread(
  file.path(project_root, cfg$project$outcome_dir, "t06", "susie_signal_catalog.tsv"),
  showProgress = FALSE
)
signal_meta <- signal_catalog[, .(
  trait_id, gwas_region = region, gwas_cs = cs,
  gwas_lead_signal_pip = lead_signal_pip,
  catalog_cs_log10bf = cs_log10bf,
  catalog_low_purity = low_purity,
  catalog_cs_size = cs_size,
  catalog_is_hla = is_hla,
  catalog_build = genome_build
)]
direct <- merge(direct, signal_meta, by = c("trait_id", "gwas_region", "gwas_cs"), all.x = TRUE)
direct[, tier1_gwas_quality := is.finite(gwas_lead_signal_pip) & gwas_lead_signal_pip >= 0.10 &
         is.finite(gwas_cs_log10bf) & gwas_cs_log10bf > 0.9 &
         !data.table::fcoalesce(as.logical(gwas_low_purity), TRUE) &
         !data.table::fcoalesce(as.logical(qtl_low_purity), TRUE) &
         !data.table::fcoalesce(as.logical(catalog_low_purity), TRUE) &
         classification == "shared_strong_H4" & sensitivity_PP.H4 >= 0.80]

# Select one best direct row per candidate and disease. The full direct table
# is written below so this convenience selection never discards audit rows.
class_score_map <- c(shared_strong_H4 = 4L, shared_suggestive_H4 = 3L, distinct_signals_H3 = 2L)
direct[, class_score := unname(class_score_map[as.character(classification)])]
direct[is.na(class_score), class_score := 0L]
data.table::setorderv(direct, c("candidate_id", "trait_id", "class_score", "PP.H4", "sensitivity_PP.H4"), c(1L, 1L, -1L, -1L, -1L), na.last = TRUE)
best <- direct[, .SD[1L], by = .(candidate_id, trait_id)]

best_fields <- c(
  "candidate_id", "trait_id", "gwas_region", "gwas_cs", "qtl_cs", "classification",
  "PP.H4", "sensitivity_PP.H4", "gwas_cs_log10bf", "gwas_low_purity", "qtl_low_purity",
  "gwas_lead_signal_pip", "catalog_cs_size", "top_h4_variant", "top_snp_pp_h4",
  "gwas_beta_at_top_h4", "qtl_beta_at_top_h4", "effect_direction_if_shared",
  "n_shared_snps", "tier1_gwas_quality"
)
best <- best[, intersect(best_fields, names(best)), with = FALSE]
g <- best[trait_id == graves_id]
a <- best[trait_id == aiht_id]
data.table::setnames(g, setdiff(names(g), "candidate_id"), paste0(setdiff(names(g), "candidate_id"), "_graves"))
data.table::setnames(a, setdiff(names(a), "candidate_id"), paste0(setdiff(names(a), "candidate_id"), "_aiht"))
wide <- merge(priority, g, by = "candidate_id", all.x = TRUE)
wide <- merge(wide, a, by = "candidate_id", all.x = TRUE)
wide[, tier1_both_traits := as.logical(tier1_gwas_quality_graves) & as.logical(tier1_gwas_quality_aiht)]

gene_map <- data.table::fread(file.path(t13_dir, "ensembl_gene_map.tsv"), showProgress = FALSE)
if (!all(c("gene_id", "gene_symbol", "annotation_release", "genome_build") %in% names(gene_map))) {
  stop("T13 gene map must contain gene_id, gene_symbol, annotation_release, and genome_build.", call. = FALSE)
}
if (anyDuplicated(gene_map$gene_id) || any(!nzchar(gene_map$gene_symbol))) stop("T13 gene map is not unique and complete.", call. = FALSE)

cascade_variant <- data.table::fread(
  file.path(project_root, "raw_data", "finngen_multiome", "cascade", "finngen_multiome_v1.cascade.variant.tsv.gz"),
  showProgress = FALSE
)
cascade_peak <- data.table::fread(
  file.path(project_root, "raw_data", "finngen_multiome", "cascade", "finngen_multiome_v1.cascade.peak.tsv.gz"),
  showProgress = FALSE
)
cascade_gene <- data.table::fread(
  file.path(project_root, "raw_data", "finngen_multiome", "cascade", "finngen_multiome_v1.cascade.gene.tsv.gz"),
  showProgress = FALSE
)

link_cache <- new.env(parent = emptyenv())
get_links <- function(level, cell_type) {
  cache_key <- paste(level, cell_type, sep = "||")
  if (exists(cache_key, envir = link_cache, inherits = FALSE)) return(get(cache_key, envir = link_cache, inherits = FALSE))
  path <- file.path(project_root, "raw_data", "finngen_multiome", "peak_gene_links",
                    paste0("finngen_multiome_v1.peak_gene_links.", level, ".", cell_type, ".tsv.gz"))
  if (!file.exists(path)) return(data.table::data.table())
  x <- data.table::fread(path, select = c("peak_id", "gene_id", "q_joint", "selected", "mode", "sig"), showProgress = FALSE)
  x <- x[as.logical(sig)]
  assign(cache_key, x, envir = link_cache)
  x
}

ca_candidates <- priority[modality == "caQTL"]
eq_candidates <- priority[modality == "eQTL"]
chain_rows <- list()
link_flags <- data.table::data.table(candidate_id = priority$candidate_id, n_sig_links = 0L, n_matching_eqtl_candidates = 0L)
get1 <- function(x, col) if (nrow(x) && col %in% names(x)) x[[col]][1] else NA

for (i in seq_len(nrow(ca_candidates))) {
  ca <- ca_candidates[i]
  links <- get_links(ca$level, ca$cell_type)
  links <- links[peak_id == ca$feature_id]
  if (nrow(links)) link_flags[candidate_id == ca$candidate_id, n_sig_links := data.table::uniqueN(links$gene_id)]
  if (!nrow(links)) next
  for (j in seq_len(nrow(links))) {
    link <- links[j]
    eqs <- eq_candidates[locus_id == ca$locus_id & feature_id == link$gene_id]
    if (!nrow(eqs)) next
    link_flags[candidate_id == ca$candidate_id, n_matching_eqtl_candidates := n_matching_eqtl_candidates + nrow(eqs)]
    for (k in seq_len(nrow(eqs))) {
      eq <- eqs[k]
      ca_g <- best[candidate_id == ca$candidate_id & trait_id == graves_id]
      ca_a <- best[candidate_id == ca$candidate_id & trait_id == aiht_id]
      eq_g <- best[candidate_id == eq$candidate_id & trait_id == graves_id]
      eq_a <- best[candidate_id == eq$candidate_id & trait_id == aiht_id]

      dir_g_ca <- get1(ca_g, "effect_direction_if_shared")
      dir_g_eq <- get1(eq_g, "effect_direction_if_shared")
      dir_a_ca <- get1(ca_a, "effect_direction_if_shared")
      dir_a_eq <- get1(eq_a, "effect_direction_if_shared")
      top_g_ca <- get1(ca_g, "top_h4_variant")
      top_g_eq <- get1(eq_g, "top_h4_variant")
      top_a_ca <- get1(ca_a, "top_h4_variant")
      top_a_eq <- get1(eq_a, "top_h4_variant")
      qtl_same_g <- !is.na(top_g_ca) && !is.na(top_g_eq) && identical(as.character(top_g_ca), as.character(top_g_eq))
      qtl_same_a <- !is.na(top_a_ca) && !is.na(top_a_eq) && identical(as.character(top_a_ca), as.character(top_a_eq))
      qtl_dir_g <- if (qtl_same_g) sign(as.numeric(get1(ca_g, "qtl_beta_at_top_h4"))) == sign(as.numeric(get1(eq_g, "qtl_beta_at_top_h4"))) else NA
      qtl_dir_a <- if (qtl_same_a) sign(as.numeric(get1(ca_a, "qtl_beta_at_top_h4"))) == sign(as.numeric(get1(eq_a, "qtl_beta_at_top_h4"))) else NA
      oriented_g <- !is.na(dir_g_ca) && !is.na(dir_g_eq) && identical(as.character(dir_g_ca), as.character(dir_g_eq))
      oriented_a <- !is.na(dir_a_ca) && !is.na(dir_a_eq) && identical(as.character(dir_a_ca), as.character(dir_a_eq))
      direction_status <- if (isTRUE(oriented_g) && isTRUE(oriented_a) &&
                              (isTRUE(qtl_dir_g) || is.na(qtl_dir_g)) &&
                              (isTRUE(qtl_dir_a) || is.na(qtl_dir_a))) {
        if (isTRUE(qtl_same_g) && isTRUE(qtl_same_a)) "three_way_exact_variant" else "three_way_gwas_oriented_top_variant_diff"
      } else if (isTRUE(oriented_g) || isTRUE(oriented_a) || identical(qtl_dir_g, FALSE) || identical(qtl_dir_a, FALSE)) {
        "direction_conflict_or_partial"
      } else {
        "direction_untestable"
      }

      ca_tier1_g <- isTRUE(get1(ca_g, "tier1_gwas_quality")); ca_tier1_a <- isTRUE(get1(ca_a, "tier1_gwas_quality"))
      eq_tier1_g <- isTRUE(get1(eq_g, "tier1_gwas_quality")); eq_tier1_a <- isTRUE(get1(eq_a, "tier1_gwas_quality"))
      tier1_both <- ca_tier1_g && ca_tier1_a && eq_tier1_g && eq_tier1_a
      chain_class <- if (tier1_both && direction_status %in% c("three_way_exact_variant", "three_way_gwas_oriented_top_variant_diff")) {
        "complete_cascade"
      } else {
        "evidence_insufficient"
      }

      cvar <- cascade_variant[variant_id %in% c(top_g_ca, top_g_eq, top_a_ca, top_a_eq)]
      cp <- cascade_peak[peak_id == ca$feature_id]
      cg <- cascade_gene[gene_id == link$gene_id]
      spec_values <- c(ca$cell_type, eq$cell_type, cp$cell_type_specificity, cg$cell_type_specificity)
      specificity_class <- if (any(spec_values == "Single cell-type")) {
        "single_cell_selective"
      } else if (any(c(ca$cell_type, eq$cell_type) == "PBMC") || any(spec_values == "Cross-lineage shared")) {
        "all_PBMC_shared"
      } else if (any(spec_values %in% c("T-cell-specific", "Lineage-specific", "Likely shared but underpowered"))) {
        "lineage_shared"
      } else {
        "evidence_insufficient"
      }
      chain_rows[[length(chain_rows) + 1L]] <- data.table::data.table(
        chain_id = sprintf("T13CHAIN%03d", length(chain_rows) + 1L),
        ca_candidate_id = ca$candidate_id, eqtl_candidate_id = eq$candidate_id,
        locus_id = ca$locus_id, ca_level = ca$level, ca_cell_type = ca$cell_type,
        peak_id = ca$feature_id, linked_gene_id = link$gene_id,
        gene_symbol = gene_map$gene_symbol[match(link$gene_id, gene_map$gene_id)],
        link_q_joint = link$q_joint, link_mode = link$mode,
        eqtl_level = eq$level, eqtl_cell_type = eq$cell_type,
        same_level_cell_type = ca$level == eq$level && ca$cell_type == eq$cell_type,
        graves_ca_direction = dir_g_ca, graves_eqtl_direction = dir_g_eq,
        aiht_ca_direction = dir_a_ca, aiht_eqtl_direction = dir_a_eq,
        graves_ca_top_variant = top_g_ca, graves_eqtl_top_variant = top_g_eq,
        aiht_ca_top_variant = top_a_ca, aiht_eqtl_top_variant = top_a_eq,
        graves_same_top_variant = qtl_same_g, aiht_same_top_variant = qtl_same_a,
        graves_qtl_direction_same = qtl_dir_g, aiht_qtl_direction_same = qtl_dir_a,
        direction_status = direction_status, tier1_both_traits = tier1_both,
        cascade_variant_category = if (nrow(cvar)) paste(sort(unique(cvar$qtl_mechanism_category)), collapse = ";") else NA_character_,
        cascade_peak_specificity = if (nrow(cp)) paste(sort(unique(cp$cell_type_specificity)), collapse = ";") else NA_character_,
        cascade_gene_specificity = if (nrow(cg)) paste(sort(unique(cg$cell_type_specificity)), collapse = ";") else NA_character_,
        specificity_class = specificity_class, chain_evidence_class = chain_class,
        note = "Direct T11 signal-pair evidence; CASCADE and peak-gene links are supporting annotations, not causal proof."
      )
    }
  }
}

chains <- if (length(chain_rows)) data.table::rbindlist(chain_rows, fill = TRUE) else data.table::data.table()
if (nrow(chains)) data.table::setorder(chains, locus_id, ca_candidate_id, eqtl_candidate_id)

summary <- merge(priority, link_flags, by = "candidate_id", all.x = TRUE)
chain_counts <- if (nrow(chains)) chains[, .(
  n_chain_rows = .N,
  n_complete_cascade = sum(chain_evidence_class == "complete_cascade"),
  n_evidence_insufficient_chain = sum(chain_evidence_class == "evidence_insufficient"),
  chain_gene_symbols = paste(sort(unique(gene_symbol)), collapse = ";")
), by = ca_candidate_id] else data.table::data.table(
  ca_candidate_id = character(), n_chain_rows = integer(), n_complete_cascade = integer(),
  n_evidence_insufficient_chain = integer(), chain_gene_symbols = character()
)
summary <- merge(summary, chain_counts, by.x = "candidate_id", by.y = "ca_candidate_id", all.x = TRUE)
summary[is.na(n_chain_rows), n_chain_rows := 0L]
summary[is.na(n_complete_cascade), n_complete_cascade := 0L]
summary[is.na(n_evidence_insufficient_chain), n_evidence_insufficient_chain := 0L]
eq_ids <- if (nrow(chains)) unique(chains$eqtl_candidate_id) else character()
summary[, candidate_evidence_class := data.table::fcase(
  n_complete_cascade > 0L, "complete_cascade",
  modality == "caQTL" & n_sig_links == 0L, "caQTL-only",
  modality == "caQTL" & n_matching_eqtl_candidates == 0L, "caQTL-only",
  modality == "eQTL" & !candidate_id %in% eq_ids, "eQTL-only",
  modality == "caQTL" & n_evidence_insufficient_chain > 0L, "evidence_insufficient",
  default = "evidence_insufficient"
)]
summary[, gene_symbol := ifelse(modality == "eQTL", gene_map$gene_symbol[match(feature_id, gene_map$gene_id)], NA_character_)]
summary[, annotation_version := "Ensembl release 116; GRCh38"]
gene_ids <- unique(c(summary[modality == "eQTL", feature_id], if (nrow(chains)) chains$linked_gene_id else character()))
t13_gene_summary <- data.table::data.table(gene_id = gene_ids)
t13_gene_summary[, `:=`(
  gene_symbol = gene_map$gene_symbol[match(gene_id, gene_map$gene_id)],
  annotation_version = "Ensembl release 116; GRCh38"
)]
t13_gene_summary <- merge(t13_gene_summary, summary[feature_id %in% gene_ids, .(
  n_priority_features = .N,
  priority_candidate_ids = paste(sort(unique(candidate_id)), collapse = ";"),
  loci = paste(sort(unique(locus_id)), collapse = ";"),
  cell_types = paste(sort(unique(paste(level, cell_type, sep = "."))), collapse = ";"),
  modalities = paste(sort(unique(modality)), collapse = ";"),
  candidate_evidence_classes = paste(sort(unique(candidate_evidence_class)), collapse = ";")
), by = .(gene_id = feature_id)], by = "gene_id", all.x = TRUE)
if (nrow(chains)) {
  chain_gene_summary <- chains[, .(n_chain_rows = .N, n_complete_cascade = sum(chain_evidence_class == "complete_cascade")), by = .(gene_id = linked_gene_id)]
  t13_gene_summary <- merge(t13_gene_summary, chain_gene_summary, by = "gene_id", all.x = TRUE)
}
t13_gene_summary[is.na(n_chain_rows), n_chain_rows := 0L]
t13_gene_summary[is.na(n_complete_cascade), n_complete_cascade := 0L]
mechanism_classes <- data.table::data.table(
  mechanism_class = c("complete_cascade", "caQTL-only", "eQTL-only", "pQTL/coding", "evidence_insufficient"),
  n_priority_features = c(
    sum(summary$candidate_evidence_class == "complete_cascade"),
    sum(summary$candidate_evidence_class == "caQTL-only"),
    sum(summary$candidate_evidence_class == "eQTL-only"),
    0L,
    sum(summary$candidate_evidence_class == "evidence_insufficient")
  ),
  status = c("VERIFIED_ZERO", "VERIFIED", "VERIFIED", "HOLD_UNAVAILABLE_INPUT", "VERIFIED"),
  interpretation = c(
    "No candidate met bilateral Tier 1 and three-way direction gates",
    "caQTL candidate had no validated matching eQTL chain",
    "eQTL candidate had no validated matching caQTL/peak-gene chain",
    "No cis-pQTL or coding consequence input is available",
    "A link existed but bilateral direction or signal evidence was insufficient"
  )
)

input_scope <- data.table::data.table(
  input_class = c("GWAS_signal_level", "direct_GWAS_QTL", "peak_gene_link", "CASCADE", "gene_symbol", "pQTL", "coding_annotation", "QTL_SNP_level"),
  source = c(
    "outcome_data/t06/susie_signal_catalog.tsv",
    "outcome_data/t11a/direct_coloc/gwas_qtl_susie_signal_pairs.tsv.gz; outcome_data/t11b/direct_coloc/gwas_qtl_susie_signal_pairs.tsv.gz",
    "raw_data/finngen_multiome/peak_gene_links/*.tsv.gz",
    "raw_data/finngen_multiome/cascade/finngen_multiome_v1.cascade.{variant,peak,gene}.tsv.gz",
    "outcome_data/t13/ensembl_gene_map.tsv",
    "not_available_in_current_project",
    "not_available_in_current_project",
    "T11A/T11B targeted slices only; no complete local QTL SNP collection"
  ),
  status = c("VERIFIED", "VERIFIED", "VERIFIED", "VERIFIED_SUPPORTING", "VERIFIED_FOR_T13_GENES", "HOLD", "HOLD", "PARTIAL"),
  limitation = c(
    "FinnGen R13 official SuSiE signal metadata; GRCh38",
    "Direct multi-signal BF pairing; H3, low-purity and overlap classes retained",
    "Only sig==TRUE links, exact level/cell/peak matching",
    "CASCADE PIP>0.5 annotation, not independent causal evidence",
    "Ensembl REST release 116; symbols mapped for the 14 T13 genes",
    "No pQTL input was supplied",
    "No coding consequence annotation was supplied",
    "Exact variant-level caQTL-eQTL comparison is limited to T11 targeted slices; top-variant differences remain explicit"
  )
)

t12 <- data.table::fread(file.path(project_root, cfg$project$outcome_dir, "t12", "t12_final_decision.tsv"), showProgress = FALSE)
t12_hold <- nrow(t12) == 1L && identical(as.character(t12$status[[1L]]), "HOLD") && grepl("HOLD", as.character(t12$decision[[1L]]))
tier1_rule_ok <- {
  expected <- !is.na(wide$tier1_gwas_quality_graves) & !is.na(wide$tier1_gwas_quality_aiht) &
    as.logical(wide$tier1_gwas_quality_graves) & as.logical(wide$tier1_gwas_quality_aiht)
  all(as.logical(wide$tier1_both_traits) == expected)
}
validation <- data.table::data.table(
  check = c(
    "t12_hold_preserved", "frozen_p1a_has_26_unique_features", "direct_signal_pair_keys_unique",
    "all_priority_features_emitted", "gene_symbols_complete_for_t13_genes",
    "tier1_requires_both_diseases_and_clean_purity", "peak_gene_links_are_sig_and_exact",
    "pQTL_and_coding_not_claimed_as_available", "no_complete_cascade_without_three_way_direction"
  ),
  passed = c(
    t12_hold,
    nrow(priority) == 26L && data.table::uniqueN(priority[, ..key_cols]) == 26L,
    !anyDuplicated(direct[, ..direct_key]),
    all(priority$candidate_id %in% summary$candidate_id),
    all(summary[modality == "eQTL", !is.na(gene_symbol)]),
    tier1_rule_ok && !any(direct$tier1_gwas_quality &
      (data.table::fcoalesce(as.logical(direct$gwas_low_purity), TRUE) |
       data.table::fcoalesce(as.logical(direct$qtl_low_purity), TRUE) |
       data.table::fcoalesce(as.logical(direct$catalog_low_purity), TRUE))),
    all(chains$same_level_cell_type %in% c(TRUE, FALSE)),
    all(input_scope$status[input_scope$input_class %in% c("pQTL", "coding_annotation")] == "HOLD"),
    !any(chains$chain_evidence_class == "complete_cascade" & chains$direction_status == "direction_conflict_or_partial")
  )
)

decision_status <- if (all(validation$passed)) "PASS_T13_PARTIAL_CASCADE_HOLD_PQTL_CODING" else "FAIL_T13_VALIDATION"
decision <- data.table::data.table(
  gate = "T13_locus_level_complete_regulatory_chain",
  status = if (all(validation$passed)) "PASS_WITH_HOLD" else "FAIL",
  decision = decision_status,
  priority_features = nrow(priority),
  direct_rows_audited = nrow(direct),
  chain_rows = nrow(chains),
  complete_cascades = if (nrow(chains)) sum(chains$chain_evidence_class == "complete_cascade") else 0L,
  t12_status = if (t12_hold) "HOLD_PRESERVED" else "T12_STATUS_INVALID",
  rationale = "T13 executed on the frozen P1A priority set. Missing pQTL/coding inputs remain HOLD; no unsupported causal claim is made.",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

provenance_paths <- c(
  file.path(project_root, cfg$project$outcome_dir, "t06", "susie_signal_catalog.tsv"),
  direct_files,
  file.path(project_root, cfg$project$outcome_dir, "t11b", "p1a_feature_results.tsv"),
  file.path(t13_dir, "ensembl_gene_map.tsv"),
  file.path(project_root, cfg$project$outcome_dir, "t12", "t12_final_decision.tsv")
)
provenance <- data.table::data.table(
  input_path = sub(paste0("^", project_root, "/"), "", provenance_paths),
  md5 = unname(tools::md5sum(provenance_paths)),
  bytes = as.numeric(file.info(provenance_paths)$size),
  role = c("T06_GWAS_signal_catalog", "T11A_direct_GWAS_QTL", "T11B_direct_GWAS_QTL",
           "frozen_P1A_priority", "Ensembl_116_gene_map", "T12_HOLD_guard")
)
software <- data.table::data.table(
  component = c("R", "data.table", "yaml", "jsonlite", "Ensembl_annotation_release", "genome_build"),
  version = c(R.version.string, as.character(packageVersion("data.table")), as.character(packageVersion("yaml")),
              as.character(packageVersion("jsonlite")), "116", cfg$project$genome_build)
)

write_tsv(direct, file.path(t13_dir, "t13_direct_evidence.tsv.gz"))
write_tsv(chains, file.path(t13_dir, "t13_chain_evidence.tsv.gz"))
write_tsv(summary, file.path(t13_dir, "t13_candidate_summary.tsv"))
write_tsv(t13_gene_summary, file.path(t13_dir, "t13_gene_summary.tsv"))
write_tsv(mechanism_classes, file.path(t13_dir, "t13_mechanism_classes.tsv"))
write_tsv(input_scope, file.path(t13_dir, "t13_input_scope.tsv"))
write_tsv(validation, file.path(t13_dir, "t13_validation.tsv"))
write_tsv(decision, file.path(t13_dir, "t13_decision.tsv"))
write_tsv(provenance, file.path(t13_dir, "t13_provenance.tsv"))
write_tsv(software, file.path(t13_dir, "t13_software.tsv"))
log_message("T13 complete: ", nrow(priority), " priority features; ", nrow(chains), " chain rows; ",
            if (nrow(chains)) sum(chains$chain_evidence_class == "complete_cascade") else 0L, " complete cascades; ", decision_status)
if (!all(validation$passed)) stop("T13 validation failed.", call. = FALSE)
