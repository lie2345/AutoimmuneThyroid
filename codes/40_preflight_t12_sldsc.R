#!/usr/bin/env Rscript

project_root <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
})
source(file.path(project_root, "codes", "lib", "common.R"))
assert_packages(c("yaml", "data.table"))
cfg <- read_analysis_config(project_root)

out_dir <- file.path(project_root, cfg$project$outcome_dir, "t12")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
peaks_path <- file.path(
  project_root, cfg$project$raw_dir, "finngen_multiome", "atac_peaks",
  "finngen_multiome_v1.consensus_peaks.bed.gz"
)
if (!file.exists(peaks_path)) stop("Missing FinnGen Multiome consensus peak BED.", call. = FALSE)

peaks <- data.table::fread(
  peaks_path, header = FALSE,
  col.names = c("chrom", "start0", "end", "source_peak_names", "read_count", "strand"),
  showProgress = FALSE
)
if (nrow(peaks) != 298224L) stop("Unexpected consensus peak count: ", nrow(peaks), call. = FALSE)
peaks[, source_tokens := strsplit(source_peak_names, ",", fixed = TRUE)]
peaks[, source_cell_types := vapply(source_tokens, function(tokens) {
  paste(sort(unique(sub("_peak_.*$", "", tokens))), collapse = ";")
}, character(1))]
peaks[, n_source_cell_types := lengths(strsplit(source_cell_types, ";", fixed = TRUE))]
peaks[, cell_type := sub("^predicted[.]celltype[.]l2[.]", "", source_cell_types)]
peaks[, width := end - start0]
peaks[, overlaps_mhc := chrom == "chr6" & start0 < as.integer(cfg$hla$end) & end > as.integer(cfg$hla$start)]

lineage_map <- c(
  B_intermediate = "B", B_memory = "B", B_naive = "B", Plasmablast = "B",
  CD4_CTL = "T", CD4_Naive = "T", CD4_Proliferating = "T", CD4_TCM = "T",
  CD4_TEM = "T", CD8_Naive = "T", CD8_TCM = "T", CD8_TEM = "T",
  MAIT = "T", Treg = "T", dnT = "T", gdT = "T",
  CD14_Mono = "Mono", CD16_Mono = "Mono",
  cDC1 = "DC", cDC2 = "DC", pDC = "DC",
  NK = "NK", NK_CD56bright = "NK", NK_Proliferating = "NK",
  HSPC = "Other", ILC = "Other", Platelet = "Other", Eryth = "Other"
)
inventory <- peaks[, .(
  n_peaks_total = .N,
  n_peaks_mhc = sum(overlaps_mhc),
  n_peaks_non_mhc = sum(!overlaps_mhc),
  bp_total = sum(width),
  bp_non_mhc = sum(width[!overlaps_mhc])
), by = cell_type]
inventory[, lineage := unname(lineage_map[cell_type])]
data.table::setcolorder(inventory, c(
  "cell_type", "lineage", "n_peaks_total", "n_peaks_mhc", "n_peaks_non_mhc", "bp_total", "bp_non_mhc"
))
data.table::setorder(inventory, lineage, cell_type)
write_tsv(inventory, file.path(out_dir, "peak_source_inventory.tsv"))

peak_checks <- data.table::data.table(
  metric = c(
    "consensus_peaks", "source_cell_types", "multi_source_cell_type_peaks",
    "mhc_overlapping_peaks", "non_mhc_peaks", "peak_genome_build"
  ),
  value = c(
    nrow(peaks), nrow(inventory), sum(peaks$n_source_cell_types != 1L),
    sum(peaks$overlaps_mhc), sum(!peaks$overlaps_mhc), cfg$project$genome_build
  )
)
write_tsv(peak_checks, file.path(out_dir, "peak_source_summary.tsv"))

graves <- cfg$primary_comparison$graves
aiht <- cfg$primary_comparison$autoimmune_hypothyroidism
phenotypes <- data.table::data.table(
  requested_analysis = c(
    "Graves_case_control", "AIHT_case_control", "AITD_combined_case_control", "Graves_vs_AIHT_direct_case_case"
  ),
  usable_input = c(
    file.exists(file.path(project_root, "outcome_data", "t05", "ldsc", "munged", paste0(graves, ".sumstats.gz"))),
    file.exists(file.path(project_root, "outcome_data", "t05", "ldsc", "munged", paste0(aiht, ".sumstats.gz"))),
    FALSE, FALSE
  ),
  input_path = c(
    file.path("outcome_data", "t05", "ldsc", "munged", paste0(graves, ".sumstats.gz")),
    file.path("outcome_data", "t05", "ldsc", "munged", paste0(aiht, ".sumstats.gz")),
    NA_character_, NA_character_
  ),
  decision = c(
    "eligible_for_T12_S_LDSC", "eligible_for_T12_S_LDSC",
    "missing_union_case_GWAS_not_derivable_from_current_summary_statistics",
    "missing_direct_case_case_GWAS_CCGWAS_is_not_a_valid_substitute_for_S_LDSC"
  )
)
write_tsv(phenotypes, file.path(out_dir, "phenotype_scope.tsv"))

baseline_dir <- file.path(project_root, "raw_data", "reference_panels", "ldsc", "baselineLD_v2.2")
plink_dir <- file.path(project_root, "raw_data", "reference_panels", "ldsc", "1000G_EUR_Phase3_plink")
frq_dir <- file.path(project_root, "raw_data", "reference_panels", "ldsc", "1000G_Phase3_frq")
weight_prefix <- file.path(project_root, cfg$ldsc$regression_weight_prefix)
hm3_path <- file.path(project_root, cfg$ldsc$hapmap3_snplist)
chain_candidates <- c(
  file.path(project_root, "raw_data", "reference_panels", "liftover", "hg38ToHg19.over.chain.gz"),
  file.path(project_root, "raw_data", "reference_panels", "liftover", "hg38ToHg19.over.chain")
)

resource_gate <- data.table::data.table(
  resource = c(
    "baselineLD_v2.2_22_chromosomes", "1000G_EUR_Phase3_PLINK_22_chromosomes",
    "1000G_Phase3_frequency_22_chromosomes", "HapMap3_regression_SNP_list",
    "HapMap3_no_MHC_regression_weights", "GRCh38_to_GRCh37_chain"
  ),
  expected_build = c("GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh38_to_GRCh37"),
  observed_files = c(
    sum(file.exists(file.path(baseline_dir, paste0("baselineLD.", 1:22, ".l2.ldscore.gz")))),
    sum(file.exists(file.path(plink_dir, paste0("1000G.EUR.QC.", 1:22, ".bim")))),
    sum(file.exists(file.path(frq_dir, paste0("1000G.EUR.QC.", 1:22, ".frq")))),
    as.integer(file.exists(hm3_path)),
    sum(file.exists(paste0(weight_prefix, 1:22, ".l2.ldscore.gz"))),
    sum(file.exists(chain_candidates))
  ),
  expected_files = c(22L, 22L, 22L, 1L, 22L, 1L)
)
resource_gate[, ready := observed_files == expected_files]
write_tsv(resource_gate, file.path(out_dir, "reference_resource_gate.tsv"))

download_dir <- file.path(project_root, "raw_data", "reference_panels", "ldsc", "t12_downloads")
provenance <- data.table::data.table(
  resource = c("baselineLD_v2.2_archive", "1000G_Phase3_PLINK_archive", "1000G_Phase3_frequency_archive"),
  source = c(
    "https://zenodo.org/records/7768714/files/1000G_Phase3_baselineLD_v2.2_ldscores.tgz",
    "https://zenodo.org/records/7768714/files/1000G_Phase3_plinkfiles.tgz",
    "https://zenodo.org/records/7768714/files/1000G_Phase3_frq.tgz"
  ),
  local_path = file.path(
    "raw_data", "reference_panels", "ldsc", "t12_downloads",
    c("1000G_Phase3_baselineLD_v2.2_ldscores.tgz", "1000G_Phase3_plinkfiles.tgz", "1000G_Phase3_frq.tgz")
  ),
  expected_md5 = c(
    "b261e0caf06a003e7522938e01b3d349",
    "a7773ab485827b533cb300c76356d76b",
    "ac29686ffd5b6378789857a522ebca77"
  )
)
provenance[, absolute_path := file.path(project_root, local_path)]
provenance[, file_present := file.exists(absolute_path)]
provenance[, observed_md5 := vapply(absolute_path, function(path) {
  if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
}, character(1))]
provenance[, md5_pass := file_present & observed_md5 == expected_md5]
provenance[, absolute_path := NULL]
write_tsv(provenance, file.path(out_dir, "reference_download_provenance.tsv"))

validation <- data.table::data.table(
  check = c(
    "consensus_peak_file_present", "consensus_peak_count_298224", "source_peak_sets_28",
    "all_peaks_have_one_source_cell_type", "MHC_exclusion_rule_frozen",
    "Graves_sumstats_ready", "AIHT_sumstats_ready", "AITD_combined_sumstats_ready",
    "direct_case_case_sumstats_ready", "baselineLD_v2.2_ready", "official_per_chromosome_PLINK_ready",
    "frequency_files_ready", "GRCh38_to_GRCh37_liftover_ready"
  ),
  passed = c(
    file.exists(peaks_path), nrow(peaks) == 298224L, nrow(inventory) == 28L,
    all(peaks$n_source_cell_types == 1L), TRUE,
    phenotypes[requested_analysis == "Graves_case_control", usable_input],
    phenotypes[requested_analysis == "AIHT_case_control", usable_input],
    phenotypes[requested_analysis == "AITD_combined_case_control", usable_input],
    phenotypes[requested_analysis == "Graves_vs_AIHT_direct_case_case", usable_input],
    resource_gate[resource == "baselineLD_v2.2_22_chromosomes", ready],
    resource_gate[resource == "1000G_EUR_Phase3_PLINK_22_chromosomes", ready],
    resource_gate[resource == "1000G_Phase3_frequency_22_chromosomes", ready],
    resource_gate[resource == "GRCh38_to_GRCh37_chain", ready]
  )
)
write_tsv(validation, file.path(out_dir, "t12_preflight_validation.tsv"))

reference_ready <- all(resource_gate$ready)
phenotype_ready <- all(phenotypes$usable_input)
decision_code <- if (!reference_ready) {
  "HOLD_T12_MISSING_REFERENCE_AND_PHENOTYPE_INPUTS"
} else if (!phenotype_ready) {
  "HOLD_T12_MISSING_REQUIRED_PHENOTYPE_INPUTS_PARTIAL_TRAITS_ELIGIBLE"
} else {
  "GO_T12_ALL_REQUIRED_INPUTS_READY"
}
rationale <- if (!reference_ready) {
  paste0(
    "One or more build-matched S-LDSC reference resources are missing and two required phenotype inputs are absent. ",
    "CC-GWAS cannot be substituted for a direct case-case GWAS in S-LDSC."
  )
} else if (!phenotype_ready) {
  paste0(
    "All build-matched S-LDSC reference resources are ready and Graves/AIHT case-control analyses are eligible, ",
    "but AITD-union and direct Graves-vs-AIHT GWAS inputs remain absent. CC-GWAS cannot be substituted for the direct case-case analysis."
  )
} else {
  "All reference and phenotype inputs required by the frozen T12 scope are ready."
}
decision <- data.table::data.table(
  gate = "T12_S_LDSC_preflight",
  status = if (reference_ready && phenotype_ready) "GO" else "HOLD",
  decision = decision_code,
  eligible_traits_now = paste(phenotypes[usable_input == TRUE, requested_analysis], collapse = ";"),
  ineligible_required_traits = paste(phenotypes[usable_input == FALSE, requested_analysis], collapse = ";"),
  failed_checks = paste(validation[passed == FALSE, check], collapse = ";"),
  rationale = rationale,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_tsv(decision, file.path(out_dir, "t12_preflight_decision.tsv"))
log_message("T12 preflight complete: ", decision$status, " (", decision$decision, ").")
