#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONFIG_PATH="${PROJECT_ROOT}/config/analysis.yml"

config_value() {
  local section="$1"
  local key="$2"
  Rscript --vanilla -e 'args <- commandArgs(trailingOnly = TRUE); cfg <- yaml::read_yaml(args[[1]]); value <- cfg[[args[[2]]]][[args[[3]]]]; if (!is.null(value)) cat(value)' \
    "${CONFIG_PATH}" "${section}" "${key}"
}

resolve_project_path() {
  local value="$1"
  if [[ "${value}" = /* ]]; then
    printf '%s' "${value}"
  else
    printf '%s/%s' "${PROJECT_ROOT}" "${value}"
  fi
}

LDSC_DIR="$(resolve_project_path "$(config_value ldsc ldsc_dir)")"
LDSC_PYTHON="$(resolve_project_path "$(config_value ldsc python)")"
LDSC_REF_PREFIX="$(resolve_project_path "$(config_value ldsc reference_ld_prefix)")"
LDSC_WEIGHT_PREFIX="$(resolve_project_path "$(config_value ldsc regression_weight_prefix)")"
HAPMAP3_SNPLIST="$(resolve_project_path "$(config_value ldsc hapmap3_snplist)")"

GRAVES="$(config_value primary_comparison graves)"
AIHT="$(config_value primary_comparison autoimmune_hypothyroidism)"
INPUT_DIR="${PROJECT_ROOT}/outcome_data/ldsc/input"
OUT_ROOT="${PROJECT_ROOT}/outcome_data/t05/ldsc"
MUNGED_DIR="${OUT_ROOT}/munged"
RESULT_DIR="${OUT_ROOT}/results"
mkdir -p "${MUNGED_DIR}" "${RESULT_DIR}"

if [[ ! -x "${LDSC_PYTHON}" ]]; then
  echo "Configured LDSC Python is not executable: ${LDSC_PYTHON}" >&2
  exit 2
fi
if [[ ! -f "${LDSC_DIR}/munge_sumstats.py" || ! -f "${LDSC_DIR}/ldsc.py" ]]; then
  echo "Missing LDSC scripts under ${LDSC_DIR}" >&2
  exit 2
fi
if [[ ! -f "${HAPMAP3_SNPLIST}" ]]; then
  echo "Missing HapMap3 allele list: ${HAPMAP3_SNPLIST}" >&2
  exit 2
fi

for chromosome in {1..22}; do
  for suffix in l2.ldscore.gz l2.M l2.M_5_50; do
    if [[ ! -f "${LDSC_REF_PREFIX}${chromosome}.${suffix}" ]]; then
      echo "Missing LDSC reference file: ${LDSC_REF_PREFIX}${chromosome}.${suffix}" >&2
      exit 2
    fi
  done
  if [[ ! -f "${LDSC_WEIGHT_PREFIX}${chromosome}.l2.ldscore.gz" ]]; then
    echo "Missing LDSC weight file: ${LDSC_WEIGHT_PREFIX}${chromosome}.l2.ldscore.gz" >&2
    exit 2
  fi
done

trait_counts() {
  local requested_trait="$1"
  awk -F '\t' -v id="${requested_trait}" '$1==id {print $4"\t"$5}' "${PROJECT_ROOT}/config/traits.tsv"
}

for trait_id in "${GRAVES}" "${AIHT}"; do
  input="${INPUT_DIR}/${trait_id}.ldsc_input.tsv.gz"
  if [[ ! -f "${input}" ]]; then
    echo "Missing non-HLA LDSC input: ${input}" >&2
    exit 2
  fi
  IFS=$'\t' read -r n_cases n_controls <<< "$(trait_counts "${trait_id}")"
  if [[ -z "${n_cases}" || -z "${n_controls}" ]]; then
    echo "Missing sample counts for ${trait_id}" >&2
    exit 2
  fi
  "${LDSC_PYTHON}" "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${input}" \
    --out "${MUNGED_DIR}/${trait_id}" \
    --merge-alleles "${HAPMAP3_SNPLIST}" \
    --snp SNP --a1 A1 --a2 A2 --p P --frq FRQ \
    --signed-sumstats BETA,0 \
    --N-cas "${n_cases}" --N-con "${n_controls}"
done

"${LDSC_PYTHON}" "${LDSC_DIR}/ldsc.py" \
  --rg "${MUNGED_DIR}/${GRAVES}.sumstats.gz,${MUNGED_DIR}/${GRAVES}.sumstats.gz" \
  --ref-ld-chr "${LDSC_REF_PREFIX}" \
  --w-ld-chr "${LDSC_WEIGHT_PREFIX}" \
  --out "${RESULT_DIR}/${GRAVES}__${GRAVES}"

"${LDSC_PYTHON}" "${LDSC_DIR}/ldsc.py" \
  --rg "${MUNGED_DIR}/${AIHT}.sumstats.gz,${MUNGED_DIR}/${AIHT}.sumstats.gz" \
  --ref-ld-chr "${LDSC_REF_PREFIX}" \
  --w-ld-chr "${LDSC_WEIGHT_PREFIX}" \
  --out "${RESULT_DIR}/${AIHT}__${AIHT}"

"${LDSC_PYTHON}" "${LDSC_DIR}/ldsc.py" \
  --rg "${MUNGED_DIR}/${GRAVES}.sumstats.gz,${MUNGED_DIR}/${AIHT}.sumstats.gz" \
  --ref-ld-chr "${LDSC_REF_PREFIX}" \
  --w-ld-chr "${LDSC_WEIGHT_PREFIX}" \
  --out "${RESULT_DIR}/${GRAVES}__${AIHT}"

Rscript "${PROJECT_ROOT}/codes/17_parse_t05_ldsc.R"
