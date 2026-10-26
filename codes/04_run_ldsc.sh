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
  if [[ -z "${value}" ]]; then
    return 0
  fi
  if [[ "${value}" = /* ]]; then
    printf '%s' "${value}"
  else
    printf '%s/%s' "${PROJECT_ROOT}" "${value}"
  fi
}

CONFIG_LDSC_DIR="$(config_value ldsc ldsc_dir)"
CONFIG_LDSC_PYTHON="$(config_value ldsc python)"
CONFIG_REF_PREFIX="$(config_value ldsc reference_ld_prefix)"
CONFIG_WEIGHT_PREFIX="$(config_value ldsc regression_weight_prefix)"
CONFIG_HAPMAP3="$(config_value ldsc hapmap3_snplist)"

LDSC_DIR="${LDSC_DIR:-$(resolve_project_path "${CONFIG_LDSC_DIR}")}"
LDSC_PYTHON="${LDSC_PYTHON:-$(resolve_project_path "${CONFIG_LDSC_PYTHON}")}"
LDSC_PYTHON="${LDSC_PYTHON:-python}"
LDSC_REF_PREFIX="${LDSC_REF_PREFIX:-$(resolve_project_path "${CONFIG_REF_PREFIX}")}"
LDSC_WEIGHT_PREFIX="${LDSC_WEIGHT_PREFIX:-$(resolve_project_path "${CONFIG_WEIGHT_PREFIX}")}"
HAPMAP3_SNPLIST="${HAPMAP3_SNPLIST:-$(resolve_project_path "${CONFIG_HAPMAP3}")}"
if [[ -z "${LDSC_DIR}" || -z "${LDSC_REF_PREFIX}" || -z "${LDSC_WEIGHT_PREFIX}" || -z "${HAPMAP3_SNPLIST}" ]]; then
  echo "Configure ldsc_dir, reference_ld_prefix, regression_weight_prefix, and hapmap3_snplist in config/analysis.yml or matching environment variables." >&2
  exit 2
fi
if [[ "${LDSC_PYTHON}" == */* ]]; then
  if [[ ! -x "${LDSC_PYTHON}" ]]; then
    echo "Configured LDSC Python is not executable: ${LDSC_PYTHON}" >&2
    exit 2
  fi
elif ! command -v "${LDSC_PYTHON}" >/dev/null 2>&1; then
  echo "Configured LDSC Python command was not found: ${LDSC_PYTHON}" >&2
  exit 2
fi
if [[ ! -f "${LDSC_DIR}/munge_sumstats.py" || ! -f "${LDSC_DIR}/ldsc.py" ]]; then
  echo "LDSC_DIR does not contain munge_sumstats.py and ldsc.py." >&2
  exit 2
fi
if [[ ! -f "${HAPMAP3_SNPLIST}" ]]; then
  echo "HAPMAP3_SNPLIST does not exist: ${HAPMAP3_SNPLIST}" >&2
  exit 2
fi
for chromosome in {1..22}; do
  for suffix in l2.ldscore.gz l2.M l2.M_5_50; do
    if [[ ! -f "${LDSC_REF_PREFIX}${chromosome}.${suffix}" ]]; then
      echo "Missing LDSC regression LD file: ${LDSC_REF_PREFIX}${chromosome}.${suffix}" >&2
      exit 2
    fi
  done
  if [[ ! -f "${LDSC_WEIGHT_PREFIX}${chromosome}.l2.ldscore.gz" ]]; then
    echo "Missing LDSC regression-weight file: ${LDSC_WEIGHT_PREFIX}${chromosome}.l2.ldscore.gz" >&2
    exit 2
  fi
done

INPUT_DIR="${PROJECT_ROOT}/outcome_data/ldsc/input"
MUNGED_DIR="${PROJECT_ROOT}/outcome_data/ldsc/munged"
RESULT_DIR="${PROJECT_ROOT}/outcome_data/ldsc/results"
mkdir -p "${MUNGED_DIR}" "${RESULT_DIR}"

while IFS=$'\t' read -r trait_id label role n_cases n_controls include_global include_locus; do
  if [[ "${trait_id}" == "trait_id" || "${include_global}" != "TRUE" ]]; then
    continue
  fi
  input="${INPUT_DIR}/${trait_id}.ldsc_input.tsv.gz"
  if [[ ! -f "${input}" ]]; then
    echo "Missing LDSC input: ${input}" >&2
    exit 2
  fi
  "${LDSC_PYTHON}" "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${input}" \
    --out "${MUNGED_DIR}/${trait_id}" \
    --merge-alleles "${HAPMAP3_SNPLIST}" \
    --snp SNP --a1 A1 --a2 A2 --p P --frq FRQ \
    --signed-sumstats BETA,0 \
    --N-cas "${n_cases}" --N-con "${n_controls}"
done < "${PROJECT_ROOT}/config/traits.tsv"

GRAVES="E4_GRAVES_STRICT"
AIHT="E4_HYTHY_AI_STRICT"
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

echo "LDSC complete. Review all three logs, then run codes/04_build_lava_overlap.R; do not assume independence for these FinnGen traits."
