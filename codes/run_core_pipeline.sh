#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:---dry-run}"

if [[ "${MODE}" != "--dry-run" && "${MODE}" != "--run" ]]; then
  echo "Usage: bash codes/run_core_pipeline.sh [--dry-run|--run]" >&2
  exit 2
fi

run_step() {
  if [[ "${MODE}" == "--dry-run" ]]; then
    printf 'DRY RUN:'
    printf ' %q' "$@"
    printf '\n'
  else
    "$@"
  fi
}

run_step Rscript "${PROJECT_ROOT}/codes/00_validate_inputs.R"
run_step Rscript "${PROJECT_ROOT}/codes/01_harmonize_gwas.R"
run_step Rscript "${PROJECT_ROOT}/codes/02_define_loci.R"
run_step Rscript "${PROJECT_ROOT}/codes/03_prepare_ld_methods.R"
run_step Rscript "${PROJECT_ROOT}/codes/06_filter_official_r12_coloc.R"
run_step Rscript "${PROJECT_ROOT}/codes/07_coloc_multiome_r13.R"

echo "LDSC and LAVA are intentionally not included in the automatic core run because their external LD resources and sample-overlap inputs must be configured and reviewed first."
echo "After those checks, run codes/04_run_ldsc.sh, codes/04_build_lava_overlap.R, codes/05_run_lava.R, then codes/08_integrate_regulatory_evidence.R."
