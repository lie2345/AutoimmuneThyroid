#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mode="${1:---dry-run}"
if [[ "$mode" != "--dry-run" && "$mode" != "--run" ]]; then
  echo "Usage: bash codes/45_run_t12_h2cts.sh [--dry-run|--run]" >&2
  exit 2
fi

python_bin="${LDSC_PYTHON:-$project_root/software/ldsc39-venv/bin/python}"
ldsc_py="$project_root/software/ldsc39/ldsc.py"
baseline="$project_root/raw_data/reference_panels/ldsc/baselineLD_v2.2/baselineLD."
weights="$project_root/raw_data/reference_panels/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
model_dir="$project_root/outcome_data/t12/models"
result_dir="$project_root/outcome_data/t12/results"
mkdir -p "$result_dir"

declare -A traits=(
  [Graves]="$project_root/outcome_data/t05/ldsc/munged/E4_GRAVES_STRICT.sumstats.gz"
  [AIHT]="$project_root/outcome_data/t05/ldsc/munged/E4_HYTHY_AI_STRICT.sumstats.gz"
)
models=(global_cell global_lineage within_lineage)

for trait in Graves AIHT; do
  for model in "${models[@]}"; do
    out="$result_dir/${trait}.${model}"
    result="${out}.cell_type_results.txt"
    if [[ -s "$result" ]]; then
      echo "$trait/$model: already present; skip"
      continue
    fi
    cmd=(
      "$python_bin" "$ldsc_py"
      --h2-cts "${traits[$trait]}"
      --ref-ld-chr "$baseline"
      --ref-ld-chr-cts "$model_dir/$model.ldcts"
      --w-ld-chr "$weights"
      --out "$out"
    )
    if [[ "$mode" == "--dry-run" ]]; then
      printf '%s/%s:' "$trait" "$model"
      printf ' %q' "${cmd[@]}"
      printf '\n'
    else
      echo "$trait/$model: running"
      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        "${cmd[@]}" > "${out}.stdout.log" 2>&1
      [[ -s "$result" ]] || { echo "$trait/$model: missing result" >&2; exit 1; }
      echo "$trait/$model: PASS"
    fi
  done
done

if [[ "$mode" == "--run" ]]; then
  echo "T12 h2-cts eligible-trait analyses: 6/6 PASS"
fi
