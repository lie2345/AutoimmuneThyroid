#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mode="${1:---dry-run}"
if [[ "$mode" != "--dry-run" && "$mode" != "--run" ]]; then
  echo "Usage: bash codes/43_run_t12_custom_ldscores.sh [--dry-run|--run]" >&2
  exit 2
fi

python_bin="${LDSC_PYTHON:-$project_root/software/ldsc39-venv/bin/python}"
ldsc_py="$project_root/software/ldsc39/ldsc.py"
plink_dir="$project_root/raw_data/reference_panels/ldsc/1000G_EUR_Phase3_plink"
ann_dir="$project_root/outcome_data/t12/annotations"
out_dir="$project_root/outcome_data/t12/ldscores/master"
hm3_source="$project_root/raw_data/reference_panels/ldsc/w_hm3.snplist"
hm3_ids="$project_root/outcome_data/t12/reference/w_hm3.ids"
mkdir -p "$out_dir" "$(dirname "$hm3_ids")"

if [[ ! -s "$hm3_ids" ]]; then
  awk 'NR > 1 {print $1}' "$hm3_source" > "$hm3_ids"
fi
[[ "$(wc -l < "$hm3_ids")" -eq 1217311 ]] || { echo "Unexpected HM3 ID count." >&2; exit 1; }

is_valid() {
  local prefix="$1"
  [[ -s "${prefix}.l2.ldscore.gz" && -s "${prefix}.l2.M" && -s "${prefix}.l2.M_5_50" ]] || return 1
  gzip -t "${prefix}.l2.ldscore.gz" || return 1
  [[ "$(zcat "${prefix}.l2.ldscore.gz" 2>/dev/null | awk 'NR==1 {print NF; exit}')" -eq 38 ]] || return 1
  [[ "$(awk 'NR==1 {print NF; exit}' "${prefix}.l2.M")" -eq 35 ]] || return 1
  [[ "$(awk 'NR==1 {print NF; exit}' "${prefix}.l2.M_5_50")" -eq 35 ]] || return 1
}

run_one() {
  local chr="$1"
  local prefix="$out_dir/t12_source_labels.$chr"
  if is_valid "$prefix"; then
    echo "chr${chr}: already valid; skip"
    return 0
  fi
  local -a cmd=(
    "$python_bin" "$ldsc_py"
    --l2
    --bfile "$plink_dir/1000G.EUR.QC.$chr"
    --ld-wind-cm 1
    --annot "$ann_dir/t12_source_labels.$chr.annot.gz"
    --thin-annot
    --print-snps "$hm3_ids"
    --out "$prefix"
  )
  if [[ "$mode" == "--dry-run" ]]; then
    printf 'chr%s:' "$chr"
    printf ' %q' "${cmd[@]}"
    printf '\n'
  else
    echo "chr${chr}: running"
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      "${cmd[@]}" > "${prefix}.stdout.log" 2>&1
    is_valid "$prefix" || { echo "chr${chr}: output validation failed" >&2; return 1; }
    echo "chr${chr}: PASS"
  fi
}

if [[ "$mode" == "--dry-run" ]]; then
  for chr in $(seq 1 22); do run_one "$chr"; done
  exit 0
fi

max_jobs=4
pids=()
labels=()
flush_jobs() {
  local failed=0
  local i
  for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
      echo "chr${labels[$i]}: FAILED; see stdout log" >&2
      failed=1
    fi
  done
  pids=()
  labels=()
  [[ "$failed" -eq 0 ]]
}

for chr in $(seq 1 22); do
  run_one "$chr" &
  pids+=("$!")
  labels+=("$chr")
  if [[ "${#pids[@]}" -ge "$max_jobs" ]]; then flush_jobs; fi
done
if [[ "${#pids[@]}" -gt 0 ]]; then flush_jobs; fi

echo "T12 custom LD-score generation: 22/22 chromosomes PASS"
