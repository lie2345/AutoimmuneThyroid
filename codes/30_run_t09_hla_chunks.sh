#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$PROJECT_ROOT/outcome_data/targeted_susie_t08/hla/chunk_transport/chunk_manifest.tsv"
MAX_PARALLEL=32

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -f "$MANIFEST" ]] || { echo "Missing manifest: $MANIFEST" >&2; exit 2; }
[[ "$MAX_PARALLEL" =~ ^[0-9]+$ ]] && [[ "$MAX_PARALLEL" -ge 1 ]] || {
  echo "--max-parallel must be a positive integer" >&2
  exit 2
}

ACQUIRE="$PROJECT_ROOT/codes/10_acquire_targeted_susie.sh"
declare -A labels_by_pid=()
active=0

wait_one() {
  local completed_pid=""
  local label=""
  if wait -n -p completed_pid; then
    label="${labels_by_pid[$completed_pid]:-$completed_pid}"
    echo "COMPLETED $label"
  else
    label="${labels_by_pid[$completed_pid]:-$completed_pid}"
    echo "FAILED $label" >&2
    unset 'labels_by_pid[$completed_pid]'
    active=$((active - 1))
    return 1
  fi
  unset 'labels_by_pid[$completed_pid]'
  active=$((active - 1))
}

failures=0
while IFS=$'\t' read -r source_id chunk_id n_chunks source_path original_region chunk_region plan_path slice_root chunk_data_path log_path; do
  [[ "$source_id" == "source_id" ]] && continue
  while [[ "$active" -ge "$MAX_PARALLEL" ]]; do
    wait_one || failures=$((failures + 1))
  done
  console_log="${log_path%.tsv}.console.log"
  (
    bash "$ACQUIRE" \
      --mode extract \
      --max-priority P0_critical \
      --plan "$plan_path" \
      --slice-root "$slice_root" \
      --log "$log_path" \
      --run \
      --keep-going
  ) >"$console_log" 2>&1 &
  pid="$!"
  labels_by_pid[$pid]="source${source_id}_chunk${chunk_id}"
  active=$((active + 1))
done < "$MANIFEST"

while [[ "$active" -gt 0 ]]; do
  wait_one || failures=$((failures + 1))
done

if [[ "$failures" -gt 0 ]]; then
  echo "$failures HLA chunk job(s) failed; inspect per-job logs." >&2
  exit 1
fi
echo "All HLA chunk transport jobs completed."
