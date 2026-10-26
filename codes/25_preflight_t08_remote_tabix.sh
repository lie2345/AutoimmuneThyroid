#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_HTSLIB_BIN="$PROJECT_ROOT/tools/htslib-1.19-curl/bin"
TABIX_BIN="${TABIX_BIN:-$PROJECT_HTSLIB_BIN/tabix}"
if [[ ! -x "$TABIX_BIN" ]]; then TABIX_BIN="$(command -v tabix || true)"; fi
TARGETS="${1:-$PROJECT_ROOT/outcome_data/t08/remote_preflight_targets.tsv}"
OUTPUT="${2:-$PROJECT_ROOT/outcome_data/t08/remote_tabix_preflight.tsv}"

if [[ ! -s "$TARGETS" ]]; then
  echo "Missing or empty preflight target table: $TARGETS" >&2
  exit 2
fi
[[ -n "$TABIX_BIN" && -x "$TABIX_BIN" ]] || { echo "tabix is required." >&2; exit 2; }
command -v timeout >/dev/null 2>&1 || { echo "timeout is required." >&2; exit 2; }

EXPECTED_HEADER=$'source_path\tmodality\tlevel\tcell_type\tlocus_id\tis_hla\tdata_url\tindex_url\tchrom\ttop_h4_variant\tPP.H4.abf\ttest_position\ttest_region'
ACTUAL_HEADER="$(head -n 1 "$TARGETS")"
if [[ "$ACTUAL_HEADER" != "$EXPECTED_HEADER" ]]; then
  echo "Unexpected preflight target header: $TARGETS" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUTPUT")"
temporary="${OUTPUT}.part.$$"
error_file="$(mktemp)"
cleanup() {
  if [[ -f "$temporary" ]]; then
    rm -f -- "$temporary"
  fi
  if [[ -f "$error_file" ]]; then
    rm -f -- "$error_file"
  fi
}
trap cleanup EXIT

printf 'source_path\tlocus_id\tis_hla\tchrom\ttest_region\tindex_status\tindex_attempts\tchromosome_status\tquery_status\tquery_attempts\tn_records\tstatus\tdetail\n' > "$temporary"

failures=0
while IFS=$'\t' read -r source_path modality level cell_type locus_id is_hla data_url index_url chrom top_h4_variant pp_h4 test_position test_region; do
  [[ "$source_path" == "source_path" ]] && continue

  index_status="not_tested"
  chromosome_status="not_tested"
  query_status="not_tested"
  n_records=0
  status="ok"
  detail="remote index, chromosome naming, and known-variant point query succeeded"
  index_attempts=0
  query_attempts=0

  query_rc=1
  records=""
  for query_attempts in 1 2 3; do
    set +e
    records="$(timeout 45 "$TABIX_BIN" -D --cache 16 "$data_url" "$test_region" 2>"$error_file")"
    query_rc=$?
    set -e
    [[ "$query_rc" -eq 0 ]] && break
    sleep "$query_attempts"
  done
  if [[ "$query_rc" -eq 0 ]]; then
    n_records="$(awk 'NF {n++} END {print n+0}' <<< "$records")"
    index_status="validated_by_query"
    chromosome_status="validated_by_query"
    if [[ "$n_records" -lt 1 ]]; then
      query_status="no_record"
      status="failed"
      detail="known Step07 top variant was not returned at its indexed position"
    else
      query_status="ok"
    fi
  else
    query_status="failed"
    status="failed"
    query_error="$(tr '\t\r\n' ' ' < "$error_file" | sed 's/  */ /g' | cut -c1-180)"

    index_rc=1
    chromosomes=""
    for index_attempts in 1 2 3; do
      set +e
      chromosomes="$(timeout 45 "$TABIX_BIN" -D -l "$data_url" 2>"$error_file")"
      index_rc=$?
      set -e
      [[ "$index_rc" -eq 0 ]] && break
      sleep "$index_attempts"
    done
    if [[ "$index_rc" -ne 0 ]]; then
      index_status="failed"
      last_error="$(tr '\t\r\n' ' ' < "$error_file" | sed 's/  */ /g' | cut -c1-240)"
      detail="point query and index read failed after retries: query=${query_error}; index=${last_error}"
    elif ! awk -v target="$chrom" '$0 == target {found=1} END {exit !found}' <<< "$chromosomes"; then
      index_status="ok"
      chromosome_status="missing"
      detail="point query failed and expected chromosome is absent from remote Tabix index: ${query_error}"
    else
      index_status="ok"
      chromosome_status="ok"
      detail="remote index and chromosome are readable, but the known-variant query failed: ${query_error}"
    fi
  fi

  if [[ "$status" != "ok" ]]; then
    failures=$((failures + 1))
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$source_path" "$locus_id" "$is_hla" "$chrom" "$test_region" "$index_status" \
    "$index_attempts" "$chromosome_status" "$query_status" "$query_attempts" "$n_records" \
    "$status" "$detail" >> "$temporary"
done < "$TARGETS"

mv -f -- "$temporary" "$OUTPUT"
total=$(( $(wc -l < "$OUTPUT") - 1 ))
echo "Remote Tabix preflight: $((total - failures))/$total passed"
echo "Audit: $OUTPUT"

if [[ "$failures" -ne 0 ]]; then
  exit 1
fi
