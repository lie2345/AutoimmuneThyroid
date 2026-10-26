#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_HTSLIB_BIN="$PROJECT_ROOT/tools/htslib-1.19-curl/bin"
TABIX_BIN="${TABIX_BIN:-$PROJECT_HTSLIB_BIN/tabix}"
BGZIP_BIN="${BGZIP_BIN:-$PROJECT_HTSLIB_BIN/bgzip}"
MANIFEST="$PROJECT_ROOT/outcome_data/targeted_susie_t08/hla/chunk_transport/chunk_manifest.tsv"
FINAL_ROOT="$PROJECT_ROOT/outcome_data/targeted_susie_t08/hla/slices"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --final-root) FINAL_ROOT="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -x "$TABIX_BIN" && -x "$BGZIP_BIN" ]] || {
  echo "Project-local tabix and bgzip are required." >&2
  exit 2
}
[[ -f "$MANIFEST" ]] || { echo "Missing manifest: $MANIFEST" >&2; exit 2; }

mapfile -t source_ids < <(awk -F '\t' 'NR>1 {print $1}' "$MANIFEST" | sort -n -u)
for source_id in "${source_ids[@]}"; do
  source_path="$(awk -F '\t' -v id="$source_id" 'NR>1 && $1==id {print $4; exit}' "$MANIFEST")"
  original_region="$(awk -F '\t' -v id="$source_id" 'NR>1 && $1==id {print $5; exit}' "$MANIFEST")"
  stage="through_P0_critical"
  relative_destination="${source_path%.SUSIE.snp.tsv.gz}.SUSIE.snp.targeted.tsv.gz"
  destination="$FINAL_ROOT/$stage/$relative_destination"
  temporary="${destination}.merge.part.$$"
  mkdir -p "$(dirname "$destination")"

  mapfile -t chunks < <(
    awk -F '\t' -v id="$source_id" 'NR>1 && $1==id {print $2 "\t" $9}' "$MANIFEST" |
      sort -n -k1,1 | cut -f2-
  )
  expected_chunks="$(awk -F '\t' -v id="$source_id" 'NR>1 && $1==id {print $3; exit}' "$MANIFEST")"
  [[ "${#chunks[@]}" -eq "$expected_chunks" ]] || {
    echo "Chunk count mismatch for $source_path" >&2
    exit 1
  }

  for chunk in "${chunks[@]}"; do
    [[ -s "$chunk" && -s "${chunk}.tbi" ]] || {
      echo "Missing chunk data/index: $chunk" >&2
      exit 1
    }
    "$BGZIP_BIN" -t "$chunk"
  done

  {
    first=1
    for chunk in "${chunks[@]}"; do
      if [[ "$first" -eq 1 ]]; then
        "$BGZIP_BIN" -cd "$chunk"
        first=0
      else
        "$BGZIP_BIN" -cd "$chunk" | awk 'NR > 1'
      fi
    done
  } | "$BGZIP_BIN" -c > "$temporary"

  "$BGZIP_BIN" -t "$temporary"
  "$TABIX_BIN" -f -s 5 -b 6 -e 6 -S 1 "$temporary"
  records="$($TABIX_BIN "$temporary" "$original_region" | awk 'NF {n++} END {print n+0}')"
  [[ "$records" -gt 0 ]] || { echo "Merged file has no planned-region records" >&2; exit 1; }
  mv -f -- "$temporary" "$destination"
  mv -f -- "${temporary}.tbi" "${destination}.tbi"
  echo "MERGED source${source_id}: $source_path ($records records)"
done

echo "All HLA chunks merged and re-indexed under: $FINAL_ROOT"
