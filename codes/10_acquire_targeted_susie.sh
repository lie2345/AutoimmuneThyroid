#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_HTSLIB_BIN="$PROJECT_ROOT/tools/htslib-1.19-curl/bin"
TABIX_BIN="${TABIX_BIN:-$PROJECT_HTSLIB_BIN/tabix}"
BGZIP_BIN="${BGZIP_BIN:-$PROJECT_HTSLIB_BIN/bgzip}"
if [[ ! -x "$TABIX_BIN" ]]; then TABIX_BIN="$(command -v tabix || true)"; fi
if [[ ! -x "$BGZIP_BIN" ]]; then BGZIP_BIN="$(command -v bgzip || true)"; fi

MODE="extract"
MAX_PRIORITY="P1_strong"
PLAN="$PROJECT_ROOT/outcome_data/targeted_susie/acquisition_plan.tsv"
CACHE_MB="64"
DO_RUN=0
FORCE=0
KEEP_GOING=0
SLICE_ROOT="$PROJECT_ROOT/outcome_data/targeted_susie/slices"
DOWNLOAD_ROOT="$PROJECT_ROOT/raw_data/finngen_multiome"
LOG_PATH="$PROJECT_ROOT/outcome_data/targeted_susie/acquisition_log.tsv"

usage() {
  cat <<'EOF'
Usage:
  bash codes/10_acquire_targeted_susie.sh [options]

Options:
  --mode extract|download       Remote Tabix locus slices (recommended) or full chromosome files.
  --max-priority LEVEL          P0_critical, P1_strong, or P2_multisignal (default: P1_strong).
  --plan FILE                   Stage-specific plan from codes/09_plan_targeted_susie.R.
  --cache-mb INT                HTSlib remote BGZF cache in MB (default: 64).
  --slice-root DIR              Destination for targeted slices.
  --download-root DIR           Destination root for full files.
  --log FILE                    Acquisition audit log.
  --run                         Execute network/file operations. Without this flag: dry-run only.
  --force                       Replace an existing destination after a new file validates.
  --keep-going                  Continue after a failed file and exit non-zero at the end.
  --help                        Show this help.

Examples:
  bash codes/10_acquire_targeted_susie.sh --mode extract --max-priority P0_critical
  bash codes/10_acquire_targeted_susie.sh --mode extract --max-priority P0_critical --run
  bash codes/10_acquire_targeted_susie.sh --mode download --max-priority P0_critical --run
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --max-priority) MAX_PRIORITY="$2"; shift 2 ;;
    --plan) PLAN="$2"; shift 2 ;;
    --cache-mb) CACHE_MB="$2"; shift 2 ;;
    --slice-root) SLICE_ROOT="$2"; shift 2 ;;
    --download-root) DOWNLOAD_ROOT="$2"; shift 2 ;;
    --log) LOG_PATH="$2"; shift 2 ;;
    --run) DO_RUN=1; shift ;;
    --force) FORCE=1; shift ;;
    --keep-going) KEEP_GOING=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$MODE" != "extract" && "$MODE" != "download" ]]; then
  echo "--mode must be extract or download." >&2
  exit 2
fi
case "$MAX_PRIORITY" in
  P0_critical|P1_strong|P2_multisignal) ;;
  *) echo "Invalid --max-priority: $MAX_PRIORITY" >&2; exit 2 ;;
esac
if [[ ! "$CACHE_MB" =~ ^[0-9]+$ ]]; then
  echo "--cache-mb must be a non-negative integer." >&2
  exit 2
fi
if [[ ! -s "$PLAN" ]]; then
  echo "Missing or empty acquisition plan: $PLAN" >&2
  echo "Run: Rscript codes/09_plan_targeted_susie.R" >&2
  exit 2
fi

EXPECTED_HEADER=$'stage\tmax_priority_rank\tmodality\tlevel\tcell_type\tchrom\tsource_path\tdata_url\tindex_url\tsize_bytes\tmd5\tindex_size_bytes\tindex_md5\tn_loci\ttabix_regions\tselection_reasons\trecommended_mode'
ACTUAL_HEADER="$(head -n 1 "$PLAN")"
if [[ "$ACTUAL_HEADER" != "$EXPECTED_HEADER" ]]; then
  echo "Unexpected acquisition plan header. Regenerate it with codes/09_plan_targeted_susie.R." >&2
  exit 2
fi

STAGE="through_${MAX_PRIORITY}"
SELECTED_PLAN="$(mktemp)"
cleanup() {
  if [[ -f "$SELECTED_PLAN" ]]; then
    rm -f -- "$SELECTED_PLAN"
  fi
}
trap cleanup EXIT

awk -F '\t' -v stage="$STAGE" 'NR == 1 || $1 == stage' "$PLAN" > "$SELECTED_PLAN"
N_FILES=$(( $(wc -l < "$SELECTED_PLAN") - 1 ))
if [[ "$N_FILES" -le 0 ]]; then
  echo "No files are planned for stage: $STAGE" >&2
  exit 2
fi

if [[ "$DO_RUN" -eq 1 ]]; then
  if [[ "$MODE" == "extract" ]]; then
    [[ -n "$TABIX_BIN" && -x "$TABIX_BIN" ]] || {
      echo "tabix is required for --mode extract. Install an HTSlib tabix CLI first." >&2
      exit 2
    }
    [[ -n "$BGZIP_BIN" && -x "$BGZIP_BIN" ]] || {
      echo "bgzip is required for --mode extract. Install an HTSlib bgzip CLI first." >&2
      exit 2
    }
    command -v curl >/dev/null 2>&1 || {
      echo "curl is required for remote header range requests." >&2
      exit 2
    }
  else
    command -v curl >/dev/null 2>&1 || { echo "curl is required." >&2; exit 2; }
    command -v md5sum >/dev/null 2>&1 || { echo "md5sum is required." >&2; exit 2; }
  fi
fi

mkdir -p "$(dirname "$LOG_PATH")"
if [[ ! -f "$LOG_PATH" ]]; then
  printf 'timestamp\tmode\tstage\tsource_path\tdestination\tstatus\tdetail\n' > "$LOG_PATH"
fi

log_status() {
  local source_path="$1"
  local destination="$2"
  local status="$3"
  local detail="$4"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -Iseconds)" "$MODE" "$STAGE" "$source_path" "$destination" "$status" "$detail" >> "$LOG_PATH"
}

fail_or_continue() {
  local message="$1"
  echo "$message" >&2
  if [[ "$KEEP_GOING" -eq 0 ]]; then
    return 1
  fi
  return 0
}

validate_file() {
  local path="$1"
  local expected_size="$2"
  local expected_md5="$3"
  [[ -s "$path" ]] || return 1
  [[ "$(stat -c '%s' "$path")" == "$expected_size" ]] || return 1
  [[ "$(md5sum "$path" | awk '{print $1}')" == "$expected_md5" ]] || return 1
}

download_one() {
  local url="$1"
  local destination="$2"
  local expected_size="$3"
  local expected_md5="$4"
  local part="${destination}.part"

  mkdir -p "$(dirname "$destination")"
  if [[ -f "$destination" ]]; then
    if validate_file "$destination" "$expected_size" "$expected_md5"; then
      return 10
    fi
    if [[ "$FORCE" -eq 0 ]]; then
      return 11
    fi
    part="${destination}.replacement.$$"
  fi

  curl --fail --location --retry 5 --retry-all-errors --continue-at - --output "$part" "$url"
  if ! validate_file "$part" "$expected_size" "$expected_md5"; then
    return 12
  fi
  mv -f -- "$part" "$destination"
  return 0
}

extract_one() {
  local data_url="$1"
  local destination="$2"
  local tabix_regions="$3"
  local temporary="${destination}.part.$$"
  local header_chunk="${temporary}.header.bgzf"
  local header_line=""
  local -a regions
  local remote_rc=1
  local attempt
  local extracted_records

  IFS=';' read -r -a regions <<< "$tabix_regions"
  if [[ "${#regions[@]}" -eq 0 || "$tabix_regions" == "NA" ]]; then
    return 20
  fi
  mkdir -p "$(dirname "$destination")"
  if [[ -e "$destination" || -e "${destination}.tbi" ]]; then
    if [[ -s "$destination" && -s "${destination}.tbi" && "$FORCE" -eq 0 ]]; then
      return 10
    fi
    if [[ "$FORCE" -eq 0 ]]; then
      return 11
    fi
  fi

  if [[ "$data_url" =~ ^https?:// ]]; then
    curl --silent --show-error --fail --location --retry 3 --retry-all-errors \
      --connect-timeout 20 --max-time 120 --range 0-65535 \
      --output "$header_chunk" "$data_url"
    if [[ -s "$header_chunk" ]]; then
      header_line="$("$BGZIP_BIN" -dc "$header_chunk" 2>/dev/null | head -n 1)"
    fi
    rm -f -- "$header_chunk"
  else
    header_line="$("$BGZIP_BIN" -dc "$data_url" 2>/dev/null | head -n 1)"
  fi
  if [[ -z "$header_line" || "$(awk -F '\t' '{print $5}' <<< "$header_line")" != "chromosome" || \
        "$(awk -F '\t' '{print $6}' <<< "$header_line")" != "position" ]]; then
    return 23
  fi

  for attempt in 1 2 3; do
    rm -f -- "$temporary" "${temporary}.tbi"
    {
      printf '%s\n' "$header_line"
      "$TABIX_BIN" -D --cache "$CACHE_MB" "$data_url" "${regions[@]}"
    } | "$BGZIP_BIN" -c > "$temporary"
    remote_rc=$?
    if [[ "$remote_rc" -eq 0 && -s "$temporary" ]]; then
      break
    fi
    if [[ "$attempt" -lt 3 ]]; then
      echo "Remote extraction attempt $attempt failed; retrying: $data_url" >&2
      sleep "$attempt"
    fi
  done

  if [[ "$remote_rc" -ne 0 || ! -s "$temporary" ]]; then
    return 21
  fi
  "$TABIX_BIN" -f -s 5 -b 6 -e 6 -S 1 "$temporary"
  extracted_records="$("$TABIX_BIN" "$temporary" "${regions[@]}" | awk 'NF {n++} END {print n+0}')"
  if [[ "$extracted_records" -lt 1 ]]; then
    rm -f -- "$temporary" "${temporary}.tbi"
    return 22
  fi
  mv -f -- "$temporary" "$destination"
  mv -f -- "${temporary}.tbi" "${destination}.tbi"
  return 0
}

echo "Mode: $MODE"
echo "Stage: $STAGE"
echo "Files: $N_FILES"
if [[ "$DO_RUN" -eq 0 ]]; then
  echo "DRY-RUN: no network access or data-file writes will occur. Add --run to execute."
fi

FAILURES=0
while IFS=$'\t' read -r stage max_rank modality level cell_type chrom source_path data_url index_url size_bytes md5 index_size_bytes index_md5 n_loci tabix_regions selection_reasons recommended_mode; do
  [[ "$stage" == "stage" ]] && continue

  if [[ "$MODE" == "extract" ]]; then
    relative_destination="${source_path%.SUSIE.snp.tsv.gz}.SUSIE.snp.targeted.tsv.gz"
    destination="$SLICE_ROOT/$STAGE/$relative_destination"
    if [[ "$DO_RUN" -eq 0 ]]; then
      echo "EXTRACT $source_path -> $destination ($n_loci loci; $tabix_regions)"
      continue
    fi
    set +e
    extract_one "$data_url" "$destination" "$tabix_regions"
    rc=$?
    set -e
    case "$rc" in
      0) log_status "$source_path" "$destination" "completed" "$n_loci loci" ;;
      10) log_status "$source_path" "$destination" "skipped_existing" "data and index already exist" ;;
      *)
        FAILURES=$((FAILURES + 1))
        log_status "$source_path" "$destination" "failed" "extract_rc=$rc"
        fail_or_continue "Failed to extract $source_path (rc=$rc)" || exit 1
        ;;
    esac
  else
    destination="$DOWNLOAD_ROOT/$source_path"
    index_destination="${destination}.tbi"
    if [[ "$DO_RUN" -eq 0 ]]; then
      gib="$(awk -v bytes="$size_bytes" 'BEGIN {printf "%.3f", bytes / 1073741824}')"
      echo "DOWNLOAD $source_path -> $destination (${gib} GiB plus index)"
      continue
    fi

    set +e
    download_one "$data_url" "$destination" "$size_bytes" "$md5"
    data_rc=$?
    set -e
    if [[ "$data_rc" -ne 0 && "$data_rc" -ne 10 ]]; then
      FAILURES=$((FAILURES + 1))
      log_status "$source_path" "$destination" "failed" "data_download_rc=$data_rc"
      fail_or_continue "Failed validation or replacement guard for $source_path (rc=$data_rc)" || exit 1
      continue
    fi

    set +e
    download_one "$index_url" "$index_destination" "$index_size_bytes" "$index_md5"
    index_rc=$?
    set -e
    if [[ "$index_rc" -ne 0 && "$index_rc" -ne 10 ]]; then
      FAILURES=$((FAILURES + 1))
      log_status "$source_path" "$destination" "failed" "index_download_rc=$index_rc"
      fail_or_continue "Failed index validation for $source_path (rc=$index_rc)" || exit 1
      continue
    fi
    if [[ "$data_rc" -eq 10 && "$index_rc" -eq 10 ]]; then
      log_status "$source_path" "$destination" "skipped_existing" "data and index already validated"
    else
      log_status "$source_path" "$destination" "completed" "MD5 and byte size validated"
    fi
  fi
done < "$SELECTED_PLAN"

if [[ "$FAILURES" -gt 0 ]]; then
  echo "Completed with $FAILURES failed file(s). See: $LOG_PATH" >&2
  exit 1
fi
if [[ "$DO_RUN" -eq 1 ]]; then
  echo "Acquisition completed. Audit log: $LOG_PATH"
else
  echo "Dry-run completed. Review the commands and estimated full-download size before adding --run."
fi
