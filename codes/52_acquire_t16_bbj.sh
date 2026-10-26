#!/usr/bin/env bash
# Fetch only the public BBJ Hashimoto source.  Promotion to .gz is owned by
# the offline QC script, which checks the published size and MD5 first.
set -euo pipefail

cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
raw_dir='raw_data/external_gwas/bbj_2021_hashimoto'
file='GCST90018635_buildGRCh37.tsv.gz'
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018635/GCST90018635_buildGRCh37.tsv.gz'
mkdir -p -- "$raw_dir"

# A final filename is never overwritten. curl resumes only the explicit .part.
if [[ ! -f "$raw_dir/$file" ]]; then
  curl -L --fail --silent --show-error --connect-timeout 20 --max-time 1800 \
    --speed-limit 10240 --speed-time 60 \
    --continue-at - "$url" -o "$raw_dir/$file.part"
fi

# On interruption rerun this script: a fresh curl invocation recalculates the
# current resume offset. Do not use curl automatic retries with this .part.

python3 codes/52_qc_t16_bbj.py
