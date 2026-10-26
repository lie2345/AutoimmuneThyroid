#!/usr/bin/env bash
# Resume the exact unrestricted official archive; never overwrite a final archive.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
dest='raw_data/external_gwas/bbj_2019_graves/hum0014.v17.GD.v1.zip'
mkdir -p raw_data/external_gwas/bbj_2019_graves
if [[ ! -f "$dest" ]]; then
  curl -L --fail --silent --show-error --connect-timeout 20 --max-time 300 \
    --speed-limit 10240 --speed-time 30 \
    --retry 4 --retry-all-errors --retry-delay 3 --continue-at - \
    'https://humandbs.dbcls.jp/files/hum0014/hum0014.v17.GD.v1.zip' \
    -o "${dest}.part"
fi
python3 codes/50_preflight_t17_bbj.py
