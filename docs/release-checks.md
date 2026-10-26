# Release preparation checks

Date: 2026-10-26. Status: **PASS_STATIC_EXPORT_CHECKS**.

- 120 source-file records checked against live source and exported SHA-256 values; code/SVG copies reproduce the originals after only the declared path replacements. One historical documentation file also has five separately recorded unavailable-link repairs.
- 41 Python files parsed with `ast.parse()`; no production modules imported.
- 50 R files parsed with installed `Rscript --vanilla` and base-R `parse()`; no source evaluated.
- 11 shell files checked with `bash -n`.
- 8 SVGs parsed as XML and checked as exact source copies; no figure was rerendered.
- Configuration JSON, local documentation links, selected code dependency closure, file-type exclusions and credential/contact/private-path patterns checked.
- 456/456 recorded static/export checks passed. No credential values are included in screen records.

Local preparation commands were `uv run --no-sync --python /usr/bin/python3 python outputs/code_release_20261026/build_export.py` and `uv run --no-sync --python /usr/bin/python3 python outputs/code_release_20261026/validate_export.py`. Release-preparation helpers and private-path audit are outside this public source tree. Checksums can be checked after download with `sha256sum -c MANIFEST.sha256` from the repository root.

Not run: production analyses, acquisition/downloads, test suites, numerical/integration checks, complete pipeline reproduction, figure rendering or independent scientific review. Static checks do not establish scientific validity or complete runtime reproducibility. The absent external diagram helper used by the historical circular builder and final Arial exporter is documented in `reproducibility.md`; data/result/baseline-artwork dependencies remain absent by design.
