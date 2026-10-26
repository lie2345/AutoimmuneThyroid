# Analysis source directory

See the repository [entrypoint map](../docs/entrypoints.md). Numeric prefixes preserve the working-project filenames and do not imply that every script should run automatically.

`lib/common.R` and `config/analysis.yml` centralize shared paths and thresholds. `run_core_pipeline.sh` prints its early-core commands by default; external LD-dependent stages remain separately gated. Scripts `00`-`65` are retained with the original numbering, including duplicate prefixes, because subsequent preparation and validation stages use those exact names. Relevant `test_*.py` files are distributed as source and were parsed, not executed, for this export.

Acquisition scripts have network and storage effects when run. The two `download_*` scripts now derive the project root from their own location. The T12 LDSC shells accept `LDSC_PYTHON` and otherwise use `software/ldsc39-venv/bin/python` within the analysis checkout. Scientific code, thresholds and algorithms are retained.
