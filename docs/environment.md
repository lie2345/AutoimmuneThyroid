# Environment specifications and evidence

These are observed stage-specific records and statically discovered requirements, not a tested portable environment or a request to install/upgrade anything.

| Stage/component | Recorded version or requirement | Local evidence path (not distributed as data) |
|---|---|---|
| T06 coloc | `coloc 5.2.3` | `outcome_data/t06/t06_provenance.tsv`, row 2 |
| CC-GWAS | package `0.1.0`; upstream commit `ce9764da946189623a0164f156ad119773bc32f5` | `software/R_library/CCGWAS/DESCRIPTION`; `documentation/analysis/T04B_CCGWAS.md`, lines 43-49 |
| T19 R | `R 4.6.1` | `outputs/t19_t24_20261017/t19/logs/sessionInfo.txt` |
| T19 R packages | data.table `1.18.4`; yaml `2.3.12`; Rsamtools `2.28.0`; GenomicRanges `1.64.0`; IRanges `2.46.0`; digest `0.6.39` | `outputs/t19_t24_20261017/t19/logs/package_versions.tsv` |
| T20 | `R 4.6.1`, base packages; Tabix, GNU gzip/awk | `outputs/t19_t24_20261017/t20/R_sessionInfo.txt`; `run_record.md` |
| Clean-v3 figure runtime | Python `3.12.3`; Matplotlib `3.6.3`; pandas `2.1.4`; NumPy `1.26.4` | clean-v3 `reviews/main_build_validation.json` and `supp_build_validation.json` runtime records |
| Display fonts | DejaVu Sans historical builds; genuine Arial Regular/Bold for final main figures | Figure source and Arial export audit; font files excluded |
| HTSlib | project-local `1.19` build with libcurl for remote Tabix; historical system build lacked HTTPS | historical pipeline guide |
| Python analysis/preparation | Python 3 standard library plus stage-specific dependencies in `python-imports.tsv`; Python 3.12 was used for recent source checks | Python source imports and release validation |

The R packages statically referenced by released source are listed in `r-packages.tsv`; base R packages are included there for transparency. Additional figure/table libraries include Pillow, Matplotlib, NumPy, pandas, ReportLab and openpyxl. The final Arial path uses Playwright/Chromium, pdfplumber, pypdfium2, pypdf and fontTools, plus `fc-list`/fontconfig. No dependency versions are invented for components absent from recorded evidence. Third-party software is obtained separately: [CC-GWAS](https://github.com/wouterpeyrot/CCGWAS), [LAVA](https://github.com/josefin-werme/LAVA), and the project's CBIIT `ldsc39` branch (its exact revision was not established by this export).

No project-wide R `renv.lock`, Python `uv.lock`, or validated standalone environment manifest was available. Installed package libraries and unrelated nested `node_modules` are excluded. The export did not install, restore or modify any environment. Existing local Python checks used `uv run --no-sync --python /usr/bin/python3 python`; R parsing used the installed `Rscript --vanilla`.

## Portable path settings

- `config/analysis.yml`: `ldsc.python` now names `software/ldsc39-venv/bin/python`; the LDSC wrappers also support their existing documented environment overrides. Supply the interpreter and LDSC source separately.
- `config/analysis.yml`: `lava.reference_prefix` now names `raw_data/reference_panels/lava/EUR/EUR`. Supply the correct build/reference files there or configure the actual location after reviewing the method requirements.
- T12 shell scripts: set `LDSC_PYTHON` to an existing compatible interpreter; otherwise the same project-relative location is used.
- Figure export scripts: set `PLAYWRIGHT_PYTHON` to an existing Playwright interpreter, otherwise the invoking interpreter is used.
- Diagram validator: an external helper is expected at `tools/diagram-design/self_check.py`; it is absent and not required to inspect the included frozen final SVGs. Its absence prevents the historical circular builder's and final Arial exporter's complete validation paths.

These changes select paths only. Algorithms, priors, thresholds, sample definitions, coefficients, labels and figure geometry were retained. The exact per-file audit is in `provenance/path_edits.tsv`.
