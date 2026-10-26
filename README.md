# AutoimmuneThyroid

Author-generated analysis, data-processing, validation and figure source for **Correspondence of fine-mapped association signals for Graves’ disease and autoimmune hypothyroidism in FinnGen**, using FinnGen R13 and complementary molecular and external association resources.

Frozen source version: **plos-one-2026-10-26**. Repository: https://github.com/lie2345/AutoimmuneThyroid. No LICENSE file is included in this release.

## What this release contains

- `codes/`: the reported analysis and preparation lineage through script 65, shared R helpers, two author-generated acquisition scripts, and relevant offline test source.
- `config/`: analysis thresholds, trait/cell labels and external-validation contracts.
- `outputs/t19_t24_20261017/t19/code/` and `t20/`: the actual descriptive HLA and exact-variant ophthalmopathy query entrypoints and their checks. Directory names retain their original provenance.
- `outputs/sci_writing_20261017/`: selected author-generated figure/table builders and their historical source dependencies; four final editable main-figure SVGs and their four exact pre-Arial input SVGs.
- `docs/`: ordered entrypoints, environment evidence, data/source layout, limitations and exclusions.
- `provenance/`: original relative paths, original/export SHA-256 values, exact path-only edits, and a separate list of documentation-only unavailable-link repairs. `MANIFEST.sha256` checks the released files except itself.

This is a curated source export, not the original working Git repository. Raw or restricted data, acquired source files, frozen numerical result tables, participant-level data, private correspondence, author metadata, manuscript files, rendered PDFs/TIFFs, external software and fonts are excluded. The TSV files in `config/`, `docs/` and `provenance/` are configuration/inventories, not numerical analysis result tables.

## Reproduction boundary

**The full workflow was not rerun or independently reproduced for this release.** Static Python AST parsing, base-R `parse()` and shell `bash -n` were used without importing or executing production scripts. Source files were compared byte-for-byte, allowing only the listed interpreter/reference-directory changes. The final SVGs are unchanged copies. Successful static checks do not validate scientific results.

An end-to-end run requires separately acquired data, reference panels, external programs, their environments and historical frozen intermediates. Several figure builders deliberately audit old result extracts and baseline artwork before drawing; those tables and binary baselines are not distributed here. A fresh clone is therefore not an independently runnable reproduction of every manuscript artifact. See [known limitations](docs/reproducibility.md) and the [input layout](docs/input-layout.tsv).

Use [entrypoints](docs/entrypoints.md) to select the correct stage. Do not execute every numbered script in sequence: some represent superseded preparation, alternative profiles or gated analyses. Existing scripts may refuse to overwrite frozen outputs; acquisition entrypoints may contact public services and download large inputs when explicitly run. Work in a separate analysis checkout with the required inputs and review each script's CLI and output paths.

## Quick orientation

Run from the repository root. The core wrapper defaults to a command-printing dry run:

```bash
bash codes/run_core_pipeline.sh --dry-run
```

This command does not check data availability and covers only the early core scripts. It was not executed during export. The broader manuscript order is documented in [entrypoints](docs/entrypoints.md). The retained [historical pipeline guide](documentation/analysis/PIPELINE.md) contains earlier-stage instructions and status statements; its five unavailable local links are retained as plain source-reference text marked not distributed. The current export scope and limits in `docs/` take precedence as release documentation.

## Environment

No new environment, package installation or dependency lock was created for this export. The source project has no study-wide `renv.lock` or `uv.lock`; a complete portable environment cannot be reconstructed from the available records. [Environment evidence](docs/environment.md) distinguishes observed stage-specific versions from unpinned requirements. Paths that identified local interpreters or another project were replaced only in export copies; configure the documented LDSC, LAVA and Playwright locations before execution.

## Interpretation retained from the analysis

T06 component-pair evidence, T07 organizational analysis units and variant-level external association checks have different denominators. The 312 released components form 304 organizational units; the 42 component comparisons are not their source denominator. CC-GWAS is a summary-statistic sensitivity analysis. BBJ selected-variant checks do not establish external component-level colocalization. HLA results are descriptive regional associations. Ophthalmopathy endpoints compare cases with population noncases and do not estimate modifiers within Graves' cases. Incomplete or untestable evidence is not negative evidence.

The Graves/AIHT T12 `h2-cts` analyses were partially executed; missing AITD-union and direct case-case inputs keep full T12 on HOLD. T13 does not establish a complete chromatin-to-expression-to-disease chain; T14 is partial contextual validation. Formal MR/PheWAS and T21-T24 target qualification were not implemented and their readiness/contact tools are excluded.
