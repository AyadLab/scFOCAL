# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```r
# Rebuild package after any change to inst/shiny/ or R/
devtools::install()          # or Build → Clean and Rebuild in RStudio

# Launch app manually
scFOCAL::runscFOCAL()

# Run all shinytest2 tests
shinytest2::test_app("inst/shiny")

# Run a single test file
testthat::test_file("inst/shiny/tests/testthat/test-upload.R")
testthat::test_file("inst/shiny/tests/testthat/test-edge-cases.R")
```

**Critical**: Tests run against the *installed* package via `system.file("shiny", package = "scFOCAL")`. Always do a Clean and Rebuild before running tests after editing anything in `inst/shiny/` or `R/`.

## Architecture

**Entry point**: `R/runApp.R` — `runscFOCAL()` calls `shiny::shinyAppDir(system.file('shiny', package = "scFOCAL"))`.

The Shiny app lives in `inst/shiny/server.R` + `inst/shiny/ui.R`.

### Key server.R reactives

| Reactive | Purpose |
|---|---|
| `rdsSeurat()` | Loads uploaded Seurat RDS; handles Seurat v5 layer detection |
| `seuratLoaded` | Boolean output; used by conditional panels and tests as upload health-check |
| `assay_to_use()` | Selects RNA vs RNA_ortho assay; guards with `req(input$seurobjRDS)` |
| `drugSignatures()` | Accepts custom drugs×genes CSV; auto-detects orientation via gene-name overlap; falls back to built-in LINCS |
| `corrMatUpload()` | Reads uploaded drugs×cells correlation matrix CSV |
| `LINCS.ResponseSigs()` / `L1000_genes` | Built-in L1000 signature data from package namespace |
| `diseaseSigs` | Triggered by `CalcDiseaseSig` action button |
| `debugRelease` | Outputs current `L1000_Release` input value; used as app alive-check in tests |

### UI tab structure

Overview → **Run scFOCAL-dev** (tabsetPanel):
1. Data Upload — Seurat RDS, assay selection
2. Exploration — dim plots, cell highlights
3. Disease Signatures — cell subsetting, heatmap/table/reversal tabs
4. Drug Response / Synergy Analysis — corrMat upload, compound selection

**Duplicate ID note**: `corrMatUpload` (input) and `referenceCompound_ui` (output) intentionally appear in two tabs — this is existing architecture. Shiny handles duplicate `uiOutput` IDs correctly at runtime. In shinytest2 tests, navigate to the correct tab *before* uploading or reading these elements to avoid "Multiple HTML elements found" warnings.

## Testing Infrastructure

Tests live in `inst/shiny/tests/testthat/`:

- `setup-chromote.R` — raises Chromote's `default_timeout` to 120s; required because the app takes ~12s to load packages on a cold start, which exceeds Chrome's default 10s `Page.navigate` timeout
- `test-upload.R` — app launch and Seurat upload health checks
- `test-edge-cases.R` — regression tests for corrMat orientation, gene case normalization, pre-upload navigation crash
- `fixtures/downsampled_seuratObj.RDS` — 134MB Seurat fixture; requires **Git LFS** before pushing to GitHub

CI: `.github/workflows/shinytest2.yml` runs on every PR to `drug_signature_update`.

## Open Bugs (as of drug_signature_update branch)

- **Bug 1** — `corrMatUpload()` in server.R double-transposes the matrix, swapping drug names and cell barcodes
- **Bug 2** — Shiny warns about duplicate `corrMatUpload`/`referenceCompound_ui` IDs in tests; workaround by navigating to the correct tab first
- **Bug 4** — `assay_to_use()` may be called before Seurat upload; needs `req(rdsSeurat())` guard
- **Bug 5** — Uploaded drug signature CSVs with lowercase gene names produce 0-gene overlap with the Seurat object (which uses uppercase), causing a silent crash; needs `toupper()` normalization
