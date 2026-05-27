library(shinytest2)

# ── Helpers ──────────────────────────────────────────────────────────────────

make_app <- function(name) {
  AppDriver$new(
    app_dir      = system.file("shiny", package = "scFOCAL"),
    name         = name,
    load_timeout = 120000,
    timeout      = 120000
  )
}

# Build a drugs x cells corrMat CSV compatible with the fixture Seurat object
make_corr_mat_csv <- function(drug_names = c("DrugA", "DrugB", "DrugC")) {
  seurat_obj <- readRDS(test_path("fixtures/downsampled_seuratObj.RDS"))
  cells <- colnames(seurat_obj)
  mat <- matrix(
    runif(length(drug_names) * length(cells)),
    nrow     = length(drug_names),
    dimnames = list(drug_names, cells)
  )
  path <- tempfile(fileext = ".csv")
  write.csv(as.data.frame(mat), path)
  path
}

# Build a custom drug signature CSV with lowercase gene names (mismatch test)
make_lowercase_drug_sig_csv <- function() {
  seurat_obj <- readRDS(test_path("fixtures/downsampled_seuratObj.RDS"))
  genes <- rownames(seurat_obj)
  drugs <- c("Drug1", "Drug2")
  mat <- matrix(
    rnorm(length(drugs) * length(genes)),
    nrow     = length(drugs),
    dimnames = list(drugs, tolower(genes))
  )
  path <- tempfile(fileext = ".csv")
  write.csv(as.data.frame(mat), path)
  path
}

# ── Bug 4: app navigated before Seurat upload should not crash ───────────────

test_that("app does not crash when accessed before uploading Seurat file", {
  # assay_to_use() calls isOrthogonAL() which calls req(input$seurobjRDS)
  # before a file is uploaded. The req() should silently pause, not crash.
  app <- make_app("pre-upload-nav")

  expect_no_error(app$get_value(output = "debugRelease"))

  app$stop()
})

# ── Bug 1: corrMat double transpose — compound dropdown must show drug names ──

test_that("corrMat upload shows drug names in compound dropdown, not cell barcodes", {
  # If the double-transpose bug is present, rownames(RDS_Final_CorrMat())
  # returns cell barcodes instead of drug names, populating the compound
  # dropdown with barcodes.
  drug_names    <- c("DrugA", "DrugB", "DrugC")
  corr_mat_path <- make_corr_mat_csv(drug_names)

  app <- make_app("corr-mat-orientation")

  app$upload_file(seurobjRDS = test_path("fixtures/downsampled_seuratObj.RDS"))
  app$wait_for_value(output = "seuratLoaded", timeout = 30000)

  app$set_inputs(uploadCorrelationMatrix = TRUE)
  app$upload_file(corrMatUpload = corr_mat_path)
  # Wait for the server to confirm the upload registered
  app$wait_for_value(output = "corrMatUploaded", timeout = 10000)

  # The Reference_L1000_or_Custom radio button lives inside a conditionalPanel
  # that is hidden until corrMat is uploaded, so set it explicitly after upload.
  app$set_inputs(Reference_L1000_or_Custom = "L1000 Derived")

  # wait_for_value blocks until the output is actually sent by Shiny,
  # avoiding the race condition with suspendWhenHidden output evaluation.
  compound_ui <- app$wait_for_value(output = "referenceCompound_ui", timeout = 15000)

  # The duplicate uiOutput ID returns a length-2 vector; collapse to one string.
  compound_html <- paste(compound_ui, collapse = "\n")

  for (drug in drug_names) {
    expect_true(
      grepl(drug, compound_html, fixed = TRUE),
      label = paste("Compound dropdown should contain", drug, "not cell barcodes")
    )
  }

  app$stop()
})

# ── Bug 2: duplicate uiOutput — compound dropdown must not be empty ───────────

test_that("reference compound dropdown is populated after corrMat upload", {
  # uiOutput('referenceCompound_ui') appears in two tabs in ui.R.
  # When the same output ID is bound twice, the second binding can
  # overwrite the first, leaving one tab's dropdown empty.
  drug_names    <- c("DrugX", "DrugY")
  corr_mat_path <- make_corr_mat_csv(drug_names)

  app <- make_app("duplicate-ui-output")

  app$upload_file(seurobjRDS = test_path("fixtures/downsampled_seuratObj.RDS"))
  app$wait_for_value(output = "seuratLoaded", timeout = 30000)

  app$set_inputs(uploadCorrelationMatrix = TRUE)
  app$upload_file(corrMatUpload = corr_mat_path)
  # Wait for the server to confirm the upload registered
  app$wait_for_value(output = "corrMatUploaded", timeout = 10000)

  # The Reference_L1000_or_Custom radio button lives inside a conditionalPanel
  # that is hidden until corrMat is uploaded, so set it explicitly after upload.
  app$set_inputs(Reference_L1000_or_Custom = "L1000 Derived")

  # wait_for_value blocks until the output is actually sent by Shiny,
  # avoiding the race condition with suspendWhenHidden output evaluation.
  compound_ui <- app$wait_for_value(output = "referenceCompound_ui", timeout = 15000)

  # The duplicate uiOutput ID returns a length-2 vector; at least one must be non-empty.
  expect_true(
    any(nchar(trimws(compound_ui)) > 0),
    label = "referenceCompound_ui should not be empty after corrMat upload"
  )

  app$stop()
})

# ── Bug 5: lowercase gene names in custom drug signature should not crash ─────

test_that("custom drug signature with lowercase gene names does not crash app", {
  # If gene names in the uploaded CSV are lowercase (e.g. tp53) but the
  # Seurat object uses uppercase (TP53), overlap detection returns 0 and
  # the server calls stop(), crashing the reactive chain silently.
  # The app should survive — either by normalising case or showing an error.
  sig_path <- make_lowercase_drug_sig_csv()

  app <- make_app("gene-case-mismatch")

  app$upload_file(seurobjRDS = test_path("fixtures/downsampled_seuratObj.RDS"))
  app$wait_for_value(output = "seuratLoaded", timeout = 30000)

  app$set_inputs(L1000_Release = "Custom Upload")
  app$upload_file(customL1000Upload = sig_path)
  app$wait_for_idle(timeout = 15000)

  # If the app crashed, reading any output would throw — this proves it's alive
  expect_no_error(app$get_value(output = "debugRelease"))

  app$stop()
})
