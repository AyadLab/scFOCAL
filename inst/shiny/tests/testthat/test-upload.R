library(shinytest2)

test_that("app loads without errors", {
  app <- AppDriver$new(
    app_dir      = system.file("shiny", package = "scFOCAL"),
    name         = "app-loads",
    load_timeout = 120000,
    timeout      = 120000
  )

  # Confirm the app started — reading any output proves it didn't crash
  expect_no_error(app$get_value(output = "debugRelease"))

  app$stop()
})

test_that("uploading a valid Seurat RDS shows success message", {
  app <- AppDriver$new(
    app_dir      = system.file("shiny", package = "scFOCAL"),
    name         = "seurat-upload",
    load_timeout = 120000,
    timeout      = 120000
  )

  app$upload_file(
    seurobjRDS = test_path("fixtures/downsampled_seuratObj.RDS")
  )

  # Wait for seuratLoaded to be TRUE. Ignore NULL (pre-upload) and any
  # transient FALSE that may arrive while the 134 MB file is being read.
  loaded <- app$wait_for_value(
    output  = "seuratLoaded",
    ignore  = list(NULL, FALSE),
    timeout = 30000
  )

  expect_true(isTRUE(loaded))

  app$stop()
})
