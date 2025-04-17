test_that("save_simulation_results and load_simulation_results round-trip correctly", {
  tmp <- tempfile(fileext = ".rds")
  obj <- list(a = 1, b = "test", c = 3.14)
  # Save
  path <- save_simulation_results(obj, tmp)
  expect_true(file.exists(path))
  # Load
  loaded <- load_simulation_results(tmp)
  expect_equal(loaded, obj)
})

test_that("load_simulation_results errors on missing file", {
  expect_error(
    load_simulation_results("nonexistent_file.rds"),
    "File non trovato"
  )
})