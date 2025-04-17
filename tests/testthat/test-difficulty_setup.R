test_that("configure_difficulty_level sets parameters for easy level", {
  cfg <- configure_difficulty_level(
    difficulty_level = "easy"
  )
  expect_equal(cfg$difficulty_level, "easy")
  expect_equal(cfg$marker_params$marker_genes_per_type, 10)
  expect_equal(cfg$marker_params$marker_expression_fold, 2.0)
  expect_equal(cfg$marker_params$marker_overlap_fold, 0.0)
  expect_equal(cfg$spatial_params$spatial_noise_intensity, 0.5)
  expect_equal(cfg$spatial_params$spatial_range, 50)
  expect_equal(cfg$dropout_params$dropout_range, c(0.1, 0.3))
  expect_equal(cfg$dropout_params$dispersion_range, c(3.0, 1.5))
  expect_equal(cfg$cell_specific_params$cell_specific_noise_sd, 0.1)
})

test_that("configure_difficulty_level warns and defaults invalid level to medium", {
  expect_warning(
    cfg <- configure_difficulty_level(
      difficulty_level = "invalid"
    ),
    "Livello di difficoltà non valido"
  )
  expect_equal(cfg$difficulty_level, "medium")
})