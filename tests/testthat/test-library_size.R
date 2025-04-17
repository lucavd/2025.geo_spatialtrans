test_that("generate_library_sizes without spatial or cell type effects", {
  # Dummy cell_df
  cell_df <- data.frame(x = 1:10, y = 1:10, intensity_cluster = factor(rep(1, 10)))
  params <- list(
    mean_library_size = 1000,
    library_size_cv = 0.1,
    spatial_effect_on_library = 0,
    cell_type_effect = FALSE
  )
  # No spatial_params needed for this test
  res <- generate_library_sizes(cell_df, library_size_params = params, spatial_params = list(), k_cell_types = 1, random_seed = 123)
  # Check length and type
  expect_equal(length(res), nrow(cell_df))
  expect_type(res, "double")
  # Check mean approximately equal to specified mean
  expect_equal(mean(res), params$mean_library_size, tolerance = params$mean_library_size * 0.1)
})