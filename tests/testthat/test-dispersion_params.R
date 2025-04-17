test_that("calculate_dispersion_params rescales mean_dist correctly without cell type effect", {
  # Prepare a dummy cell_df without boundary_dist
  cell_df <- data.frame(intensity_cluster = factor(c(1,2,3)))
  mean_dist <- c(1, 2, 3)
  # No cell type effect to simplify
  dp <- list(dispersion_range = c(2, 1), cell_type_dispersion_effect = 0)
  res <- calculate_dispersion_params(cell_df, mean_dist, dropout_params = dp, k_cell_types = 3, random_seed = 1)
  # Rescaled from [1,3] to [2,1]: expect c(2,1.5,1)
  expect_equal(res, c(2, 1.5, 1), tolerance = 1e-8)
})

test_that("calculate_dispersion_params uses boundary_dist when provided", {
  # Create cell_df with boundary_dist
  cell_df <- data.frame(
    intensity_cluster = factor(1:3),
    boundary_dist = c(0, 0.5, 1)
  )
  dp <- list(dispersion_range = c(2, 1), cell_type_dispersion_effect = 0)
  res <- calculate_dispersion_params(cell_df, NULL, dropout_params = dp, k_cell_types = 3, random_seed = 1)
  # dispersion = 1 + boundary_dist * (2 - 1) = 1 + bd
  expect_equal(res, c(1, 1.5, 2), tolerance = 1e-8)
})

test_that("calculate_dispersion_params warns on invalid dispersion_range and uses default", {
  cell_df <- data.frame(intensity_cluster = factor(1:4))
  bad_dp <- list(dispersion_range = c(1), cell_type_dispersion_effect = 0)
  expect_warning(
    res <- calculate_dispersion_params(cell_df, mean_dist = NULL, dropout_params = bad_dp, k_cell_types = 2, random_seed = 1),
    "dispersion_range deve essere un vettore numerico di lunghezza 2"
  )
  # Should use default c(2,1) -> mean = 1.5 for all
  expect_equal(res, rep(1.5, 4), tolerance = 1e-8)
})