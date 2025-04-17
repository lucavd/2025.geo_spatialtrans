test_that("calculate_dropout_probabilities rescales mean_dist correctly", {
  cell_df <- data.frame(boundary_dist = rep(NA, 3))
  mean_dist <- c(1, 2, 3)
  dp <- list(dropout_range = c(0.2, 0.5))
  res <- calculate_dropout_probabilities(cell_df, mean_dist, dropout_params = dp, spatial_params = list(gradient_regions = FALSE))
  # Rescaled from [1,3] to [0.2,0.5]
  expect_equal(res, c(0.2, 0.35, 0.5), tolerance = 1e-8)
})

test_that("calculate_dropout_probabilities uses boundary_dist when gradient_regions TRUE", {
  cell_df <- data.frame(boundary_dist = c(0, 0.5, 1))
  mean_dist <- c(NA, NA, NA)
  dp <- list(dropout_range = c(0.1, 0.4))
  sp <- list(gradient_regions = TRUE)
  res <- calculate_dropout_probabilities(cell_df, mean_dist, dropout_params = dp, spatial_params = sp)
  # dropout = 0.1 + (1 - bd) * 0.3
  expect_equal(res, c(0.4, 0.25, 0.1), tolerance = 1e-8)
})

test_that("apply_dropout without expression_dependent_dropout leaves values unchanged when base_dropout zero", {
  expr <- c(5, 0, 10, 2)
  base_dropout <- rep(0, length(expr))
  params <- list(expression_dependent_dropout = FALSE)
  res <- apply_dropout(expr, base_dropout, dropout_params = params, random_seed = 42)
  expect_equal(res, expr)
})