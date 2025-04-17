test_that("calculate_spatial_distances returns correct distances and densities", {
  # Simple square grid with two clusters
  cell_df <- data.frame(
    x = c(0, 1, 0, 1),
    y = c(0, 0, 1, 1),
    intensity_cluster = factor(c(1, 1, 2, 2))
  )
  # Compute distances and densities
  res <- calculate_spatial_distances(cell_df, chunk_size = 2, random_seed = 42)
  # Check output structure
  expect_type(res, "list")
  expect_equal(dim(res$dist_mat), c(4L, 4L))
  # All same-cluster distances mean should be (0 + 1)/2 = 0.5
  expect_equal(unname(res$mean_dist), rep(0.5, 4), tolerance = 1e-8)
  # Local density: only self-distance < quantile, so 1/4
  expect_equal(unname(res$local_density), rep(0.25, 4), tolerance = 1e-8)
})