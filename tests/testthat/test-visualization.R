test_that("create_simulation_plots returns a ggplot object", {
  cell_df <- data.frame(
    x = 1:4,
    y = 4:1,
    intensity_cluster = factor(c("1", "2", "1", "2"))
  )
  config <- list(
    grid_mode = TRUE,
    grid_resolution = 2,
    n_genes = 3,
    threshold_value = 0.5
  )
  difficulty <- list(
    difficulty_level = "easy",
    marker_params = list(
      marker_genes_per_type = 1,
      marker_expression_fold = 2.0
    )
  )
  p <- create_simulation_plots(cell_df, config, difficulty)
  # Check that p is a ggplot object
  expect_s3_class(p, "ggplot")
  # Should contain tile layer for grid_mode=TRUE
  layers <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomTile" %in% layers)
})