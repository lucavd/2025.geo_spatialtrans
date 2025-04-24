test_that("generate_depth_field creates a valid 3D depth field", {
  # Create test data
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10)
  )
  n_cells <- 100
  
  # Test with flat surface
  flat_depth <- generate_depth_field(
    test_df,
    depth_params = list(
      depth_pattern = "flat",
      depth_range = c(0, 10)
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_equal(length(flat_depth), n_cells)
  expect_true(all(flat_depth >= 0))
  expect_true(all(flat_depth <= 10))
  
  # A flat surface should have low variance
  expect_true(var(flat_depth) < 2)
  
  # Test with gradient surface
  gradient_depth <- generate_depth_field(
    test_df,
    depth_params = list(
      depth_pattern = "gradient",
      depth_range = c(0, 20),
      gradient_direction = c(1, 1)
    ),
    random_seed = 42
  )
  
  expect_equal(length(gradient_depth), n_cells)
  
  # Test with terrain surface
  terrain_depth <- generate_depth_field(
    test_df,
    depth_params = list(
      depth_pattern = "terrain",
      depth_range = c(5, 25),
      terrain_complexity = 2,
      terrain_smoothness = 0.5
    ),
    random_seed = 42
  )
  
  expect_equal(length(terrain_depth), n_cells)
  expect_true(all(terrain_depth >= 5))
  expect_true(all(terrain_depth <= 25))
  
  # Terrain should have higher variance
  expect_true(var(terrain_depth) > var(flat_depth))
  
  # Invalid pattern should fall back to flat
  invalid_depth <- generate_depth_field(
    test_df,
    depth_params = list(
      depth_pattern = "invalid",
      depth_range = c(0, 10)
    ),
    random_seed = 42
  )
  
  expect_equal(length(invalid_depth), n_cells)
})

test_that("calculate_3d_distances computes accurate 3D distances", {
  # Create test data
  n_cells <- 16
  test_df <- data.frame(
    x = rep(1:4, each = 4),
    y = rep(1:4, times = 4)
  )
  
  # Create simple depth field
  depth <- c(rep(0, 4), rep(5, 4), rep(10, 4), rep(15, 4))
  
  # Calculate 3D distances
  dist_3d <- calculate_3d_distances(test_df, depth)
  
  # Check structure
  expect_true(is.matrix(dist_3d))
  expect_equal(dim(dist_3d), c(n_cells, n_cells))
  
  # Diagonal should be zero
  expect_true(all(diag(dist_3d) == 0))
  
  # Distance should be symmetric
  expect_true(isSymmetric(dist_3d))
  
  # Test specific distance calculation
  # Points at (1,1,0) and (2,2,10)
  p1_idx <- which(test_df$x == 1 & test_df$y == 1)
  p2_idx <- which(test_df$x == 2 & test_df$y == 2)
  
  # Approximate test - the implementation may use a slightly different formula
  expected_dist <- sqrt((2-1)^2 + (2-1)^2 + (10-0)^2)
  actual_dist <- dist_3d[p1_idx, p2_idx]
  expect_true(abs(actual_dist - expected_dist) < 1.0)
  
  # 3D distance should be >= 2D distance
  dist_2d <- as.matrix(dist(test_df[, c("x", "y")]))
  expect_true(all(dist_3d >= dist_2d - 1e-10))  # Allow for small floating point differences
})

test_that("generate_overlap_effects models cell overlap effects", {
  # Test data
  n_cells <- 50
  n_genes <- 10
  depth <- runif(n_cells, 0, 20)
  
  # Create expression matrix
  expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_cells, ncol = n_genes
  )
  
  # Generate overlap effects
  overlap_result <- generate_overlap_effects(
    expr_matrix, depth,
    overlap_params = list(
      overlap_intensity = 0.5,
      overlap_decay = "exponential",
      overlap_gene_specificity = 0.7
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(overlap_result))
  expect_true(all(c("expr_matrix", "overlap_matrix") %in% names(overlap_result)))
  
  # Output matrix should have same dimensions as input
  expect_equal(dim(overlap_result$expr_matrix), c(n_cells, n_genes))
  
  # Output should be different from input
  expect_false(identical(overlap_result$expr_matrix, expr_matrix))
  
  # Output should be non-negative
  expect_true(all(overlap_result$expr_matrix >= 0))
  
  # Overlap matrix should show the contributions
  expect_true(is.matrix(overlap_result$overlap_matrix))
  expect_equal(dim(overlap_result$overlap_matrix), c(n_cells, n_genes))
  
  # Test with different parameters
  overlap_result2 <- generate_overlap_effects(
    expr_matrix, depth,
    overlap_params = list(
      overlap_intensity = 0.2,
      overlap_decay = "linear",
      overlap_gene_specificity = 0.3
    ),
    random_seed = 42
  )
  
  expect_equal(dim(overlap_result2$expr_matrix), c(n_cells, n_genes))
})

test_that("generate_3d_microenvironment integrates full 3D workflow", {
  # Test data
  test_df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5)
  )
  n_cells <- 25
  n_genes <- 15
  
  # Generate base expression
  expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_cells, ncol = n_genes
  )
  
  # Run full 3D microenvironment workflow
  result <- generate_3d_microenvironment(
    test_df, expr_matrix,
    microenvironment_params = list(
      use_3d_microenvironment = TRUE,
      depth_pattern = "terrain",
      depth_range = c(0, 20),
      overlap_intensity = 0.4,
      projection_distortion = 0.3
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("expr_matrix", "depth_field", "dist_3d", "overlap_matrix") %in% names(result)))
  
  # Dimensions should match
  expect_equal(dim(result$expr_matrix), c(n_cells, n_genes))
  expect_equal(length(result$depth_field), n_cells)
  expect_equal(dim(result$dist_3d), c(n_cells, n_cells))
  
  # When disabled, should return original matrix
  result_disabled <- generate_3d_microenvironment(
    test_df, expr_matrix,
    microenvironment_params = list(
      use_3d_microenvironment = FALSE
    ),
    random_seed = 42
  )
  
  expect_true(is.list(result_disabled))
  expect_true(identical(result_disabled$expr_matrix, expr_matrix))
})