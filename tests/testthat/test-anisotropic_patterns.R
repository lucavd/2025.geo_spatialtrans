test_that("generate_backbone_structures creates valid spatial structures", {
  # Create test data
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10)
  )
  n_cells <- 100
  
  # Test with linear structure
  linear_result <- generate_backbone_structures(
    test_df,
    structure_params = list(
      n_structures = 1,
      structure_type = "linear",
      structure_width = 2,
      structure_length_factor = 0.8
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(linear_result))
  expect_true(all(c("distance_matrices", "structure_mask") %in% names(linear_result)))
  
  # Should have one distance matrix per structure
  expect_length(linear_result$distance_matrices, 1)
  expect_equal(dim(linear_result$distance_matrices[[1]]), c(n_cells, n_cells))
  
  # Structure mask should identify cells inside structures
  expect_equal(length(linear_result$structure_mask), n_cells)
  expect_true(any(linear_result$structure_mask > 0))
  
  # Test with branched structures
  branched_result <- generate_backbone_structures(
    test_df,
    structure_params = list(
      n_structures = 2,
      structure_type = "branched",
      n_branches = 3,
      structure_width = 1.5
    ),
    random_seed = 42
  )
  
  expect_length(branched_result$distance_matrices, 2)
  expect_true(any(branched_result$structure_mask > 0))
  
  # Test with network structures
  network_result <- generate_backbone_structures(
    test_df,
    structure_params = list(
      n_structures = 1,
      structure_type = "network",
      network_density = 0.1,
      structure_width = 1
    ),
    random_seed = 42
  )
  
  expect_length(network_result$distance_matrices, 1)
  expect_true(any(network_result$structure_mask > 0))
})

test_that("generate_anisotropic_expression creates expression patterns along structures", {
  # Create test data
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10)
  )
  n_cells <- 100
  n_genes <- 10
  n_structures <- 2
  
  # Create mock distance matrices and structure mask
  distance_matrices <- list()
  for (i in 1:n_structures) {
    distance_matrices[[i]] <- matrix(runif(n_cells * n_cells, 0, 10), nrow = n_cells)
  }
  
  structure_mask <- rep(0, n_cells)
  structure_mask[sample(1:n_cells, 30)] <- sample(1:n_structures, 30, replace = TRUE)
  
  # Test gradient pattern
  gradient_result <- generate_anisotropic_expression(
    test_df,
    distance_matrices,
    structure_mask,
    n_genes,
    expression_params = list(
      anisotropic_pattern = "gradient",
      anisotropic_gene_fraction = 0.5,
      pattern_smoothness = 0.2
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.matrix(gradient_result))
  expect_equal(dim(gradient_result), c(n_cells, n_genes))
  
  # Test oscillating pattern
  oscillating_result <- generate_anisotropic_expression(
    test_df,
    distance_matrices,
    structure_mask,
    n_genes,
    expression_params = list(
      anisotropic_pattern = "oscillating",
      anisotropic_gene_fraction = 0.8,
      oscillation_frequency = 0.2
    ),
    random_seed = 42
  )
  
  expect_true(is.matrix(oscillating_result))
  expect_equal(dim(oscillating_result), c(n_cells, n_genes))
  
  # Test hotspot pattern
  hotspot_result <- generate_anisotropic_expression(
    test_df,
    distance_matrices,
    structure_mask,
    n_genes,
    expression_params = list(
      anisotropic_pattern = "hotspot",
      n_hotspots_per_structure = 2,
      hotspot_radius = 5
    ),
    random_seed = 42
  )
  
  expect_true(is.matrix(hotspot_result))
  expect_equal(dim(hotspot_result), c(n_cells, n_genes))
})

test_that("apply_anisotropic_effects modifies expression based on anisotropic patterns", {
  # Test data
  n_cells <- 50
  n_genes <- 20
  
  # Create base expression
  expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_cells, ncol = n_genes
  )
  
  # Create anisotropic pattern matrix
  aniso_matrix <- matrix(
    runif(n_cells * n_genes, -1, 1),
    nrow = n_cells, ncol = n_genes
  )
  
  # Create structure mask
  structure_mask <- rep(0, n_cells)
  structure_mask[sample(1:n_cells, 20)] <- 1
  
  # Apply effects - multiplicative
  result_mult <- apply_anisotropic_effects(
    expr_matrix,
    aniso_matrix,
    structure_mask,
    effect_params = list(
      anisotropic_effect_type = "multiplicative",
      anisotropic_effect_strength = 0.8,
      background_effect_fraction = 0.2
    )
  )
  
  # Check structure
  expect_true(is.matrix(result_mult))
  expect_equal(dim(result_mult), c(n_cells, n_genes))
  
  # Should be different from original
  expect_false(identical(result_mult, expr_matrix))
  
  # Apply effects - additive
  result_add <- apply_anisotropic_effects(
    expr_matrix,
    aniso_matrix,
    structure_mask,
    effect_params = list(
      anisotropic_effect_type = "additive",
      anisotropic_effect_strength = 0.5,
      background_effect_fraction = 0.1
    )
  )
  
  expect_true(is.matrix(result_add))
  expect_equal(dim(result_add), c(n_cells, n_genes))
  
  # Result should be non-negative
  expect_true(all(result_add >= 0))
})

test_that("generate_anisotropic_patterns integrates full anisotropic workflow", {
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
  
  # Run full anisotropic workflow
  result <- generate_anisotropic_patterns(
    test_df, expr_matrix,
    anisotropic_params = list(
      use_anisotropic_patterns = TRUE,
      n_structures = 2,
      structure_type = "linear",
      anisotropic_pattern = "gradient",
      anisotropic_gene_fraction = 0.6,
      anisotropic_effect_strength = 0.8
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("expr_matrix", "structure_mask", "distance_matrices") %in% names(result)))
  
  # Dimensions should match
  expect_equal(dim(result$expr_matrix), c(n_cells, n_genes))
  expect_equal(length(result$structure_mask), n_cells)
  
  # There should be structures identified
  expect_true(any(result$structure_mask > 0))
  
  # When disabled, should return original matrix
  result_disabled <- generate_anisotropic_patterns(
    test_df, expr_matrix,
    anisotropic_params = list(
      use_anisotropic_patterns = FALSE
    ),
    random_seed = 42
  )
  
  expect_true(is.list(result_disabled))
  expect_true(identical(result_disabled$expr_matrix, expr_matrix))
  expect_true(all(result_disabled$structure_mask == 0))
})