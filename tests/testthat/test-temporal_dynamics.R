test_that("generate_pseudotime_field creates valid pseudotime gradients", {
  # Create test data
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10)
  )
  
  # Test gradient mode
  pt_gradient <- generate_pseudotime_field(
    test_df, 
    temporal_params = list(
      pseudotime_mode = "gradient",
      pseudotime_origin = c(1, 1)
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_equal(length(pt_gradient), 100)
  expect_true(all(pt_gradient >= 0 & pt_gradient <= 1))
  
  # Values should have a gradient structure
  # Further from origin should have higher pseudotime
  origin_idx <- which(test_df$x == 1 & test_df$y == 1)
  far_idx <- which(test_df$x == 10 & test_df$y == 10)
  expect_true(pt_gradient[far_idx] > pt_gradient[origin_idx])
  
  # Test focal mode
  pt_focal <- generate_pseudotime_field(
    test_df, 
    temporal_params = list(
      pseudotime_mode = "focal",
      n_foci = 3,
      focal_radius = 0.2
    ),
    random_seed = 42
  )
  
  expect_equal(length(pt_focal), 100)
  expect_true(all(pt_focal >= 0 & pt_focal <= 1))
  
  # Test bifurcation mode
  pt_bifurcation <- generate_pseudotime_field(
    test_df, 
    temporal_params = list(
      pseudotime_mode = "bifurcation",
      bifurcation_point = c(5, 5),
      branch_angles = c(45, 135)
    ),
    random_seed = 42
  )
  
  expect_equal(length(pt_bifurcation), 100)
  expect_true(all(pt_bifurcation >= 0 & pt_bifurcation <= 1))
  
  # Invalid mode should fall back to gradient
  pt_invalid <- generate_pseudotime_field(
    test_df, 
    temporal_params = list(
      pseudotime_mode = "invalid"
    ),
    random_seed = 42
  )
  
  expect_equal(length(pt_invalid), 100)
})

test_that("generate_gene_trajectories creates valid gene-pseudotime patterns", {
  # Create pseudotime vector
  n_cells <- 100
  pseudotime <- seq(0, 1, length.out = n_cells)
  n_genes <- 20
  
  # Test transient pattern
  traj_transient <- generate_gene_trajectories(
    pseudotime, n_genes,
    trajectory_params = list(
      pattern_distribution = c(transient = 1),
      trajectory_smoothness = 0.1
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.matrix(traj_transient))
  expect_equal(dim(traj_transient), c(n_cells, n_genes))
  
  # Test mixed patterns
  traj_mixed <- generate_gene_trajectories(
    pseudotime, n_genes,
    trajectory_params = list(
      pattern_distribution = c(
        monotonic = 0.25,
        transient = 0.25, 
        cyclic = 0.25,
        bifurcating = 0.25
      ),
      trajectory_smoothness = 0.1
    ),
    random_seed = 42
  )
  
  expect_true(is.matrix(traj_mixed))
  expect_equal(dim(traj_mixed), c(n_cells, n_genes))
  
  # Test with custom gene modules
  gene_modules <- list(
    module1 = 1:5,
    module2 = 6:10,
    module3 = 11:15,
    module4 = 16:20
  )
  
  traj_modules <- generate_gene_trajectories(
    pseudotime, n_genes,
    trajectory_params = list(
      pattern_distribution = c(monotonic = 0.5, transient = 0.5),
      trajectory_smoothness = 0.1,
      use_gene_modules = TRUE
    ),
    gene_modules = gene_modules,
    random_seed = 42
  )
  
  expect_true(is.matrix(traj_modules))
  expect_equal(dim(traj_modules), c(n_cells, n_genes))
  
  # Check that genes within same module have correlated patterns
  module_cor <- cor(traj_modules[, 1:5])
  expect_true(all(module_cor > 0.5))
})

test_that("generate_rna_velocity creates valid velocity vectors", {
  # Create test data
  n_cells <- 50
  n_genes <- 10
  pseudotime <- seq(0, 1, length.out = n_cells)
  
  # Create gene expression trajectory
  expr_matrix <- matrix(
    runif(n_cells * n_genes, 0, 10),
    nrow = n_cells, ncol = n_genes
  )
  
  # Generate velocities
  velocity_result <- generate_rna_velocity(
    expr_matrix, pseudotime,
    velocity_params = list(
      velocity_strength = 1.0,
      velocity_noise = 0.2
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(velocity_result))
  expect_true(all(c("velocity", "unspliced") %in% names(velocity_result)))
  
  expect_true(is.matrix(velocity_result$velocity))
  expect_true(is.matrix(velocity_result$unspliced))
  
  expect_equal(dim(velocity_result$velocity), c(n_cells, n_genes))
  expect_equal(dim(velocity_result$unspliced), c(n_cells, n_genes))
  
  # Test with velocity disabled
  velocity_disabled <- generate_rna_velocity(
    expr_matrix, pseudotime,
    velocity_params = list(
      velocity_strength = 0
    ),
    random_seed = 42
  )
  
  # Should return all zeros for velocity
  expect_true(all(velocity_disabled$velocity == 0))
})

test_that("generate_temporal_dynamics integrates full temporal workflow", {
  # Create test data
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
  
  # Run full temporal workflow
  result <- generate_temporal_dynamics(
    test_df, expr_matrix,
    temporal_params = list(
      use_temporal_dynamics = TRUE,
      pseudotime_mode = "gradient",
      temporal_gene_fraction = 0.8,
      pattern_distribution = c(monotonic = 0.7, transient = 0.3),
      trajectory_strength = 0.8,
      include_velocity = TRUE
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("expr_matrix", "pseudotime", "velocity", "unspliced") %in% names(result)))
  
  # Dimensions should match
  expect_equal(dim(result$expr_matrix), c(n_cells, n_genes))
  expect_equal(length(result$pseudotime), n_cells)
  expect_equal(dim(result$velocity), c(n_cells, n_genes))
  
  # Values should be realistic
  expect_true(all(result$pseudotime >= 0 & result$pseudotime <= 1))
  expect_true(all(result$expr_matrix >= 0))
  
  # When disabled, should return original matrix
  result_disabled <- generate_temporal_dynamics(
    test_df, expr_matrix,
    temporal_params = list(
      use_temporal_dynamics = FALSE
    ),
    random_seed = 42
  )
  
  expect_true(is.list(result_disabled))
  expect_true(identical(result_disabled$expr_matrix, expr_matrix))
  expect_true(is.numeric(result_disabled$pseudotime))
  expect_true(all(result_disabled$velocity == 0))
})