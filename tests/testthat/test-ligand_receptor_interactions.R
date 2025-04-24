test_that("define_ligand_receptor_db creates a valid database", {
  # Test with default parameters
  n_genes <- 100
  lr_db <- define_ligand_receptor_db(n_genes = n_genes, random_seed = 42)
  
  # Check structure and content
  expect_true(is.list(lr_db))
  expect_true("interactions" %in% names(lr_db))
  expect_true("gene_roles" %in% names(lr_db))
  
  interactions <- lr_db$interactions
  expect_true(is.list(interactions))
  expect_equal(length(interactions), 20)  # Default n_interactions
  
  # Check first interaction structure
  first_interaction <- interactions[[1]]
  expect_true(is.list(first_interaction))
  expect_true(all(c("ligand_id", "receptor_id", "target_ids", "target_effects") %in% names(first_interaction)))
  
  # Check gene roles
  gene_roles <- lr_db$gene_roles
  expect_equal(length(gene_roles), n_genes)
  expect_true(all(gene_roles %in% c("ligand", "receptor", "other")))
  
  # Test with custom parameters
  custom_params <- list(
    n_interactions = 10,
    min_target_genes = 2,
    max_target_genes = 8,
    effect_strength = c(0.2, 1.5),
    inhibitory_prob = 0.4
  )
  
  lr_db2 <- define_ligand_receptor_db(
    n_genes = n_genes,
    interaction_params = custom_params,
    random_seed = 42
  )
  
  expect_equal(length(lr_db2$interactions), 10)
})

test_that("compute_lr_signaling_effects calculates distance-based effects", {
  # Create test data
  n_cells <- 20
  n_genes <- 30
  
  test_df <- data.frame(
    x = rep(1:5, each = 4),
    y = rep(1:4, times = 5)
  )
  
  # Create expression matrix
  expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_cells, ncol = n_genes
  )
  
  # Create a simplified interaction database
  interactions <- list()
  interactions[[1]] <- list(
    ligand_id = 1,
    ligand_name = "LIG_1",
    receptor_id = 2,
    receptor_name = "REC_2",
    target_ids = 3:7,
    target_effects = c(0.8, -0.6, 0.5, 0.7, -0.4),
    decay_factor = 20,
    activation_threshold = 0.1,
    interaction_type = "activating"
  )
  
  interactions[[2]] <- list(
    ligand_id = 10,
    ligand_name = "LIG_10",
    receptor_id = 15,
    receptor_name = "REC_15",
    target_ids = 20:25,
    target_effects = c(0.9, 0.8, -0.7, 0.5, 0.6, -0.5),
    decay_factor = 30,
    activation_threshold = 0.2,
    interaction_type = "inhibitory"
  )
  
  lr_db <- list(
    interactions = interactions,
    gene_roles = rep("other", n_genes)
  )
  lr_db$gene_roles[c(1, 10)] <- "ligand"
  lr_db$gene_roles[c(2, 15)] <- "receptor"
  
  # Compute distances
  dist_mat <- as.matrix(dist(test_df[, c("x", "y")]))
  
  # Run function with default parameters
  result <- compute_lr_signaling_effects(
    test_df,
    expr_matrix,
    lr_db,
    dist_mat,
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n_cells, n_genes))
  
  # Test with custom parameters
  custom_params <- list(
    signal_propagation_mode = "threshold",
    max_signaling_distance = 50,
    signal_amplification = 1.5
  )
  
  result2 <- compute_lr_signaling_effects(
    test_df,
    expr_matrix,
    lr_db,
    dist_mat,
    interaction_params = custom_params,
    random_seed = 42
  )
  
  expect_true(is.matrix(result2))
  expect_equal(dim(result2), c(n_cells, n_genes))
})

test_that("apply_lr_signaling applies signaling effects to expression", {
  # Create test data
  n_cells <- 15
  n_genes <- 10
  
  # Base expression matrix
  expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_cells, ncol = n_genes
  )
  
  # Signaling effects matrix
  effect_matrix <- matrix(
    runif(n_cells * n_genes, -1, 1),
    nrow = n_cells, ncol = n_genes
  )
  
  # Apply effects with default parameters
  result <- apply_lr_signaling(
    expr_matrix, 
    effect_matrix,
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n_cells, n_genes))
  
  # Should be different from original
  expect_false(identical(result, expr_matrix))
  
  # Values should be non-negative
  expect_true(all(result >= 0))
  
  # Apply effects with custom parameters
  custom_params <- list(
    signaling_weight = 0.8,
    adjust_method = "multiplicative",
    min_effect = 0.05,
    max_effect = 2.0
  )
  
  result2 <- apply_lr_signaling(
    expr_matrix, 
    effect_matrix,
    integration_params = custom_params,
    random_seed = 42
  )
  
  expect_true(is.matrix(result2))
  expect_equal(dim(result2), c(n_cells, n_genes))
  
  # Test with different integration method
  custom_params2 <- list(
    adjust_method = "additive"
  )
  
  result3 <- apply_lr_signaling(
    expr_matrix, 
    effect_matrix,
    integration_params = custom_params2,
    random_seed = 42
  )
  
  expect_true(is.matrix(result3))
})

test_that("generate_lr_interactions integrates full L-R interaction workflow", {
  # Create test data
  n_cells <- 16
  n_genes <- 20
  
  test_df <- data.frame(
    x = rep(1:4, each = 4),
    y = rep(1:4, times = 4)
  )
  
  expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_cells, ncol = n_genes
  )
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(test_df[, c("x", "y")]))
  
  # Run full LR interaction workflow
  result <- generate_lr_interactions(
    test_df, 
    expr_matrix,
    dist_mat,
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("expression", "original_expression", "signaling_effects", "interaction_db") %in% names(result)))
  
  # Output matrix should have same dimensions as input
  expect_equal(dim(result$expression), c(n_cells, n_genes))
  
  # Interaction DB should have expected structure
  expect_true(is.list(result$interaction_db))
  expect_true("interactions" %in% names(result$interaction_db))
  
  # Test with custom parameters
  custom_params <- list(
    n_interactions = 5,
    max_signaling_distance = 30,
    signal_propagation_mode = "threshold",
    adjust_method = "additive"
  )
  
  result2 <- generate_lr_interactions(
    test_df, 
    expr_matrix,
    dist_mat,
    lr_params = custom_params,
    random_seed = 42
  )
  
  expect_true(is.list(result2))
  expect_equal(dim(result2$expression), c(n_cells, n_genes))
})