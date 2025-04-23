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

test_that("generate_gene_specific_dropout_factors creates correct structure", {
  n_genes <- 50
  n_cells <- 10
  dropout_params <- list(
    gene_dropout_variability = 0.3,
    use_gene_specific_dropout = TRUE,
    gc_content_effect = 0.5,
    length_effect = 0.3,
    sequence_effect = 0.4
  )
  
  res <- generate_gene_specific_dropout_factors(n_genes, n_cells, dropout_params, random_seed = 42)
  
  # Check returned structure
  expect_type(res, "list")
  expect_named(res, c("gene_dropout_factors", "dropout_baseline_shift", "gc_content", "gene_length"))
  
  # Check dimensions
  expect_length(res$gene_dropout_factors, n_genes)
  expect_length(res$dropout_baseline_shift, n_genes)
  expect_length(res$gc_content, n_genes)
  expect_length(res$gene_length, n_genes)
  
  # Check ranges
  expect_true(all(res$gene_dropout_factors >= 0.5 & res$gene_dropout_factors <= 2.0))
  expect_true(all(res$gc_content >= 0 & res$gc_content <= 1))
  expect_true(all(res$gene_length > 0))
})

test_that("apply_dropout with gene_specific_factors applies correct dropout", {
  set.seed(42) # For reproducibility
  expr <- matrix(rep(10, 4*5), nrow = 4) # 4 cells, 5 genes, all expression = 10
  base_dropout <- rep(0.2, 4)
  
  # Create gene-specific factors that make the first gene more prone to dropout (2x)
  # and the second gene less prone (0.5x)
  gene_specific_factors <- list(
    gene_dropout_factors = c(2.0, 0.5, 1.0, 1.0, 1.0),
    dropout_baseline_shift = c(0.1, -0.1, 0, 0, 0)
  )
  
  # Apply dropout with gene-specific factors
  params <- list(
    expression_dependent_dropout = FALSE,
    use_gene_specific_dropout = TRUE,
    gene_effect_weight = 1.0  # Use only gene-specific effect for clearer testing
  )
  
  # Fixed seed to make test deterministic
  set.seed(123)
  
  # Apply dropout
  res <- apply_dropout(expr, base_dropout, gene_specific_factors, dropout_params = params, random_seed = 123)
  
  # Check that the first gene (factor 2.0) has more zeros than others
  zeros_per_gene <- colSums(res == 0)
  expect_true(zeros_per_gene[1] >= zeros_per_gene[3])
  
  # Check that the second gene (factor 0.5) has fewer zeros than baseline
  expect_true(zeros_per_gene[2] <= zeros_per_gene[3])
})

test_that("generate_ambient_rna creates ambient contamination matrix of correct size", {
  # Set up test data
  n_cells <- 5
  n_genes <- 3
  k_cell_types <- 2
  
  cell_df <- data.frame(
    x = 1:n_cells, 
    y = 1:n_cells,
    intensity_cluster = factor(rep(1:2, length.out = n_cells))
  )
  
  dist_mat <- as.matrix(dist(cbind(cell_df$x, cell_df$y)))
  
  # Create a simple mean expression list
  mean_expression_list <- list(
    rep(10, n_genes),
    rep(5, n_genes)
  )
  
  # Test basic configuration
  ambient_params <- list(
    ambient_contamination_rate = 0.1,
    ambient_diffusion_distance = 10,
    tissue_leakage_factor = 0.5,
    background_noise = 0.05
  )
  
  # Generate ambient RNA
  res <- generate_ambient_rna(cell_df, dist_mat, mean_expression_list, ambient_params, random_seed = 42)
  
  # Check dimensions
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(n_cells, n_genes))
  
  # Check values are reasonable (should be small contamination values)
  expect_true(all(res >= 0))
  expect_true(mean(res) > 0 && mean(res) < 1)
})

test_that("apply_ambient_rna_contamination adds contamination correctly", {
  # Create test matrices
  expr_matrix <- matrix(10, nrow = 3, ncol = 2)
  ambient_matrix <- matrix(1, nrow = 3, ncol = 2)
  
  # Apply contamination with scale factor 2.0
  res <- apply_ambient_rna_contamination(
    expr_matrix, 
    ambient_matrix,
    ambient_params = list(ambient_scale_factor = 2.0)
  )
  
  # Check dimensions unchanged
  expect_equal(dim(res), dim(expr_matrix))
  
  # Check values increased by correct amount
  expect_equal(res, matrix(10 + 2*1, nrow = 3, ncol = 2))
})