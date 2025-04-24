test_that("generate_splicing_variants creates valid splice variant assignments", {
  # Test parameters
  n_genes <- 50
  
  # Test with default parameters
  variants <- generate_splicing_variants(
    n_genes,
    splicing_params = list(
      splicing_fraction = 0.6,
      n_splicing_variants = 2
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(variants))
  expect_true(all(c("genes_with_variants", "variant_counts") %in% names(variants)))
  
  # Count of genes with splicing should match parameter
  expect_length(variants$genes_with_variants, round(n_genes * 0.6))
  
  # Each gene should have specified number of variants
  expect_true(all(variants$variant_counts[variants$genes_with_variants] == 2))
  
  # Other genes should have 1 variant (the canonical form)
  expect_true(all(variants$variant_counts[-variants$genes_with_variants] == 1))
  
  # Test with custom parameters
  variants2 <- generate_splicing_variants(
    n_genes,
    splicing_params = list(
      splicing_fraction = 0.3,
      n_splicing_variants = 3
    ),
    random_seed = 42
  )
  
  expect_length(variants2$genes_with_variants, round(n_genes * 0.3))
  expect_true(all(variants2$variant_counts[variants2$genes_with_variants] == 3))
})

test_that("calculate_splicing_propensity creates spatial propensity patterns", {
  # Create test data
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    cluster = factor(rep(1:2, each = 50))
  )
  n_cells <- 100
  
  # Genes with variants
  genes_with_variants <- c(1, 3, 5, 7, 9)
  n_variants <- 2
  n_genes_with_variants <- length(genes_with_variants)
  
  # Test with gradient pattern
  propensity_gradient <- calculate_splicing_propensity(
    test_df,
    genes_with_variants,
    n_variants,
    splicing_params = list(
      splicing_spatial_pattern = "gradient",
      splicing_cluster_specific = FALSE
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.array(propensity_gradient))
  expect_equal(dim(propensity_gradient), c(n_cells, n_genes_with_variants, n_variants))
  
  # Should be probabilities (sum to 1 across variants)
  variant_sums <- apply(propensity_gradient, c(1, 2), sum)
  expect_true(all(abs(variant_sums - 1) < 1e-10))
  
  # Test with cluster-specific pattern
  propensity_cluster <- calculate_splicing_propensity(
    test_df,
    genes_with_variants,
    n_variants,
    splicing_params = list(
      splicing_spatial_pattern = "random",
      splicing_cluster_specific = TRUE
    ),
    random_seed = 42
  )
  
  expect_true(is.array(propensity_cluster))
  expect_equal(dim(propensity_cluster), c(n_cells, n_genes_with_variants, n_variants))
  
  # Should observe cluster-specific patterns
  cluster1_avg <- apply(propensity_cluster[1:50, , ], c(2, 3), mean)
  cluster2_avg <- apply(propensity_cluster[51:100, , ], c(2, 3), mean)
  # At least some genes should show different splicing patterns between clusters
  expect_true(any(abs(cluster1_avg - cluster2_avg) > 0.1))
  
  # Test with focal pattern
  propensity_focal <- calculate_splicing_propensity(
    test_df,
    genes_with_variants,
    n_variants,
    splicing_params = list(
      splicing_spatial_pattern = "focal",
      n_splicing_foci = 3
    ),
    random_seed = 42
  )
  
  expect_true(is.array(propensity_focal))
  expect_equal(dim(propensity_focal), c(n_cells, n_genes_with_variants, n_variants))
})

test_that("generate_splicing_expression creates variant-specific expression", {
  # Test data
  n_cells <- 50
  n_genes <- 10
  base_expr <- matrix(rpois(n_cells * n_genes, lambda = 10), nrow = n_cells, ncol = n_genes)
  
  genes_with_variants <- c(1, 3, 5, 7, 9)
  n_variants <- 2
  variant_counts <- rep(1, n_genes)
  variant_counts[genes_with_variants] <- n_variants
  
  # Generate propensity array
  propensity <- array(0, dim = c(n_cells, length(genes_with_variants), n_variants))
  # Set random probabilities that sum to 1 for each cell-gene combination
  for (i in 1:n_cells) {
    for (j in 1:length(genes_with_variants)) {
      probs <- runif(n_variants)
      propensity[i, j, ] <- probs / sum(probs)
    }
  }
  
  # Generate splicing expression
  result <- generate_splicing_expression(
    base_expr, 
    genes_with_variants, 
    variant_counts,
    propensity,
    splicing_params = list(
      splicing_strength = 1.0,
      variant_expression_ratios = c(1.0, 0.8)
    )
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("expr_matrix", "variant_matrices") %in% names(result)))
  
  # Main matrix should have same dimensions
  expect_equal(dim(result$expr_matrix), dim(base_expr))
  
  # Should have one matrix per gene with variants
  expect_equal(length(result$variant_matrices), length(genes_with_variants))
  
  # Each variant matrix should have dimensions [n_cells, n_variants]
  expect_equal(dim(result$variant_matrices[[1]]), c(n_cells, n_variants))
  
  # Sum of variant expression should approximately equal base expression
  # (allowing for rounding and stochasticity)
  for (i in 1:length(genes_with_variants)) {
    gene_idx <- genes_with_variants[i]
    total_expr <- rowSums(result$variant_matrices[[i]])
    expect_true(all(abs(total_expr - base_expr[, gene_idx]) < 1))
  }
})

test_that("generate_alternative_splicing integrates full splicing workflow", {
  # Test data
  test_df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5),
    cluster = factor(rep(1:5, each = 5))
  )
  n_cells <- 25
  n_genes <- 20
  
  # Base expression matrix
  expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 10),
    nrow = n_cells, ncol = n_genes
  )
  
  # Run full splicing workflow
  result <- generate_alternative_splicing(
    test_df, expr_matrix,
    splicing_params = list(
      use_alternative_splicing = TRUE,
      splicing_fraction = 0.4,
      n_splicing_variants = 2,
      splicing_spatial_pattern = "gradient",
      splicing_strength = 0.8
    ),
    random_seed = 42
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("expr_matrix", "genes_with_variants", "variant_matrices") %in% names(result)))
  
  # Output matrix should have same dimensions as input
  expect_equal(dim(result$expr_matrix), dim(expr_matrix))
  
  # Should have expected number of genes with variants
  expect_length(result$genes_with_variants, round(n_genes * 0.4))
  
  # When disabled, should return original matrix
  result_disabled <- generate_alternative_splicing(
    test_df, expr_matrix,
    splicing_params = list(
      use_alternative_splicing = FALSE
    ),
    random_seed = 42
  )
  
  expect_true(is.list(result_disabled))
  expect_true(identical(result_disabled$expr_matrix, expr_matrix))
  expect_length(result_disabled$genes_with_variants, 0)
  expect_length(result_disabled$variant_matrices, 0)
})