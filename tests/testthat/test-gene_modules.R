test_that("generate_gene_modules produces correct module sizes and no overlap", {
  n_genes <- 20
  n_cells <- 5
  params <- list(
    use_gene_modules = TRUE,
    n_gene_modules = 4,
    module_hierarchical = FALSE,
    module_overlap = 0,
    module_size_distribution = "uniform",
    n_latent_factors = 0,
    module_network_density = 0,
    module_correlation = 0.7
  )
  res <- generate_gene_modules(n_genes, n_cells, cell_specific_params = params, random_seed = 42)
  # gene_modules list
  expect_type(res, "list")
  expect_true("gene_modules" %in% names(res) || is.list(res))
  gm <- res$gene_modules
  # Should have 4 modules of equal size (20/4)
  expect_length(gm, 4)
  sizes <- vapply(gm, length, integer(1))
  expect_true(all(sizes == 5))
  # Modules should be disjoint and cover all genes
  all_genes <- sort(unlist(gm))
  expect_equal(length(unique(all_genes)), n_genes)
  expect_equal(sort(unique(all_genes)), 1:n_genes)
  # module_noise and module_network dimensions
  expect_true(is.matrix(res$module_noise))
  expect_equal(dim(res$module_noise), c(n_cells, n_genes))
  expect_true(is.matrix(res$module_network))
  expect_equal(dim(res$module_network), c(4, 4))
  # latent_factors should be NULL for n_latent_factors = 0
  expect_null(res$latent_factors)
})