test_that("generate_expression_profiles crea una matrice di espressione valida", {
  # Crea un dataset di test semplice
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Chiamata alla funzione con parametri minimi
  set.seed(123)
  result <- generate_expression_profiles(
    cell_df = test_df,
    n_genes = 10,
    k_cell_types = 2,
    use_spatial_correlation = FALSE,  # Disabilita per semplificare il test
    hybrid_params = list(
      use_hybrid_cells = FALSE
    ),
    cell_specific_params = list(
      use_gene_modules = FALSE,
      cell_specific_noise_sd = 0.1
    ),
    random_seed = 123
  )
  
  # Verifiche sulla struttura del risultato
  expect_type(result, "list")
  expect_true("expression" %in% names(result))
  expect_true("library_size" %in% names(result))
  expect_true("dispersion_param" %in% names(result))
  
  # Verifica dimensioni della matrice di espressione
  expect_equal(dim(result$expression), c(100, 10))  # 100 cellule, 10 geni
  
  # Verifica che i valori siano numeri interi non negativi (conteggi)
  expect_true(all(result$expression >= 0))
  expect_true(all(result$expression == round(result$expression)))
  
  # Verifica che i parametri library_size e dispersion siano corretti
  expect_equal(length(result$library_size), 100)
  expect_equal(length(result$dispersion_param), 100)
  expect_true(all(result$library_size > 0))
  expect_true(all(result$dispersion_param > 0))
})

test_that("generate_expression_profiles rispetta il seed per la riproducibilità", {
  # Crea un dataset di test semplice
  test_df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5),
    intensity_cluster = factor(c(rep(1, 12), rep(2, 13)))
  )
  
  # Parametri semplificati
  common_params <- list(
    cell_df = test_df,
    n_genes = 5,
    k_cell_types = 2,
    use_spatial_correlation = FALSE,
    hybrid_params = list(use_hybrid_cells = FALSE),
    cell_specific_params = list(use_gene_modules = FALSE, cell_specific_noise_sd = 0.1)
  )
  
  # Prima esecuzione
  set.seed(42)
  result1 <- do.call(generate_expression_profiles, c(common_params, list(random_seed = 42)))
  
  # Seconda esecuzione con lo stesso seed
  set.seed(99)  # Impostiamo un seed diverso per R, ma la funzione dovrebbe usare il suo
  result2 <- do.call(generate_expression_profiles, c(common_params, list(random_seed = 42)))
  
  # Terza esecuzione con seed diverso
  set.seed(42)  # Ripristiniamo il seed di R
  result3 <- do.call(generate_expression_profiles, c(common_params, list(random_seed = 99)))
  
  # Verifica la riproducibilità
  expect_identical(result1$expression, result2$expression)
  expect_false(identical(result1$expression, result3$expression))
})
