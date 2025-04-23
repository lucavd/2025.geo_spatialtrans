# Test per moduli di correlazione spaziale

test_that("generate_spatial_correlation crea pattern spaziali validi", {
  # Crea un dataset semplice di test
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Test con metodo GRF
  grf_noise <- generate_spatial_correlation(
    test_df,
    spatial_params = list(
      spatial_noise_intensity = 1.0,
      spatial_range = 20
    ),
    correlation_method = "grf",
    random_seed = 123
  )
  
  # Verifica risultato GRF
  expect_equal(length(grf_noise), 100)
  expect_true(all(!is.na(grf_noise)))
  
  # Test con metodo CAR
  car_noise <- generate_spatial_correlation(
    test_df,
    spatial_params = list(
      spatial_noise_intensity = 1.0,
      spatial_range = 20
    ),
    correlation_method = "car",
    random_seed = 123
  )
  
  # Verifica risultato CAR
  expect_equal(length(car_noise), 100)
  expect_true(all(!is.na(car_noise)))
  
  # Test con metodo non supportato
  null_noise <- generate_spatial_correlation(
    test_df,
    spatial_params = list(
      spatial_noise_intensity = 1.0,
      spatial_range = 20
    ),
    correlation_method = "invalid",
    random_seed = 123
  )
  
  # Dovrebbe restituire NULL
  expect_null(null_noise)
})

test_that("generate_multiscale_correlation crea pattern multi-scala", {
  # Crea un dataset semplice di test
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Test con parametri di base
  ms_result <- generate_multiscale_correlation(
    test_df,
    spatial_params = list(
      use_multiscale = TRUE,
      global_contribution = 0.6,
      macro_range = 50,
      micro_range = 10,
      n_hierarchical_levels = 2
    ),
    random_seed = 42
  )
  
  # Verifica risultato
  expect_equal(length(ms_result), 100)
  expect_true(all(!is.na(ms_result)))
  
  # Test disabilitato
  null_result <- generate_multiscale_correlation(
    test_df,
    spatial_params = list(use_multiscale = FALSE),
    random_seed = 42
  )
  
  # Dovrebbe restituire NULL
  expect_null(null_result)
})

test_that("generate_nonstationary_correlation produce risultati validi", {
  # Crea un dataset semplice di test
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Test con parametri di base
  ns_result <- generate_nonstationary_correlation(
    test_df,
    spatial_params = list(
      use_nonstationary = TRUE,
      nonstationary_type = "gradient",
      range_min = 10,
      range_max = 50
    ),
    random_seed = 42
  )
  
  # Verifica risultato
  expect_true(is.list(ns_result))
  expect_true("noise" %in% names(ns_result))
  expect_equal(length(ns_result$noise), 100)
  expect_true(all(!is.na(ns_result$noise)))
  
  # Verifica che restituisca anche le mappe di parametri
  expect_true("range_map" %in% names(ns_result))
  expect_true("intensity_map" %in% names(ns_result))
  
  # Test disabilitato
  null_result <- generate_nonstationary_correlation(
    test_df,
    spatial_params = list(use_nonstationary = FALSE),
    random_seed = 42
  )
  
  # Dovrebbe restituire NULL
  expect_null(null_result)
})

test_that("generate_composite_spatial_patterns combina correttamente i pattern", {
  # Crea un dataset semplice di test
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Test con parametri di base - usiamo solo multiscale per evitare problemi con dataset piccoli
  cp_result <- generate_composite_spatial_patterns(
    test_df,
    spatial_params = list(
      use_composite = TRUE,
      composite_methods = c("multiscale"),
      method_weights = c(1.0),
      multiscale_params = list(
        use_multiscale = TRUE,
        global_contribution = 0.6,
        macro_range = 50,
        micro_range = 10
      )
    ),
    random_seed = 42
  )
  
  # Verifica risultato
  expect_equal(length(cp_result), 100)
  expect_true(all(!is.na(cp_result)))
  
  # Test disabilitato
  null_result <- generate_composite_spatial_patterns(
    test_df,
    spatial_params = list(use_composite = FALSE),
    random_seed = 42
  )
  
  # Dovrebbe restituire NULL
  expect_null(null_result)
})

test_that("generate_hybrid_cells identifica cellule ai confini", {
  # Crea un dataset test con confini chiari
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Calcola la matrice di distanza
  dist_mat <- as.matrix(dist(test_df[, c("x", "y")]))
  
  # Genera cellule ibride
  hybrid_matrix <- generate_hybrid_cells(
    test_df,
    dist_mat,
    k_cell_types = 2,
    hybrid_params = list(
      use_hybrid_cells = TRUE,
      max_hybrid_pairs = 20,
      hybrid_intensity_range = c(0.2, 0.5)
    ),
    random_seed = 123
  )
  
  # Verifica struttura risultato
  expect_equal(dim(hybrid_matrix), c(100, 2))
  
  # Verifico che ci siano cellule ibride
  hybrid_cells <- which(rowSums(hybrid_matrix) > 0)
  expect_true(length(hybrid_cells) > 0)
  expect_true(length(hybrid_cells) <= 40)  # max 20 coppie = max 40 cellule
  
  # Le intensità dovrebbero essere nel range specificato
  nonzero_intensities <- hybrid_matrix[hybrid_matrix > 0]
  expect_true(all(nonzero_intensities >= 0.2 & nonzero_intensities <= 0.5))
  
  # Test con cellule ibride disabilitate
  hybrid_matrix_disabled <- generate_hybrid_cells(
    test_df,
    dist_mat,
    k_cell_types = 2,
    hybrid_params = list(
      use_hybrid_cells = FALSE,
      max_hybrid_pairs = 20,
      hybrid_intensity_range = c(0.2, 0.5)
    ),
    random_seed = 123
  )
  
  # Dovrebbe essere una matrice di zeri
  expect_equal(sum(hybrid_matrix_disabled), 0)
})