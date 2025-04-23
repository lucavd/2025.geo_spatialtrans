# Test per le funzioni wrapper

test_that("generate_expression_profiles_wrapper coordina correttamente i moduli", {
  # Creo un dataset di test semplice
  test_df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5),
    intensity_cluster = factor(rep(1:2, each = 12, length.out = 25))
  )
  
  # Chiama la funzione wrapper con parametri semplificati
  set.seed(123)
  result <- generate_expression_profiles(
    cell_df = test_df,
    n_genes = 10,
    k_cell_types = 2,
    use_spatial_correlation = FALSE,  # Disabilito per semplicità
    hybrid_params = list(
      use_hybrid_cells = FALSE
    ),
    cell_specific_params = list(
      use_gene_modules = FALSE,
      cell_specific_noise_sd = 0.1
    ),
    random_seed = 123
  )
  
  # Verifica la struttura del risultato
  expect_type(result, "list")
  expect_true(all(c("expression", "library_size", "dispersion_param", "mean_expression_list") %in% names(result)))
  
  # Verifica la dimensione della matrice di espressione
  expect_equal(dim(result$expression), c(25, 10))  # 25 celle, 10 geni
  
  # Verifica che i conteggi siano non negativi e interi
  expect_true(all(result$expression >= 0))
  expect_true(all(result$expression == as.integer(result$expression)))
  
  # Verifica dimensioni libreria
  expect_equal(length(result$library_size), 25)
  expect_true(all(result$library_size > 0))
  
  # Verifica parametri di dispersione
  expect_equal(length(result$dispersion_param), 25)
  expect_true(all(result$dispersion_param > 0))
})

test_that("simulate_spatial_transcriptomics_wrapper integra correttamente i componenti", {
  # Usa l'immagine synthetic_tissue1.png dall'archivio
  image_path <- file.path("..", "..", "images", "synthetic_tissue1.png")
  
  # Verifica che l'immagine esista
  if (!file.exists(image_path)) {
    skip("Immagine di test non trovata")
  }
  
  # Crea un percorso temporaneo per l'output
  output_path <- tempfile(fileext = ".rds")
  
  # Esegui la simulazione con parametri minimi per velocizzare il test
  result <- simulate_spatial_transcriptomics(
    image_path = image_path,
    output_path = output_path,
    n_cells = 100,        # Numero ridotto di celle
    n_genes = 10,          # Numero ridotto di geni
    k_cell_types = 2,      # Solo 2 tipi cellulari
    grid_mode = TRUE,
    grid_resolution = 10,  # Risoluzione più bassa per velocità
    random_seed = 123,
    # Evitare calcoli pesanti per il test
    use_spatial_correlation = FALSE,
    hybrid_params = list(use_hybrid_cells = FALSE),
    cell_specific_params = list(use_gene_modules = FALSE)
  )
  
  # Verifiche sul risultato
  expect_true(file.exists(output_path))
  expect_type(result, "list")
  expect_true(all(c("coordinates", "intensity_cluster", "expression", "parameters") %in% names(result)))
  
  # Verifica dimensioni matrice espressione
  # La dimensione effettiva dipende dall'immagine e dal grid_mode
  # Controlliamo che ci siano delle celle e che il numero corrisponda a coordinates
  expect_equal(nrow(result$expression), nrow(result$coordinates))
  expect_equal(ncol(result$expression), 10)   # 10 geni
  
  # Verifica che i conteggi siano non negativi e interi
  expect_true(all(result$expression >= 0))
  expect_true(all(result$expression == as.integer(result$expression)))
  
  # Pulizia del file temporaneo
  if (file.exists(output_path)) {
    file.remove(output_path)
  }
})