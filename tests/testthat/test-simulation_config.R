# Test per i moduli di configurazione della simulazione

test_that("initialize_simulation_config valida i parametri", {
  # Crea un file temporaneo per i test
  img_path <- tempfile(fileext = ".png")
  file.create(img_path)
  on.exit(unlink(img_path))
  
  # Test con parametri validi
  config <- initialize_simulation_config(
    image_path = img_path,
    n_cells = 100,
    n_genes = 50,
    random_seed = 123
  )
  
  # Verifica struttura output
  expect_type(config, "list")
  expect_equal(config$image_path, img_path)
  expect_equal(config$n_cells, 100)
  expect_equal(config$n_genes, 50)
  
  # Test errore con immagine non esistente
  expect_error(
    initialize_simulation_config(image_path = "file_non_esistente.png"),
    "È necessario specificare un percorso di immagine valido"
  )
  
  # Test con parametri non validi
  expect_error(
    initialize_simulation_config(image_path = img_path, n_cells = -10),
    "n_cells deve essere positivo"
  )
  
  expect_warning(
    config_warning <- initialize_simulation_config(
      image_path = img_path,
      threshold_value = 1.5
    ),
    "threshold_value dovrebbe essere tra 0 e 1"
  )
  
  expect_equal(config_warning$threshold_value, 0.7)
})

test_that("configure_difficulty_level imposta correttamente i parametri", {
  # Test livello easy
  easy_config <- configure_difficulty_level(
    difficulty_level = "easy",
    marker_params = list(
      marker_genes_per_type = NULL,
      marker_expression_fold = NULL
    ),
    dropout_params = list(dropout_range = NULL)
  )
  
  # Verifica parametri easy
  expect_equal(easy_config$marker_params$marker_genes_per_type, 10)
  expect_equal(easy_config$marker_params$marker_expression_fold, 2.0)
  expect_equal(easy_config$marker_params$marker_overlap_fold, 0.0)
  expect_equal(easy_config$dropout_params$dropout_range, c(0.1, 0.3))
  
  # Test livello medium
  medium_config <- configure_difficulty_level(
    difficulty_level = "medium",
    marker_params = list(marker_genes_per_type = 15),  # Override
    dropout_params = list(dropout_range = NULL)
  )
  
  # Verifica parametri medium con override
  expect_equal(medium_config$marker_params$marker_genes_per_type, 15)
  expect_equal(medium_config$marker_params$marker_expression_fold, 1.2)
  expect_equal(medium_config$dropout_params$dropout_range, c(0.2, 0.5))
  
  # Test livello non valido
  expect_warning(
    invalid_config <- configure_difficulty_level(
      difficulty_level = "invalid",
      marker_params = list(marker_genes_per_type = NULL)
    ),
    "Livello di difficoltà non valido"
  )
  
  # Dovrebbe usare medium come default
  expect_equal(invalid_config$difficulty_level, "medium")
  expect_equal(invalid_config$marker_params$marker_genes_per_type, 7)
})