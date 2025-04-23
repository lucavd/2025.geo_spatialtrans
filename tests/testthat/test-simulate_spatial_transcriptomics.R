test_that("simulate_spatial_transcriptomics verifica l'input immagine", {
  # Test con image_path NULL
  expect_error(
    simulate_spatial_transcriptomics(image_path = NULL),
    "È necessario specificare un percorso di immagine valido"
  )
})

# Questo test richiede un'immagine di test - in una versione reale del pacchetto,
# potremmo includere un'immagine di esempio nelle fixtures di test

test_that("simulate_spatial_transcriptomics produce risultati con la struttura corretta", {
  # Setup
  image_path <- file.path("..", "..", "images", "synthetic_tissue1.png")
  tmp_dir <- tempdir()
  tmp_output <- file.path(tmp_dir, "test_output.rds")

  # Chiamata alla funzione con parametri semplificati per velocizzare il test
  result <- simulate_spatial_transcriptomics(
    image_path = image_path,
    output_path = tmp_output,
    n_genes = 10,
    k_cell_types = 2,
    grid_mode = TRUE,
    grid_resolution = 5,
    n_cells = 100,
    random_seed = 123,
    use_spatial_correlation = FALSE,
    hybrid_params = list(use_hybrid_cells = FALSE),
    cell_specific_params = list(use_gene_modules = FALSE)
  )

  # Verifica della struttura dell'output
  expect_type(result, "list")
  expect_true(all(c("coordinates", "intensity_cluster", "expression", "parameters") %in% names(result)))

  # Verifica che il file sia stato creato
  expect_true(file.exists(tmp_output))

  # Cleanup
  if (file.exists(tmp_output)) file.remove(tmp_output)
})

# In un ambiente di sviluppo reale, dobbiamo decidere come gestire
# la dipendenza dall'immagine per i test (includere nelle fixtures o
# simulare i dati di immagine). Per ora, il test commentato sopra
# serve come esempio di come potrebbe essere strutturato.