# Script per verificare i nuovi moduli spacchettati per geo.spatialtrans
# Questo script dimostra come utilizzare i nuovi moduli e come sono interconnessi

# Carica i file dei moduli in ordine
library(tictoc)
library(dplyr)
library(ggplot2)

tic("Caricamento moduli")
# Carica tutti i file dei moduli
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (file in sort(files)) {
  cat("Caricamento:", basename(file), "\n")
  source(file)
}
toc()

# Configurazione del reporter di test personalizzato per report più dettagliati
detailed_reporter <- testthat::ListReporter$new()

# Verifica i test
tic("Esecuzione tests")
test_results <- testthat::test_dir("tests/testthat/", reporter = detailed_reporter)
toc()

# Report dettagliato
cat("\n=== REPORT DETTAGLIATO DEI TEST ===\n")
cat("Numero totale di test eseguiti:", sum(sapply(detailed_reporter$get_results(), function(x) length(x$results))), "\n")

# Conta test passati, falliti, con avvisi e saltati
test_data <- sapply(detailed_reporter$get_results(), function(x) {
  results <- x$results
  passed <- sum(sapply(results, function(r) r$passed && !inherits(r, "expectation_warning")))
  warnings <- sum(sapply(results, function(r) inherits(r, "expectation_warning")))
  errors <- sum(sapply(results, function(r) inherits(r, "expectation_error")))
  failures <- sum(sapply(results, function(r) inherits(r, "expectation_failure")))
  skipped <- sum(sapply(results, function(r) inherits(r, "expectation_skip")))
  c(passed = passed, warnings = warnings, errors = errors, failures = failures, skipped = skipped, total = length(results))
})

# Crea un riassunto per file
test_summary <- data.frame(
  file = names(detailed_reporter$get_results()),
  passed = test_data["passed", ],
  warnings = test_data["warnings", ],
  errors = test_data["errors", ],
  failures = test_data["failures", ],
  skipped = test_data["skipped", ],
  total = test_data["total", ]
)

# Stampa il riassunto
for (i in 1:nrow(test_summary)) {
  cat(sprintf("%s: %d/%d test passati", 
              test_summary$file[i], 
              test_summary$passed[i],
              test_summary$total[i]))
  
  if (test_summary$warnings[i] > 0) 
    cat(sprintf(", %d avvisi", test_summary$warnings[i]))
  if (test_summary$errors[i] > 0) 
    cat(sprintf(", %d errori", test_summary$errors[i]))
  if (test_summary$failures[i] > 0) 
    cat(sprintf(", %d fallimenti", test_summary$failures[i]))
  if (test_summary$skipped[i] > 0) 
    cat(sprintf(", %d saltati", test_summary$skipped[i]))
  
  cat("\n")
}

# Riepilogo
total_passed <- sum(test_summary$passed)
total_tests <- sum(test_summary$total)
coverage <- round(total_passed / total_tests * 100, 1)

cat("\nRIEPILOGO:\n")
cat(sprintf("Test passati: %d/%d (%.1f%%)\n", total_passed, total_tests, coverage))
cat(sprintf("Avvisi: %d\n", sum(test_summary$warnings)))
cat(sprintf("Errori: %d\n", sum(test_summary$errors)))
cat(sprintf("Fallimenti: %d\n", sum(test_summary$failures)))
cat(sprintf("Test saltati: %d\n", sum(test_summary$skipped)))

# Esempio 1: Generazione di un profilo di espressione semplice
cat("\n\n=== ESEMPIO 1: Profilo di espressione semplice ===\n")
tic("Esempio 1")

# Crea un dataframe di celle semplice
simple_cells <- data.frame(
  x = rep(1:10, each = 10),
  y = rep(1:10, times = 10),
  intensity_cluster = factor(rep(c(1, 2, 3), each = 33, length.out = 100))
)

# Genera profili di espressione
simple_expression <- generate_expression_profiles(
  cell_df = simple_cells,
  n_genes = 15,
  k_cell_types = 3,
  use_spatial_correlation = TRUE,
  correlation_method = "grf",
  random_seed = 123
)

# Visualizza le dimensioni della matrice
cat("\nDimensioni matrice espressione:", dim(simple_expression$expression), "\n")
toc()

# Esempio 2: Utilizzo dei singoli moduli
cat("\n\n=== ESEMPIO 2: Utilizzo moduli separati ===\n")
tic("Esempio 2")

# 1. Inizializza parametri
params <- initialize_expression_params(
  marker_params = list(
    marker_genes_per_type = 5,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  ),
  dropout_params = list(
    dropout_range = c(0.1, 0.3),
    dispersion_range = c(2.0, 1.0)
  ),
  random_seed = 456
)

# 2. Genera profili baseline
baseline_expr <- generate_baseline_expression(
  n_genes = 15,
  k_cell_types = 3,
  marker_params = params$marker_params,
  random_seed = 456
)

# 3. Calcola distanze
spatial_dist <- calculate_spatial_distances(
  simple_cells,
  random_seed = 456
)

# 4. Calcola dispersione
dispersion <- calculate_dispersion_params(
  simple_cells,
  spatial_dist$mean_dist,
  params$dropout_params,
  k_cell_types = 3,
  random_seed = 456
)

# Visualizza risultati
cat("\nMedia dispersione:", mean(dispersion), "\n")
toc()

# Esempio 3: Configurazione delle difficoltà
cat("\n\n=== ESEMPIO 3: Configurazione difficoltà ===\n")
tic("Esempio 3")

# Ottieni parametri per diversi livelli di difficoltà
easy_config <- configure_difficulty_level("easy")
medium_config <- configure_difficulty_level("medium")
hard_config <- configure_difficulty_level("hard")

# Stampa differenze nei marker genes
cat("\nMarker genes per tipo (easy):", easy_config$marker_params$marker_genes_per_type, "\n")
cat("Marker genes per tipo (medium):", medium_config$marker_params$marker_genes_per_type, "\n")
cat("Marker genes per tipo (hard):", hard_config$marker_params$marker_genes_per_type, "\n")

# Stampa differenze nel fold change
cat("\nFold change (easy):", easy_config$marker_params$marker_expression_fold, "\n")
cat("Fold change (medium):", medium_config$marker_params$marker_expression_fold, "\n")
cat("Fold change (hard):", hard_config$marker_params$marker_expression_fold, "\n")

# Stampa differenze nel dropout
cat("\nDropout range (easy):", paste(easy_config$dropout_params$dropout_range, collapse="-"), "\n")
cat("Dropout range (medium):", paste(medium_config$dropout_params$dropout_range, collapse="-"), "\n")
cat("Dropout range (hard):", paste(hard_config$dropout_params$dropout_range, collapse="-"), "\n")
toc()

# Esempio 4: Simulazione completa con immagine di test
cat("\n\n=== ESEMPIO 4: Simulazione completa con immagine ===\n")
tic("Esempio 4")

# Percorso immagine di test
image_path <- "images/synthetic_tissue1.png"
if (!file.exists(image_path)) {
  cat("ATTENZIONE: Immagine di test non trovata:", image_path, "\n")
} else {
  # Esegui una simulazione molto semplificata per velocità
  result <- simulate_spatial_transcriptomics(
    image_path = image_path,
    output_path = tempfile(fileext = ".rds"),
    n_cells = 50,
    n_genes = 8,
    k_cell_types = 2,
    grid_mode = TRUE,
    grid_resolution = 15,
    random_seed = 789,
    use_spatial_correlation = FALSE,
    hybrid_params = list(use_hybrid_cells = FALSE),
    cell_specific_params = list(use_gene_modules = FALSE)
  )
  
  # Mostra informazioni
  cat("\nDimensioni matrice espressione:", dim(result$expression), "\n")
  cat("Numero di tipi cellulari:", length(unique(result$intensity_cluster)), "\n")
}
toc()

cat("\n\nVerifica completata!\n")