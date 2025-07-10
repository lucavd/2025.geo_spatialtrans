#!/usr/bin/env Rscript
# Script di debug per verificare parametri library size

cat("=== DEBUG PARAMETRI LIBRARY SIZE ===\n")

# 1. Caricamento delle funzioni e librerie
cat("Caricamento librerie...\n")
suppressPackageStartupMessages({
  library(png)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(ClusterR)
  library(sp)
  library(gstat)
})

cat("Caricamento funzioni...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# 2. Configurazione parametri ottimizzati
cat("\n=== CONFIGURAZIONE PARAMETRI ===\n")
diff_cfg <- configure_difficulty_level("medium")

# Parametri ottimizzati
optimized_params <- list(
  mean_library_size = 7000,
  library_size_cv = 0.28,
  spatial_effect_on_library = 0.1,
  cell_type_effect = TRUE
)

cat("Parametri ottimizzati:\n")
str(optimized_params)

# Applicazione ai due livelli
diff_cfg$cell_specific_params$library_size_params <- optimized_params
diff_cfg$library_size_params <- optimized_params

cat("\nParametri applicati a diff_cfg:\n")
cat("- cell_specific_params$library_size_params$mean_library_size:", 
    diff_cfg$cell_specific_params$library_size_params$mean_library_size, "\n")
cat("- library_size_params$mean_library_size:", 
    diff_cfg$library_size_params$mean_library_size, "\n")

# 3. Test su dataset ridotto
cat("\n=== TEST SU DATASET RIDOTTO ===\n")
cfg <- initialize_simulation_config(
  image_path = "images/granuloma.png",
  output_path = "results/debug_test.rds",
  output_plot = "results/debug_test.png",
  n_genes = 1000,  # Ridotto per test rapido
  k_cell_types = 5,
  threshold_value = 0.7,
  random_seed = 42,
  pixel_size_um = 10,
  grid_mode = TRUE,
  grid_resolution = 50,  # Griglia più rada
  grid_spacing = 0,
  use_fixed_grid = TRUE,
  fixed_grid_width_mm = 2.0,  # Area ridotta
  fixed_grid_height_mm = 2.0
)

# 4. Simulazione con debug usando parametri diretti
cat("\nEsecuzione simulazione con debug...\n")
simulation_results <- simulate_spatial_transcriptomics(
  image_path = "images/granuloma.png",
  output_path = "results/debug_test.rds",
  output_plot = "results/debug_test.png",
  n_genes = 1000,  # Ridotto per test rapido
  k_cell_types = 5,
  threshold_value = 0.7,
  random_seed = 42,
  pixel_size_um = 10,
  grid_mode = TRUE,
  grid_resolution = 50,  # Griglia più rada
  grid_spacing = 0,
  use_fixed_grid = TRUE,
  fixed_grid_width_mm = 2.0,  # Area ridotta
  fixed_grid_height_mm = 2.0,
  difficulty_level = "medium",
  library_size_params = optimized_params  # PARAMETRI OTTIMIZZATI DIRETTI
)

# 5. Verifica risultati
cat("\n=== VERIFICA RISULTATI ===\n")
if (!is.null(simulation_results$expression)) {
  total_umi <- Matrix::colSums(simulation_results$expression)
  
  cat("Statistiche UMI generate:\n")
  cat("- Media UMI:", round(mean(total_umi)), "\n")
  cat("- Mediana UMI:", round(median(total_umi)), "\n")
  cat("- Range UMI:", paste(range(total_umi), collapse=" - "), "\n")
  cat("- SD UMI:", round(sd(total_umi)), "\n")
  cat("- CV UMI:", round(sd(total_umi)/mean(total_umi), 3), "\n")
  
  # Verifica target
  target_mean <- optimized_params$mean_library_size
  actual_mean <- mean(total_umi)
  diff_percent <- abs(actual_mean - target_mean) / target_mean * 100
  
  cat("\nVERIFICA TARGET:\n")
  cat("- Target mean:", target_mean, "\n")
  cat("- Actual mean:", round(actual_mean), "\n")
  cat("- Differenza:", round(diff_percent, 1), "%\n")
  
  if (diff_percent < 20) {
    cat("✓ PARAMETRI APPLICATI CORRETTAMENTE\n")
  } else {
    cat("✗ PARAMETRI NON APPLICATI - ERRORE\n")
  }
  
  # Salva risultati per ispezione
  saveRDS(simulation_results, "results/debug_test.rds")
  cat("\nRisultati salvati in: results/debug_test.rds\n")
  
} else {
  cat("✗ ERRORE: Simulazione fallita\n")
}

cat("\n=== DEBUG COMPLETATO ===\n")