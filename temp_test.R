# Script di test per debug della generazione dei profili di espressione

# Carica tutte le librerie necessarie
library(tictoc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(imager)
library(MASS)
library(cluster)
library(ClusterR)
library(gstat)
library(sp)
library(scales)
library(spatstat)
library(spdep)
library(fields)

# Carica tutti i file dei moduli
tic("Caricamento moduli")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (file in sort(files)) {
  cat("Caricamento:", basename(file), "\n")
  source(file)
}
toc()

# Crea un set di dati di test
set.seed(123)
test_df <- data.frame(
  x = rep(1:20, each = 20),
  y = rep(1:20, times = 20),
  intensity_cluster = factor(sample(1:2, 400, replace=TRUE))
)

# Prova a generare i profili di espressione
tic("Test generazione profili")
tryCatch({
  result <- generate_expression_profiles(
    cell_df = test_df,
    n_genes = 20,
    k_cell_types = 2,
    use_spatial_correlation = TRUE,
    marker_params = list(
      marker_genes_per_type = 5,
      marker_expression_fold = 1.5
    ),
    dropout_params = list(
      dropout_range = c(0.2, 0.5),
      expression_dependent_dropout = TRUE
    ),
    cell_specific_params = list(
      use_gene_modules = TRUE,
      module_hierarchical = FALSE,
      module_overlap = 0.1,
      cell_specific_noise_sd = 0.2
    ),
    hybrid_params = list(
      use_hybrid_cells = TRUE,
      hybrid_intensity_range = c(0.2, 0.5)
    ),
    random_seed = 123
  )
  
  # Verifica se il risultato è valido
  cat("\nVerifica risultato:\n")
  cat("Dimensioni della matrice di espressione:", dim(result$expression), "\n")
  cat("Conteggio medio:", mean(rowSums(result$expression)), "\n")
  cat("Percentuale dropout:", 100 * sum(result$expression == 0) / prod(dim(result$expression)), "%\n")
  
}, error = function(e) {
  cat("\nErrore durante la generazione dei profili:\n")
  print(e)
  
  # Verifica quali funzioni dipendenti potrebbero causare il problema
  cat("\nVerifica funzioni principali:\n")
  tryCatch({
    set.seed(123)
    params <- initialize_expression_params(
      marker_params = list(marker_genes_per_type = 5, marker_expression_fold = 1.5),
      spatial_params = list(),
      dropout_params = list(dropout_range = c(0.2, 0.5)),
      library_size_params = list(),
      cell_specific_params = list(use_gene_modules = TRUE),
      hybrid_params = list(use_hybrid_cells = TRUE),
      random_seed = 123
    )
    cat("initialize_expression_params: OK\n")
  }, error = function(e) {
    cat("initialize_expression_params: ERRORE -", e$message, "\n")
  })
  
  tryCatch({
    set.seed(123)
    mean_expression_list <- generate_baseline_expression(20, 2, list(marker_genes_per_type = 5), 123)
    cat("generate_baseline_expression: OK\n")
  }, error = function(e) {
    cat("generate_baseline_expression: ERRORE -", e$message, "\n")
  })
  
  tryCatch({
    set.seed(123)
    spatial_distances <- calculate_spatial_distances(test_df, chunk_size = NULL, 123)
    cat("calculate_spatial_distances: OK\n")
  }, error = function(e) {
    cat("calculate_spatial_distances: ERRORE -", e$message, "\n")
  })
  
  tryCatch({
    set.seed(123)
    if (exists("spatial_distances")) {
      dispersion_param <- calculate_dispersion_params(
        test_df, spatial_distances$mean_dist, 
        list(dropout_range = c(0.2, 0.5)), 2, 123
      )
      cat("calculate_dispersion_params: OK\n")
    } else {
      cat("calculate_dispersion_params: SKIP (dipende da spatial_distances)\n")
    }
  }, error = function(e) {
    cat("calculate_dispersion_params: ERRORE -", e$message, "\n")
  })
  
  cat("\nGli errori possono verificarsi in altre funzioni.\n")
})
toc()