#' Test Full Simulation Script
#'
#' Questo script esegue una simulazione completa di trascrizione spaziale
#' utilizzando le funzioni del pacchetto senza dipendere dal framework targets.

# Carica tutte le funzioni in ordine
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (file in sort(files)) { source(file) }

# Carica le librerie necessarie
library(imager)
library(ggplot2)
library(dplyr)
library(ClusterR)
library(sp)
library(gstat)

# Aumenta il limite di memoria per future
options(future.globals.maxSize = 100 * 1024^2) # 100 GB

# Imposta il seed per riproducibilità
set.seed(42)

# Definisci parametri di simulazione
params <- list(
  image_path = "images/colon.png",
  output_path = "results/test_simulation.rds",
  output_plot = "results/test_simulation_plot.png",
  n_genes = 50,        # Ridotto per test
  n_cells = 5000,      # Ridotto per test
  k_cell_types = 3,    # Numero di tipi cellulari
  difficulty_level = "medium"
)

# Esegui la simulazione
cat("Avvio simulazione con", params$n_genes, "geni e", params$n_cells, "celle...\n")

risultato <- simulate_spatial_transcriptomics(
  image_path = params$image_path,
  output_path = params$output_path,
  output_plot = params$output_plot,
  n_genes = params$n_genes,
  n_cells = params$n_cells,
  k_cell_types = params$k_cell_types,
  difficulty_level = params$difficulty_level
)

# Stampa statistiche principali
cat("\nSimulazione completata!\n")
cat("Dimensioni matrice di espressione:", dim(risultato$expression)[1], "celle x", 
    dim(risultato$expression)[2], "geni\n")
cat("Distribuzione dei cluster:\n")
print(table(risultato$intensity_cluster))
cat("\nPercorso del report:", params$output_report, "\n")
cat("Percorso del plot:", params$output_plot, "\n")

# Visualizza il plot se in sessione interattiva
if (interactive()) {
  cat("Apertura visualizzazione...\n")
  
  # Ottieni configurazione dalla simulazione
  config <- get_simulation_config(
    image_path = params$image_path,
    n_cells = params$n_cells,
    k_cell_types = params$k_cell_types,
    clustering_method = params$clustering_method
  )
  
  difficulty_config <- setup_difficulty_parameters(
    difficulty_level = params$difficulty_level,
    n_genes = params$n_genes,
    k_cell_types = params$k_cell_types
  )
  
  # Crea visualizzazione
  plot <- create_simulation_plots(
    cell_df = data.frame(
      x = risultato$coordinates[,1],
      y = risultato$coordinates[,2],
      intensity_cluster = as.factor(risultato$intensity_cluster)
    ),
    config = config,
    difficulty_config = difficulty_config,
    expression_results = risultato$expression
  )
  
  print(plot)
}
