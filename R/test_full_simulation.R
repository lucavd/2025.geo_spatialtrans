#' Test Full Simulation Script
#'
#' Questo script esegue una simulazione completa di trascrizione spaziale
#' utilizzando le funzioni del pacchetto senza dipendere dal framework targets.

# Carica tutte le funzioni in ordine
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (file in sort(files)) { source(file) }

# Carica le librerie necessarie
library(png)
library(ggplot2)
library(dplyr)
library(ClusterR)
library(sp)
library(gstat)
library(Matrix)
# Profiling and memory monitoring:
# • Use profmem::profmem() or lineprof to pinpoint the biggest temporary allocations.
# • Use pryr::mem_used() at key steps to watch your code's memory footprint.

# Aumenta il limite di memoria per future
options(future.globals.maxSize = 100 * 1024^2) # 100 GB

# Imposta il seed per riproducibilità
set.seed(42)

# Definisci parametri di simulazione
params <- list(
  image_path = "images/generated2.png",
  output_path = "results/test_simulation_gen2.rds",
  output_plot = "results/test_simulation_gen2_plot.png",
  n_genes = 50,        # Ridotto per test
  n_cells = 5000,      # Ridotto per test
  k_cell_types = 3,    # Numero di tipi cellulari
  difficulty_level = "medium"
)

# Chunked simulation: process spots in blocks and build sparse expression matrix
config <- initialize_simulation_config(
  image_path   = params$image_path,
  output_path  = params$output_path,
  output_plot  = params$output_plot,
  n_cells      = params$n_cells,
  n_genes      = params$n_genes,
  k_cell_types = params$k_cell_types
)
difficulty_config <- configure_difficulty_level(
  difficulty_level = params$difficulty_level
)
# Setup library size parameters (default values)
library_size_params <- list(
  mean_library_size       = 10000,
  library_size_cv         = 0.3,
  spatial_effect_on_library = 0.5,
  cell_type_effect        = TRUE
)
# Setup hybrid cell parameters (default values)
hybrid_params <- list(
  use_hybrid_cells       = TRUE,
  max_hybrid_pairs       = 1000,
  hybrid_intensity_range = c(0.2, 0.5)
)
# Prepare image and cluster labels
image_data <- prepare_image(config$image_path, threshold_value = config$threshold_value)
clustered_data <- cluster_image(
  img_df_thresh = image_data$img_df_thresh,
  k_cell_types  = config$k_cell_types,
  random_seed   = config$random_seed
)
# Create sampling grid of spots/cells
cell_df <- create_sampling_grid(
  img_df_thresh = clustered_data,
  img_array      = image_data$img_array,
  img_width      = image_data$width,
  img_height     = image_data$height,
  grid_mode      = config$grid_mode,
  n_cells        = config$n_cells,
  k_cell_types   = config$k_cell_types,
  grid_resolution = config$grid_resolution,
  grid_spacing   = config$grid_spacing,
  use_fixed_grid = config$use_fixed_grid,
  fixed_grid_width_mm  = config$fixed_grid_width_mm,
  fixed_grid_height_mm = config$fixed_grid_height_mm,
  pixel_size_um  = config$pixel_size_um,
  threshold_value = config$threshold_value,
  random_seed    = config$random_seed
)
# Chunk settings: number of spots per block
chunk_size <- 1000
spot_indices <- seq_len(nrow(cell_df))
chunks <- split(spot_indices, ceiling(seq_along(spot_indices) / chunk_size))
# Allocate list to hold sparse blocks
sparse_blocks <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  idx <- chunks[[i]]
  block_df <- cell_df[idx, ]
  # Generate expression for this block (disable full spatial correlation)
  expr_res <- generate_expression_profiles(
    cell_df               = block_df,
    n_genes               = config$n_genes,
    k_cell_types          = config$k_cell_types,
    marker_params         = difficulty_config$marker_params,
    spatial_params        = difficulty_config$spatial_params,
    dropout_params        = difficulty_config$dropout_params,
    library_size_params   = library_size_params,
    cell_specific_params  = difficulty_config$cell_specific_params,
    hybrid_params         = hybrid_params,
    use_spatial_correlation = FALSE,
    correlation_method    = "grf",
    random_seed           = config$random_seed
  )
  # Convert dense block to sparse
  sparse_blocks[[i]] <- Matrix(expr_res$expression, sparse = TRUE)
  # Clean up intermediate objects
  rm(expr_res)
  gc()
  message(sprintf("Finished chunk %d/%d", i, length(chunks)))
}
# Combine sparse blocks into full expression matrix by stacking rows (cells)
full_expression <- do.call(rbind, sparse_blocks)
rm(sparse_blocks)
gc()
# Assemble result list
risultato <- list(
  expression        = full_expression,
  coordinates       = as.matrix(cell_df[, c("x", "y")]),
  intensity_cluster = cell_df$intensity_cluster
)
# Save simulated results
saveRDS(risultato, file = config$output_path)

# Stampa statistiche principali
cat("\nSimulazione completata!\n")
cat("Dimensioni matrice di espressione:", dim(risultato$expression)[1], "celle x",
    dim(risultato$expression)[2], "geni\n")
cat("Distribuzione dei cluster:\n")
print(table(risultato$intensity_cluster))
cat("\nPercorso del report: results/test_simulation_gen2_report.txt\n")
cat("Percorso del plot:", params$output_plot, "\n")

# Visualizza il plot se in sessione interattiva
if (interactive()) {
  cat("Apertura visualizzazione...\n")

  # Prepara dataframe per il plot
  cell_plot_df <- data.frame(
    x = risultato$coordinates[,1],
    y = risultato$coordinates[,2],
    intensity_cluster = as.factor(risultato$intensity_cluster)
  )

  # Crea e salva il plot
  generate_and_save_plots(
    cell_df = cell_plot_df,
    config = config,
    difficulty_config = difficulty_config,
    output_plot = params$output_plot,
    expression_results = risultato$expression
  )
}
