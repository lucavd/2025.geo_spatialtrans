#!/usr/bin/env Rscript
# Simplified simulation runner for spatial transcriptomics

cat("Loading functions...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

cat("Loading libraries...\n")
library(png)
library(ggplot2)
library(dplyr)
library(Matrix)
library(sp)
library(gstat)
library(tictoc)

# Increase memory limit for large simulations
options(future.globals.maxSize = 100 * 1024^2)

cat("Configuring simulation parameters...\n")
image_path <- "images/granuloma.png"
output_rds <- "results/simple_simulation.rds"
output_plot <- "results/simple_simulation.png"

# Ensure output directory exists
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)

cat("Starting simulation...\n")
tic("Total simulation time")
sim_results <- simulate_spatial_transcriptomics(
  image_path = image_path,
  output_path = output_rds,
  output_plot = output_plot,
  n_genes = 2000,
  k_cell_types = 10,
  threshold_value = 0.7,
  random_seed = 42,
  pixel_size_um = 10,
  grid_mode = TRUE,
  grid_resolution = 30,
  grid_spacing = 0,
  use_fixed_grid = TRUE,
  fixed_grid_width_mm = 6.5,
  fixed_grid_height_mm = 6.5,
  difficulty_level = "medium",
  use_spatial_correlation = TRUE,
  correlation_method = "grf",
  # Parametri per ottenere valori biologicamente realistici
  library_size_params = list(
    mean_library_size = 8000,      # Biologicamente realistico per Visium HD
    library_size_cv = 0.3,         # 30% CV è tipico per spatial transcriptomics
    spatial_effect_on_library = 0.1,  # Leggero effetto spaziale
    cell_type_effect = TRUE         # Abilitato per realismo biologico
  )
)
toc()
## Ensure expression matrix is genes x cells with proper rownames
expr_mat <- sim_results$expression
# Transpose if rows (cells) exceed columns (genes)
if (nrow(expr_mat) > ncol(expr_mat)) {
  expr_mat <- t(expr_mat)
}
# Set gene rownames if missing
if (is.null(rownames(expr_mat))) {
  rownames(expr_mat) <- paste0("gene_", seq_len(nrow(expr_mat)))
}
sim_results$expression <- expr_mat

cat("Generating validation plots...\n")
validation_dir <- "plots/validazione_simple"
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

  tic("Validation plots time")
  # Identify top marker genes per cluster for spatial plotting
  marker_list <- identify_marker_genes(sim_results, n_markers = 3, min_ratio = 1.5)
  marker_genes <- unique(unlist(marker_list))
  cat("Marker genes selected for validation:", paste(marker_genes, collapse = ", "), "\n")
  generate_validation_plots(
    sim_results = sim_results,
    marker_genes = marker_genes,
    output_dir = validation_dir,
    file_prefix = "simple",
    skip_dim_reduction = TRUE
  )
  toc()

cat("Simulation and validation complete.\n")
cat("Results saved to:", output_rds, "\n")
cat("Main plot saved to:", output_plot, "\n")
cat("Validation plots saved to:", validation_dir, "\n")
