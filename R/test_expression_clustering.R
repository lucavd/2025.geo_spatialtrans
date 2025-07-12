#!/usr/bin/env Rscript
# Test script per clustering basato su espressione

# 1. Caricamento delle funzioni
cat("Caricamento delle funzioni...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# Carica funzione di validazione
source("R/validation/validate_expression_clustering.R")

# 2. Caricamento librerie necessarie
cat("Caricamento librerie...\n")
library(png)
library(ggplot2)
library(dplyr)
library(Matrix)
library(tictoc)
library(sp)
library(gstat)

# Configurazione per la memoria
options(future.globals.maxSize = 50 * 1024^2)  # 50 GB
options(future.rng.onMisuse = "ignore")

# 3. Parametri di test (più piccoli per velocità)
cat("Configurazione test expression-based clustering...\n")
pixel_size_um <- 3000 / 300  # 10 µm/px

cfg <- initialize_simulation_config(
  image_path          = "images/granuloma.png",
  output_path         = "results/test_expression_clustering.rds",
  output_plot         = "results/test_expression_clustering.png",
  n_genes             = 500,      # Molto ridotto per test
  k_cell_types        = 5,
  threshold_value     = 0.7,
  random_seed         = 42,
  pixel_size_um       = pixel_size_um,
  grid_mode           = TRUE,
  grid_resolution     = 50,        # Ridotto per test
  grid_spacing        = 0,
  use_fixed_grid      = TRUE,
  fixed_grid_width_mm = 2.0,       # Molto ridotto per test
  fixed_grid_height_mm= 2.0
)

# Configurazione difficoltà
diff_cfg <- configure_difficulty_level("medium")

# Parametri biologicamente realistici (semplificati)
diff_cfg$cell_specific_params$library_size_params <- list(
  mean_library_size = 5000,       # Ridotto per test
  library_size_cv = 0.25,
  spatial_effect_on_library = 0.1,
  cell_type_effect = TRUE
)

diff_cfg$library_size_params <- diff_cfg$cell_specific_params$library_size_params
diff_cfg$dropout_params$dropout_range <- c(0.3, 0.5)
diff_cfg$dispersion_params$dispersion_range <- c(8.0, 4.0)

# Flag per immagine sintetica
use_synthetic_image <- TRUE
synthetic_complexity <- 2

# 4. Preparazione immagine
if (use_synthetic_image) {
  tic("Generazione immagine sintetica")
  syn <- generate_synthetic_tissue(
    width_px  = 800,      # Molto ridotto per test
    height_px = 800,
    complexity = synthetic_complexity,
    seed = cfg$random_seed,
    output_path = "results/test_synthetic_tissue.png"
  )
  
  img_array_norm <- syn$img_matrix / 255
  img_df <- syn$img_df
  img_df$value <- img_df$intensity / 255
  img_df_thresh <- img_df[img_df$value < cfg$threshold_value, c("x", "y", "value")]
  
  img_dat <- list(
    width = 800,
    height = 800,
    img_df_thresh = img_df_thresh,
    img_array = img_array_norm
  )
  toc()
} else {
  tic("Preparazione immagine")
  img_dat <- prepare_image(cfg$image_path, cfg$threshold_value)
  toc()
}

# 5. Creazione griglia SENZA clustering preliminare
tic("Creazione griglia")

# Aggiungi una colonna intensity_cluster all'img_df_thresh per compatibilità
img_df_thresh_with_cluster <- img_dat$img_df_thresh
img_df_thresh_with_cluster$intensity_cluster <- 1

# Creiamo una griglia diretta dall'immagine thresholded
cell_df <- create_sampling_grid(
  img_df_thresh         = img_df_thresh_with_cluster,
  img_array             = img_dat$img_array,
  img_width             = img_dat$width,
  img_height            = img_dat$height,
  grid_mode             = cfg$grid_mode,
  grid_resolution       = cfg$grid_resolution,
  grid_spacing          = cfg$grid_spacing,
  use_fixed_grid        = cfg$use_fixed_grid,
  fixed_grid_width_mm   = cfg$fixed_grid_width_mm,
  fixed_grid_height_mm  = cfg$fixed_grid_height_mm,
  pixel_size_um         = cfg$pixel_size_um,
  threshold_value       = cfg$threshold_value,
  random_seed           = cfg$random_seed
)

# Aggiungi una colonna intensity_cluster temporanea per compatibilità
cell_df$intensity_cluster <- 1

cat("Punti generati:", nrow(cell_df), "\n")
toc()

# 6. NUOVO: Expression-based clustering
tic("Expression-based clustering")

# Prepara parametri per expression clustering
expression_params <- list(
  marker_params = diff_cfg$marker_params,
  spatial_params = diff_cfg$spatial_params,
  dropout_params = diff_cfg$dropout_params,
  cell_specific_params = diff_cfg$cell_specific_params
)

# Applica clustering basato su espressione
clustering_result <- expression_based_clustering(
  cell_df = cell_df,
  n_genes = cfg$n_genes,
  k_cell_types = cfg$k_cell_types,
  expression_params = expression_params,
  random_seed = cfg$random_seed,
  clustering_method = "louvain",
  n_pcs = 20,              # Ridotto per test
  k_neighbors = 8,         # Ridotto per test
  resolution = 1.5         # Aumentato per più cluster
)

# Aggiorna cell_df con i nuovi cluster
cell_df <- clustering_result$cell_df
expr_matrix <- clustering_result$expression_matrix

cat("Clustering completato!\n")
toc()

# 7. Validazione del clustering
tic("Validazione clustering")
validation_results <- validate_expression_clustering(
  clustering_result = clustering_result,
  cell_df = cell_df,
  output_dir = "R/validation",
  prefix = "test_expression_clustering"
)
toc()

# 8. Salvataggio risultati
tic("Salvataggio risultati")
if (is.null(rownames(expr_matrix))) {
  rownames(expr_matrix) <- paste0("gene_", seq_len(nrow(expr_matrix)))
}

risultato <- list(
  expression = expr_matrix,
  coordinates = cell_df[, c("x", "y")],
  intensity_cluster = factor(cell_df$intensity_cluster),
  clustering_method = "expression_based",
  pca_coords = clustering_result$pca_coords,
  validation_results = validation_results,
  parameters = list(
    marker_params = expression_params$marker_params,
    spatial_params = expression_params$spatial_params,
    dropout_params = expression_params$dropout_params,
    cell_specific_params = expression_params$cell_specific_params,
    n_genes = cfg$n_genes,
    clustering_method = "expression_based"
  )
)

dir.create(dirname(cfg$output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(risultato, cfg$output_path)
toc()

# 9. Plot principale
tic("Generazione plot principale")
generate_and_save_plots(
  cell_df = cell_df,
  config = cfg,
  difficulty_config = diff_cfg,
  output_plot = cfg$output_plot
)
toc()

# 10. Statistiche riassuntive
cat("\n=== CLUSTERING BASATO SU ESPRESSIONE COMPLETATO ===\n")
cat("Risultati salvati in:", cfg$output_path, "\n")
cat("- Numero di geni:", cfg$n_genes, "\n")
cat("- Numero di celle:", nrow(cell_df), "\n")
cat("- Numero di cluster trovati:", length(unique(cell_df$intensity_cluster)), "\n")
cat("- Silhouette score medio:", round(validation_results$silhouette_score$average, 3), "\n")
cat("- Coerenza spaziale media:", round(validation_results$spatial_coherence$average, 3), "\n")
cat("- Ratio separazione espressione:", round(validation_results$expression_separation$separation_ratio, 3), "\n")

cat("\n=== VANTAGGI CLUSTERING BASATO SU ESPRESSIONE ===\n")
cat("1. Cluster basati su somiglianze biologiche reali\n")
cat("2. Forme irregolari naturalmente derivate dall'espressione\n")
cat("3. Validazione quantitativa della qualità biologica\n")
cat("4. Maggiore realismo per benchmark di metodi\n")

cat("\nValidazione salvata in: R/validation/\n")
cat("- File di validazione: test_expression_clustering_validation_report.md\n")
cat("- Visualizzazioni PCA e spaziali generate\n")