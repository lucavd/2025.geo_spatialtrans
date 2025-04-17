#!/usr/bin/env Rscript
# Script per simulazione full‑size Visium HD su granuloma.png

# 1. Caricamento delle funzioni
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# 2. Caricamento librerie necessarie
library(png)
library(ggplot2)
library(dplyr)
library(Matrix)
# ClusterR for fast kmeans++ implementation (provides KMeans_rcpp)
# ClusterR for fast kmeans++ implementation (provides KMeans_rcpp)
library(ClusterR)
library(sp)
library(gstat)

## 3. Parametri di simulazione
pixel_size_um <- 3000 / 300  # 300 px = 3 mm ⇒ 10 µm/px
cfg <- initialize_simulation_config(
  image_path          = "images/granuloma.png",
  output_path         = "results/visiumHD_full.rds",
  output_plot         = "results/visiumHD_full.png",
  n_genes             = 200,
  k_cell_types        = 10,
  threshold_value     = 0.7,
  random_seed         = 42,
  pixel_size_um       = pixel_size_um,
  grid_mode           = TRUE,
  grid_resolution     = 10,
  grid_spacing        = 0,
  use_fixed_grid      = TRUE,
  fixed_grid_width_mm = 6.5,
  fixed_grid_height_mm= 6.5
)
diff_cfg <- configure_difficulty_level("medium")

# 4. Preparazione immagine e clustering
img_dat <- prepare_image(cfg$image_path, cfg$threshold_value)
clust  <- cluster_image(
  img_df_thresh = img_dat$img_df_thresh,
  k_cell_types  = cfg$k_cell_types,
  random_seed   = cfg$random_seed
)

# 5. Creazione griglia fissa centrata
cell_df <- create_sampling_grid(
  img_df_thresh         = clust,
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

# 6. Generazione chunked dei profili di espressione
chunk_size <- 5000
spots      <- seq_len(nrow(cell_df))
chunks     <- split(spots, ceiling(seq_along(spots) / chunk_size))
sparse_list <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  idx     <- chunks[[i]]
  sub_df  <- cell_df[idx, ]
  expr_r  <- generate_expression_profiles(
    cell_df               = sub_df,
    n_genes               = cfg$n_genes,
    k_cell_types          = cfg$k_cell_types,
    marker_params         = diff_cfg$marker_params,
    spatial_params        = diff_cfg$spatial_params,
    dropout_params        = diff_cfg$dropout_params,
    cell_specific_params  = diff_cfg$cell_specific_params,
    use_spatial_correlation = TRUE,
    correlation_method    = "grf",
    random_seed           = cfg$random_seed
  )
  sparse_list[[i]] <- Matrix(expr_r$expression, sparse = TRUE)
  message(sprintf("Chunk %d/%d completato", i, length(chunks)))
}
full_expr <- do.call(rbind, sparse_list)

# 7. Salvataggio dei risultati
risultato <- list(
  expression        = full_expr,
  coordinates       = as.matrix(cell_df[, c("x", "y")]),
  intensity_cluster = cell_df$intensity_cluster
)
dir.create(dirname(cfg$output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(risultato, cfg$output_path)

# 8. Visualizzazione e salvataggio del plot
generate_and_save_plots(
  cell_df            = cell_df,
  config             = cfg,
  difficulty_config  = diff_cfg,
  output_plot        = cfg$output_plot,
  expression_results = full_expr
)
message("Simulazione full‑size completata. Risultati in: ", cfg$output_path)
