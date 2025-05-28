#!/usr/bin/env Rscript
# Test rapido della simulazione full-size con parametri ridotti

cat("=== TEST RAPIDO SIMULAZIONE FULL-SIZE ===\n")
cat("Questo test usa parametri ridotti per verificare che tutto funzioni.\n\n")

# 1. Caricamento delle funzioni
cat("Caricamento delle funzioni...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# 2. Caricamento librerie
suppressPackageStartupMessages({
  library(png)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(ClusterR)
  library(sp)
  library(gstat)
  library(tictoc)
})

# Configurazione memoria
options(future.globals.maxSize = 50 * 1024^2)  # 50 GB per test
options(future.rng.onMisuse = "ignore")

# 3. Parametri TEST (ridotti)
cat("\nConfigurazione test con parametri ridotti...\n")
pixel_size_um <- 3000 / 300

cfg <- initialize_simulation_config(
  image_path          = "images/granuloma.png",
  output_path         = "results/test_fullsize.rds",
  output_plot         = "results/test_fullsize.png",
  n_genes             = 5000,      # Ridotto per test (vs 20k full)
  k_cell_types        = 5,         # Ridotto per test (vs 10 full)
  threshold_value     = 0.7,
  random_seed         = 42,
  pixel_size_um       = pixel_size_um,
  grid_mode           = TRUE,
  grid_resolution     = 10,        # Ridotto per test (vs 30 full)
  grid_spacing        = 0,
  use_fixed_grid      = TRUE,
  fixed_grid_width_mm = 6.5,
  fixed_grid_height_mm= 6.5
)

# Usa configurazione medium con parametri biologici
diff_cfg <- configure_difficulty_level("medium")
diff_cfg$cell_specific_params$library_size_params <- list(
  mean_library_size = 8000,
  library_size_cv = 0.3,
  spatial_effect_on_library = 0.1,
  cell_type_effect = TRUE
)

# 4. Pipeline rapida
tic("Test completo")

# Immagine e clustering
img_dat <- prepare_image(cfg$image_path, cfg$threshold_value)
clust <- cluster_image(
  img_df_thresh = img_dat$img_df_thresh,
  k_cell_types  = cfg$k_cell_types,
  random_seed   = cfg$random_seed
)

# Griglia
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

cat("\nPunti generati:", nrow(cell_df), "\n")
cat("Generazione espressione per", cfg$n_genes, "geni...\n")

# Test su subset ancora più piccolo
n_test_cells <- min(500, nrow(cell_df))
test_df <- cell_df[1:n_test_cells, ]

# Genera espressione
expr_result <- generate_expression_profiles(
  cell_df               = test_df,
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

# Verifica orientamento matrice
if (nrow(expr_result$expression) != cfg$n_genes) {
  expr_result$expression <- t(expr_result$expression)
}

# Statistiche rapide
umi_counts <- colSums(expr_result$expression)
cat("\n=== RISULTATI TEST ===\n")
cat("Dimensioni matrice:", dim(expr_result$expression), "\n")
cat("UMI medio:", round(mean(umi_counts)), "\n")
cat("UMI mediano:", round(median(umi_counts)), "\n")
cat("Range UMI:", round(min(umi_counts)), "-", round(max(umi_counts)), "\n")
cat("Sparsità:", round(mean(expr_result$expression == 0) * 100, 1), "%\n")
cat("Espressione mediana (non-zero):", round(median(expr_result$expression[expr_result$expression > 0]), 2), "\n")

toc()

# Validazione biologica rapida
cat("\n=== VALIDAZIONE BIOLOGICA RAPIDA ===\n")
realistic_cells <- sum(umi_counts >= 1000 & umi_counts <= 15000)
cat("Celle con UMI realistici (1k-15k):", realistic_cells, "/", length(umi_counts),
    "(", round(realistic_cells/length(umi_counts)*100, 1), "%)\n")

# Controllo distribuzione genica
gene_means <- rowMeans(expr_result$expression)
gene_zeros <- rowMeans(expr_result$expression == 0)

cat("\nDistribuzione espressione genica:\n")
cat("- Geni con espressione media < 0.1:", sum(gene_means < 0.1), 
    "(", round(sum(gene_means < 0.1)/cfg$n_genes*100, 1), "%)\n")
cat("- Geni con >90% zeri:", sum(gene_zeros > 0.9),
    "(", round(sum(gene_zeros > 0.9)/cfg$n_genes*100, 1), "%)\n")

# Test superato?
test_passed <- (mean(umi_counts) > 3000 && mean(umi_counts) < 12000) &&
               (median(expr_result$expression[expr_result$expression > 0]) > 0.1) &&
               (realistic_cells/length(umi_counts) > 0.5)

cat("\n=== TEST ", ifelse(test_passed, "SUPERATO ✓", "FALLITO ✗"), " ===\n")

if (test_passed) {
  cat("\nIl test è superato! Parametri biologici sembrano corretti.\n")
  cat("Puoi procedere con la simulazione full-size usando:\n")
  cat("  Rscript R/run_full_size_optimized.R\n")
} else {
  cat("\nIl test non è superato. Verifica i parametri biologici.\n")
  cat("Problemi rilevati:\n")
  if (mean(umi_counts) <= 3000 || mean(umi_counts) >= 12000) {
    cat("- UMI medio fuori range (atteso 3k-12k, ottenuto", round(mean(umi_counts)), ")\n")
  }
  if (median(expr_result$expression[expr_result$expression > 0]) <= 0.1) {
    cat("- Espressione mediana troppo bassa (atteso >0.1, ottenuto", 
        round(median(expr_result$expression[expr_result$expression > 0]), 3), ")\n")
  }
  if (realistic_cells/length(umi_counts) <= 0.5) {
    cat("- Troppe poche celle con UMI realistici (atteso >50%, ottenuto",
        round(realistic_cells/length(umi_counts)*100, 1), "%)\n")
  }
}

# Cleanup
rm(expr_result, umi_counts)
gc()

cat("\nMemoria utilizzata (MB):", round(sum(gc()[,2]), 1), "\n")