#!/usr/bin/env Rscript
# Script ottimizzato per simulazione full-size Visium HD biologicamente realistica

# 1. Caricamento delle funzioni
cat("Caricamento delle funzioni...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# 2. Caricamento librerie necessarie
cat("Caricamento librerie...\n")
library(png)
library(ggplot2)
library(dplyr)
library(Matrix)
library(ClusterR)
library(sp)
library(gstat)
library(tictoc)

# Configurazione per la memoria
options(future.globals.maxSize = 200 * 1024^2)  # 200 GB per simulazioni full-size
options(future.rng.onMisuse = "ignore")

## 3. Parametri di simulazione biologicamente realistici
cat("Configurazione simulazione full-size biologicamente realistica...\n")
pixel_size_um <- 3000 / 300  # 300 px = 3 mm ⇒ 10 µm/px

cfg <- initialize_simulation_config(
  image_path          = "images/granuloma.png",
  output_path         = "results/visiumHD_biological.rds",
  output_plot         = "results/visiumHD_biological.png",
  n_genes             = 20000,     # Numero realistico di geni per Visium HD
  k_cell_types        = 10,
  threshold_value     = 0.7,
  random_seed         = 42,
  pixel_size_um       = pixel_size_um,
  grid_mode           = TRUE,
  grid_resolution     = 30,        # Griglia fitta per Visium HD
  grid_spacing        = 0,
  use_fixed_grid      = TRUE,
  fixed_grid_width_mm = 6.5,
  fixed_grid_height_mm= 6.5
)

# Configurazione difficoltà con parametri biologicamente validati
diff_cfg <- configure_difficulty_level("medium")

# Sovrascrivi con parametri biologicamente realistici
# Library size parameters
diff_cfg$cell_specific_params$library_size_params <- list(
  mean_library_size = 8000,      # Validato per Visium HD
  library_size_cv = 0.3,         # 30% CV tipico
  spatial_effect_on_library = 0.1,  # Leggero effetto spaziale
  cell_type_effect = TRUE          # Variazione per tipo cellulare
)

# Dropout parameters ottimizzati
diff_cfg$dropout_params$dropout_range <- c(0.4, 0.6)  # Range realistico
diff_cfg$dropout_params$k <- 1.5  # Forma della curva di dropout

# Dispersion parameters
diff_cfg$dispersion_params$dispersion_range <- c(10.0, 5.0)  # Validati

# Spatial correlation più realistica
diff_cfg$spatial_params$correlation_length <- diff_cfg$spatial_params$correlation_length * 2.5
diff_cfg$spatial_params$noise_level <- diff_cfg$spatial_params$noise_level * 1.5

# Marker gene parameters per maggiore variabilità biologica
diff_cfg$marker_params$strength_mean <- diff_cfg$marker_params$strength_mean * 1.2
diff_cfg$marker_params$strength_sd <- diff_cfg$marker_params$strength_sd * 1.8

# Flag per validazione (può essere disabilitato per risparmiare memoria)
validate <- TRUE

# --- NUOVO: flag per immagine sintetica ------------------------------------
use_synthetic_image <- TRUE   # Imposta a FALSE per usare immagine reale
synthetic_complexity <- 2     # 1=blob, 2=Voronoi, 3=mix

# 4. Preparazione immagine e clustering
if (use_synthetic_image) {
  tic("Generazione immagine sintetica")
  syn <- generate_synthetic_tissue(
    width_px  = 6800,
    height_px = 6500,
    complexity = synthetic_complexity,
    seed = cfg$random_seed,
    output_path = "results/synthetic_tissue.png"
  )
  # Normalizza valori 0-255 -> 0-1 e prepara dataframe come atteso
  img_array_norm <- syn$img_matrix / 255
  img_df <- syn$img_df
  img_df$value <- img_df$intensity / 255
  img_df_thresh <- img_df[img_df$value < cfg$threshold_value, c("x", "y", "value")]
  
  img_dat <- list(
    width = 6800,
    height = 6500,
    img_df_thresh = img_df_thresh,
    img_array = img_array_norm
  )
  cat("Immagine sintetica salvata in results/synthetic_tissue.png\n")
  toc()
} else {
  tic("Preparazione immagine")
  img_dat <- prepare_image(cfg$image_path, cfg$threshold_value)
  cat("Dimensioni immagine:", img_dat$width, "x", img_dat$height, "pixel\n")
  toc()
}

clust <- cluster_image(
  img_df_thresh = img_dat$img_df_thresh,
  k_cell_types  = cfg$k_cell_types,
  random_seed   = cfg$random_seed,
  clustering_method = "dbscan_graph"  # Usa pipeline DBSCAN + Graph clustering
)
# Rimuovi punti senza cluster assegnato (NA)
clust <- clust[!is.na(clust$intensity_cluster), ]

# 5. Creazione griglia
tic("Creazione griglia")
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
cat("Punti generati:", nrow(cell_df), "\n")
toc()

# --- CLUSTER MORPHOLOGY PLOT (minimal) ---
cat("Generazione figura di morfologia dei cluster...\n")
library(ggplot2)
library(dplyr)

# Calcola distanza media locale per ogni punto (rispetto agli altri punti dello stesso cluster)
calc_mean_distance <- function(df, n_neighbors = 15) {
  # Per grandi dataset: calcola solo sui 15 vicini più vicini
  if (!requireNamespace("FNN", quietly = TRUE)) install.packages("FNN")
  library(FNN)
  df <- df %>% mutate(cluster = as.factor(intensity_cluster))
  df$mean_distance <- NA
  for (cl in levels(df$cluster)) {
    idx <- which(df$cluster == cl)
    if (length(idx) > n_neighbors) {
      nn <- get.knn(df[idx, c("x", "y")], k = n_neighbors)
      df$mean_distance[idx] <- rowMeans(nn$nn.dist)
    } else {
      df$mean_distance[idx] <- NA
    }
  }
  df
}

cell_df_plot <- calc_mean_distance(cell_df)

p <- ggplot(cell_df_plot, aes(x = x, y = y, color = mean_distance)) +
  geom_point(size = 1) +
  facet_wrap(~ intensity_cluster, ncol = 4) +
  scale_color_viridis_c(option = "inferno", direction = -1) +
  theme_minimal() +
  labs(title = "Morfologia dei cluster (distanza media locale)",
       x = "x", y = "y", color = "Distanza media") +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

dir.create("plots", showWarnings = FALSE)
ggsave("plots/cluster_morphology.png", p, width = 10, height = 4)
# --- FINE CLUSTER MORPHOLOGY PLOT ---

# 6. Generazione ottimizzata dei profili di espressione
tic("Generazione espressione con parametri biologici")

# Dimensione chunk ottimizzata per 20k geni
chunk_size <- 2000  # Ridotto per gestire più geni
spots <- seq_len(nrow(cell_df))
chunks <- split(spots, ceiling(seq_along(spots) / chunk_size))
cat("Divisione in", length(chunks), "chunks per gestire", cfg$n_genes, "geni\n")

# Lista per risultati
sparse_list <- vector("list", length(chunks))

# Progress tracking
pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)

for (i in seq_along(chunks)) {
  setTxtProgressBar(pb, i)
  
  idx <- chunks[[i]]
  sub_df <- cell_df[idx, ]
  
  # Garbage collection
  if (i %% 5 == 0) gc()
  
  # Perturbazione delle posizioni
  set.seed(cfg$random_seed + i*100)
  perturbed_df <- sub_df
  grid_cell_size <- 1 / cfg$grid_resolution
  jitter_scale <- grid_cell_size * 0.45
  
  perturbed_df$x <- perturbed_df$x + runif(nrow(perturbed_df), -jitter_scale, jitter_scale)
  perturbed_df$y <- perturbed_df$y + runif(nrow(perturbed_df), -jitter_scale, jitter_scale)
  
  # Variazione parametri spaziali per chunk
  chunk_spatial_params <- diff_cfg$spatial_params
  chunk_spatial_params$correlation_length <- chunk_spatial_params$correlation_length * runif(1, 0.9, 1.1)
  
  # Genera espressione con tutti i parametri biologici
  expr_r <- generate_expression_profiles(
    cell_df               = perturbed_df,
    n_genes               = cfg$n_genes,
    k_cell_types          = cfg$k_cell_types,
    marker_params         = diff_cfg$marker_params,
    spatial_params        = chunk_spatial_params,
    dropout_params        = diff_cfg$dropout_params,
    cell_specific_params  = diff_cfg$cell_specific_params,
    use_spatial_correlation = TRUE,
    correlation_method    = "grf",
    random_seed           = cfg$random_seed + i
  )
  
  # Converti in matrice sparsa
  if (nrow(expr_r$expression) == cfg$n_genes) {
    sparse_list[[i]] <- Matrix(expr_r$expression, sparse = TRUE)
  } else if (ncol(expr_r$expression) == cfg$n_genes) {
    sparse_list[[i]] <- Matrix(t(expr_r$expression), sparse = TRUE)
  }
}

close(pb)
cat("\n")

# Combina chunks
cat("Combinazione delle matrici sparse...\n")
if (nrow(sparse_list[[1]]) == cfg$n_genes) {
  full_expr <- Matrix(0, nrow = cfg$n_genes, ncol = nrow(cell_df), sparse = TRUE)
  col_offset <- 0
  for (i in seq_along(sparse_list)) {
    chunk_size <- ncol(sparse_list[[i]])
    col_indices <- col_offset + (1:chunk_size)
    full_expr[, col_indices] <- sparse_list[[i]]
    col_offset <- col_offset + chunk_size
  }
} else {
  full_expr <- do.call(rbind, sparse_list)
}

cat("Dimensione matrice finale:", dim(full_expr), "\n")
toc()

# 7. Validazione biologica della matrice finale
cat("\n=== VALIDAZIONE BIOLOGICA MATRICE FINALE ===\n")
tic("Validazione biologica")
# Carica la funzione se non già disponibile
if (!exists("validate_biological_plausibility")) {
  source("R/functions/06j_expression_generation.R")
}
full_expr <- validate_biological_plausibility(full_expr)
toc()

# 8. Statistiche rapide post-validazione
cat("\n=== STATISTICHE RAPIDE POST-VALIDAZIONE ===\n")
umi_sample <- colSums(full_expr[, sample(ncol(full_expr), min(1000, ncol(full_expr)))])
cat("UMI medio (campione):", round(mean(umi_sample)), "\n")
cat("UMI mediano (campione):", round(median(umi_sample)), "\n")
cat("Sparsità:", round(mean(full_expr == 0) * 100, 1), "%\n")

# 9. Salvataggio risultati
tic("Salvataggio risultati")
if (is.null(rownames(full_expr))) {
  rownames(full_expr) <- paste0("gene_", seq_len(nrow(full_expr)))
}

risultato <- list(
  expression = full_expr,
  coordinates = cell_df[, c("x", "y")],
  intensity_cluster = factor(cell_df$intensity_cluster,
                           levels = sort(unique(as.numeric(cell_df$intensity_cluster)))),
  parameters = list(
    marker_params = diff_cfg$marker_params,
    spatial_params = diff_cfg$spatial_params,
    dropout_params = diff_cfg$dropout_params,
    cell_specific_params = diff_cfg$cell_specific_params,
    n_genes = cfg$n_genes,
    library_size_mean = diff_cfg$cell_specific_params$library_size_params$mean_library_size
  )
)

dir.create(dirname(cfg$output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(risultato, cfg$output_path)
toc()

# 10. Plot principale
tic("Generazione plot principale")
generate_and_save_plots(
  cell_df = cell_df,
  config = cfg,
  difficulty_config = diff_cfg,
  output_plot = cfg$output_plot
)
toc()

# 11. Validazione biologica report
cat("\n=== SIMULAZIONE FULL-SIZE COMPLETATA ===\n")
cat("Risultati salvati in:", cfg$output_path, "\n")
cat("- Numero di geni:", cfg$n_genes, "\n")
cat("- Numero di celle:", nrow(cell_df), "\n")

# Validazione automatica post-simulazione
cat("\n=== VALIDAZIONE BIOLOGICA AUTOMATICA ===\n")
source("R/biological_validation_report.R")
generate_biological_validation_report(cfg$output_path, output_dir = "R/validation", report_name = "visiumHD_biological")
cat("\nReport di validazione salvato in: R/validation\n")
cat("- Library size medio configurato:", diff_cfg$cell_specific_params$library_size_params$mean_library_size, "\n")

cat("\n=== NOTE BIOLOGICHE ===\n")
cat("Questa è una simulazione full-size con parametri biologicamente realistici:\n")
cat("- 20,000 geni (tipico per esperimenti Visium HD)\n")
cat("- Library size ~8,000 UMI/cella (validato per Visium HD)\n")
cat("- Distribuzione genica: 40% non-espressi, 35% low, 20% medium, 5% high expression\n")
cat("- Dropout modeling realistico basato su espressione media\n")
cat("- Correlazione spaziale e variabilità biologica incluse\n")

if (validate) {
  cat("\nLa validazione biologica verificherà che i risultati siano plausibili.\n")
}