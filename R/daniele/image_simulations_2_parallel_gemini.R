#!/usr/bin/env Rscript

# Script: simulazione di trascrittomica spaziale
#  dove le cellule sono raggruppate in base a:
#   - intensità di grigio (cluster "cell type")
#   - distanza spaziale (maggiore correlazione se vicine)
#   - soglia per escludere i pixel troppo chiari/scuri.

# ========================================
# Caricamento librerie
# ========================================
library(imager)
library(tidyverse)
library(MASS)
library(cluster)
library(parallel)

# Imposta numero di core da utilizzare (ottimizzato per la stabilità)
num_cores <- min(8, max(1, detectCores() - 2)) # Limite più conservativo
cat(sprintf("Utilizzo %d core per il calcolo parallelo con %d core disponibili...\n", num_cores, detectCores()))

# ========================================
# 1) Parametri principali
# ========================================
image_path            <- here::here("R/daniele/colon.png")
n_cells               <- 30000
n_genes               <- 10
k_cell_types          <- 3
use_spatial_correlation <- TRUE
threshold_value       <- 0.7

# ========================================
# 2) Caricamento e conversione immagine
# ========================================
img <- load.image(image_path)
if (spectrum(img) == 4) {
  img <- rm.alpha(img)
}
if (spectrum(img) == 3) {
  img <- grayscale(img)
}
img <- squeeze(img)
cat("Dimensioni dell'immagine cimg:\n")
print(dim(img))

# ========================================
# 3) Estraggo intensità e coordinate di tutti i pixel
# ========================================
img_array <- as.array(img)
w <- nrow(img_array)
h <- ncol(img_array)
img_cimg <- as.cimg(img_array, dims = c(w, h, 1, 1))
dim(img_cimg)
spectrum(img_cimg)
img_df <- as.data.frame(img_cimg)
img_df <- img_df %>%
  dplyr::select(x, y, value)
cat("Range intensità immagine:\n")
print(summary(img_df$value))

# ========================================
# 4) Applico la soglia per escludere pixel indesiderati
# ========================================
img_df_thresh <- img_df %>%
  filter(value < threshold_value)
cat(sprintf("Pixel totali: %d\n", nrow(img_df)))
cat(sprintf("Pixel dopo la soglia (< %.2f): %d\n", threshold_value, nrow(img_df_thresh)))

# ========================================
# 5) Cluster dei pixel in base all'intensità
#    (solo i pixel rimanenti dopo la soglia)
# ========================================
set.seed(123)
km_intensity <- kmeans(img_df_thresh$value, centers = k_cell_types, nstart = 10)
img_df_thresh$intensity_cluster <- factor(km_intensity$cluster)

# ========================================
# 6) Campiono n_cells pixel (implementazione vettorizzata)
#    proporzionale alla dimensione di ogni cluster
# ========================================
clust_counts <- table(img_df_thresh$intensity_cluster)
clust_freq   <- clust_counts / sum(clust_counts)
cells_per_cluster <- round(clust_freq * n_cells)

set.seed(123)
cell_df_list <- lapply(levels(img_df_thresh$intensity_cluster), function(clust) {
  n_sub <- cells_per_cluster[clust]
  df_sub <- img_df_thresh[img_df_thresh$intensity_cluster == clust, ]
  if (nrow(df_sub) > 0) {
    idx_sub <- sample(seq_len(nrow(df_sub)), min(n_sub, nrow(df_sub)), replace = FALSE)
    return(df_sub[idx_sub, ])
  } else {
    return(NULL)
  }
})
cell_df <- do.call(rbind, cell_df_list)
cat(sprintf("Celle totali campionate: %d\n", nrow(cell_df)))

# ========================================
# 7) Generazione dei profili di espressione
# ========================================
# 7.1) Definisco un "profilo medio" per ciascun cluster (vettorizzato)
mean_expression_list <- lapply(1:k_cell_types, function(k) {
  mu <- rep(2, n_genes)
  start_idx <- (k - 1) * 10 + 1
  end_idx   <- min(k * 10, n_genes)
  if (start_idx <= end_idx) {
    mu[start_idx:end_idx] <- mu[start_idx:end_idx] + 2
  }
  return(mu)
})

# 7.2) Costruisco matrice di correlazione iniziale (intra vs inter cluster)
cat("Inizio costruzione matrice di correlazione iniziale...\n")
rho_intra <- 0.6
rho_inter <- 0.1
N <- nrow(cell_df)
cluster_labels <- cell_df$intensity_cluster
cluster_indices <- as.integer(cluster_labels)
equal_clusters <- outer(cluster_indices, cluster_indices, "==")
cell_cor <- matrix(0, nrow=N, ncol=N)
cell_cor[equal_clusters] <- rho_intra
cell_cor[!equal_clusters] <- rho_inter
diag(cell_cor) <- 1
cat("Matrice di correlazione iniziale costruita.\n")

# 7.3) Opzionale: aggiungo contributo dalla distanza spaziale
if (use_spatial_correlation) {
  cat("Inizio calcolo matrice di distanza...\n")
  coords <- cell_df[, c("x", "y")]
  # Tenta il metodo vettoriale, dovrebbe essere più efficiente se funziona
  tryCatch({
    dist_mat <- as.matrix(dist(coords))
    cat("Matrice di distanza calcolata con metodo vettoriale.\n")
  }, error = function(e) {
    cat("Errore nel calcolo vettoriale, considera di ridurre n_cells se persiste.\n")
    stop(e) # Interrompi se il metodo vettoriale fallisce, per ora
    # Se proprio necessario, si potrebbe implementare una versione a blocchi non parallela
    # per evitare i problemi di memoria della parallelizzazione manuale.
  })
  cat("Calcolo della similarità spaziale...\n")
  range_param <- 30
  spatial_sim <- exp(- dist_mat / range_param)
  alpha <- 0.5
  beta  <- 0.5
  cell_cor <- alpha * cell_cor + beta * spatial_sim
  diag(cell_cor) <- 1
  cat("Correlazione spaziale completata.\n")
}

# 7.4) Genero la matrice di espressione (parallelizzato e con vettorizzazione interna)
cat("Generazione dei profili di espressione (parallela)...\n")
N <- nrow(cell_df)
cluster_labels <- cell_df$intensity_cluster

expression_data <- t(mclapply(seq_len(n_genes), function(g) {
  # Vettore base con correlazione cell_cor
  base_vec <- mvrnorm(n=1, mu=rep(0,N), Sigma=cell_cor)

  # Aggiungo la media di cluster (completamente vettorizzato)
  cluster_indices_numeric <- as.integer(cluster_labels)
  mean_expression_vector <- sapply(cluster_indices_numeric, function(cluster_id) mean_expression_list[[cluster_id]][g])
  gene_values <- base_vec + mean_expression_vector
  return(gene_values)
}, mc.cores = num_cores))

# Converto in conteggi
expression_data <- exp(expression_data)
lambda_mat      <- expression_data
expression_data <- matrix(
  rpois(length(lambda_mat), lambda=lambda_mat),
  nrow=N, ncol=n_genes
)

# Dropout (vettorizzato)
dropout_rate <- 0.3
num_zero     <- round(length(expression_data) * dropout_rate)
set.seed(123)
zero_idx <- sample.int(length(expression_data), num_zero)
expression_data[zero_idx] <- 0

cat("Generazione dei profili completata (parallela).\n")

# Risultato finale
result <- list(
  coordinates        = cell_df[, c("x","y")],
  intensity_cluster  = cell_df$intensity_cluster,
  expression         = expression_data,
  threshold_used     = threshold_value
)

# ========================================
# 8) Visualizzazione rapida
# ========================================
p <- ggplot(cell_df, aes(x=x, y=y, color=intensity_cluster)) +
  geom_point(size=2, alpha=0.7) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  labs(title = sprintf("Distribuzione spaziale (threshold=%.2f)", threshold_value),
       color = "Cell Type")

# p

ggplot2::ggsave(here::here("R/daniele/output.png"), plot = p, device = "png", dpi = 600 )

levels(result$intensity_cluster) <- paste0("cells_", letters[1:k_cell_types])
cat("Tabelle finali di cluster:\n")
table(result$intensity_cluster)
saveRDS(result, file = "data/simulated_image_correlation.rds")
cat("Script completato con successo.\n")
# Fine script
