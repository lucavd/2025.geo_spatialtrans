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

# Imposta numero di core da utilizzare
num_cores <- max(1, detectCores() - 10)
cat(sprintf("Utilizzo %d core per il calcolo parallelo...\n", num_cores))

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
# 7.1) Definisco un "profilo medio" per ciascun cluster
mean_expression_list <- list()
for (k in 1:k_cell_types) {
  mu <- rep(2, n_genes)
  start_idx <- (k-1)*10 + 1
  end_idx   <- min(k*10, n_genes)
  if (start_idx <= end_idx) {
    mu[start_idx:end_idx] <- mu[start_idx:end_idx] + 2
  }
  mean_expression_list[[k]] <- mu
}

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
  tryCatch({
    dist_mat <- as.matrix(dist(coords))
    cat("Matrice di distanza calcolata con metodo vettoriale.\n")
  }, error = function(e) {
    cat("Errore nel calcolo vettoriale, uso metodo parallelo...\n")
    chunk_size <- ceiling(N / min(100, num_cores * 4))
    row_indices <- 1:N
    row_chunks <- split(row_indices, ceiling(seq_along(row_indices) / chunk_size))
    cat(sprintf("Diviso il lavoro in %d blocchi\n", length(row_chunks)))
    calc_chunk <- function(chunk_rows) {
      result <- matrix(0, nrow=length(chunk_rows), ncol=N)
      x_coords <- coords$x
      y_coords <- coords$y
      for (i in seq_along(chunk_rows)) {
        row_idx <- chunk_rows[i]
        x_diff <- x_coords[row_idx] - x_coords
        y_diff <- y_coords[row_idx] - y_coords
        result[i,] <- sqrt(x_diff^2 + y_diff^2)
      }
      return(result)
    }
    dist_chunks <- mclapply(row_chunks, calc_chunk, mc.cores=num_cores)
    dist_mat <- matrix(0, nrow=N, ncol=N)
    for (i in seq_along(row_chunks)) {
      dist_mat[row_chunks[[i]],] <- dist_chunks[[i]]
    }
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

# Prova a eseguire senza parallelizzazione per il debugging
# expression_list <- lapply(seq_len(n_genes), function(g) {
expression_list <- mclapply(seq_len(n_genes), function(g) {
  # Vettore base con correlazione cell_cor
  base_vec <- mvrnorm(n=1, mu=rep(0,N), Sigma=cell_cor)

  # Aggiungo la media di cluster (vettorizzato)
  cluster_indices_numeric <- as.integer(cluster_labels)
  mu_values <- sapply(cluster_indices_numeric, function(cluster_id) mean_expression_list[[cluster_id]][g])
  gene_values <- base_vec + mu_values

  # Aggiungi un controllo per assicurarti che gene_values sia numerico
  if (!is.numeric(gene_values)) {
    cat(sprintf("Attenzione: gene_values non è numerico per il gene %d\n", g))
    return(rep(NA, N)) # Restituisci un vettore di NA per indicare il problema
  }
  return(gene_values)
}, mc.cores = num_cores)

# Verifica se ci sono elementi NULL nella lista
if (any(sapply(expression_list, is.null))) {
  stop("Errore: La lista di espressione contiene elementi NULL. Controlla la funzione interna di mclapply.")
}

# Combina la lista in una matrice
expression_data <- do.call(cbind, expression_list)

# Converto in conteggi
cat("Verifica classe di expression_data prima di exp(): ", class(expression_data), "\n")
if (!is.numeric(expression_data)) {
  stop("Errore: expression_data non è numerica prima di exp().")
}
expression_data <- exp(expression_data)
lambda_mat      <- expression_data
expression_data <- matrix(
  rpois(length(lambda_mat), lambda=lambda_mat),
  nrow=N, ncol=n_genes
)

# Dropout
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

ggplot2::ggsave(here::here("R/daniele/output.png"), plot = p, device = "png", dpi = 300 )

levels(result$intensity_cluster) <- paste0("cells_", letters[1:k_cell_types])
cat("Tabelle finali di cluster:\n")
table(result$intensity_cluster)
saveRDS(result, file = "data/simulated_image_correlation.rds")
cat("Script completato con successo.\n")
# Fine script
