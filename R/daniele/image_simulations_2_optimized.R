#!/usr/bin/env Rscript

# Script: simulazione di trascrittomica spaziale (ottimizzato)
# 1) Caricamento immagine e threshold
# 2) Clustering con kmeans++
# 3) Sampling celle
# 4) Generazione dei profili
# 5) Visualizzazione

library(imager)
library(tidyverse)
library(MASS)
library(cluster)
library(ClusterR)
library(future)
library(future.apply)
library(tictoc)

tictoc::tic()

# =======================
# 1) Parametri principali
# =======================
image_path            <- here::here("R/daniele/colon.png")
n_cells               <- 3000
n_genes               <- 10
k_cell_types          <- 3
use_spatial_correlation <- TRUE
threshold_value       <- 0.7

# ========================
# 2) Caricamento immagine e conversione in data frame sogliato
# ========================
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

# Converto manualmente l'immagine in array e poi in cimg (x,y,cc,t)
img_array <- as.array(img)
w <- nrow(img_array)  # ad esempio 350
h <- ncol(img_array)  # ad esempio 271
img_cimg <- as.cimg(img_array, dims = c(w, h, 1, 1))

cat("Dimensioni img_cimg:\n")
print(dim(img_cimg))
cat("Numero canali:\n")
print(spectrum(img_cimg))

# Converto in data frame (x, y, value)
img_df <- as.data.frame(img_cimg) %>%
  dplyr::select(x, y, value)

cat("Range intensità immagine:\n")
print(summary(img_df$value))

# Applico soglia
img_df_thresh <- img_df %>%
  filter(value < threshold_value)

cat(sprintf("Pixel totali: %d\n", nrow(img_df)))
cat(sprintf("Pixel dopo soglia (< %.2f): %d\n", threshold_value, nrow(img_df_thresh)))

# ========================
# 3) Cluster dei pixel in base all'intensità, usando kmeans++
# ========================
set.seed(123)
km_intensity <- KMeans_rcpp(
  as.matrix(img_df_thresh$value),
  clusters    = k_cell_types,
  num_init    = 5,
  initializer = 'kmeans++',
  seed        = 123
)

img_df_thresh <- img_df_thresh %>%
  mutate(intensity_cluster = factor(km_intensity$clusters))

# ========================
# 4) Campiono n_cells pixel in proporzione al cluster
# ========================
clust_counts <- table(img_df_thresh$intensity_cluster)
clust_freq   <- clust_counts / sum(clust_counts)
cells_per_cluster <- round(clust_freq * n_cells)

cell_list <- vector("list", length(levels(img_df_thresh$intensity_cluster)))
set.seed(123)
for (i in seq_along(levels(img_df_thresh$intensity_cluster))) {
  clust_name <- levels(img_df_thresh$intensity_cluster)[i]
  n_sub      <- cells_per_cluster[clust_name]
  df_sub     <- img_df_thresh %>% filter(intensity_cluster == clust_name)
  idx_sub    <- sample(seq_len(nrow(df_sub)), min(n_sub, nrow(df_sub)), replace = FALSE)
  cell_list[[i]] <- df_sub[idx_sub, ]
}

cell_df <- bind_rows(cell_list)
cat(sprintf("Celle totali campionate: %d\n", nrow(cell_df)))

# ========================
# 5) Generazione dei profili di espressione
# ========================
# 5a) Medie di espressione per cluster
mean_expression_list <- list()
for (k in seq_len(k_cell_types)) {
  mu <- rep(2, n_genes)  # baseline log(7) ~ 2
  start_idx <- (k - 1) * 10 + 1
  end_idx   <- min(k * 10, n_genes)
  if (start_idx <= end_idx) {
    mu[start_idx:end_idx] <- mu[start_idx:end_idx] + 2
  }
  mean_expression_list[[k]] <- mu
}

# 5b) Matrice di correlazione intra/inter cluster
rho_intra <- 0.6
rho_inter <- 0.1
N <- nrow(cell_df)
cluster_labels <- cell_df$intensity_cluster

same_cluster <- outer(cluster_labels, cluster_labels, FUN = `==`)
cell_cor <- ifelse(same_cluster, rho_intra, rho_inter)
diag(cell_cor) <- 1

# 5c) Correlazione spaziale opzionale
if (use_spatial_correlation) {
  coords   <- cell_df[, c("x", "y")]
  dist_mat <- as.matrix(dist(coords))
  range_param <- 30
  spatial_sim <- exp(-dist_mat / range_param)
  alpha <- 0.5
  beta  <- 0.5
  cell_cor <- alpha * cell_cor + beta * spatial_sim
  diag(cell_cor) <- 1
}

# 5d) Generazione espressione (vettoriale)
expression_data <- matrix(0, nrow = N, ncol = n_genes)
for (g in seq_len(n_genes)) {
  base_vec <- mvrnorm(n = 1, mu = rep(0, N), Sigma = cell_cor)
  cluster_means <- sapply(cluster_labels, function(cl) mean_expression_list[[cl]][g])
  expression_data[, g] <- base_vec + cluster_means
}

# Converto in conteggi e applico dropout
expression_data <- exp(expression_data)
expression_data <- matrix(
  rpois(length(expression_data), lambda = expression_data),
  nrow = N, ncol = n_genes
)

dropout_rate <- 0.3
num_zero     <- round(length(expression_data) * dropout_rate)
set.seed(123)
zero_idx <- sample.int(length(expression_data), num_zero)
expression_data[zero_idx] <- 0

result <- list(
  coordinates       = cell_df[, c("x", "y")],
  intensity_cluster = cell_df$intensity_cluster,
  expression        = expression_data,
  threshold_used    = threshold_value
)

# ========================
# 6) Visualizzazione
# ========================
p <- ggplot(cell_df, aes(x = x, y = y, color = intensity_cluster)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  labs(title = sprintf("Distribuzione spaziale (threshold=%.2f)", threshold_value),
       color = "Cell Type")

print(p)

ggplot2::ggsave(here::here("R/daniele/output.png"), plot = p, device = "png", dpi = 300)

# Rinomino i cluster
levels(result$intensity_cluster) <- paste0("cells_", letters[1:k_cell_types])

cat("Tabelle finali di cluster:\n")
print(table(result$intensity_cluster))

saveRDS(result, file = "data/simulated_image_correlation.rds")

tictoc::toc()
