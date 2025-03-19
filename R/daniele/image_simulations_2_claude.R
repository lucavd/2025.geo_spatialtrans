#!/usr/bin/env Rscript

# Script: simulazione di trascrittomica spaziale con miglioramenti
# 1) Caricamento immagine e threshold
# 2) Clustering con kmeans++
# 3) Sampling celle
# 4) Generazione dei profili con modello Negative Binomial
# 5) Simulazione di dropout spaziale
# 6) Visualizzazione

library(imager)
library(tidyverse)
library(MASS)
library(cluster)
library(ClusterR)
library(future)
library(future.apply)
library(tictoc)
library(gstat)
library(sp)
library(scales)

tictoc::tic()

# =======================
# 1) Parametri principali
# =======================
image_path            <- here::here("R/daniele/granuloma.png")
n_cells               <- 3000
n_genes               <- 10
k_cell_types          <- 10
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

# 5b) Calcolo distanze e densità locale per il modello di dropout
N <- nrow(cell_df)
cluster_labels <- cell_df$intensity_cluster
coords <- cell_df %>% dplyr::select(x, y)
dist_mat <- as.matrix(dist(coords))

# Calcola la distanza media di ciascuna cellula rispetto alle altre del proprio cluster
mean_dist <- numeric(N)
for (i in 1:N) {
  mean_dist[i] <- mean(dist_mat[i, cluster_labels == cluster_labels[i]])
}

# Calcola la densità locale (per il modello di dropout)
local_density <- apply(dist_mat, 1, function(row) mean(row < quantile(row, 0.1)))

# 5c) Dispersione crescente ai bordi dei cluster (minore al centro)
dispersion_param <- rescale(mean_dist, to = c(5, 0.5))  # Valori: centro (alto) -> bordi (basso)

# 5d) Dropout probability crescente ai bordi
dropout_prob <- rescale(mean_dist, to = c(0.1, 0.5))  # Maggior dropout ai bordi

# 5e) Generazione espressione usando Negative Binomial
expression_data <- matrix(0, nrow = N, ncol = n_genes)

# Identificazione geni stabili (sub-Poissoniani)
stable_genes <- sample(n_genes, max(1, round(n_genes * 0.1)))  # 10% geni stabili

# Simulazione con Gaussian Process per correlazione spaziale continua
if (use_spatial_correlation) {
  # Converti cell_df in oggetto spatial
  sp_df <- cell_df
  coordinates(sp_df) <- ~ x + y

  # Crea un modello gstat per il GP
  gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                  beta = 0, model = vgm(psill=1, range=30, model="Exp"), nmax=20)

  # Genera il noise spaziale
  set.seed(123)
  gp_noise <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1

  # Normalizza il noise
  gp_noise <- scale(gp_noise)
}

# Genera l'espressione genica
for (g in seq_len(n_genes)) {
  cl <- as.integer(cluster_labels)
  mu_vals <- sapply(cl, function(x) mean_expression_list[[x]][g])

  # Aggiungi correlazione spaziale se richiesta
  if (use_spatial_correlation) {
    mu_vals <- mu_vals + 0.5 * gp_noise
  }

  if (g %in% stable_genes) {
    # Modello sub-Poisson: Binomiale con p alto e n moderato
    p <- 0.9
    n_trial <- round(exp(mu_vals)/(1-p))
    expression_data[, g] <- rbinom(N, n_trial, p)
  } else {
    # Negative Binomial con dispersione variabile spazialmente
    expression_data[, g] <- rnbinom(N, mu = exp(mu_vals), size = dispersion_param)
  }

  # Applica dropout spazialmente variabile
  zero_idx <- rbinom(N, 1, dropout_prob) == 1
  expression_data[zero_idx, g] <- 0
}

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
  geom_point(size = 0.1, alpha = 0.7) +
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
