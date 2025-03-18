#!/usr/bin/env Rscript

# Script: simulazione di trascrittomica spaziale
#  dove le cellule sono raggruppate in base a:
#   - intensità di grigio (cluster "cell type")
#   - distanza spaziale (maggiore correlazione se vicine)
#   - soglia per escludere i pixel troppo chiari/scuri.

# ========================================
# Caricamento librerie
# ========================================
library(imager)    # Per caricamento e gestione immagine
library(tidyverse) # Per funzioni d'aiuto (dplyr, ggplot2, ecc.)
library(MASS)      # Per mvrnorm
library(cluster)   # Per silhouette/clustering o altri metodi
library(parallel)  # Per parallelizzazione di base con mclapply

# Imposta numero di core da utilizzare
num_cores <- max(1, detectCores() - 10)
cat(sprintf("Utilizzo %d core per il calcolo parallelo...\n", num_cores))

# ========================================
# 1) Parametri principali
# ========================================
image_path            <- here::here("R/daniele/colon.png") # Path dell'immagine
n_cells               <- 30000      # Numero di celle da simulare
n_genes               <- 10           # Numero di geni
k_cell_types          <- 3             # Numero di "tipi" cellulari in base all'intensità
use_spatial_correlation <- TRUE        # Se vuoi usare anche la distanza (x,y)

# Parametro di soglia per escludere pixel troppo chiari (o scuri).
# Ad esempio, threshold_value = 0.9 esclude pixel con intensità >= 0.9
threshold_value       <- 0.7

# ========================================
# 2) Caricamento e conversione immagine
# ========================================
img <- load.image(image_path)

# Se c'è un canale alpha (RGBA=4 canali), lo rimuovo
if (spectrum(img) == 4) {
  img <- rm.alpha(img)
}

# Se è a 3 canali (RGB), converto in grigio
if (spectrum(img) == 3) {
  img <- grayscale(img)
}

# Schiaccia dimensioni unitarie (se l'immagine è [W,H,1,1])
img <- squeeze(img)

# A questo punto, img potrebbe essere un oggetto "cimg" solo 2D (W,H).
# Verifico dimensioni
cat("Dimensioni dell'immagine cimg:\n")
print(dim(img))

# ========================================
# 3) Estraggo intensità e coordinate di tutti i pixel
# ========================================
# Se dim(img) = c(W, H), 'imager::as.data.frame.cimg()' si aspetta 4 dimensioni.
# Convertiamo manualmente in array 2D e poi in cimg a 4 dimensioni (x,y,cc,t).

img_array <- as.array(img)
# Creo un cimg con dims=(width, height, 1, 1).
# NOTA: a volte width e height possono essere invertiti, dipende dall'ordine interno.
# Per sicurezza, puoi scambiare i valori se noti un'inversione.
w <- nrow(img_array)     # es. 350
h <- ncol(img_array)     # es. 271
img_cimg <- as.cimg(img_array, dims = c(w, h, 1, 1))

dim(img_cimg)
spectrum(img_cimg)

# Converto in data frame: x, y, cc, value
img_df <- as.data.frame(img_cimg)

# Seleziono solo le colonne di interesse (x,y,value)
img_df <- img_df %>%
  dplyr::select(x, y, value)

# Controllo range intensità
cat("Range intensità immagine:\n")
print(summary(img_df$value))

# ========================================
# 4) Applico la soglia per escludere pixel indesiderati
# ========================================
# Esempio: voglio escludere i pixel con intensità >= threshold_value
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
# 6) Campiono n_cells pixel
#    proporzionale alla dimensione di ogni cluster
# ========================================
clust_counts <- table(img_df_thresh$intensity_cluster)
clust_freq   <- clust_counts / sum(clust_counts)
cells_per_cluster <- round(clust_freq * n_cells)

cell_df <- data.frame()
set.seed(123)
for (clust in levels(img_df_thresh$intensity_cluster)) {
  n_sub  <- cells_per_cluster[clust]
  df_sub <- img_df_thresh %>% filter(intensity_cluster == clust)

  # Campiono n_sub righe
  # Attenzione: se n_sub > nrow(df_sub), sostituisci con un min(...) o rimani con un warning
  idx_sub <- sample(seq_len(nrow(df_sub)), n_sub, replace = FALSE)
  sampled_sub <- df_sub[idx_sub, ]

  cell_df <- rbind(cell_df, sampled_sub)
}

cat(sprintf("Celle totali campionate: %d\n", nrow(cell_df)))

# ========================================
# 7) Generazione dei profili di espressione
# ========================================
# 7.1) Definisco un "profilo medio" per ciascun cluster
mean_expression_list <- list()
for (k in 1:k_cell_types) {
  # baseline ~ log(7) => 2
  mu <- rep(2, n_genes)

  start_idx <- (k-1)*10 + 1
  end_idx   <- min(k*10, n_genes)
  if (start_idx <= end_idx) {
    # up-reg di 2 => e^2 ~ 7x fold change
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

# Metodo OTTIMIZZATO per costruire la matrice di correlazione
# Evita il doppio ciclo for usando matrici e operazioni vettoriali
cat("Usando metodo ottimizzato per la costruzione della matrice di correlazione...\n")

# Creo vettori di indici cluster
cluster_indices <- as.integer(cluster_labels)

# Creo una matrice di confronto (TRUE quando i cluster sono uguali)
# Questo è molto più veloce del doppio ciclo for
equal_clusters <- outer(cluster_indices, cluster_indices, "==")

# Assegno i valori di correlazione in base al confronto dei cluster
cell_cor <- matrix(0, nrow=N, ncol=N)
cell_cor[equal_clusters] <- rho_intra
cell_cor[!equal_clusters] <- rho_inter
diag(cell_cor) <- 1

cat("Matrice di correlazione iniziale costruita.\n")

# 7.3) Opzionale: aggiungo contributo dalla distanza spaziale
if (use_spatial_correlation) {
  cat("Inizio calcolo matrice di distanza...\n")
  coords <- cell_df[, c("x", "y")]

  # Calcolo distanza con metodo vettoriale ottimizzato se possibile
  tryCatch({
    # Metodo vettoriale (utilizza la funzione dist() che è ottimizzata in C)
    dist_mat <- as.matrix(dist(coords))
    cat("Matrice di distanza calcolata con metodo vettoriale.\n")
  }, error = function(e) {
    cat("Errore nel calcolo vettoriale, uso metodo parallelo...\n")

    # Dividi il lavoro in blocchi
    chunk_size <- ceiling(N / min(100, num_cores * 4))  # Dimensione chunk
    row_indices <- 1:N
    row_chunks <- split(row_indices, ceiling(seq_along(row_indices) / chunk_size))

    cat(sprintf("Diviso il lavoro in %d blocchi\n", length(row_chunks)))

    # Funzione per calcolare le distanze per un blocco di righe
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

    # Prova a usare mclapply (parallelizzazione su UNIX/Linux)
    dist_chunks <- mclapply(row_chunks, calc_chunk, mc.cores=num_cores)

    # Combina i risultati
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

# 7.4) Genero la matrice di espressione
cat("Generazione dei profili di espressione...\n")
expression_data <- matrix(0, nrow=N, ncol=n_genes)

for (g in seq_len(n_genes)) {
  # Vettore base con correlazione cell_cor
  base_vec <- mvrnorm(n=1, mu=rep(0,N), Sigma=cell_cor)

  # Aggiungo la media di cluster
  gene_values <- numeric(N)
  for (i in seq_len(N)) {
    cl <- as.integer(cluster_labels[i])
    mu_val <- mean_expression_list[[cl]][g]
    gene_values[i] <- base_vec[i] + mu_val
  }

  expression_data[,g] <- gene_values
}

# Converto in conteggi
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

cat("Generazione dei profili completata.\n")

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
# Visualizzo la distribuzione spaziale colorando per "intensity_cluster"
p <- ggplot(cell_df, aes(x=x, y=y, color=intensity_cluster)) +
  geom_point(size=2, alpha=0.7) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  labs(title = sprintf("Distribuzione spaziale (threshold=%.2f)", threshold_value),
       color = "Cell Type")

# p

ggplot2::ggsave(here::here("R/daniele/output.png"), plot = p, device = "png", dpi = 600 )

# (Opzionale) Rinomino i cluster con "cells_a", "cells_b", ...
levels(result$intensity_cluster) <- paste0("cells_", letters[1:k_cell_types])

cat("Tabelle finali di cluster:\n")
table(result$intensity_cluster)

saveRDS(result, file = "data/simulated_image_correlation.rds")

cat("Script completato con successo.\n")
# Fine script
