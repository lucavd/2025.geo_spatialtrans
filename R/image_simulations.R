#!/usr/bin/env Rscript

# Script: simulazione di trascrittomica spaziale parametrizzata
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

simulate_spatial_transcriptomics <- function(
  # Parametri generali
  image_path = NULL,            # Path dell'immagine
  output_path = "data/simulated_image_correlation.rds", # Path per salvare i risultati
  output_plot = NULL,           # Path per salvare il plot (opzionale)
  n_cells = 20000,              # Numero di celle da campionare
  n_genes = 100,                # Numero di geni totali
  k_cell_types = 5,             # Numero di tipi cellulari
  threshold_value = 0.7,        # Soglia per il thresholding dell'immagine
  random_seed = 123,            # Seed per riproducibilità
  
  # Parametri per la simulazione dell'espressione genica
  difficulty_level = "hard",    # "easy", "medium", o "hard"
  use_spatial_correlation = TRUE, # Usare correlazione spaziale?
  
  # Parametri avanzati (personalizzabili)
  marker_params = list(
    marker_genes_per_type = NULL,  # Numero di geni marker per tipo cellulare
    marker_expression_fold = NULL, # Intensità dei geni marker (incremento log-fold)
    marker_overlap_fold = NULL     # Intensità dell'overlapping tra tipi adiacenti
  ),
  
  spatial_params = list(
    spatial_noise_intensity = NULL, # Intensità del rumore spaziale
    spatial_range = NULL,           # Range di correlazione spaziale
    random_noise_sd = NULL          # Deviazione standard del rumore casuale
  ),
  
  dropout_params = list(
    dropout_range = NULL,           # Range di probabilità di dropout
    dispersion_range = NULL         # Range del parametro di dispersione
  ),
  
  hybrid_params = list(
    use_hybrid_cells = TRUE,        # Creare cellule ibride ai confini?
    max_hybrid_pairs = 1000,        # Numero massimo di coppie ibride
    hybrid_intensity_range = c(0.2, 0.5) # Range di intensità dell'effetto ibrido
  ),
  
  cell_specific_params = list(
    cell_specific_noise_sd = NULL   # Deviazione standard del rumore cellula-specifico
  )
) {
  tictoc::tic()
  set.seed(random_seed)
  
  # Imposta i parametri in base al livello di difficoltà se non specificati direttamente
  if (is.null(marker_params$marker_genes_per_type)) {
    marker_params$marker_genes_per_type <- switch(
      difficulty_level,
      "easy" = 10,    # Facile: 10 geni marker per tipo
      "medium" = 7,   # Medio: 7 geni marker per tipo
      "hard" = 5      # Difficile: 5 geni marker per tipo
    )
  }
  
  if (is.null(marker_params$marker_expression_fold)) {
    marker_params$marker_expression_fold <- switch(
      difficulty_level,
      "easy" = 2.0,    # Facile: forte differenziale (+2.0 log fold)
      "medium" = 1.2,  # Medio: differenziale moderato (+1.2 log fold)
      "hard" = 0.8     # Difficile: differenziale basso (+0.8 log fold)
    )
  }
  
  if (is.null(marker_params$marker_overlap_fold)) {
    marker_params$marker_overlap_fold <- switch(
      difficulty_level,
      "easy" = 0.0,    # Facile: nessun overlapping
      "medium" = 0.2,  # Medio: overlapping basso (+0.2 log fold)
      "hard" = 0.4     # Difficile: overlapping alto (+0.4 log fold)
    )
  }
  
  if (is.null(spatial_params$spatial_noise_intensity)) {
    spatial_params$spatial_noise_intensity <- switch(
      difficulty_level,
      "easy" = 0.5,    # Facile: rumore spaziale basso
      "medium" = 1.0,  # Medio: rumore spaziale moderato
      "hard" = 1.5     # Difficile: rumore spaziale alto
    )
  }
  
  if (is.null(spatial_params$spatial_range)) {
    spatial_params$spatial_range <- switch(
      difficulty_level,
      "easy" = 50,     # Facile: correlazione a lungo range (più smooth)
      "medium" = 30,   # Medio: correlazione a medio range
      "hard" = 15      # Difficile: correlazione a corto range (più irregolare)
    )
  }
  
  if (is.null(spatial_params$random_noise_sd)) {
    spatial_params$random_noise_sd <- switch(
      difficulty_level,
      "easy" = 0.1,    # Facile: poco rumore casuale
      "medium" = 0.2,  # Medio: rumore casuale moderato
      "hard" = 0.4     # Difficile: molto rumore casuale
    )
  }
  
  if (is.null(dropout_params$dropout_range)) {
    dropout_params$dropout_range <- switch(
      difficulty_level,
      "easy" = c(0.1, 0.3),    # Facile: poco dropout
      "medium" = c(0.2, 0.5),  # Medio: dropout moderato
      "hard" = c(0.3, 0.7)     # Difficile: molto dropout
    )
  }
  
  if (is.null(dropout_params$dispersion_range)) {
    dropout_params$dispersion_range <- switch(
      difficulty_level,
      "easy" = c(3.0, 1.5),    # Facile: bassa dispersione
      "medium" = c(2.0, 1.0),  # Medio: dispersione moderata
      "hard" = c(1.5, 0.8)     # Difficile: alta dispersione
    )
  }
  
  if (is.null(cell_specific_params$cell_specific_noise_sd)) {
    cell_specific_params$cell_specific_noise_sd <- switch(
      difficulty_level,
      "easy" = 0.1,    # Facile: poca variabilità tra cellule
      "medium" = 0.2,  # Medio: variabilità moderata
      "hard" = 0.3     # Difficile: alta variabilità
    )
  }
  
  # Verifica che l'immagine sia specificata
  if (is.null(image_path)) {
    stop("È necessario specificare il percorso dell'immagine (image_path)")
  }
  
  # ========================
  # 2) Caricamento immagine e conversione in data frame sogliato
  # ========================
  cat(sprintf("Carico immagine: %s\n", image_path))
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
  set.seed(random_seed)
  km_intensity <- KMeans_rcpp(
    as.matrix(img_df_thresh$value),
    clusters    = k_cell_types,
    num_init    = 5,
    initializer = 'kmeans++',
    seed        = random_seed
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
  set.seed(random_seed)
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
  # 5a) Medie di espressione per cluster - parametrizzate
  mean_expression_list <- list()
  for (k in seq_len(k_cell_types)) {
    mu <- rep(2, n_genes)  # baseline log(7) ~ 2
    
    # Applicazione dei marker specifici con parametri personalizzati
    start_idx <- (k - 1) * marker_params$marker_genes_per_type + 1
    end_idx   <- min(k * marker_params$marker_genes_per_type, n_genes)
    if (start_idx <= end_idx) {
      mu[start_idx:end_idx] <- mu[start_idx:end_idx] + marker_params$marker_expression_fold
    }
    
    # Aggiungi espressione parziale nei cluster adiacenti (overlapping)
    if (k > 1 && marker_params$marker_overlap_fold > 0) {
      prev_markers <- ((k-2) * marker_params$marker_genes_per_type + 1):min((k-1) * marker_params$marker_genes_per_type, n_genes)
      if (length(prev_markers) > 0) {
        mu[prev_markers] <- mu[prev_markers] + marker_params$marker_overlap_fold
      }
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
  
  # 5c) Impostazione parametri di dispersione in base ai parametri
  dispersion_param <- rescale(mean_dist, to = dropout_params$dispersion_range)
  
  # 5d) Impostazione probabilità di dropout in base ai parametri
  dropout_prob <- rescale(mean_dist, to = dropout_params$dropout_range)
  
  # 5e) Generazione espressione usando Negative Binomial
  expression_data <- matrix(0, nrow = N, ncol = n_genes)
  
  # Identificazione geni stabili (sub-Poissoniani)
  stable_genes <- sample(n_genes, max(1, round(n_genes * 0.1)))  # 10% geni stabili
  
  # Simulazione con Gaussian Process per correlazione spaziale continua
  if (use_spatial_correlation) {
    # Converti cell_df in oggetto spatial
    sp_df <- cell_df
    coordinates(sp_df) <- ~ x + y
    
    # Crea un modello gstat per il GP con parametri personalizzati
    gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                    beta = 0, model = vgm(psill=spatial_params$spatial_noise_intensity, 
                                         range=spatial_params$spatial_range, 
                                         model="Exp"), 
                    nmax=10)
    
    # Genera il noise spaziale
    set.seed(random_seed)
    gp_noise <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
    
    # Normalizza il noise
    gp_noise <- scale(gp_noise)
  }
  
  # Crea cellule ibride ai confini tra cluster, se richiesto
  hybrid_matrix <- matrix(0, nrow = N, ncol = k_cell_types)
  
  if (hybrid_params$use_hybrid_cells) {
    # Identifica coppie di cellule vicine appartenenti a cluster diversi
    hybrid_pairs <- list()
    for (i in 1:N) {
      # Trova cellule vicine (tra le 20 più vicine)
      neighbors <- order(dist_mat[i,])[2:20]
      diff_cluster_neighbors <- neighbors[cluster_labels[neighbors] != cluster_labels[i]]
      
      # Se ci sono vicini di cluster diversi, aggiungi alla lista
      if (length(diff_cluster_neighbors) > 0) {
        hybrid_pairs[[length(hybrid_pairs) + 1]] <- c(i, diff_cluster_neighbors[1])
      }
    }
    
    # Limita a max_hybrid_pairs coppie casuali per efficienza
    if (length(hybrid_pairs) > hybrid_params$max_hybrid_pairs) {
      set.seed(random_seed + 1)
      hybrid_pairs <- hybrid_pairs[sample(length(hybrid_pairs), hybrid_params$max_hybrid_pairs)]
    }
    
    # Crea una matrice di ibridazione (quanto ogni cellula è ibrida di un altro cluster)
    for (pair in hybrid_pairs) {
      cell1 <- pair[1]
      cell2 <- pair[2]
      
      # Prendi i cluster delle due cellule
      cluster1 <- as.integer(cluster_labels[cell1])
      cluster2 <- as.integer(cluster_labels[cell2])
      
      # La cellula 1 è in parte del cluster 2
      hybrid_matrix[cell1, cluster2] <- runif(1, 
                                             hybrid_params$hybrid_intensity_range[1], 
                                             hybrid_params$hybrid_intensity_range[2])
      
      # La cellula 2 è in parte del cluster 1
      hybrid_matrix[cell2, cluster1] <- runif(1, 
                                             hybrid_params$hybrid_intensity_range[1], 
                                             hybrid_params$hybrid_intensity_range[2])
    }
  }
  
  # Genera l'espressione genica con variabilità aggiuntiva tra cellule dello stesso tipo
  for (g in seq_len(n_genes)) {
    cl <- as.integer(cluster_labels)
    
    # Aggiungi variabilità cellula-specifica indipendente dal cluster
    cell_specific_effect <- rnorm(N, 0, cell_specific_params$cell_specific_noise_sd)
    
    # Calcola medie di espressione di base
    base_expr <- sapply(cl, function(x) mean_expression_list[[x]][g])
    
    # Applica effetto di ibridazione tra cluster (cellule ai confini)
    if (hybrid_params$use_hybrid_cells) {
      for (k in 1:k_cell_types) {
        # Trova cellule che hanno componente del cluster k
        hybrid_cells <- which(hybrid_matrix[, k] > 0)
        if (length(hybrid_cells) > 0) {
          # Applica l'effetto ibrido usando il profilo del cluster k
          hybrid_effect <- mean_expression_list[[k]][g] * hybrid_matrix[hybrid_cells, k]
          base_expr[hybrid_cells] <- base_expr[hybrid_cells] * (1 - hybrid_matrix[hybrid_cells, k]) + hybrid_effect
        }
      }
    }
    
    # Combina con l'effetto cellula-specifico
    mu_vals <- base_expr + cell_specific_effect
    
    # Aggiungi correlazione spaziale se richiesta
    if (use_spatial_correlation) {
      mu_vals <- mu_vals + spatial_params$spatial_noise_intensity * gp_noise
      
      # Aggiungi noise casuale addizionale per confondere i pattern
      random_noise <- rnorm(length(mu_vals), 0, spatial_params$random_noise_sd)
      mu_vals <- mu_vals + random_noise
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
  
  # Prepara l'output
  result <- list(
    coordinates       = cell_df[, c("x", "y")],
    intensity_cluster = cell_df$intensity_cluster,
    expression        = expression_data,
    threshold_used    = threshold_value,
    parameters        = list(
      image_path = image_path,
      n_cells = n_cells,
      n_genes = n_genes,
      k_cell_types = k_cell_types,
      difficulty_level = difficulty_level,
      use_spatial_correlation = use_spatial_correlation,
      marker_params = marker_params,
      spatial_params = spatial_params,
      dropout_params = dropout_params,
      hybrid_params = hybrid_params,
      cell_specific_params = cell_specific_params
    )
  )
  
  # ========================
  # 6) Visualizzazione
  # ========================
  p <- ggplot(cell_df, aes(x = x, y = y, color = intensity_cluster)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(title = sprintf("Distribuzione spaziale (livello difficoltà: %s)", difficulty_level),
         subtitle = sprintf("Threshold: %.2f, Geni marker: %d per tipo, Fold-change: %.1f", 
                           threshold_value, 
                           marker_params$marker_genes_per_type,
                           marker_params$marker_expression_fold),
         color = "Cell Type")
  
  print(p)
  
  # Salva il plot se richiesto
  if (!is.null(output_plot)) {
    ggplot2::ggsave(output_plot, plot = p, device = "png", dpi = 300)
  }
  
  # Rinomino i cluster
  levels(result$intensity_cluster) <- paste0("cells_", letters[1:k_cell_types])
  
  cat("Tabelle finali di cluster:\n")
  print(table(result$intensity_cluster))
  
  # Salva il risultato
  saveRDS(result, file = output_path)
  
  tictoc::toc()
  
  return(result)
}

# ========================
# Esempi d'uso
# ========================

# Esempio 1: dataset "facile" (ben separato, con marker forti)
# simulate_spatial_transcriptomics(
#   image_path = here::here("images/colon.png"),
#   output_path = "data/simulated_easy_correlation.rds",
#   output_plot = "results/simulated_easy_distribution.png",
#   difficulty_level = "easy",
#   n_cells = 20000,
#   n_genes = 100
# )

# Esempio 2: dataset "difficile" (confini sfumati, marker deboli)
# simulate_spatial_transcriptomics(
#   image_path = here::here("images/granuloma.png"),
#   output_path = "data/simulated_hard_correlation.rds",
#   output_plot = "results/simulated_hard_distribution.png",
#   difficulty_level = "hard",
#   n_cells = 30000,
#   n_genes = 100
# )

# Esempio 3: dataset con parametri personalizzati
# simulate_spatial_transcriptomics(
#   image_path = here::here("images/colon.png"),
#   output_path = "data/simulated_custom_correlation.rds",
#   output_plot = "results/simulated_custom_distribution.png",
#   difficulty_level = "medium",  # base di partenza
#   n_cells = 25000,
#   n_genes = 150,
#   k_cell_types = 7,
#   marker_params = list(
#     marker_genes_per_type = 8,
#     marker_expression_fold = 1.0,
#     marker_overlap_fold = 0.3
#   ),
#   spatial_params = list(
#     spatial_noise_intensity = 1.2,
#     spatial_range = 25,
#     random_noise_sd = 0.3
#   ),
#   dropout_params = list(
#     dropout_range = c(0.25, 0.6),
#     dispersion_range = c(1.8, 0.9)
#   )
# )