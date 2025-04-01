#!/usr/bin/env Rscript

# Script: simulazione di trascrittomica spaziale parametrizzata per Visium HD
# 1) Caricamento immagine e threshold
# 2) Creazione di una griglia regolare ad alta risoluzione (tipo Visium HD)
# 3) Clustering con kmeans++
# 4) Generazione dei profili con modello Negative Binomial
# 5) Simulazione di dropout spaziale e correlazione
# 6) Visualizzazione

# Aumenta il limite massimo di dimensione dei dati per la parallelizzazione
options(future.globals.maxSize = 10000 * 1024^2)  # 10GB

# Imposta il numero di core disponibili e disabilita i limiti di parallelizzazione 
options(mc.cores = 110)
options(future.fork.enable = TRUE)
options(future.rng.onMisuse = "ignore")

# Carica il pacchetto parallelly se disponibile e imposta opzioni
if (requireNamespace("parallelly", quietly = TRUE)) {
  try({
    options(parallelly.availableCores.fallback = 110)
    options(parallelly.validateCores.limit = 999)
  }, silent = TRUE)
}

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
library(spatstat)    # Per funzioni spaziali avanzate
library(spdep)       # Per analisi di autocorrelazione spaziale
library(fields)      # Per simulazione GRF (Gaussian Random Fields)

simulate_spatial_transcriptomics <- function(
  # Parametri generali
  image_path = NULL,            # Path dell'immagine
  output_path = "data/simulated_image_correlation.rds", # Path per salvare i risultati
  output_plot = NULL,           # Path per salvare il plot (opzionale)
  n_cells = 20000,              # Numero di celle da campionare (se grid_mode=FALSE)
  n_genes = 100,                # Numero di geni totali
  k_cell_types = 5,             # Numero di tipi cellulari
  threshold_value = 0.7,        # Soglia per il thresholding dell'immagine
  random_seed = 123,            # Seed per riproducibilità
  pixel_size_um = 1,            # Dimensione di ogni pixel dell'immagine in μm (valore predefinito: 1μm per pixel)

  # Parametri della griglia Visium HD
  grid_mode = TRUE,             # Usare modalità griglia (Visium HD) invece di sampling casuale
  grid_resolution = 2,          # Dimensione della griglia in μm (es. 2μm x 2μm per Visium HD)
  grid_spacing = 0,             # Spazio tra bin della griglia (0 per griglia continua)
  use_fixed_grid = FALSE,       # Usare una griglia fissa con dimensioni predefinte (Visium HD-like)
  fixed_grid_width_mm = 6.5,    # Larghezza della griglia fissa in mm (standard Visium HD)
  fixed_grid_height_mm = 6.5,   # Altezza della griglia fissa in mm (standard Visium HD)

  # Parametri per la simulazione dell'espressione genica
  difficulty_level = "hard",    # "easy", "medium", o "hard"
  use_spatial_correlation = TRUE, # Usare correlazione spaziale?
  correlation_method = "grf",   # Metodo di correlazione spaziale: "grf" (Gaussian Random Field) o "car" (Conditional Autoregressive)

  # Parametri avanzati (personalizzabili)
  marker_params = list(
    marker_genes_per_type = NULL,  # Numero di geni marker per tipo cellulare
    marker_expression_fold = NULL, # Intensità dei geni marker (incremento log-fold)
    marker_overlap_fold = NULL     # Intensità dell'overlapping tra tipi adiacenti
  ),

  spatial_params = list(
    spatial_noise_intensity = NULL, # Intensità del rumore spaziale
    spatial_range = NULL,           # Range di correlazione spaziale
    random_noise_sd = NULL,         # Deviazione standard del rumore casuale
    gradient_regions = FALSE,       # Usare gradienti invece di confini netti tra regioni
    gradient_width = 5,             # Ampiezza del gradiente (in unità della griglia)
    gradient_exponent = 1.5         # Esponente per la transizione del gradiente (1=lineare, >1=più graduale al centro)
  ),

  dropout_params = list(
    dropout_range = NULL,           # Range di probabilità di dropout
    dispersion_range = NULL,        # Range del parametro di dispersione
    cell_type_dispersion_effect = 0.2, # Effetto del tipo cellulare sulla dispersione
    expression_dependent_dropout = TRUE, # Dropout dipendente dal livello di espressione
    dropout_curve_midpoint = 0.5,   # Punto medio della curva di dropout (in scala di espressione normalizzata)
    dropout_curve_steepness = 5     # Ripidità della curva di dropout
  ),

  library_size_params = list(
    mean_library_size = 10000,      # Dimensione media della libreria per spot
    library_size_cv = 0.3,          # Coefficiente di variazione della dimensione della libreria
    spatial_effect_on_library = 0.5, # Effetto spaziale sulla dimensione della libreria (0-1)
    cell_type_effect = TRUE         # Effetto del tipo cellulare sulla dimensione libreria
  ),

  hybrid_params = list(
    use_hybrid_cells = TRUE,        # Creare cellule ibride ai confini?
    max_hybrid_pairs = 1000,        # Numero massimo di coppie ibride
    hybrid_intensity_range = c(0.2, 0.5) # Range di intensità dell'effetto ibrido
  ),

  cell_specific_params = list(
    cell_specific_noise_sd = NULL,  # Deviazione standard del rumore cellula-specifico
    use_gene_modules = TRUE,        # Usare moduli di geni co-espressi
    n_gene_modules = 5,             # Numero di moduli genici
    module_correlation = 0.7        # Correlazione tra geni nello stesso modulo
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
  # 4a) Creazione di griglia Visium HD o sampling casuale
  # ========================
  if (grid_mode) {
    cat("Modalità griglia attiva - Creazione griglia Visium HD-like\n")

    # Calcola le dimensioni dell'immagine in μm
    img_width_um <- w * pixel_size_um
    img_height_um <- h * pixel_size_um
    
    # Determina le dimensioni della griglia e le coordinate
    if (use_fixed_grid) {
      # Usa una griglia fissa con dimensioni standard Visium HD (6.5mm x 6.5mm)
      # Converti da mm a μm
      grid_width_um <- fixed_grid_width_mm * 1000  # 6.5mm = 6500μm
      grid_height_um <- fixed_grid_height_mm * 1000  # 6.5mm = 6500μm
      
      # Calcola il numero di bin necessari per la griglia fissa
      n_bins_x <- ceiling(grid_width_um / grid_resolution)
      n_bins_y <- ceiling(grid_height_um / grid_resolution)
      
      cat(sprintf("Dimensione immagine originale: %d x %d μm\n", img_width_um, img_height_um))
      cat(sprintf("Usando griglia fissa Visium HD: %.1f x %.1f mm (%d x %d μm)\n", 
                  fixed_grid_width_mm, fixed_grid_height_mm, grid_width_um, grid_height_um))
      cat(sprintf("Risoluzione griglia: %d μm, risultato: %d x %d bin = %d bin totali\n",
                  grid_resolution, n_bins_x, n_bins_y, n_bins_x * n_bins_y))
      
      # Crea le coordinate della griglia fissa in μm
      x_coords <- seq(0, grid_width_um - grid_resolution, by = grid_resolution + grid_spacing)
      y_coords <- seq(0, grid_height_um - grid_resolution, by = grid_resolution + grid_spacing)
    } else {
      # Usa una griglia che si adatta all'immagine (comportamento originale)
      grid_width_um <- img_width_um
      grid_height_um <- img_height_um
      
      # Calcola il numero di bin necessari
      n_bins_x <- ceiling(img_width_um / grid_resolution)
      n_bins_y <- ceiling(img_height_um / grid_resolution)
      
      cat(sprintf("Dimensione immagine: %d x %d μm\n", img_width_um, img_height_um))
      cat(sprintf("Risoluzione griglia: %d μm, risultato: %d x %d bin = %d bin totali\n",
                  grid_resolution, n_bins_x, n_bins_y, n_bins_x * n_bins_y))
      
      # Crea le coordinate della griglia adattata all'immagine in μm
      x_coords <- seq(0, img_width_um - grid_resolution, by = grid_resolution + grid_spacing)
      y_coords <- seq(0, img_height_um - grid_resolution, by = grid_resolution + grid_spacing)
    }

    # Crea un dataframe con tutte le coordinate possibili
    grid_points <- expand.grid(x = x_coords, y = y_coords)

    # Per la griglia fissa, dobbiamo assicurarci che i punti corrispondano all'immagine
    if (use_fixed_grid) {
      # Calcola l'offset per centrare l'immagine nella griglia, se necessario
      if (grid_width_um > img_width_um) {
        offset_x <- (grid_width_um - img_width_um) / 2
      } else {
        offset_x <- 0
      }
      
      if (grid_height_um > img_height_um) {
        offset_y <- (grid_height_um - img_height_um) / 2
      } else {
        offset_y <- 0
      }
      
      # Assegna a ciascun punto della griglia il valore dell'immagine, se il punto è all'interno dell'immagine
      grid_df <- grid_points %>%
        rowwise() %>%
        mutate(
          # Calcola le coordinate relative all'immagine, considerando l'offset
          rel_x = x - offset_x,
          rel_y = y - offset_y,
          
          # Controlla se il punto è all'interno dell'immagine
          is_in_image = (rel_x >= 0 && rel_x < img_width_um && rel_y >= 0 && rel_y < img_height_um),
          
          # Se il punto è fuori dall'immagine, assegna un valore superiore alla soglia
          # altrimenti prendi il valore dall'immagine
          value = if (is_in_image) {
            # Converti da coordinate μm a indici di pixel
            img_x = min(max(round(rel_x / pixel_size_um), 1), w)
            img_y = min(max(round(rel_y / pixel_size_um), 1), h)
            img_array[img_x, img_y]
          } else {
            1.0  # Valore superiore alla soglia, sarà filtrato
          }
        ) %>%
        filter(value < threshold_value) %>%  # Applica la soglia
        ungroup()
    } else {
      # Comportamento originale per la griglia adattata all'immagine
      grid_df <- grid_points %>%
        rowwise() %>%
        mutate(
          # Converti da μm a indici di pixel nell'immagine originale
          img_x = min(max(round(x / pixel_size_um), 1), w),
          img_y = min(max(round(y / pixel_size_um), 1), h),
          value = img_array[img_x, img_y]
        ) %>%
        filter(value < threshold_value) %>%  # Applica la stessa soglia
        ungroup()
    }

    # Assegna cluster usando i centroidi k-means dell'immagine originale
    # Trova i centroidi dei cluster
    cluster_centroids <- sapply(levels(img_df_thresh$intensity_cluster), function(cl) {
      mean(img_df_thresh$value[img_df_thresh$intensity_cluster == cl])
    })

    # Assegna ogni punto griglia al cluster più vicino
    grid_df <- grid_df %>%
      mutate(intensity_cluster = factor(apply(outer(value, cluster_centroids, FUN = function(x, y) abs(x - y)),
                                            1, which.min)))

    # Questa è la nostra "cell_df" finale
    cell_df <- grid_df

    # Se vogliamo applicare gradienti ai confini tra regioni
    if (spatial_params$gradient_regions) {
      cat("Applicazione gradienti tra regioni\n")

      # Crea una matrice di distanza per calcolare quanto ogni punto è vicino al confine
      cell_coords <- cell_df %>% dplyr::select(x, y)

      # Per ogni cluster, identifica i punti di bordo
      all_clusters <- levels(cell_df$intensity_cluster)
      boundary_dists <- matrix(Inf, nrow = nrow(cell_df), ncol = length(all_clusters))
      
      # Utilizzare parallelizzazione per i cluster, forza il numero di workers
      future::plan(future::multisession, workers = min(110, 15))
      # Ignora riferimenti per ridurre la memoria necessaria
      options(future.globals.onReference = "ignore")
      
      # Pre-calcola indici x,y di matrice per tutti i punti in una sola volta
      x_idx <- ceiling((cell_df$x - min(cell_df$x)) / grid_resolution) + 1
      y_idx <- ceiling((cell_df$y - min(cell_df$y)) / grid_resolution) + 1
      x_idx <- pmin(pmax(x_idx, 1), n_bins_x)
      y_idx <- pmin(pmax(y_idx, 1), n_bins_y)
      
      # Processa i cluster in parallelo
      cluster_results <- future_lapply(seq_along(all_clusters), function(i) {
        cl <- all_clusters[i]
        # Punti in questo cluster
        in_cluster <- cell_df$intensity_cluster == cl

        # Crea una matrice binaria per il cluster e calcola la distanza dal bordo
        cluster_mat <- matrix(0, nrow = n_bins_x, ncol = n_bins_y)
        
        # Utilizziamo gli indici x,y pre-calcolati

        # Segna i punti in questo cluster - vettorizzato
        idx <- in_cluster
        if (sum(idx) > 0) {
          # Converti indici logici in indici di matrice in una singola operazione
          cluster_mat[cbind(x_idx[idx], y_idx[idx])] <- 1

          # Crea una matrice di bordo vettorizzando la logica del bordo
          # Modo più efficiente per trovare i bordi
          border_mat <- matrix(0, nrow = n_bins_x, ncol = n_bins_y)

          # Applica filtro per rilevare i bordi
          # Prima identifica tutti i punti interni (circondati da 1)
          interior <- matrix(0, nrow = n_bins_x, ncol = n_bins_y)
          interior[2:(n_bins_x-1), 2:(n_bins_y-1)] <-
            cluster_mat[2:(n_bins_x-1), 2:(n_bins_y-1)] *
            cluster_mat[1:(n_bins_x-2), 2:(n_bins_y-1)] *
            cluster_mat[3:n_bins_x, 2:(n_bins_y-1)] *
            cluster_mat[2:(n_bins_x-1), 1:(n_bins_y-2)] *
            cluster_mat[2:(n_bins_x-1), 3:n_bins_y]

          # I bordi sono i punti che sono nel cluster ma non interni
          border_mat <- cluster_mat * (1 - (interior > 0))

          # Aggiungi anche tutti i punti di bordo esterno (i confini della griglia)
          border_mat[1,] <- border_mat[1,] | (cluster_mat[1,] > 0)
          border_mat[n_bins_x,] <- border_mat[n_bins_x,] | (cluster_mat[n_bins_x,] > 0)
          border_mat[,1] <- border_mat[,1] | (cluster_mat[,1] > 0)
          border_mat[,n_bins_y] <- border_mat[,n_bins_y] | (cluster_mat[,n_bins_y] > 0)

          # Trova le coordinate dei punti di bordo
          border_indices <- which(border_mat > 0, arr.ind = TRUE)

          if (nrow(border_indices) > 0) {
            # Converti indici in coordinate reali in μm
            border_coords <- matrix(0, nrow = nrow(border_indices), ncol = 2)
            border_coords[,1] <- min(cell_df$x) + (border_indices[,1] - 1) * grid_resolution
            border_coords[,2] <- min(cell_df$y) + (border_indices[,2] - 1) * grid_resolution

            # Per ogni punto nella griglia, calcola la distanza minima da un punto di bordo
            # in modo vettorizzato
            cell_coords_mat <- as.matrix(cell_coords)

            # Questa operazione calcola tutte le distanze punto-bordo in una volta
            # ed estrae il minimo per ogni punto
            # La funzione pdist calcola la distanza euclidea tra tutte le coppie di punti
            # in due matrici
            # Utilizza l'elaborazione parallela a blocchi per evitare problemi di memoria
            boundary_dists[, i] <- Inf
            block_size <- 100  # Blocchi più grandi per maggiore efficienza parallela
            n_blocks <- ceiling(nrow(border_coords) / block_size)
            
            # Elabora anche i punti della griglia a blocchi per risparmiare memoria
            grid_block_size <- 500  # Aumentato per bilanciare efficienza/memoria
            grid_blocks <- ceiling(nrow(cell_coords_mat) / grid_block_size)
            
            # Configura parallelizzazione per i calcoli di distanza
            # Limitiamo il numero di workers per evitare sovraccarico di memoria
            future::plan(future::multisession, workers = min(110, 20))
            
            # Definiamo le coppie di blocchi da processare in un'unica lista
            block_pairs <- expand.grid(b = 1:n_blocks, gb = 1:grid_blocks)
            
            # Parallelizziamo il calcolo delle distanze
            results <- future_lapply(1:nrow(block_pairs), function(pair_idx) {
              b <- block_pairs$b[pair_idx]
              gb <- block_pairs$gb[pair_idx]
              
              # Calcola gli indici per questo blocco di punti bordo
              start_idx <- (b-1) * block_size + 1
              end_idx <- min(b * block_size, nrow(border_coords))
              block_coords <- border_coords[start_idx:end_idx, , drop = FALSE]
              
              # Calcola gli indici per questo blocco di punti griglia
              grid_start <- (gb-1) * grid_block_size + 1
              grid_end <- min(gb * grid_block_size, nrow(cell_coords_mat))
              
              # Calcola distanze solo per il sottoinsieme corrente
              block_dists <- fields::rdist(cell_coords_mat[grid_start:grid_end, , drop = FALSE], block_coords)
              
              # Calcola i minimi per riga e restituisci insieme agli indici della griglia
              list(
                grid_range = c(grid_start, grid_end),
                min_dists = apply(block_dists, 1, min)
              )
            })
            
            # Unisci i risultati
            for (r in results) {
              grid_start <- r$grid_range[1]
              grid_end <- r$grid_range[2]
              boundary_dists[grid_start:grid_end, i] <- pmin(
                boundary_dists[grid_start:grid_end, i],
                r$min_dists
              )
            }
          }
        }
      })

      # Normalizza le distanze per ogni cluster
      for (i in seq_along(all_clusters)) {
        if (any(is.finite(boundary_dists[, i]))) {
          boundary_dists[, i] <- boundary_dists[, i] / max(boundary_dists[, i], na.rm = TRUE)
        }
      }

      # Salva le distanze dal bordo per uso successivo
      cell_df$boundary_dist <- apply(boundary_dists, 1, min, na.rm = TRUE)

      # Se dist < gradient_width, consideriamo il punto in zona di transizione
      cell_df$in_gradient <- cell_df$boundary_dist < (spatial_params$gradient_width * grid_resolution /
                                                     max(img_width_um, img_height_um))
    }

  } else {
    # Modalità campionamento originale
    cat("Modalità campionamento casuale attiva\n")

    # ========================
    # 4b) Campiono n_cells pixel in proporzione al cluster
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
  }

  cat(sprintf("Celle/bin totali: %d\n", nrow(cell_df)))

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
  # Versione vettorizzata e ottimizzata
  mean_dist <- numeric(N)
  unique_clusters <- unique(cluster_labels)
  
  # future_lapply per parallelizzare
  # Utilizziamo parallelizzazione bilanciata
  future::plan(future::multisession, workers = min(110, 40))
  # Ignora riferimenti per ridurre la memoria necessaria
  options(future.globals.onReference = "ignore")
  
  # Chunk più grandi per migliorare l'efficienza
  chunk_size <- max(1, ceiling(N/1000))
  chunks <- split(1:N, ceiling(seq_along(1:N)/chunk_size))
  
  mean_dist <- future_lapply(chunks, function(chunk_idx) {
    result <- numeric(length(chunk_idx))
    for (j in seq_along(chunk_idx)) {
      i <- chunk_idx[j]
      cl <- cluster_labels[i]
      same_cluster <- which(cluster_labels == cl)
      result[j] <- mean(dist_mat[i, same_cluster])
    }
    return(result)
  }) %>% unlist()
  
  # Calcola la densità locale (per il modello di dropout)
  # Versione ottimizzata con chunking
  local_density <- future_lapply(chunks, function(chunk_idx) {
    result <- numeric(length(chunk_idx))
    for (j in seq_along(chunk_idx)) {
      i <- chunk_idx[j]
      row <- dist_mat[i,]
      q <- quantile(row, 0.1)
      result[j] <- mean(row < q)
    }
    return(result)
  }) %>% unlist()

  # 5c) Impostazione parametri di dispersione in base ai parametri
  if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df)) {
    # Se utilizziamo gradienti, facciamo variare la dispersione in base alla distanza dal confine
    # Più vicino al confine = più variabilità (parametro dispersione più basso)
    # boundary_dist è già normalizzato da 0 a 1
    dispersion_param <- dropout_params$dispersion_range[2] +
      cell_df$boundary_dist * (dropout_params$dispersion_range[1] - dropout_params$dispersion_range[2])
  } else {
    # Altrimenti usiamo il metodo originale basato sulla distanza media
    dispersion_param <- rescale(mean_dist, to = dropout_params$dispersion_range)
  }

  # Aggiungi effetto del tipo cellulare sulla dispersione
  # Alcuni tipi cellulari sono intrinsecamente più variabili di altri
  set.seed(random_seed + 4)
  # Crea effetti casuali per ogni tipo cellulare
  type_dispersion_effects <- runif(k_cell_types,
                                 min = 1 - dropout_params$cell_type_dispersion_effect,
                                 max = 1 + dropout_params$cell_type_dispersion_effect)

  # Applica effetto moltiplicativo per tipo cellulare
  for (k in 1:k_cell_types) {
    dispersion_param[cluster_labels == k] <- dispersion_param[cluster_labels == k] * type_dispersion_effects[k]
  }

  # 5d) Simulazione delle dimensioni delle librerie
  # La dimensione della libreria influisce sui conteggi finali
  set.seed(random_seed + 1)

  # Genera dimensioni libreria con distribuzione log-normale
  library_size_sd <- library_size_params$mean_library_size * library_size_params$library_size_cv
  log_mean <- log(library_size_params$mean_library_size^2 /
                 sqrt(library_size_params$mean_library_size^2 + library_size_sd^2))
  log_sd <- sqrt(log(1 + (library_size_sd^2 / library_size_params$mean_library_size^2)))

  library_size <- rlnorm(N, meanlog = log_mean, sdlog = log_sd)

  # Aggiungi effetto spaziale sulla dimensione libreria se richiesto
  if (library_size_params$spatial_effect_on_library > 0) {
    # Converti cell_df in oggetto spatial per il GP
    sp_df_lib <- cell_df
    coordinates(sp_df_lib) <- ~ x + y

    # Crea un GP per l'effetto spaziale sulla dimensione libreria
    # Usa una correlazione spaziale a raggio più ampio per la dimensione libreria
    lib_gp <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                   beta = 0, model = vgm(psill = 1.0,
                                        range = spatial_params$spatial_range * 1.5,
                                        model = "Exp"),
                   nmax = 15)

    # Genera il noise spaziale
    set.seed(random_seed + 2)
    lib_noise <- predict(lib_gp, newdata = sp_df_lib, nsim = 1)$sim1

    # Normalizza e scala il noise
    lib_noise <- scale(lib_noise)
    lib_effect <- library_size_params$spatial_effect_on_library * lib_noise

    # Applica alla dimensione libreria (effetto moltiplicativo)
    library_size <- library_size * exp(lib_effect)
  }

  # Aggiungi effetto del tipo cellulare sulla dimensione della libreria
  if (library_size_params$cell_type_effect) {
    # Diversi tipi cellulari hanno diversi contenuti di RNA
    cell_type_effect <- numeric(N)

    # Crea effetti diversi per diversi tipi cellulari
    # Alcuni tipi cellulari hanno librerie sistematicamente più grandi
    set.seed(random_seed + 3)
    type_effects <- rnorm(k_cell_types, mean = 0, sd = 0.2)  # Effetti casuali per tipo

    # Assegna effetto in base al tipo cellulare
    for (k in 1:k_cell_types) {
      cell_type_effect[cluster_labels == k] <- type_effects[k]
    }

    # Applica effetto moltiplicativo
    library_size <- library_size * exp(cell_type_effect)
  }

  # 5e) Impostazione probabilità di dropout
  if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df)) {
    # Più dropout vicino al confine
    base_dropout <- dropout_params$dropout_range[1] +
      (1 - cell_df$boundary_dist) * (dropout_params$dropout_range[2] - dropout_params$dropout_range[1])
  } else {
    # Metodo originale
    base_dropout <- rescale(mean_dist, to = dropout_params$dropout_range)
  }

  # 5f) Generazione espressione usando Negative Binomial
  expression_data <- matrix(0, nrow = N, ncol = n_genes)

  # Identificazione geni stabili (sub-Poissoniani)
  stable_genes <- sample(n_genes, max(1, round(n_genes * 0.1)))  # 10% geni stabili

  # Crea moduli di geni co-espressi
  gene_modules <- NULL
  module_noise <- NULL

  if (cell_specific_params$use_gene_modules) {
    # Calcola il numero di geni per modulo
    genes_per_module <- ceiling(n_genes / cell_specific_params$n_gene_modules)

    # Assegna geni ai moduli
    gene_modules <- list()
    for (m in 1:cell_specific_params$n_gene_modules) {
      start_idx <- (m-1) * genes_per_module + 1
      end_idx <- min(m * genes_per_module, n_genes)
      gene_modules[[m]] <- start_idx:end_idx
    }

    # Crea rumore correlato per ogni modulo
    # Questo rumore sarà applicato a tutti i geni dello stesso modulo
    module_noise <- matrix(0, nrow = N, ncol = n_genes)

    # Genera rumore base per ogni modulo
    base_module_noise <- matrix(rnorm(cell_specific_params$n_gene_modules * N),
                               nrow = N, ncol = cell_specific_params$n_gene_modules)

    # Applica il rumore del modulo a ciascun gene appartenente al modulo
    for (m in 1:length(gene_modules)) {
      module_genes <- gene_modules[[m]]
      # Assegna lo stesso rumore base a tutti i geni del modulo, scalato per la correlazione
      module_noise[, module_genes] <- base_module_noise[, m] * cell_specific_params$module_correlation
    }
  }

  # Simulazione della correlazione spaziale continua - versione ottimizzata
  if (use_spatial_correlation) {
    if (correlation_method == "grf") {
      # Metodo GRF (Gaussian Random Field) con gstat
      # Utilizza blocchi più piccoli per ridurre la memoria
      sp_df <- cell_df
      coordinates(sp_df) <- ~ x + y
      
      # Usa parallelizzazione per accelerare il calcolo - limitiamo per evitare problemi di memoria
      future::plan(future::multisession, workers = min(110, 20))
      # Ignora riferimenti per ridurre la memoria necessaria
      options(future.globals.onReference = "ignore")
      
      # Calcola il numero ottimale di blocchi in base al numero di punti
      n_points <- nrow(sp_df)
      # Blocchi più piccoli per efficienza con parallelizzazione massiccia
      block_size <- min(1000, max(200, floor(n_points/20)))
      n_blocks <- ceiling(n_points / block_size)
      
      cat(sprintf("Processando correlazione spaziale in %d blocchi di %d punti\n", 
                 n_blocks, block_size))
      
      # Inizializza il vettore del rumore
      gp_noise <- numeric(n_points)
      
      if (n_blocks > 1) {
        # Processo a blocchi per risparmiare memoria
        for (b in 1:n_blocks) {
          start_idx <- (b-1) * block_size + 1
          end_idx <- min(b * block_size, n_points)
          
          # Crea un modello gstat solo per questo blocco
          sp_block <- sp_df[start_idx:end_idx,]
          
          gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                         beta = 0, model = vgm(psill = spatial_params$spatial_noise_intensity * 2,
                                              range = spatial_params$spatial_range * 0.8,
                                              model = "Exp"),
                         nmax = 20)
                         
          # Genera il noise spaziale per questo blocco
          set.seed(random_seed + b) # Seed diverso per ogni blocco ma riproducibile
          noise_block <- predict(gp_sim, newdata = sp_block, nsim = 1)$sim1
          
          # Assegna al vettore completo
          gp_noise[start_idx:end_idx] <- noise_block
        }
      } else {
        # Metodo originale se abbiamo pochi punti
        gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                       beta = 0, model = vgm(psill = spatial_params$spatial_noise_intensity * 2,
                                            range = spatial_params$spatial_range * 0.8,
                                            model = "Exp"),
                       nmax = 20)
        
        # Genera il noise spaziale
        set.seed(random_seed)
        gp_noise <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
      }
      
      # Normalizza il noise ma mantieni maggiore varianza
      gp_noise <- scale(gp_noise) * 1.5

    } else if (correlation_method == "car") {
      # Metodo CAR (Conditional Autoregressive) con spdep
      # Crea una griglia di vicinanza per i punti
      coords_matrix <- as.matrix(cell_df[, c("x", "y")])

      # Definiamo vicinanza come punti entro una certa distanza
      nb <- dnearneigh(coords_matrix, 0, spatial_params$spatial_range)

      # Se abbiamo troppi pochi vicini, aumentiamo la distanza
      if (min(card(nb)) < 3) {
        cat("Aumentiamo il raggio per includere più vicini\n")
        nb <- dnearneigh(coords_matrix, 0, spatial_params$spatial_range * 2)
      }

      # Creiamo la matrice dei pesi spaziali
      lw <- nb2listw(nb, style = "W")

      # Simula un processo spaziale autocorrelato
      set.seed(random_seed)
      gp_noise <- spam.spGauss(lw, sigma2 = spatial_params$spatial_noise_intensity, n = 1)

      # Normalizza
      gp_noise <- scale(gp_noise)
    } else {
      # Metodo di fallback
      cat("Metodo di correlazione spaziale non riconosciuto, uso GRF come fallback\n")

      # Usa un metodo più semplice basato su fields
      set.seed(random_seed)
      coords_matrix <- as.matrix(cell_df[, c("x", "y")])

      # Normalizza coordinate per evitare problemi numerici
      coords_norm <- scale(coords_matrix)

      # Crea la matrice di covarianza spaziale
      cov_mat <- exp(-rdist(coords_norm) / (spatial_params$spatial_range / 100))

      # Simuliamo un campo casuale gaussiano
      gp_noise <- t(chol(cov_mat)) %*% rnorm(nrow(coords_norm))

      # Normalizza
      gp_noise <- as.vector(scale(gp_noise))
    }
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
  # Preallochiamo la matrice di espressione - migliora efficienza
  expression_data <- matrix(0, nrow = N, ncol = n_genes)
  
  # Converti cluster_labels in interi una sola volta
  cl <- as.integer(cluster_labels)
  
  # Parallelizzazione massiccia per la generazione dell'espressione genica
  # Dividi i geni in chunk per parallelizzare
  gene_chunks <- split(seq_len(n_genes), ceiling(seq_len(n_genes)/min(50, ceiling(n_genes/10))))
  
  # Pre-calcola alcune strutture di dati comuni a tutti i geni
  # Crea una matrice di medie di espressione per tipo di cellula e gene
  # Organizzato come matrice [gene, cluster]
  all_mean_expr <- matrix(0, nrow = n_genes, ncol = k_cell_types)
  for (g in seq_len(n_genes)) {
    for (k in 1:k_cell_types) {
      all_mean_expr[g, k] <- mean_expression_list[[k]][g]
    }
  }
  
  # Configura multi-sessione con parallelizzazione bilanciata
  future::plan(future::multisession, workers = min(110, 40))
  # Ignora riferimenti per ridurre la memoria necessaria
  options(future.globals.onReference = "ignore")
  
  # Genera l'espressione genica in parallelo per chunk di geni
  expression_chunks <- future_lapply(gene_chunks, function(genes_subset) {
    # Alloca lo storage per l'espressione di questo chunk
    chunk_expression <- matrix(0, nrow = N, ncol = length(genes_subset))
    
    for (i in seq_along(genes_subset)) {
      g <- genes_subset[i]
      
      # Aggiungi variabilità cellula-specifica indipendente dal cluster
      cell_specific_effect <- rnorm(N, 0, cell_specific_params$cell_specific_noise_sd)

    # Calcola medie di espressione di base - versione molto più veloce usando indexing
    base_expr <- all_mean_expr[g, cl]

    # Applica effetto di ibridazione o gradienti
    if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df) &&
        exists("in_gradient", where = cell_df)) {
      # Usa l'approccio gradiente per i confini - versione vettorizzata

      # Identifica tutti i punti in zona gradiente
      gradient_points <- which(cell_df$in_gradient)

      if (length(gradient_points) > 0) {
        # Ottieni i cluster originali per tutti i punti gradiente contemporaneamente
        orig_clusters <- as.integer(cluster_labels[gradient_points])

        # Pre-calcola i pesi del gradiente per tutti i punti di confine
        # Usa esponente personalizzabile per controllo più fine della forma del gradiente
        gradient_weights <- (1 - cell_df$boundary_dist[gradient_points])^spatial_params$gradient_exponent

        # Crea una matrice di medie di espressione per tipo di cellula
        # Questo è più rapido che cercare ripetutamente nel mean_expression_list
        cluster_expr_matrix <- sapply(1:k_cell_types, function(k) mean_expression_list[[k]][g])

        # Prepara una matrice di pesi per tutti i tipi cellulari
        cell_type_weights <- matrix(0, nrow = length(gradient_points), ncol = k_cell_types)

        # Identifica, per ogni punto gradiente, il cluster più vicino diverso da quello originale
        for (i in seq_along(gradient_points)) {
          point_idx <- gradient_points[i]
          orig_cl <- orig_clusters[i]

          # Trova cluster vicini con soglia dinamica basata sull'ampiezza del gradiente
          # Usa una soglia maggiore per trovare più cluster vicini
          dist_threshold <- spatial_params$gradient_width * 1.5 * grid_resolution
          nearby_points <- which(dist_mat[point_idx, ] < dist_threshold)

          if (length(nearby_points) > 0) {
            nearby_cls <- as.integer(cluster_labels[nearby_points])
            nearby_cls <- nearby_cls[nearby_cls != orig_cl]

            if (length(nearby_cls) > 0) {
              # Conta le occorrenze di ciascun cluster vicino
              cls_counts <- table(nearby_cls)

              # Normalizza i conteggi per ottenere pesi relativi
              cls_weights <- as.numeric(cls_counts) / sum(cls_counts)

              # Per ogni tipo cellulare vicino, aggiungi il peso appropriato
              for (cl_idx in seq_along(cls_counts)) {
                cl_type <- as.integer(names(cls_counts)[cl_idx])
                cell_type_weights[i, cl_type] <- cls_weights[cl_idx]
              }
            }
          }
        }

        # Normalizza i pesi per tipo cellulare, preservando il peso totale del gradiente
        row_sums <- rowSums(cell_type_weights)
        valid_rows <- row_sums > 0
        if (any(valid_rows)) {
          cell_type_weights[valid_rows, ] <- cell_type_weights[valid_rows, ] / row_sums[valid_rows]

          # Calcola l'espressione ponderata di tutti i tipi cellulari in una singola operazione
          weighted_expr <- cell_type_weights %*% cluster_expr_matrix

          # Applica il mix con il peso del gradiente
          # Usa versione vettorizzata per tutti i punti in una volta
          base_expr[gradient_points] <- base_expr[gradient_points] * (1 - gradient_weights) +
                                        weighted_expr * gradient_weights
        }
      }
    } else if (hybrid_params$use_hybrid_cells) {
      # Approccio con cellule ibride - vettorizzato
      # Calcola l'effetto di tutti i cluster contemporaneamente

      # Crea una matrice dove ogni riga è la media di espressione per ogni tipo di cellula
      all_cluster_expr <- sapply(1:k_cell_types, function(k) mean_expression_list[[k]][g])

      # Applica l'effetto ibrido per tutti i punti e cluster in un'operazione matriciale
      # base_expr finale = base_expr attuale * (1-peso_ibrido) + espressione_altro_cluster * peso_ibrido

      # Calcola l'effetto ibrido complessivo
      hybrid_effect <- hybrid_matrix %*% all_cluster_expr

      # Calcola il peso complessivo dell'effetto ibrido su ogni cellula
      hybrid_weight <- rowSums(hybrid_matrix)

      # Applica solo a cellule che hanno un effetto ibrido
      hybrid_cells <- which(hybrid_weight > 0)
      if (length(hybrid_cells) > 0) {
        # Effetto ibrido totale per ogni cellula
        base_expr[hybrid_cells] <- base_expr[hybrid_cells] * (1 - hybrid_weight[hybrid_cells]) +
                                 hybrid_effect[hybrid_cells]
      }
    }

    # Combina con l'effetto cellula-specifico
    mu_vals <- base_expr + cell_specific_effect

    # Aggiungi effetto dei moduli genici se abilitato
    if (cell_specific_params$use_gene_modules && !is.null(module_noise)) {
      # Aggiungi il rumore correlato del modulo genico a cui appartiene questo gene
      mu_vals <- mu_vals + module_noise[, g]
    }

    # Aggiungi correlazione spaziale se richiesta
    if (use_spatial_correlation) {
      mu_vals <- mu_vals + spatial_params$spatial_noise_intensity * gp_noise

      # Aggiungi noise casuale addizionale per confondere i pattern
      random_noise <- rnorm(length(mu_vals), 0, spatial_params$random_noise_sd)
      mu_vals <- mu_vals + random_noise
    }

    # Genera conteggi di espressione
    if (g %in% stable_genes) {
      # Modello sub-Poisson: Binomiale con p alto e n moderato
      p <- 0.9
      n_trial <- round(exp(mu_vals)/(1-p))
      raw_counts <- rbinom(N, n_trial, p)
    } else {
      # Negative Binomial con dispersione variabile spazialmente
      raw_counts <- rnbinom(N, mu = exp(mu_vals), size = dispersion_param)
    }

    # Applica l'effetto della dimensione della libreria
    # Scala i conteggi in base alla dimensione relativa della libreria
    # (mantenendo la somma totale constante per evitare bias nelle medie)
    scaled_counts <- raw_counts * (library_size / mean(library_size))
    # Arrotonda a numeri interi (conteggi)
    chunk_expression[, i] <- round(scaled_counts)

    # Applica dropout in base al modello specificato - versione vettorizzata
    if (dropout_params$expression_dependent_dropout) {
      # Dropout dipendente dal livello di espressione - usa una funzione logistica
      # Più il gene è espresso, meno probabile è il dropout

      # Definisci una funzione vettorizzata per normalizzare tra 0 e 1
      scale01_vec <- function(x) {
        if (all(x == x[1])) return(rep(0.5, length(x)))
        (x - min(x)) / (max(x) - min(x))
      }

      # Normalizza l'espressione del gene corrente
      norm_expr <- scale01_vec(chunk_expression[, i])

      # Calcola la probabilità di dropout con una funzione logistica vettorizzata
      # Formula: p(dropout) = 1 / (1 + exp((expr - midpoint) * steepness))
      dropout_prob_expr <- 1 / (1 + exp((norm_expr - dropout_params$dropout_curve_midpoint) *
                                       dropout_params$dropout_curve_steepness))

      # Combina con il dropout spaziale base (media pesata) - operazione vettorizzata
      dropout_prob <- 0.7 * dropout_prob_expr + 0.3 * base_dropout

      # Tronca i valori al range [0,1] in un'unica operazione
      dropout_prob <- pmin(pmax(dropout_prob, 0), 1)

      # Applica dropout in modo vettorizzato
      zero_idx <- runif(N) < dropout_prob  # Più efficiente di rbinom() per vettori lunghi
      chunk_expression[zero_idx, i] <- 0
    } else {
      # Modello di dropout originale (solo spaziale) - vettorizzato
      zero_idx <- runif(N) < base_dropout  # Più efficiente di rbinom() per vettori lunghi
      chunk_expression[zero_idx, i] <- 0
    }
    
    } # Fine del ciclo for sui geni di questo chunk
    
    return(chunk_expression)
  }) # Fine del future_lapply
  
  # Combina i risultati dei chunk in una singola matrice di espressione
  for (i in seq_along(gene_chunks)) {
    genes_subset <- gene_chunks[[i]]
    expression_data[, genes_subset] <- expression_chunks[[i]]
  }

  # Funzione helper per normalizzare tra 0 e 1
  scale01 <- function(x) {
    if (max(x) == min(x)) return(rep(0.5, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }

  # Prepara l'output
  result <- list(
    coordinates       = cell_df[, c("x", "y")],
    intensity_cluster = cell_df$intensity_cluster,
    expression        = expression_data,
    threshold_used    = threshold_value,
    library_size      = library_size,  # Aggiungiamo libreria per ogni spot
    dispersion_param  = dispersion_param,  # Parametro di dispersione per ogni spot
    gene_modules      = gene_modules,  # Informazioni sui moduli di geni
    parameters        = list(
      image_path = image_path,
      pixel_size_um = pixel_size_um,   # Aggiungiamo il parametro della dimensione del pixel
      n_cells = n_cells,
      n_genes = n_genes,
      k_cell_types = k_cell_types,
      difficulty_level = difficulty_level,
      grid_mode = grid_mode,
      grid_resolution = grid_resolution,
      use_fixed_grid = use_fixed_grid,  # Usiamo griglia fissa?
      fixed_grid_width_mm = fixed_grid_width_mm,  # Dimensioni della griglia in mm
      fixed_grid_height_mm = fixed_grid_height_mm,
      use_spatial_correlation = use_spatial_correlation,
      correlation_method = correlation_method,
      marker_params = marker_params,
      spatial_params = spatial_params,
      dropout_params = dropout_params,
      library_size_params = library_size_params,
      hybrid_params = hybrid_params,
      cell_specific_params = cell_specific_params
    )
  )

  # Aggiungi informazioni aggiuntive se disponibili
  if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df)) {
    result$boundary_dist = cell_df$boundary_dist
  }

  # Calcola e aggiungi l'autocorrelazione spaziale globale (Moran's I) per alcuni geni
  if (requireNamespace("spdep", quietly = TRUE)) {
    coords_matrix <- as.matrix(cell_df[, c("x", "y")])
    nb <- spdep::dnearneigh(coords_matrix, 0, spatial_params$spatial_range)
    lw <- spdep::nb2listw(nb, style = "W")

    # Calcola Moran's I per i primi 10 geni (o tutti se < 10)
    moran_genes <- min(10, ncol(expression_data))
    moran_i <- numeric(moran_genes)

    for (g in 1:moran_genes) {
      # Gestisci casi degeneri (tutti zeri)
      if (all(expression_data[, g] == 0)) {
        moran_i[g] <- NA
      } else {
        tryCatch({
          moran_i[g] <- spdep::moran.test(expression_data[, g], lw)$estimate[1]
        }, error = function(e) {
          moran_i[g] <- NA
        })
      }
    }

    result$spatial_autocorrelation <- list(
      moran_i = moran_i,
      genes_tested = 1:moran_genes
    )
  }

  # ========================
  # 6) Visualizzazione
  # ========================
  # Plot della distribuzione spaziale dei cluster
  if (grid_mode) {
    # Per griglia, usiamo geom_tile
    p <- ggplot(cell_df, aes(x = x, y = y, fill = intensity_cluster)) +
      geom_tile(width = grid_resolution, height = grid_resolution) +
      scale_y_reverse() +
      coord_fixed() +
      theme_minimal() +
      labs(title = sprintf("Visium HD (2µm) - Livello difficoltà: %s", difficulty_level),
           subtitle = sprintf("Griglia %dµm, %d bin, %d geni (marker/tipo: %d, fold: %.1f)",
                             grid_resolution, nrow(cell_df), n_genes,
                             marker_params$marker_genes_per_type,
                             marker_params$marker_expression_fold),
           fill = "Cell Type")
  } else {
    # Per sampling casuale, usiamo punti
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
  }

  print(p)

  # Plot della dimensione libreria
  if (exists("library_size")) {
    p_lib <- ggplot(data.frame(x = cell_df$x, y = cell_df$y, lib_size = library_size),
                   aes(x = x, y = y, fill = lib_size)) +
      geom_tile(width = grid_resolution, height = grid_resolution) +
      scale_fill_viridis_c(option = "plasma") +
      scale_y_reverse() +
      coord_fixed() +
      theme_minimal() +
      labs(title = "Distribuzione della dimensione della libreria",
           fill = "Library size")

    print(p_lib)
  }

  # Plot dell'autocorrelazione spaziale se calcolata
  if (exists("spatial_autocorrelation", where = result)) {
    p_auto <- ggplot(data.frame(gene = result$spatial_autocorrelation$genes_tested,
                               moran_i = result$spatial_autocorrelation$moran_i),
                    aes(x = gene, y = moran_i)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = "Indice di Moran per geni selezionati",
           x = "Gene index", y = "Moran's I")

    print(p_auto)
  }

  # Salva i plot se richiesto
  if (!is.null(output_plot)) {
    # Salva il plot principale
    ggplot2::ggsave(output_plot, plot = p, device = "png", dpi = 300)

    # Salva anche i plot aggiuntivi se esistono
    if (exists("p_lib")) {
      ggplot2::ggsave(sub(".png$", "_library_size.png", output_plot),
                     plot = p_lib, device = "png", dpi = 300)
    }

    if (exists("p_auto")) {
      ggplot2::ggsave(sub(".png$", "_autocorrelation.png", output_plot),
                     plot = p_auto, device = "png", dpi = 300)
    }
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

# Esempio 1: Visium HD-like con griglia a 2µm
# simulate_spatial_transcriptomics(
#   image_path = "images/colon.png",
#   output_path = "data/simulated_visiumhd_simple.rds",
#   output_plot = "results/simulated_visiumhd_simple.png",
#   grid_mode = TRUE,                    # Attiva la modalità griglia
#   grid_resolution = 2,                 # Risoluzione 2µm come Visium HD
#   difficulty_level = "medium",
#   n_genes = 100,
#   k_cell_types = 5,
#   correlation_method = "grf",          # Gaussian Random Field
#   spatial_params = list(
#     gradient_regions = TRUE,           # Attiva gradienti tra regioni
#     gradient_width = 5,                # Larghezza del gradiente (in unità griglia)
#     gradient_exponent = 1.5,           # Forma del gradiente
#     spatial_noise_intensity = 1.0,
#     spatial_range = 30,
#     random_noise_sd = 0.2
#   ),
#   dropout_params = list(
#     expression_dependent_dropout = TRUE,  # Dropout dipendente dal livello di espressione
#     dropout_curve_midpoint = 0.5,
#     dropout_curve_steepness = 5,
#     cell_type_dispersion_effect = 0.2    # Effetto del tipo cellulare sulla dispersione
#   ),
#   library_size_params = list(
#     mean_library_size = 10000,         # Media conteggi per spot
#     library_size_cv = 0.3,             # Coefficiente di variazione
#     spatial_effect_on_library = 0.5,   # Effetto spaziale sulla dimensione libreria
#     cell_type_effect = TRUE            # Effetto del tipo cellulare
#   ),
#   cell_specific_params = list(
#     cell_specific_noise_sd = 0.2,
#     use_gene_modules = TRUE,
#     n_gene_modules = 5,                # Moduli di co-espressione genica
#     module_correlation = 0.6           # Correlazione tra geni dello stesso modulo
#   )
# )

# Esempio 2: dataset "facile" con metodo originale (pre-HD)
# simulate_spatial_transcriptomics(
#   image_path = "images/colon.png",
#   output_path = "data/simulated_basic.rds",
#   output_plot = "results/simulated_basic.png",
#   grid_mode = FALSE,                    # Disattiva la modalità griglia
#   difficulty_level = "easy",
#   n_cells = 20000,
#   n_genes = 100
# )

# Esempio 3: dataset difficile con pattern CAR (Conditional Autoregressive)
# simulate_spatial_transcriptomics(
#   image_path = "images/granuloma.png",
#   output_path = "data/simulated_visiumhd_car.rds",
#   output_plot = "results/simulated_visiumhd_car.png",
#   grid_mode = TRUE,
#   grid_resolution = 2,
#   difficulty_level = "hard",
#   n_genes = 100,
#   correlation_method = "car",           # Conditional Autoregressive model
#   spatial_params = list(
#     gradient_regions = TRUE,
#     gradient_width = 8,
#     gradient_exponent = 1.5,
#     spatial_noise_intensity = 1.8,
#     spatial_range = 20,
#     random_noise_sd = 0.15
#   ),
#   cell_specific_params = list(
#     use_gene_modules = TRUE,
#     n_gene_modules = 4,
#     module_correlation = 0.8
#   )
# )

# Esempio 4: dataset Visium HD con migliorie biologiche
simulate_spatial_transcriptomics(
  image_path = "images/granuloma.png",
  output_path = "data/simulated_visiumhd.rds",
  output_plot = "results/simulated_visiumhd.png",
  pixel_size_um = 10,           # Ogni pixel vale 10μm nell'immagine granuloma.png
  grid_mode = TRUE,
  grid_resolution = 2,          # 2μm x 2μm per bin (standard Visium HD)
  use_fixed_grid = TRUE,        # Usa griglia fissa standard Visium HD
  fixed_grid_width_mm = 6.5,    # 6.5mm x 6.5mm (standard Visium HD)
  fixed_grid_height_mm = 6.5,
  threshold_value = 0.6,
  difficulty_level = "medium",
  n_genes = 150,
  k_cell_types = 7,
  # Correlazione forte ma con range corto
  correlation_method = "grf",
  spatial_params = list(
    spatial_noise_intensity = 0.8,     # Intensità più alta
    spatial_range = 25,                # Range più corto
    random_noise_sd = 0.1,            # Ridotto per favorire pattern spaziali
    gradient_regions = TRUE,
    gradient_width = 15,               # Gradiente più ampio
    gradient_exponent = 1.5            # Transizione non lineare del gradiente
  ),
  # Marker con forte sovrapposizione
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 2.0,      # Differenza moderatamente alta
    marker_overlap_fold = 0.1          # Moderata sovrapposizione
  ),
  # Dropout e dispersione
  dropout_params = list(
    dropout_range = c(0.05, 0.3),
    dispersion_range = c(3.0, 1.2),
    cell_type_dispersion_effect = 0.3,  # Effetto specifico di tipo cellulare
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.3,
    dropout_curve_steepness = 8         # Moderatamente ripida
  ),
  # Effetti sulla dimensione libreria
  library_size_params = list(
    mean_library_size = 8000,
    library_size_cv = 0.4,              # Variabilità moderata-alta
    spatial_effect_on_library = 0.7,    # Moderato effetto spaziale
    cell_type_effect = TRUE             # Effetto del tipo cellulare
  ),
  # Parametri di correlazione di geni
  cell_specific_params = list(
    cell_specific_noise_sd = 0.15,
    use_gene_modules = TRUE,
    n_gene_modules = 6,                 # Programmi di espressione ben definiti
    module_correlation = 0.8            # Correlazione moderata-alta
  )
)
