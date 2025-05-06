#' Esegue clustering sui pixel dell'immagine
#'
#' @param img_df_thresh Dataframe dei pixel filtrati
#' @param k_cell_types Numero di tipi cellulari
#' @param random_seed Seed per riproducibilità
#' @param clustering_method Metodo di clustering: "spatial_kmeans", "kmeans++", o "slic"
#' @param spatial_weight Peso della componente spaziale (0-1+)
#' @param estimate_k Se stimare automaticamente k
#' @param k_estimation_method Metodo per stimare k: "silhouette" o "elbow"
#' @return Dataframe con informazioni di cluster aggiunte
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom dplyr mutate
#' @export
cluster_image <- function(
  img_df_thresh,
  k_cell_types,
  random_seed = 123,
  clustering_method = "spatial_kmeans",
  spatial_weight = 0.5,
  estimate_k = FALSE,
  k_estimation_method = "silhouette"
) {
  set.seed(random_seed)
  
  # Se necessario, stima automaticamente il numero di cluster
  if (estimate_k) {
    # Seleziona le caratteristiche in base al metodo di clustering
    if (clustering_method == "spatial_kmeans") {
      # Per spatial_kmeans, considera intensità e posizione
      spatial_coords <- scale(as.matrix(img_df_thresh[, c("x", "y")]))
      intensity_vals <- scale(as.matrix(img_df_thresh$value))
      features <- cbind(intensity_vals, spatial_coords * spatial_weight)
    } else {
      # Per kmeans++, considera solo l'intensità
      features <- as.matrix(img_df_thresh$value)
    }
    
    # Stima k con il metodo scelto
    if (k_estimation_method == "silhouette") {
      k_cell_types <- estimate_k_silhouette(features, max_k = min(20, floor(sqrt(nrow(features)))), random_seed = random_seed)
    } else if (k_estimation_method == "elbow") {
      k_cell_types <- estimate_k_elbow(features, max_k = min(20, floor(sqrt(nrow(features)))), random_seed = random_seed)
    } else {
      warning("Metodo di stima k non riconosciuto. Usando k_cell_types fornito.")
    }
  }
  
  # Seleziona il metodo di clustering
  if (clustering_method == "kmeans++") {
    # Metodo kmeans++ originale (solo intensità)
    km_result <- KMeans_rcpp(
      as.matrix(img_df_thresh$value),
      clusters    = k_cell_types,
      num_init    = 5,
      initializer = 'kmeans++',
      seed        = random_seed
    )
  } else if (clustering_method == "spatial_kmeans") {
    # Nuovo metodo spatial_kmeans (intensità + posizione)
    km_result <- spatial_kmeans(
      img_df_thresh, 
      k_cell_types, 
      spatial_weight, 
      random_seed
    )
  } else if (clustering_method == "slic") {
    # Clustering con SLIC (superpixel)
    km_result <- slic_clustering(
      img_df_thresh, 
      k_cell_types, 
      random_seed
    )
  } else {
    stop("Metodo di clustering non riconosciuto. Scegliere tra 'spatial_kmeans', 'kmeans++' o 'slic'.")
  }
  
  # Aggiungo il cluster al dataframe
  img_df_thresh <- img_df_thresh %>%
    mutate(intensity_cluster = factor(km_result$clusters))
  
  return(img_df_thresh)
}

#' Esegue clustering spatial_kmeans che considera sia intensità che posizione
#'
#' @param img_df_thresh Dataframe dei pixel filtrati
#' @param k_cell_types Numero di tipi cellulari
#' @param spatial_weight Peso della componente spaziale
#' @param random_seed Seed per riproducibilità
#' @return Lista con i risultati del clustering
#' @importFrom ClusterR KMeans_rcpp
spatial_kmeans <- function(img_df_thresh, k_cell_types, spatial_weight = 0.5, random_seed = 123) {
  # 1. Normalizzazione delle caratteristiche
  spatial_coords <- scale(as.matrix(img_df_thresh[, c("x", "y")]))  # Normalizza coordinate
  intensity_vals <- scale(as.matrix(img_df_thresh$value))           # Normalizza intensità
  
  # 2. Combinazione con peso
  combined_features <- cbind(intensity_vals, spatial_coords * spatial_weight)
  
  # 3. Clustering k-means++
  km_combined <- KMeans_rcpp(
    combined_features,
    clusters    = k_cell_types,
    num_init    = 5,
    initializer = 'kmeans++',
    seed        = random_seed
  )
  
  return(km_combined)
}

#' Esegue clustering con SLIC (Simple Linear Iterative Clustering)
#'
#' @param img_df_thresh Dataframe dei pixel filtrati
#' @param k_cell_types Numero di tipi cellulari finali
#' @param random_seed Seed per riproducibilità
#' @return Lista con i risultati del clustering
#' @importFrom stats aggregate kmeans
slic_clustering <- function(img_df_thresh, k_cell_types, random_seed = 123) {
  set.seed(random_seed)
  
  # 1. Verifica disponibilità dei pacchetti necessari
  if (!requireNamespace("supercells", quietly = TRUE) || 
      !requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    warning("I pacchetti 'supercells', 'terra' o 'sf' non sono disponibili. Utilizzando spatial_kmeans come fallback.")
    # Fallback a spatial_kmeans
    return(spatial_kmeans(img_df_thresh, k_cell_types, spatial_weight = 0.3, random_seed = random_seed))
  } 
  
  # 2. Conversione del dataframe in una griglia raster
  tryCatch({
    # Determina le dimensioni dell'immagine
    width <- max(img_df_thresh$x)
    height <- max(img_df_thresh$y)
    
    # Converti i dati in una matrice
    # Inizializza una matrice vuota
    img_matrix <- matrix(NA, nrow = height, ncol = width)
    
    # Riempie la matrice con i valori di intensità
    for (i in 1:nrow(img_df_thresh)) {
      x <- img_df_thresh$x[i]
      y <- img_df_thresh$y[i]
      if (x >= 1 && x <= width && y >= 1 && y <= height) {
        img_matrix[y, x] <- img_df_thresh$value[i]
      }
    }
    
    # Converti la matrice in un oggetto SpatRaster (richiesto da supercells)
    rast <- terra::rast(img_matrix)
    terra::ext(rast) <- c(0, width, 0, height)  # Imposta l'estensione
    
    # Calcola il numero ottimale di superpixel
    n_superpixels <- min(k_cell_types * 15, floor(nrow(img_df_thresh) / 10))
    
    # Usa supercells (SLIC)
    sc <- supercells::supercells(
      x = rast,           # Input raster data
      k = n_superpixels,  # Number of superpixels
      compactness = 10,   # Controlla la compattezza dei superpixel
      dist_fun = "euclidean"
    )
    
    # 5. Estrai i valori medi per ogni superpixel
    # Converti la geometria sf in un dataframe semplice
    superpixel_data <- sf::st_drop_geometry(sc)
    
    # 6. Clustering dei superpixel basato sui valori medi
    # Estrai la prima banda se ci sono più bande
    if (ncol(superpixel_data) > 2 && "value" %in% colnames(superpixel_data)) {
      superpixel_clusters <- kmeans(superpixel_data$value, centers = k_cell_types, nstart = 5)
    } else {
      # Usa il primo valore numerico disponibile (oltre a cell)
      numeric_cols <- sapply(superpixel_data, is.numeric)
      numeric_cols["cell"] <- FALSE  # Escludi la colonna cell
      if (sum(numeric_cols) > 0) {
        first_value_col <- names(numeric_cols)[which(numeric_cols)[1]]
        superpixel_clusters <- kmeans(superpixel_data[[first_value_col]], centers = k_cell_types, nstart = 5)
      } else {
        # Fallback se non ci sono colonne numeriche
        warning("Nessuna colonna numerica trovata in superpixel_data, usando spatial_kmeans come fallback")
        result <- spatial_kmeans(img_df_thresh, k_cell_types, spatial_weight = 0.3, random_seed = random_seed)
        return(result)
      }
    }
    
    # 7. Crea una mappa per i cluster
    cluster_lookup <- superpixel_clusters$cluster
    names(cluster_lookup) <- superpixel_data$cell  # Usa l'ID della cella come chiave
    
    # Estrai i centroidi dei superpixel per il nearest neighbor
    superpixel_centers <- sf::st_coordinates(sf::st_centroid(sc))
    
    # Inizializza il vettore dei cluster finali
    final_clusters <- numeric(nrow(img_df_thresh))
    
    # Mappiamo i cluster dai superpixel alle celle originali usando nearest neighbor
    for (i in 1:nrow(img_df_thresh)) {
      # Trova il superpixel più vicino a questo punto
      distances <- sqrt((superpixel_centers[, 1] - img_df_thresh$x[i])^2 + 
                       (superpixel_centers[, 2] - img_df_thresh$y[i])^2)
      nearest_superpixel <- which.min(distances)
      
      # Assegna il cluster direttamente
      final_clusters[i] <- superpixel_clusters$cluster[nearest_superpixel]
    }
  }, error = function(e) {
    message("Errore nell'uso di supercells: ", e$message)
    message("Utilizzando spatial_kmeans come fallback...")
    result <- spatial_kmeans(img_df_thresh, k_cell_types, spatial_weight = 0.3, random_seed = random_seed)
    final_clusters <- result$clusters
  })
  
  # 8. Creazione di un oggetto risultato compatibile con KMeans_rcpp
  result <- list(
    clusters = final_clusters,
    WCSS_per_cluster = rep(0, k_cell_types),  # Placeholder per withiness
    centroids = matrix(0, nrow = k_cell_types, ncol = 1)  # Placeholder per centroidi
  )
  
  return(result)
}

#' Stima il numero ottimale di cluster usando il metodo elbow
#'
#' @param features Matrice delle caratteristiche
#' @param max_k Numero massimo di cluster da testare
#' @param random_seed Seed per riproducibilità
#' @return Numero ottimale di cluster stimato
#' @importFrom ClusterR KMeans_rcpp
estimate_k_elbow <- function(features, max_k = 20, random_seed = 123) {
  set.seed(random_seed)
  
  # Limita max_k in base alla dimensione dei dati
  max_k <- min(max_k, floor(sqrt(nrow(features))))
  max_k <- max(3, max_k)  # Almeno 3 cluster
  
  # Calcola WCSS (Within-Cluster Sum of Squares) per diversi valori di k
  wcss <- numeric(max_k - 1)
  for (k in 2:max_k) {
    kmeans_result <- KMeans_rcpp(features, clusters = k, seed = random_seed)
    wcss[k-1] <- sum(kmeans_result$WCSS_per_cluster)
  }
  
  # Trova il punto di gomito calcolando il cambio nella derivata seconda
  d1 <- diff(wcss)
  d2 <- diff(d1)
  k_optimal <- which.max(abs(d2)) + 2
  
  # Fallback a un valore mediano se il metodo non trova un chiaro punto di gomito
  if (is.na(k_optimal) || k_optimal < 2) {
    k_optimal <- ceiling(max_k / 2)
  }
  
  return(k_optimal)
}

#' Stima il numero ottimale di cluster usando il metodo silhouette
#'
#' @param features Matrice delle caratteristiche
#' @param max_k Numero massimo di cluster da testare
#' @param random_seed Seed per riproducibilità
#' @return Numero ottimale di cluster stimato
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom cluster silhouette
#' @importFrom stats dist
estimate_k_silhouette <- function(features, max_k = 20, random_seed = 123) {
  set.seed(random_seed)
  
  # Limita max_k in base alla dimensione dei dati
  max_k <- min(max_k, floor(sqrt(nrow(features))))
  max_k <- max(3, max_k)  # Almeno 3 cluster
  
  # Calcola il silhouette score per diversi valori di k
  sil_scores <- numeric(max_k - 1)
  
  # Prima calcola la matrice di distanza (una sola volta per efficienza)
  dist_matrix <- dist(features)
  
  for (k in 2:max_k) {
    kmeans_result <- KMeans_rcpp(features, clusters = k, seed = random_seed)
    sil <- silhouette(kmeans_result$clusters, dist_matrix)
    sil_scores[k-1] <- mean(sil[, "sil_width"])
  }
  
  # Trova il k con il miglior punteggio di silhouette
  k_optimal <- which.max(sil_scores) + 1
  
  # Fallback a un valore mediano se il metodo non trova un chiaro massimo
  if (is.na(k_optimal) || k_optimal < 2) {
    k_optimal <- ceiling(max_k / 2)
  }
  
  return(k_optimal)
}