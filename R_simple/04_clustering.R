#' Esegue clustering sui pixel dell'immagine
#'
#' @param img_df_thresh Dataframe dei pixel filtrati
#' @param k_cell_types Numero di tipi cellulari
#' @param random_seed Seed per riproducibilità
#' @param spatial_weight Peso della componente spaziale (0-1+)
#' @return Dataframe con informazioni di cluster aggiunte
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom dplyr mutate
#' @export
cluster_image <- function(
  img_df_thresh,
  k_cell_types,
  random_seed = 123,
  spatial_weight = 0.5
) {
  set.seed(random_seed)
  
  # Esegui spatial k-means
  km_result <- spatial_kmeans(
    img_df_thresh, 
    k_cell_types, 
    spatial_weight, 
    random_seed
  )
  
  # Aggiungi il cluster al dataframe
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
  if (!requireNamespace("ClusterR", quietly = TRUE)) {
    stop("Il pacchetto 'ClusterR' è necessario per questa funzione")
  }
  
  # 1. Normalizzazione delle caratteristiche
  spatial_coords <- scale(as.matrix(img_df_thresh[, c("x", "y")]))  # Normalizza coordinate
  intensity_vals <- scale(as.matrix(img_df_thresh$value))           # Normalizza intensità
  
  # 2. Combinazione con peso
  combined_features <- cbind(intensity_vals, spatial_coords * spatial_weight)
  
  # 3. Clustering k-means++
  km_combined <- ClusterR::KMeans_rcpp(
    combined_features,
    clusters    = k_cell_types,
    num_init    = 5,
    initializer = 'kmeans++',
    seed        = random_seed
  )
  
  return(km_combined)
}