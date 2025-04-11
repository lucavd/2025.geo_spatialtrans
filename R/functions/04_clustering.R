#' Esegue clustering kmeans++ sui pixel dell'immagine
#'
#' @param img_df_thresh Dataframe dei pixel filtrati
#' @param k_cell_types Numero di tipi cellulari
#' @param random_seed Seed per riproducibilità
#' @return Dataframe con informazioni di cluster aggiunte
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom dplyr mutate
#' @export
cluster_image <- function(img_df_thresh, k_cell_types, random_seed = 123) {
  set.seed(random_seed)
  
  # Eseguo kmeans++ sull'intensità dei pixel
  km_intensity <- KMeans_rcpp(
    as.matrix(img_df_thresh$value),
    clusters    = k_cell_types,
    num_init    = 5,
    initializer = 'kmeans++',
    seed        = random_seed
  )
  
  # Aggiungo il cluster al dataframe
  img_df_thresh <- img_df_thresh %>%
    mutate(intensity_cluster = factor(km_intensity$clusters))
  
  return(img_df_thresh)
}