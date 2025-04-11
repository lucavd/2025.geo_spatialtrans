#' Calcola le distanze e le densità locali per ogni cella
#'
#' Calcola le distanze tra celle/spot e metriche di densità locale
#' utilizzate poi per i modelli di dropout e dispersione.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param chunk_size Dimensione del chunk per parallelizzazione
#' @param random_seed Seed per riproducibilità
#' @return Lista con distanze e densità calcolate
#' @importFrom future.apply future_lapply
#' @export
calculate_spatial_distances <- function(
  cell_df,
  chunk_size = NULL,
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai i dati necessari
  N <- nrow(cell_df)
  cluster_labels <- cell_df$intensity_cluster
  coords <- cell_df %>% dplyr::select(x, y)
  dist_mat <- as.matrix(dist(coords))
  
  # Calcola la dimensione ottimale del chunk se non fornita
  if (is.null(chunk_size)) {
    chunk_size <- max(1, ceiling(N/500))
  }
  
  # Dividi in chunks per ottimizzare l'elaborazione
  chunks <- split(1:N, ceiling(seq_along(1:N)/chunk_size))
  
  # Calcola la distanza media per ciascuna cellula rispetto alle altre del proprio cluster
  mean_dist <- future_lapply(chunks, function(chunk_idx) {
    result <- numeric(length(chunk_idx))
    for (j in seq_along(chunk_idx)) {
      i <- chunk_idx[j]
      cl <- cluster_labels[i]
      same_cluster <- which(cluster_labels == cl)
      result[j] <- mean(dist_mat[i, same_cluster])
    }
    return(result)
  }, future.scheduling = 1, future.chunk.size = NULL, future.seed = TRUE) %>% unlist()
  
  # Calcola la densità locale (per il modello di dropout)
  local_density <- future_lapply(chunks, function(chunk_idx) {
    result <- numeric(length(chunk_idx))
    for (j in seq_along(chunk_idx)) {
      i <- chunk_idx[j]
      row <- dist_mat[i,]
      q <- quantile(row, 0.1)
      result[j] <- mean(row < q)
    }
    return(result)
  }, future.scheduling = 1, future.chunk.size = NULL, future.seed = TRUE) %>% unlist()
  
  # Restituisci i risultati
  return(list(
    dist_mat = dist_mat,
    mean_dist = mean_dist,
    local_density = local_density
  ))
}