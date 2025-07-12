#' Calcola le distanze e le densità locali per ogni cella
#'
#' Calcola le distanze tra celle/spot e metriche di densità locale
#' utilizzate poi per i modelli di dropout e dispersione.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param chunk_size Dimensione del chunk per parallelizzazione
#' @param random_seed Seed per riproducibilità
#' @return Lista con distanze e densità calcolate
# Utilizziamo lapply standard
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
  # Calcolo chunked di distanze per evitare allocazione densa
  # NOTA: dist_mat completo evitato per risparmio memoria
  mean_dist <- numeric(N)
  local_density <- numeric(N)
  for (i in 1:N) {
    same_cluster <- which(cluster_labels == cluster_labels[i])
    dists <- sqrt((coords$x[i] - coords$x[same_cluster])^2 + (coords$y[i] - coords$y[same_cluster])^2)
    mean_dist[i] <- mean(dists)
    q <- quantile(dists, 0.1)
    local_density[i] <- mean(dists < q)
  }
  
  # Calcola la dimensione ottimale del chunk se non fornita
  if (is.null(chunk_size)) {
    chunk_size <- max(1, ceiling(N/500))
  }
  
  # Dividi in chunks per ottimizzare l'elaborazione
  chunks <- split(1:N, ceiling(seq_along(1:N)/chunk_size))
  
  # Calcola la distanza media e densità locale SENZA dist_mat
  mean_dist <- numeric(N)
  local_density <- numeric(N)
  for (i in 1:N) {
    same_cluster <- which(cluster_labels == cluster_labels[i])
    dists <- sqrt((coords$x[i] - coords$x[same_cluster])^2 + (coords$y[i] - coords$y[same_cluster])^2)
    mean_dist[i] <- mean(dists)
    q <- quantile(dists, 0.1)
    local_density[i] <- mean(dists < q)
  }
  # Restituisci solo le statistiche
  return(list(
    mean_dist = mean_dist,
    local_density = local_density
  ))
}