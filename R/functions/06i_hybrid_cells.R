#' Identifica e genera cellule ibride tra tipi cellulari
#'
#' Trova celle di confine e crea profili di espressione ibridi
#' mescolando caratteristiche di tipi cellulari adiacenti.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param dist_mat Matrice delle distanze tra celle
#' @param k_cell_types Numero di tipi cellulari
#' @param hybrid_params Parametri per le cellule ibride
#' @param random_seed Seed per riproducibilità
#' @return Matrice di ibridazione per ogni cella e tipo cellulare
#' @export
generate_hybrid_cells <- function(
  cell_df,
  dist_mat,
  k_cell_types,
  hybrid_params = list(
    use_hybrid_cells = TRUE,
    max_hybrid_pairs = 1000,
    hybrid_intensity_range = c(0.2, 0.5)
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Numero di celle
  N <- nrow(cell_df)
  
  # Crea una matrice di ibridazione inizializzata a zero
  hybrid_matrix <- matrix(0, nrow = N, ncol = k_cell_types)
  
  # Procedi solo se è richiesta l'ibridazione
  if (hybrid_params$use_hybrid_cells) {
    # Estrai i cluster labels
    cluster_labels <- cell_df$intensity_cluster
    
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
      hybrid_pairs <- hybrid_pairs[sample(length(hybrid_pairs), hybrid_params$max_hybrid_pairs)]
    }
    
    # Crea una matrice di ibridazione 
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
  
  return(hybrid_matrix)
}