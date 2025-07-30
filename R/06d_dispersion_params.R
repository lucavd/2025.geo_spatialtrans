#' Calcola i parametri di dispersione per la simulazione
#'
#' Genera parametri di dispersione variabili, tenendo conto della distanza
#' dai confini, della posizione spaziale e del tipo cellulare.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param mean_dist Vettore con le distanze medie per ogni cella
#' @param dropout_params Parametri di dropout, incluso il range di dispersione
#' @param k_cell_types Numero di tipi cellulari
#' @param random_seed Seed per riproducibilità
#' @return Vettore con parametri di dispersione per ogni cella
#' @importFrom scales rescale
#' @export
calculate_dispersion_params <- function(
  cell_df,
  mean_dist,
  dropout_params = list(
    dispersion_range = c(2.0, 1.0),
    cell_type_dispersion_effect = 0.2
  ),
  k_cell_types,
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai i cluster labels
  cluster_labels <- cell_df$intensity_cluster
  
  # Verifica che dispersion_range sia un vettore numerico di lunghezza 2
  if (!is.numeric(dropout_params$dispersion_range) || length(dropout_params$dispersion_range) != 2) {
    warning("dispersion_range deve essere un vettore numerico di lunghezza 2, impostato a c(2.0, 1.0)")
    dropout_params$dispersion_range <- c(2.0, 1.0)
  }
  
  # Calcola la dispersione base in base alla distanza dai confini o alla distanza media
  if ("boundary_dist" %in% colnames(cell_df) && !all(is.na(cell_df$boundary_dist))) {
    # Se utilizziamo gradienti, facciamo variare la dispersione in base alla distanza dal confine
    dispersion_param <- dropout_params$dispersion_range[2] +
      cell_df$boundary_dist * (dropout_params$dispersion_range[1] - dropout_params$dispersion_range[2])
  } else {
    # Altrimenti usiamo il metodo originale basato sulla distanza media
    # Gestisci il caso in cui mean_dist sia NULL o contenga NA
    if (is.null(mean_dist) || any(is.na(mean_dist))) {
      # Se mean_dist non è valido, usa un valore costante
      dispersion_param <- rep(mean(dropout_params$dispersion_range), nrow(cell_df))
    } else {
      dispersion_param <- scales::rescale(mean_dist, to = dropout_params$dispersion_range)
    }
  }
  
  # Aggiungi effetto del tipo cellulare sulla dispersione
  if (!is.null(dropout_params$cell_type_dispersion_effect) && 
      is.numeric(dropout_params$cell_type_dispersion_effect)) {
    type_dispersion_effects <- runif(
      k_cell_types,
      min = 1 - dropout_params$cell_type_dispersion_effect,
      max = 1 + dropout_params$cell_type_dispersion_effect
    )
    
    # Applica effetto moltiplicativo per tipo cellulare
    for (k in 1:k_cell_types) {
      idx <- which(cluster_labels == k)
      if (length(idx) > 0) {
        dispersion_param[idx] <- dispersion_param[idx] * type_dispersion_effects[k]
      }
    }
  }
  
  return(dispersion_param)
}