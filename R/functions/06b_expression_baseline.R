#' Genera profili di espressione baseline per ciascun tipo cellulare
#'
#' Crea i profili di espressione di base per ciascun tipo cellulare,
#' includendo effetti di marker genici specifici e sovrapposizione.
#'
#' @param n_genes Numero di geni da simulare
#' @param k_cell_types Numero di tipi cellulari
#' @param marker_params Parametri dei marker genici
#' @param random_seed Seed per riproducibilità
#' @return Lista con medie di espressione per ogni tipo cellulare
#' @export
generate_baseline_expression <- function(
  n_genes,
  k_cell_types,
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Crea le medie di espressione per ogni cluster
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
    if (k > 1 && !is.null(marker_params$marker_overlap_fold) && marker_params$marker_overlap_fold > 0) {
      prev_start_idx <- (k-2) * marker_params$marker_genes_per_type + 1
      prev_end_idx <- min((k-1) * marker_params$marker_genes_per_type, n_genes)
      
      # Verifica che l'intervallo sia valido prima di procedere
      if (prev_start_idx <= prev_end_idx && prev_start_idx <= n_genes) {
        prev_markers <- prev_start_idx:prev_end_idx
        mu[prev_markers] <- mu[prev_markers] + marker_params$marker_overlap_fold
      }
    }
    
    mean_expression_list[[k]] <- mu
  }
  
  # Restituisci la lista di medie di espressione per cluster
  return(mean_expression_list)
}