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
    # Imposta seed specifico per tipo cellulare per riproducibilità
    set.seed(random_seed + k)
    
    # Genera distribuzione biologicamente realistica dell'espressione genica
    # Adattata per simulazioni full-size (20k geni)
    if (n_genes > 10000) {
      # Distribuzione per simulazioni full-size
      # 60% geni non espressi, 30% lowly, 8% medium, 2% highly expressed
      n_zero <- round(n_genes * 0.60)
      n_low <- round(n_genes * 0.30)
      n_med <- round(n_genes * 0.08)
      n_high <- n_genes - n_zero - n_low - n_med
      
      # Genera valori mu per ciascuna categoria
      mu_zero <- rep(-20, n_zero)                      # Praticamente zero
      mu_low <- rnorm(n_low, mean = -3, sd = 0.5)      # Bassa espressione
      mu_med <- rnorm(n_med, mean = 0, sd = 0.4)       # Media espressione  
      mu_high <- rnorm(n_high, mean = 2, sd = 0.3)     # Alta espressione (housekeeping)
      
      # Applica bounds
      mu_low <- pmax(mu_low, -7)
      mu_low <- pmin(mu_low, -2)
      mu_med <- pmax(mu_med, -3)
      mu_med <- pmin(mu_med, 0.5)
      mu_high <- pmax(mu_high, 0)
      mu_high <- pmin(mu_high, 3)
      
    } else {
      # Distribuzione originale per simulazioni semplificate (2k geni)
      # 85% geni lowly expressed, 10% medium, 5% highly expressed
      n_zero <- 0  # Nessun gene completamente non espresso
      n_low <- round(n_genes * 0.85)
      n_med <- round(n_genes * 0.10)
      n_high <- n_genes - n_low - n_med
      
      mu_zero <- numeric(0)  # Array vuoto
      mu_low <- rnorm(n_low, mean = -4.5, sd = 0.4)
      mu_med <- rnorm(n_med, mean = -1.5, sd = 0.4)
      mu_high <- rnorm(n_high, mean = 0.5, sd = 0.3)
      
      # Applica bounds originali
      mu_low <- pmax(mu_low, -7)
      mu_med <- pmax(mu_med, -4)
      mu_med <- pmin(mu_med, -1)
      mu_high <- pmax(mu_high, -1)
      mu_high <- pmin(mu_high, 1.5)
    }
    
    # Combina e mescola per distribuzione casuale
    mu <- c(mu_zero, mu_low, mu_med, mu_high)
    mu <- sample(mu)  # randomizza l'ordine dei geni
    
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