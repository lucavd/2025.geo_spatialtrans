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
    # 85% geni lowly expressed, 10% medium, 5% highly expressed
    n_low <- round(n_genes * 0.85)
    n_med <- round(n_genes * 0.10)
    n_high <- n_genes - n_low - n_med
    
    # Genera valori mu per ciascuna categoria
    # Ultimo fine-tuning per aumentare mediana espressione mantenendo UMI realistici
    mu_low <- rnorm(n_low, mean = -4.5, sd = 0.4)   # exp(-4.5) ≈ 0.011
    mu_med <- rnorm(n_med, mean = -1.5, sd = 0.4)   # exp(-1.5) ≈ 0.22
    mu_high <- rnorm(n_high, mean = 0.5, sd = 0.3)  # exp(0.5) ≈ 1.65
    
    # Applica bounds per evitare valori estremi
    mu_low <- pmax(mu_low, -7)    # Floor per low expressed
    mu_med <- pmax(mu_med, -4)    # Floor per medium
    mu_med <- pmin(mu_med, -1)    # Ceiling per medium
    mu_high <- pmax(mu_high, -1)  # Floor per high
    mu_high <- pmin(mu_high, 1.5) # Ceiling per high (evita outliers)
    
    # Combina e mescola per distribuzione casuale
    mu <- c(mu_low, mu_med, mu_high)
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