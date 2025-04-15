#' Genera la matrice finale di espressione genica
#'
#' Combina tutti i componenti per generare la matrice finale di espressione
#' con tutti gli effetti biologici e tecnici.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param mean_expression_list Lista con le medie di espressione per tipo cellulare
#' @param n_genes Numero di geni
#' @param library_size Vettore con dimensioni di libreria
#' @param dispersion_param Vettore con parametri di dispersione
#' @param base_dropout Vettore con probabilità base di dropout
#' @param hybrid_matrix Matrice di ibridazione
#' @param module_noise Matrice con rumore correlato dei moduli
#' @param gp_noise Vettore con rumore di correlazione spaziale
#' @param latent_factors Fattori latenti che influenzano i moduli di espressione
#' @param module_network Rete di interazioni tra i diversi moduli di geni
#' @param spatial_params Parametri spaziali
#' @param dropout_params Parametri di dropout
#' @param cell_specific_params Parametri cellula-specifici
#' @param use_spatial_correlation Se usare la correlazione spaziale
#' @param random_seed Seed per riproducibilità
#' @return Matrice di espressione finale
#' @importFrom MASS rnbinom
#' @importFrom future.apply future_lapply
#' @export
generate_expression_matrix <- function(
  cell_df,
  mean_expression_list,
  n_genes,
  library_size,
  dispersion_param,
  base_dropout,
  hybrid_matrix,
  module_noise,
  gp_noise,
  latent_factors = NULL,
  module_network = NULL,
  spatial_params = list(
    spatial_noise_intensity = 1.0,
    random_noise_sd = 0.2
  ),
  dropout_params = list(
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5
  ),
  cell_specific_params = list(
    cell_specific_noise_sd = 0.2,
    use_gene_modules = TRUE,
    module_hierarchical = FALSE,
    module_overlap = 0.1,
    module_size_distribution = "exponential",
    n_latent_factors = 3,
    module_network_density = 0.2,
    latent_factor_strength = 0.8
  ),
  use_spatial_correlation = TRUE,
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai variabili
  N <- nrow(cell_df)
  cluster_labels <- cell_df$intensity_cluster
  k_cell_types <- length(unique(cluster_labels))
  
  # Identificazione geni stabili (sub-Poissoniani)
  stable_genes <- sample(n_genes, max(1, round(n_genes * 0.1)))  # 10% geni stabili
  
  # Pre-calcola la matrice di medie di espressione per tipo di cellula e gene
  all_mean_expr <- matrix(0, nrow = n_genes, ncol = k_cell_types)
  for (g in seq_len(n_genes)) {
    for (k in 1:k_cell_types) {
      all_mean_expr[g, k] <- mean_expression_list[[k]][g]
    }
  }
  
  # Converti cluster_labels in interi una sola volta
  cl <- as.integer(cluster_labels)
  
  # Genera l'espressione genica in chunk
  chunk_size <- max(5, ceiling(n_genes/32))
  gene_chunks <- split(seq_len(n_genes), ceiling(seq_len(n_genes)/chunk_size))
  
  expression_chunks <- future_lapply(gene_chunks, function(genes_subset) {
    # Alloca lo storage per l'espressione di questo chunk
    chunk_expression <- matrix(0, nrow = N, ncol = length(genes_subset))
    
    # Genera rumore cellula-specifico una volta sola per tutto il chunk
    all_cell_specific_effects <- matrix(
      rnorm(N * length(genes_subset), 0, cell_specific_params$cell_specific_noise_sd),
      nrow = N, ncol = length(genes_subset)
    )
    
    for (i in seq_along(genes_subset)) {
      g <- genes_subset[i]
      
      # Usa il rumore pre-generato
      cell_specific_effect <- all_cell_specific_effects[, i]
      
      # Calcola medie di espressione di base - versione molto più veloce usando indexing
      base_expr <- all_mean_expr[g, cl]
      
      # Applica effetto di ibridazione
      if (!is.null(hybrid_matrix) && sum(hybrid_matrix) > 0) {
        # Calcola l'effetto di tutti i cluster contemporaneamente
        all_cluster_expr <- sapply(1:k_cell_types, function(k) mean_expression_list[[k]][g])
        
        # Calcola l'effetto ibrido complessivo
        hybrid_effect <- hybrid_matrix %*% all_cluster_expr
        
        # Calcola il peso complessivo dell'effetto ibrido su ogni cellula
        hybrid_weight <- rowSums(hybrid_matrix)
        
        # Applica solo a cellule che hanno un effetto ibrido
        hybrid_cells <- which(hybrid_weight > 0)
        if (length(hybrid_cells) > 0) {
          # Effetto ibrido totale per ogni cellula
          base_expr[hybrid_cells] <- base_expr[hybrid_cells] * (1 - hybrid_weight[hybrid_cells]) +
                                   hybrid_effect[hybrid_cells]
        }
      }
      
      # Combina con l'effetto cellula-specifico
      mu_vals <- base_expr + cell_specific_effect
      
      # Aggiungi effetto dei moduli genici se abilitato
      if (cell_specific_params$use_gene_modules && !is.null(module_noise)) {
        # Aggiungi il rumore correlato del modulo genico a cui appartiene questo gene
        mu_vals <- mu_vals + module_noise[, g]
      }
      
      # Aggiungi correlazione spaziale se richiesta
      if (use_spatial_correlation && !is.null(gp_noise)) {
        mu_vals <- mu_vals + spatial_params$spatial_noise_intensity * gp_noise
        
        # Aggiungi noise casuale addizionale per confondere i pattern
        random_noise <- rnorm(length(mu_vals), 0, spatial_params$random_noise_sd)
        mu_vals <- mu_vals + random_noise
      }
      
      # Genera conteggi di espressione
      if (g %in% stable_genes) {
        # Modello sub-Poisson: Binomiale con p alto e n moderato
        p <- 0.9
        n_trial <- round(exp(mu_vals)/(1-p))
        raw_counts <- rbinom(N, n_trial, p)
      } else {
        # Negative Binomial con dispersione variabile spazialmente
        raw_counts <- rnbinom(N, mu = exp(mu_vals), size = dispersion_param)
      }
      
      # Applica l'effetto della dimensione della libreria
      scaled_counts <- raw_counts * (library_size / mean(library_size))
      # Arrotonda a numeri interi (conteggi)
      chunk_expression[, i] <- round(scaled_counts)
      
      # Applica dropout in base al modello specificato
      # Definisci una funzione vettorizzata per normalizzare tra 0 e 1
      scale01_vec <- function(x) {
        if (all(x == x[1])) return(rep(0.5, length(x)))
        (x - min(x)) / (max(x) - min(x))
      }
      
      if (dropout_params$expression_dependent_dropout) {
        # Normalizza l'espressione del gene corrente
        norm_expr <- scale01_vec(chunk_expression[, i])
        
        # Calcola la probabilità di dropout con una funzione logistica
        dropout_prob_expr <- 1 / (1 + exp((norm_expr - dropout_params$dropout_curve_midpoint) *
                                     dropout_params$dropout_curve_steepness))
        
        # Combina con il dropout spaziale base
        dropout_prob <- 0.7 * dropout_prob_expr + 0.3 * base_dropout
        
        # Tronca i valori al range [0,1]
        dropout_prob <- pmin(pmax(dropout_prob, 0), 1)
        
        # Applica dropout
        zero_idx <- runif(N) < dropout_prob
        chunk_expression[zero_idx, i] <- 0
      } else {
        # Modello di dropout originale (solo spaziale)
        zero_idx <- runif(N) < base_dropout
        chunk_expression[zero_idx, i] <- 0
      }
    }
    
    return(chunk_expression)
  }, future.scheduling = 1, future.seed = TRUE)
  
  # Combina i risultati dei chunk in una singola matrice di espressione
  expression_data <- matrix(0, nrow = N, ncol = n_genes)
  for (i in seq_along(gene_chunks)) {
    genes_subset <- gene_chunks[[i]]
    expression_data[, genes_subset] <- expression_chunks[[i]]
  }
  
  return(expression_data)
}