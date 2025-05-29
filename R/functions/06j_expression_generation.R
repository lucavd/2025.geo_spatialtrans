#' Genera la matrice finale di espressione genica
#'
#' Combina tutti i componenti per generare la matrice finale di espressione
#' con tutti gli effetti biologici e tecnici, includendo ambient RNA contamination 
#' e dropout gene-specifico.
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
#' @param ambient_params Parametri per RNA ambientale
#' @param cell_specific_params Parametri cellula-specifici
#' @param use_spatial_correlation Se usare la correlazione spaziale
#' @param dist_mat Matrice di distanza tra celle (per ambient RNA)
#' @param random_seed Seed per riproducibilità
#' @return Matrice di espressione finale
#' @importFrom MASS rnbinom
#' @importFrom future.apply future_lapply
#' @importFrom stats rbinom rnorm runif
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
    dropout_curve_steepness = 5,
    use_gene_specific_dropout = TRUE,
    gene_dropout_variability = 0.3,
    gc_content_effect = 0.5,
    length_effect = 0.3,
    sequence_effect = 0.4,
    gene_effect_weight = 0.3
  ),
  ambient_params = list(
    use_ambient_rna = FALSE,
    ambient_contamination_rate = 0.05,
    ambient_diffusion_distance = 30,
    tissue_leakage_factor = 0.7,
    background_noise = 0.1
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
  dist_mat = NULL,
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai variabili
  N <- nrow(cell_df)
  # Usa i livelli originali per mappare le medie di espressione
  cluster_labels <- cell_df$intensity_cluster
  # Il numero di tipi cellulari corrisponde alla lunghezza di mean_expression_list
  k_cell_types <- length(mean_expression_list)
  
  # Identificazione geni stabili (sub-Poissoniani)
  stable_genes <- sample(n_genes, max(1, round(n_genes * 0.1)))  # 10% geni stabili
  
  # Per simulazioni full-size, conta i geni effettivamente espressi
  if (n_genes > 10000) {
    # Trova geni con espressione > -10 in almeno un tipo cellulare
    expressed_genes_mask <- sapply(1:n_genes, function(g) {
      any(sapply(mean_expression_list, function(ml) ml[g] > -10))
    })
    n_expressed_genes <- sum(expressed_genes_mask)
    cat("DEBUG - Geni effettivamente espressi:", n_expressed_genes, "su", n_genes, "\n")
  } else {
    # Per simulazioni piccole, tutti i geni sono considerati espressi
    n_expressed_genes <- n_genes
  }
  
  # Pre-calcola la matrice di medie di espressione per tipo di cellula e gene
  all_mean_expr <- matrix(0, nrow = n_genes, ncol = k_cell_types)
  for (g in seq_len(n_genes)) {
    for (k in 1:k_cell_types) {
      all_mean_expr[g, k] <- mean_expression_list[[k]][g]
    }
  }
  
  # Converti cluster_labels in interi una sola volta
  cl <- as.integer(cluster_labels)
  
  # Genera fattori di dropout gene-specifici se richiesto
  use_gene_specific_dropout <- ifelse(is.null(dropout_params$use_gene_specific_dropout), 
                                    TRUE, dropout_params$use_gene_specific_dropout)
  
  gene_specific_factors <- NULL
  if (use_gene_specific_dropout) {
    gene_specific_factors <- generate_gene_specific_dropout_factors(
      n_genes, N, dropout_params, random_seed
    )
  }
  
  # Genera l'espressione genica in chunk
  chunk_size <- max(5, ceiling(n_genes/32))
  gene_chunks <- split(seq_len(n_genes), ceiling(seq_len(n_genes)/chunk_size))
  
  # Usa parallelizzazione con future_lapply se disponibile
  expression_chunks <- lapply(gene_chunks, function(genes_subset) {
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
        intensity <- ifelse(!is.null(spatial_params$spatial_noise_intensity), 
                           spatial_params$spatial_noise_intensity, 1.0)
        mu_vals <- mu_vals + intensity * gp_noise
        
        # Aggiungi noise casuale addizionale per confondere i pattern
        random_sd <- ifelse(!is.null(spatial_params$random_noise_sd), 
                           spatial_params$random_noise_sd, 0.2)
        random_noise <- rnorm(length(mu_vals), 0, random_sd)
        mu_vals <- mu_vals + random_noise
      }
      
      # Genera conteggi di espressione
      # Usa un approccio che mantiene scale biologiche realistiche
      
      if (g %in% stable_genes) {
        # Modello sub-Poisson: Binomiale con p alto e n moderato
        p <- 0.9
        # Scala exp(mu) per frazione di library size dedicata a questo gene
        lambda <- exp(mu_vals) * library_size / n_expressed_genes
        # NUOVO: Limita lambda per evitare valori estremi
        lambda <- pmin(lambda, 5000)  # Max 5000 UMI per gene (realistico per Visium HD)
        n_trial <- round(lambda/(1-p))
        raw_counts <- rbinom(N, n_trial, p)
      } else {
        # Negative Binomial con dispersione variabile spazialmente
        # Scala exp(mu) per frazione di library size dedicata a questo gene
        lambda <- exp(mu_vals) * library_size / n_expressed_genes
        # NUOVO: Limita lambda per evitare valori estremi
        lambda <- pmin(lambda, 5000)  # Max 5000 UMI per gene (realistico per Visium HD)
        raw_counts <- rnbinom(N, mu = lambda, size = dispersion_param)
      }
      
      # I conteggi sono già scalati per library size
      scaled_counts <- raw_counts
      # NUOVO: Applica cap finale per sicurezza
      scaled_counts <- pmin(scaled_counts, 10000)  # Assolutamente nessun gene > 10k UMI
      # Arrotonda a numeri interi (conteggi)
      chunk_expression[, i] <- round(scaled_counts)
      
      # Estrai fattori gene-specifici per il dropout se disponibili
      gene_factors <- NULL
      if (use_gene_specific_dropout && !is.null(gene_specific_factors)) {
        gene_factors <- list(
          gene_dropout_factors = gene_specific_factors$gene_dropout_factors[g],
          dropout_baseline_shift = gene_specific_factors$dropout_baseline_shift[g]
        )
      }
      
      # Applica dropout usando la versione migliorata con supporto gene-specifico
      chunk_expression[, i] <- apply_dropout(
        chunk_expression[, i], 
        base_dropout, 
        gene_factors,
        dropout_params, 
        random_seed + g  # Usa un seed diverso per ogni gene
      )
    }
    
    return(chunk_expression)
  })
  
  # Combina i risultati dei chunk in una singola matrice di espressione
  expression_data <- matrix(0, nrow = N, ncol = n_genes)
  for (i in seq_along(gene_chunks)) {
    genes_subset <- gene_chunks[[i]]
    expression_data[, genes_subset] <- expression_chunks[[i]]
  }
  
  # Aggiungi contaminazione da RNA ambientale se richiesto
  use_ambient_rna <- ifelse(is.null(ambient_params$use_ambient_rna), 
                          FALSE, ambient_params$use_ambient_rna)
  
  if (use_ambient_rna) {
    # Genera la contaminazione da RNA ambientale
    ambient_contamination <- generate_ambient_rna(
      cell_df, dist_mat, mean_expression_list, ambient_params, random_seed
    )
    
    # Applica la contaminazione alla matrice di espressione
    expression_data <- apply_ambient_rna_contamination(
      expression_data, ambient_contamination, ambient_params
    )
  }
  
  # NUOVO: Validazione biologica finale
  # Controlla e corregge eventuali valori biologicamente irrealistici
  expression_data <- validate_biological_plausibility(expression_data)
  
  return(expression_data)
}

#' Valida e corregge la plausibilità biologica dell'espressione
#'
#' @param expression_matrix Matrice di espressione sparsa
#' @return Matrice corretta
#' @noRd
validate_biological_plausibility <- function(expression_matrix) {
  # Parametri biologici realistici
  max_umi_per_cell <- 50000  # Massimo per Visium HD
  max_umi_per_gene_per_cell <- 5000  # Nessun gene dovrebbe superare questo
  
  # Converti in matrice sparsa se necessario
  if (!inherits(expression_matrix, "sparseMatrix")) {
    expression_matrix <- Matrix(expression_matrix, sparse = TRUE)
  }
  
  # 1. Controlla UMI totali per cella
  cell_totals <- Matrix::colSums(expression_matrix)
  cells_over_limit <- which(cell_totals > max_umi_per_cell)
  
  if (length(cells_over_limit) > 0) {
    cat(sprintf("Correzione biologica: %d celle con >%d UMI\n", 
                length(cells_over_limit), max_umi_per_cell))
    
    # Processa ogni cella che supera il limite
    for (i in seq_along(cells_over_limit)) {
      cell <- cells_over_limit[i]
      scale_factor <- max_umi_per_cell / cell_totals[cell]
      
      # Ottieni gli indici dei valori non-zero per questa cella
      cell_indices <- which(expression_matrix[, cell] != 0)
      
      # Scala i valori
      for (gene in cell_indices) {
        new_val <- round(expression_matrix[gene, cell] * scale_factor)
        expression_matrix[gene, cell] <- max(new_val, 1)  # Mantieni almeno 1
      }
    }
    
    # Verifica che il cap sia stato applicato
    new_totals <- Matrix::colSums(expression_matrix[, cells_over_limit])
    if (any(new_totals > max_umi_per_cell * 1.1)) {
      cat("ATTENZIONE: Alcune celle ancora sopra il limite dopo correzione\n")
    }
  }
  
  # 2. Controlla valori estremi per gene
  extreme_values <- which(expression_matrix@x > max_umi_per_gene_per_cell)
  if (length(extreme_values) > 0) {
    cat(sprintf("Correzione biologica: %d valori >%d UMI per gene\n", 
                length(extreme_values), max_umi_per_gene_per_cell))
    expression_matrix@x[extreme_values] <- max_umi_per_gene_per_cell
  }
  
  # 3. Rimuovi valori < 1 dopo le correzioni
  expression_matrix@x[expression_matrix@x < 1] <- 0
  expression_matrix <- Matrix::drop0(expression_matrix)
  
  # Report finale
  final_cell_totals <- Matrix::colSums(expression_matrix)
  cat("\nStatistiche post-validazione biologica:\n")
  cat(sprintf("- Max UMI per valore: %d\n", max(expression_matrix@x)))
  cat(sprintf("- Max UMI per cella: %d\n", max(final_cell_totals)))
  cat(sprintf("- Mediana UMI per cella: %.0f\n", median(final_cell_totals)))
  cat(sprintf("- Media espressione (non-zero): %.2f\n", mean(expression_matrix@x)))
  cat(sprintf("- Celle con >%dk UMI: %d\n", max_umi_per_cell/1000, sum(final_cell_totals > max_umi_per_cell)))
  
  return(expression_matrix)
}