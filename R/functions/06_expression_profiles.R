#' Genera profili di espressione genica
#'
#' Crea profili di espressione genica per ogni cellula/spot,
#' incorporando effetti biologici come correlazione spaziale,
#' dropout, dimensioni diverse delle librerie e clustering di geni.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param n_genes Numero di geni da simulare
#' @param k_cell_types Numero di tipi cellulari
#' @param marker_params Parametri dei marker genici
#' @param spatial_params Parametri spaziali
#' @param dropout_params Parametri di dropout
#' @param library_size_params Parametri dimensione libreria
#' @param cell_specific_params Parametri cellula-specifici
#' @param hybrid_params Parametri cellule ibride
#' @param use_spatial_correlation Se TRUE usa correlazione spaziale
#' @param correlation_method Metodo di correlazione ("grf" o "car")
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrice di espressione e metadati associati
#' @importFrom MASS rnbinom
#' @importFrom future.apply future_lapply
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @importFrom scales rescale
#' @export
generate_expression_profiles <- function(
  cell_df,
  n_genes,
  k_cell_types,
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  ),
  spatial_params = list(
    spatial_noise_intensity = 1.0,
    spatial_range = 30,
    random_noise_sd = 0.2,
    gradient_regions = FALSE,
    gradient_width = 5,
    gradient_exponent = 1.5
  ),
  dropout_params = list(
    dropout_range = c(0.2, 0.5),
    dispersion_range = c(2.0, 1.0),
    cell_type_dispersion_effect = 0.2,
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5
  ),
  library_size_params = list(
    mean_library_size = 10000,
    library_size_cv = 0.3,
    spatial_effect_on_library = 0.5,
    cell_type_effect = TRUE
  ),
  cell_specific_params = list(
    cell_specific_noise_sd = 0.2,
    use_gene_modules = TRUE,
    n_gene_modules = 5,
    module_correlation = 0.7
  ),
  hybrid_params = list(
    use_hybrid_cells = TRUE,
    max_hybrid_pairs = 1000,
    hybrid_intensity_range = c(0.2, 0.5)
  ),
  use_spatial_correlation = TRUE,
  correlation_method = "grf",
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # 5a) Medie di espressione per cluster - parametrizzate
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
    if (k > 1 && marker_params$marker_overlap_fold > 0) {
      prev_markers <- ((k-2) * marker_params$marker_genes_per_type + 1):min((k-1) * marker_params$marker_genes_per_type, n_genes)
      if (length(prev_markers) > 0) {
        mu[prev_markers] <- mu[prev_markers] + marker_params$marker_overlap_fold
      }
    }
    
    mean_expression_list[[k]] <- mu
  }
  
  # 5b) Calcolo distanze e densità locale per il modello di dropout
  N <- nrow(cell_df)
  cluster_labels <- cell_df$intensity_cluster
  coords <- cell_df %>% dplyr::select(x, y)
  dist_mat <- as.matrix(dist(coords))
  
  # Calcola la distanza media di ciascuna cellula rispetto alle altre del proprio cluster
  # Versione vettorizzata e ottimizzata
  mean_dist <- numeric(N)
  unique_clusters <- unique(cluster_labels)
  
  # Calcolo più efficiente con chunking ottimizzato
  chunk_size <- max(1, ceiling(N/500))
  chunks <- split(1:N, ceiling(seq_along(1:N)/chunk_size))
  
  # Utilizziamo future_lapply con scheduling migliorato
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
  
  # 5c) Impostazione parametri di dispersione in base ai parametri
  if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df)) {
    # Se utilizziamo gradienti, facciamo variare la dispersione in base alla distanza dal confine
    dispersion_param <- dropout_params$dispersion_range[2] +
      cell_df$boundary_dist * (dropout_params$dispersion_range[1] - dropout_params$dispersion_range[2])
  } else {
    # Altrimenti usiamo il metodo originale basato sulla distanza media
    dispersion_param <- scales::rescale(mean_dist, to = dropout_params$dispersion_range)
  }
  
  # Aggiungi effetto del tipo cellulare sulla dispersione
  set.seed(random_seed + 4)
  type_dispersion_effects <- runif(k_cell_types,
                                 min = 1 - dropout_params$cell_type_dispersion_effect,
                                 max = 1 + dropout_params$cell_type_dispersion_effect)
  
  # Applica effetto moltiplicativo per tipo cellulare
  for (k in 1:k_cell_types) {
    dispersion_param[cluster_labels == k] <- dispersion_param[cluster_labels == k] * type_dispersion_effects[k]
  }
  
  # 5d) Simulazione delle dimensioni delle librerie
  set.seed(random_seed + 1)
  
  # Genera dimensioni libreria con distribuzione log-normale
  library_size_sd <- library_size_params$mean_library_size * library_size_params$library_size_cv
  log_mean <- log(library_size_params$mean_library_size^2 /
                 sqrt(library_size_params$mean_library_size^2 + library_size_sd^2))
  log_sd <- sqrt(log(1 + (library_size_sd^2 / library_size_params$mean_library_size^2)))
  
  library_size <- rlnorm(N, meanlog = log_mean, sdlog = log_sd)
  
  # Aggiungi effetto spaziale sulla dimensione libreria se richiesto
  if (library_size_params$spatial_effect_on_library > 0) {
    # Converti cell_df in oggetto spatial per il GP
    sp_df_lib <- cell_df
    coordinates(sp_df_lib) <- ~ x + y
    
    # Crea un GP per l'effetto spaziale sulla dimensione libreria
    lib_gp <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                   beta = 0, model = vgm(psill = 1.0,
                                        range = spatial_params$spatial_range * 1.5,
                                        model = "Exp"),
                   nmax = 15)
    
    # Genera il noise spaziale
    set.seed(random_seed + 2)
    lib_noise <- predict(lib_gp, newdata = sp_df_lib, nsim = 1)$sim1
    
    # Normalizza e scala il noise
    lib_noise <- scale(lib_noise)
    lib_effect <- library_size_params$spatial_effect_on_library * lib_noise
    
    # Applica alla dimensione libreria (effetto moltiplicativo)
    library_size <- library_size * exp(lib_effect)
  }
  
  # Aggiungi effetto del tipo cellulare sulla dimensione della libreria
  if (library_size_params$cell_type_effect) {
    # Diversi tipi cellulari hanno diversi contenuti di RNA
    cell_type_effect <- numeric(N)
    
    # Crea effetti diversi per diversi tipi cellulari
    set.seed(random_seed + 3)
    type_effects <- rnorm(k_cell_types, mean = 0, sd = 0.2)  # Effetti casuali per tipo
    
    # Assegna effetto in base al tipo cellulare
    for (k in 1:k_cell_types) {
      cell_type_effect[cluster_labels == k] <- type_effects[k]
    }
    
    # Applica effetto moltiplicativo
    library_size <- library_size * exp(cell_type_effect)
  }
  
  # 5e) Impostazione probabilità di dropout
  if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df)) {
    # Più dropout vicino al confine
    base_dropout <- dropout_params$dropout_range[1] +
      (1 - cell_df$boundary_dist) * (dropout_params$dropout_range[2] - dropout_params$dropout_range[1])
  } else {
    # Metodo originale
    base_dropout <- scales::rescale(mean_dist, to = dropout_params$dropout_range)
  }
  
  # 5f) Generazione espressione usando Negative Binomial
  expression_data <- matrix(0, nrow = N, ncol = n_genes)
  
  # Identificazione geni stabili (sub-Poissoniani)
  stable_genes <- sample(n_genes, max(1, round(n_genes * 0.1)))  # 10% geni stabili
  
  # Crea moduli di geni co-espressi
  gene_modules <- NULL
  module_noise <- NULL
  
  if (cell_specific_params$use_gene_modules) {
    # Calcola il numero di geni per modulo
    genes_per_module <- ceiling(n_genes / cell_specific_params$n_gene_modules)
    
    # Assegna geni ai moduli
    gene_modules <- list()
    for (m in 1:cell_specific_params$n_gene_modules) {
      start_idx <- (m-1) * genes_per_module + 1
      end_idx <- min(m * genes_per_module, n_genes)
      gene_modules[[m]] <- start_idx:end_idx
    }
    
    # Crea rumore correlato per ogni modulo
    module_noise <- matrix(0, nrow = N, ncol = n_genes)
    
    # Genera rumore base per ogni modulo
    base_module_noise <- matrix(rnorm(cell_specific_params$n_gene_modules * N),
                               nrow = N, ncol = cell_specific_params$n_gene_modules)
    
    # Applica il rumore del modulo a ciascun gene appartenente al modulo
    for (m in 1:length(gene_modules)) {
      module_genes <- gene_modules[[m]]
      # Assegna lo stesso rumore base a tutti i geni del modulo, scalato per la correlazione
      module_noise[, module_genes] <- base_module_noise[, m] * cell_specific_params$module_correlation
    }
  }
  
  # Simulazione della correlazione spaziale continua - versione semplificata per il pacchetto
  gp_noise <- NULL
  if (use_spatial_correlation) {
    # Crea oggetto spatial points
    sp_df <- cell_df
    coordinates(sp_df) <- ~ x + y
    
    # Parametri per gstat
    gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                   beta = 0, model = vgm(psill = spatial_params$spatial_noise_intensity * 2,
                                        range = spatial_params$spatial_range * 0.8,
                                        model = "Exp"),
                   nmax = 20)
    
    # Genera il noise spaziale
    set.seed(random_seed)
    gp_noise <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
    
    # Normalizza il noise
    gp_noise <- scale(gp_noise) * 1.5
  }
  
  # Crea cellule ibride ai confini tra cluster, se richiesto
  hybrid_matrix <- matrix(0, nrow = N, ncol = k_cell_types)
  
  if (hybrid_params$use_hybrid_cells) {
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
      set.seed(random_seed + 1)
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
  
  # Pre-calcola alcune strutture di dati comuni a tutti i geni
  # Crea una matrice di medie di espressione per tipo di cellula e gene
  all_mean_expr <- matrix(0, nrow = n_genes, ncol = k_cell_types)
  for (g in seq_len(n_genes)) {
    for (k in 1:k_cell_types) {
      all_mean_expr[g, k] <- mean_expression_list[[k]][g]
    }
  }
  
  # Converti cluster_labels in interi una sola volta
  cl <- as.integer(cluster_labels)
  
  # Genera l'espressione genica in chunk - versione semplificata
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
      if (hybrid_params$use_hybrid_cells) {
        # Approccio con cellule ibride - vettorizzato
        # Calcola l'effetto di tutti i cluster contemporaneamente
        
        # Crea una matrice dove ogni riga è la media di espressione per ogni tipo di cellula
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
      
      # Applica dropout in base al modello specificato - versione vettorizzata
      if (dropout_params$expression_dependent_dropout) {
        # Definisci una funzione vettorizzata per normalizzare tra 0 e 1
        scale01_vec <- function(x) {
          if (all(x == x[1])) return(rep(0.5, length(x)))
          (x - min(x)) / (max(x) - min(x))
        }
        
        # Normalizza l'espressione del gene corrente
        norm_expr <- scale01_vec(chunk_expression[, i])
        
        # Calcola la probabilità di dropout con una funzione logistica vettorizzata
        dropout_prob_expr <- 1 / (1 + exp((norm_expr - dropout_params$dropout_curve_midpoint) *
                                       dropout_params$dropout_curve_steepness))
        
        # Combina con il dropout spaziale base (media pesata) - operazione vettorizzata
        dropout_prob <- 0.7 * dropout_prob_expr + 0.3 * base_dropout
        
        # Tronca i valori al range [0,1] in un'unica operazione
        dropout_prob <- pmin(pmax(dropout_prob, 0), 1)
        
        # Applica dropout in modo vettorizzato
        zero_idx <- runif(N) < dropout_prob
        chunk_expression[zero_idx, i] <- 0
      } else {
        # Modello di dropout originale (solo spaziale) - vettorizzato
        zero_idx <- runif(N) < base_dropout
        chunk_expression[zero_idx, i] <- 0
      }
    } # Fine del ciclo for sui geni di questo chunk
    
    return(chunk_expression)
  }, future.scheduling = 1, future.seed = TRUE) # Fine del future_lapply
  
  # Combina i risultati dei chunk in una singola matrice di espressione
  expression_data <- matrix(0, nrow = N, ncol = n_genes)
  for (i in seq_along(gene_chunks)) {
    genes_subset <- gene_chunks[[i]]
    expression_data[, genes_subset] <- expression_chunks[[i]]
  }
  
  # Prepara l'output
  result <- list(
    expression = expression_data,
    library_size = library_size,
    dispersion_param = dispersion_param,
    gene_modules = gene_modules,
    mean_expression_list = mean_expression_list
  )
  
  return(result)
}