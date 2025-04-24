#' Modellazione semplificata di splicing alternativo per testing
#'
#' Versioni semplificate delle funzioni di splicing alternativo adattate per testing
#' 
#' @import stats
#' @import sp
#' @import gstat

#' Genera varianti di splicing per i geni (versione semplificata per testing)
#'
#' @param n_genes Numero totale di geni
#' @param splicing_params Parametri per lo splicing alternativo
#' @param random_seed Seed per riproducibilità
#' @return Lista con definizione delle varianti di splicing
#' @export
generate_splicing_variants <- function(
  n_genes,
  splicing_params = list(
    splicing_fraction = 0.6,
    n_splicing_variants = 2
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  splicing_fraction <- splicing_params$splicing_fraction
  n_variants <- splicing_params$n_splicing_variants
  
  # Valori predefiniti se non specificati
  if (is.null(splicing_fraction)) splicing_fraction <- 0.6
  if (is.null(n_variants)) n_variants <- 2
  
  # Seleziona geni con varianti di splicing
  n_genes_with_variants <- round(n_genes * splicing_fraction)
  genes_with_variants <- sort(sample(1:n_genes, n_genes_with_variants))
  
  # Crea vettore di conteggio varianti
  variant_counts <- rep(1, n_genes)
  variant_counts[genes_with_variants] <- n_variants
  
  # Restituisci risultati
  return(list(
    genes_with_variants = genes_with_variants,
    variant_counts = variant_counts
  ))
}

#' Calcola propensione di splicing (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param genes_with_variants Vettore di ID geni con varianti di splicing
#' @param n_variants Numero di varianti per gene
#' @param splicing_params Parametri per il calcolo delle propensioni
#' @param random_seed Seed per riproducibilità
#' @return Array 3D di propensioni di splicing [cellule, geni con varianti, varianti]
#' @export
calculate_splicing_propensity <- function(
  cell_df,
  genes_with_variants,
  n_variants,
  splicing_params = list(
    splicing_spatial_pattern = "gradient",
    splicing_cluster_specific = FALSE,
    n_splicing_foci = 3
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  spatial_pattern <- splicing_params$splicing_spatial_pattern
  cluster_specific <- splicing_params$splicing_cluster_specific
  n_foci <- splicing_params$n_splicing_foci
  
  # Valori predefiniti se non specificati
  if (is.null(spatial_pattern)) spatial_pattern <- "gradient"
  if (is.null(cluster_specific)) cluster_specific <- FALSE
  if (is.null(n_foci)) n_foci <- 3
  
  # Dimensioni
  n_cells <- nrow(cell_df)
  n_genes_with_variants <- length(genes_with_variants)
  
  # Crea array di propensioni
  propensity <- array(0, dim = c(n_cells, n_genes_with_variants, n_variants))
  
  # Crea pattern spaziale in base al tipo specificato
  if (spatial_pattern == "gradient") {
    # Gradiente dall'angolo in alto a sinistra
    x_norm <- (cell_df$x - min(cell_df$x)) / (max(cell_df$x) - min(cell_df$x))
    y_norm <- (cell_df$y - min(cell_df$y)) / (max(cell_df$y) - min(cell_df$y))
    
    # Per ogni gene con varianti
    for (i in 1:n_genes_with_variants) {
      # Angolo casuale per il gradiente
      angle <- runif(1, 0, 2*pi)
      gradient <- x_norm * cos(angle) + y_norm * sin(angle)
      
      # Trasforma in propensioni per tutte le varianti
      for (j in 1:n_variants) {
        if (j == 1) {
          propensity[, i, j] <- gradient
        } else {
          propensity[, i, j] <- 1 - gradient
        }
      }
    }
  } else if (spatial_pattern == "focal") {
    # Pattern con foci puntuali
    for (i in 1:n_genes_with_variants) {
      # Crea foci casuali
      n_gene_foci <- sample(1:n_foci, 1)
      focal_x <- runif(n_gene_foci, min(cell_df$x), max(cell_df$x))
      focal_y <- runif(n_gene_foci, min(cell_df$y), max(cell_df$y))
      
      # Per ogni cellula, calcola distanza dal focus più vicino
      focal_pattern <- rep(1, n_cells)
      for (c in 1:n_cells) {
        min_dist <- Inf
        for (f in 1:n_gene_foci) {
          dist <- sqrt((cell_df$x[c] - focal_x[f])^2 + (cell_df$y[c] - focal_y[f])^2)
          min_dist <- min(min_dist, dist)
        }
        focal_pattern[c] <- exp(-min_dist / (0.2 * (max(cell_df$x) - min(cell_df$x))))
      }
      
      # Normalizza il pattern
      focal_pattern <- (focal_pattern - min(focal_pattern)) / (max(focal_pattern) - min(focal_pattern))
      
      # Trasforma in propensioni per tutte le varianti
      for (j in 1:n_variants) {
        if (j == 1) {
          propensity[, i, j] <- focal_pattern
        } else {
          propensity[, i, j] <- 1 - focal_pattern
        }
      }
    }
  } else {
    # Pattern casuale
    for (i in 1:n_genes_with_variants) {
      for (j in 1:n_variants) {
        propensity[, i, j] <- runif(n_cells)
      }
    }
  }
  
  # Applica specifici pattern per cluster se richiesto
  if (cluster_specific && "cluster" %in% colnames(cell_df)) {
    clusters <- unique(cell_df$cluster)
    n_clusters <- length(clusters)
    
    # Per ogni cluster, modifica le propensioni
    for (cl in 1:n_clusters) {
      cluster_id <- clusters[cl]
      cells_in_cluster <- which(cell_df$cluster == cluster_id)
      
      if (length(cells_in_cluster) > 0) {
        # Per ogni gene con varianti, genera preferenze di cluster
        for (i in 1:n_genes_with_variants) {
          # Preferenza casuale per varianti in questo cluster
          cluster_prefs <- runif(n_variants)
          cluster_prefs <- cluster_prefs / sum(cluster_prefs)
          
          # Modifica propensioni delle cellule in questo cluster
          for (j in 1:n_variants) {
            propensity[cells_in_cluster, i, j] <- propensity[cells_in_cluster, i, j] * cluster_prefs[j] * 2
          }
        }
      }
    }
  }
  
  # Normalizza le propensioni per ogni cellula-gene in modo che sommino a 1 tra le varianti
  for (c in 1:n_cells) {
    for (g in 1:n_genes_with_variants) {
      sum_prop <- sum(propensity[c, g, ])
      if (sum_prop > 0) {
        propensity[c, g, ] <- propensity[c, g, ] / sum_prop
      } else {
        # Se tutte le propensioni sono 0, assegna pari probabilità
        propensity[c, g, ] <- rep(1/n_variants, n_variants)
      }
    }
  }
  
  return(propensity)
}

#' Genera espressione di varianti di splicing (versione semplificata per testing)
#'
#' @param expr_matrix Matrice di espressione originale
#' @param genes_with_variants Vettore di ID geni con varianti di splicing
#' @param variant_counts Vettore di conteggio varianti per gene
#' @param propensity Array 3D di propensioni di splicing
#' @param splicing_params Parametri per generazione di espressione
#' @return Lista con matrice di espressione modificata e matrici di varianti
#' @export
generate_splicing_expression <- function(
  expr_matrix,
  genes_with_variants,
  variant_counts,
  propensity,
  splicing_params = list(
    splicing_strength = 1.0,
    variant_expression_ratios = c(1.0, 0.8)
  )
) {
  # Estrai parametri
  splicing_strength <- splicing_params$splicing_strength
  variant_ratios <- splicing_params$variant_expression_ratios
  
  # Valori predefiniti se non specificati
  if (is.null(splicing_strength)) splicing_strength <- 1.0
  if (is.null(variant_ratios)) variant_ratios <- c(1.0, 0.8)
  
  # Dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  n_genes_with_variants <- length(genes_with_variants)
  
  # Inizializza matrice di espressione modificata
  modified_expr <- expr_matrix
  
  # Prepara una lista per le matrici di espressione delle varianti
  variant_matrices <- list()
  
  # Per ogni gene con splicing
  for (i in 1:n_genes_with_variants) {
    gene_id <- genes_with_variants[i]
    n_variants <- variant_counts[gene_id]
    
    # Crea matrice di espressione per le varianti di questo gene
    gene_variants <- matrix(0, nrow = n_cells, ncol = n_variants)
    
    # Estrai espressione di base del gene
    base_expr <- expr_matrix[, gene_id]
    
    # Estrai propensioni per questo gene
    gene_propensity <- propensity[, i, ]
    
    # Per ogni cellula, distribuisci l'espressione tra le varianti
    for (c in 1:n_cells) {
      if (base_expr[c] > 0) {
        # Applica le proporzioni di splicing
        for (v in 1:n_variants) {
          # Effetto variante-specifico: applica il rapporto di espressione
          effect <- variant_ratios[min(v, length(variant_ratios))]
          gene_variants[c, v] <- base_expr[c] * gene_propensity[c, v] * effect * splicing_strength
        }
      }
    }
    
    # Memorizza la matrice delle varianti
    variant_matrices[[i]] <- gene_variants
    
    # Aggiorna l'espressione totale del gene (somma delle varianti)
    modified_expr[, gene_id] <- rowSums(gene_variants)
  }
  
  return(list(
    expr_matrix = modified_expr,
    variant_matrices = variant_matrices
  ))
}

#' Genera modello completo di splicing alternativo (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expr_matrix Matrice di espressione originale
#' @param splicing_params Parametri per lo splicing alternativo
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrice modificata e informazioni sulle varianti
#' @export
generate_alternative_splicing <- function(
  cell_df,
  expr_matrix,
  splicing_params = list(
    use_alternative_splicing = TRUE,
    splicing_fraction = 0.4,
    n_splicing_variants = 2,
    splicing_spatial_pattern = "gradient",
    splicing_strength = 0.8
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Verifica se lo splicing alternativo è abilitato
  use_splicing <- splicing_params$use_alternative_splicing
  if (is.null(use_splicing)) use_splicing <- TRUE
  
  # Se lo splicing è disabilitato, restituisci i valori originali
  if (!use_splicing) {
    return(list(
      expr_matrix = expr_matrix,
      genes_with_variants = integer(0),
      variant_matrices = list()
    ))
  }
  
  # Estrai dimensioni
  n_genes <- ncol(expr_matrix)
  
  # 1. Genera informazioni sulle varianti di splicing
  variants <- generate_splicing_variants(
    n_genes = n_genes,
    splicing_params = splicing_params,
    random_seed = random_seed
  )
  
  # Estrai i geni con varianti e conteggi
  genes_with_variants <- variants$genes_with_variants
  variant_counts <- variants$variant_counts
  
  # Se nessun gene ha varianti, restituisci i valori originali
  if (length(genes_with_variants) == 0) {
    return(list(
      expr_matrix = expr_matrix,
      genes_with_variants = integer(0),
      variant_matrices = list()
    ))
  }
  
  # 2. Calcola propensioni di splicing
  propensity <- calculate_splicing_propensity(
    cell_df = cell_df,
    genes_with_variants = genes_with_variants,
    n_variants = splicing_params$n_splicing_variants,
    splicing_params = splicing_params,
    random_seed = random_seed
  )
  
  # 3. Genera espressione delle varianti
  result <- generate_splicing_expression(
    expr_matrix = expr_matrix,
    genes_with_variants = genes_with_variants,
    variant_counts = variant_counts,
    propensity = propensity,
    splicing_params = splicing_params
  )
  
  # Restituisci i risultati
  return(list(
    expr_matrix = result$expr_matrix,
    genes_with_variants = genes_with_variants,
    variant_matrices = result$variant_matrices
  ))
}