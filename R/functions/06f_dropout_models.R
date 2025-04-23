#' Calcola probabilità di dropout per ogni cella
#'
#' Genera le probabilità di dropout in base alla distanza dai confini
#' o alla distanza media tra celle dello stesso tipo.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param mean_dist Vettore con le distanze medie per ogni cella
#' @param dropout_params Parametri di dropout
#' @param spatial_params Parametri spaziali
#' @importFrom scales rescale
#' @return Vettore con le probabilità base di dropout per ogni cella
#' @export
calculate_dropout_probabilities <- function(
  cell_df,
  mean_dist,
  dropout_params = list(
    dropout_range = c(0.2, 0.5)
  ),
  spatial_params = list(
    gradient_regions = FALSE
  )
) {
  # Calcola le probabilità di dropout base
  if (!is.null(spatial_params$gradient_regions) && 
      spatial_params$gradient_regions && 
      "boundary_dist" %in% colnames(cell_df) && 
      !all(is.na(cell_df$boundary_dist))) {
    # Più dropout vicino al confine
    base_dropout <- dropout_params$dropout_range[1] +
      (1 - cell_df$boundary_dist) * (dropout_params$dropout_range[2] - dropout_params$dropout_range[1])
  } else {
    # Metodo originale basato sulla distanza media
    base_dropout <- scales::rescale(mean_dist, to = dropout_params$dropout_range)
  }
  
  return(base_dropout)
}

#' Genera fattori di dropout specifici per gene
#'
#' Crea una lista di fattori che influenzano il dropout diversamente per ciascun gene.
#' Simula l'effetto di sequenza, lunghezza genica, contenuto GC, e altre proprietà
#' chimico-fisiche che rendono alcuni geni più suscettibili al dropout.
#'
#' @param n_genes Numero di geni
#' @param n_cells Numero di celle
#' @param dropout_params Parametri di dropout
#' @param random_seed Seed per riproducibilità
#' @return Lista contenente i fattori di dropout gene-specifici
#' @importFrom stats rnorm rbeta
#' @export
generate_gene_specific_dropout_factors <- function(
  n_genes,
  n_cells,
  dropout_params = list(
    gene_dropout_variability = 0.3,
    use_gene_specific_dropout = TRUE,
    gc_content_effect = 0.5,
    length_effect = 0.3,
    sequence_effect = 0.4
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Verifica se utilizzare il dropout gene-specifico
  use_gene_specific <- ifelse(is.null(dropout_params$use_gene_specific_dropout), 
                             TRUE, dropout_params$use_gene_specific_dropout)
  
  if (!use_gene_specific) {
    # Restituisci valori uniformi se il dropout gene-specifico non è richiesto
    return(list(
      gene_dropout_factors = rep(1, n_genes),
      dropout_baseline_shift = rep(0, n_genes)
    ))
  }
  
  # Determina la variabilità del dropout gene-specifico
  gene_dropout_variability <- ifelse(is.null(dropout_params$gene_dropout_variability), 
                                    0.3, dropout_params$gene_dropout_variability)
  
  # Intensità di ciascun effetto
  gc_effect <- ifelse(is.null(dropout_params$gc_content_effect), 0.5, dropout_params$gc_content_effect)
  length_effect <- ifelse(is.null(dropout_params$length_effect), 0.3, dropout_params$length_effect)  
  seq_effect <- ifelse(is.null(dropout_params$sequence_effect), 0.4, dropout_params$sequence_effect)
  
  # Genera simulazione del contenuto GC (distribuzione beta)
  gc_content <- rbeta(n_genes, 5, 5)
  
  # Genera simulazione della lunghezza del gene (log-normale)
  gene_length <- exp(rnorm(n_genes, log(2000), 0.7))
  
  # Genera effetto di sequenza (es. strutture secondarie, ripetizioni)
  sequence_factor <- rnorm(n_genes, 0, 1)
  
  # Combina gli effetti per creare il fattore di dropout gene-specifico
  # Alto valore = più suscettibile al dropout
  gene_dropout_factors <- 1 + 
    gc_effect * (gc_content - 0.5) +
    length_effect * (scale(log(gene_length))) +
    seq_effect * sequence_factor
  
  # Normalizza i fattori
  gene_dropout_factors <- scale(gene_dropout_factors)
  
  # Scala i fattori in base alla variabilità desiderata
  gene_dropout_factors <- 1 + gene_dropout_variability * gene_dropout_factors
  
  # Tronca i valori a un range ragionevole [0.5, 2.0]
  gene_dropout_factors <- pmin(pmax(gene_dropout_factors, 0.5), 2.0)
  
  # Crea anche uno shift di base per il dropout (alcuni geni sono intrinsecamente
  # più suscettibili al dropout indipendentemente dal livello di espressione)
  dropout_baseline_shift <- gene_dropout_variability * rnorm(n_genes, 0, 0.1)
  
  return(list(
    gene_dropout_factors = gene_dropout_factors,
    dropout_baseline_shift = dropout_baseline_shift,
    gc_content = gc_content,
    gene_length = gene_length
  ))
}

#' Genera contaminazione da RNA ambientale
#'
#' Simula la contaminazione da RNA ambientale (ambient RNA) che può diffondersi 
#' tra gli spot in tecnologie di trascrittomica spaziale.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param dist_mat Matrice di distanza tra celle (opzionale)
#' @param mean_expression_list Lista con le medie di espressione per tipo cellulare
#' @param ambient_params Parametri per la contaminazione ambientale
#' @param random_seed Seed per riproducibilità
#' @return Matrice di contaminazione ambientale
#' @importFrom stats rnorm runif
#' @export
generate_ambient_rna <- function(
  cell_df,
  dist_mat = NULL,
  mean_expression_list,
  ambient_params = list(
    ambient_contamination_rate = 0.05,
    ambient_diffusion_distance = 30,
    tissue_leakage_factor = 0.7,
    background_noise = 0.1
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  n_cells <- nrow(cell_df)
  n_genes <- length(mean_expression_list[[1]])
  k_cell_types <- length(mean_expression_list)
  
  # Determina il tasso di contaminazione ambientale
  ambient_rate <- ifelse(is.null(ambient_params$ambient_contamination_rate), 
                        0.05, ambient_params$ambient_contamination_rate)
  
  # Crea una "zuppa" globale di RNA ambientale come media pesata di tutti i tipi cellulari
  ambient_pool <- rep(0, n_genes)
  
  # Pesi differenti per tipi cellulari (alcuni possono contribuire più di altri)
  cell_type_weights <- runif(k_cell_types, 0.5, 1.5)
  cell_type_weights <- cell_type_weights / sum(cell_type_weights)
  
  # Calcola la media pesata per creare la zuppa di RNA ambientale
  for (k in 1:k_cell_types) {
    ambient_pool <- ambient_pool + cell_type_weights[k] * mean_expression_list[[k]]
  }
  
  # Normalizza
  ambient_pool <- ambient_pool / sum(ambient_pool) * n_genes
  
  # Applica un fattore di rumore casuale per ogni gene (alcuni geni si diffondono più facilmente)
  gene_diffusion_factors <- runif(n_genes, 0.5, 1.5)
  
  # Crea la matrice base di contaminazione ambientale
  ambient_contamination <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Contributo di background uniforme
  background_noise <- ifelse(is.null(ambient_params$background_noise), 
                           0.1, ambient_params$background_noise)
  
  # Aggiungi contributo base di ambient RNA (uguale per tutte le celle)
  for (g in 1:n_genes) {
    ambient_contamination[, g] <- ambient_pool[g] * gene_diffusion_factors[g]
  }
  
  # Aggiungi variazione spaziale se richiesta e se la matrice di distanza è disponibile
  if (!is.null(dist_mat) && !is.null(ambient_params$ambient_diffusion_distance) && 
      ambient_params$ambient_diffusion_distance > 0) {
    
    diffusion_distance <- ambient_params$ambient_diffusion_distance
    tissue_leakage <- ifelse(is.null(ambient_params$tissue_leakage_factor), 
                           0.7, ambient_params$tissue_leakage_factor)
    
    # Calcola l'influenza di ogni cellula sulla contaminazione delle altre
    # basandosi sulla distanza spaziale (computazionalmente costoso, usa solo per n_cells piccolo)
    if (n_cells <= 1000) {
      # Per ogni cella, aggiungi contaminazione dalle cellule vicine
      for (i in 1:n_cells) {
        # Calcola il peso di ogni cella j sulla cella i
        weights <- exp(-dist_mat[i, ] / diffusion_distance)
        weights[i] <- 0  # Escludi auto-contaminazione
        weights <- weights / (sum(weights) + 1e-10)
        
        # Per ogni gene, aggiungi contaminazione proporzionale alla distanza
        for (g in 1:n_genes) {
          local_ambient <- sum(weights * ambient_contamination[, g])
          ambient_contamination[i, g] <- ambient_contamination[i, g] * (1 - tissue_leakage) + 
                                      local_ambient * tissue_leakage
        }
      }
    } else {
      # Versione semplificata per dataset grandi
      # Aggiungi variabilità spaziale basata su cluster
      cluster_ids <- as.numeric(as.factor(cell_df$intensity_cluster))
      n_clusters <- max(cluster_ids)
      
      # Crea un profilo ambientale specifico per cluster
      for (cl in 1:n_clusters) {
        cells_in_cluster <- which(cluster_ids == cl)
        n_cells_in_cluster <- length(cells_in_cluster)
        
        if (n_cells_in_cluster > 0) {
          # Crea una variazione cluster-specifica
          cluster_factor <- runif(n_genes, 0.8, 1.2)
          
          # Applica alle celle di questo cluster
          for (g in 1:n_genes) {
            ambient_contamination[cells_in_cluster, g] <- 
              ambient_contamination[cells_in_cluster, g] * cluster_factor[g]
          }
        }
      }
    }
  }
  
  # Aggiungi una componente di rumore di background
  if (background_noise > 0) {
    background <- matrix(
      rnorm(n_cells * n_genes, mean = 0, sd = background_noise),
      nrow = n_cells, ncol = n_genes
    )
    # Assicurati che il rumore non renda negative le contaminazioni
    ambient_contamination <- pmax(ambient_contamination + background, 0)
  }
  
  # Scala la contaminazione in base al tasso desiderato
  ambient_contamination <- ambient_contamination * ambient_rate
  
  return(ambient_contamination)
}

#' Applica contaminazione da RNA ambientale alla matrice di espressione
#'
#' Aggiunge la contaminazione da RNA ambientale alla matrice di espressione.
#'
#' @param expression_matrix Matrice di espressione originale
#' @param ambient_contamination Matrice di contaminazione ambientale
#' @param ambient_params Parametri per la contaminazione ambientale
#' @return Matrice di espressione con contaminazione ambientale applicata
#' @export
apply_ambient_rna_contamination <- function(
  expression_matrix,
  ambient_contamination,
  ambient_params = list(
    ambient_scale_factor = 1.0
  )
) {
  # Verifica dimensioni compatibili
  if (!identical(dim(expression_matrix), dim(ambient_contamination))) {
    stop("Le dimensioni della matrice di espressione e della contaminazione ambientale devono corrispondere")
  }
  
  # Estrai fattore di scala
  scale_factor <- ifelse(is.null(ambient_params$ambient_scale_factor), 
                        1.0, ambient_params$ambient_scale_factor)
  
  # Applica la contaminazione: matrice originale + contaminazione scalata
  contaminated_matrix <- expression_matrix + ambient_contamination * scale_factor
  
  # Arrotonda a numeri interi (conteggi)
  contaminated_matrix <- round(contaminated_matrix)
  
  return(contaminated_matrix)
}

#' Applica dropout all'espressione genica
#'
#' Applica dropout all'espressione genica in modo espressione-dipendente e/o spaziale,
#' considerando anche fattori gene-specifici.
#'
#' @param expression Matrice o vettore di espressione
#' @param base_dropout Probabilità base di dropout
#' @param gene_specific_factors Lista con fattori di dropout gene-specifici
#' @param dropout_params Parametri di dropout
#' @param random_seed Seed per riproducibilità
#' @return Matrice o vettore di espressione con dropout applicato
#' @export
apply_dropout <- function(
  expression,
  base_dropout,
  gene_specific_factors = NULL,
  dropout_params = list(
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5,
    use_gene_specific_dropout = TRUE,
    gene_effect_weight = 0.3
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Determina se l'input è una matrice o un vettore
  is_matrix <- is.matrix(expression)
  
  # Se è un vettore, trattalo come una matrice con una colonna
  if (!is_matrix) {
    expression <- matrix(expression, ncol = 1)
  }
  
  # Definisci una funzione vettorizzata per normalizzare tra 0 e 1
  scale01_vec <- function(x) {
    if (all(x == x[1])) return(rep(0.5, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }
  
  # Numero di celle e geni
  N <- nrow(expression)
  G <- ncol(expression)
  
  # Verifica se utilizzare il dropout gene-specifico
  use_gene_specific <- ifelse(is.null(dropout_params$use_gene_specific_dropout), 
                             TRUE, dropout_params$use_gene_specific_dropout)
  
  # Peso dell'effetto gene-specifico
  gene_effect_weight <- ifelse(is.null(dropout_params$gene_effect_weight), 
                              0.3, dropout_params$gene_effect_weight)
  
  # Verifica disponibilità dei fattori gene-specifici
  if (use_gene_specific && is.null(gene_specific_factors)) {
    # Se richiesti ma non forniti, crea valori di default
    gene_specific_factors <- list(
      gene_dropout_factors = rep(1, G),
      dropout_baseline_shift = rep(0, G)
    )
    use_gene_specific <- FALSE  # Disabilita se non sono disponibili fattori reali
  }
  
  # Applica dropout a ogni colonna (gene)
  for (g in 1:G) {
    # Applica fattori gene-specifici se richiesto
    gene_factor <- 1
    gene_shift <- 0
    
    if (use_gene_specific) {
      gene_factor <- gene_specific_factors$gene_dropout_factors[g]
      gene_shift <- gene_specific_factors$dropout_baseline_shift[g]
    }
    
    if (dropout_params$expression_dependent_dropout) {
      # Normalizza l'espressione del gene corrente
      norm_expr <- scale01_vec(expression[, g])
      
      # Calcola la probabilità di dropout con una funzione logistica
      # Aggiusta il midpoint e la pendenza in base ai fattori gene-specifici
      adjusted_midpoint <- dropout_params$dropout_curve_midpoint
      adjusted_steepness <- dropout_params$dropout_curve_steepness / gene_factor
      
      if (use_gene_specific) {
        adjusted_midpoint <- adjusted_midpoint + gene_shift
      }
      
      dropout_prob_expr <- 1 / (1 + exp((norm_expr - adjusted_midpoint) * adjusted_steepness))
      
      # Combina con il dropout spaziale base
      # I pesi sono ora: espressione-dipendente (0.6) + spaziale (0.3) + gene-specifico (0.1)
      if (use_gene_specific) {
        dropout_prob <- (1 - gene_effect_weight) * 
          (0.7 * dropout_prob_expr + 0.3 * base_dropout) + 
          gene_effect_weight * (gene_factor * base_dropout + gene_shift)
      } else {
        dropout_prob <- 0.7 * dropout_prob_expr + 0.3 * base_dropout
      }
      
      # Tronca i valori al range [0,1]
      dropout_prob <- pmin(pmax(dropout_prob, 0), 1)
    } else {
      # Usa solo il modello di dropout spaziale con influenze gene-specifiche
      if (use_gene_specific) {
        dropout_prob <- base_dropout * gene_factor + gene_shift
        dropout_prob <- pmin(pmax(dropout_prob, 0), 1)
      } else {
        dropout_prob <- base_dropout
      }
    }
    
    # Applica dropout
    zero_idx <- runif(N) < dropout_prob
    expression[zero_idx, g] <- 0
  }
  
  # Restituisci nella forma originale
  if (!is_matrix) {
    expression <- as.vector(expression)
  }
  
  return(expression)
}