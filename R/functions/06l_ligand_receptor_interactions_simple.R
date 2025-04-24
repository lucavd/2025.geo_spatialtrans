#' Modellazione semplificata di interazioni ligando-recettore per testing
#'
#' Versioni semplificate delle funzioni di interazioni ligando-recettore adattate per testing
#' 
#' @import stats

#' Definisce un database di interazioni ligando-recettore (versione semplificata per testing)
#'
#' @param n_interactions Numero di interazioni L-R
#' @param gene_pool Vettore di ID geni tra cui scegliere (optional)
#' @param strength_range Range di intensità di interazione
#' @param random_seed Seed per riproducibilità
#' @return Dataframe con interazioni L-R
#' @export
define_ligand_receptor_db <- function(
  n_interactions = 10,
  gene_pool = NULL,
  strength_range = c(0.1, 1.0),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Crea pool di geni se non specificato
  if (is.null(gene_pool)) {
    gene_pool <- 1:100  # Default: 100 geni
  }
  
  n_genes <- length(gene_pool)
  
  # Crea database di interazioni L-R
  lr_db <- data.frame(
    ligand = sample(gene_pool, n_interactions),
    receptor = sample(gene_pool, n_interactions),
    strength = runif(n_interactions, strength_range[1], strength_range[2]),
    range = runif(n_interactions, 5, 20)
  )
  
  # Assicura che ligando e recettore siano diversi
  for (i in 1:n_interactions) {
    while (lr_db$ligand[i] == lr_db$receptor[i]) {
      lr_db$receptor[i] <- sample(gene_pool, 1)
    }
  }
  
  return(lr_db)
}

#' Calcola effetti di segnalazione ligando-recettore (versione semplificata per testing)
#'
#' @param expr_matrix Matrice di espressione
#' @param lr_db Database delle interazioni L-R
#' @param dist_mat Matrice delle distanze tra cellule
#' @param params Parametri per il calcolo degli effetti
#' @param random_seed Seed per riproducibilità
#' @return Matrice degli effetti di segnalazione
#' @export
compute_lr_signaling_effects <- function(
  expr_matrix,
  lr_db,
  dist_mat,
  params = list(
    lr_distance_decay = "exponential",
    lr_effect_type = "multiplicative"
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  decay_type <- params$lr_distance_decay
  effect_type <- params$lr_effect_type
  
  # Valori predefiniti se non specificati
  if (is.null(decay_type)) decay_type <- "exponential"
  if (is.null(effect_type)) effect_type <- "multiplicative"
  
  # Estrai dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  n_interactions <- nrow(lr_db)
  
  # Inizializza matrice di effetti
  effects_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Per ogni interazione L-R
  for (i in 1:n_interactions) {
    ligand_id <- lr_db$ligand[i]
    receptor_id <- lr_db$receptor[i]
    interaction_strength <- lr_db$strength[i]
    interaction_range <- lr_db$range[i]
    
    # Salta se i geni non sono nella matrice di espressione
    if (ligand_id > n_genes || receptor_id > n_genes) {
      next
    }
    
    # Estrai espressione di ligando e recettore
    ligand_expr <- expr_matrix[, ligand_id]
    receptor_expr <- expr_matrix[, receptor_id]
    
    # Per ogni cellula emittente
    for (sender in 1:n_cells) {
      # Salta se non esprime il ligando
      if (ligand_expr[sender] <= 0) {
        next
      }
      
      # Per ogni cellula ricevente
      for (receiver in 1:n_cells) {
        # Salta se è la stessa cellula o non esprime il recettore
        if (sender == receiver || receptor_expr[receiver] <= 0) {
          next
        }
        
        # Distanza tra cellule
        distance <- dist_mat[sender, receiver]
        
        # Calcola fattore di decadimento con la distanza
        decay_factor <- 0
        if (decay_type == "exponential") {
          decay_factor <- exp(-distance / interaction_range)
        } else if (decay_type == "linear") {
          decay_factor <- max(0, 1 - distance / interaction_range)
        } else {
          # Quadratico come fallback
          decay_factor <- max(0, (1 - distance / interaction_range)^2)
        }
        
        # Calcola intensità del segnale
        signal_strength <- ligand_expr[sender] * receptor_expr[receiver] * 
                          interaction_strength * decay_factor
        
        # Aggiungi effetto alla cellula ricevente
        if (signal_strength > 0) {
          # Target i geni a caso con l'effetto, per semplicità
          n_target_genes <- max(2, round(n_genes * 0.05))
          target_genes <- sample(1:n_genes, n_target_genes)
          
          for (target in target_genes) {
            effects_matrix[receiver, target] <- effects_matrix[receiver, target] + 
                                               signal_strength / n_target_genes
          }
        }
      }
    }
  }
  
  # Normalizza gli effetti
  if (max(effects_matrix) > 0) {
    effects_matrix <- effects_matrix / max(effects_matrix)
  }
  
  return(effects_matrix)
}

#' Applica effetti di segnalazione all'espressione (versione semplificata per testing)
#'
#' @param expr_matrix Matrice di espressione originale
#' @param effect_matrix Matrice degli effetti di segnalazione
#' @param params Parametri per l'applicazione degli effetti
#' @param random_seed Seed per riproducibilità
#' @return Matrice di espressione modificata
#' @export
apply_lr_signaling <- function(
  expr_matrix,
  effect_matrix,
  params = list(
    lr_effect_type = "multiplicative",
    lr_effect_strength = 1.0
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  effect_type <- params$lr_effect_type
  effect_strength <- params$lr_effect_strength
  
  # Valori predefiniti se non specificati
  if (is.null(effect_type)) effect_type <- "multiplicative"
  if (is.null(effect_strength)) effect_strength <- 1.0
  
  # Scala gli effetti per l'intensità desiderata
  scaled_effects <- effect_matrix * effect_strength
  
  # Applica gli effetti in base al tipo
  if (effect_type == "multiplicative") {
    # Effetto moltiplicativo: 1 + effetto (quindi 1 = nessun effetto)
    factors <- 1 + scaled_effects
    modified_expr <- expr_matrix * factors
    
  } else if (effect_type == "additive") {
    # Effetto additivo
    modified_expr <- expr_matrix + scaled_effects * mean(expr_matrix)
    
  } else {
    # Approccio misto come fallback
    factors <- 1 + scaled_effects * 0.5
    add_component <- scaled_effects * 0.5 * mean(expr_matrix)
    modified_expr <- expr_matrix * factors + add_component
  }
  
  # Assicura che l'espressione rimanga non-negativa
  modified_expr[modified_expr < 0] <- 0
  
  return(modified_expr)
}

#' Genera modello completo di interazioni ligando-recettore (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expr_matrix Matrice di espressione originale
#' @param params Parametri per le interazioni L-R
#' @param random_seed Seed per riproducibilità
#' @return Lista con risultati delle interazioni L-R
#' @export
generate_lr_interactions <- function(
  cell_df,
  expr_matrix,
  dist_mat = NULL,
  params = list(
    use_lr_interactions = TRUE,
    n_lr_interactions = 5,
    lr_distance_decay = "exponential",
    lr_effect_type = "multiplicative",
    lr_effect_strength = 0.8
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Verifica se le interazioni L-R sono abilitate
  use_lr <- params$use_lr_interactions
  if (is.null(use_lr)) use_lr <- TRUE
  
  # Se disabilitate, restituisci valori originali
  if (!use_lr) {
    return(list(
      expr_matrix = expr_matrix,
      lr_database = data.frame(),
      lr_effects = matrix(0, nrow(expr_matrix), ncol(expr_matrix))
    ))
  }
  
  # Numero di interazioni
  n_interactions <- params$n_lr_interactions
  if (is.null(n_interactions)) n_interactions <- 5
  
  # Estrai dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  
  # Calcola matrice delle distanze se non fornita
  if (is.null(dist_mat)) {
    dist_mat <- as.matrix(dist(cell_df[, c("x", "y")]))
  }
  
  # 1. Definisci database di interazioni L-R
  lr_db <- define_ligand_receptor_db(
    n_interactions = n_interactions,
    gene_pool = 1:n_genes,
    random_seed = random_seed
  )
  
  # 2. Calcola effetti di segnalazione
  lr_effects <- compute_lr_signaling_effects(
    expr_matrix = expr_matrix,
    lr_db = lr_db,
    dist_mat = dist_mat,
    params = list(
      lr_distance_decay = params$lr_distance_decay,
      lr_effect_type = params$lr_effect_type
    ),
    random_seed = random_seed
  )
  
  # 3. Applica effetti all'espressione
  modified_expr <- apply_lr_signaling(
    expr_matrix = expr_matrix,
    effect_matrix = lr_effects,
    params = list(
      lr_effect_type = params$lr_effect_type,
      lr_effect_strength = params$lr_effect_strength
    ),
    random_seed = random_seed
  )
  
  # Restituisci risultati
  return(list(
    expr_matrix = modified_expr,
    lr_database = lr_db,
    lr_effects = lr_effects
  ))
}