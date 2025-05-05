#' Modellazione semplificata di interazioni ligando-recettore per testing
#'
#' Versioni semplificate delle funzioni di interazioni ligando-recettore adattate per testing
#' 
#' @import stats

#' Definisce un database di interazioni ligando-recettore (versione semplificata per testing)
#'
#' @param n_genes Numero totale di geni nella simulazione (o numero di interazioni se non specificato)
#' @param gene_pool Vettore di ID geni tra cui scegliere (optional)
#' @param strength_range Range di intensità di interazione
#' @param random_seed Seed per riproducibilità
#' @return Dataframe con interazioni L-R
#' @export
define_ligand_receptor_db <- function(
  n_genes = 100,
  n_interactions = NULL,
  gene_pool = NULL,
  strength_range = c(0.1, 1.0),
  interaction_params = NULL,
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Gestisci parametri aggiuntivi
  if (!is.null(interaction_params)) {
    if (!is.null(interaction_params$n_interactions)) {
      n_interactions <- interaction_params$n_interactions
    }
  }
  
  # Determina il numero di interazioni
  if (is.null(n_interactions)) {
    n_interactions = 20  # Valore predefinito è 20
  }
  
  # Crea pool di geni se non specificato
  if (is.null(gene_pool)) {
    gene_pool <- 1:n_genes
  }
  
  # Crea lista di interazioni L-R
  interactions <- list()
  
  for (i in 1:n_interactions) {
    # Seleziona geni ligando e recettore
    ligand_id <- sample(gene_pool, 1)
    receptor_id <- sample(gene_pool, 1)
    
    # Assicura che siano diversi
    while (ligand_id == receptor_id) {
      receptor_id <- sample(gene_pool, 1)
    }
    
    # Seleziona geni target (da 2 a 10)
    n_targets <- sample(2:10, 1)
    target_ids <- sample(gene_pool, n_targets)
    
    # Genera effetti per ogni target (positivi o negativi)
    target_effects <- runif(n_targets, -1, 1)
    
    # Crea voce per questa interazione
    interactions[[i]] <- list(
      ligand_id = ligand_id,
      receptor_id = receptor_id,
      target_ids = target_ids,
      target_effects = target_effects,
      decay_factor = runif(1, 5, 30),
      activation_threshold = runif(1, 0.1, 0.3),
      interaction_type = sample(c("activating", "inhibitory"), 1)
    )
  }
  
  # Crea ruoli dei geni
  gene_roles <- rep("other", n_genes)
  
  # Assegna ruoli basati sulle interazioni
  for (interaction in interactions) {
    gene_roles[interaction$ligand_id] <- "ligand"
    gene_roles[interaction$receptor_id] <- "receptor"
  }
  
  # Costruisci e restituisci il database completo
  return(list(
    interactions = interactions,
    gene_roles = gene_roles
  ))
}

#' Calcola effetti di segnalazione ligando-recettore (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expr_matrix Matrice di espressione
#' @param lr_db Database delle interazioni L-R
#' @param dist_mat Matrice delle distanze tra cellule
#' @param params Parametri per il calcolo degli effetti
#' @param interaction_params Alias alternativo per params (per compatibilità con test)
#' @param random_seed Seed per riproducibilità
#' @return Matrice degli effetti di segnalazione
#' @export
compute_lr_signaling_effects <- function(
  cell_df,
  expr_matrix,
  lr_db,
  dist_mat,
  params = list(
    signal_propagation_mode = "exponential",
    max_signaling_distance = 40,
    signal_amplification = 1.0
  ),
  interaction_params = NULL,
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Gestisci parametri da entrambe le fonti possibili
  if (!is.null(interaction_params)) {
    params <- interaction_params
  }
  
  # Estrai parametri
  propagation_mode <- "exponential"
  max_distance <- 40
  amplification <- 1.0
  
  # Estrai dai parametri se disponibili
  if (!is.null(params) && is.list(params)) {
    if (!is.null(params$signal_propagation_mode)) {
      propagation_mode <- params$signal_propagation_mode
    }
    if (!is.null(params$max_signaling_distance)) {
      max_distance <- params$max_signaling_distance
    }
    if (!is.null(params$signal_amplification)) {
      amplification <- params$signal_amplification
    }
  }
  
  # Estrai dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  
  # Inizializza matrice di effetti
  effects_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Accesso alle interazioni tramite la struttura list
  interactions <- lr_db$interactions
  
  # Per ogni interazione L-R
  for (interaction in interactions) {
    ligand_id <- interaction$ligand_id
    receptor_id <- interaction$receptor_id
    target_ids <- interaction$target_ids
    target_effects <- interaction$target_effects
    decay_factor <- interaction$decay_factor
    
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
        
        # Salta se la distanza supera il massimo
        if (distance > max_distance) {
          next
        }
        
        # Calcola fattore di decadimento con la distanza
        dist_factor <- 0
        if (propagation_mode == "exponential") {
          dist_factor <- exp(-distance / decay_factor)
        } else if (propagation_mode == "linear") {
          dist_factor <- max(0, 1 - distance / decay_factor)
        } else if (propagation_mode == "threshold") {
          dist_factor <- ifelse(distance <= decay_factor, 1, 0)
        } else {
          # Quadratico come fallback
          dist_factor <- max(0, (1 - distance / decay_factor)^2)
        }
        
        # Calcola intensità del segnale
        signal_strength <- ligand_expr[sender] * receptor_expr[receiver] * dist_factor * amplification
        
        # Aggiungi effetto alla cellula ricevente per i geni target
        if (signal_strength > 0) {
          for (t in 1:length(target_ids)) {
            target_id <- target_ids[t]
            effect <- target_effects[t]
            
            if (target_id <= n_genes) {
              effects_matrix[receiver, target_id] <- effects_matrix[receiver, target_id] + 
                                                   signal_strength * effect
            }
          }
        }
      }
    }
  }
  
  # Normalizza gli effetti se necessario
  effects_range <- range(effects_matrix)
  if (effects_range[2] > effects_range[1]) {
    effects_matrix <- (effects_matrix - effects_range[1]) / (effects_range[2] - effects_range[1]) * 2 - 1
  }
  
  return(effects_matrix)
}

#' Applica effetti di segnalazione all'espressione (versione semplificata per testing)
#'
#' @param expr_matrix Matrice di espressione originale
#' @param effect_matrix Matrice degli effetti di segnalazione
#' @param params Parametri per l'applicazione degli effetti
#' @param integration_params Alias alternativo per params (per compatibilità con test)
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
  integration_params = NULL,
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Gestisci parametri da entrambe le fonti possible
  if (!is.null(integration_params)) {
    params <- integration_params
  }
  
  # Estrai parametri
  effect_type <- "multiplicative"  # Default
  effect_strength <- 1.0  # Default
  
  # Estrai dai parametri se disponibili
  if (!is.null(params) && is.list(params)) {
    if (!is.null(params$lr_effect_type)) {
      effect_type <- params$lr_effect_type
    }
    if (!is.null(params$lr_effect_strength)) {
      effect_strength <- params$lr_effect_strength
    }
  }
  
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
#' @param dist_mat Matrice delle distanze tra cellule (opzionale)
#' @param params Parametri per le interazioni L-R
#' @param lr_params Alias alternativo per params (per compatibilità con test)
#' @param random_seed Seed per riproducibilità
#' @return Lista con risultati delle interazioni L-R
#' @export
generate_lr_interactions <- function(
  cell_df,
  expr_matrix,
  dist_mat = NULL,
  params = list(
    use_lr_interactions = TRUE,
    n_interactions = 20,
    signal_propagation_mode = "exponential",
    max_signaling_distance = 40,
    adjust_method = "multiplicative"
  ),
  lr_params = NULL,
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Gestisci parametri da entrambe le fonti possibili
  if (!is.null(lr_params)) {
    params <- lr_params
  }
  
  # Verifica se le interazioni L-R sono abilitate
  use_lr <- TRUE  # Default
  if (!is.null(params) && is.list(params) && !is.null(params$use_lr_interactions)) {
    use_lr <- params$use_lr_interactions
  }
  
  # Se disabilitate, restituisci valori originali
  if (!use_lr) {
    return(list(
      expression = expr_matrix,
      original_expression = expr_matrix,
      signaling_effects = matrix(0, nrow(expr_matrix), ncol(expr_matrix)),
      interaction_db = list(interactions = list())
    ))
  }
  
  # Numero di interazioni
  n_interactions <- 20  # Default
  if (!is.null(params) && is.list(params)) {
    if (!is.null(params$n_interactions)) {
      n_interactions <- params$n_interactions
    } else if (!is.null(params$n_lr_interactions)) {
      n_interactions <- params$n_lr_interactions
    }
  }
  
  # Estrai dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  
  # Calcola matrice delle distanze se non fornita
  if (is.null(dist_mat)) {
    dist_mat <- as.matrix(dist(cell_df[, c("x", "y")]))
  }
  
  # 1. Definisci database di interazioni L-R
  lr_db <- define_ligand_receptor_db(
    n_genes = n_genes,
    n_interactions = n_interactions,
    interaction_params = params,
    random_seed = random_seed
  )
  
  # 2. Calcola effetti di segnalazione
  lr_effects <- compute_lr_signaling_effects(
    cell_df = cell_df,
    expr_matrix = expr_matrix,
    lr_db = lr_db,
    dist_mat = dist_mat,
    params = params,
    random_seed = random_seed
  )
  
  # 3. Applica effetti all'espressione
  modified_expr <- apply_lr_signaling(
    expr_matrix = expr_matrix,
    effect_matrix = lr_effects,
    params = params,
    random_seed = random_seed
  )
  
  # Restituisci risultati con struttura e nomi compatibili con i test
  return(list(
    expression = modified_expr,
    original_expression = expr_matrix,
    signaling_effects = lr_effects,
    interaction_db = list(interactions = lr_db$interactions),
    
    # Alias aggiuntivi per altri usi possibili
    expr_matrix = modified_expr,
    lr_database = lr_db,
    lr_effects = lr_effects
  ))
}