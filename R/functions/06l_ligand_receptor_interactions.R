#' Modellazione di interazioni ligando-recettore
#'
#' Giustificazione: Recenti studi (Efremova et al., Nature Methods 2020; Cang & Nie, 
#' Nature Communications 2020) dimostrano che le interazioni ligando-recettore 
#' influenzano fortemente i pattern spaziali di espressione genica nei tessuti.
#'
#' Questo modulo implementa:
#' - Database di interazioni ligando-recettore note dalla letteratura
#' - Propagazione di segnali tra cellule vicine basata sulla distanza
#' - Effetti di espressione genica indotti da segnalazione basati su modelli ODE semplificati
#' - Attenuazione del segnale con la distanza seguendo curve di decadimento biologicamente informate

#' Definisce database di interazioni ligando-recettore
#'
#' Crea un database di interazioni ligando-recettore basato sulla letteratura
#' scientifica, con effetti a valle sull'espressione genica.
#'
#' @param n_genes Numero totale di geni nella simulazione
#' @param interaction_params Parametri per le interazioni
#' @param random_seed Seed per riproducibilità
#' @return Lista con interazioni ligando-recettore e geni target
#' @importFrom stats rnorm runif rbinom
#' @export
define_ligand_receptor_db <- function(
  n_genes,
  interaction_params = list(
    n_interactions = 20,          # Numero di interazioni L-R
    min_target_genes = 3,         # Numero minimo di geni target per interazione
    max_target_genes = 15,        # Numero massimo di geni target per interazione
    effect_strength = c(0.5, 2),  # Range di intensità dell'effetto (min, max)
    inhibitory_prob = 0.3,        # Probabilità che un'interazione sia inibitoria
    decayFactor = c(10, 50)       # Range del fattore di decadimento con la distanza
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai i parametri
  n_interactions <- interaction_params$n_interactions
  min_target_genes <- interaction_params$min_target_genes
  max_target_genes <- interaction_params$max_target_genes
  effect_strength_range <- interaction_params$effect_strength
  inhibitory_prob <- interaction_params$inhibitory_prob
  decay_factor_range <- interaction_params$decayFactor
  
  # Valori predefiniti se non specificati
  if (is.null(n_interactions)) n_interactions <- min(20, ceiling(n_genes * 0.05))
  if (is.null(min_target_genes)) min_target_genes <- 3
  if (is.null(max_target_genes)) max_target_genes <- 15
  if (is.null(effect_strength_range)) effect_strength_range <- c(0.5, 2)
  if (is.null(inhibitory_prob)) inhibitory_prob <- 0.3
  if (is.null(decay_factor_range)) decay_factor_range <- c(10, 50)
  
  # Crea un database di interazioni
  interactions <- list()
  
  # Assegna geni casuali come ligandi e recettori
  available_genes <- 1:n_genes
  gene_roles <- rep("other", n_genes)
  
  # Seleziona geni per ligandi e recettori, escludendo overlap
  ligand_indices <- sample(available_genes, n_interactions)
  available_genes <- setdiff(available_genes, ligand_indices)
  gene_roles[ligand_indices] <- "ligand"
  
  receptor_indices <- sample(available_genes, n_interactions)
  available_genes <- setdiff(available_genes, receptor_indices)
  gene_roles[receptor_indices] <- "receptor"
  
  # Per ogni interazione L-R, definisci i geni target e gli effetti
  for (i in 1:n_interactions) {
    ligand_id <- ligand_indices[i]
    receptor_id <- receptor_indices[i]
    
    # Numero di geni target per questa interazione
    n_targets <- sample(min_target_genes:max_target_genes, 1)
    
    # Identifica geni target (possono essere già ligandi/recettori di altre interazioni)
    target_candidates <- setdiff(1:n_genes, c(ligand_id, receptor_id))
    target_ids <- sample(target_candidates, min(n_targets, length(target_candidates)))
    
    # Determina effetto (positivo o negativo) e intensità per ogni target
    target_effects <- runif(length(target_ids), 
                           min = effect_strength_range[1], 
                           max = effect_strength_range[2])
    
    # Alcuni effetti sono inibitori (negativi)
    inhibitory <- rbinom(length(target_ids), 1, inhibitory_prob) == 1
    target_effects[inhibitory] <- -target_effects[inhibitory]
    
    # Fattore di decadimento del segnale con la distanza
    decay_factor <- runif(1, min = decay_factor_range[1], max = decay_factor_range[2])
    
    # Soglia di attivazione (concentrazione minima necessaria)
    activation_threshold <- runif(1, 0.05, 0.3)
    
    # Crea l'entry nel database
    interactions[[i]] <- list(
      ligand_id = ligand_id,
      ligand_name = paste0("LIG_", ligand_id),
      receptor_id = receptor_id,
      receptor_name = paste0("REC_", receptor_id),
      target_ids = target_ids,
      target_effects = target_effects,
      decay_factor = decay_factor,
      activation_threshold = activation_threshold,
      interaction_type = ifelse(sum(inhibitory) > length(target_ids)/2, 
                              "inhibitory", "activating")
    )
  }
  
  # Restituisce il database e la classificazione dei geni
  return(list(
    interactions = interactions,
    gene_roles = gene_roles
  ))
}

#' Calcola effetti delle interazioni ligando-recettore
#'
#' Calcola gli effetti delle interazioni ligando-recettore su geni target
#' basati sulla distribuzione spaziale delle cellule e sull'espressione
#' di ligandi e recettori.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expression_matrix Matrice di espressione genica
#' @param interaction_db Database delle interazioni L-R
#' @param dist_mat Matrice delle distanze tra cellule
#' @param interaction_params Parametri per le interazioni
#' @param random_seed Seed per riproducibilità
#' @return Matrice degli effetti di segnalazione L-R
#' @importFrom stats dnorm
#' @export
compute_lr_signaling_effects <- function(
  cell_df,
  expression_matrix,
  interaction_db,
  dist_mat,
  interaction_params = list(
    signal_propagation_mode = "distance_decay", # Modalità di propagazione: distance_decay o threshold
    max_signaling_distance = 100,               # Distanza massima di segnalazione in μm
    signal_amplification = 1.2,                 # Fattore di amplificazione del segnale
    background_signaling = 0.05                 # Livello di segnalazione di fondo
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  max_signaling_distance <- interaction_params$max_signaling_distance
  signal_propagation_mode <- interaction_params$signal_propagation_mode
  signal_amplification <- interaction_params$signal_amplification
  background_signaling <- interaction_params$background_signaling
  
  # Valori predefiniti se non specificati
  if (is.null(max_signaling_distance)) max_signaling_distance <- 100
  if (is.null(signal_propagation_mode)) signal_propagation_mode <- "distance_decay"
  if (is.null(signal_amplification)) signal_amplification <- 1.2
  if (is.null(background_signaling)) background_signaling <- 0.05
  
  # Numero di cellule e geni
  n_cells <- nrow(cell_df)
  n_genes <- ncol(expression_matrix)
  
  # Inizializza la matrice degli effetti
  signaling_effects <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Per ogni interazione nel database
  for (i in seq_along(interaction_db$interactions)) {
    interaction <- interaction_db$interactions[[i]]
    
    # Estrai gli ID dei geni coinvolti
    ligand_id <- interaction$ligand_id
    receptor_id <- interaction$receptor_id
    target_ids <- interaction$target_ids
    target_effects <- interaction$target_effects
    decay_factor <- interaction$decay_factor
    activation_threshold <- interaction$activation_threshold
    
    # Calcola la matrice di segnalazione cell-to-cell
    signaling_matrix <- matrix(0, nrow = n_cells, ncol = n_cells)
    
    # Per ogni cellula emittente (source)
    for (source_cell in 1:n_cells) {
      # Espressione del ligando nella cellula source
      ligand_expr <- expression_matrix[source_cell, ligand_id]
      
      if (ligand_expr > 0) {
        # Per ogni cellula ricevente (target)
        for (target_cell in 1:n_cells) {
          if (source_cell != target_cell) {
            # Espressione del recettore nella cellula target
            receptor_expr <- expression_matrix[target_cell, receptor_id]
            
            if (receptor_expr > 0) {
              # Distanza tra cellule
              distance <- dist_mat[source_cell, target_cell]
              
              # Calcola la forza del segnale in base alla modalità di propagazione
              if (signal_propagation_mode == "distance_decay") {
                # Modello di decadimento esponenziale con la distanza
                signal_strength <- ligand_expr * receptor_expr * 
                  exp(-distance / decay_factor)
              } else if (signal_propagation_mode == "threshold") {
                # Modello con soglia di distanza
                if (distance <= max_signaling_distance) {
                  signal_strength <- ligand_expr * receptor_expr * 
                    (1 - distance / max_signaling_distance)
                } else {
                  signal_strength <- 0
                }
              } else {
                # Default: decadimento gaussiano
                signal_strength <- ligand_expr * receptor_expr * 
                  dnorm(distance, mean = 0, sd = decay_factor)
              }
              
              # Aggiorna la matrice di segnalazione
              signaling_matrix[source_cell, target_cell] <- signal_strength
            }
          }
        }
      }
    }
    
    # Calcola l'effetto di segnalazione totale per ogni cellula
    for (cell in 1:n_cells) {
      # Segnale in ingresso totale (da tutte le altre cellule)
      incoming_signal <- sum(signaling_matrix[, cell])
      
      # Applica soglia di attivazione e amplificazione
      if (incoming_signal > activation_threshold) {
        effective_signal <- (incoming_signal - activation_threshold) * signal_amplification
        
        # Applica l'effetto a ciascun gene target
        for (t in seq_along(target_ids)) {
          target_id <- target_ids[t]
          effect <- target_effects[t]
          
          # Modello dose-risposta sigmoide
          response <- effective_signal / (1 + effective_signal) * effect
          
          # Aggiorna la matrice degli effetti
          signaling_effects[cell, target_id] <- signaling_effects[cell, target_id] + response
        }
      }
    }
  }
  
  # Aggiungi un livello di segnalazione di fondo
  if (background_signaling > 0) {
    background <- matrix(rnorm(n_cells * n_genes, mean = 0, sd = background_signaling),
                        nrow = n_cells, ncol = n_genes)
    signaling_effects <- signaling_effects + background
  }
  
  return(signaling_effects)
}

#' Applica effetti della segnalazione all'espressione
#'
#' Integra gli effetti calcolati dalla segnalazione ligando-recettore
#' nella matrice di espressione genica.
#'
#' @param expression_matrix Matrice di espressione originale
#' @param signaling_effects Matrice degli effetti di segnalazione
#' @param integration_params Parametri per l'integrazione
#' @param random_seed Seed per riproducibilità
#' @return Matrice di espressione modificata con effetti di segnalazione
#' @importFrom stats plogis
#' @export
apply_lr_signaling <- function(
  expression_matrix,
  signaling_effects,
  integration_params = list(
    signaling_weight = 0.6,      # Peso degli effetti di segnalazione
    adjust_method = "adaptive",  # Metodo di integrazione: additive, multiplicative, o adaptive
    min_effect = 0.01,           # Soglia minima per applicare effetti
    max_effect = 3.0             # Limite massimo dell'effetto
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  signaling_weight <- integration_params$signaling_weight
  adjust_method <- integration_params$adjust_method
  min_effect <- integration_params$min_effect
  max_effect <- integration_params$max_effect
  
  # Valori predefiniti se non specificati
  if (is.null(signaling_weight)) signaling_weight <- 0.6
  if (is.null(adjust_method)) adjust_method <- "adaptive"
  if (is.null(min_effect)) min_effect <- 0.01
  if (is.null(max_effect)) max_effect <- 3.0
  
  # Limita gli effetti alla soglia specificata
  signaling_effects[signaling_effects < min_effect & signaling_effects > -min_effect] <- 0
  signaling_effects[signaling_effects > max_effect] <- max_effect
  signaling_effects[signaling_effects < -max_effect] <- -max_effect
  
  # Applica gli effetti in base al metodo selezionato
  if (adjust_method == "additive") {
    # Metodo additivo (baseline + effetto)
    adjusted_matrix <- expression_matrix + signaling_weight * signaling_effects
    
  } else if (adjust_method == "multiplicative") {
    # Metodo moltiplicativo (baseline * fattore)
    # Converte gli effetti in fattori moltiplicativi (0.5 = dimezzato, 2 = raddoppiato)
    factors <- exp(signaling_weight * signaling_effects)
    adjusted_matrix <- expression_matrix * factors
    
  } else if (adjust_method == "adaptive") {
    # Metodo adattivo: additivo per effetti piccoli, moltiplicativo per grandi
    # Per ogni cella ed elemento della matrice
    adjusted_matrix <- expression_matrix
    
    for (i in 1:nrow(expression_matrix)) {
      for (j in 1:ncol(expression_matrix)) {
        effect <- signaling_effects[i, j]
        if (abs(effect) >= min_effect) {
          base_expr <- expression_matrix[i, j]
          
          if (abs(effect) < 0.5) {
            # Effetto piccolo: additivo
            adj_expr <- base_expr + signaling_weight * effect
          } else {
            # Effetto grande: moltiplicativo
            factor <- exp(signaling_weight * effect)
            adj_expr <- base_expr * factor
          }
          
          # Assicura che l'espressione rimanga non-negativa
          adjusted_matrix[i, j] <- max(0, adj_expr)
        }
      }
    }
  } else {
    # Default: applicazione sigmoide (più biologicamente plausibile)
    # Trasforma gli effetti usando una funzione sigmoide
    sigmoid_effects <- 2 * plogis(signaling_weight * signaling_effects) - 1
    adjusted_matrix <- expression_matrix * (1 + sigmoid_effects)
  }
  
  # Assicura che tutti i valori siano non-negativi
  adjusted_matrix[adjusted_matrix < 0] <- 0
  
  return(adjusted_matrix)
}

#' Genera effetti completi di interazioni ligando-recettore
#'
#' Funzione wrapper che esegue l'intero processo di modellazione delle
#' interazioni ligando-recettore e la loro integrazione nell'espressione genica.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expression_matrix Matrice di espressione originale
#' @param dist_mat Matrice delle distanze tra cellule
#' @param lr_params Parametri per le interazioni L-R
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrice di espressione modificata e metadati
#' @export
generate_lr_interactions <- function(
  cell_df,
  expression_matrix,
  dist_mat,
  lr_params = list(
    n_interactions = 20,                 # Numero di interazioni L-R
    max_signaling_distance = 100,        # Distanza massima di segnalazione
    signal_propagation_mode = "distance_decay", # Modalità di propagazione
    signaling_weight = 0.6,              # Peso degli effetti di segnalazione
    min_target_genes = 3,                # Min geni target per interazione
    max_target_genes = 15,               # Max geni target per interazione
    effect_strength = c(0.5, 2.0),       # Range intensità effetto
    inhibitory_prob = 0.3,               # Prob. interazione inibitoria
    decayFactor = c(10, 50),             # Range fattore decadimento
    signal_amplification = 1.2,          # Amplificazione del segnale
    background_signaling = 0.05,         # Segnalazione di fondo
    adjust_method = "adaptive"           # Metodo integrazione
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Numero di geni
  n_genes <- ncol(expression_matrix)
  
  # 1. Definisci il database delle interazioni L-R
  interaction_params <- list(
    n_interactions = lr_params$n_interactions,
    min_target_genes = lr_params$min_target_genes,
    max_target_genes = lr_params$max_target_genes,
    effect_strength = lr_params$effect_strength,
    inhibitory_prob = lr_params$inhibitory_prob,
    decayFactor = lr_params$decayFactor
  )
  
  interaction_db <- define_ligand_receptor_db(
    n_genes = n_genes,
    interaction_params = interaction_params,
    random_seed = random_seed
  )
  
  # 2. Calcola gli effetti della segnalazione
  signaling_params <- list(
    signal_propagation_mode = lr_params$signal_propagation_mode,
    max_signaling_distance = lr_params$max_signaling_distance,
    signal_amplification = lr_params$signal_amplification,
    background_signaling = lr_params$background_signaling
  )
  
  signaling_effects <- compute_lr_signaling_effects(
    cell_df = cell_df,
    expression_matrix = expression_matrix,
    interaction_db = interaction_db,
    dist_mat = dist_mat,
    interaction_params = signaling_params,
    random_seed = random_seed
  )
  
  # 3. Applica gli effetti all'espressione
  integration_params <- list(
    signaling_weight = lr_params$signaling_weight,
    adjust_method = lr_params$adjust_method,
    min_effect = 0.01,
    max_effect = 3.0
  )
  
  adjusted_expression <- apply_lr_signaling(
    expression_matrix = expression_matrix,
    signaling_effects = signaling_effects,
    integration_params = integration_params,
    random_seed = random_seed
  )
  
  # Restituisci i risultati
  return(list(
    expression = adjusted_expression,
    original_expression = expression_matrix,
    signaling_effects = signaling_effects,
    interaction_db = interaction_db
  ))
}