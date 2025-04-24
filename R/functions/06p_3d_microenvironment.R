#' Effetti del microambiente tridimensionale
#'
#' Giustificazione: Lavori recenti (Dries et al., Nature Methods 2021; 
#' Littman et al., Nature Biotechnology 2021) evidenziano l'importanza della 
#' struttura 3D nei pattern di espressione, anche in sezioni 2D.
#'
#' Questo modulo implementa:
#' - Simulazione dell'effetto di proiezione 2D di strutture 3D (effetti di sovrapposizione)
#' - Modelli della profondità variabile di campionamento nelle tecnologie spaziali
#' - Incorporamento di effetti di prossimità 3D non catturati dalla distanza 2D

#' Genera livelli di profondità per simulare la terza dimensione
#'
#' Crea un campo di profondità che rappresenta la posizione delle celle
#' nella terza dimensione (z), per simulare effetti dovuti a proiezioni 2D.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param depth_params Parametri per la generazione della profondità
#' @param random_seed Seed per riproducibilità
#' @return Vettore di valori di profondità per ogni cellula
#' @importFrom stats runif rnorm
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_depth_field <- function(
  cell_df,
  depth_params = list(
    depth_pattern = "smooth",     # "smooth", "layered", "complex" o "random"
    depth_range = c(-10, 10),     # Range di profondità (μm)
    n_layers = 3,                 # Numero di strati per modalità "layered"
    layer_fuzziness = 0.3,        # Transizione graduale tra strati (0-1)
    spatial_coherence = 0.8,      # Coerenza spaziale della profondità (0-1)
    depth_noise = 0.2,            # Rumore casuale nella profondità
    z_resolution = 1              # Risoluzione in z (μm) per discretizzazione
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  depth_pattern <- depth_params$depth_pattern
  depth_range <- depth_params$depth_range
  n_layers <- depth_params$n_layers
  layer_fuzziness <- depth_params$layer_fuzziness
  spatial_coherence <- depth_params$spatial_coherence
  depth_noise <- depth_params$depth_noise
  z_resolution <- depth_params$z_resolution
  
  # Valori predefiniti se non specificati
  if (is.null(depth_pattern)) depth_pattern <- "smooth"
  if (is.null(depth_range)) depth_range <- c(-10, 10)
  if (is.null(n_layers)) n_layers <- 3
  if (is.null(layer_fuzziness)) layer_fuzziness <- 0.3
  if (is.null(spatial_coherence)) spatial_coherence <- 0.8
  if (is.null(depth_noise)) depth_noise <- 0.2
  if (is.null(z_resolution)) z_resolution <- 1
  
  # Numero di cellule
  n_cells <- nrow(cell_df)
  
  # Inizializza campo di profondità
  depth_field <- numeric(n_cells)
  
  # Genera campi di profondità in base al pattern specificato
  if (depth_pattern == "smooth") {
    # Campo di profondità spazialmente coerente
    # Converti cell_df in oggetto spatial
    sp_df <- cell_df
    sp::coordinates(sp_df) <- ~ x + y
    
    # Crea un campo gaussiano spazialmente correlato
    range_param <- mean(c(diff(range(cell_df$x)), diff(range(cell_df$y)))) * 0.3
    
    gp_sim <- gstat::gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                    beta = 0, model = gstat::vgm(psill = 1.0, 
                                           range = range_param, 
                                           model = "Exp"), 
                    nmax = 30)
    
    # Predizione del campo gaussiano
    depth_field <- gstat::predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
    
    # Normalizza e scala al range specificato
    depth_field <- scale01(depth_field)
    depth_field <- depth_range[1] + depth_field * (depth_range[2] - depth_range[1])
    
  } else if (depth_pattern == "layered") {
    # Struttura a strati con transizioni
    # Genera campo di base spazialmente correlato per determinare l'appartenenza a strati
    sp_df <- cell_df
    sp::coordinates(sp_df) <- ~ x + y
    
    # Campo di base per determinare l'appartenenza a strati
    gp_sim <- gstat::gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                    beta = 0, model = gstat::vgm(psill = 1.0, 
                                           range = mean(c(diff(range(cell_df$x)), 
                                                       diff(range(cell_df$y)))) * 0.3, 
                                           model = "Exp"), 
                    nmax = 30)
    
    base_field <- gstat::predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
    
    # Normalizza
    base_field <- scale01(base_field)
    
    # Genera valori discreti per gli strati
    layer_assignments <- cut(base_field, breaks = n_layers, 
                            labels = FALSE, include.lowest = TRUE)
    
    # Valore medio di ogni strato
    layer_values <- seq(depth_range[1], depth_range[2], length.out = n_layers)
    
    # Assegna profondità in base allo strato, con transizione graduale se richiesto
    if (layer_fuzziness > 0) {
      # Calcola la distanza dal centro dello strato per ogni cella
      for (i in 1:n_cells) {
        layer <- layer_assignments[i]
        
        # Distanza normalizzata dal centro dello strato (0 = centro, 1 = confine)
        layer_pos <- base_field[i] * n_layers - layer + 0.5
        dist_from_center <- abs(layer_pos - 0.5) * 2  # 0 al centro, 1 ai confini
        
        # Aggiungi rumore proporzionale alla distanza dal centro e alla fuzziness
        layer_noise <- rnorm(1, 0, layer_fuzziness * dist_from_center)
        
        # Interpola tra strati adiacenti se siamo vicini ai confini
        if (dist_from_center > 0.8) {
          # Determina strato adiacente
          if (layer_pos > 0.5) {
            next_layer <- min(layer + 1, n_layers)
          } else {
            next_layer <- max(layer - 1, 1)
          }
          
          # Calcola peso per interpolazione (0.5-1.0 in base a quanto siamo vicini al confine)
          blend_weight <- (dist_from_center - 0.8) / 0.2 * layer_fuzziness
          
          # Interpola
          depth_field[i] <- (1 - blend_weight) * layer_values[layer] + 
                          blend_weight * layer_values[next_layer] + 
                          layer_noise * (depth_range[2] - depth_range[1]) / n_layers
        } else {
          depth_field[i] <- layer_values[layer] + 
                         layer_noise * (depth_range[2] - depth_range[1]) / n_layers
        }
      }
    } else {
      # Assegnazione diretta senza fuzzyness
      for (i in 1:n_cells) {
        depth_field[i] <- layer_values[layer_assignments[i]]
      }
    }
    
  } else if (depth_pattern == "complex") {
    # Combinazione di pattern diversi per strutture tissutali complesse
    # Genera più campi e li combina
    
    # Campo base a larga scala
    sp_df <- cell_df
    sp::coordinates(sp_df) <- ~ x + y
    
    # Parametri del tessuto
    width <- diff(range(cell_df$x))
    height <- diff(range(cell_df$y))
    
    # Campo 1: Pattern a larga scala (es. curvatura globale del tessuto)
    gp_sim1 <- gstat::gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                     beta = 0, model = gstat::vgm(psill = 1.0, 
                                            range = max(width, height) * 0.7, 
                                            model = "Sph"), 
                     nmax = 40)
    
    field1 <- gstat::predict(gp_sim1, newdata = sp_df, nsim = 1)$sim1
    
    # Campo 2: Strutture a media scala (es. domini tissutali)
    gp_sim2 <- gstat::gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                     beta = 0, model = gstat::vgm(psill = 1.0, 
                                            range = mean(c(width, height)) * 0.2, 
                                            model = "Exp"), 
                     nmax = 30)
    
    field2 <- gstat::predict(gp_sim2, newdata = sp_df, nsim = 1)$sim1
    
    # Campo 3: Dettagli a piccola scala (es. variazioni locali)
    gp_sim3 <- gstat::gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                     beta = 0, model = gstat::vgm(psill = 1.0, 
                                            range = min(width, height) * 0.05, 
                                            model = "Gau"), 
                     nmax = 20)
    
    field3 <- gstat::predict(gp_sim3, newdata = sp_df, nsim = 1)$sim1
    
    # Combina i campi con pesi differenti
    depth_field <- 0.5 * scale01(field1) + 0.3 * scale01(field2) + 0.2 * scale01(field3)
    
    # Normalizza e scala
    depth_field <- scale01(depth_field)
    depth_field <- depth_range[1] + depth_field * (depth_range[2] - depth_range[1])
    
  } else {
    # Default: random con coerenza spaziale
    random_depth <- rnorm(n_cells, mean = mean(depth_range), sd = diff(depth_range) / 4)
    
    if (spatial_coherence > 0) {
      # Aggiungi coerenza spaziale
      sp_df <- cell_df
      sp::coordinates(sp_df) <- ~ x + y
      
      # Crea un campo gaussiano spazialmente correlato
      range_param <- mean(c(diff(range(cell_df$x)), diff(range(cell_df$y)))) * 0.2
      
      gp_sim <- gstat::gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                      beta = 0, model = gstat::vgm(psill = 1.0, 
                                             range = range_param, 
                                             model = "Exp"), 
                      nmax = 30)
      
      spatial_comp <- gstat::predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
      
      # Normalizza
      spatial_comp <- scale01(spatial_comp)
      
      # Combina componente casuale con componente spaziale
      depth_field <- (1 - spatial_coherence) * random_depth + 
                   spatial_coherence * (depth_range[1] + spatial_comp * (depth_range[2] - depth_range[1]))
    } else {
      depth_field <- random_depth
    }
  }
  
  # Aggiungi rumore casuale
  if (depth_noise > 0) {
    noise_scale <- diff(depth_range) * depth_noise
    depth_field <- depth_field + rnorm(n_cells, 0, noise_scale)
  }
  
  # Limita al range specificato
  depth_field <- pmax(depth_range[1], pmin(depth_range[2], depth_field))
  
  # Discretizza in base alla risoluzione in z se specificato
  if (z_resolution > 0) {
    depth_field <- round(depth_field / z_resolution) * z_resolution
  }
  
  return(depth_field)
}

#' Funzione di utilità per normalizzare vettori o matrici nell'intervallo [0,1]
#'
#' @param x Vettore o matrice da normalizzare
#' @return Vettore o matrice normalizzato
#' @keywords internal
scale01 <- function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    return(rep(0.5, length(x)))
  }
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

#' Calcola distanze 3D tra cellule
#'
#' Genera una matrice di distanze 3D tra cellule, utilizzando i valori
#' di profondità per la coordinata z.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param depth_field Vettore di valori di profondità
#' @param dist_mat Matrice di distanze 2D esistente (opzionale)
#' @param dist_params Parametri per il calcolo delle distanze
#' @param random_seed Seed per riproducibilità
#' @return Matrice di distanze 3D
#' @export
calculate_3d_distances <- function(
  cell_df,
  depth_field,
  dist_mat = NULL,
  dist_params = list(
    distance_metric = "euclidean",  # "euclidean", "manhattan" o "weighted"
    z_weight = 0.7,                # Peso della dimensione z (per distance_metric="weighted")
    use_sparse = TRUE,             # Utilizza matrice sparsa per grandi dataset
    max_distance = NULL            # Distanza massima da memorizzare
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  distance_metric <- dist_params$distance_metric
  z_weight <- dist_params$z_weight
  use_sparse <- dist_params$use_sparse
  max_distance <- dist_params$max_distance
  
  # Valori predefiniti se non specificati
  if (is.null(distance_metric)) distance_metric <- "euclidean"
  if (is.null(z_weight)) z_weight <- 0.7
  if (is.null(use_sparse)) use_sparse <- TRUE
  
  # Numero di cellule
  n_cells <- nrow(cell_df)
  
  # Se la matrice di distanze 2D è già calcolata, la riutilizziamo
  if (is.null(dist_mat)) {
    # Calcola matrice di distanze 2D
    dist_mat_2d <- matrix(0, nrow = n_cells, ncol = n_cells)
    
    # Calcola distanze euclidee 2D
    for (i in 1:n_cells) {
      for (j in i:n_cells) {
        if (i == j) {
          dist_mat_2d[i, j] <- 0
        } else {
          dist_2d <- sqrt((cell_df$x[i] - cell_df$x[j])^2 + 
                         (cell_df$y[i] - cell_df$y[j])^2)
          dist_mat_2d[i, j] <- dist_mat_2d[j, i] <- dist_2d
        }
      }
    }
  } else {
    # Usa la matrice esistente
    dist_mat_2d <- dist_mat
  }
  
  # Calcola matrice di distanze 3D
  if (distance_metric == "euclidean") {
    # Distanza euclidea 3D standard
    dist_mat_3d <- matrix(0, nrow = n_cells, ncol = n_cells)
    
    for (i in 1:n_cells) {
      for (j in i:n_cells) {
        if (i == j) {
          dist_mat_3d[i, j] <- 0
        } else {
          z_diff <- depth_field[i] - depth_field[j]
          dist_3d <- sqrt(dist_mat_2d[i, j]^2 + z_diff^2)
          dist_mat_3d[i, j] <- dist_mat_3d[j, i] <- dist_3d
        }
      }
    }
  } else if (distance_metric == "manhattan") {
    # Distanza Manhattan 3D (somma delle differenze assolute)
    dist_mat_3d <- matrix(0, nrow = n_cells, ncol = n_cells)
    
    # Per Manhattan, abbiamo bisogno di scomporre le distanze 2D
    for (i in 1:n_cells) {
      for (j in i:n_cells) {
        if (i == j) {
          dist_mat_3d[i, j] <- 0
        } else {
          x_diff <- abs(cell_df$x[i] - cell_df$x[j])
          y_diff <- abs(cell_df$y[i] - cell_df$y[j])
          z_diff <- abs(depth_field[i] - depth_field[j])
          dist_3d <- x_diff + y_diff + z_diff
          dist_mat_3d[i, j] <- dist_mat_3d[j, i] <- dist_3d
        }
      }
    }
  } else if (distance_metric == "weighted") {
    # Distanza euclidea con peso diverso per la dimensione z
    dist_mat_3d <- matrix(0, nrow = n_cells, ncol = n_cells)
    
    for (i in 1:n_cells) {
      for (j in i:n_cells) {
        if (i == j) {
          dist_mat_3d[i, j] <- 0
        } else {
          z_diff <- depth_field[i] - depth_field[j]
          dist_3d <- sqrt(dist_mat_2d[i, j]^2 + (z_weight * z_diff)^2)
          dist_mat_3d[i, j] <- dist_mat_3d[j, i] <- dist_3d
        }
      }
    }
  } else {
    # Default: euclidea
    dist_mat_3d <- matrix(0, nrow = n_cells, ncol = n_cells)
    
    for (i in 1:n_cells) {
      for (j in i:n_cells) {
        if (i == j) {
          dist_mat_3d[i, j] <- 0
        } else {
          z_diff <- depth_field[i] - depth_field[j]
          dist_3d <- sqrt(dist_mat_2d[i, j]^2 + z_diff^2)
          dist_mat_3d[i, j] <- dist_mat_3d[j, i] <- dist_3d
        }
      }
    }
  }
  
  # Applica la distanza massima se specificata
  if (!is.null(max_distance)) {
    dist_mat_3d[dist_mat_3d > max_distance] <- Inf
  }
  
  # Restituisci matrice sparsa se richiesto e disponibile
  if (use_sparse && requireNamespace("Matrix", quietly = TRUE) && n_cells > 1000) {
    result <- Matrix::Matrix(dist_mat_3d, sparse = TRUE)
  } else {
    result <- dist_mat_3d
  }
  
  return(result)
}

#' Simula effetti di sovrapposizione cellulare
#'
#' Modella l'effetto della proiezione di dati 3D su una superficie 2D,
#' che causa sovrapposizioni di cellule ed effetti di contaminazione.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param depth_field Vettore di valori di profondità
#' @param expression_matrix Matrice di espressione originale
#' @param overlap_params Parametri per gli effetti di sovrapposizione
#' @param random_seed Seed per riproducibilità
#' @return Matrice di espressione modificata con effetti di sovrapposizione
#' @export
generate_overlap_effects <- function(
  cell_df,
  depth_field,
  expression_matrix,
  overlap_params = list(
    overlap_radius = 5,           # Raggio di sovrapposizione in μm
    depth_threshold = 3,          # Differenza minima di profondità per considerare sovrapposizione
    max_overlap_effect = 0.3,     # Effetto massimo di sovrapposizione (0-1)
    effect_decay = "exponential", # "linear", "exponential" o "threshold"
    depth_visibility_mask = NULL, # Maschera di visibilità per limitare gli effetti alle celle visibili
    weight_by_area = 0.5          # Ponderazione per area di sovrapposizione
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  overlap_radius <- overlap_params$overlap_radius
  depth_threshold <- overlap_params$depth_threshold
  max_overlap_effect <- overlap_params$max_overlap_effect
  effect_decay <- overlap_params$effect_decay
  depth_visibility_mask <- overlap_params$depth_visibility_mask
  weight_by_area <- overlap_params$weight_by_area
  
  # Valori predefiniti se non specificati
  if (is.null(overlap_radius)) overlap_radius <- 5
  if (is.null(depth_threshold)) depth_threshold <- 3
  if (is.null(max_overlap_effect)) max_overlap_effect <- 0.3
  if (is.null(effect_decay)) effect_decay <- "exponential"
  if (is.null(weight_by_area)) weight_by_area <- 0.5
  
  # Numero di cellule e geni
  n_cells <- nrow(cell_df)
  n_genes <- ncol(expression_matrix)
  
  # Inizializza la matrice per gli effetti di sovrapposizione
  overlap_contribution <- matrix(0, nrow = n_cells, ncol = n_genes)
  overlap_weights <- numeric(n_cells)
  
  # Calcola le coordinate 2D (x, y) e la profondità (z)
  coords_2d <- cbind(cell_df$x, cell_df$y)
  
  # Se non è specificata una maschera di visibilità, tutte le celle sono considerate visibili
  if (is.null(depth_visibility_mask)) {
    visibility_mask <- rep(TRUE, n_cells)
  } else {
    visibility_mask <- depth_visibility_mask
  }
  
  # Per ogni cellula, calcola l'effetto di sovrapposizione
  for (i in 1:n_cells) {
    # Salta se la cella non è visibile
    if (!visibility_mask[i]) next
    
    # Trova celle vicine in 2D
    for (j in 1:n_cells) {
      # Una cella non si sovrappone a se stessa
      if (i == j) next
      
      # Calcola distanza 2D e differenza di profondità
      dist_2d <- sqrt(sum((coords_2d[i, ] - coords_2d[j, ])^2))
      z_diff <- depth_field[j] - depth_field[i]
      
      # Consideriamo solo celle che sono:
      # 1. Vicine in 2D (entro il raggio di sovrapposizione)
      # 2. Più in alto (z maggiore) rispetto alla cella corrente
      # 3. Con differenza di profondità superiore alla soglia
      if (dist_2d <= overlap_radius && z_diff >= depth_threshold) {
        # Calcola l'effetto di sovrapposizione basato sulla distanza 2D e la differenza di profondità
        
        if (effect_decay == "linear") {
          # Decadimento lineare con la distanza
          distance_effect <- 1 - (dist_2d / overlap_radius)
          
        } else if (effect_decay == "exponential") {
          # Decadimento esponenziale
          distance_effect <- exp(-dist_2d / (overlap_radius / 2))
          
        } else if (effect_decay == "threshold") {
          # Effetto a soglia (costante entro il raggio)
          distance_effect <- 1.0
          
        } else {
          # Default: decadimento gaussiano
          distance_effect <- exp(-(dist_2d^2) / (2 * (overlap_radius/2)^2))
        }
        
        # Calcola l'effetto di sovrapposizione
        if (weight_by_area > 0) {
          # Calcoliamo l'area di sovrapposizione (approssimata come cerchio)
          # L'area di sovrapposizione è approssimata dalla formula di intersezione di due cerchi
          d <- dist_2d
          R <- overlap_radius
          
          if (d < 2 * R) {
            # Area di intersezione di due cerchi di raggio R a distanza d
            area_overlap <- 2 * R^2 * acos(d / (2 * R)) - 0.5 * d * sqrt(4 * R^2 - d^2)
            area_max <- pi * R^2  # Area massima (cerchio completo)
            
            # Normalizza tra 0 e 1
            area_effect <- area_overlap / area_max
          } else {
            area_effect <- 0
          }
          
          # Combina effetto distanza e area
          overlap_effect <- (1 - weight_by_area) * distance_effect + weight_by_area * area_effect
        } else {
          overlap_effect <- distance_effect
        }
        
        # Modula l'effetto in base alla differenza di profondità
        depth_effect <- min(1, z_diff / (2 * depth_threshold))
        
        # L'effetto finale è proporzionale sia alla distanza 2D che alla differenza di profondità
        final_effect <- max_overlap_effect * overlap_effect * depth_effect
        
        # Aggiungi il contributo di questa cella sovrapposta
        overlap_contribution[i, ] <- overlap_contribution[i, ] + final_effect * expression_matrix[j, ]
        overlap_weights[i] <- overlap_weights[i] + final_effect
      }
    }
  }
  
  # Calcola la matrice finale con gli effetti di sovrapposizione
  # Per ogni cella, combiniamo l'espressione originale con il contributo di sovrapposizione
  modified_expression <- expression_matrix + overlap_contribution
  
  # Normalizziamo per evitare che l'effetto di sovrapposizione aumenti eccessivamente l'espressione totale
  # Ma solo per celle che hanno un effetto di sovrapposizione
  for (i in 1:n_cells) {
    if (overlap_weights[i] > 0) {
      normalization_factor <- 1 + 0.5 * overlap_weights[i]  # Fattore di normalizzazione (aggiustabile)
      modified_expression[i, ] <- modified_expression[i, ] / normalization_factor
    }
  }
  
  return(modified_expression)
}

#' Genera un modello di proiezioni 3D complete
#'
#' Funzione wrapper che implementa l'intero modello di effetti di microambiente 3D,
#' inclusi profondità variabile, distanze 3D, e effetti di sovrapposizione.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expression_matrix Matrice di espressione originale
#' @param dist_mat Matrice di distanze 2D (opzionale)
#' @param microenv_params Parametri per il microambiente 3D
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici di espressione e metadati 3D
#' @export
generate_3d_microenvironment <- function(
  cell_df,
  expression_matrix,
  dist_mat = NULL,
  microenv_params = list(
    # Parametri per la generazione della profondità
    depth_pattern = "smooth",
    depth_range = c(-10, 10),
    spatial_coherence = 0.8,
    depth_noise = 0.2,
    
    # Parametri per le distanze 3D
    use_3d_distances = TRUE,
    distance_metric = "weighted",
    z_weight = 0.7,
    
    # Parametri per gli effetti di sovrapposizione
    model_overlap = TRUE,
    overlap_radius = 5,
    depth_threshold = 3,
    max_overlap_effect = 0.3,
    effect_decay = "exponential"
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  depth_pattern <- microenv_params$depth_pattern
  depth_range <- microenv_params$depth_range
  spatial_coherence <- microenv_params$spatial_coherence
  depth_noise <- microenv_params$depth_noise
  
  use_3d_distances <- microenv_params$use_3d_distances
  distance_metric <- microenv_params$distance_metric
  z_weight <- microenv_params$z_weight
  
  model_overlap <- microenv_params$model_overlap
  overlap_radius <- microenv_params$overlap_radius
  depth_threshold <- microenv_params$depth_threshold
  max_overlap_effect <- microenv_params$max_overlap_effect
  effect_decay <- microenv_params$effect_decay
  
  # Valori predefiniti se non specificati
  if (is.null(depth_pattern)) depth_pattern <- "smooth"
  if (is.null(depth_range)) depth_range <- c(-10, 10)
  if (is.null(spatial_coherence)) spatial_coherence <- 0.8
  if (is.null(depth_noise)) depth_noise <- 0.2
  
  if (is.null(use_3d_distances)) use_3d_distances <- TRUE
  if (is.null(distance_metric)) distance_metric <- "weighted"
  if (is.null(z_weight)) z_weight <- 0.7
  
  if (is.null(model_overlap)) model_overlap <- TRUE
  if (is.null(overlap_radius)) overlap_radius <- 5
  if (is.null(depth_threshold)) depth_threshold <- 3
  if (is.null(max_overlap_effect)) max_overlap_effect <- 0.3
  if (is.null(effect_decay)) effect_decay <- "exponential"
  
  # 1. Genera campo di profondità
  depth_params <- list(
    depth_pattern = depth_pattern,
    depth_range = depth_range,
    spatial_coherence = spatial_coherence,
    depth_noise = depth_noise
  )
  
  depth_field <- generate_depth_field(
    cell_df = cell_df,
    depth_params = depth_params,
    random_seed = random_seed
  )
  
  # 2. Calcola distanze 3D se richiesto
  dist_mat_3d <- NULL
  if (use_3d_distances) {
    dist_params <- list(
      distance_metric = distance_metric,
      z_weight = z_weight,
      use_sparse = TRUE
    )
    
    dist_mat_3d <- calculate_3d_distances(
      cell_df = cell_df,
      depth_field = depth_field,
      dist_mat = dist_mat,
      dist_params = dist_params,
      random_seed = random_seed
    )
  }
  
  # 3. Genera effetti di sovrapposizione se richiesto
  modified_expression <- expression_matrix
  if (model_overlap) {
    overlap_params <- list(
      overlap_radius = overlap_radius,
      depth_threshold = depth_threshold,
      max_overlap_effect = max_overlap_effect,
      effect_decay = effect_decay
    )
    
    modified_expression <- generate_overlap_effects(
      cell_df = cell_df,
      depth_field = depth_field,
      expression_matrix = expression_matrix,
      overlap_params = overlap_params,
      random_seed = random_seed
    )
  }
  
  # 4. Arricchisci cell_df con i valori di profondità
  cell_df_3d <- cell_df
  cell_df_3d$depth <- depth_field
  
  # Restituisci i risultati
  return(list(
    expression = modified_expression,      # Matrice di espressione modificata
    original_expression = expression_matrix, # Matrice originale
    cell_df_3d = cell_df_3d,               # Dataframe celle con profondità
    depth_field = depth_field,             # Campo di profondità
    dist_mat_3d = dist_mat_3d              # Matrice distanze 3D
  ))
}