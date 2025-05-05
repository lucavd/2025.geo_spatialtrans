#' Modellazione semplificata di pattern anisotropici per testing
#'
#' Versioni semplificate delle funzioni di pattern anisotropici adattate per testing
#' 
#' @import stats

#' Genera strutture backbone (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param structure_params Parametri per le strutture
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici di distanza e maschera delle strutture
#' @export
generate_backbone_structures <- function(
  cell_df,
  structure_params = list(
    n_structures = 1,
    structure_type = "linear",
    structure_width = 2,
    structure_length_factor = 0.8,
    n_branches = 3,
    network_density = 0.1
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  n_structures <- structure_params$n_structures
  structure_type <- structure_params$structure_type
  width <- structure_params$structure_width
  length_factor <- structure_params$structure_length_factor
  n_branches <- structure_params$n_branches
  network_density <- structure_params$network_density
  
  # Valori predefiniti se non specificati
  if (is.null(n_structures)) n_structures <- 1
  if (is.null(structure_type)) structure_type <- "linear"
  if (is.null(width)) width <- 2
  if (is.null(length_factor)) length_factor <- 0.8
  if (is.null(n_branches)) n_branches <- 3
  if (is.null(network_density)) network_density <- 0.1
  
  # Estrai dimensioni
  n_cells <- nrow(cell_df)
  
  # Range della griglia spaziale
  x_range <- range(cell_df$x)
  y_range <- range(cell_df$y)
  grid_width <- x_range[2] - x_range[1]
  grid_height <- y_range[2] - y_range[1]
  
  # Converti width da unità relative a assolute
  width_abs <- width * min(grid_width, grid_height) / 20
  
  # Inizializza liste e vettori per risultati
  distance_matrices <- list()
  structure_mask <- rep(0, n_cells)
  
  # Genera ogni struttura
  for (s in 1:n_structures) {
    # Crea una matrice di distanza per questa struttura
    dist_mat <- matrix(Inf, nrow = n_cells, ncol = n_cells)
    diag(dist_mat) <- 0
    
    # Vettore di appartenenza delle cellule a questa struttura
    structure_cells <- rep(FALSE, n_cells)
    
    if (structure_type == "linear") {
      # Genera una struttura lineare
      # Punto di inizio e direzione
      start_x <- runif(1, x_range[1], x_range[2])
      start_y <- runif(1, y_range[1], y_range[2])
      angle <- runif(1, 0, 2*pi)
      
      # Lunghezza come frazione della diagonale
      diagonal <- sqrt(grid_width^2 + grid_height^2)
      length <- diagonal * length_factor
      
      # Punto finale
      end_x <- start_x + cos(angle) * length
      end_y <- start_y + sin(angle) * length
      
      # Per ogni cellula, calcola la distanza dalla linea
      for (i in 1:n_cells) {
        cell_x <- cell_df$x[i]
        cell_y <- cell_df$y[i]
        
        # Distanza al segmento di linea
        dist_to_line <- distanceToSegment(cell_x, cell_y, start_x, start_y, end_x, end_y)
        
        # Se la cellula è vicina alla linea, è parte della struttura
        if (dist_to_line <= width_abs) {
          structure_cells[i] <- TRUE
          
          # Calcola la posizione relativa lungo la linea (parametro t)
          t <- projectPointOnLine(cell_x, cell_y, start_x, start_y, end_x, end_y)
          
          # Per ogni altra cellula nella struttura, calcola distanza lungo la linea
          for (j in 1:n_cells) {
            if (structure_cells[j]) {
              cell_x2 <- cell_df$x[j]
              cell_y2 <- cell_df$y[j]
              t2 <- projectPointOnLine(cell_x2, cell_y2, start_x, start_y, end_x, end_y)
              
              # Distanza lungo la struttura
              line_dist <- abs(t - t2) * length
              
              # Aggiorna la matrice di distanza
              dist_mat[i, j] <- line_dist
              dist_mat[j, i] <- line_dist
            }
          }
        }
      }
      
    } else if (structure_type == "branched") {
      # Genera una struttura ramificata
      # Punto centrale
      center_x <- runif(1, x_range[1] + grid_width*0.2, x_range[2] - grid_width*0.2)
      center_y <- runif(1, y_range[1] + grid_height*0.2, y_range[2] - grid_height*0.2)
      
      # Punti finali dei rami
      branch_ends_x <- c()
      branch_ends_y <- c()
      
      diagonal <- sqrt(grid_width^2 + grid_height^2)
      branch_length <- diagonal * length_factor * 0.6  # Rami più corti della struttura lineare
      
      # Genera rami in direzioni casuali
      branch_angles <- seq(0, 2*pi * (1 - 1/n_branches), length.out = n_branches)
      branch_angles <- branch_angles + runif(1, 0, 2*pi)  # Rotazione casuale
      
      for (b in 1:n_branches) {
        angle <- branch_angles[b]
        end_x <- center_x + cos(angle) * branch_length
        end_y <- center_y + sin(angle) * branch_length
        branch_ends_x <- c(branch_ends_x, end_x)
        branch_ends_y <- c(branch_ends_y, end_y)
      }
      
      # Per ogni cellula, calcola distanza minima da un ramo
      for (i in 1:n_cells) {
        cell_x <- cell_df$x[i]
        cell_y <- cell_df$y[i]
        
        min_dist <- Inf
        closest_branch <- 0
        closest_t <- 0
        
        # Trova il ramo più vicino
        for (b in 1:n_branches) {
          dist_to_branch <- distanceToSegment(cell_x, cell_y, center_x, center_y, 
                                             branch_ends_x[b], branch_ends_y[b])
          
          if (dist_to_branch < min_dist) {
            min_dist <- dist_to_branch
            closest_branch <- b
            closest_t <- projectPointOnLine(cell_x, cell_y, center_x, center_y, 
                                           branch_ends_x[b], branch_ends_y[b])
          }
        }
        
        # Se la cellula è vicina a un ramo, è parte della struttura
        if (min_dist <= width_abs) {
          structure_cells[i] <- TRUE
          
          # Per calcolare distanze tra cellule, considera il percorso attraverso il centro
          # se le cellule sono su rami diversi
          for (j in 1:n_cells) {
            if (structure_cells[j]) {
              cell_x2 <- cell_df$x[j]
              cell_y2 <- cell_df$y[j]
              
              min_dist2 <- Inf
              closest_branch2 <- 0
              closest_t2 <- 0
              
              for (b in 1:n_branches) {
                dist_to_branch <- distanceToSegment(cell_x2, cell_y2, center_x, center_y, 
                                                  branch_ends_x[b], branch_ends_y[b])
                
                if (dist_to_branch < min_dist2) {
                  min_dist2 <- dist_to_branch
                  closest_branch2 <- b
                  closest_t2 <- projectPointOnLine(cell_x2, cell_y2, center_x, center_y, 
                                                 branch_ends_x[b], branch_ends_y[b])
                }
              }
              
              # Calcola distanza
              if (closest_branch == closest_branch2) {
                # Stessa ramificazione, distanza diretta lungo il ramo
                branch_length_i <- branch_length * closest_t
                branch_length_j <- branch_length * closest_t2
                branch_dist <- abs(branch_length_i - branch_length_j)
                
                dist_mat[i, j] <- branch_dist
                dist_mat[j, i] <- branch_dist
              } else {
                # Ramificazioni diverse, distanza attraverso il centro
                branch_length_i <- branch_length * closest_t
                branch_length_j <- branch_length * closest_t2
                
                # Distanza da i al centro più distanza dal centro a j
                branch_dist <- branch_length_i + branch_length_j
                
                dist_mat[i, j] <- branch_dist
                dist_mat[j, i] <- branch_dist
              }
            }
          }
        }
      }
      
    } else if (structure_type == "network") {
      # Genera una rete di strutture interconnesse
      # Punti nodali casuali
      n_nodes <- max(3, round(n_cells * network_density * 0.05))
      node_x <- runif(n_nodes, x_range[1], x_range[2])
      node_y <- runif(n_nodes, y_range[1], y_range[2])
      
      # Costruisci un grafo connesso
      edges <- c()
      
      # Minimo spanning tree per garantire connessione
      mst <- buildMST(node_x, node_y)
      edges <- mst
      
      # Aggiungi archi casuali addizionali
      n_extra_edges <- round(n_nodes * network_density)
      for (e in 1:n_extra_edges) {
        i <- sample(1:n_nodes, 1)
        j <- sample(1:n_nodes, 1)
        if (i != j) {
          edges <- rbind(edges, c(i, j))
        }
      }
      
      # Per ogni cellula, trova la distanza minima da un arco
      for (i in 1:n_cells) {
        cell_x <- cell_df$x[i]
        cell_y <- cell_df$y[i]
        
        min_dist <- Inf
        closest_edge <- c(0, 0)
        closest_t <- 0
        
        for (e in 1:nrow(edges)) {
          n1 <- edges[e, 1]
          n2 <- edges[e, 2]
          
          dist_to_edge <- distanceToSegment(cell_x, cell_y, 
                                          node_x[n1], node_y[n1], 
                                          node_x[n2], node_y[n2])
          
          if (dist_to_edge < min_dist) {
            min_dist <- dist_to_edge
            closest_edge <- c(n1, n2)
            closest_t <- projectPointOnLine(cell_x, cell_y, 
                                          node_x[n1], node_y[n1], 
                                          node_x[n2], node_y[n2])
          }
        }
        
        # Se la cellula è vicina a un arco, è parte della struttura
        if (min_dist <= width_abs) {
          structure_cells[i] <- TRUE
          
          # Salva nodo più vicino per calcoli di distanza
          # Inizializza l'attributo node_info se non esiste ancora
          if (is.null(attr(structure_cells, "node_info"))) {
            attr(structure_cells, "node_info") <- vector("list", n_cells)
          }
          attr(structure_cells, "node_info")[[i]] <- list(
            edge = closest_edge,
            t = closest_t
          )
          
          # Le distanze sulla rete richiederebbero un algoritmo di percorso minimo completo
          # Per semplicità, usiamo distanze euclidee ma con un fattore moltiplicativo
          # che riflette la natura tortuosa della rete
          for (j in 1:n_cells) {
            if (structure_cells[j] && i != j) {
              euclidean_dist <- sqrt((cell_df$x[i] - cell_df$x[j])^2 + 
                                     (cell_df$y[i] - cell_df$y[j])^2)
              
              # Fattore di tortuosità: distanza maggiore lungo la rete rispetto alla linea retta
              network_factor <- 1.5
              dist_mat[i, j] <- euclidean_dist * network_factor
              dist_mat[j, i] <- euclidean_dist * network_factor
            }
          }
        }
      }
    }
    
    # Aggiorna maschera globale delle strutture
    structure_mask[structure_cells] <- s
    
    # Aggiungi matrice di distanza alla lista
    distance_matrices[[s]] <- dist_mat
  }
  
  return(list(
    distance_matrices = distance_matrices,
    structure_mask = structure_mask
  ))
}

#' Genera espressione genica anisotropica (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param distance_matrices Lista di matrici di distanza per strutture
#' @param structure_mask Vettore di appartenenza delle cellule alle strutture
#' @param n_genes Numero di geni da modellare
#' @param expression_params Parametri per l'espressione anisotropica
#' @param random_seed Seed per riproducibilità
#' @return Matrice di pattern di espressione anisotropica
#' @export
generate_anisotropic_expression <- function(
  cell_df,
  distance_matrices,
  structure_mask,
  n_genes = NULL,
  expression_params = list(
    anisotropic_pattern = "gradient", 
    anisotropic_gene_fraction = 0.5,
    pattern_smoothness = 0.2,
    oscillation_frequency = 0.2,
    n_hotspots_per_structure = 2,
    hotspot_radius = 5
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  pattern_type <- expression_params$anisotropic_pattern
  gene_fraction <- expression_params$anisotropic_gene_fraction
  smoothness <- expression_params$pattern_smoothness
  oscillation_freq <- expression_params$oscillation_frequency
  n_hotspots <- expression_params$n_hotspots_per_structure
  hotspot_radius <- expression_params$hotspot_radius
  
  # Valori predefiniti se non specificati
  if (is.null(pattern_type)) pattern_type <- "gradient"
  if (is.null(gene_fraction)) gene_fraction <- 0.5
  if (is.null(smoothness)) smoothness <- 0.2
  if (is.null(oscillation_freq)) oscillation_freq <- 0.2
  if (is.null(n_hotspots)) n_hotspots <- 2
  if (is.null(hotspot_radius)) hotspot_radius <- 5
  
  # Estrai dimensioni
  n_cells <- nrow(cell_df)
  n_structures <- length(distance_matrices)
  
  # Seleziona i geni che seguiranno pattern anisotropici
  n_aniso_genes <- round(n_genes * gene_fraction)
  aniso_genes <- sample(1:n_genes, n_aniso_genes)
  
  # Inizializza matrice di espressione
  expr_patterns <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Assegna geni alle strutture (distribuzione uniforme)
  gene_structures <- sample(1:n_structures, n_aniso_genes, replace = TRUE)
  
  # Per ogni gene anisotropico
  for (g_idx in 1:n_aniso_genes) {
    g <- aniso_genes[g_idx]
    struct_id <- gene_structures[g_idx]
    
    # Cellule in questa struttura
    struct_cells <- which(structure_mask == struct_id)
    
    if (length(struct_cells) == 0) {
      next  # Salta se non ci sono cellule in questa struttura
    }
    
    # Estrai la matrice di distanza per questa struttura
    dist_mat <- distance_matrices[[struct_id]]
    
    # Genera pattern in base al tipo specificato
    if (pattern_type == "gradient") {
      # Scegli una cellula come punto di inizio del gradiente
      start_cell <- sample(struct_cells, 1)
      
      # Per ogni cellula nella struttura, calcola la distanza dalla cellula iniziale
      for (i in struct_cells) {
        # Distanza lungo la struttura
        distance <- dist_mat[start_cell, i]
        
        if (is.infinite(distance)) {
          # Fallback a distanza euclidea se non definita nella struttura
          distance <- sqrt((cell_df$x[i] - cell_df$x[start_cell])^2 + 
                          (cell_df$y[i] - cell_df$y[start_cell])^2)
        }
        
        # Normalizza la distanza
        max_dist <- max(dist_mat[start_cell, struct_cells])
        if (is.infinite(max_dist)) max_dist <- sqrt((max(cell_df$x) - min(cell_df$x))^2 + 
                                                  (max(cell_df$y) - min(cell_df$y))^2)
        
        norm_dist <- distance / max_dist
        
        # Trasforma in pattern di espressione (decrescente con la distanza)
        expr_patterns[i, g] <- 1 - norm_dist
      }
      
      # Aggiungi rumore al pattern
      if (smoothness > 0) {
        noise <- rnorm(length(struct_cells), 0, smoothness)
        expr_patterns[struct_cells, g] <- expr_patterns[struct_cells, g] + noise
      }
      
    } else if (pattern_type == "oscillating") {
      # Scegli una cellula come punto di inizio
      start_cell <- sample(struct_cells, 1)
      
      # Per ogni cellula nella struttura, calcola pattern oscillante basato sulla distanza
      for (i in struct_cells) {
        # Distanza lungo la struttura
        distance <- dist_mat[start_cell, i]
        
        if (is.infinite(distance)) {
          # Fallback a distanza euclidea
          distance <- sqrt((cell_df$x[i] - cell_df$x[start_cell])^2 + 
                          (cell_df$y[i] - cell_df$y[start_cell])^2)
        }
        
        # Pattern oscillante
        freq <- oscillation_freq * 2 * pi
        expr_patterns[i, g] <- sin(distance * freq)
      }
      
      # Aggiungi rumore al pattern
      if (smoothness > 0) {
        noise <- rnorm(length(struct_cells), 0, smoothness)
        expr_patterns[struct_cells, g] <- expr_patterns[struct_cells, g] + noise
      }
      
    } else if (pattern_type == "hotspot") {
      # Genera hotspot casuali lungo la struttura
      n_spots <- sample(1:n_hotspots, 1)
      hotspot_cells <- sample(struct_cells, n_spots)
      
      # Per ogni cellula nella struttura, calcola l'influenza dai hotspot
      for (i in struct_cells) {
        hotspot_influence <- 0
        
        for (h in hotspot_cells) {
          # Distanza dal hotspot
          distance <- dist_mat[h, i]
          
          if (is.infinite(distance)) {
            # Fallback a distanza euclidea
            distance <- sqrt((cell_df$x[i] - cell_df$x[h])^2 + 
                            (cell_df$y[i] - cell_df$y[h])^2)
          }
          
          # Calcola influenza del hotspot (decadimento esponenziale)
          radius <- hotspot_radius
          influence <- exp(-distance / radius)
          hotspot_influence <- max(hotspot_influence, influence)
        }
        
        expr_patterns[i, g] <- hotspot_influence
      }
      
      # Aggiungi rumore al pattern
      if (smoothness > 0) {
        noise <- rnorm(length(struct_cells), 0, smoothness)
        expr_patterns[struct_cells, g] <- expr_patterns[struct_cells, g] + noise
      }
    } else {
      # Pattern casuale come fallback
      expr_patterns[struct_cells, g] <- runif(length(struct_cells))
    }
    
    # Normalizza il pattern per questo gene
    if (length(struct_cells) > 0) {
      # Verifica se ci sono NA o NaN nei valori
      pattern_values <- expr_patterns[struct_cells, g]
      has_na <- any(is.na(pattern_values))
      
      if (has_na) {
        # Ripara valori NA/NaN
        expr_patterns[struct_cells, g][is.na(pattern_values)] <- 0.5
        pattern_values <- expr_patterns[struct_cells, g]
      }
      
      # Procedi con la normalizzazione
      min_val <- min(pattern_values)
      max_val <- max(pattern_values)
      
      if (!is.na(min_val) && !is.na(max_val) && max_val > min_val) {
        expr_patterns[struct_cells, g] <- (expr_patterns[struct_cells, g] - min_val) / (max_val - min_val)
      } else {
        # Fallback in caso di problemi con min/max
        expr_patterns[struct_cells, g] <- 0.5
      }
    }
  }
  
  # Per i geni non anisotropici, genera pattern casuali
  non_aniso_genes <- setdiff(1:n_genes, aniso_genes)
  expr_patterns[, non_aniso_genes] <- matrix(runif(n_cells * length(non_aniso_genes)), 
                                           nrow = n_cells, ncol = length(non_aniso_genes))
  
  return(expr_patterns)
}

#' Applica effetti anisotropici all'espressione genica (versione semplificata per testing)
#'
#' @param expr_matrix Matrice di espressione originale
#' @param aniso_matrix Matrice di pattern anisotropici
#' @param structure_mask Vettore di appartenenza delle cellule alle strutture
#' @param effect_params Parametri per gli effetti anisotropici
#' @return Matrice di espressione modificata con effetti anisotropici
#' @export
apply_anisotropic_effects <- function(
  expr_matrix,
  aniso_matrix,
  structure_mask,
  effect_params = NULL
) {
  # Estrai parametri con gestione di NULL
  if (is.null(effect_params)) {
    effect_type <- "multiplicative"
    effect_strength <- 0.8
    background_fraction <- 0.2
  } else {
    effect_type <- effect_params$anisotropic_effect_type
    effect_strength <- effect_params$anisotropic_effect_strength
    background_fraction <- effect_params$background_effect_fraction
    
    # Valori predefiniti se non specificati
    if (is.null(effect_type)) effect_type <- "multiplicative"
    if (is.null(effect_strength)) effect_strength <- 0.8
    if (is.null(background_fraction)) background_fraction <- 0.2
  }
  
  # Estrai dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  
  # Crea una matrice di effetti finali
  effect_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Normalizza la matrice anisotropica per ottenere effetti tra -1 e 1
  aniso_scaled <- aniso_matrix * 2 - 1
  
  # Applica effetto completo alle cellule nelle strutture e effetto ridotto alle altre
  for (i in 1:n_cells) {
    if (structure_mask[i] > 0) {
      # Cellula in una struttura
      effect_matrix[i, ] <- aniso_scaled[i, ] * effect_strength
    } else {
      # Cellula fuori dalle strutture - effetto ridotto
      effect_matrix[i, ] <- aniso_scaled[i, ] * effect_strength * background_fraction
    }
  }
  
  # Applica gli effetti in base al tipo
  if (effect_type == "multiplicative") {
    # Trasforma effetti in fattori moltiplicativi (range: 0.5 - 2.0 per effect_strength = 1.0)
    factors <- exp(effect_matrix)
    modified_expr <- expr_matrix * factors
    
  } else if (effect_type == "additive") {
    # Effetti additivi (aggiunge o sottrae)
    # Scala gli effetti in base al livello di espressione per realismo
    modified_expr <- expr_matrix + effect_matrix * expr_matrix * 0.5
    
  } else {
    # Approccio misto come fallback
    additive_component <- effect_matrix * expr_matrix * 0.3
    factors <- exp(effect_matrix * 0.7)
    modified_expr <- expr_matrix * factors + additive_component
  }
  
  # Assicura che l'espressione rimanga non-negativa
  modified_expr[modified_expr < 0] <- 0
  
  return(modified_expr)
}

#' Genera modello completo di pattern anisotropici (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expr_matrix Matrice di espressione originale
#' @param anisotropic_params Parametri per i pattern anisotropici
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrice modificata e informazioni sulle strutture
#' @export
generate_anisotropic_patterns <- function(
  cell_df,
  expr_matrix,
  anisotropic_params = NULL,
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Gestione parametri
  if (is.null(anisotropic_params)) {
    # Valori predefiniti
    use_anisotropic <- TRUE
    n_structures <- 2
    structure_type <- "linear"
    anisotropic_pattern <- "gradient"
    anisotropic_gene_fraction <- 0.6
    anisotropic_effect_strength <- 0.8
  } else {
    # Estrai dai parametri forniti
    use_anisotropic <- anisotropic_params$use_anisotropic_patterns
    if (is.null(use_anisotropic)) use_anisotropic <- TRUE
    
    n_structures <- anisotropic_params$n_structures
    if (is.null(n_structures)) n_structures <- 2
    
    structure_type <- anisotropic_params$structure_type
    if (is.null(structure_type)) structure_type <- "linear"
    
    anisotropic_pattern <- anisotropic_params$anisotropic_pattern
    if (is.null(anisotropic_pattern)) anisotropic_pattern <- "gradient"
    
    anisotropic_gene_fraction <- anisotropic_params$anisotropic_gene_fraction
    if (is.null(anisotropic_gene_fraction)) anisotropic_gene_fraction <- 0.6
    
    anisotropic_effect_strength <- anisotropic_params$anisotropic_effect_strength
    if (is.null(anisotropic_effect_strength)) anisotropic_effect_strength <- 0.8
  }
  
  # Se disabilitati, restituisci valori originali
  if (!use_anisotropic) {
    return(list(
      expr_matrix = expr_matrix,
      structure_mask = rep(0, nrow(cell_df)),
      distance_matrices = list()
    ))
  }
  
  # Estrai dimensioni
  n_cells <- nrow(cell_df)
  n_genes <- ncol(expr_matrix)
  
  # 1. Genera strutture backbone
  structures <- generate_backbone_structures(
    cell_df = cell_df,
    structure_params = list(
      n_structures = n_structures,
      structure_type = structure_type,
      structure_width = 2,
      structure_length_factor = 0.8
    ),
    random_seed = random_seed
  )
  
  # Estrai informazioni sulle strutture
  distance_matrices <- structures$distance_matrices
  structure_mask <- structures$structure_mask
  
  # 2. Genera pattern di espressione anisotropica
  aniso_patterns <- generate_anisotropic_expression(
    cell_df = cell_df,
    distance_matrices = distance_matrices,
    structure_mask = structure_mask,
    n_genes = n_genes,
    expression_params = list(
      anisotropic_pattern = anisotropic_pattern,
      anisotropic_gene_fraction = anisotropic_gene_fraction,
      pattern_smoothness = 0.2
    ),
    random_seed = random_seed
  )
  
  # 3. Applica gli effetti all'espressione
  modified_expr <- apply_anisotropic_effects(
    expr_matrix = expr_matrix,
    aniso_matrix = aniso_patterns,
    structure_mask = structure_mask,
    effect_params = list(
      anisotropic_effect_type = "multiplicative",
      anisotropic_effect_strength = anisotropic_effect_strength,
      background_effect_fraction = 0.2
    )
  )
  
  # Restituisci risultati
  return(list(
    expr_matrix = modified_expr,
    structure_mask = structure_mask,
    distance_matrices = distance_matrices,
    aniso_patterns = aniso_patterns
  ))
}

# Funzioni di utilità

#' Calcola la distanza da un punto a un segmento di linea
#'
#' @param px,py Coordinate del punto
#' @param x1,y1,x2,y2 Coordinate degli estremi del segmento
#' @return Distanza minima dal punto al segmento
distanceToSegment <- function(px, py, x1, y1, x2, y2) {
  # Vettore segmento
  A <- c(x2 - x1, y2 - y1)
  # Vettore da punto 1 al punto di interesse
  B <- c(px - x1, py - y1)
  
  # Lunghezza del segmento al quadrato
  len_sq <- sum(A^2)
  
  if (len_sq == 0) {
    # Il segmento è un punto
    return(sqrt(sum(B^2)))
  }
  
  # Proiezione di B su A (parametro t lungo il segmento)
  t <- max(0, min(1, sum(A * B) / len_sq))
  
  # Punto più vicino sul segmento
  projection <- c(x1, y1) + t * A
  
  # Distanza dal punto alla proiezione
  return(sqrt(sum((c(px, py) - projection)^2)))
}

#' Calcola il parametro t della proiezione di un punto su una linea
#'
#' @param px,py Coordinate del punto
#' @param x1,y1,x2,y2 Coordinate degli estremi del segmento
#' @return Parametro t della proiezione (0-1 se sul segmento)
projectPointOnLine <- function(px, py, x1, y1, x2, y2) {
  # Vettore segmento
  A <- c(x2 - x1, y2 - y1)
  # Vettore da punto 1 al punto di interesse
  B <- c(px - x1, py - y1)
  
  # Lunghezza del segmento al quadrato
  len_sq <- sum(A^2)
  
  if (len_sq == 0) {
    # Il segmento è un punto
    return(0)
  }
  
  # Proiezione di B su A (parametro t lungo il segmento)
  t <- sum(A * B) / len_sq
  
  # Limita t all'intervallo [0,1] per il segmento
  return(max(0, min(1, t)))
}

#' Costruisce un Minimum Spanning Tree per un insieme di punti
#'
#' @param x,y Vettori di coordinate dei punti
#' @return Matrice di indici che definiscono gli archi del MST
buildMST <- function(x, y) {
  n <- length(x)
  if (n <= 1) return(matrix(ncol = 2, nrow = 0))
  
  # Calcola distanze tra tutti i punti
  dist_mat <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d
    }
  }
  
  # Algoritmo di Prim per MST
  in_tree <- rep(FALSE, n)
  in_tree[1] <- TRUE
  edges <- matrix(ncol = 2, nrow = 0)
  
  for (t in 1:(n-1)) {
    min_dist <- Inf
    min_i <- 0
    min_j <- 0
    
    for (i in which(in_tree)) {
      for (j in which(!in_tree)) {
        if (dist_mat[i, j] < min_dist) {
          min_dist <- dist_mat[i, j]
          min_i <- i
          min_j <- j
        }
      }
    }
    
    edges <- rbind(edges, c(min_i, min_j))
    in_tree[min_j] <- TRUE
  }
  
  return(edges)
}