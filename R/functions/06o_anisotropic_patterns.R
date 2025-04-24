#' Pattern anisotropici dipendenti da strutture tissutali
#'
#' Giustificazione: Studi recenti (Burgess et al., Nature Methods 2022; Rood et al., 
#' Nature Biotechnology 2022) mostrano che molti pattern di espressione seguono 
#' strutture biologiche con orientamento specifico (vasi, fibre nervose).
#'
#' Questo modulo implementa:
#' - Simulazione di strutture vascolari o nervose con pattern di espressione associati
#' - Direzionalità variabile dell'anisotropia basata su un "backbone" strutturale
#' - Gradiente di espressione ortogonale alle strutture principali

#' Genera strutture lineari anisotropiche
#'
#' Crea backbone strutturali (es. vasi, fibre nervose) con direzionalità
#' specifica, sui quali possono essere ancorati pattern di espressione.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param structure_params Parametri per le strutture anisotropiche
#' @param random_seed Seed per riproducibilità
#' @return Lista con strutture lineari e campi direzionali
#' @importFrom stats runif rbinom rnorm
#' @importFrom sp coordinates SpatialPointsDataFrame SpatialLines SpatialLinesDataFrame Line Lines
#' @importFrom gstat gstat predict vgm
#' @export
generate_backbone_structures <- function(
  cell_df,
  structure_params = list(
    n_structures = 5,              # Numero di strutture lineari (es. vasi)
    structure_type = "vessel",     # Tipo di struttura: "vessel", "nerve", "boundary", o "mixed"
    curvature = 0.3,               # Livello di curvatura (0 = linee rette, 1 = molto curve)
    width_range = c(5, 15),        # Range delle larghezze delle strutture
    bifurcation_prob = 0.3,        # Probabilità di biforcazione
    edge_avoidance = 0.8,          # Tendenza a evitare i margini del tessuto
    attraction_to_clusters = 0.4,  # Tendenza a seguire/attraversare cluster specifici
    padding = 20                   # Padding attorno all'area di campionamento
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  n_structures <- structure_params$n_structures
  structure_type <- structure_params$structure_type
  curvature <- structure_params$curvature
  width_range <- structure_params$width_range
  bifurcation_prob <- structure_params$bifurcation_prob
  edge_avoidance <- structure_params$edge_avoidance
  attraction_to_clusters <- structure_params$attraction_to_clusters
  padding <- structure_params$padding
  
  # Valori predefiniti se non specificati
  if (is.null(n_structures)) n_structures <- 5
  if (is.null(structure_type)) structure_type <- "vessel"
  if (is.null(curvature)) curvature <- 0.3
  if (is.null(width_range)) width_range <- c(5, 15)
  if (is.null(bifurcation_prob)) bifurcation_prob <- 0.3
  if (is.null(edge_avoidance)) edge_avoidance <- 0.8
  if (is.null(attraction_to_clusters)) attraction_to_clusters <- 0.4
  if (is.null(padding)) padding <- 20
  
  # Limiti dell'area di campionamento
  x_range <- range(cell_df$x)
  y_range <- range(cell_df$y)
  width <- diff(x_range)
  height <- diff(y_range)
  
  # Area estesa con padding
  x_min <- x_range[1] - padding
  x_max <- x_range[2] + padding
  y_min <- y_range[1] - padding
  y_max <- y_range[2] + padding
  
  # Determina parametri specifici per tipo di struttura
  if (structure_type == "vessel") {
    # Vasi sanguigni tendono ad essere più curvi e con biforcazioni
    curvature <- max(curvature, 0.4)
    segment_length <- width / 15
    step_range <- c(width / 30, width / 15)
    width_factor <- 1.0  # Larghezza relativa
  } else if (structure_type == "nerve") {
    # Fibre nervose tendono ad essere più lineari con meno biforcazioni
    curvature <- min(curvature, 0.2)
    bifurcation_prob <- min(bifurcation_prob, 0.2)
    segment_length <- width / 12
    step_range <- c(width / 25, width / 12)
    width_factor <- 0.8  # Leggermente più sottili
  } else if (structure_type == "boundary") {
    # Confini di tessuto tendono a formare curve chiuse
    segment_length <- width / 10
    step_range <- c(width / 20, width / 10)
    width_factor <- 1.2  # Leggermente più larghi
  } else {
    # Valori di default per "mixed" o altro
    segment_length <- width / 15
    step_range <- c(width / 30, width / 15)
    width_factor <- 1.0
  }
  
  # Lista per salvare strutture
  structures <- list()
  
  # Lista di matrici per memorizzare campi di distanza
  distance_fields <- list()
  
  # Lista per memorizzare campi direzionali
  direction_fields <- list()
  
  # Contatore globale per le strutture (incluse le biforcazioni)
  structure_counter <- 0
  
  # Per ogni struttura richiesta
  for (s in 1:n_structures) {
    # Decidi se generare un vaso lineare o una struttura più complessa
    if (structure_type == "boundary" && s <= min(2, n_structures)) {
      # Per i confini tissutali, crea curve chiuse/circolari
      structure_data <- generate_boundary_structure(
        x_min, x_max, y_min, y_max, 
        curvature = curvature,
        n_points = ceiling(20 + runif(1, -5, 10))
      )
    } else {
      # Struttura standard (vaso o fibra)
      
      # Punto di partenza: inizia dai margini con possibilità di regolazione
      if (runif(1) < edge_avoidance) {
        # Inizia da una posizione interna
        start_x <- runif(1, x_min + padding, x_max - padding)
        start_y <- runif(1, y_min + padding, y_max - padding)
        
        # Direzione iniziale casuale
        angle <- runif(1, 0, 2*pi)
        direction <- c(cos(angle), sin(angle))
      } else {
        # Inizia dal margine
        edge <- sample(1:4, 1)  # 1=superiore, 2=destro, 3=inferiore, 4=sinistro
        
        if (edge == 1) {
          start_x <- runif(1, x_min, x_max)
          start_y <- y_min
          direction <- c(runif(1, -0.5, 0.5), 1)  # Principalmente verso l'alto
        } else if (edge == 2) {
          start_x <- x_max
          start_y <- runif(1, y_min, y_max)
          direction <- c(-1, runif(1, -0.5, 0.5))  # Principalmente verso sinistra
        } else if (edge == 3) {
          start_x <- runif(1, x_min, x_max)
          start_y <- y_max
          direction <- c(runif(1, -0.5, 0.5), -1)  # Principalmente verso il basso
        } else {
          start_x <- x_min
          start_y <- runif(1, y_min, y_max)
          direction <- c(1, runif(1, -0.5, 0.5))  # Principalmente verso destra
        }
        
        # Normalizza la direzione
        direction <- direction / sqrt(sum(direction^2))
      }
      
      # Genera percorso con random walk direzionato
      structure_data <- generate_directed_path(
        start_x, start_y, direction,
        x_min, x_max, y_min, y_max,
        curvature = curvature,
        segment_length = segment_length,
        n_steps = ceiling(min(width, height) / segment_length * runif(1, 0.5, 1.5)),
        step_range = step_range,
        bifurcation_prob = bifurcation_prob,
        padding = padding
      )
    }
    
    # Per ogni segmento della struttura
    for (i in 1:length(structure_data$paths)) {
      # Incrementa il contatore della struttura
      structure_counter <- structure_counter + 1
      
      # Estrai il percorso attuale
      current_path <- structure_data$paths[[i]]
      current_directions <- structure_data$directions[[i]]
      
      # Larghezza di questa struttura
      if (length(width_range) == 2) {
        structure_width <- runif(1, width_range[1], width_range[2]) * width_factor
      } else {
        structure_width <- width_range * width_factor
      }
      
      # Memorizza i dettagli della struttura
      structures[[structure_counter]] <- list(
        id = structure_counter,
        parent_id = structure_data$parent_ids[i],
        type = structure_type,
        path = current_path,
        directions = current_directions,
        width = structure_width
      )
      
      # Crea un campo di distanza per questa struttura
      distance_field <- calculate_distance_to_path(
        cell_df, current_path, current_directions, structure_width
      )
      
      distance_fields[[structure_counter]] <- distance_field
      
      # Crea un campo direzionale (vettori normalizzati che indicano la direzione della struttura)
      direction_field <- calculate_direction_field(
        cell_df, current_path, current_directions, distance_field, structure_width
      )
      
      direction_fields[[structure_counter]] <- direction_field
    }
  }
  
  # Crea campi di distanza aggregati per tutte le strutture
  combined_distance <- matrix(Inf, nrow = nrow(cell_df), ncol = 1)
  closest_structure <- rep(0, nrow(cell_df))
  
  # Per ogni struttura, aggiorna la distanza minima e l'ID della struttura più vicina
  for (s in 1:length(distance_fields)) {
    distances <- distance_fields[[s]]
    
    # Aggiorna la struttura più vicina per ogni cellula
    update_mask <- distances < combined_distance
    combined_distance[update_mask] <- distances[update_mask]
    closest_structure[update_mask] <- s
  }
  
  # Crea la matrice per il campo direzionale totale (2 colonne: x, y)
  combined_directions <- matrix(0, nrow = nrow(cell_df), ncol = 2)
  
  # Utilizza il campo direzionale della struttura più vicina
  for (i in 1:nrow(cell_df)) {
    if (closest_structure[i] > 0) {
      s <- closest_structure[i]
      combined_directions[i, ] <- direction_fields[[s]][i, ]
    }
  }
  
  # Restituisci tutti i dati
  return(list(
    structures = structures,
    distance_fields = distance_fields,
    direction_fields = direction_fields,
    combined_distance = combined_distance,
    closest_structure = closest_structure,
    combined_directions = combined_directions
  ))
}

#' Genera un percorso direzionato con possibili biforcazioni
#'
#' Funzione interna per generare un percorso principale con 
#' potenziali biforcazioni, usando un random walk direzionato.
#'
#' @param start_x Coordinata x iniziale
#' @param start_y Coordinata y iniziale
#' @param initial_direction Vettore direzione iniziale [x, y]
#' @param x_min Limite inferiore x
#' @param x_max Limite superiore x
#' @param y_min Limite inferiore y
#' @param y_max Limite superiore y
#' @param curvature Livello di curvatura del percorso (0-1)
#' @param segment_length Lunghezza di un segmento
#' @param n_steps Numero di passi per il percorso principale
#' @param step_range Range della lunghezza del passo
#' @param bifurcation_prob Probabilità di biforcazione in ogni passo
#' @param padding Padding dai margini
#' @param max_depth Profondità massima di biforcazione, per evitare ricorsione eccessiva
#' @return Lista con percorsi, direzioni e parent_ids
#' @keywords internal
generate_directed_path <- function(
  start_x, start_y, initial_direction,
  x_min, x_max, y_min, y_max,
  curvature = 0.3,
  segment_length = 10,
  n_steps = 20,
  step_range = c(5, 15),
  bifurcation_prob = 0.3,
  padding = 20,
  max_depth = 3
) {
  # Inizializza vettori per il percorso
  path <- matrix(0, nrow = n_steps + 1, ncol = 2)
  path[1, ] <- c(start_x, start_y)
  
  # Vettore per le direzioni
  directions <- matrix(0, nrow = n_steps + 1, ncol = 2)
  directions[1, ] <- initial_direction
  
  # Lista di percorsi e direzioni per gestire le biforcazioni
  paths <- list()
  all_directions <- list()
  parent_ids <- numeric(0)
  
  # Mantieni traccia delle biforcazioni
  bifurcation_points <- list()
  
  # Genera punti per il percorso principale
  valid_steps <- 1  # Conta solo i passi validi
  current_dir <- initial_direction
  
  for (i in 2:(n_steps + 1)) {
    # Variazione casuale della direzione (curvatura)
    angle_variation <- rnorm(1, 0, curvature * pi/4)  # Variazione più forte = più curvatura
    
    # Ruota il vettore direzione
    cos_val <- cos(angle_variation)
    sin_val <- sin(angle_variation)
    new_dir_x <- current_dir[1] * cos_val - current_dir[2] * sin_val
    new_dir_y <- current_dir[1] * sin_val + current_dir[2] * cos_val
    
    # Normalizza il nuovo vettore direzione
    new_dir <- c(new_dir_x, new_dir_y)
    new_dir <- new_dir / sqrt(sum(new_dir^2))
    
    # Lunghezza del passo
    step_length <- runif(1, step_range[1], step_range[2])
    
    # Calcola la nuova posizione
    new_x <- path[valid_steps, 1] + new_dir[1] * step_length
    new_y <- path[valid_steps, 2] + new_dir[2] * step_length
    
    # Controlla se la nuova posizione è all'interno dei limiti (con padding)
    if (new_x >= x_min && new_x <= x_max && new_y >= y_min && new_y <= y_max) {
      valid_steps <- valid_steps + 1
      path[valid_steps, ] <- c(new_x, new_y)
      directions[valid_steps, ] <- new_dir
      current_dir <- new_dir
      
      # Considera la possibilità di biforcazione
      if (runif(1) < bifurcation_prob && max_depth > 0) {
        bifurcation_points[[length(bifurcation_points) + 1]] <- list(
          step = valid_steps,
          position = c(new_x, new_y),
          direction = new_dir,
          depth = max_depth
        )
      }
    }
    
    # Esci se raggiungiamo il margine
    if (new_x < x_min + padding || new_x > x_max - padding || 
        new_y < y_min + padding || new_y > y_max - padding) {
      break
    }
  }
  
  # Trimming di punti non utilizzati
  path <- path[1:valid_steps, , drop = FALSE]
  directions <- directions[1:valid_steps, , drop = FALSE]
  
  # Aggiungi il percorso principale
  paths[[1]] <- path
  all_directions[[1]] <- directions
  parent_ids[1] <- 0  # Il percorso principale non ha parent
  
  # Processa le biforcazioni
  for (b in seq_along(bifurcation_points)) {
    bp <- bifurcation_points[[b]]
    
    # Crea una nuova direzione con una variazione maggiore
    angle_variation <- runif(1, pi/6, pi/3) * sample(c(-1, 1), 1)
    
    cos_val <- cos(angle_variation)
    sin_val <- sin(angle_variation)
    branch_dir_x <- bp$direction[1] * cos_val - bp$direction[2] * sin_val
    branch_dir_y <- bp$direction[1] * sin_val + bp$direction[2] * cos_val
    
    branch_dir <- c(branch_dir_x, branch_dir_y)
    branch_dir <- branch_dir / sqrt(sum(branch_dir^2))
    
    # Genera un ramo più corto
    branch_steps <- ceiling(n_steps * runif(1, 0.3, 0.7))
    
    # Crea un nuovo percorso per la biforcazione
    branch_result <- generate_directed_path(
      bp$position[1], bp$position[2], branch_dir,
      x_min, x_max, y_min, y_max,
      curvature = curvature * 1.2,  # Rami tendono ad essere più curvi
      segment_length = segment_length * 0.8,  # Rami più corti
      n_steps = branch_steps,
      step_range = step_range * 0.8,
      bifurcation_prob = bifurcation_prob * 0.5,  # Meno probabilità di ulteriori biforcazioni
      padding = padding,
      max_depth = bp$depth - 1
    )
    
    # Aggiungi rami alla lista di percorsi
    n_existing_paths <- length(paths)
    for (p in seq_along(branch_result$paths)) {
      path_idx <- n_existing_paths + p
      paths[[path_idx]] <- branch_result$paths[[p]]
      all_directions[[path_idx]] <- branch_result$directions[[p]]
      
      # Aggiusta l'ID del parent
      if (branch_result$parent_ids[p] == 0) {
        parent_ids[path_idx] <- 1  # Collega al percorso principale
      } else {
        parent_ids[path_idx] <- n_existing_paths + branch_result$parent_ids[p]
      }
    }
  }
  
  return(list(
    paths = paths,
    directions = all_directions,
    parent_ids = parent_ids
  ))
}

#' Genera struttura di confine tissutale
#'
#' Funzione interna per generare un confine di tessuto come curva chiusa.
#'
#' @param x_min Limite inferiore x
#' @param x_max Limite superiore x
#' @param y_min Limite inferiore y
#' @param y_max Limite superiore y
#' @param curvature Livello di curvatura
#' @param n_points Numero di punti nel percorso
#' @return Lista con percorsi, direzioni e parent_ids
#' @keywords internal
generate_boundary_structure <- function(
  x_min, x_max, y_min, y_max,
  curvature = 0.3,
  n_points = 20
) {
  width <- x_max - x_min
  height <- y_max - y_min
  center_x <- (x_min + x_max) / 2
  center_y <- (y_min + y_max) / 2
  
  # Scegli tra diverse forme di confine
  boundary_type <- sample(c("ellipse", "complex"), 1)
  
  if (boundary_type == "ellipse") {
    # Genera ellisse
    a <- width * (0.4 + runif(1, -0.1, 0.1))  # Semiasse maggiore
    b <- height * (0.4 + runif(1, -0.1, 0.1))  # Semiasse minore
    
    # Angoli per i punti dell'ellisse
    angles <- seq(0, 2*pi, length.out = n_points)
    
    # Crea l'ellisse
    path <- matrix(0, nrow = n_points, ncol = 2)
    directions <- matrix(0, nrow = n_points, ncol = 2)
    
    for (i in 1:n_points) {
      # Punto sull'ellisse
      path[i, 1] <- center_x + a * cos(angles[i])
      path[i, 2] <- center_y + b * sin(angles[i])
      
      # Direzione tangente all'ellisse
      directions[i, 1] <- -sin(angles[i])
      directions[i, 2] <- cos(angles[i])
      
      # Normalizza la direzione
      directions[i, ] <- directions[i, ] / sqrt(sum(directions[i, ]^2))
    }
  } else {
    # Genera forma complessa tramite una serie di punti di controllo
    n_control_points <- sample(5:8, 1)
    
    # Genera punti di controllo attorno a un cerchio
    control_angles <- seq(0, 2*pi, length.out = n_control_points)
    control_radius <- min(width, height) * 0.4
    control_points <- matrix(0, nrow = n_control_points, ncol = 2)
    
    # Distorci leggermente i raggi per creare forme irregolari
    radius_variations <- runif(n_control_points, 0.7, 1.3)
    
    for (i in 1:n_control_points) {
      control_points[i, 1] <- center_x + control_radius * radius_variations[i] * cos(control_angles[i])
      control_points[i, 2] <- center_y + control_radius * radius_variations[i] * sin(control_angles[i])
    }
    
    # Interpola tra punti di controllo
    path <- matrix(0, nrow = n_points, ncol = 2)
    directions <- matrix(0, nrow = n_points, ncol = 2)
    
    # Assicura che il percorso si chiuda tornando al primo punto
    control_points <- rbind(control_points, control_points[1, ])
    
    # Crea una interpolazione con spline per un contorno più liscio
    for (i in 1:n_points) {
      t <- (i - 1) / (n_points - 1) * n_control_points
      seg_idx <- floor(t) + 1
      t_frac <- t - floor(t)
      
      if (seg_idx < n_control_points + 1) {
        p1 <- control_points[seg_idx, ]
        p2 <- control_points[seg_idx + 1, ]
        
        # Interpolazione lineare con un po' di rumore
        noise_x <- rnorm(1, 0, min(width, height) * 0.02 * curvature)
        noise_y <- rnorm(1, 0, min(width, height) * 0.02 * curvature)
        
        path[i, 1] <- p1[1] * (1 - t_frac) + p2[1] * t_frac + noise_x
        path[i, 2] <- p1[2] * (1 - t_frac) + p2[2] * t_frac + noise_y
        
        # Calcola la direzione approssimativa (tangente)
        dx <- p2[1] - p1[1]
        dy <- p2[2] - p1[2]
        
        # Normalizza
        dir_norm <- sqrt(dx^2 + dy^2)
        if (dir_norm > 0) {
          directions[i, 1] <- dx / dir_norm
          directions[i, 2] <- dy / dir_norm
        } else {
          # Direzione predefinita se la norma è vicina a zero
          directions[i, ] <- c(1, 0)
        }
      }
    }
  }
  
  # Crea la struttura di output
  paths <- list(path)
  all_directions <- list(directions)
  parent_ids <- 0  # Nessun parent per il percorso principale
  
  return(list(
    paths = paths,
    directions = all_directions,
    parent_ids = parent_ids
  ))
}

#' Calcola la distanza di ogni cellula dal percorso
#'
#' Funzione interna per calcolare la distanza minima da ogni cellula al percorso.
#'
#' @param cell_df Dataframe delle cellule
#' @param path Matrice delle coordinate del percorso
#' @param directions Matrice delle direzioni lungo il percorso
#' @param structure_width Larghezza della struttura
#' @return Vettore di distanze per ogni cellula
#' @keywords internal
calculate_distance_to_path <- function(
  cell_df,
  path,
  directions,
  structure_width
) {
  n_cells <- nrow(cell_df)
  n_path_points <- nrow(path)
  
  # Vettore per le distanze
  min_distances <- rep(Inf, n_cells)
  
  # Per ogni cellula, calcola la distanza minima dal percorso
  for (i in 1:n_cells) {
    cell_pos <- c(cell_df$x[i], cell_df$y[i])
    
    # Calcola distanza per ogni segmento del percorso
    for (p in 1:(n_path_points - 1)) {
      segment_start <- path[p, ]
      segment_end <- path[p + 1, ]
      
      # Vettore del segmento
      segment_vec <- segment_end - segment_start
      segment_length <- sqrt(sum(segment_vec^2))
      
      if (segment_length > 0) {
        # Vettore normalizzato
        segment_dir <- segment_vec / segment_length
        
        # Vettore dalla cellula al punto iniziale del segmento
        cell_to_start <- cell_pos - segment_start
        
        # Proiezione sulla direzione del segmento
        projection <- sum(cell_to_start * segment_dir)
        
        # Punto più vicino sul segmento
        if (projection < 0) {
          # Prima del punto iniziale
          closest_point <- segment_start
        } else if (projection > segment_length) {
          # Dopo il punto finale
          closest_point <- segment_end
        } else {
          # Sul segmento
          closest_point <- segment_start + projection * segment_dir
        }
        
        # Calcola distanza dal punto più vicino
        distance <- sqrt(sum((cell_pos - closest_point)^2))
        
        # Aggiorna la distanza minima
        min_distances[i] <- min(min_distances[i], distance)
      }
    }
  }
  
  # Normalizza le distanze rispetto alla larghezza della struttura
  # (distanza/larghezza = 1 al bordo della struttura)
  normalized_distances <- min_distances / structure_width
  
  return(normalized_distances)
}

#' Calcola il campo direzionale per ogni cellula
#'
#' Funzione interna per calcolare la direzione della struttura più vicina ad ogni cellula.
#'
#' @param cell_df Dataframe delle cellule
#' @param path Matrice delle coordinate del percorso
#' @param directions Matrice delle direzioni lungo il percorso
#' @param distance_field Campo di distanza calcolato
#' @param structure_width Larghezza della struttura
#' @return Matrice di direzioni (x, y) per ogni cellula
#' @keywords internal
calculate_direction_field <- function(
  cell_df,
  path,
  directions,
  distance_field,
  structure_width
) {
  n_cells <- nrow(cell_df)
  n_path_points <- nrow(path)
  
  # Matrice per il campo direzionale
  direction_field <- matrix(0, nrow = n_cells, ncol = 2)
  
  # Per ogni cellula
  for (i in 1:n_cells) {
    cell_pos <- c(cell_df$x[i], cell_df$y[i])
    closest_point_idx <- 1
    min_distance <- Inf
    
    # Trova il punto più vicino sul percorso
    for (p in 1:n_path_points) {
      path_point <- path[p, ]
      distance <- sqrt(sum((cell_pos - path_point)^2))
      
      if (distance < min_distance) {
        min_distance <- distance
        closest_point_idx <- p
      }
    }
    
    # Assegna la direzione del punto più vicino
    direction_field[i, ] <- directions[closest_point_idx, ]
    
    # Per punti molto lontani, la direzione ha meno significato
    # Possiamo attenuare la direzione con la distanza
    if (distance_field[i] > 2) {  # Oltre il doppio della larghezza della struttura
      attenuation <- 1 / (distance_field[i] - 1)  # Attenua gradualmente
      direction_field[i, ] <- direction_field[i, ] * attenuation
    }
  }
  
  return(direction_field)
}

#' Genera campi di espressione anisotropici
#'
#' Genera campi di espressione genica che seguono strutture anisotropiche
#' come vasi sanguigni o fibre nervose.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param backbone_structures Risultato della funzione generate_backbone_structures
#' @param n_genes Numero totale di geni
#' @param expression_params Parametri per i pattern di espressione
#' @param random_seed Seed per riproducibilità
#' @return Matrice di effetti di espressione anisotropici
#' @importFrom stats runif rbinom
#' @export
generate_anisotropic_expression <- function(
  cell_df,
  backbone_structures,
  n_genes,
  expression_params = list(
    fraction_anisotropic_genes = 0.1,  # Frazione di geni con pattern anisotropici
    expression_decay_range = c(1, 3),  # Parametri di decadimento dell'espressione con la distanza
    gradient_strength = 0.7,           # Forza del gradiente perpendicolare alle strutture
    effect_size_range = c(0.5, 2.0),   # Range dell'effetto di espressione (min, max)
    corridor_effect = 0.4,             # Effetto "corridoio" per espressione lungo la struttura
    structure_type_specificity = 0.8   # Specificità per tipo di struttura (vessel vs nerve vs boundary)
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  fraction_anisotropic_genes <- expression_params$fraction_anisotropic_genes
  expression_decay_range <- expression_params$expression_decay_range
  gradient_strength <- expression_params$gradient_strength
  effect_size_range <- expression_params$effect_size_range
  corridor_effect <- expression_params$corridor_effect
  structure_type_specificity <- expression_params$structure_type_specificity
  
  # Valori predefiniti se non specificati
  if (is.null(fraction_anisotropic_genes)) fraction_anisotropic_genes <- 0.1
  if (is.null(expression_decay_range)) expression_decay_range <- c(1, 3)
  if (is.null(gradient_strength)) gradient_strength <- 0.7
  if (is.null(effect_size_range)) effect_size_range <- c(0.5, 2.0)
  if (is.null(corridor_effect)) corridor_effect <- 0.4
  if (is.null(structure_type_specificity)) structure_type_specificity <- 0.8
  
  # Estrai dati delle strutture
  structures <- backbone_structures$structures
  combined_distance <- backbone_structures$combined_distance
  closest_structure <- backbone_structures$closest_structure
  combined_directions <- backbone_structures$combined_directions
  
  # Se non ci sono strutture, restituisci una matrice vuota
  if (length(structures) == 0) {
    return(matrix(0, nrow = nrow(cell_df), ncol = n_genes))
  }
  
  # Determina i tipi di strutture presenti
  structure_types <- sapply(structures, function(s) s$type)
  unique_types <- unique(structure_types)
  
  # Numero di celle e geni anisotropici
  n_cells <- nrow(cell_df)
  n_anisotropic_genes <- round(n_genes * fraction_anisotropic_genes)
  
  # Inizializza matrice di effetti
  anisotropic_effects <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Seleziona geni anisotropici
  anisotropic_genes <- sample(1:n_genes, n_anisotropic_genes)
  
  # Assegna parametri e genera effetti per ciascun gene anisotropico
  for (g_idx in 1:length(anisotropic_genes)) {
    g <- anisotropic_genes[g_idx]
    
    # Tipo di espressione anisotropica per questo gene
    expression_type <- sample(c("parallel", "perpendicular", "mixed"), 1,
                            prob = c(0.4, 0.4, 0.2))
    
    # Seleziona un tipo di struttura se c'è specificità
    if (runif(1) < structure_type_specificity && length(unique_types) > 1) {
      preferred_type <- sample(unique_types, 1)
      
      # Trova le strutture del tipo selezionato
      type_structures <- which(structure_types == preferred_type)
      
      # Per le celle che hanno una struttura del tipo preferito come più vicina
      type_cells <- which(closest_structure %in% type_structures)
      
      # Crea maschera di specificità (1 per strutture preferite, valore basso per altre)
      type_specificity <- rep(0.2, n_cells)  # Valore base basso per altre strutture
      type_specificity[type_cells] <- 1      # Valore pieno per strutture preferite
    } else {
      # Nessuna specificità per tipo, tutte le strutture influenzano allo stesso modo
      type_specificity <- rep(1, n_cells)
    }
    
    # Parametri di espressione per questo gene
    decay_param <- runif(1, expression_decay_range[1], expression_decay_range[2])
    effect_size <- runif(1, effect_size_range[1], effect_size_range[2])
    
    # Calcola l'effetto per ogni cellula
    for (i in 1:n_cells) {
      # Distanza normalizzata dalla struttura più vicina
      dist <- combined_distance[i]
      
      # Direzioni della struttura
      dir_x <- combined_directions[i, 1]
      dir_y <- combined_directions[i, 2]
      
      # Calcola componente di base (distanza)
      if (dist < Inf) {
        # Decadimento esponenziale con la distanza
        distance_effect <- exp(-dist * decay_param)
        
        # Modello di espressione in base al tipo
        if (expression_type == "parallel") {
          # Espressione che segue la direzione della struttura
          # Con attenuazione per creare un "corridoio" di espressione
          if (corridor_effect > 0 && dist <= 1) {
            # Crea un effetto "corridoio" - espressione più forte lungo l'asse della struttura
            corridor_factor <- 1 - corridor_effect * (dist * 2 - 1)^2  # 1 al centro, decresce verso i bordi
            distance_effect <- distance_effect * corridor_factor
          }
          
        } else if (expression_type == "perpendicular") {
          # Espressione perpendicolare - crea un gradiente ortogonale alla struttura
          # Calcola il segno (quale lato della struttura)
          cell_pos <- c(cell_df$x[i], cell_df$y[i])
          
          # Struttura più vicina
          s_idx <- closest_structure[i]
          if (s_idx > 0) {
            # Calcola un vettore perpendicolare alla direzione della struttura
            perp_dir <- c(-dir_y, dir_x)  # Rotazione di 90 gradi
            
            # Punto sulla struttura
            closest_point_idx <- which.min(sapply(1:nrow(structures[[s_idx]]$path), function(p) {
              point <- structures[[s_idx]]$path[p, ]
              return(sum((cell_pos - point)^2))
            }))
            closest_point <- structures[[s_idx]]$path[closest_point_idx, ]
            
            # Vettore dalla struttura alla cellula
            struct_to_cell <- cell_pos - closest_point
            
            # Proiezione sul vettore perpendicolare per determinare il lato
            side_value <- sum(struct_to_cell * perp_dir)
            
            # Modifica l'effetto in base al lato (gradiente) - un lato positivo, un lato negativo
            side_factor <- tanh(side_value * gradient_strength)
            
            # Applica il fattore di lato
            distance_effect <- distance_effect * (1 + side_factor)
          }
          
        } else if (expression_type == "mixed") {
          # Combinazione di effetti paralleli e perpendicolari
          # Applica corridor_effect come nel caso parallelo
          if (corridor_effect > 0 && dist <= 1) {
            corridor_factor <- 1 - corridor_effect * (dist * 2 - 1)^2
            distance_effect <- distance_effect * corridor_factor
          }
          
          # Aggiungi anche effetto perpendicolare più debole
          cell_pos <- c(cell_df$x[i], cell_df$y[i])
          s_idx <- closest_structure[i]
          if (s_idx > 0) {
            perp_dir <- c(-dir_y, dir_x)
            closest_point_idx <- which.min(sapply(1:nrow(structures[[s_idx]]$path), function(p) {
              point <- structures[[s_idx]]$path[p, ]
              return(sum((cell_pos - point)^2))
            }))
            closest_point <- structures[[s_idx]]$path[closest_point_idx, ]
            struct_to_cell <- cell_pos - closest_point
            side_value <- sum(struct_to_cell * perp_dir)
            side_factor <- tanh(side_value * gradient_strength * 0.5)  # Effetto dimezzato
            distance_effect <- distance_effect * (1 + side_factor * 0.5)
          }
        }
        
        # Applica specificità per tipo e effetto finale
        anisotropic_effects[i, g] <- effect_size * distance_effect * type_specificity[i]
      }
    }
  }
  
  # Crea una tabella di metadati per i geni anisotropici
  anisotropic_metadata <- data.frame(
    gene_id = anisotropic_genes,
    is_anisotropic = TRUE,
    structure_type = rep(NA, length(anisotropic_genes)),
    expression_type = rep(NA, length(anisotropic_genes)),
    effect_size = rep(NA, length(anisotropic_genes)),
    stringsAsFactors = FALSE
  )
  
  # Restituisci gli effetti e i metadati
  return(list(
    anisotropic_effects = anisotropic_effects,
    anisotropic_metadata = anisotropic_metadata
  ))
}

#' Applica effetti di espressione anisotropici
#'
#' Integra pattern di espressione anisotropici nella matrice di espressione.
#'
#' @param expression_matrix Matrice di espressione originale
#' @param anisotropic_effects Effetti anisotropici dalla funzione generate_anisotropic_expression
#' @param integration_params Parametri per l'integrazione
#' @param random_seed Seed per riproducibilità
#' @return Matrice di espressione modificata con effetti anisotropici
#' @export
apply_anisotropic_effects <- function(
  expression_matrix,
  anisotropic_effects,
  integration_params = list(
    integration_mode = "multiplicative",  # "additive" o "multiplicative"
    integration_weight = 0.8              # Peso per l'integrazione
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  integration_mode <- integration_params$integration_mode
  integration_weight <- integration_params$integration_weight
  
  # Valori predefiniti se non specificati
  if (is.null(integration_mode)) integration_mode <- "multiplicative"
  if (is.null(integration_weight)) integration_weight <- 0.8
  
  # Estrai gli effetti anisotropici
  if (is.list(anisotropic_effects) && "anisotropic_effects" %in% names(anisotropic_effects)) {
    effects <- anisotropic_effects$anisotropic_effects
  } else {
    effects <- anisotropic_effects
  }
  
  # Verifica le dimensioni
  if (nrow(expression_matrix) != nrow(effects) || ncol(expression_matrix) != ncol(effects)) {
    stop("Le dimensioni della matrice di espressione e degli effetti anisotropici non corrispondono")
  }
  
  # Applica gli effetti in base al modo di integrazione
  if (integration_mode == "additive") {
    # Modo additivo
    modified_expression <- expression_matrix + integration_weight * effects
    
  } else if (integration_mode == "multiplicative") {
    # Modo moltiplicativo: converte gli effetti in fattori moltiplicativi
    factors <- exp(integration_weight * effects)
    modified_expression <- expression_matrix * factors
    
  } else {
    # Modo per difetto
    warning("Modalità di integrazione non riconosciuta, verrà usato il modo moltiplicativo")
    factors <- exp(integration_weight * effects)
    modified_expression <- expression_matrix * factors
  }
  
  # Assicura che non ci siano valori negativi
  modified_expression[modified_expression < 0] <- 0
  
  return(modified_expression)
}

#' Genera pattern di espressione anisotropici completi
#'
#' Funzione wrapper che esegue l'intero processo di generazione di pattern
#' anisotropici e la loro integrazione nell'espressione genica.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expression_matrix Matrice di espressione originale
#' @param anisotropy_params Parametri per l'anisotropia
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici di espressione e strutture anisotropiche
#' @export
generate_anisotropic_patterns <- function(
  cell_df,
  expression_matrix,
  anisotropy_params = list(
    # Parametri per le strutture
    n_structures = 5,
    structure_type = "mixed",  # "vessel", "nerve", "boundary", o "mixed"
    curvature = 0.3,
    width_range = c(5, 15),
    bifurcation_prob = 0.3,
    
    # Parametri per l'espressione
    fraction_anisotropic_genes = 0.1,
    expression_decay_range = c(1, 3),
    gradient_strength = 0.7,
    effect_size_range = c(0.5, 2.0),
    corridor_effect = 0.4,
    structure_type_specificity = 0.8,
    
    # Parametri per l'integrazione
    integration_mode = "multiplicative",
    integration_weight = 0.8
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Dimensioni
  n_genes <- ncol(expression_matrix)
  
  # Assegna tipo di struttura in caso di "mixed"
  if (anisotropy_params$structure_type == "mixed") {
    structure_types <- c("vessel", "nerve", "boundary")
    structure_params <- anisotropy_params
    
    n_structures <- anisotropy_params$n_structures
    n_types <- length(structure_types)
    
    # Dividi le strutture tra i vari tipi
    structures_per_type <- numeric(n_types)
    
    # Assicura almeno una struttura per tipo, poi distribuzione casuale
    structures_per_type[1:min(n_types, n_structures)] <- 1
    if (n_structures > n_types) {
      remaining <- n_structures - n_types
      distribution <- rmultinom(1, remaining, rep(1, n_types))[, 1]
      structures_per_type <- structures_per_type + distribution
    }
    
    # Genera strutture per ciascun tipo
    all_structures <- list()
    for (t in 1:n_types) {
      if (structures_per_type[t] > 0) {
        structure_params$structure_type <- structure_types[t]
        structure_params$n_structures <- structures_per_type[t]
        
        type_structures <- generate_backbone_structures(
          cell_df = cell_df,
          structure_params = structure_params,
          random_seed = random_seed + t
        )
        
        if (t == 1) {
          all_structures <- type_structures
        } else {
          # Unisci strutture
          n_existing <- length(all_structures$structures)
          for (s in 1:length(type_structures$structures)) {
            all_structures$structures[[n_existing + s]] <- type_structures$structures[[s]]
            type_structures$structures[[s]]$id <- n_existing + s
          }
          
          # Unisci campi di distanza
          for (s in 1:length(type_structures$distance_fields)) {
            all_structures$distance_fields[[n_existing + s]] <- type_structures$distance_fields[[s]]
          }
          
          # Unisci campi direzionali
          for (s in 1:length(type_structures$direction_fields)) {
            all_structures$direction_fields[[n_existing + s]] <- type_structures$direction_fields[[s]]
          }
          
          # Aggiorna distanza combinata e struttura più vicina
          update_mask <- type_structures$combined_distance < all_structures$combined_distance
          all_structures$combined_distance[update_mask] <- type_structures$combined_distance[update_mask]
          
          # Aggiorna indici delle strutture più vicine
          all_structures$closest_structure[update_mask] <- type_structures$closest_structure[update_mask] + n_existing
          
          # Aggiorna direzioni combinate
          direction_update <- update_mask
          all_structures$combined_directions[direction_update, ] <- type_structures$combined_directions[direction_update, ]
        }
      }
    }
    
    backbone_structures <- all_structures
  } else {
    # Crea strutture di un singolo tipo
    structure_params <- list(
      n_structures = anisotropy_params$n_structures,
      structure_type = anisotropy_params$structure_type,
      curvature = anisotropy_params$curvature,
      width_range = anisotropy_params$width_range,
      bifurcation_prob = anisotropy_params$bifurcation_prob
    )
    
    backbone_structures <- generate_backbone_structures(
      cell_df = cell_df,
      structure_params = structure_params,
      random_seed = random_seed
    )
  }
  
  # Genera pattern di espressione anisotropici
  expression_params <- list(
    fraction_anisotropic_genes = anisotropy_params$fraction_anisotropic_genes,
    expression_decay_range = anisotropy_params$expression_decay_range,
    gradient_strength = anisotropy_params$gradient_strength,
    effect_size_range = anisotropy_params$effect_size_range,
    corridor_effect = anisotropy_params$corridor_effect,
    structure_type_specificity = anisotropy_params$structure_type_specificity
  )
  
  anisotropic_effects <- generate_anisotropic_expression(
    cell_df = cell_df,
    backbone_structures = backbone_structures,
    n_genes = n_genes,
    expression_params = expression_params,
    random_seed = random_seed
  )
  
  # Applica gli effetti all'espressione
  integration_params <- list(
    integration_mode = anisotropy_params$integration_mode,
    integration_weight = anisotropy_params$integration_weight
  )
  
  modified_expression <- apply_anisotropic_effects(
    expression_matrix = expression_matrix,
    anisotropic_effects = anisotropic_effects,
    integration_params = integration_params,
    random_seed = random_seed
  )
  
  # Restituisci i risultati
  return(list(
    expression = modified_expression,
    original_expression = expression_matrix,
    anisotropic_effects = anisotropic_effects$anisotropic_effects,
    anisotropic_metadata = anisotropic_effects$anisotropic_metadata,
    backbone_structures = backbone_structures
  ))
}