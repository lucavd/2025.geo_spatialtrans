#' Modellazione semplificata di dinamiche temporali per testing
#'
#' Versioni semplificate delle funzioni di dinamiche temporali adattate per testing
#' 
#' @import stats

#' Genera campo di pseudotempo (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param temporal_params Parametri per il campo di pseudotempo
#' @param random_seed Seed per riproducibilità
#' @return Vettore di pseudotempo per ogni cellula
#' @export
generate_pseudotime_field <- function(
  cell_df,
  temporal_params = list(
    pseudotime_mode = "gradient",
    pseudotime_origin = c(0, 0),
    n_foci = 3,
    focal_radius = 0.2,
    bifurcation_point = c(0.5, 0.5),
    branch_angles = c(45, 135)
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  mode <- temporal_params$pseudotime_mode
  origin <- temporal_params$pseudotime_origin
  n_foci <- temporal_params$n_foci
  focal_radius <- temporal_params$focal_radius
  bifurcation_point <- temporal_params$bifurcation_point
  branch_angles <- temporal_params$branch_angles
  
  # Valori predefiniti se non specificati
  if (is.null(mode)) mode <- "gradient"
  if (is.null(origin)) origin <- c(0, 0)
  if (is.null(n_foci)) n_foci <- 3
  if (is.null(focal_radius)) focal_radius <- 0.2
  if (is.null(bifurcation_point)) bifurcation_point <- c(0.5, 0.5)
  if (is.null(branch_angles)) branch_angles <- c(45, 135)
  
  # Estrai dimensioni
  n_cells <- nrow(cell_df)
  
  # Normalizza coordinate per facilitare calcoli
  x_range <- range(cell_df$x)
  y_range <- range(cell_df$y)
  
  x_norm <- (cell_df$x - x_range[1]) / (x_range[2] - x_range[1])
  y_norm <- (cell_df$y - y_range[1]) / (y_range[2] - y_range[1])
  
  # Inizializza vettore di pseudotempo
  pseudotime <- rep(0, n_cells)
  
  # Genera pattern di pseudotempo in base al tipo
  if (mode == "gradient") {
    # Gradiente lineare dall'origine
    # Normalizza origine se non già in [0,1]
    if (any(origin < 0) || any(origin > 1)) {
      origin_x <- (origin[1] - x_range[1]) / (x_range[2] - x_range[1])
      origin_y <- (origin[2] - y_range[1]) / (y_range[2] - y_range[1])
    } else {
      origin_x <- origin[1]
      origin_y <- origin[2]
    }
    
    # Calcola distanza dall'origine
    for (i in 1:n_cells) {
      pseudotime[i] <- sqrt((x_norm[i] - origin_x)^2 + (y_norm[i] - origin_y)^2)
    }
    
    # Normalizza a [0,1]
    pseudotime <- pseudotime / max(pseudotime)
    
  } else if (mode == "focal") {
    # Pattern con foci puntuali
    # Genera foci casuali
    focal_x <- runif(n_foci)
    focal_y <- runif(n_foci)
    
    # Per ogni cellula, calcola distanza minima dai foci
    for (i in 1:n_cells) {
      min_dist <- Inf
      for (f in 1:n_foci) {
        dist <- sqrt((x_norm[i] - focal_x[f])^2 + (y_norm[i] - focal_y[f])^2)
        min_dist <- min(min_dist, dist)
      }
      
      # Inverti la distanza per ottenere pseudotempo alto vicino ai foci
      pseudotime[i] <- exp(-min_dist / focal_radius)
    }
    
    # Normalizza a [0,1]
    pseudotime <- (pseudotime - min(pseudotime)) / (max(pseudotime) - min(pseudotime))
    
  } else if (mode == "bifurcation") {
    # Pattern con biforcazione
    # Normalizza punto di biforcazione se non già in [0,1]
    if (any(bifurcation_point < 0) || any(bifurcation_point > 1)) {
      bif_x <- (bifurcation_point[1] - x_range[1]) / (x_range[2] - x_range[1])
      bif_y <- (bifurcation_point[2] - y_range[1]) / (y_range[2] - y_range[1])
    } else {
      bif_x <- bifurcation_point[1]
      bif_y <- bifurcation_point[2]
    }
    
    # Converti angoli in radianti
    angles_rad <- branch_angles * pi / 180
    
    # Calcola pseudotempo come minima distanza da uno dei rami
    for (i in 1:n_cells) {
      # Distanza al punto di biforcazione
      dist_to_bif <- sqrt((x_norm[i] - bif_x)^2 + (y_norm[i] - bif_y)^2)
      
      if (dist_to_bif < 0.1) {
        # Vicino alla biforcazione, pseudotempo proporzionale alla distanza
        pseudotime[i] <- dist_to_bif * 5  # Scala per dare valori ragionevoli
      } else {
        # Calcola distanza minima dai rami
        branch_dist <- rep(Inf, length(angles_rad))
        
        for (b in 1:length(angles_rad)) {
          # Vettore unitario lungo il ramo
          vx <- cos(angles_rad[b])
          vy <- sin(angles_rad[b])
          
          # Vettore dalla biforcazione alla cellula
          dx <- x_norm[i] - bif_x
          dy <- y_norm[i] - bif_y
          
          # Proiezione ortogonale
          proj <- dx * vx + dy * vy
          
          if (proj <= 0) {
            # Cellula "prima" della biforcazione
            branch_dist[b] <- dist_to_bif
          } else {
            # Distanza ortogonale al ramo
            orth_dist <- abs(dx * vy - dy * vx)
            
            # Distanza totale considerando sia proiezione che ortogonalità
            branch_dist[b] <- sqrt(orth_dist^2 + (dist_to_bif - proj)^2)
          }
        }
        
        # Pseudotempo basato sulla minima distanza da un ramo
        min_branch_dist <- min(branch_dist)
        pseudotime[i] <- 0.2 + min_branch_dist * 0.8  # Scala per dare valori ragionevoli
      }
    }
    
    # Normalizza a [0,1]
    pseudotime <- (pseudotime - min(pseudotime)) / (max(pseudotime) - min(pseudotime))
  } else {
    # Gradiente come fallback
    for (i in 1:n_cells) {
      pseudotime[i] <- sqrt((x_norm[i] - 0.5)^2 + (y_norm[i] - 0.5)^2)
    }
    
    # Normalizza a [0,1]
    pseudotime <- pseudotime / max(pseudotime)
  }
  
  return(pseudotime)
}

#' Genera traiettorie geniche lungo pseudotempo (versione semplificata per testing)
#'
#' @param pseudotime Vettore di pseudotempo per ogni cellula
#' @param n_genes Numero di geni da modellare
#' @param trajectory_params Parametri per le traiettorie di espressione
#' @param gene_modules Lista di moduli genici (opzionale)
#' @param random_seed Seed per riproducibilità
#' @return Matrice di traiettorie di espressione
#' @export
generate_gene_trajectories <- function(
  pseudotime,
  n_genes,
  trajectory_params = list(
    pattern_distribution = c(monotonic = 0.4, transient = 0.3, cyclic = 0.2, bifurcating = 0.1),
    trajectory_smoothness = 0.1,
    use_gene_modules = FALSE
  ),
  gene_modules = NULL,
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  pattern_dist <- trajectory_params$pattern_distribution
  smoothness <- trajectory_params$trajectory_smoothness
  use_modules <- trajectory_params$use_gene_modules
  
  # Valori predefiniti se non specificati
  if (is.null(pattern_dist)) {
    pattern_dist <- c(monotonic = 0.4, transient = 0.3, cyclic = 0.2, bifurcating = 0.1)
  }
  if (is.null(smoothness)) smoothness <- 0.1
  if (is.null(use_modules)) use_modules <- FALSE
  
  # Normalizza la distribuzione dei pattern
  pattern_dist <- pattern_dist / sum(pattern_dist)
  
  # Estrai dimensioni
  n_cells <- length(pseudotime)
  
  # Ordina pseudotempo per facilitare la generazione di traiettorie
  pt_order <- order(pseudotime)
  
  # Inizializza matrice di traiettorie
  trajectories <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Determina tipo di pattern per ogni gene
  pattern_types <- sample(names(pattern_dist), n_genes, replace = TRUE, 
                         prob = pattern_dist)
  
  # Se usiamo moduli, raggruppa geni simili
  gene_groups <- 1:n_genes  # Default: ogni gene è nel suo gruppo
  
  if (use_modules && !is.null(gene_modules)) {
    # Usa i moduli forniti
    gene_groups <- rep(0, n_genes)
    
    for (m in 1:length(gene_modules)) {
      module_genes <- gene_modules[[m]]
      gene_groups[module_genes] <- m
    }
  } else if (use_modules) {
    # Crea moduli casuali se non forniti
    n_modules <- max(1, round(n_genes / 10))
    gene_groups <- sample(1:n_modules, n_genes, replace = TRUE)
  }
  
  # Genera parametri per ogni gruppo
  unique_groups <- unique(gene_groups)
  group_params <- list()
  
  for (g in unique_groups) {
    # Parametri comuni per tutti i geni nel gruppo
    group_params[[as.character(g)]] <- list(
      amplitude = runif(1, 0.5, 2.0),
      phase = runif(1, 0, 2*pi),
      peak_time = runif(1),
      peak_width = runif(1, 0.1, 0.3),
      bifurcation_point = runif(1, 0.3, 0.7),
      bifurcation_strength = runif(1, 0.5, 2.0)
    )
  }
  
  # Genera traiettoria per ogni gene
  for (g in 1:n_genes) {
    pattern <- pattern_types[g]
    group <- gene_groups[g]
    params <- group_params[[as.character(group)]]
    
    # Applica variazione per distinguere geni nello stesso gruppo
    gene_params <- list(
      amplitude = params$amplitude * runif(1, 0.8, 1.2),
      phase = params$phase + runif(1, -0.5, 0.5),
      peak_time = params$peak_time + runif(1, -0.1, 0.1),
      peak_width = params$peak_width * runif(1, 0.8, 1.2),
      bifurcation_point = params$bifurcation_point + runif(1, -0.05, 0.05),
      bifurcation_strength = params$bifurcation_strength * runif(1, 0.9, 1.1)
    )
    
    # Genera traiettoria in base al tipo di pattern
    trajectory <- rep(0, n_cells)
    
    if (pattern == "monotonic") {
      # Pattern monotono (crescente o decrescente)
      is_increasing <- sample(c(TRUE, FALSE), 1)
      
      for (i in 1:n_cells) {
        pt <- pseudotime[i]
        if (is_increasing) {
          trajectory[i] <- gene_params$amplitude * (pt^1.5)  # Crescita non lineare
        } else {
          trajectory[i] <- gene_params$amplitude * (1 - pt^1.5)  # Decrescita non lineare
        }
      }
      
    } else if (pattern == "transient") {
      # Pattern transiente (picco o valle)
      is_peak <- sample(c(TRUE, FALSE), 1)
      
      for (i in 1:n_cells) {
        pt <- pseudotime[i]
        dist_to_peak <- abs(pt - gene_params$peak_time)
        
        if (is_peak) {
          # Picco gaussiano
          trajectory[i] <- gene_params$amplitude * exp(-(dist_to_peak^2) / (2 * gene_params$peak_width^2))
        } else {
          # Valle gaussiana
          trajectory[i] <- gene_params$amplitude * (1 - exp(-(dist_to_peak^2) / (2 * gene_params$peak_width^2)))
        }
      }
      
    } else if (pattern == "cyclic") {
      # Pattern ciclico (oscillante)
      frequency <- runif(1, 2, 5)  # Numero di cicli
      
      for (i in 1:n_cells) {
        pt <- pseudotime[i]
        trajectory[i] <- gene_params$amplitude * 0.5 * (1 + sin(2 * pi * frequency * pt + gene_params$phase))
      }
      
    } else if (pattern == "bifurcating") {
      # Pattern con biforcazione
      branch <- sample(c(-1, 1), n_cells, replace = TRUE)  # Assegna cellule casualmente ai rami
      
      for (i in 1:n_cells) {
        pt <- pseudotime[i]
        
        if (pt < gene_params$bifurcation_point) {
          # Prima della biforcazione
          trajectory[i] <- gene_params$amplitude * pt / gene_params$bifurcation_point
        } else {
          # Dopo la biforcazione
          branch_effect <- (pt - gene_params$bifurcation_point) / (1 - gene_params$bifurcation_point)
          branch_effect <- branch_effect * gene_params$bifurcation_strength * branch[i]
          
          # Base + effetto specifico del ramo
          trajectory[i] <- gene_params$amplitude + branch_effect
        }
      }
    } else {
      # Pattern lineare come fallback
      for (i in 1:n_cells) {
        trajectory[i] <- gene_params$amplitude * pseudotime[i]
      }
    }
    
    # Aggiungi rumore alla traiettoria
    if (smoothness > 0) {
      noise <- rnorm(n_cells, 0, smoothness * mean(abs(trajectory)))
      trajectory <- trajectory + noise
    }
    
    # Assicura valori non negativi
    trajectory[trajectory < 0] <- 0
    
    # Aggiungi alla matrice
    trajectories[, g] <- trajectory
  }
  
  return(trajectories)
}

#' Genera vettori di velocità RNA (versione semplificata per testing)
#'
#' @param expr_matrix Matrice di espressione
#' @param pseudotime Vettore di pseudotempo
#' @param velocity_params Parametri per i vettori di velocità
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici di velocità e RNA unspliced
#' @export
generate_rna_velocity <- function(
  expr_matrix,
  pseudotime,
  velocity_params = list(
    velocity_strength = 1.0,
    velocity_noise = 0.2
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  strength <- velocity_params$velocity_strength
  noise <- velocity_params$velocity_noise
  
  # Valori predefiniti se non specificati
  if (is.null(strength)) strength <- 1.0
  if (is.null(noise)) noise <- 0.2
  
  # Se la forza è 0, restituisci matrici di 0
  if (strength <= 0) {
    return(list(
      velocity = matrix(0, nrow = nrow(expr_matrix), ncol = ncol(expr_matrix)),
      unspliced = matrix(0, nrow = nrow(expr_matrix), ncol = ncol(expr_matrix))
    ))
  }
  
  # Estrai dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  
  # Ordina cellule per pseudotempo
  pt_order <- order(pseudotime)
  
  # Inizializza matrici
  velocity_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  unspliced_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Per ogni gene
  for (g in 1:n_genes) {
    # Estrai espressione del gene
    gene_expr <- expr_matrix[, g]
    
    # Calcola derivata dell'espressione rispetto al pseudotempo
    # Metodo: differenze finite
    
    # Ordina espressione per pseudotempo
    expr_ordered <- gene_expr[pt_order]
    
    # Calcola gradiente usando differenze finite
    gradient <- diff(c(expr_ordered[1], expr_ordered)) / diff(c(0, pseudotime[pt_order]))
    
    # Riordina nella sequenza originale
    gene_velocity <- rep(0, n_cells)
    gene_velocity[pt_order] <- gradient
    
    # Scala la velocità per la forza desiderata
    gene_velocity <- gene_velocity * strength
    
    # Aggiungi rumore alla velocità
    if (noise > 0) {
      velocity_noise <- rnorm(n_cells, 0, noise * mean(abs(gene_velocity)))
      gene_velocity <- gene_velocity + velocity_noise
    }
    
    # Modello semplice di RNA unspliced: espressione + velocità positiva
    gene_unspliced <- gene_expr + pmax(0, gene_velocity)
    
    # Aggiungi alle matrici
    velocity_matrix[, g] <- gene_velocity
    unspliced_matrix[, g] <- gene_unspliced
  }
  
  # Assicura che l'RNA unspliced sia non negativo
  unspliced_matrix[unspliced_matrix < 0] <- 0
  
  return(list(
    velocity = velocity_matrix,
    unspliced = unspliced_matrix
  ))
}

#' Genera modello completo di dinamiche temporali (versione semplificata per testing)
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expr_matrix Matrice di espressione originale
#' @param temporal_params Parametri per le dinamiche temporali
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrice modificata e informazioni temporali
#' @export
generate_temporal_dynamics <- function(
  cell_df,
  expr_matrix,
  temporal_params = list(
    use_temporal_dynamics = TRUE,
    pseudotime_mode = "gradient",
    temporal_gene_fraction = 0.8,
    pattern_distribution = c(monotonic = 0.7, transient = 0.3),
    trajectory_strength = 0.8,
    include_velocity = TRUE
  ),
  random_seed = 42
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Verifica se le dinamiche temporali sono abilitate
  use_temporal <- temporal_params$use_temporal_dynamics
  if (is.null(use_temporal)) use_temporal <- TRUE
  
  # Se disabilitate, restituisci valori originali
  if (!use_temporal) {
    return(list(
      expr_matrix = expr_matrix,
      pseudotime = rep(0, nrow(cell_df)),
      velocity = matrix(0, nrow = nrow(expr_matrix), ncol = ncol(expr_matrix)),
      unspliced = matrix(0, nrow = nrow(expr_matrix), ncol = ncol(expr_matrix))
    ))
  }
  
  # Estrai dimensioni
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  
  # 1. Genera campo di pseudotempo
  pseudotime <- generate_pseudotime_field(
    cell_df = cell_df,
    temporal_params = temporal_params,
    random_seed = random_seed
  )
  
  # Frazione di geni temporali
  temporal_fraction <- temporal_params$temporal_gene_fraction
  if (is.null(temporal_fraction)) temporal_fraction <- 0.8
  
  # Numero di geni che seguono traiettorie temporali
  n_temporal_genes <- round(n_genes * temporal_fraction)
  temporal_genes <- sample(1:n_genes, n_temporal_genes)
  
  # 2. Genera traiettorie geniche
  gene_trajectories <- generate_gene_trajectories(
    pseudotime = pseudotime,
    n_genes = n_temporal_genes,
    trajectory_params = list(
      pattern_distribution = temporal_params$pattern_distribution,
      trajectory_smoothness = 0.1
    ),
    random_seed = random_seed
  )
  
  # 3. Costruisci espressione modificata integrando le traiettorie
  modified_expr <- expr_matrix
  
  # Forza delle traiettorie
  traj_strength <- temporal_params$trajectory_strength
  if (is.null(traj_strength)) traj_strength <- 0.8
  
  # Combina espressione originale con traiettorie
  for (i in 1:n_temporal_genes) {
    g <- temporal_genes[i]
    
    # Espressione originale con un peso di 1 - traj_strength
    base_contrib <- expr_matrix[, g] * (1 - traj_strength)
    
    # Contributo della traiettoria con un peso di traj_strength
    traj_contrib <- gene_trajectories[, i] * traj_strength * mean(expr_matrix[, g])
    
    # Combina i contributi
    modified_expr[, g] <- base_contrib + traj_contrib
  }
  
  # 4. Genera velocità RNA se richiesto
  velocity <- matrix(0, nrow = n_cells, ncol = n_genes)
  unspliced <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Controlla se generare velocità
  include_velocity <- temporal_params$include_velocity
  if (is.null(include_velocity)) include_velocity <- TRUE
  
  if (include_velocity) {
    velocity_result <- generate_rna_velocity(
      expr_matrix = modified_expr,
      pseudotime = pseudotime,
      velocity_params = list(
        velocity_strength = 1.0,
        velocity_noise = 0.2
      ),
      random_seed = random_seed
    )
    
    velocity <- velocity_result$velocity
    unspliced <- velocity_result$unspliced
  }
  
  # Restituisci risultati
  return(list(
    expr_matrix = modified_expr,
    pseudotime = pseudotime,
    velocity = velocity,
    unspliced = unspliced
  ))
}