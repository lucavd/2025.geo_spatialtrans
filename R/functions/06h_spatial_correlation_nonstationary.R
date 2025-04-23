#' Genera pattern di correlazione spaziale non-stazionari
#'
#' Crea pattern di correlazione spaziale con parametri che variano 
#' in base alla posizione, creando regioni con struttura di correlazione diversa.
#'
#' @param cell_df Dataframe delle celle con coordinate
#' @param spatial_params Parametri spaziali non-stazionari
#' @param random_seed Seed per riproducibilità
#' @return Vettore con il rumore spaziale non-stazionario generato
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_nonstationary_correlation <- function(
  cell_df,
  spatial_params = list(
    # Parametri generali
    use_nonstationary = TRUE,
    nonstationary_type = "gradient",  # "gradient", "patch", o "adaptive"
    
    # Parametri per il metodo gradient
    range_min = 10,                  # Range minimo (piccola scala)
    range_max = 100,                 # Range massimo (grande scala)
    intensity_min = 0.5,             # Intensità minima
    intensity_max = 2.0,             # Intensità massima
    gradient_direction = c(1, 1),    # Direzione del gradiente
    gradient_strength = 1.0,         # Forza del gradiente
    
    # Parametri per il metodo patch
    n_patches = 5,                   # Numero di patch con parametri diversi
    patch_size_range = c(20, 50),    # Range delle dimensioni delle patch
    blend_patches = TRUE,            # Mescolare i confini delle patch
    blend_width = 10,                # Larghezza della zona di transizione
    
    # Parametri per il metodo adaptive
    adapt_to_feature = "density",    # "density" o "cluster"
    feature_impact = 0.8,            # Impatto della feature sul range/intensità
    feature_scaling = "linear"       # "linear" o "exponential"
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Controlla se usare correlazione non-stazionaria
  use_nonstationary <- ifelse(!is.null(spatial_params$use_nonstationary), 
                             spatial_params$use_nonstationary, TRUE)
  
  if (!use_nonstationary) {
    return(NULL)
  }
  
  # Converti cell_df in oggetto spatial
  sp_df <- cell_df
  coordinates(sp_df) <- ~ x + y
  coords <- coordinates(sp_df)
  N <- nrow(coords)
  
  # Determina il tipo di modello non-stazionario
  nonstationary_type <- ifelse(!is.null(spatial_params$nonstationary_type), 
                              spatial_params$nonstationary_type, "gradient")
  
  # Inizializza i risultati
  range_values <- rep(0, N)
  intensity_values <- rep(0, N)
  
  # ---- 1. Calcola parametri spazialmente variabili in base al metodo ----
  if (nonstationary_type == "gradient") {
    
    # Estrai parametri del gradiente
    range_min <- ifelse(!is.null(spatial_params$range_min), 
                       spatial_params$range_min, 10)
    
    range_max <- ifelse(!is.null(spatial_params$range_max), 
                       spatial_params$range_max, 100)
    
    intensity_min <- ifelse(!is.null(spatial_params$intensity_min), 
                           spatial_params$intensity_min, 0.5)
    
    intensity_max <- ifelse(!is.null(spatial_params$intensity_max), 
                           spatial_params$intensity_max, 2.0)
    
    # Direzione del gradiente (vettore normalizzato)
    grad_dir <- ifelse(!is.null(spatial_params$gradient_direction), 
                      list(spatial_params$gradient_direction), 
                      list(c(1, 1)))[[1]]
    
    grad_dir <- grad_dir / sqrt(sum(grad_dir^2))
    
    gradient_strength <- ifelse(!is.null(spatial_params$gradient_strength), 
                               spatial_params$gradient_strength, 1.0)
    
    # Calcola la proiezione di ogni punto sulla direzione del gradiente
    proj_values <- coords %*% grad_dir
    
    # Normalizza tra 0 e 1
    norm_proj <- (proj_values - min(proj_values)) / (max(proj_values) - min(proj_values))
    
    # Calcola i parametri usando il gradiente
    norm_proj_scaled <- norm_proj^gradient_strength  # Permette di rendere non-lineare
    range_values <- range_min + (range_max - range_min) * norm_proj_scaled
    intensity_values <- intensity_min + (intensity_max - intensity_min) * norm_proj_scaled
    
  } else if (nonstationary_type == "patch") {
    
    # Estrai parametri delle patch
    n_patches <- ifelse(!is.null(spatial_params$n_patches), 
                       spatial_params$n_patches, 5)
    
    patch_size_range <- ifelse(!is.null(spatial_params$patch_size_range), 
                              list(spatial_params$patch_size_range), 
                              list(c(20, 50)))[[1]]
    
    blend_patches <- ifelse(!is.null(spatial_params$blend_patches), 
                           spatial_params$blend_patches, TRUE)
    
    blend_width <- ifelse(!is.null(spatial_params$blend_width), 
                         spatial_params$blend_width, 10)
    
    # 1. Genera centri casuali delle patch
    # Determina il range delle coordinate
    x_range <- range(coords[, 1])
    y_range <- range(coords[, 2])
    
    # Genera i centri delle patch
    patch_centers <- matrix(
      c(
        runif(n_patches, x_range[1], x_range[2]),
        runif(n_patches, y_range[1], y_range[2])
      ),
      ncol = 2
    )
    
    # 2. Genera parametri casuali per ogni patch
    patch_params <- data.frame(
      range = runif(n_patches, 
                   spatial_params$range_min, 
                   spatial_params$range_max),
      intensity = runif(n_patches, 
                       spatial_params$intensity_min, 
                       spatial_params$intensity_max),
      size = runif(n_patches, 
                  patch_size_range[1], 
                  patch_size_range[2])
    )
    
    # 3. Assegna parametri a ogni punto in base alla patch più vicina
    if (blend_patches) {
      # Calcola la distanza di ogni punto da ogni patch
      dist_to_patches <- matrix(0, nrow = N, ncol = n_patches)
      for (p in 1:n_patches) {
        dist_to_patches[, p] <- sqrt(
          (coords[, 1] - patch_centers[p, 1])^2 + 
          (coords[, 2] - patch_centers[p, 2])^2
        )
      }
      
      # Calcola pesi inversamente proporzionali alla distanza e dimensione patch
      patch_weights <- matrix(0, nrow = N, ncol = n_patches)
      for (p in 1:n_patches) {
        # Attenuazione in base alla dimensione della patch
        patch_weights[, p] <- exp(-(dist_to_patches[, p] / patch_params$size[p])^2)
        
        # Rimuovi influenza per punti troppo lontani
        too_far <- dist_to_patches[, p] > patch_params$size[p] + blend_width
        patch_weights[too_far, p] <- 0
      }
      
      # Normalizza i pesi a 1 per ogni punto
      row_sums <- rowSums(patch_weights)
      row_sums[row_sums == 0] <- 1  # Evita divisioni per zero
      patch_weights <- patch_weights / row_sums
      
      # Calcola i parametri finali come media pesata dei parametri delle patch
      for (i in 1:N) {
        range_values[i] <- sum(patch_weights[i, ] * patch_params$range)
        intensity_values[i] <- sum(patch_weights[i, ] * patch_params$intensity)
      }
      
    } else {
      # Versione semplice: assegna ogni punto alla patch più vicina
      for (i in 1:N) {
        dist_to_patches <- apply(patch_centers, 1, function(center) {
          sqrt(sum((coords[i, ] - center)^2))
        })
        
        closest_patch <- which.min(dist_to_patches)
        range_values[i] <- patch_params$range[closest_patch]
        intensity_values[i] <- patch_params$intensity[closest_patch]
      }
    }
    
  } else if (nonstationary_type == "adaptive") {
    
    # Estrai parametri adattivi
    adapt_to_feature <- ifelse(!is.null(spatial_params$adapt_to_feature), 
                              spatial_params$adapt_to_feature, "density")
    
    feature_impact <- ifelse(!is.null(spatial_params$feature_impact), 
                            spatial_params$feature_impact, 0.8)
    
    feature_scaling <- ifelse(!is.null(spatial_params$feature_scaling), 
                             spatial_params$feature_scaling, "linear")
    
    range_min <- ifelse(!is.null(spatial_params$range_min), 
                       spatial_params$range_min, 10)
    
    range_max <- ifelse(!is.null(spatial_params$range_max), 
                       spatial_params$range_max, 100)
    
    intensity_min <- ifelse(!is.null(spatial_params$intensity_min), 
                           spatial_params$intensity_min, 0.5)
    
    intensity_max <- ifelse(!is.null(spatial_params$intensity_max), 
                           spatial_params$intensity_max, 2.0)
    
    # Calcola la feature a cui adattarsi
    feature_values <- numeric(N)
    if (adapt_to_feature == "density") {
      # Calcola la densità locale (numero di vicini entro una certa distanza)
      density_radius <- mean(c(range_min, range_max)) * 0.2
      
      for (i in 1:N) {
        dist_i <- sqrt(rowSums((coords - matrix(coords[i, ], nrow = N, ncol = 2, byrow = TRUE))^2))
        feature_values[i] <- sum(dist_i < density_radius) / N
      }
    } else if (adapt_to_feature == "cluster") {
      # Usa i cluster ID se disponibili, altrimenti crea cluster temporanei
      if ("intensity_cluster" %in% colnames(cell_df)) {
        cluster_ids <- as.integer(cell_df$intensity_cluster)
        
        # Calcola la distanza dal confine del cluster
        for (i in 1:N) {
          own_cluster <- cluster_ids[i]
          other_clusters <- which(cluster_ids != own_cluster)
          
          if (length(other_clusters) > 0) {
            # Trova la distanza minima a un punto di un altro cluster
            min_dist <- min(sqrt(rowSums((coords[other_clusters, ] - 
                                         matrix(coords[i, ], 
                                               nrow = length(other_clusters), 
                                               ncol = 2, 
                                               byrow = TRUE))^2)))
            
            # Normalizza la distanza
            feature_values[i] <- min_dist / max(min_dist)
          } else {
            feature_values[i] <- 1  # Se c'è un solo cluster
          }
        }
      } else {
        # Crea cluster temporanei usando k-means
        temp_clusters <- kmeans(coords, centers = 3)$cluster
        
        # Calcola la distanza dal confine del cluster
        for (i in 1:N) {
          own_cluster <- temp_clusters[i]
          other_clusters <- which(temp_clusters != own_cluster)
          
          if (length(other_clusters) > 0) {
            min_dist <- min(sqrt(rowSums((coords[other_clusters, ] - 
                                         matrix(coords[i, ], 
                                               nrow = length(other_clusters), 
                                               ncol = 2, 
                                               byrow = TRUE))^2)))
            
            # Usa un valore che diminuisce con la distanza
            feature_values[i] <- exp(-min_dist / 30)
          } else {
            feature_values[i] <- 0  # Improbabile con k-means
          }
        }
      }
    }
    
    # Normalizza la feature tra 0 e 1
    feature_values <- (feature_values - min(feature_values)) / 
      (max(feature_values) - min(feature_values) + 1e-10)
    
    # Applica scaling alla feature
    if (feature_scaling == "exponential") {
      feature_values <- feature_values^2
    }
    
    # Calcola i parametri adattandoli alla feature
    # Effetto sulla scala: alta densità -> range minore
    range_values <- range_max - (range_max - range_min) * 
      (feature_values * feature_impact + (1 - feature_impact) * 0.5)
    
    # Effetto sull'intensità: alta densità -> intensità maggiore
    intensity_values <- intensity_min + (intensity_max - intensity_min) * 
      (feature_values * feature_impact + (1 - feature_impact) * 0.5)
  }
  
  # ---- 2. Genera il noise non-stazionario usando processi gaussiani localizzati ----
  # Idea: dividiamo i punti in regioni simili e generiamo un GP per ciascuna
  
  # Discretizza i parametri in bin
  n_bins <- ifelse(N <= 1000, 5, 10)  # Meno bin per dataset piccoli
  range_bins <- cut(range_values, breaks = n_bins, labels = FALSE)
  intensity_bins <- cut(intensity_values, breaks = n_bins, labels = FALSE)
  
  # Crea un identificatore combinato per ogni regione
  region_ids <- (range_bins - 1) * n_bins + intensity_bins
  unique_regions <- sort(unique(region_ids))
  n_regions <- length(unique_regions)
  
  # Per ogni regione, genera un processo gaussiano con i parametri appropriati
  region_noise <- matrix(0, nrow = N, ncol = n_regions)
  for (r in 1:n_regions) {
    region <- unique_regions[r]
    region_points <- which(region_ids == region)
    
    if (length(region_points) > 1) {
      # Parametri medi per questa regione
      avg_range <- mean(range_values[region_points])
      avg_intensity <- mean(intensity_values[region_points])
      
      # Crea un GP per questa regione usando parametri locali
      # Usiamo un subset locale per efficienza computazionale nei dataset grandi
      if (length(region_points) > 200) {
        sample_points <- sample(region_points, 200)
      } else {
        sample_points <- region_points
      }
      
      # Oggetto spatial per la simulazione
      region_sp_df <- sp_df[sample_points, ]
      
      # Simula il campo gaussiano per questa regione
      region_gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                            beta = 0, model = vgm(psill = avg_intensity,
                                                 range = avg_range,
                                                 model = "Exp"),
                            nmax = 30)
      
      # Predice i valori per tutti i punti
      region_predict <- predict(region_gp_sim, newdata = sp_df, nsim = 1)$sim1
      
      # Memorizza il risultato
      region_noise[, r] <- region_predict
    }
  }
  
  # ---- 3. Combina il noise da tutte le regioni usando un approccio di miscelazione ----
  combined_noise <- numeric(N)
  
  # Distanze da ogni regione (usate per la miscelazione)
  for (i in 1:N) {
    # Trova la regione del punto corrente
    point_region <- which(unique_regions == region_ids[i])
    
    # Di default, usa principalmente il GP della propria regione
    region_weights <- rep(0, n_regions)
    region_weights[point_region] <- 0.7
    
    # Aggiunge contributi dalle regioni vicine
    # Trova le regioni con parametri simili
    similar_regions <- which(abs(unique_regions - region_ids[i]) <= n_bins + 1)
    similar_regions <- setdiff(similar_regions, point_region)
    
    # Distribuisce il restante 30% tra le regioni vicine
    if (length(similar_regions) > 0) {
      region_weights[similar_regions] <- 0.3 / length(similar_regions)
    } else {
      # Se non ci sono regioni simili, aumenta il peso della propria
      region_weights[point_region] <- 1.0
    }
    
    # Calcola il valore combinato
    combined_noise[i] <- sum(region_weights * region_noise[i, ])
  }
  
  # ---- 4. Finalizza il risultato ----
  # Normalizza il rumore finale
  final_noise <- scale(combined_noise)
  
  # Restituisci sia il rumore che i parametri spazialmente variabili
  return(list(
    noise = as.vector(final_noise),
    range_map = range_values,
    intensity_map = intensity_values,
    region_map = region_ids
  ))
}

#' Genera pattern di correlazione spaziale con anisotropia variabile
#'
#' Crea pattern di correlazione spaziale con direzioni preferenziali
#' che variano nello spazio, ad esempio per simulare flussi o strutture tissutali.
#'
#' @param cell_df Dataframe delle celle con coordinate
#' @param spatial_params Parametri spaziali per l'anisotropia
#' @param random_seed Seed per riproducibilità
#' @return Vettore con il rumore spaziale anisotropico generato
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_anisotropic_correlation <- function(
  cell_df,
  spatial_params = list(
    # Parametri generali
    use_anisotropy = TRUE,
    anisotropy_type = "flow",     # "flow", "radial", o "custom"
    
    # Parametri generali di correlazione
    range_main = 50,              # Range nella direzione principale
    range_secondary = 15,         # Range nella direzione secondaria (perpendicolare)
    intensity = 1.2,              # Intensità complessiva del pattern
    
    # Parametri per il tipo "flow"
    flow_field_type = "uniform",  # "uniform", "vortex", o "gradient"
    flow_direction = c(1, 0),     # Direzione del flusso uniforme
    flow_center = NULL,           # Centro del vortice (default: centro dell'area)
    flow_strength = 1.0,          # Intensità del campo di flusso
    
    # Parametri per il tipo "radial"
    radial_center = NULL,         # Centro del pattern radiale
    radial_type = "outward",      # "outward", "inward", o "circular"
    
    # Parametri per il tipo "custom"
    custom_direction_field = NULL # Matrice personalizzata di direzioni
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Controlla se usare correlazione anisotropica
  use_anisotropy <- ifelse(!is.null(spatial_params$use_anisotropy), 
                          spatial_params$use_anisotropy, TRUE)
  
  if (!use_anisotropy) {
    return(NULL)
  }
  
  # Converti cell_df in oggetto spatial
  sp_df <- cell_df
  coordinates(sp_df) <- ~ x + y
  coords <- coordinates(sp_df)
  N <- nrow(coords)
  
  # Estrai i parametri
  anisotropy_type <- ifelse(!is.null(spatial_params$anisotropy_type), 
                           spatial_params$anisotropy_type, "flow")
  
  range_main <- ifelse(!is.null(spatial_params$range_main), 
                      spatial_params$range_main, 50)
  
  range_secondary <- ifelse(!is.null(spatial_params$range_secondary), 
                           spatial_params$range_secondary, 15)
  
  intensity <- ifelse(!is.null(spatial_params$intensity), 
                     spatial_params$intensity, 1.2)
  
  # 1. Genera un campo di direzioni anisotropiche per ogni posizione
  direction_field <- matrix(0, nrow = N, ncol = 2)
  
  if (anisotropy_type == "flow") {
    # Estrai parametri del flusso
    flow_field_type <- ifelse(!is.null(spatial_params$flow_field_type), 
                             spatial_params$flow_field_type, "uniform")
    
    flow_strength <- ifelse(!is.null(spatial_params$flow_strength), 
                           spatial_params$flow_strength, 1.0)
    
    # Genera il campo di flusso
    if (flow_field_type == "uniform") {
      # Estrai direzione del flusso
      flow_dir <- ifelse(!is.null(spatial_params$flow_direction), 
                        list(spatial_params$flow_direction), 
                        list(c(1, 0)))[[1]]
      
      # Normalizza il vettore direzione
      flow_dir <- flow_dir / sqrt(sum(flow_dir^2))
      
      # Direzione uniforme per tutti i punti
      direction_field[, 1] <- flow_dir[1]
      direction_field[, 2] <- flow_dir[2]
      
    } else if (flow_field_type == "vortex") {
      # Determina il centro del vortice
      if (is.null(spatial_params$flow_center)) {
        # Default è il centro della griglia
        center_x <- mean(range(coords[, 1]))
        center_y <- mean(range(coords[, 2]))
      } else {
        center_x <- spatial_params$flow_center[1]
        center_y <- spatial_params$flow_center[2]
      }
      
      # Calcola il vettore dai punti al centro, poi ruota di 90 gradi
      for (i in 1:N) {
        # Vettore dal centro al punto
        dx <- coords[i, 1] - center_x
        dy <- coords[i, 2] - center_y
        
        # Distanza dal centro
        dist <- sqrt(dx^2 + dy^2)
        
        # Normalizza e ruota di 90 gradi (per ottenere vettori tangenziali)
        if (dist > 0) {
          # Vettore tangenziale normalizzato
          direction_field[i, 1] <- -dy / dist
          direction_field[i, 2] <- dx / dist
        } else {
          # Al centro, usa una direzione casuale
          random_angle <- runif(1, 0, 2*pi)
          direction_field[i, 1] <- cos(random_angle)
          direction_field[i, 2] <- sin(random_angle)
        }
      }
    } else if (flow_field_type == "gradient") {
      # Genera un campo gaussiano per simulare un gradiente
      gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                     beta = 0, model = vgm(psill = 1.0,
                                          range = max(range(coords[, 1]), 
                                                     range(coords[, 2])) / 3,
                                          model = "Sph"),
                     nmax = 40)
      
      grad_field <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
      
      # Calcola il gradiente approssimato discretizzando in una griglia
      grid_size <- ceiling(sqrt(N) / 10)
      if (grid_size < 5) grid_size <- 5
      
      x_range <- range(coords[, 1])
      y_range <- range(coords[, 2])
      
      x_grid <- seq(x_range[1], x_range[2], length.out = grid_size)
      y_grid <- seq(y_range[1], y_range[2], length.out = grid_size)
      
      grid_values <- matrix(NA, nrow = grid_size, ncol = grid_size)
      
      # Interpola valori sulla griglia
      for (i in 1:grid_size) {
        for (j in 1:grid_size) {
          # Trova i punti vicini al nodo della griglia
          x_node <- x_grid[i]
          y_node <- y_grid[j]
          
          dist_to_node <- sqrt((coords[, 1] - x_node)^2 + (coords[, 2] - y_node)^2)
          closest_points <- order(dist_to_node)[1:min(5, N)]
          
          # Media pesata in base alla distanza
          weights <- 1 / (dist_to_node[closest_points] + 1e-5)
          weights <- weights / sum(weights)
          
          grid_values[i, j] <- sum(grad_field[closest_points] * weights)
        }
      }
      
      # Calcola gradiente sulla griglia
      grad_x <- matrix(0, nrow = grid_size, ncol = grid_size)
      grad_y <- matrix(0, nrow = grid_size, ncol = grid_size)
      
      for (i in 2:(grid_size-1)) {
        for (j in 2:(grid_size-1)) {
          grad_x[i, j] <- (grid_values[i+1, j] - grid_values[i-1, j]) / 
                          (x_grid[i+1] - x_grid[i-1])
          
          grad_y[i, j] <- (grid_values[i, j+1] - grid_values[i, j-1]) / 
                          (y_grid[j+1] - y_grid[j-1])
        }
      }
      
      # Estrapola il gradiente a ciascun punto originale
      for (i in 1:N) {
        # Trova la cella della griglia contenente il punto
        x_idx <- findInterval(coords[i, 1], x_grid)
        y_idx <- findInterval(coords[i, 2], y_grid)
        
        # Limita agli indici validi
        x_idx <- min(max(x_idx, 2), grid_size-1)
        y_idx <- min(max(y_idx, 2), grid_size-1)
        
        # Prendi il gradiente dalla griglia
        gx <- grad_x[x_idx, y_idx]
        gy <- grad_y[x_idx, y_idx]
        
        # Normalizza
        grad_norm <- sqrt(gx^2 + gy^2)
        if (grad_norm > 0) {
          direction_field[i, 1] <- gx / grad_norm
          direction_field[i, 2] <- gy / grad_norm
        } else {
          # Se gradiente nullo, usa direzione casuale
          random_angle <- runif(1, 0, 2*pi)
          direction_field[i, 1] <- cos(random_angle)
          direction_field[i, 2] <- sin(random_angle)
        }
      }
    }
    
  } else if (anisotropy_type == "radial") {
    # Determina il centro del pattern radiale
    if (is.null(spatial_params$radial_center)) {
      # Default è il centro della griglia
      center_x <- mean(range(coords[, 1]))
      center_y <- mean(range(coords[, 2]))
    } else {
      center_x <- spatial_params$radial_center[1]
      center_y <- spatial_params$radial_center[2]
    }
    
    # Tipo di pattern radiale
    radial_type <- ifelse(!is.null(spatial_params$radial_type), 
                         spatial_params$radial_type, "outward")
    
    # Calcola vettori direzionali
    for (i in 1:N) {
      # Vettore dal centro al punto
      dx <- coords[i, 1] - center_x
      dy <- coords[i, 2] - center_y
      
      # Distanza dal centro
      dist <- sqrt(dx^2 + dy^2)
      
      if (dist > 0) {
        if (radial_type == "outward") {
          # Direzioni radiali uscenti dal centro
          direction_field[i, 1] <- dx / dist
          direction_field[i, 2] <- dy / dist
        } else if (radial_type == "inward") {
          # Direzioni radiali verso il centro
          direction_field[i, 1] <- -dx / dist
          direction_field[i, 2] <- -dy / dist
        } else if (radial_type == "circular") {
          # Direzioni circolari intorno al centro (tangenziali)
          direction_field[i, 1] <- -dy / dist
          direction_field[i, 2] <- dx / dist
        }
      } else {
        # Al centro, usa una direzione casuale
        random_angle <- runif(1, 0, 2*pi)
        direction_field[i, 1] <- cos(random_angle)
        direction_field[i, 2] <- sin(random_angle)
      }
    }
    
  } else if (anisotropy_type == "custom" && !is.null(spatial_params$custom_direction_field)) {
    # Usa il campo di direzioni personalizzato fornito dall'utente
    direction_field <- spatial_params$custom_direction_field
    
    # Assicurati che abbia le dimensioni corrette
    if (nrow(direction_field) != N || ncol(direction_field) != 2) {
      warning("Campo di direzioni personalizzato ha dimensioni errate. Usando direzioni casuali.")
      
      # Genera direzioni casuali
      random_angles <- runif(N, 0, 2*pi)
      direction_field <- cbind(cos(random_angles), sin(random_angles))
    }
  } else {
    # Default: direzioni casuali
    random_angles <- runif(N, 0, 2*pi)
    direction_field <- cbind(cos(random_angles), sin(random_angles))
  }
  
  # 2. Genera processo gaussiano anisotropico
  
  # Metodo di approssimazione: dividi lo spazio in regioni con direzione simile
  # e genera processi gaussiani locali con anisotropia allineata alla direzione
  
  # Raggruppa punti con direzioni simili (metodo k-means sulle direzioni)
  n_direction_groups <- min(ifelse(N <= 500, 5, 10), max(2, N / 5))
  direction_clusters <- kmeans(direction_field, centers = n_direction_groups)$cluster
  
  # Per ogni gruppo direzionale, genera un processo gaussiano anisotropico
  group_noise <- matrix(0, nrow = N, ncol = n_direction_groups)
  
  for (g in 1:n_direction_groups) {
    # Punti in questo gruppo
    group_points <- which(direction_clusters == g)
    
    if (length(group_points) > 0) {
      # Direzione media del gruppo
      mean_direction <- colMeans(direction_field[group_points, , drop = FALSE])
      mean_direction <- mean_direction / sqrt(sum(mean_direction^2))
      
      # Definisci matrice di anisotropia
      # La matrice di rotazione per allineare alla direzione principale
      angle <- atan2(mean_direction[2], mean_direction[1])
      
      # Per efficienza computazionale nei dataset grandi, usa un subset
      if (length(group_points) > 200) {
        subset_points <- sample(group_points, 200)
      } else {
        subset_points <- group_points
      }
      
      # Creiamo un nuovo dataframe per simulare con anisotropia
      local_coords <- coords[subset_points, ]
      
      # Ruota le coordinate per allineare con gli assi
      rotated_coords <- cbind(
        local_coords[, 1] * cos(-angle) - local_coords[, 2] * sin(-angle),
        local_coords[, 1] * sin(-angle) + local_coords[, 2] * cos(-angle)
      )
      
      # Scala le coordinate per creare l'anisotropia
      aniso_factor <- range_secondary / range_main
      scaled_coords <- cbind(
        rotated_coords[, 1],
        rotated_coords[, 2] * aniso_factor
      )
      
      # Crea un oggetto spatial con le coordinate trasformate
      group_sp_df <- sp_df[subset_points, ]
      coordinates(group_sp_df) <- scaled_coords
      
      # Simula il campo gaussiano sulle coordinate trasformate
      gp_sim <- gstat(formula = z ~ 1, locations = coordinates(group_sp_df), dummy = TRUE,
                     beta = 0, model = vgm(psill = intensity,
                                          range = range_main,
                                          model = "Exp"),
                     nmax = 30)
      
      # Crea un oggetto spatial per la predizione con coordinate trasformate
      predict_coords <- coords
      rotated_predict <- cbind(
        predict_coords[, 1] * cos(-angle) - predict_coords[, 2] * sin(-angle),
        predict_coords[, 1] * sin(-angle) + predict_coords[, 2] * cos(-angle)
      )
      
      scaled_predict <- cbind(
        rotated_predict[, 1],
        rotated_predict[, 2] * aniso_factor
      )
      
      predict_sp_df <- sp_df
      coordinates(predict_sp_df) <- scaled_predict
      
      # Predici i valori per tutti i punti
      group_field <- predict(gp_sim, newdata = predict_sp_df, nsim = 1)$sim1
      
      # Memorizza i risultati
      group_noise[, g] <- group_field
    }
  }
  
  # 3. Miscela i campi locali in base alla similarità di direzione
  
  # Calcola pesi per la miscelazione in base alla similarità di direzione
  final_noise <- numeric(N)
  
  for (i in 1:N) {
    # Direzione di questo punto
    point_dir <- direction_field[i, ]
    
    # Calcola similarità con le direzioni medie dei gruppi
    group_weights <- numeric(n_direction_groups)
    
    for (g in 1:n_direction_groups) {
      # Punti in questo gruppo
      group_points <- which(direction_clusters == g)
      
      if (length(group_points) > 0) {
        # Direzione media del gruppo
        mean_direction <- colMeans(direction_field[group_points, , drop = FALSE])
        mean_direction <- mean_direction / sqrt(sum(mean_direction^2))
        
        # Calcola similarità (prodotto scalare delle direzioni)
        similarity <- sum(point_dir * mean_direction)
        
        # Trasforma in peso (più alto per direzioni simili)
        group_weights[g] <- max(0, similarity)^2
      }
    }
    
    # Normalizza i pesi
    group_weights <- group_weights / sum(group_weights)
    
    # Calcola il valore finale come media pesata
    final_noise[i] <- sum(group_weights * group_noise[i, ])
  }
  
  # 4. Finalizza il risultato
  
  # Normalizza il rumore
  final_noise <- scale(final_noise) * intensity
  
  # Restituisci sia il rumore che il campo di direzioni
  return(list(
    noise = as.vector(final_noise),
    direction_field = direction_field
  ))
}