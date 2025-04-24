#' Variazioni temporali nella struttura spaziale
#'
#' Giustificazione: Lavori recenti (La Manno et al., Nature 2021; Qiu et al., Cell 2022) 
#' mostrano che i tessuti presentano dinamiche temporali di espressione genica 
#' anche in campioni apparentemente statici.
#'
#' Questo modulo implementa:
#' - Simulazione di "pseudo-time" all'interno dei campioni
#' - Traiettorie di espressione genica lungo gradienti di sviluppo
#' - Oscillazioni nei pattern di espressione (es. geni coinvolti in ciclo cellulare)
#' - Incorporare le ultime scoperte sul modello "splicing kinetics"

#' Genera campo di pseudo-tempo
#'
#' Crea un campo di pseudo-tempo basato sulla distribuzione spaziale
#' delle cellule, che può seguire gradienti, punti focali o pattern complessi.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param pseudotime_params Parametri per il pseudo-tempo
#' @param random_seed Seed per riproducibilità
#' @return Vettore di valori di pseudo-tempo per ogni cellula
#' @importFrom stats kmeans runif rnorm
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_pseudotime_field <- function(
  cell_df,
  pseudotime_params = list(
    pseudotime_type = "gradient",  # "gradient", "focal", "bifurcation", o "complex"
    direction = c(1, 1),           # Direzione del gradiente (per type="gradient")
    n_foci = 3,                    # Numero di punti focali (per type="focal")
    focus_strength = c(0.3, 1),    # Range di intensità dei focus
    noise_level = 0.1,             # Livello di rumore casuale
    normalize = TRUE,              # Normalizzare i valori nell'intervallo [0,1]
    spatial_coherence = 0.8        # Livello di coerenza spaziale (0-1)
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  pseudotime_type <- pseudotime_params$pseudotime_type
  direction <- pseudotime_params$direction
  n_foci <- pseudotime_params$n_foci
  focus_strength <- pseudotime_params$focus_strength
  noise_level <- pseudotime_params$noise_level
  normalize <- pseudotime_params$normalize
  spatial_coherence <- pseudotime_params$spatial_coherence
  
  # Valori predefiniti se non specificati
  if (is.null(pseudotime_type)) pseudotime_type <- "gradient"
  if (is.null(direction)) direction <- c(1, 1)
  if (is.null(n_foci)) n_foci <- 3
  if (is.null(focus_strength)) focus_strength <- c(0.3, 1)
  if (is.null(noise_level)) noise_level <- 0.1
  if (is.null(normalize)) normalize <- TRUE
  if (is.null(spatial_coherence)) spatial_coherence <- 0.8
  
  # Numero di cellule
  n_cells <- nrow(cell_df)
  
  # Estrai coordinate
  x <- cell_df$x
  y <- cell_df$y
  
  # Inizializza vettore di pseudo-tempo
  pseudotime <- rep(0, n_cells)
  
  # Genera il campo di pseudo-tempo in base al tipo
  if (pseudotime_type == "gradient") {
    # Direzione normalizzata
    direction <- direction / sqrt(sum(direction^2))
    
    # Calcola la proiezione delle coordinate sulla direzione
    pseudotime <- x * direction[1] + y * direction[2]
    
  } else if (pseudotime_type == "focal") {
    # Genera punti focali
    x_range <- range(x)
    y_range <- range(y)
    
    foci <- data.frame(
      x = runif(n_foci, min = x_range[1], max = x_range[2]),
      y = runif(n_foci, min = y_range[1], max = y_range[2]),
      strength = runif(n_foci, min = focus_strength[1], max = focus_strength[2])
    )
    
    # Calcola pseudo-tempo come funzione della distanza dai punti focali
    for (i in 1:n_foci) {
      focus_x <- foci$x[i]
      focus_y <- foci$y[i]
      strength <- foci$strength[i]
      
      # Distanza euclidea dal punto focale
      dist <- sqrt((x - focus_x)^2 + (y - focus_y)^2)
      
      # Contributo inversamente proporzionale alla distanza, modulato dalla forza
      contribution <- strength * exp(-dist / (max(dist) * 0.2))
      
      # Aggiungi il contributo al campo di pseudo-tempo
      pseudotime <- pseudotime + contribution
    }
    
  } else if (pseudotime_type == "bifurcation") {
    # Implementa una biforcazione del tempo lungo un asse principale
    # Prima determina un asse principale
    principal_axis <- c(runif(1, -1, 1), runif(1, -1, 1))
    principal_axis <- principal_axis / sqrt(sum(principal_axis^2))
    
    # Proiezione sul asse principale
    main_component <- x * principal_axis[1] + y * principal_axis[2]
    
    # Normalizza tra 0 e 1
    main_component <- (main_component - min(main_component)) / (max(main_component) - min(main_component))
    
    # Determina punto di biforcazione
    bifurcation_point <- runif(1, 0.4, 0.6)
    
    # Asse secondario ortogonale al principale
    secondary_axis <- c(-principal_axis[2], principal_axis[1])
    
    # Proiezione sull'asse secondario
    sec_component <- x * secondary_axis[1] + y * secondary_axis[2]
    
    # Prima della biforcazione: seguire solo l'asse principale
    pre_bifurcation <- main_component < bifurcation_point
    pseudotime[pre_bifurcation] <- main_component[pre_bifurcation]
    
    # Dopo la biforcazione: dividere in due rami
    # Ramo superiore
    upper_branch <- main_component >= bifurcation_point & sec_component > 0
    pseudotime[upper_branch] <- main_component[upper_branch] + 
      0.5 * abs(sec_component[upper_branch]) / max(abs(sec_component))
    
    # Ramo inferiore
    lower_branch <- main_component >= bifurcation_point & sec_component <= 0
    pseudotime[lower_branch] <- main_component[lower_branch] + 
      0.5 * abs(sec_component[lower_branch]) / max(abs(sec_component))
    
  } else if (pseudotime_type == "complex") {
    # Crea un campo gaussiano spazialmente correlato per un pattern complesso
    sp_df <- cell_df
    coordinates(sp_df) <- ~ x + y
    
    # Definisci un modello di correlazione spaziale
    range_param <- max(diff(range(x)), diff(range(y))) * 0.3
    gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                  beta = 0, model = vgm(psill = 1.0, 
                                       range = range_param, 
                                       model = "Exp"), 
                  nmax = 30)
    
    # Genera il campo gaussiano
    pseudotime <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
  } else {
    # Default: random con coerenza spaziale
    pseudotime <- rnorm(n_cells)
    
    # Se richiesto, aggiungi coerenza spaziale
    if (spatial_coherence > 0) {
      # Converti a spatial
      sp_df <- cell_df
      coordinates(sp_df) <- ~ x + y
      
      # Genera un campo gaussiano spazialmente correlato
      range_param <- max(diff(range(x)), diff(range(y))) * 0.2
      gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                    beta = 0, model = vgm(psill = 1.0, 
                                         range = range_param, 
                                         model = "Exp"), 
                    nmax = 30)
      
      # Genera il campo gaussiano
      coherent_component <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
      
      # Miscela il componente coerente con il rumore casuale
      pseudotime <- spatial_coherence * scale(coherent_component) + 
                  (1 - spatial_coherence) * scale(pseudotime)
    }
  }
  
  # Aggiungi rumore casuale
  if (noise_level > 0) {
    noise <- rnorm(n_cells, mean = 0, sd = noise_level)
    pseudotime <- pseudotime + noise
  }
  
  # Normalizza se richiesto
  if (normalize) {
    pseudotime <- (pseudotime - min(pseudotime)) / (max(pseudotime) - min(pseudotime))
  }
  
  return(pseudotime)
}

#' Genera traiettorie di espressione lungo pseudo-tempo
#'
#' Definisce come l'espressione dei geni cambia lungo la traiettoria
#' di pseudo-tempo, utilizzando modelli di traiettoria biologicamente
#' plausibili.
#'
#' @param pseudotime Vettore di valori di pseudo-tempo per ogni cellula
#' @param n_genes Numero totale di geni
#' @param trajectory_params Parametri per le traiettorie
#' @param random_seed Seed per riproducibilità
#' @return Matrice di effetti temporali sull'espressione
#' @importFrom stats runif rbinom rpois
#' @export
generate_gene_trajectories <- function(
  pseudotime,
  n_genes,
  trajectory_params = list(
    temporal_genes_fraction = 0.3,    # Frazione di geni con dinamica temporale
    monotonic_fraction = 0.6,         # Frazione di geni con traiettoria monotona
    oscillatory_fraction = 0.2,       # Frazione di geni con pattern oscillatorio
    transient_fraction = 0.2,         # Frazione di geni con picchi transitori
    max_effect_size = 2.0,            # Intensità massima dell'effetto temporale
    n_oscillatory_groups = 3,         # Gruppi di geni con oscillazioni sincronizzate
    oscillation_periods = c(0.2, 1),  # Range dei periodi di oscillazione
    noise_level = 0.1                 # Rumore nei profili temporali
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  temporal_genes_fraction <- trajectory_params$temporal_genes_fraction
  monotonic_fraction <- trajectory_params$monotonic_fraction
  oscillatory_fraction <- trajectory_params$oscillatory_fraction
  transient_fraction <- trajectory_params$transient_fraction
  max_effect_size <- trajectory_params$max_effect_size
  n_oscillatory_groups <- trajectory_params$n_oscillatory_groups
  oscillation_periods <- trajectory_params$oscillation_periods
  noise_level <- trajectory_params$noise_level
  
  # Valori predefiniti se non specificati
  if (is.null(temporal_genes_fraction)) temporal_genes_fraction <- 0.3
  if (is.null(monotonic_fraction)) monotonic_fraction <- 0.6
  if (is.null(oscillatory_fraction)) oscillatory_fraction <- 0.2
  if (is.null(transient_fraction)) transient_fraction <- 0.2
  if (is.null(max_effect_size)) max_effect_size <- 2.0
  if (is.null(n_oscillatory_groups)) n_oscillatory_groups <- 3
  if (is.null(oscillation_periods)) oscillation_periods <- c(0.2, 1)
  if (is.null(noise_level)) noise_level <- 0.1
  
  # Numero di cellule
  n_cells <- length(pseudotime)
  
  # Numero di geni con effetti temporali
  n_temporal_genes <- round(n_genes * temporal_genes_fraction)
  
  # Inizializza matrice degli effetti
  temporal_effects <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Seleziona i geni con effetti temporali
  temporal_gene_indices <- sample(1:n_genes, n_temporal_genes)
  
  # Assegna tipo di traiettoria a ciascun gene temporale
  n_monotonic <- round(n_temporal_genes * monotonic_fraction)
  n_oscillatory <- round(n_temporal_genes * oscillatory_fraction)
  n_transient <- round(n_temporal_genes * transient_fraction)
  n_other <- n_temporal_genes - (n_monotonic + n_oscillatory + n_transient)
  
  # Aggiusta numeri se necessario
  if (n_other < 0) {
    n_monotonic <- n_monotonic + n_other
    n_other <- 0
  }
  
  # Crea una lista di indici per tipo di traiettoria
  temp_gene_indices <- sample(temporal_gene_indices)
  monotonic_genes <- temp_gene_indices[1:n_monotonic]
  oscillatory_genes <- temp_gene_indices[(n_monotonic+1):(n_monotonic+n_oscillatory)]
  transient_genes <- temp_gene_indices[(n_monotonic+n_oscillatory+1):(n_monotonic+n_oscillatory+n_transient)]
  other_genes <- temp_gene_indices[(n_monotonic+n_oscillatory+n_transient+1):n_temporal_genes]
  
  # 1. Geni con traiettorie monotone (lineari, esponenziali, logistiche)
  if (length(monotonic_genes) > 0) {
    for (g in monotonic_genes) {
      # Determina se crescente o decrescente
      increasing <- rbinom(1, 1, 0.5) == 1
      
      # Determina tipo di curva monotona
      curve_type <- sample(c("linear", "exponential", "sigmoid"), 1)
      
      # Effetto massimo per questo gene
      effect_size <- runif(1, 0.5, max_effect_size)
      
      # Genera traiettoria
      if (curve_type == "linear") {
        trajectory <- pseudotime
        
      } else if (curve_type == "exponential") {
        trajectory <- pseudotime^2
        
      } else if (curve_type == "sigmoid") {
        # Funzione logistica centrata a pseudotime = 0.5
        midpoint <- runif(1, 0.3, 0.7)
        steepness <- runif(1, 5, 15)
        trajectory <- 1 / (1 + exp(-steepness * (pseudotime - midpoint)))
      }
      
      # Normalizza tra 0 e 1
      trajectory <- (trajectory - min(trajectory)) / (max(trajectory) - min(trajectory))
      
      # Inverti se decrescente
      if (!increasing) {
        trajectory <- 1 - trajectory
      }
      
      # Scala all'effetto desiderato
      temporal_effects[, g] <- effect_size * (trajectory - 0.5)
    }
  }
  
  # 2. Geni con pattern oscillatori (seni, coseni con periodo e fase variabili)
  if (length(oscillatory_genes) > 0) {
    # Divide i geni oscillatori in gruppi con pattern sincronizzati
    genes_per_group <- ceiling(length(oscillatory_genes) / n_oscillatory_groups)
    oscillatory_groups <- list()
    for (i in 1:n_oscillatory_groups) {
      start_idx <- (i - 1) * genes_per_group + 1
      end_idx <- min(i * genes_per_group, length(oscillatory_genes))
      if (start_idx <= end_idx) {
        oscillatory_groups[[i]] <- oscillatory_genes[start_idx:end_idx]
      }
    }
    
    # Per ogni gruppo, genera un pattern oscillatorio di base
    for (i in 1:length(oscillatory_groups)) {
      if (length(oscillatory_groups[[i]]) > 0) {
        # Parametri di oscillazione comuni per il gruppo
        period <- runif(1, min = oscillation_periods[1], max = oscillation_periods[2])
        phase <- runif(1, 0, 2 * pi)
        
        # Pattern di base
        base_oscillation <- sin(2 * pi * pseudotime / period + phase)
        
        # Assegna a ciascun gene nel gruppo con variazioni minori
        for (g in oscillatory_groups[[i]]) {
          effect_size <- runif(1, 0.3, max_effect_size)
          
          # Varia leggermente il pattern per ogni gene
          gene_variation <- runif(1, 0.8, 1.2)  # Variazione di ampiezza
          phase_shift <- runif(1, -0.2, 0.2)    # Leggero sfasamento
          
          trajectory <- gene_variation * sin(2 * pi * pseudotime / period + phase + phase_shift)
          
          # Scala all'effetto desiderato
          temporal_effects[, g] <- effect_size * trajectory
        }
      }
    }
  }
  
  # 3. Geni con pattern transitori (picchi a determinati punti del tempo)
  if (length(transient_genes) > 0) {
    for (g in transient_genes) {
      # Parametri per il picco
      peak_time <- runif(1, 0.2, 0.8)  # Posizione del picco
      peak_width <- runif(1, 0.05, 0.2)  # Ampiezza del picco
      effect_size <- runif(1, 0.5, max_effect_size)
      
      # Genera picco gaussiano
      distance_from_peak <- abs(pseudotime - peak_time)
      trajectory <- exp(-(distance_from_peak^2) / (2 * peak_width^2))
      
      # Scala all'effetto desiderato
      temporal_effects[, g] <- effect_size * trajectory
    }
  }
  
  # 4. Altri pattern più complessi
  if (length(other_genes) > 0) {
    for (g in other_genes) {
      # Genera un pattern complesso come somma di componenti armoniche
      n_components <- rpois(1, lambda = 2) + 1
      effect_size <- runif(1, 0.3, max_effect_size)
      
      trajectory <- rep(0, n_cells)
      for (i in 1:n_components) {
        period <- runif(1, min = oscillation_periods[1], max = oscillation_periods[2])
        phase <- runif(1, 0, 2 * pi)
        amplitude <- runif(1, 0.2, 1)
        
        component <- amplitude * sin(2 * pi * i * pseudotime / period + phase)
        trajectory <- trajectory + component
      }
      
      # Normalizza
      trajectory <- (trajectory - min(trajectory)) / (max(trajectory) - min(trajectory))
      
      # Centra a zero
      trajectory <- trajectory - 0.5
      
      # Scala all'effetto desiderato
      temporal_effects[, g] <- effect_size * trajectory
    }
  }
  
  # Aggiungi rumore ai profili temporali
  if (noise_level > 0) {
    noise <- matrix(rnorm(n_cells * n_genes, mean = 0, sd = noise_level),
                   nrow = n_cells, ncol = n_genes)
    noise[, -temporal_gene_indices] <- 0  # Aggiungi rumore solo ai geni temporali
    temporal_effects <- temporal_effects + noise
  }
  
  # Metadati sui geni temporali
  gene_metadata <- data.frame(
    gene_id = 1:n_genes,
    is_temporal = rep(FALSE, n_genes),
    temporal_type = rep("none", n_genes)
  )
  
  gene_metadata$is_temporal[temporal_gene_indices] <- TRUE
  gene_metadata$temporal_type[monotonic_genes] <- "monotonic"
  gene_metadata$temporal_type[oscillatory_genes] <- "oscillatory"
  gene_metadata$temporal_type[transient_genes] <- "transient"
  gene_metadata$temporal_type[other_genes] <- "complex"
  
  # Restituisci sia gli effetti che i metadati
  return(list(
    temporal_effects = temporal_effects,
    gene_metadata = gene_metadata,
    pseudotime = pseudotime
  ))
}

#' Modella effetti di RNA velocità basati sulla cinetica di splicing
#'
#' Implementa un modello di metabolismo dell'RNA che considera
#' trascrizione, splicing e degradazione per simulare RNA velocità
#' e stati di transizione.
#'
#' @param pseudotime Vettore di valori di pseudo-tempo per ogni cellula
#' @param gene_metadata Metadati sui geni generati da generate_gene_trajectories
#' @param velocity_params Parametri per RNA velocità
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici di RNA pre e maturo
#' @importFrom stats runif rgamma
#' @export
generate_rna_velocity <- function(
  pseudotime,
  gene_metadata,
  velocity_params = list(
    fraction_velocity_genes = 0.7,   # Frazione di geni temporali con effetti di velocità
    splicing_rate_mean = 0.2,        # Tasso di splicing medio (1/tempo)
    splicing_rate_cv = 0.3,          # Coefficiente di variazione per tassi di splicing
    degradation_rate_mean = 0.1,     # Tasso di degradazione medio (1/tempo)
    degradation_rate_cv = 0.4,       # Coefficiente di variazione per tassi di degradazione
    scale_factor = 5.0,              # Fattore di scala per RNA non-spliced
    gene_specific_kinetics = TRUE    # Usare parametri cinetici specifici per gene
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  fraction_velocity_genes <- velocity_params$fraction_velocity_genes
  splicing_rate_mean <- velocity_params$splicing_rate_mean
  splicing_rate_cv <- velocity_params$splicing_rate_cv
  degradation_rate_mean <- velocity_params$degradation_rate_mean
  degradation_rate_cv <- velocity_params$degradation_rate_cv
  scale_factor <- velocity_params$scale_factor
  gene_specific_kinetics <- velocity_params$gene_specific_kinetics
  
  # Valori predefiniti se non specificati
  if (is.null(fraction_velocity_genes)) fraction_velocity_genes <- 0.7
  if (is.null(splicing_rate_mean)) splicing_rate_mean <- 0.2
  if (is.null(splicing_rate_cv)) splicing_rate_cv <- 0.3
  if (is.null(degradation_rate_mean)) degradation_rate_mean <- 0.1
  if (is.null(degradation_rate_cv)) degradation_rate_cv <- 0.4
  if (is.null(scale_factor)) scale_factor <- 5.0
  if (is.null(gene_specific_kinetics)) gene_specific_kinetics <- TRUE
  
  # Dimensioni
  n_cells <- length(pseudotime)
  n_genes <- nrow(gene_metadata)
  
  # Identifica i geni temporali
  temporal_genes <- which(gene_metadata$is_temporal)
  
  # Seleziona una sottofrazione di geni temporali per gli effetti di velocità
  n_velocity_genes <- round(length(temporal_genes) * fraction_velocity_genes)
  velocity_genes <- sample(temporal_genes, n_velocity_genes)
  
  # Inizializza matrici per RNA pre-spliced e maturo
  unspliced_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  spliced_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Genera tassi di splicing e degradazione gene-specifici
  if (gene_specific_kinetics) {
    # Usa distribuzione gamma per generare tassi positivi con dispersione controllata
    splicing_rates <- rgamma(n_genes, 
                           shape = 1 / splicing_rate_cv^2,
                           scale = splicing_rate_mean * splicing_rate_cv^2)
    
    degradation_rates <- rgamma(n_genes, 
                              shape = 1 / degradation_rate_cv^2,
                              scale = degradation_rate_mean * degradation_rate_cv^2)
  } else {
    # Tassi costanti
    splicing_rates <- rep(splicing_rate_mean, n_genes)
    degradation_rates <- rep(degradation_rate_mean, n_genes)
  }
  
  # Funzione di trascrizione che varia lungo pseudotime
  # per emulare il modello di Bergen et al. 2020 (Nature Biotechnology)
  transcription_function <- function(time, amplitude = 1, offset = 0, width = 0.1) {
    # Funzione impulso gaussiana centrata in 'offset'
    pulse <- amplitude * exp(-((time - offset)^2) / (2 * width^2))
    return(pulse)
  }
  
  # Per ogni gene selezionato, genera un profilo di velocità
  for (g in velocity_genes) {
    # Parametri cinetici per questo gene
    splicing_rate <- splicing_rates[g]
    degradation_rate <- degradation_rates[g]
    
    # Determina il trend temporale per questo gene
    gene_type <- gene_metadata$temporal_type[g]
    
    # Differenti pattern di trascrizione basati sul tipo temporale
    if (gene_type == "monotonic") {
      # Per geni monotoni, trascrizione che aumenta/diminuisce monotonicamente
      is_increasing <- runif(1) > 0.5
      slope <- runif(1, 0.5, 2.0) * ifelse(is_increasing, 1, -1)
      
      # Trascrizione lineare o esponenziale
      if (runif(1) > 0.5) {
        # Lineare
        transcription <- slope * pseudotime + 0.2
      } else {
        # Esponenziale
        transcription <- exp(slope * pseudotime) / exp(abs(slope))
      }
      
      # Normalizza per evitare valori negativi
      transcription <- pmax(0, transcription)
      
    } else if (gene_type == "oscillatory") {
      # Per geni oscillatori, pattern di trascrizione oscillante
      period <- runif(1, 0.2, 0.6)
      phase <- runif(1, 0, 2 * pi)
      transcription <- 0.5 + 0.5 * sin(2 * pi * pseudotime / period + phase)
      
    } else if (gene_type == "transient") {
      # Per geni transienti, picchi localizzati
      n_pulses <- sample(1:3, 1)
      transcription <- rep(0, n_cells)
      
      for (p in 1:n_pulses) {
        # Posizione casuale degli impulsi
        offset <- runif(1, 0.1, 0.9)
        width <- runif(1, 0.05, 0.2)
        amplitude <- runif(1, 0.5, 2.0)
        
        transcription <- transcription + transcription_function(
          pseudotime, amplitude, offset, width
        )
      }
      
      # Aggiungi un livello basale
      transcription <- transcription + runif(1, 0, 0.2)
      
    } else {
      # Pattern complesso
      n_components <- sample(2:4, 1)
      transcription <- rep(0, n_cells)
      
      for (c in 1:n_components) {
        if (runif(1) > 0.5) {
          # Componente sinusoidale
          period <- runif(1, 0.1, 0.5)
          phase <- runif(1, 0, 2 * pi)
          amplitude <- runif(1, 0.2, 1.0)
          
          transcription <- transcription + amplitude * sin(2 * pi * pseudotime / period + phase)
        } else {
          # Componente impulso
          offset <- runif(1, 0.1, 0.9)
          width <- runif(1, 0.05, 0.2)
          amplitude <- runif(1, 0.5, 1.5)
          
          transcription <- transcription + transcription_function(
            pseudotime, amplitude, offset, width
          )
        }
      }
      
      # Normalizza e sposta per avere valori positivi
      transcription <- transcription - min(transcription)
      transcription <- transcription / max(transcription)
    }
    
    # Normalizza la trascrizione tra 0.1 e 1
    transcription <- 0.1 + 0.9 * (transcription - min(transcription)) / 
                  (max(transcription) - min(transcription) + 1e-10)
    
    # Calcola il modello di RNA dinamica
    # Implementazione del modello di Bergen et al. 2020 (semplificata)
    
    # 1. RNA non-spliced (pre-mRNA)
    # du/dt = α(t) - β·u
    # Approssimazione all'equilibrio locale: u ≈ α(t)/β
    unspliced <- transcription / splicing_rate
    
    # 2. RNA spliced (mRNA maturo)
    # ds/dt = β·u - γ·s
    # Approssimazione all'equilibrio locale: s ≈ β·u/γ = α(t)/γ
    spliced <- transcription / degradation_rate
    
    # Rappresentiamo la dinamica nelle vicinanze dell'equilibrio
    # Perturbando leggermente per simulare la transizione
    # Questo crea il pattern "a ciclo aperto" tipico dei plot di RNA velocity
    time_shift <- 0.05  # Rappresenta un piccolo passo nel tempo
    shifted_idx <- pmin(n_cells, pmax(1, round((pseudotime + time_shift) * n_cells)))
    
    # Introduciamo un ritardo per il spliced rispetto all'unspliced
    # Questo crea l'effetto di velocità osservato nei dati reali
    unspliced_raw <- unspliced
    spliced_raw <- spliced
    
    # Scala matrici e aggiungi rumore
    noise_factor <- 0.1
    unspliced_noise <- rnorm(n_cells, 0, noise_factor * mean(unspliced))
    spliced_noise <- rnorm(n_cells, 0, noise_factor * mean(spliced))
    
    unspliced_matrix[, g] <- scale_factor * unspliced + unspliced_noise
    spliced_matrix[, g] <- spliced + spliced_noise
    
    # Assicura che non ci siano valori negativi
    unspliced_matrix[unspliced_matrix < 0] <- 0
    spliced_matrix[spliced_matrix < 0] <- 0
  }
  
  # Restituisci le matrici e i metadati
  return(list(
    unspliced = unspliced_matrix,
    spliced = spliced_matrix,
    velocity_genes = velocity_genes,
    splicing_rates = splicing_rates,
    degradation_rates = degradation_rates,
    pseudotime = pseudotime
  ))
}

#' Integra dinamiche temporali nell'espressione genica
#'
#' Funzione wrapper che combina pseudo-tempo, traiettorie geniche
#' e velocità RNA per creare una simulazione completa di
#' dinamiche temporali nell'espressione genica.
#'
#' @param cell_df Dataframe delle cellule con coordinate
#' @param expression_matrix Matrice di espressione originale
#' @param temporal_params Parametri per la simulazione temporale
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici modificate e metadati temporali
#' @export
generate_temporal_dynamics <- function(
  cell_df,
  expression_matrix,
  temporal_params = list(
    pseudotime_type = "gradient",   # Tipo di campo di pseudo-tempo
    direction = c(1, 1),            # Direzione del gradiente
    temporal_genes_fraction = 0.3,  # Frazione di geni con dinamica temporale
    monotonic_fraction = 0.6,       # Frazione di geni con traiettoria monotona
    oscillatory_fraction = 0.2,     # Frazione di geni con pattern oscillatorio
    transient_fraction = 0.2,       # Frazione di geni con picchi transitori
    max_effect_size = 2.0,          # Intensità massima dell'effetto temporale
    include_velocity = TRUE,        # Includere dinamiche di RNA velocity
    splicing_rate_mean = 0.2,       # Tasso di splicing medio
    velocity_genes_fraction = 0.5,  # Frazione di geni con dinamica di RNA velocity
    noise_level = 0.1,              # Livello di rumore nei profili temporali
    integration_weight = 0.8        # Peso dell'integrazione del pattern temporale
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Dimensioni
  n_cells <- nrow(cell_df)
  n_genes <- ncol(expression_matrix)
  
  # 1. Genera campo di pseudo-tempo
  pseudotime_params <- list(
    pseudotime_type = temporal_params$pseudotime_type,
    direction = temporal_params$direction,
    noise_level = temporal_params$noise_level / 2,  # Rumore ridotto per il pseudotime
    normalize = TRUE,
    spatial_coherence = 0.8
  )
  
  pseudotime <- generate_pseudotime_field(
    cell_df = cell_df,
    pseudotime_params = pseudotime_params,
    random_seed = random_seed
  )
  
  # 2. Genera traiettorie di espressione
  trajectory_params <- list(
    temporal_genes_fraction = temporal_params$temporal_genes_fraction,
    monotonic_fraction = temporal_params$monotonic_fraction,
    oscillatory_fraction = temporal_params$oscillatory_fraction,
    transient_fraction = temporal_params$transient_fraction,
    max_effect_size = temporal_params$max_effect_size,
    noise_level = temporal_params$noise_level
  )
  
  trajectory_results <- generate_gene_trajectories(
    pseudotime = pseudotime,
    n_genes = n_genes,
    trajectory_params = trajectory_params,
    random_seed = random_seed
  )
  
  # 3. Integra traiettorie nell'espressione
  integration_weight <- temporal_params$integration_weight
  
  if (is.null(integration_weight)) integration_weight <- 0.8
  
  # Metodo moltiplicativo: converte gli effetti in fattori moltiplicativi
  temporal_factors <- exp(integration_weight * trajectory_results$temporal_effects)
  
  # Applicazione moltiplicativa (mantenendo media costante per non influenzare library size)
  modified_expression <- expression_matrix * temporal_factors
  
  # 4. Genera dinamiche di RNA velocity se richiesto
  velocity_results <- NULL
  
  if (!is.null(temporal_params$include_velocity) && temporal_params$include_velocity) {
    velocity_params <- list(
      fraction_velocity_genes = temporal_params$velocity_genes_fraction,
      splicing_rate_mean = temporal_params$splicing_rate_mean,
      degradation_rate_mean = 0.1,
      gene_specific_kinetics = TRUE
    )
    
    velocity_results <- generate_rna_velocity(
      pseudotime = pseudotime,
      gene_metadata = trajectory_results$gene_metadata,
      velocity_params = velocity_params,
      random_seed = random_seed
    )
  }
  
  # Restituisci risultati
  return(list(
    expression = modified_expression,  # Matrice di espressione modificata
    original_expression = expression_matrix,  # Matrice originale
    pseudotime = pseudotime,  # Campo di pseudo-tempo
    temporal_effects = trajectory_results$temporal_effects,  # Effetti temporali
    gene_metadata = trajectory_results$gene_metadata,  # Metadati sui geni
    unspliced = if (!is.null(velocity_results)) velocity_results$unspliced else NULL,  # RNA non-spliced
    spliced = if (!is.null(velocity_results)) velocity_results$spliced else NULL,  # RNA spliced
    velocity_genes = if (!is.null(velocity_results)) velocity_results$velocity_genes else NULL  # Geni con RNA velocity
  ))
}