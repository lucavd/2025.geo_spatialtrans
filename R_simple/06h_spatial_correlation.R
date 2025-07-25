#' Genera effetti di correlazione spaziale
#'
#' Crea pattern di correlazione spaziale usando processi gaussiani
#' o altri metodi, con supporto per modelli avanzati multi-scala e non-stazionari.
#'
#' @param cell_df Dataframe delle celle con coordinate
#' @param spatial_params Parametri spaziali
#' @param correlation_method Metodo di correlazione ("grf", "car", "multiscale", "nonstationary", "anisotropic", "multidomain")
#' @param random_seed Seed per riproducibilità
#' @return Vettore o lista con il rumore spaziale generato e metadati
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_spatial_correlation <- function(
  cell_df,
  spatial_params = list(
    spatial_noise_intensity = 1.0,
    spatial_range = 30,
    
    # Parametri per modelli avanzati
    use_multiscale = FALSE,
    use_nonstationary = FALSE,
    use_anisotropy = FALSE,
    use_multidomain = FALSE
  ),
  correlation_method = "grf",
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Metodo di correlazione avanzato richiesto?
  advanced_methods <- c("multiscale", "nonstationary", "anisotropic", "multidomain")
  
  # Se è richiesto un metodo avanzato, usa le funzioni specifiche
  if (correlation_method %in% advanced_methods) {
    result <- switch(
      correlation_method,
      "multiscale" = generate_multiscale_correlation(cell_df, spatial_params, random_seed),
      "nonstationary" = generate_nonstationary_correlation(cell_df, spatial_params, random_seed),
      "anisotropic" = generate_anisotropic_correlation(cell_df, spatial_params, random_seed),
      "multidomain" = generate_multidomain_correlation(cell_df, spatial_params, random_seed)
    )
    
    # Se è stato restituito un oggetto lista, estrai solo il noise
    if (is.list(result) && "noise" %in% names(result)) {
      return(result$noise)
    } else {
      return(result)
    }
  }
  
  # Per metodi tradizionali (grf, car) continua con l'implementazione originale
  
  # Inizializza il risultato
  gp_noise <- NULL
  
  # Procedi solo se abbiamo un metodo di correlazione
  if (correlation_method %in% c("grf", "car")) {
    # Converti cell_df in oggetto spatial
    sp_df <- cell_df
    coordinates(sp_df) <- ~ x + y
    
    if (correlation_method == "grf") {
      # Verifica se usare multi-scala
      if (!is.null(spatial_params$use_multiscale) && spatial_params$use_multiscale) {
        # Usa pattern multi-scala
        result <- generate_multiscale_correlation(cell_df, spatial_params, random_seed)
        if (is.list(result)) {
          gp_noise <- result$noise
        } else {
          gp_noise <- result
        }
      } else {
        # Gaussian Random Field standard
        gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                       beta = 0, model = vgm(psill = spatial_params$spatial_noise_intensity * 2,
                                            range = spatial_params$spatial_range * 0.8,
                                            model = "Exp"),
                       nmax = 20)
        
        # Genera il noise spaziale
        gp_noise <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
        
        # Normalizza il noise
        gp_noise <- scale(gp_noise) * 1.5
      }
    } else if (correlation_method == "car") {
      # Conditional Autoregressive Model - implementazione semplificata
      
      # Verifica se usare anisotropia o non-stazionarietà
      if ((!is.null(spatial_params$use_anisotropy) && spatial_params$use_anisotropy) ||
          (!is.null(spatial_params$use_nonstationary) && spatial_params$use_nonstationary)) {
        # Usa versione avanzata basata su metodo appropriato
        if (!is.null(spatial_params$use_anisotropy) && spatial_params$use_anisotropy) {
          result <- generate_anisotropic_correlation(cell_df, spatial_params, random_seed)
        } else {
          result <- generate_nonstationary_correlation(cell_df, spatial_params, random_seed)
        }
        
        if (is.list(result)) {
          gp_noise <- result$noise
        } else {
          gp_noise <- result
        }
      } else {
        # Versione originale del CAR
        coords <- coordinates(sp_df)
        dist_mat <- as.matrix(dist(coords))
        
        # Crea una matrice di pesi per i vicini
        threshold_dist <- spatial_params$spatial_range * 0.5
        W <- (dist_mat <= threshold_dist) * (1 - dist_mat/threshold_dist)
        diag(W) <- 0
        
        # Normalizza i pesi
        W <- sweep(W, 1, rowSums(W) + 1e-10, "/")
        
        # Genera rumore base
        base_noise <- rnorm(nrow(coords))
        
        # Applica effetto CAR (versione semplificata)
        spatial_effect <- 0.8  # Forza dell'effetto spaziale
        
        # Iterazioni per convergenza
        noise <- base_noise
        for (i in 1:5) {
          noise <- (1 - spatial_effect) * base_noise + spatial_effect * (W %*% noise)
        }
        
        # Normalizza
        gp_noise <- as.vector(scale(noise)) * spatial_params$spatial_noise_intensity
      }
    }
  }
  
  return(gp_noise)
}

#' Genera pattern spaziali compositi combinando diversi metodi
#'
#' Crea pattern di correlazione spaziale complessi sovrapponendo
#' diversi tipi di correlazione spaziale (ad es. multiscala + anisotropia).
#'
#' @param cell_df Dataframe delle celle con coordinate
#' @param spatial_params Parametri spaziali compositi
#' @param random_seed Seed per riproducibilità
#' @return Vettore con il rumore spaziale composito generato
#' @importFrom sp coordinates
#' @export
generate_composite_spatial_patterns <- function(
  cell_df,
  spatial_params = list(
    # Parametri globali
    use_composite = TRUE,
    composite_methods = c("multiscale", "anisotropic"),  # Metodi da combinare
    method_weights = c(0.6, 0.4),                        # Pesi relativi
    
    # Parametri specifici per ogni metodo
    multiscale_params = list(
      use_multiscale = TRUE,
      global_contribution = 0.7,
      macro_range = 80,
      micro_range = 15,
      n_hierarchical_levels = 3
    ),
    
    anisotropic_params = list(
      use_anisotropy = TRUE,
      anisotropy_type = "flow",
      flow_field_type = "gradient",
      range_main = 40,
      range_secondary = 10
    ),
    
    nonstationary_params = list(
      use_nonstationary = TRUE,
      nonstationary_type = "patch",
      n_patches = 4,
      blend_patches = TRUE
    ),
    
    multidomain_params = list(
      use_multidomain = TRUE,
      n_domains = 3,
      blend_regions = TRUE
    )
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Controlla se usare pattern compositi
  use_composite <- ifelse(!is.null(spatial_params$use_composite), spatial_params$use_composite, TRUE)
  
  if (!use_composite) {
    return(NULL)
  }
  
  # Definisci i metodi da combinare
  composite_methods <- ifelse(!is.null(spatial_params$composite_methods),
                             list(spatial_params$composite_methods),
                             list(c("multiscale", "anisotropic")))[[1]]
  
  # Definisci i pesi per ogni metodo
  method_weights <- ifelse(!is.null(spatial_params$method_weights),
                          list(spatial_params$method_weights),
                          list(rep(1/length(composite_methods), length(composite_methods))))[[1]]
  
  # Assicurati che i pesi sommino a 1
  method_weights <- method_weights / sum(method_weights)
  
  # Genera pattern per ogni metodo richiesto
  pattern_list <- list()
  for (method in composite_methods) {
    # Ottieni i parametri specifici per il metodo
    method_params <- switch(
      method,
      "multiscale" = spatial_params$multiscale_params,
      "anisotropic" = spatial_params$anisotropic_params,
      "nonstationary" = spatial_params$nonstationary_params,
      "multidomain" = spatial_params$multidomain_params,
      NULL
    )
    
    # Se i parametri non sono specificati, usa un set di default
    if (is.null(method_params)) {
      method_params <- list()
    }
    
    # Attiva l'uso del metodo
    method_params[[paste0("use_", method)]] <- TRUE
    
    # Genera il pattern usando il metodo appropriato
    pattern <- switch(
      method,
      "multiscale" = generate_multiscale_correlation(cell_df, method_params, random_seed + 1),
      "anisotropic" = generate_anisotropic_correlation(cell_df, method_params, random_seed + 2),
      "nonstationary" = generate_nonstationary_correlation(cell_df, method_params, random_seed + 3),
      "multidomain" = generate_multidomain_correlation(cell_df, method_params, random_seed + 4),
      NULL
    )
    
    # Estrai il rumore dalla lista se necessario
    if (is.list(pattern) && "noise" %in% names(pattern)) {
      pattern_list[[method]] <- pattern$noise
    } else if (!is.null(pattern)) {
      pattern_list[[method]] <- pattern
    }
  }
  
  # Combina i pattern generati usando i pesi specificati
  composite_pattern <- numeric(nrow(cell_df))
  
  # Pesi attuali basati sui metodi disponibili
  available_methods <- names(pattern_list)
  available_weights <- method_weights[match(available_methods, composite_methods)]
  available_weights <- available_weights / sum(available_weights)
  
  # Applica la combinazione pesata
  for (i in seq_along(available_methods)) {
    method <- available_methods[i]
    weight <- available_weights[i]
    
    # Normalizza il pattern prima di combinarlo
    normalized_pattern <- scale(pattern_list[[method]])
    
    # Aggiungi al pattern composito
    composite_pattern <- composite_pattern + weight * normalized_pattern
  }
  
  # Normalizza il pattern finale
  composite_pattern <- scale(composite_pattern)
  
  return(as.vector(composite_pattern))
}