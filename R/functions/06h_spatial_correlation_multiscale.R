#' Genera pattern di correlazione spaziale multi-scala
#'
#' Crea pattern di correlazione spaziale gerarchici combinando
#' processi gaussiani a scale diverse (micro e macro domini).
#'
#' @param cell_df Dataframe delle celle con coordinate
#' @param spatial_params Parametri spaziali multi-scala
#' @param random_seed Seed per riproducibilità
#' @return Vettore con il rumore spaziale multi-scala generato
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_multiscale_correlation <- function(
  cell_df,
  spatial_params = list(
    # Parametri globali
    use_multiscale = TRUE,
    global_contribution = 0.5,   # Proporzione dell'effetto globale vs locale
    
    # Parametri macro-scala
    macro_range = 100,           # Range del processo a larga scala (macro-domini)
    macro_intensity = 1.0,       # Intensità dell'effetto a larga scala
    macro_model = "Exp",         # Modello di variogramma per macro-scala
    
    # Parametri micro-scala
    micro_range = 20,            # Range del processo a fine scala (micro-domini)
    micro_intensity = 0.8,       # Intensità dell'effetto a piccola scala
    micro_model = "Sph",         # Modello di variogramma per micro-scala
    
    # Parametri per annidamento
    n_hierarchical_levels = 2,   # Numero di livelli gerarchici (default: 2)
    hierarchy_scaling = 0.5      # Fattore di scala tra livelli gerarchici
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Controlla se usare correlazione multi-scala
  use_multiscale <- ifelse(!is.null(spatial_params$use_multiscale), 
                           spatial_params$use_multiscale, TRUE)
  
  if (!use_multiscale) {
    return(NULL)
  }
  
  # Converti cell_df in oggetto spatial
  sp_df <- cell_df
  coordinates(sp_df) <- ~ x + y
  
  # Estrai parametri per macro e micro scala
  global_contribution <- ifelse(!is.null(spatial_params$global_contribution),
                               spatial_params$global_contribution, 0.5)
  
  macro_range <- ifelse(!is.null(spatial_params$macro_range),
                       spatial_params$macro_range, 100)
  
  macro_intensity <- ifelse(!is.null(spatial_params$macro_intensity),
                           spatial_params$macro_intensity, 1.0)
  
  macro_model <- ifelse(!is.null(spatial_params$macro_model),
                       spatial_params$macro_model, "Exp")
  
  micro_range <- ifelse(!is.null(spatial_params$micro_range),
                       spatial_params$micro_range, 20)
  
  micro_intensity <- ifelse(!is.null(spatial_params$micro_intensity),
                           spatial_params$micro_intensity, 0.8)
  
  micro_model <- ifelse(!is.null(spatial_params$micro_model),
                       spatial_params$micro_model, "Sph")
  
  # Parametri per livelli gerarchici
  n_levels <- ifelse(!is.null(spatial_params$n_hierarchical_levels),
                    spatial_params$n_hierarchical_levels, 2)
  
  hierarchy_scaling <- ifelse(!is.null(spatial_params$hierarchy_scaling),
                             spatial_params$hierarchy_scaling, 0.5)
  
  # Limita n_levels a valori ragionevoli (1-5)
  n_levels <- min(max(1, n_levels), 5)
  
  # Genera noise per ogni scala gerarchica
  hierarchical_noise <- list()
  
  # 1. Genera processo a larga scala (macro-domini)
  macro_gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                      beta = 0, model = vgm(psill = macro_intensity,
                                           range = macro_range,
                                           model = macro_model),
                      nmax = 40)
  
  macro_noise <- predict(macro_gp_sim, newdata = sp_df, nsim = 1)$sim1
  hierarchical_noise[[1]] <- scale(macro_noise)
  
  # 2. Genera processi a scala più fine (micro-domini)
  if (n_levels > 1) {
    for (level in 2:n_levels) {
      # Scala il range in base al livello gerarchico
      level_range <- micro_range * (hierarchy_scaling^(level-2))
      level_intensity <- micro_intensity * (1.2^(level-2))
      
      # Genera processo a scala più fine
      micro_gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                          beta = 0, model = vgm(psill = level_intensity,
                                               range = level_range,
                                               model = micro_model),
                          nmax = 20)
      
      micro_noise <- predict(micro_gp_sim, newdata = sp_df, nsim = 1)$sim1
      hierarchical_noise[[level]] <- scale(micro_noise)
    }
  }
  
  # 3. Combina i processi a diverse scale gerarchiche
  combined_noise <- hierarchical_noise[[1]] * global_contribution
  
  if (n_levels > 1) {
    # Calcola i pesi per ogni livello gerarchico
    remaining_weight <- 1 - global_contribution
    level_weights <- remaining_weight * (1 - hierarchy_scaling)^(0:(n_levels-2))
    level_weights <- level_weights / sum(level_weights)
    
    # Combina i livelli gerarchici con i rispettivi pesi
    for (level in 2:n_levels) {
      combined_noise <- combined_noise + 
        hierarchical_noise[[level]] * level_weights[level-1]
    }
  }
  
  # Normalizza il risultato finale
  combined_noise <- scale(combined_noise) * 
    (global_contribution * macro_intensity + (1 - global_contribution) * micro_intensity)
  
  return(as.vector(combined_noise))
}

#' Genera pattern di correlazione spaziale multi-dominio
#'
#' Crea pattern di correlazione spaziale che variaro tra domini
#' utilizzando diverse impostazioni del processo gaussiano.
#'
#' @param cell_df Dataframe delle celle con coordinate
#' @param domain_df Dataframe con la definizione dei domini
#' @param spatial_params Parametri spaziali multi-dominio
#' @param random_seed Seed per riproducibilità
#' @return Vettore con il rumore spaziale combinato per tutti i domini
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_multidomain_correlation <- function(
  cell_df,
  domain_df = NULL,
  spatial_params = list(
    # Parametri globali
    use_multidomain = TRUE,
    
    # Parametri dominio
    n_domains = 3,              # Numero di domini separati (se domain_df non fornito)
    domain_separation = 0.7,    # Chiarezza della separazione tra domini (0-1)
    blend_regions = TRUE,       # Usare transizioni graduali tra domini
    blend_width = 15,           # Ampiezza della zona di transizione
    
    # Parametri GRF per ogni dominio
    domain_params = list(
      list(
        range = 80,            # Primo dominio: grande range
        intensity = 1.2,
        model = "Exp"
      ),
      list(
        range = 20,            # Secondo dominio: piccolo range, alta intensità
        intensity = 2.0,
        model = "Sph"
      ),
      list(
        range = 40,            # Terzo dominio: range medio, bassa intensità
        intensity = 0.6,
        model = "Gau"
      )
    )
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Controlla se usare correlazione multi-dominio
  use_multidomain <- ifelse(!is.null(spatial_params$use_multidomain), 
                           spatial_params$use_multidomain, TRUE)
  
  if (!use_multidomain) {
    return(NULL)
  }
  
  # Converti cell_df in oggetto spatial
  sp_df <- cell_df
  coordinates(sp_df) <- ~ x + y
  coords <- coordinates(sp_df)
  N <- nrow(coords)
  
  # Estrai i parametri
  n_domains <- ifelse(!is.null(spatial_params$n_domains), 
                     spatial_params$n_domains, 3)
  
  domain_separation <- ifelse(!is.null(spatial_params$domain_separation), 
                             spatial_params$domain_separation, 0.7)
  
  blend_regions <- ifelse(!is.null(spatial_params$blend_regions), 
                         spatial_params$blend_regions, TRUE)
  
  blend_width <- ifelse(!is.null(spatial_params$blend_width), 
                       spatial_params$blend_width, 15)
  
  # 1. Determina l'appartenenza di ogni cella a un dominio
  if (is.null(domain_df)) {
    # Genera domini artificiali se non forniti
    
    # Opzione 1: Genera cluster spaziali usando k-means sui dati spaziali
    if (domain_separation >= 0.5) {
      domains <- kmeans(coords, centers = n_domains, nstart = 10)$cluster
    } 
    # Opzione 2: Genera domini da un processo gaussiano per confini più naturali
    else {
      # Genera GP a larga scala
      domain_gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                           beta = 0, model = vgm(psill = 1.0,
                                                range = max(coords[,1]) / 2,
                                                model = "Exp"),
                           nmax = 40)
      
      domain_field <- predict(domain_gp_sim, newdata = sp_df, nsim = 1)$sim1
      domain_field_norm <- domain_field - min(domain_field)
      domain_field_norm <- domain_field_norm / max(domain_field_norm)
      
      # Discretizza il campo per creare domini
      domains <- cut(domain_field_norm, breaks = n_domains, labels = FALSE)
    }
  } else {
    # Usa domini forniti dall'utente
    domains <- domain_df$domain
  }
  
  # 2. Per ogni dominio, genera un campo gaussiano con parametri specifici
  domain_fields <- list()
  for (d in 1:n_domains) {
    # Ottieni parametri per questo dominio
    if (!is.null(spatial_params$domain_params) && 
        length(spatial_params$domain_params) >= d) {
      
      domain_range <- spatial_params$domain_params[[d]]$range
      domain_intensity <- spatial_params$domain_params[[d]]$intensity
      domain_model <- spatial_params$domain_params[[d]]$model
      
    } else {
      # Parametri di default se non specificati
      domain_range <- 40 * (1 + d/5)
      domain_intensity <- 1.0 * (1 + (d %% 3) / 2)
      domain_model <- c("Exp", "Sph", "Gau")[(d %% 3) + 1]
    }
    
    # Genera campo gaussiano per questo dominio
    domain_gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                         beta = 0, model = vgm(psill = domain_intensity,
                                              range = domain_range,
                                              model = domain_model),
                         nmax = 30)
    
    domain_noise <- predict(domain_gp_sim, newdata = sp_df, nsim = 1)$sim1
    domain_fields[[d]] <- scale(domain_noise) * domain_intensity
  }
  
  # 3. Combinare i campi in base all'appartenenza ai domini
  combined_noise <- numeric(N)
  
  if (blend_regions) {
    # Calcola distanze dai confini dei domini
    domain_dists <- matrix(Inf, nrow = N, ncol = n_domains)
    for (d in 1:n_domains) {
      domain_cells <- which(domains == d)
      non_domain_cells <- which(domains != d)
      
      # Per ogni cella nel dominio, calcola la distanza minima alle celle fuori dominio
      if (length(domain_cells) > 0 && length(non_domain_cells) > 0) {
        for (i in domain_cells) {
          domain_dists[i, d] <- min(sqrt(rowSums((coords[non_domain_cells,] - 
                                                 matrix(coords[i,], 
                                                       nrow = length(non_domain_cells), 
                                                       ncol = 2, 
                                                       byrow = TRUE))^2)))
        }
      }
    }
    
    # Normalizza le distanze e calcola i pesi per ogni dominio
    for (i in 1:N) {
      # Calcola i pesi di ogni dominio per questa cella
      weights <- rep(0, n_domains)
      own_domain <- domains[i]
      
      # Distanza dai confini
      distance_to_border <- domain_dists[i, own_domain]
      
      # Confini bruschi (peso 1 per il proprio dominio)
      if (distance_to_border >= blend_width || !blend_regions) {
        weights[own_domain] <- 1
      } 
      # Transizione graduale (peso proporzionale alla distanza dal confine)
      else {
        # Peso base per il proprio dominio
        weights[own_domain] <- 0.5 + 0.5 * (distance_to_border / blend_width)
        
        # Cerca domini vicini
        for (d in 1:n_domains) {
          if (d != own_domain) {
            # Calcola distanza al dominio d
            cells_in_d <- which(domains == d)
            if (length(cells_in_d) > 0) {
              dist_to_d <- min(sqrt(rowSums((coords[cells_in_d,] - 
                                             matrix(coords[i,], 
                                                    nrow = length(cells_in_d), 
                                                    ncol = 2, 
                                                    byrow = TRUE))^2)))
              
              # Se abbastanza vicino, aggiunge influenza
              if (dist_to_d < blend_width) {
                weights[d] <- 0.5 * (1 - dist_to_d / blend_width)
              }
            }
          }
        }
        
        # Normalizza pesi
        weights <- weights / sum(weights)
      }
      
      # Combina i campi in base ai pesi
      for (d in 1:n_domains) {
        if (weights[d] > 0) {
          combined_noise[i] <- combined_noise[i] + domain_fields[[d]][i] * weights[d]
        }
      }
    }
  } else {
    # Versione semplice: assegna ogni cella completamente al suo dominio
    for (i in 1:N) {
      d <- domains[i]
      combined_noise[i] <- domain_fields[[d]][i]
    }
  }
  
  # Normalizza il risultato finale
  combined_noise <- scale(combined_noise)
  
  # Restituisci sia il noise che le etichette di dominio
  return(list(
    noise = as.vector(combined_noise),
    domains = domains
  ))
}