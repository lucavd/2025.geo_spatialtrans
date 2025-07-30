#' Genera dimensioni di libreria per ogni cella
#'
#' Simula le dimensioni delle librerie per ciascuna cella/spot,
#' includendo effetti spaziali e specifici del tipo cellulare.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param library_size_params Parametri per la simulazione della dimensione libreria
#' @param spatial_params Parametri spaziali
#' @param k_cell_types Numero di tipi cellulari
#' @param random_seed Seed per riproducibilità
#' @return Vettore con le dimensioni delle librerie per ogni cella
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_library_sizes <- function(
  cell_df,
  library_size_params = list(
    mean_library_size = 10000,
    library_size_cv = 0.3,
    spatial_effect_on_library = 0.5,
    cell_type_effect = TRUE
  ),
  spatial_params = list(
    spatial_range = 30
  ),
  k_cell_types,
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai variabili
  N <- nrow(cell_df)
  cluster_labels <- cell_df$intensity_cluster
  
  # Genera dimensioni libreria con distribuzione log-normale
  library_size_sd <- library_size_params$mean_library_size * library_size_params$library_size_cv
  log_mean <- log(library_size_params$mean_library_size^2 /
                 sqrt(library_size_params$mean_library_size^2 + library_size_sd^2))
  log_sd <- sqrt(log(1 + (library_size_sd^2 / library_size_params$mean_library_size^2)))
  
  library_size <- rlnorm(N, meanlog = log_mean, sdlog = log_sd)
  
  # Aggiungi effetto spaziale sulla dimensione libreria se richiesto
  if (!is.null(library_size_params$spatial_effect_on_library) && 
      library_size_params$spatial_effect_on_library > 0) {
    # Converti cell_df in oggetto spatial per il GP
    sp_df_lib <- cell_df
    coordinates(sp_df_lib) <- ~ x + y
    
    # Crea un GP per l'effetto spaziale sulla dimensione libreria
    lib_gp <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                   beta = 0, model = vgm(psill = 1.0,
                                        range = spatial_params$spatial_range * 1.5,
                                        model = "Exp"),
                   nmax = 15)
    
    # Genera il noise spaziale
    lib_noise <- predict(lib_gp, newdata = sp_df_lib, nsim = 1)$sim1
    
    # Normalizza e scala il noise
    lib_noise <- scale(lib_noise)
    lib_effect <- library_size_params$spatial_effect_on_library * lib_noise
    
    # Applica alla dimensione libreria (effetto moltiplicativo)
    library_size <- library_size * exp(lib_effect)
  }
  
  # Aggiungi effetto del tipo cellulare sulla dimensione della libreria
  if (!is.null(library_size_params$cell_type_effect) && library_size_params$cell_type_effect) {
    # Diversi tipi cellulari hanno diversi contenuti di RNA
    cell_type_effect <- numeric(N)
    
    # Crea effetti diversi per diversi tipi cellulari
    type_effects <- rnorm(k_cell_types, mean = 0, sd = 0.2)  # Effetti casuali per tipo
    
    # Assegna effetto in base al tipo cellulare
    for (k in 1:k_cell_types) {
      cell_type_effect[cluster_labels == k] <- type_effects[k]
    }
    
    # Applica effetto moltiplicativo
    library_size <- library_size * exp(cell_type_effect)
  }
  
  # Applica cap ai valori di library_size per evitare outlier estremi
  library_size <- pmax(library_size, 1000)
  library_size <- pmin(library_size, 30000)

  # Debug: stampa statistiche library size dopo il cap
  cat("DEBUG - Library size (capped) stats: mean =", round(mean(library_size)), 
      ", median =", round(median(library_size)), 
      ", range = [", round(min(library_size)), ",", round(max(library_size)), "]\n")
  print(summary(library_size))

  return(library_size)
}