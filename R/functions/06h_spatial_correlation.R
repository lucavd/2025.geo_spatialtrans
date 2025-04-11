#' Genera effetti di correlazione spaziale
#'
#' Crea pattern di correlazione spaziale usando processi gaussiani
#' o altri metodi.
#'
#' @param cell_df Dataframe delle celle con coordinate
#' @param spatial_params Parametri spaziali
#' @param correlation_method Metodo di correlazione ("grf" o "car")
#' @param random_seed Seed per riproducibilità
#' @return Vettore con il rumore spaziale generato
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_spatial_correlation <- function(
  cell_df,
  spatial_params = list(
    spatial_noise_intensity = 1.0,
    spatial_range = 30
  ),
  correlation_method = "grf",
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Inizializza il risultato
  gp_noise <- NULL
  
  # Procedi solo se abbiamo un metodo di correlazione
  if (correlation_method %in% c("grf", "car")) {
    # Converti cell_df in oggetto spatial
    sp_df <- cell_df
    coordinates(sp_df) <- ~ x + y
    
    if (correlation_method == "grf") {
      # Gaussian Random Field
      gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                     beta = 0, model = vgm(psill = spatial_params$spatial_noise_intensity * 2,
                                          range = spatial_params$spatial_range * 0.8,
                                          model = "Exp"),
                     nmax = 20)
      
      # Genera il noise spaziale
      gp_noise <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
      
      # Normalizza il noise
      gp_noise <- scale(gp_noise) * 1.5
    } else if (correlation_method == "car") {
      # Conditional Autoregressive Model - implementazione semplificata
      # In una versione completa, usare spdep o altre librerie per CAR
      
      # Per ora, usiamo una versione semplificata basata su vicinato
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
  
  return(gp_noise)
}