#' Normalizza un vettore tra 0 e 1
#'
#' @param x Vettore numerico da normalizzare
#' @return Vettore normalizzato tra 0 e 1
#' @examples
#' scale01(c(1, 2, 3, 4, 5))
#' @export
scale01 <- function(x) {
  if (max(x) == min(x)) return(rep(0.5, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

#' Normalizza un vettore tra 0 e 1 (versione vettorizzata)
#'
#' @param x Vettore numerico da normalizzare
#' @return Vettore normalizzato tra 0 e 1
#' @noRd
scale01_vec <- function(x) {
  if (all(x == x[1])) return(rep(0.5, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

#' Calcola Moran's I per dati point-based (spatial transcriptomics)
#'
#' @param spatial_data Data frame con colonne x, y e variable per analisi
#' @param variable Nome della variabile da analizzare (default: "intensity")
#' @param max_points Numero massimo di punti per efficiency (default: 5000)
#' @param distance_percentile Percentile distanza per soglia vicinanza (default: 0.1)
#' @return Lista con Moran's I, p-value (se calcolabile), e dettagli
#' @examples
#' df <- data.frame(x = runif(100), y = runif(100), intensity = runif(100))
#' moran_result <- compute_morans_i(df, "intensity")
#' @export
compute_morans_i <- function(spatial_data, 
                            variable = "intensity", 
                            max_points = 5000,
                            distance_percentile = 0.1) {
  
  # Controlli input
  if (!all(c("x", "y") %in% colnames(spatial_data))) {
    stop("spatial_data deve contenere colonne 'x' e 'y'")
  }
  
  if (!variable %in% colnames(spatial_data)) {
    stop(paste("Variabile", variable, "non trovata in spatial_data"))
  }
  
  # Rimuovi NA e campiona se necessario
  data_clean <- spatial_data[complete.cases(spatial_data[, c("x", "y", variable)]), ]
  
  if (nrow(data_clean) == 0) {
    return(list(moran_i = NA, interpretation = "No valid data", n_points = 0))
  }
  
  # Campionamento per efficiency su dataset grandi
  if (nrow(data_clean) > max_points) {
    sample_idx <- sample(nrow(data_clean), max_points)
    data_clean <- data_clean[sample_idx, ]
    sampled <- TRUE
  } else {
    sampled <- FALSE
  }
  
  n <- nrow(data_clean)
  if (n < 4) {
    return(list(moran_i = NA, interpretation = "Insufficient data points", n_points = n))
  }
  
  # Estrai coordinate e variabile
  coords <- as.matrix(data_clean[, c("x", "y")])
  values <- data_clean[[variable]]
  
  # Normalizza valori per Moran's I
  values_norm <- scale(values)[, 1]
  
  # Calcola matrice distanze
  dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
  
  # Soglia distanza basata su percentile per definire vicinanza
  # Usa distanze dal centroide per soglia più robusta
  centroid <- c(mean(coords[, 1]), mean(coords[, 2]))
  distances_from_center <- sqrt((coords[, 1] - centroid[1])^2 + 
                               (coords[, 2] - centroid[2])^2)
  distance_threshold <- quantile(distances_from_center, distance_percentile)
  
  # Crea matrice pesi spaziali
  # W[i,j] = 1 se distance <= threshold, 0 altrimenti  
  # Diagonale = 0 (un punto non è vicino a se stesso)
  W <- (dist_matrix <= distance_threshold) * 1.0
  diag(W) <- 0
  
  # Calcola totale pesi
  W_sum <- sum(W)
  if (W_sum == 0) {
    return(list(moran_i = NA, 
               interpretation = "No spatial neighbors found", 
               n_points = n, 
               distance_threshold = distance_threshold))
  }
  
  # Calcola Moran's I
  # Formula: I = (n/W) * (Σ w_ij * (x_i - x̄) * (x_j - x̄)) / Σ(x_i - x̄)²
  
  # Numeratore: somma prodotti pesati delle deviazioni
  numerator <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (W[i, j] > 0) {
        numerator <- numerator + W[i, j] * values_norm[i] * values_norm[j]
      }
    }
  }
  
  # Denominatore: somma quadrati deviazioni  
  denominator <- sum(values_norm^2)
  
  # Moran's I finale
  moran_i <- (n / W_sum) * (numerator / denominator)
  
  # Interpretazione biologica
  interpretation <- if (moran_i > 0.4) {
    "Strong positive spatial autocorrelation (excellent biological realism)"
  } else if (moran_i > 0.2) {
    "Moderate positive spatial autocorrelation (good spatial patterns)"
  } else if (moran_i > 0.05) {
    "Weak positive spatial autocorrelation (marginal spatial structure)"
  } else if (moran_i > -0.05) {
    "No significant spatial autocorrelation (random pattern)"
  } else {
    "Negative spatial autocorrelation (dispersed pattern)"
  }
  
  # Risultato completo
  result <- list(
    moran_i = moran_i,
    interpretation = interpretation,
    n_points = n,
    n_neighbors = W_sum / 2,  # Diviso 2 perché matrice simmetrica
    distance_threshold = distance_threshold,
    sampled = sampled,
    spatial_randomness = abs(moran_i) < 0.05
  )
  
  return(result)
}