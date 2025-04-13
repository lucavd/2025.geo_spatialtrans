#' Calcola probabilità di dropout per ogni cella
#'
#' Genera le probabilità di dropout in base alla distanza dai confini
#' o alla distanza media tra celle dello stesso tipo.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param mean_dist Vettore con le distanze medie per ogni cella
#' @param dropout_params Parametri di dropout
#' @param spatial_params Parametri spaziali
#' @importFrom scales rescale
#' @return Vettore con le probabilità base di dropout per ogni cella
#' @export
calculate_dropout_probabilities <- function(
  cell_df,
  mean_dist,
  dropout_params = list(
    dropout_range = c(0.2, 0.5)
  ),
  spatial_params = list(
    gradient_regions = FALSE
  )
) {
  # Calcola le probabilità di dropout base
  if (spatial_params$gradient_regions && "boundary_dist" %in% colnames(cell_df) && !all(is.na(cell_df$boundary_dist))) {
    # Più dropout vicino al confine
    base_dropout <- dropout_params$dropout_range[1] +
      (1 - cell_df$boundary_dist) * (dropout_params$dropout_range[2] - dropout_params$dropout_range[1])
  } else {
    # Metodo originale basato sulla distanza media
    base_dropout <- scales::rescale(mean_dist, to = dropout_params$dropout_range)
  }
  
  return(base_dropout)
}

#' Applica dropout all'espressione genica
#'
#' Applica dropout all'espressione genica in modo espressione-dipendente e/o spaziale.
#'
#' @param expression Matrice o vettore di espressione
#' @param base_dropout Probabilità base di dropout
#' @param dropout_params Parametri di dropout
#' @param random_seed Seed per riproducibilità
#' @return Matrice o vettore di espressione con dropout applicato
#' @export
apply_dropout <- function(
  expression,
  base_dropout,
  dropout_params = list(
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Determina se l'input è una matrice o un vettore
  is_matrix <- is.matrix(expression)
  
  # Se è un vettore, trattalo come una matrice con una colonna
  if (!is_matrix) {
    expression <- matrix(expression, ncol = 1)
  }
  
  # Definisci una funzione vettorizzata per normalizzare tra 0 e 1
  scale01_vec <- function(x) {
    if (all(x == x[1])) return(rep(0.5, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }
  
  # Numero di celle
  N <- nrow(expression)
  
  # Applica dropout a ogni colonna (gene)
  for (g in 1:ncol(expression)) {
    if (dropout_params$expression_dependent_dropout) {
      # Normalizza l'espressione del gene corrente
      norm_expr <- scale01_vec(expression[, g])
      
      # Calcola la probabilità di dropout con una funzione logistica
      dropout_prob_expr <- 1 / (1 + exp((norm_expr - dropout_params$dropout_curve_midpoint) *
                                     dropout_params$dropout_curve_steepness))
      
      # Combina con il dropout spaziale base
      dropout_prob <- 0.7 * dropout_prob_expr + 0.3 * base_dropout
      
      # Tronca i valori al range [0,1]
      dropout_prob <- pmin(pmax(dropout_prob, 0), 1)
    } else {
      # Usa solo il modello di dropout spaziale
      dropout_prob <- base_dropout
    }
    
    # Applica dropout
    zero_idx <- runif(N) < dropout_prob
    expression[zero_idx, g] <- 0
  }
  
  # Restituisci nella forma originale
  if (!is_matrix) {
    expression <- as.vector(expression)
  }
  
  return(expression)
}