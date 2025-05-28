#' Configura i parametri in base al livello di difficoltà
#'
#' Imposta i parametri della simulazione in base al livello di difficoltà scelto.
#'
#' @param difficulty_level Livello di difficoltà ("easy", "medium", "hard")
#' @param marker_params Lista di parametri per i geni marker
#' @param spatial_params Lista di parametri spaziali
#' @param dropout_params Lista di parametri di dropout
#' @param cell_specific_params Lista di parametri cellula-specifici
#' @return Lista con parametri aggiornati in base al livello di difficoltà
#' @export
configure_difficulty_level <- function(
  difficulty_level = "hard",
  marker_params = list(
    marker_genes_per_type = NULL,
    marker_expression_fold = NULL,
    marker_overlap_fold = NULL
  ),
  spatial_params = list(
    spatial_noise_intensity = NULL,
    spatial_range = NULL,
    random_noise_sd = NULL,
    gradient_regions = FALSE,
    gradient_width = 5,
    gradient_exponent = 1.5
  ),
  dropout_params = list(
    dropout_range = NULL,
    dispersion_range = NULL,
    cell_type_dispersion_effect = 0.2,
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5
  ),
  cell_specific_params = list(
    cell_specific_noise_sd = NULL,
    use_gene_modules = TRUE,
    n_gene_modules = 5,
    module_correlation = 0.7
  )
) {
  # Valida il livello di difficoltà
  if (!difficulty_level %in% c("easy", "medium", "hard")) {
    warning("Livello di difficoltà non valido, impostato su 'medium'")
    difficulty_level <- "medium"
  }
  
  # Imposta i parametri in base al livello di difficoltà se non specificati direttamente
  # Per marker_params
  if (is.null(marker_params$marker_genes_per_type)) {
    marker_params$marker_genes_per_type <- switch(
      difficulty_level,
      "easy" = 10,
      "medium" = 7,
      "hard" = 5
    )
  }
  
  if (is.null(marker_params$marker_expression_fold)) {
    marker_params$marker_expression_fold <- switch(
      difficulty_level,
      "easy" = 2.0,
      "medium" = 1.2,
      "hard" = 0.8
    )
  }
  
  if (is.null(marker_params$marker_overlap_fold)) {
    marker_params$marker_overlap_fold <- switch(
      difficulty_level,
      "easy" = 0.0,
      "medium" = 0.2,
      "hard" = 0.4
    )
  }
  
  # Per spatial_params
  if (is.null(spatial_params$spatial_noise_intensity)) {
    spatial_params$spatial_noise_intensity <- switch(
      difficulty_level,
      "easy" = 0.5,
      "medium" = 1.0,
      "hard" = 1.5
    )
  }
  
  if (is.null(spatial_params$spatial_range)) {
    spatial_params$spatial_range <- switch(
      difficulty_level,
      "easy" = 50,
      "medium" = 30,
      "hard" = 15
    )
  }
  
  if (is.null(spatial_params$random_noise_sd)) {
    spatial_params$random_noise_sd <- switch(
      difficulty_level,
      "easy" = 0.1,
      "medium" = 0.2,
      "hard" = 0.4
    )
  }
  
  # Per dropout_params
  if (is.null(dropout_params$dropout_range)) {
    dropout_params$dropout_range <- switch(
      difficulty_level,
      "easy" = c(0.5, 0.7),
      "medium" = c(0.4, 0.6),
      "hard" = c(0.7, 0.9)
    )
  }
  
  if (is.null(dropout_params$dispersion_range)) {
    dropout_params$dispersion_range <- switch(
      difficulty_level,
      "easy" = c(15.0, 10.0),
      "medium" = c(10.0, 5.0),
      "hard" = c(5.0, 2.0)
    )
  }
  
  if (is.null(cell_specific_params$cell_specific_noise_sd)) {
    cell_specific_params$cell_specific_noise_sd <- switch(
      difficulty_level,
      "easy" = 0.1,
      "medium" = 0.2,
      "hard" = 0.3
    )
  }
  
  # Restituisci i parametri aggiornati
  return(list(
    difficulty_level = difficulty_level,
    marker_params = marker_params,
    spatial_params = spatial_params,
    dropout_params = dropout_params,
    cell_specific_params = cell_specific_params
  ))
}