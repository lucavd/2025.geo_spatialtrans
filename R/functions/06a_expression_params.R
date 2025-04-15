#' Inizializza parametri per la generazione di profili di espressione
#'
#' Questa funzione inizializza e valida i parametri necessari per la generazione
#' di profili di espressione genica.
#'
#' @param marker_params Parametri relativi ai geni marker
#' @param spatial_params Parametri spaziali
#' @param dropout_params Parametri di dropout
#' @param library_size_params Parametri dimensione libreria
#' @param cell_specific_params Parametri cellula-specifici
#' @param hybrid_params Parametri per cellule ibride
#' @param random_seed Seed per riproducibilità
#' @return Lista con parametri validati
#' @export
initialize_expression_params <- function(
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  ),
  spatial_params = list(
    spatial_noise_intensity = 1.0,
    spatial_range = 30,
    random_noise_sd = 0.2,
    gradient_regions = FALSE,
    gradient_width = 5,
    gradient_exponent = 1.5
  ),
  dropout_params = list(
    dropout_range = c(0.2, 0.5),
    dispersion_range = c(2.0, 1.0),
    cell_type_dispersion_effect = 0.2,
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5
  ),
  library_size_params = list(
    mean_library_size = 10000,
    library_size_cv = 0.3,
    spatial_effect_on_library = 0.5,
    cell_type_effect = TRUE
  ),
  cell_specific_params = list(
    cell_specific_noise_sd = 0.2,
    use_gene_modules = TRUE,
    n_gene_modules = 5,
    module_correlation = 0.7,
    module_hierarchical = FALSE,
    module_overlap = 0.1,
    module_size_distribution = "exponential",
    n_latent_factors = 3,
    module_network_density = 0.2,
    latent_factor_strength = 0.8
  ),
  hybrid_params = list(
    use_hybrid_cells = TRUE,
    max_hybrid_pairs = 1000,
    hybrid_intensity_range = c(0.2, 0.5)
  ),
  random_seed = 123
) {
  # Validazione parametri
  
  # Validazione marker_params
  if (is.null(marker_params$marker_genes_per_type) || 
      !is.numeric(marker_params$marker_genes_per_type) || 
      marker_params$marker_genes_per_type < 0) {
    warning("marker_genes_per_type deve essere un numero non negativo, impostato a 10")
    marker_params$marker_genes_per_type <- 10
  }
  
  if (is.null(marker_params$marker_expression_fold) || 
      !is.numeric(marker_params$marker_expression_fold)) {
    warning("marker_expression_fold deve essere numerico, impostato a 1.5")
    marker_params$marker_expression_fold <- 1.5
  }
  
  if (is.null(marker_params$marker_overlap_fold) || 
      !is.numeric(marker_params$marker_overlap_fold) || 
      marker_params$marker_overlap_fold < 0) {
    warning("marker_overlap_fold deve essere un numero non negativo, impostato a 0.2")
    marker_params$marker_overlap_fold <- 0.2
  }
  
  # Validazione dropout_params
  if (!is.numeric(dropout_params$dropout_range) || length(dropout_params$dropout_range) != 2) {
    warning("dropout_range deve essere un vettore numerico di lunghezza 2, impostato a c(0.2, 0.5)")
    dropout_params$dropout_range <- c(0.2, 0.5)
  }
  
  if (!is.numeric(dropout_params$dispersion_range) || length(dropout_params$dispersion_range) != 2) {
    warning("dispersion_range deve essere un vettore numerico di lunghezza 2, impostato a c(2.0, 1.0)")
    dropout_params$dispersion_range <- c(2.0, 1.0)
  }
  
  # Imposta il seed
  set.seed(random_seed)
  
  # Restituisci tutti i parametri validati
  return(list(
    marker_params = marker_params,
    spatial_params = spatial_params,
    dropout_params = dropout_params,
    library_size_params = library_size_params,
    cell_specific_params = cell_specific_params,
    hybrid_params = hybrid_params,
    random_seed = random_seed
  ))
}