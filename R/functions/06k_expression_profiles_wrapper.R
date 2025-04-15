#' Genera profili di espressione genica
#'
#' Crea profili di espressione genica per ogni cellula/spot,
#' incorporando effetti biologici come correlazione spaziale,
#' dropout, dimensioni diverse delle librerie e clustering di geni.
#' Questa è una funzione wrapper che coordina tutti i passaggi.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param n_genes Numero di geni da simulare
#' @param k_cell_types Numero di tipi cellulari
#' @param marker_params Parametri dei marker genici
#' @param spatial_params Parametri spaziali
#' @param dropout_params Parametri di dropout
#' @param library_size_params Parametri dimensione libreria
#' @param cell_specific_params Parametri cellula-specifici
#' @param hybrid_params Parametri cellule ibride
#' @param use_spatial_correlation Se TRUE usa correlazione spaziale
#' @param correlation_method Metodo di correlazione ("grf" o "car")
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrice di espressione e metadati associati
#' @export
generate_expression_profiles <- function(
  cell_df,
  n_genes,
  k_cell_types,
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
  use_spatial_correlation = TRUE,
  correlation_method = "grf",
  random_seed = 123
) {
  # Inizializza e valida i parametri
  params <- initialize_expression_params(
    marker_params, spatial_params, dropout_params,
    library_size_params, cell_specific_params, hybrid_params,
    random_seed
  )
  
  # Genera i profili di espressione baseline per ogni tipo cellulare
  mean_expression_list <- generate_baseline_expression(
    n_genes, k_cell_types, params$marker_params, random_seed
  )
  
  # Calcola le distanze spaziali e le densità locali
  spatial_distances <- calculate_spatial_distances(
    cell_df, chunk_size = NULL, random_seed
  )
  
  # Calcola i parametri di dispersione
  dispersion_param <- calculate_dispersion_params(
    cell_df, spatial_distances$mean_dist, 
    params$dropout_params, k_cell_types, random_seed
  )
  
  # Genera le dimensioni delle librerie
  library_size <- generate_library_sizes(
    cell_df, params$library_size_params, 
    params$spatial_params, k_cell_types, random_seed
  )
  
  # Calcola le probabilità di dropout
  base_dropout <- calculate_dropout_probabilities(
    cell_df, spatial_distances$mean_dist, 
    params$dropout_params, params$spatial_params
  )
  
  # Genera i moduli genici
  gene_modules_result <- generate_gene_modules(
    n_genes, nrow(cell_df), params$cell_specific_params, random_seed
  )
  
  # Genera la correlazione spaziale
  gp_noise <- NULL
  if (use_spatial_correlation) {
    gp_noise <- generate_spatial_correlation(
      cell_df, params$spatial_params, correlation_method, random_seed
    )
  }
  
  # Genera le cellule ibride
  hybrid_matrix <- generate_hybrid_cells(
    cell_df, spatial_distances$dist_mat, 
    k_cell_types, params$hybrid_params, random_seed
  )
  
  # Genera la matrice di espressione finale
  expression_data <- generate_expression_matrix(
    cell_df, mean_expression_list, n_genes, library_size, dispersion_param,
    base_dropout, hybrid_matrix, gene_modules_result$module_noise, gp_noise,
    gene_modules_result$latent_factors, gene_modules_result$module_network,
    params$spatial_params, params$dropout_params, params$cell_specific_params,
    use_spatial_correlation, random_seed
  )
  
  # Prepara l'output
  result <- list(
    expression = expression_data,
    library_size = library_size,
    dispersion_param = dispersion_param,
    gene_modules = gene_modules_result$gene_modules,
    latent_factors = gene_modules_result$latent_factors,
    module_network = gene_modules_result$module_network,
    mean_expression_list = mean_expression_list
  )
  
  return(result)
}