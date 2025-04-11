#' Esegue la pipeline di simulazione di dati di trascrittomica spaziale
#'
#' Coordina i vari step della simulazione, dall'elaborazione dell'immagine
#' alla generazione dei profili di espressione.
#'
#' @param config Configurazione di base della simulazione
#' @param difficulty_config Configurazione di difficoltà
#' @param use_spatial_correlation Se usare correlazione spaziale
#' @param correlation_method Metodo di correlazione ("grf" o "car")
#' @param library_size_params Parametri dimensione libreria
#' @param hybrid_params Parametri cellule ibride
#' @return Lista con risultati della simulazione: celle, espressione, metadati
#' @importFrom tictoc tic toc
#' @export
run_simulation_pipeline <- function(
  config,
  difficulty_config,
  use_spatial_correlation = TRUE,
  correlation_method = "grf",
  library_size_params = list(
    mean_library_size = 10000,
    library_size_cv = 0.3,
    spatial_effect_on_library = 0.5,
    cell_type_effect = TRUE
  ),
  hybrid_params = list(
    use_hybrid_cells = TRUE,
    max_hybrid_pairs = 1000,
    hybrid_intensity_range = c(0.2, 0.5)
  )
) {
  # Avvia timer
  tictoc::tic()
  
  # Imposta il seed
  set.seed(config$random_seed)

  # 1. Preparazione immagine
  image_data <- prepare_image(config$image_path, config$threshold_value)
  
  # 2. Clustering dell'immagine
  clustered_data <- cluster_image(
    img_df_thresh = image_data$img_df_thresh,
    k_cell_types = config$k_cell_types,
    random_seed = config$random_seed
  )
  
  # 3. Creazione griglia di campionamento
  cell_df <- create_sampling_grid(
    img_df_thresh = clustered_data,
    img_array = image_data$img_array,
    img_width = image_data$width,
    img_height = image_data$height,
    grid_mode = config$grid_mode,
    n_cells = config$n_cells,
    k_cell_types = config$k_cell_types,
    grid_resolution = config$grid_resolution,
    grid_spacing = config$grid_spacing,
    use_fixed_grid = config$use_fixed_grid,
    fixed_grid_width_mm = config$fixed_grid_width_mm,
    fixed_grid_height_mm = config$fixed_grid_height_mm,
    pixel_size_um = config$pixel_size_um,
    threshold_value = config$threshold_value,
    random_seed = config$random_seed
  )
  
  # 4. Generazione profili di espressione
  expression_results <- generate_expression_profiles(
    cell_df = cell_df,
    n_genes = config$n_genes,
    k_cell_types = config$k_cell_types,
    marker_params = difficulty_config$marker_params,
    spatial_params = difficulty_config$spatial_params,
    dropout_params = difficulty_config$dropout_params,
    library_size_params = library_size_params,
    hybrid_params = hybrid_params,
    cell_specific_params = difficulty_config$cell_specific_params,
    use_spatial_correlation = use_spatial_correlation,
    correlation_method = correlation_method,
    random_seed = config$random_seed
  )
  
  # Rinomina i cluster
  levels(cell_df$intensity_cluster) <- paste0("cells_", letters[1:config$k_cell_types])
  
  # 5. Prepara il risultato
  result <- list(
    coordinates = cell_df[, c("x", "y")],
    intensity_cluster = cell_df$intensity_cluster,
    expression = expression_results$expression,
    threshold_used = config$threshold_value,
    library_size = expression_results$library_size,
    dispersion_param = expression_results$dispersion_param,
    gene_modules = expression_results$gene_modules,
    parameters = list(
      image_path = config$image_path,
      pixel_size_um = config$pixel_size_um,
      n_cells = config$n_cells,
      n_genes = config$n_genes,
      k_cell_types = config$k_cell_types,
      difficulty_level = difficulty_config$difficulty_level,
      grid_mode = config$grid_mode,
      grid_resolution = config$grid_resolution,
      use_fixed_grid = config$use_fixed_grid,
      fixed_grid_width_mm = config$fixed_grid_width_mm,
      fixed_grid_height_mm = config$fixed_grid_height_mm,
      use_spatial_correlation = use_spatial_correlation,
      correlation_method = correlation_method,
      marker_params = difficulty_config$marker_params,
      spatial_params = difficulty_config$spatial_params,
      dropout_params = difficulty_config$dropout_params,
      library_size_params = library_size_params,
      hybrid_params = hybrid_params,
      cell_specific_params = difficulty_config$cell_specific_params
    )
  )
  
  # Registra tempo di esecuzione
  execution_time <- tictoc::toc(quiet = TRUE)
  result$execution_time <- execution_time$toc - execution_time$tic
  
  return(result)
}