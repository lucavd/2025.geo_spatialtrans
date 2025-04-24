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
  
  # 5. Applicazione di moduli biologici aggiuntivi
  
  # Inizializza una lista per metadati dei moduli
  module_metadata <- list()
  
  # Calcola la matrice di distanza una volta sola se necessaria
  dist_mat <- NULL
  if (config$lr_params$use_lr_interactions || 
      config$anisotropic_params$use_anisotropic_patterns ||
      config$micro3d_params$use_3d_microenvironment) {
    dist_mat <- as.matrix(dist(cell_df[, c("x", "y")]))
  }
  
  # 5.1. Ligand-Receptor Interactions
  if (config$lr_params$use_lr_interactions) {
    cat("Applicazione di interazioni ligando-recettore...\n")
    lr_result <- generate_lr_interactions(
      cell_df = cell_df,
      expr_matrix = expression_results$expression,
      dist_mat = dist_mat,
      lr_params = config$lr_params,
      random_seed = config$random_seed
    )
    # Aggiorna la matrice di espressione
    expression_results$expression <- lr_result$expression
    # Salva metadati
    module_metadata$lr_interactions <- list(
      interaction_db = lr_result$interaction_db,
      signaling_effects = lr_result$signaling_effects
    )
  }
  
  # 5.2. Temporal Dynamics
  if (config$temporal_params$use_temporal_dynamics) {
    cat("Applicazione di dinamiche temporali...\n")
    temporal_result <- generate_temporal_dynamics(
      cell_df = cell_df,
      expr_matrix = expression_results$expression,
      temporal_params = config$temporal_params,
      random_seed = config$random_seed
    )
    # Aggiorna la matrice di espressione
    expression_results$expression <- temporal_result$expr_matrix
    # Salva metadati
    module_metadata$temporal_dynamics <- list(
      pseudotime = temporal_result$pseudotime,
      velocity = temporal_result$velocity,
      unspliced = temporal_result$unspliced
    )
  }
  
  # 5.3. Alternative Splicing
  if (config$splicing_params$use_alternative_splicing) {
    cat("Applicazione di splicing alternativo...\n")
    splicing_result <- generate_alternative_splicing(
      cell_df = cell_df,
      expr_matrix = expression_results$expression,
      splicing_params = config$splicing_params,
      random_seed = config$random_seed
    )
    # Aggiorna la matrice di espressione
    expression_results$expression <- splicing_result$expr_matrix
    # Salva metadati
    module_metadata$alternative_splicing <- list(
      genes_with_variants = splicing_result$genes_with_variants,
      variant_matrices = splicing_result$variant_matrices
    )
  }
  
  # 5.4. Anisotropic Patterns
  if (config$anisotropic_params$use_anisotropic_patterns) {
    cat("Applicazione di pattern anisotropici...\n")
    aniso_result <- generate_anisotropic_patterns(
      cell_df = cell_df,
      expr_matrix = expression_results$expression,
      anisotropic_params = config$anisotropic_params,
      random_seed = config$random_seed
    )
    # Aggiorna la matrice di espressione
    expression_results$expression <- aniso_result$expr_matrix
    # Salva metadati
    module_metadata$anisotropic_patterns <- list(
      structure_mask = aniso_result$structure_mask,
      distance_matrices = aniso_result$distance_matrices
    )
  }
  
  # 5.5. 3D Microenvironment
  if (config$micro3d_params$use_3d_microenvironment) {
    cat("Applicazione di effetti microambiente 3D...\n")
    micro3d_result <- generate_3d_microenvironment(
      cell_df = cell_df,
      expr_matrix = expression_results$expression,
      dist_mat = dist_mat,
      micro3d_params = config$micro3d_params,
      random_seed = config$random_seed
    )
    # Aggiorna la matrice di espressione
    expression_results$expression <- micro3d_result$expr_matrix
    # Salva metadati
    module_metadata$microenvironment_3d <- list(
      z_positions = micro3d_result$z_positions,
      layer_assignments = micro3d_result$layer_assignments
    )
  }
  
  # 6. Prepara il risultato
  result <- list(
    coordinates = cell_df[, c("x", "y")],
    intensity_cluster = cell_df$intensity_cluster,
    expression = expression_results$expression,
    threshold_used = config$threshold_value,
    library_size = expression_results$library_size,
    dispersion_param = expression_results$dispersion_param,
    gene_modules = expression_results$gene_modules,
    
    # Aggiungi i metadati dei moduli biologici aggiuntivi
    module_data = module_metadata,
    
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
      cell_specific_params = difficulty_config$cell_specific_params,
      
      # Parametri dei nuovi moduli biologici
      lr_params = config$lr_params,
      temporal_params = config$temporal_params,
      splicing_params = config$splicing_params,
      anisotropic_params = config$anisotropic_params,
      micro3d_params = config$micro3d_params
    )
  )
  
  # Registra tempo di esecuzione
  execution_time <- tictoc::toc(quiet = TRUE)
  result$execution_time <- execution_time$toc - execution_time$tic
  
  return(result)
}