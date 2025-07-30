#' Simula dati di trascrittomica spaziale
#' 
#' Funzione principale per la simulazione completa di dati di trascrittomica spaziale
#' parametrizzati per piattaforme come Visium HD, con controllo completo su
#' aspetti biologici e tecnici.
#'
#' @param image_path Path dell'immagine
#' @param output_path Path per salvare i risultati (default: "data/simulated_image_correlation.rds")
#' @param output_plot Path per salvare il plot (opzionale, default: NULL)
#' @param n_cells Numero di celle da campionare (se grid_mode=FALSE) (default: 20000)
#' @param n_genes Numero di geni totali (default: 100)
#' @param k_cell_types Numero di tipi cellulari (default: 5)
#' @param threshold_value Soglia per il thresholding dell'immagine (default: 0.7)
#' @param random_seed Seed per riproducibilità (default: 123)
#' @param pixel_size_um Dimensione di ogni pixel in um (default: 1)
#' @param grid_mode Usare modalità griglia (Visium HD) invece di sampling casuale (default: TRUE)
#' @param grid_resolution Dimensione della griglia in um (default: 2)
#' @param grid_spacing Spazio tra bin della griglia (default: 0)
#' @param use_fixed_grid Usare una griglia fissa con dimensioni predefinite (default: FALSE)
#' @param fixed_grid_width_mm Larghezza della griglia fissa in mm (default: 6.5)
#' @param fixed_grid_height_mm Altezza della griglia fissa in mm (default: 6.5)
#' @param difficulty_level Livello di difficoltà ("easy", "medium", "hard") (default: "hard")
#' @param use_spatial_correlation Usare correlazione spaziale (default: TRUE)
#' @param correlation_method Metodo di correlazione spaziale ("grf" o "car") (default: "grf")
#' @param marker_params Lista di parametri per i geni marker
#' @param spatial_params Lista di parametri spaziali
#' @param dropout_params Lista di parametri di dropout
#' @param library_size_params Lista di parametri per dimensione libreria
#' @param hybrid_params Lista di parametri per cellule ibride
#' @param cell_specific_params Lista di parametri cellula-specifici
#' @return Lista con i risultati della simulazione
#' @importFrom tictoc tic toc
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_y_reverse coord_fixed theme_minimal labs
#' @export
simulate_spatial_transcriptomics <- function(
  # Parametri generali
  image_path = NULL,
  output_path = "data/simulated_image_correlation.rds",
  output_plot = NULL,
  n_cells = 20000,
  n_genes = 100,
  k_cell_types = 5,
  threshold_value = 0.7,
  random_seed = 123,
  pixel_size_um = 1,
  
  # Parametri della griglia Visium HD
  grid_mode = TRUE,
  grid_resolution = 2,
  grid_spacing = 0,
  use_fixed_grid = FALSE,
  fixed_grid_width_mm = 6.5,
  fixed_grid_height_mm = 6.5,
  
  # Parametri per la simulazione dell'espressione genica
  difficulty_level = "hard",
  use_spatial_correlation = TRUE,
  correlation_method = "grf",
  
  # Parametri avanzati (personalizzabili)
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
  ),
  
  cell_specific_params = list(
    cell_specific_noise_sd = NULL,
    use_gene_modules = TRUE,
    n_gene_modules = 5,
    module_correlation = 0.7
  )
) {
  tictoc::tic()
  set.seed(random_seed)
  
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
      "easy" = c(0.1, 0.3),
      "medium" = c(0.2, 0.5),
      "hard" = c(0.3, 0.7)
    )
  }
  
  if (is.null(dropout_params$dispersion_range)) {
    dropout_params$dispersion_range <- switch(
      difficulty_level,
      "easy" = c(3.0, 1.5),
      "medium" = c(2.0, 1.0),
      "hard" = c(1.5, 0.8)
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
  
  # 1. Preparazione immagine
  image_data <- prepare_image(image_path, threshold_value)
  
  # 2. Clustering dell'immagine
  clustered_data <- cluster_image(
    img_df_thresh = image_data$img_df_thresh,
    k_cell_types = k_cell_types,
    random_seed = random_seed
  )
  
  # 3. Creazione griglia di campionamento
  cell_df <- create_sampling_grid(
    img_df_thresh = clustered_data,
    img_array = image_data$img_array,
    img_width = image_data$width,
    img_height = image_data$height,
    grid_mode = grid_mode,
    n_cells = n_cells,
    k_cell_types = k_cell_types,
    grid_resolution = grid_resolution,
    grid_spacing = grid_spacing,
    use_fixed_grid = use_fixed_grid,
    fixed_grid_width_mm = fixed_grid_width_mm,
    fixed_grid_height_mm = fixed_grid_height_mm,
    pixel_size_um = pixel_size_um,
    threshold_value = threshold_value,
    random_seed = random_seed
  )
  
  # 4. Generazione profili di espressione
  expression_results <- generate_expression_profiles(
    cell_df = cell_df,
    n_genes = n_genes,
    k_cell_types = k_cell_types,
    marker_params = marker_params,
    spatial_params = spatial_params,
    dropout_params = dropout_params,
    library_size_params = library_size_params,
    hybrid_params = hybrid_params,
    cell_specific_params = cell_specific_params,
    use_spatial_correlation = use_spatial_correlation,
    correlation_method = correlation_method,
    random_seed = random_seed
  )
  
  # 5. Visualizzazione dei risultati
  if (grid_mode) {
    # Per griglia, usiamo geom_tile
    p <- ggplot(cell_df, aes(x = x, y = y, fill = intensity_cluster)) +
      geom_tile(width = grid_resolution, height = grid_resolution) +
      scale_y_reverse() +
      coord_fixed() +
      theme_minimal() +
      labs(title = sprintf("Visium HD (2μm) - Livello difficoltà: %s", difficulty_level),
           subtitle = sprintf("Griglia %dμm, %d bin, %d geni (marker/tipo: %d, fold: %.1f)",
                             grid_resolution, nrow(cell_df), n_genes,
                             marker_params$marker_genes_per_type,
                             marker_params$marker_expression_fold),
           fill = "Cell Type")
  } else {
    # Per sampling casuale, usiamo punti
    p <- ggplot(cell_df, aes(x = x, y = y, color = intensity_cluster)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_y_reverse() +
      coord_fixed() +
      theme_minimal() +
      labs(title = sprintf("Distribuzione spaziale (livello difficoltà: %s)", difficulty_level),
           subtitle = sprintf("Threshold: %.2f, Geni marker: %d per tipo, Fold-change: %.1f",
                             threshold_value,
                             marker_params$marker_genes_per_type,
                             marker_params$marker_expression_fold),
           color = "Cell Type")
  }
  
  print(p)
  
  # Salva i plot se richiesto
  if (!is.null(output_plot)) {
    # Salva il plot principale
    ggplot2::ggsave(output_plot, plot = p, device = "png", dpi = 300)
  }
  
  # Rinomino i cluster
  levels(cell_df$intensity_cluster) <- paste0("cells_", letters[1:k_cell_types])
  
  # 6. Prepara e salva il risultato
  result <- list(
    coordinates = cell_df[, c("x", "y")],
    intensity_cluster = cell_df$intensity_cluster,
    expression = expression_results$expression,
    threshold_used = threshold_value,
    library_size = expression_results$library_size,
    dispersion_param = expression_results$dispersion_param,
    gene_modules = expression_results$gene_modules,
    parameters = list(
      image_path = image_path,
      pixel_size_um = pixel_size_um,
      n_cells = n_cells,
      n_genes = n_genes,
      k_cell_types = k_cell_types,
      difficulty_level = difficulty_level,
      grid_mode = grid_mode,
      grid_resolution = grid_resolution,
      use_fixed_grid = use_fixed_grid,
      fixed_grid_width_mm = fixed_grid_width_mm,
      fixed_grid_height_mm = fixed_grid_height_mm,
      use_spatial_correlation = use_spatial_correlation,
      correlation_method = correlation_method,
      marker_params = marker_params,
      spatial_params = spatial_params,
      dropout_params = dropout_params,
      library_size_params = library_size_params,
      hybrid_params = hybrid_params,
      cell_specific_params = cell_specific_params
    )
  )
  
  # Crea directory se non esiste
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  
  # Salva il risultato
  saveRDS(result, file = output_path)
  
  tictoc::toc()
  
  return(result)
}