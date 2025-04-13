#' Simula dati di trascrittomica spaziale
#' 
#' Funzione principale per la simulazione completa di dati di trascrittomica spaziale
#' parametrizzati per piattaforme come Visium HD, con controllo completo su
#' aspetti biologici e tecnici. Questa è una funzione wrapper che coordina tutti i componenti.
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
  # Inizializza la configurazione
  config <- initialize_simulation_config(
    image_path, output_path, output_plot, n_cells, n_genes, k_cell_types,
    threshold_value, random_seed, pixel_size_um, grid_mode, grid_resolution,
    grid_spacing, use_fixed_grid, fixed_grid_width_mm, fixed_grid_height_mm
  )
  
  # Configura il livello di difficoltà
  difficulty_config <- configure_difficulty_level(
    difficulty_level, marker_params, spatial_params, 
    dropout_params, cell_specific_params
  )
  
  # Esegui la pipeline di simulazione
  result <- run_simulation_pipeline(
    config, difficulty_config, use_spatial_correlation,
    correlation_method, library_size_params, hybrid_params
  )
  
  # Preleva i dati delle celle per la visualizzazione
  # (aggiungiamo questa estrazione poiché cell_df non è disponibile direttamente in questo scope)
  cell_df <- data.frame(
    x = result$coordinates$x,
    y = result$coordinates$y,
    intensity_cluster = result$intensity_cluster
  )
  
  # Genera e salva i plot
  generate_and_save_plots(
    cell_df, config, difficulty_config, output_plot
  )
  
  # Salva i risultati
  save_simulation_results(result, config$output_path)
  
  return(result)
}