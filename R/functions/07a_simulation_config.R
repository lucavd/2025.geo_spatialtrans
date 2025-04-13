#' Inizializza la configurazione della simulazione
#'
#' Configura i parametri di base per la simulazione di dati di trascrittomica spaziale.
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
#' @return Lista con la configurazione della simulazione
#' @export
initialize_simulation_config <- function(
  image_path = NULL,
  output_path = "data/simulated_image_correlation.rds",
  output_plot = NULL,
  n_cells = 20000,
  n_genes = 100,
  k_cell_types = 5,
  threshold_value = 0.7,
  random_seed = 123,
  pixel_size_um = 1,
  grid_mode = TRUE,
  grid_resolution = 2,
  grid_spacing = 0,
  use_fixed_grid = FALSE,
  fixed_grid_width_mm = 6.5,
  fixed_grid_height_mm = 6.5
) {
  # Validazione di base
  if (is.null(image_path) || !file.exists(image_path)) {
    stop("È necessario specificare un percorso di immagine valido")
  }
  
  # Verifica parametri numerici
  if (n_cells <= 0) stop("n_cells deve essere positivo")
  if (n_genes <= 0) stop("n_genes deve essere positivo")
  if (k_cell_types <= 0) stop("k_cell_types deve essere positivo")
  if (threshold_value <= 0 || threshold_value >= 1) {
    warning("threshold_value dovrebbe essere tra 0 e 1, impostato a 0.7")
    threshold_value <- 0.7
  }
  
  # Controlla directory di output
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Restituisce la configurazione
  config <- list(
    image_path = image_path,
    output_path = output_path,
    output_plot = output_plot,
    n_cells = n_cells,
    n_genes = n_genes,
    k_cell_types = k_cell_types,
    threshold_value = threshold_value,
    random_seed = random_seed,
    pixel_size_um = pixel_size_um,
    grid_mode = grid_mode,
    grid_resolution = grid_resolution,
    grid_spacing = grid_spacing,
    use_fixed_grid = use_fixed_grid,
    fixed_grid_width_mm = fixed_grid_width_mm,
    fixed_grid_height_mm = fixed_grid_height_mm
  )
  
  return(config)
}