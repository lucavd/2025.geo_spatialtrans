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
#' @param use_ambient_rna Abilitare la contaminazione da RNA ambientale (default: FALSE)
#' @param ambient_contamination_rate Tasso di contaminazione RNA ambientale (default: 0.05)
#' @param use_gene_specific_dropout Abilitare il dropout gene-specifico (default: TRUE)
#' @param gene_dropout_variability Variabilità del dropout gene-specifico (default: 0.3)
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
  fixed_grid_height_mm = 6.5,
  use_ambient_rna = FALSE,
  ambient_contamination_rate = 0.05,
  use_gene_specific_dropout = TRUE,
  gene_dropout_variability = 0.3
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
  
  # Verifica parametri di ambient RNA e dropout
  if (ambient_contamination_rate < 0 || ambient_contamination_rate > 1) {
    warning("ambient_contamination_rate dovrebbe essere tra 0 e 1, impostato a 0.05")
    ambient_contamination_rate <- 0.05
  }
  
  if (gene_dropout_variability < 0 || gene_dropout_variability > 1) {
    warning("gene_dropout_variability dovrebbe essere tra 0 e 1, impostato a 0.3")
    gene_dropout_variability <- 0.3
  }
  
  # Controlla directory di output
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Crea parametri per ambient RNA
  ambient_params <- list(
    use_ambient_rna = use_ambient_rna,
    ambient_contamination_rate = ambient_contamination_rate,
    ambient_diffusion_distance = 30,
    tissue_leakage_factor = 0.7,
    background_noise = 0.1
  )
  
  # Crea parametri per dropout gene-specifico
  dropout_gene_specific_params <- list(
    use_gene_specific_dropout = use_gene_specific_dropout,
    gene_dropout_variability = gene_dropout_variability,
    gc_content_effect = 0.5,
    length_effect = 0.3,
    sequence_effect = 0.4,
    gene_effect_weight = 0.3
  )
  
  # Crea parametri per i nuovi moduli biologici
  
  # 1. Ligand-Receptor Interactions
  lr_params <- list(
    use_lr_interactions = FALSE,   # Disabilitato di default
    n_interactions = 20,
    signal_propagation_mode = "exponential",
    max_signaling_distance = 40,
    adjust_method = "multiplicative",
    signal_amplification = 1.0
  )
  
  # 2. Temporal Dynamics
  temporal_params <- list(
    use_temporal_dynamics = FALSE,  # Disabilitato di default
    pseudotime_mode = "gradient",
    pseudotime_origin = c(0, 0),
    temporal_gene_fraction = 0.6,
    pattern_distribution = c(monotonic = 0.4, transient = 0.3, cyclic = 0.2, bifurcating = 0.1),
    trajectory_strength = 0.8,
    include_velocity = TRUE
  )
  
  # 3. Alternative Splicing
  splicing_params <- list(
    use_alternative_splicing = FALSE,  # Disabilitato di default
    splicing_fraction = 0.3,
    n_splicing_variants = 2,
    splicing_spatial_pattern = "gradient",
    splicing_cluster_specific = FALSE,
    splicing_strength = 0.7
  )
  
  # 4. Anisotropic Patterns
  anisotropic_params <- list(
    use_anisotropic_patterns = FALSE,  # Disabilitato di default
    n_structures = 2,
    structure_type = "linear",
    anisotropic_pattern = "gradient",
    anisotropic_gene_fraction = 0.5,
    anisotropic_effect_strength = 0.8
  )
  
  # 5. 3D Microenvironment
  micro3d_params <- list(
    use_3d_microenvironment = FALSE,  # Disabilitato di default
    n_layers = 5,
    layer_specificity = 0.7,
    projection_noise = 0.2,
    z_decay_factor = 0.5
  )
  
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
    fixed_grid_height_mm = fixed_grid_height_mm,
    ambient_params = ambient_params,
    dropout_gene_specific_params = dropout_gene_specific_params,
    
    # Aggiungi nuovi moduli biologici
    lr_params = lr_params,
    temporal_params = temporal_params,
    splicing_params = splicing_params,
    anisotropic_params = anisotropic_params,
    micro3d_params = micro3d_params
  )
  
  return(config)
}