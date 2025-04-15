# Esempio 1: dataset "facile" (ben separato, con marker forti)
simulate_spatial_transcriptomics(
  image_path = here::here("images/colon.png"),
  output_path = "data/simulated_easy_correlation.rds",
  output_plot = "results/simulated_easy_distribution.png",
  difficulty_level = "easy",
  n_cells = 20000,
  n_genes = 100
)


# Esempio: genera grafici diagnostici per un dataset simulato
generate_diagnostic_plots(
  input_rds_path = "data/simulated_easy_correlation.rds",
  output_dir = "images/simulation_plots",
  base_name = "sim_analysis"  # opzionale, altrimenti usa il nome del file
)


analyze_and_compare_clusters(rds_path = "data/simulated_improved_bio.rds")

run_simulation_testing("data/simulated_easy_correlation.rds", "results/simulation_tests/")


# Esempio 1: Visium HD-like con griglia a 2µm
simulate_spatial_transcriptomics(
  image_path = here::here("images/colon.png"),
  output_path = "data/simulated_visium_hd.rds",
  output_plot = "results/simulated_visium_hd.png",
  grid_mode = TRUE,                    # Attiva la modalità griglia
  grid_resolution = 2,                 # Risoluzione 2µm come Visium HD
  difficulty_level = "medium",
  n_genes = 100,
  k_cell_types = 5,
  correlation_method = "grf",          # Gaussian Random Field
  spatial_params = list(
    gradient_regions = TRUE,           # Attiva gradienti tra regioni
    gradient_width = 5,                # Larghezza del gradiente (in unità griglia)
    spatial_noise_intensity = 1.0,
    spatial_range = 30,
    random_noise_sd = 0.2
  ),
  dropout_params = list(
    expression_dependent_dropout = TRUE,  # Dropout dipendente dal livello di espressione
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5
  ),
  library_size_params = list(
    mean_library_size = 10000,         # Media conteggi per spot
    library_size_cv = 0.3,             # Coefficiente di variazione
    spatial_effect_on_library = 0.5    # Effetto spaziale sulla dimensione libreria
  )
)
