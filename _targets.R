library(targets)
library(tidyr)
library(dplyr)
library(imager)
library(ggplot2)
library(ClusterR)

# Configura crew per parallelizzazione
library(crew)
controller <- crew::crew_controller_local(
  name = "spatial_transcriptomics",
  workers = 2,  # Numero di worker in parallelo
  seconds_idle = 120,  # Tempo di inattività prima della chiusura automatica
  reset_globals = FALSE  # Preserva variabili globali tra esecuzioni
)

# Configura targets con crew
tar_option_set(
  packages = c(
    "tidyr", "dplyr", "imager", "ggplot2", "ClusterR", 
    "MASS", "cluster", "gstat", "sp", "scales", "spatstat", "spdep", "fields"
  ),
  controller = controller,
  storage = "worker",
  retrieval = "worker",
  memory = "transient",
  garbage_collection = TRUE
)

# Carica funzioni principali
source("R/functions/01_configuration.R")
source("R/functions/02_helper_functions.R")
source("R/functions/07a_simulation_config.R")
source("R/functions/07b_difficulty_setup.R")
source("R/functions/03_image_processing.R")
source("R/functions/04_clustering.R")
source("R/functions/05_grid_sampling.R")
source("R/functions/06a_expression_params.R") 
source("R/functions/06b_expression_baseline.R")
source("R/functions/06c_spatial_distances.R")
source("R/functions/06d_dispersion_params.R")
source("R/functions/06e_library_size.R")
source("R/functions/06f_dropout_models.R")
source("R/functions/06g_gene_modules.R")
source("R/functions/06h_spatial_correlation.R")
source("R/functions/06i_hybrid_cells.R")
source("R/functions/06j_expression_generation.R")
source("R/functions/06k_expression_profiles_wrapper.R")
source("R/functions/07c_simulation_pipeline.R")
source("R/functions/07d_visualization.R")
source("R/functions/07e_results_handling.R")
source("R/functions/07f_simulate_spatial_transcriptomics_wrapper.R")

# Crea directory dei risultati
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# Parametri base (ridotti per test)
params <- list(
  image_path = "images/colon.png",
  n_cells = 500,      # Ridotto per velocità
  n_genes = 30,       # Ridotto per velocità
  k_cell_types = 2,   # Ridotto per semplicità
  difficulty_level = "easy",
  grid_mode = TRUE,
  random_seed = 123
)

# Pipeline
list(
  tar_target(
    config,
    initialize_simulation_config(
      image_path = params$image_path,
      output_path = "results/test_simulation.rds",
      n_cells = params$n_cells,
      n_genes = params$n_genes,
      k_cell_types = params$k_cell_types,
      grid_mode = params$grid_mode,
      random_seed = params$random_seed
    )
  ),
  
  tar_target(
    difficulty,
    configure_difficulty_level(
      difficulty_level = params$difficulty_level
    )
  ),
  
  tar_target(
    image_data,
    prepare_image(config$image_path, config$threshold_value)
  ),
  
  tar_target(
    clustered_data,
    cluster_image(
      img_df_thresh = image_data$img_df_thresh,
      k_cell_types = config$k_cell_types,
      random_seed = config$random_seed
    )
  ),
  
  tar_target(
    cell_df,
    create_sampling_grid(
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
  ),
  
  tar_target(
    expression_profiles,
    generate_expression_profiles(
      cell_df = cell_df,
      n_genes = config$n_genes,
      k_cell_types = config$k_cell_types,
      use_spatial_correlation = TRUE,
      marker_params = list(
        marker_genes_per_type = difficulty$marker_params$marker_genes_per_type,
        marker_expression_fold = difficulty$marker_params$marker_expression_fold
      ),
      dropout_params = list(
        dropout_range = c(0.2, 0.5),
        expression_dependent_dropout = TRUE
      ),
      cell_specific_params = list(
        use_gene_modules = TRUE,
        module_hierarchical = FALSE,
        module_overlap = 0.1,
        cell_specific_noise_sd = 0.2
      ),
      hybrid_params = list(
        use_hybrid_cells = TRUE,
        hybrid_intensity_range = c(0.2, 0.5)
      ),
      random_seed = config$random_seed
    )
  ),
  
  tar_target(
    plots,
    create_simulation_plots(cell_df, config, difficulty, expression_profiles$expression)
  ),
  
  tar_target(
    final_report,
    {
      # Assembla i risultati finali
      result <- list(
        coordinates = cell_df,
        intensity_cluster = cell_df$intensity_cluster,
        expression = expression_profiles$expression,
        parameters = list(
          config = config,
          difficulty = difficulty
        )
      )
      
      # Crea un report con le informazioni principali
      report_content <- c(
        "# Rapporto di Simulazione",
        "",
        paste("Data:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Immagine:", config$image_path),
        paste("Numero di celle:", nrow(cell_df)),
        paste("Tipi cellulari:", config$k_cell_types),
        paste("Geni:", config$n_genes),
        paste("Livello difficoltà:", difficulty$difficulty_level),
        "",
        "## Parametri",
        paste("- Grid mode:", ifelse(config$grid_mode, "SI", "NO")),
        paste("- Risoluzione griglia:", config$grid_resolution, "um"),
        paste("- Random seed:", config$random_seed),
        paste("- Correlazione spaziale:", "SI"),
        paste("- Moduli genici:", "SI"),
        "",
        "## Statistiche celle",
        paste("- Celle totali:", nrow(cell_df)),
        paste("- Distribuzione tipi cellulari:"),
        paste("  -", names(table(cell_df$intensity_cluster)), ":", table(cell_df$intensity_cluster)),
        "",
        "## Statistiche espressione",
        paste("- Dimensione matrice espressione:", paste(dim(expression_profiles$expression), collapse=" x ")),
        paste("- Conteggio medio per cella:", round(mean(rowSums(expression_profiles$expression)), 2)),
        paste("- Percentuale dropout:", round(100 * sum(expression_profiles$expression == 0) / prod(dim(expression_profiles$expression)), 2), "%")
      )
      
      # Salva il report in un file di testo
      report_file <- "results/simulation_report.txt"
      writeLines(report_content, report_file)
      
      # Salva il risultato completo come RDS
      result_file <- "results/simulation_report.rds"
      saveRDS(result, result_file)
      
      # Salva anche il plot
      plot_file <- "results/simulation_plot.png"
      ggplot2::ggsave(plot_file, plot = plots, width = 8, height = 6, dpi = 300)
      
      # Stampa nel terminale
      cat("Simulazione completata con successo!\n")
      cat(sprintf("- Numero di celle: %d\n", nrow(cell_df)))
      cat(sprintf("- Tipi cellulari: %d\n", config$k_cell_types))
      cat(sprintf("- Geni: %d\n", config$n_genes))
      cat(sprintf("- Dimensione matrice espressione: %s\n", paste(dim(expression_profiles$expression), collapse=" x ")))
      cat(sprintf("- Immagine: %s\n", config$image_path))
      cat(sprintf("- Report salvato: %s\n", report_file))
      cat(sprintf("- Dati completi salvati: %s\n", result_file))
      cat(sprintf("- Plot salvato: %s\n", plot_file))
      
      # Restituisci il percorso del report
      list(report_file = report_file, result_file = result_file, plot_file = plot_file)
    }
  )
)