library(targets)
library(tidyr)
library(dplyr)
library(imager)
library(ggplot2)
library(ClusterR)

# Carica funzioni principali
source("R/functions/01_configuration.R")
source("R/functions/02_helper_functions.R")
source("R/functions/07a_simulation_config.R")
source("R/functions/07b_difficulty_setup.R")
source("R/functions/03_image_processing.R")
source("R/functions/04_clustering.R")
source("R/functions/05_grid_sampling.R")
source("R/functions/07d_visualization.R")

# Crea directory dei risultati
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# Parametri base
params <- list(
  image_path = "images/generated2.png",
  n_cells = 1000,
  n_genes = 50,
  k_cell_types = 3,
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
    plots,
    create_simulation_plots(cell_df, config, difficulty)
  ),
  
  tar_target(
    final_report,
    {
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
        "",
        "## Statistiche celle",
        paste("- Celle totali:", nrow(cell_df)),
        paste("- Distribuzione tipi cellulari:"),
        paste("  -", names(table(cell_df$intensity_cluster)), ":", table(cell_df$intensity_cluster))
      )
      
      # Salva il report in un file di testo
      report_file <- "results/simulation_report.txt"
      writeLines(report_content, report_file)
      
      # Salva anche il plot
      plot_file <- "results/simulation_plot.png"
      ggplot2::ggsave(plot_file, plot = plots, width = 8, height = 6, dpi = 300)
      
      # Stampa nel terminale
      cat("Simulazione completata con successo!\n")
      cat(sprintf("- Numero di celle: %d\n", nrow(cell_df)))
      cat(sprintf("- Tipi cellulari: %d\n", config$k_cell_types))
      cat(sprintf("- Immagine: %s\n", config$image_path))
      cat(sprintf("- Report salvato: %s\n", report_file))
      cat(sprintf("- Plot salvato: %s\n", plot_file))
      
      # Restituisci il percorso del report
      list(report_file = report_file, plot_file = plot_file)
    }
  )
)