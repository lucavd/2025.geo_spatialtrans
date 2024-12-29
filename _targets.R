# Carica i pacchetti necessari
library(targets)
library(tarchetypes)
library(crew)

# Configura il controller per il calcolo parallelo
controller <- crew_controller_local(
  name = "spatialgeo",
  workers = min(110, parallel::detectCores() - 1),
  seconds_idle = 120
)

# Imposta le opzioni per targets
tar_option_set(
  controller = controller,
  storage = "worker",
  memory = "transient",
  error = "continue"
)

# Definizione dei target
list(
  # Target per creare le directory necessarie
  tar_target(
    name = setup_dirs,
    command = {
      for (dir in c("results", "data")) {
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
      }
      TRUE
    }
  ),

  # File di input come dipendenze
  tar_target(
    name = simulation_file,
    command = "R/01_simulate_data.R",
    format = "file"
  ),
  tar_target(
    name = seurat_file,
    command = "R/02_bio_gauss_seurat.R",
    format = "file"
  ),
  tar_target(
    name = tweedie_file,
    command = "R/03_bio_nongauss_tweedie.R",
    format = "file"
  ),
  tar_target(
    name = gam_file,
    command = "R/04_geo_gauss_gam.R",
    format = "file"
  ),
  tar_target(
    name = lgcp_file,
    command = "R/05_geo_nongauss_lgcp.R",
    format = "file"
  ),
  tar_target(
    name = rf_file,
    command = "R/06_geo_nongauss_rf.R",
    format = "file"
  ),

  tar_target(
    name = gnn_file,
    command = "R/07_GNN.R",
    format = "file"
  ),

  # Target per la simulazione dei dati
  tar_target(
    name = spatial_simulation,
    command = {
      if (setup_dirs) {  # Use setup_dirs as a dependency
        source(simulation_file)
        list(
          high = readRDS("data/simulated_high_correlation.rds"),
          medium = readRDS("data/simulated_medium_correlation.rds"),
          low = readRDS("data/simulated_low_correlation.rds")
        )
      }
    },
    format = "rds"
  ),

  tar_target(
    name = seurat_analysis,
    command = {
      source(seurat_file)
      spatial_simulation  # Use as dependency to ensure it runs after simulation
      list(
        high = readRDS("results/seurat_results_simulated_high_correlation.rds"),
        medium = readRDS("results/seurat_results_simulated_medium_correlation.rds"),
        low = readRDS("results/seurat_results_simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  tar_target(
    name = tweedie_analysis,
    command = {
      source(tweedie_file)
      spatial_simulation  # Use as dependency
      list(
        high = readRDS("results/tweedie_results_simulated_high_correlation.rds"),
        medium = readRDS("results/tweedie_results_simulated_medium_correlation.rds"),
        low = readRDS("results/tweedie_results_simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  tar_target(
    name = gam_analysis,
    command = {
      source(gam_file)
      spatial_simulation  # Use as dependency
      list(
        high = readRDS("results/gam_geospatial_results_simulated_high_correlation.rds"),
        medium = readRDS("results/gam_geospatial_results_simulated_medium_correlation.rds"),
        low = readRDS("results/gam_geospatial_results_simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  tar_target(
    name = lgcp_analysis,
    command = {
      source(lgcp_file)
      spatial_simulation  # Use as dependency
      list(
        high = readRDS("results/lgcp_results_simulated_high_correlation.rds"),
        medium = readRDS("results/lgcp_results_simulated_medium_correlation.rds"),
        low = readRDS("results/lgcp_results_simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  tar_target(
    name = rf_analysis,
    command = {
      source(rf_file)
      spatial_simulation  # Use as dependency
      list(
        high = readRDS("results/ranger_results_simulated_high_correlation.rds"),
        medium = readRDS("results/ranger_results_simulated_medium_correlation.rds"),
        low = readRDS("results/ranger_results_simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  tar_target(
    name = gnn_analysis,
    command = {
      source(gnn_file)
      spatial_simulation  # Use as dependency
      list(
        high = readRDS("results/gnn_results_simulated_high_correlation.rds"),
        medium = readRDS("results/gnn_results_simulated_medium_correlation.rds"),
        low = readRDS("results/gnn_results_simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  # Target per l'analisi finale dei risultati
  tar_target(
    name = final_results,
    command = {
      # Include all analysis results as dependencies
      seurat_analysis
      tweedie_analysis
      gam_analysis
      lgcp_analysis
      rf_analysis
      gnn_analysis

      source("R/results.R")
      list(
        metrics = read.csv("results/all_correlation_levels_metrics.csv"),
        plots = list(
          metrics = "results/all_correlation_levels_metrics_plot.png",
          errors = "results/all_correlation_levels_fp_fn_errors.png",
          f1_comparison = "results/f1_score_comparison.png"
        )
      )
    },
    format = "rds"
  )
)
