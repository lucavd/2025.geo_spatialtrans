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

  # Target per la simulazione dei dati
  tar_target(
    name = spatial_simulation,
    command = {
      source(simulation_file)
      list(
        high = readRDS("data/simulated_high_correlation.rds"),
        medium = readRDS("data/simulated_medium_correlation.rds"),
        low = readRDS("data/simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  tar_target(
    name = seurat_analysis,
    command = {
      source(seurat_file)
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
      list(
        high = readRDS("results/rf_spatial_results_simulated_high_correlation.rds"),
        medium = readRDS("results/rf_spatial_results_simulated_medium_correlation.rds"),
        low = readRDS("results/rf_spatial_results_simulated_low_correlation.rds")
      )
    },
    format = "rds"
  ),

  # Target per l'analisi finale dei risultati
  tar_target(
    name = final_results,
    command = {
      source("R/results.R")
      list(
        metrics = read.csv("results/all_correlation_levels_metrics.csv"),
        plots = list(
          metrics = "results/all_correlation_levels_metrics_plot.png",
          errors = "results/all_correlation_levels_fp_fn_errors.png"
        )
      )
    },
    format = "rds"
  )
)
