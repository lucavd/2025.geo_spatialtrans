# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Structure

This R package for spatial transcriptomics simulation has been fully modularized with these key components:

- `R/functions/`: Contains all the modularized functions, numbered for dependency order
  - `01_configuration.R`: Configuration settings
  - `02_helper_functions.R`: General utilities
  - `03_image_processing.R`: Image loading and preprocessing
  - `04_clustering.R`: Clustering algorithms
  - `05_grid_sampling.R`: Grid creation and sampling
  - `06*_expression_profiles*.R`: Expression profile generation (11 modules)
  - `07*_simulation*.R`: Main simulation pipeline (6 modules)
- `tests/testthat/`: Contains unit tests for all functions
- `DESCRIPTION`: Package metadata and dependencies
- `run_tests.R`: Script to run all tests with proper library loading
- `_targets.R`: Targets pipeline definition

## Build/Test Commands

- Load all functions for testing: 
  ```r
  # Load all functions in order
  files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
  for (file in sort(files)) { source(file) }
  ```

- Run tests with proper library loading:
  ```r
  source("run_tests.R")
  ```

- Run tests:
  ```r
  testthat::test_dir("tests/testthat/")
  ```

- Run a specific test:
  ```r
  testthat::test_file("tests/testthat/test-expression_modules.R")
  ```

- Build package:
  ```r
  devtools::document()
  devtools::build()
  ```

## Targets Integration

The package now includes a functional targets pipeline for spatial transcriptomics simulation.

1. Run the complete pipeline:
   ```r
   library(targets)
   library(tidyr)
   library(dplyr)
   library(imager)
   library(ggplot2)
   library(ClusterR)
   tar_make()
   ```

2. Visualize the pipeline:
   ```r
   tar_visnetwork()
   ```

3. View generated results in the `results/` directory:
   - `simulation_report.txt`: Textual report with simulation details
   - `simulation_plot.png`: Visualization of the simulated data

## Key Features of the Targets Pipeline

- **Modular design**: Each step of the simulation is a separate target
- **Reproducibility**: Fixed random seeds ensure reproducible results 
- **Visualization**: Automatic plotting and report generation
- **Error handling**: Built-in workspace preservation for debugging

## Next Development Steps

1. **Complete remaining modularization**:
   - Finish modularizing `analyze_and_compare_clusters.R`
   - Finish modularizing `generate_synthetic_tissue.R`

2. **Extend targets integration**:
   - Add parameter-based branching for multiple simulations
   - Set up distributed computing options
   - Add more comprehensive reporting

3. **Add example data and vignettes**:
   - Create demo datasets
   - Write tutorial vignettes

## Modular Structure Benefits

The package has been modularized to support efficient pipeline execution:

1. **Expression Profiles (06*.R files)**:
   - `06a_expression_params.R`: Parameter initialization
   - `06b_expression_baseline.R`: Baseline profiles
   - `06c_spatial_distances.R`: Distance calculations
   - `06d_dispersion_params.R`: Dispersion parameters
   - `06e_library_size.R`: Library size generation
   - `06f_dropout_models.R`: Dropout modeling
   - `06g_gene_modules.R`: Gene module generation
   - `06h_spatial_correlation.R`: Spatial correlation
   - `06i_hybrid_cells.R`: Hybrid cell handling
   - `06j_expression_generation.R`: Matrix generation
   - `06k_expression_profiles_wrapper.R`: Wrapper function

2. **Simulation Pipeline (07*.R files)**:
   - `07a_simulation_config.R`: Configuration
   - `07b_difficulty_setup.R`: Difficulty parameters
   - `07c_simulation_pipeline.R`: Main pipeline
   - `07d_visualization.R`: Visualization
   - `07e_results_handling.R`: Results management
   - `07f_simulate_spatial_transcriptomics_wrapper.R`: Wrapper function

## Code Style Guidelines

- **Indentation**: 2 spaces
- **Function Names**: snake_case (e.g., `calculate_dispersion_params`)
- **Parameter Structure**: Group related parameters in lists
- **Documentation**: All function parameters documented with roxygen style
- **Error Handling**: Use `stop()` for errors, `warning()` for warnings, `tryCatch()` for exceptions
- **Testing**: Each module has corresponding tests