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
- `verify_modules.R`: Script to verify modules function correctly
- `verify_targets.R`: Example targets integration

## Build/Test Commands

- Load all functions for testing: 
  ```r
  # Load all functions in order
  files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
  for (file in sort(files)) { source(file) }
  ```

- Run tests with detailed report:
  ```r
  source("verify_modules.R")
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

To use with targets:

1. Examine the example in `verify_targets.R`
2. Run the pipeline:
   ```r
   library(targets)
   tar_make()
   ```
3. Visualize pipeline:
   ```r
   tar_visnetwork()
   ```

## Next Development Steps

1. **Complete remaining modularization**:
   - Finish modularizing `analyze_and_compare_clusters.R`
   - Finish modularizing `generate_synthetic_tissue.R`

2. **Implement full targets integration**:
   - Refine the `_targets.R` file
   - Add persistent caching
   - Set up distributed computing options

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