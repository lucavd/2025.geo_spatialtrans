# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Structure

This R package for spatial transcriptomics simulation has been modularized with these key components:

- `R/functions/`: Contains all the modularized functions, numbered for dependency order
- `tests/testthat/`: Contains unit tests for all functions
- `DESCRIPTION`: Package metadata and dependencies

## Build/Test Commands

- Load all functions for testing: 
  ```r
  # Load all functions in order
  files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
  for (file in sort(files)) { source(file) }
  ```

- Run tests:
  ```r
  testthat::test_dir("tests/testthat/")
  ```

- Run a specific test:
  ```r
  testthat::test_file("tests/testthat/test-helper_functions.R")
  ```

- Future package build (when ready):
  ```r
  devtools::document()
  devtools::build()
  ```

## Next Development Steps

1. **Complete the modularization**:
   - Finish modularizing `analyze_and_compare_clusters.R`
   - Finish modularizing `generate_synthetic_tissue.R`

2. **Set up renv for dependency management**:
   ```r
   install.packages("renv")
   renv::init()
   ```

3. **Set up targets for pipeline workflows**:
   ```r
   install.packages("targets")
   # Create _targets.R configuration file
   ```

4. **Add example data and vignettes**

## Code Style Guidelines

- **Indentation**: 2 spaces (NumSpacesForTab: 2)
- **Function Names**: snake_case (e.g., `simulate_spatial_transcriptomics`)
- **Parameter Structure**: Group related parameters in lists (e.g., `marker_params`, `spatial_params`)
- **Documentation**: All function parameters documented with roxygen style
- **Error Handling**: Use `stop()` for errors, `warning()` for warnings, `tryCatch()` for exceptions
- **Testing**: Each function has corresponding tests in `tests/testthat/`

## Major Components

- **Image Loading & Preprocessing**: `03_image_processing.R`
- **Clustering**: `04_clustering.R`
- **Grid Sampling**: `05_grid_sampling.R`
- **Expression Profile Generation**: `06_expression_profiles.R`
- **Main Simulation Function**: `09_simulate_spatial_transcriptomics.R`