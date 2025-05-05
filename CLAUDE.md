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
  
- Note: The test suite requires increasing the future.globals.maxSize limit to 100 GB:
  ```r
  options(future.globals.maxSize = 100 * 1024^2) # 100 GB
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


## Next Development Steps

1. **Complete remaining modularization**:
   - Finish modularizing `analyze_and_compare_clusters.R`
   - Finish modularizing `generate_synthetic_tissue.R`

2. **Implement package infrastructure**:
   - Complete roxygen documentation
   - Set up proper package namespace
   - Create package installation workflows

3. **Add example data and vignettes**:
   - Create demo datasets
   - Write tutorial vignettes

4. **Explore advanced spatial modeling**:
   - Test multi-scale and hierarchical spatial models
   - Implement anisotropic and non-stationary patterns
   - Create composite spatial patterns for realistic tissues

## Modular Structure Benefits

The package has been modularized to support efficient pipeline execution:

1. **Expression Profiles (06*.R files)**:
   - `06a_expression_params.R`: Parameter initialization
   - `06b_expression_baseline.R`: Baseline profiles
   - `06c_spatial_distances.R`: Distance calculations
   - `06d_dispersion_params.R`: Dispersion parameters
   - `06e_library_size.R`: Library size generation
   - `06f_dropout_models.R`: Dropout modeling with ambient RNA
   - `06g_gene_modules.R`: Gene module generation
   - `06h_spatial_correlation.R`: Basic spatial correlation
   - `06h_spatial_correlation_multiscale.R`: Multi-scale correlations
   - `06h_spatial_correlation_nonstationary.R`: Non-stationary patterns
   - `06i_hybrid_cells.R`: Hybrid cell handling
   - `06j_expression_generation.R`: Matrix generation
   - `06k_expression_profiles_wrapper.R`: Wrapper function
   - `06l_ligand_receptor_interactions_simple.R`: Ligand-receptor interactions (simplified)
   - `06m_temporal_dynamics_simple.R`: Temporal dynamics (simplified)
   - `06n_alternative_splicing_simple.R`: Alternative splicing (simplified)
   - `06o_anisotropic_patterns_simple.R`: Anisotropic patterns (simplified)
   - `06p_3d_microenvironment_simple.R`: 3D microenvironment (simplified)

2. **Simulation Pipeline (07*.R files)**:
   - `07a_simulation_config.R`: Configuration
   - `07b_difficulty_setup.R`: Difficulty parameters
   - `07c_simulation_pipeline.R`: Main pipeline
   - `07d_visualization.R`: Visualization
   - `07e_results_handling.R`: Results management
   - `07f_simulate_spatial_transcriptomics_wrapper.R`: Wrapper function

3. **Validation and Analysis (08*.R files)**:
   - `08_validation_plots.R`: Comprehensive validation plots generation

## Code Style Guidelines

- **Indentation**: 2 spaces
- **Function Names**: snake_case (e.g., `calculate_dispersion_params`)
- **Parameter Structure**: Group related parameters in lists
- **Documentation**: All function parameters documented with roxygen style
- **Error Handling**: Use `stop()` for errors, `warning()` for warnings, `tryCatch()` for exceptions
- **Testing**: Each module has corresponding tests

## Technical Improvement Process

When implementing technical improvements to the codebase:

1. **Clearly define the improvement idea**
2. **Implement changes methodically** without shortcuts or data fabrication
3. **Update/add tests** to validate the improvements
4. **Run and fix tests** until they pass successfully
5. **Update documentation** to reflect the changes
6. **Update CLAUDE.md** with any new processes or guidelines