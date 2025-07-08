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


## ✅ Clustering Improvements (Gennaio 2025)

### New Clustering Methods Implemented

Il sistema di clustering è stato esteso con metodi biologicamente realistici per sostituire k-means:

#### Pipeline DBSCAN + Graph Clustering
- **File**: `R/functions/04_clustering.R`
- **Metodo**: `clustering_method = "dbscan_graph"`
- **Strategia**: Pipeline consecutiva per forme irregolari biologicamente plausibili

```r
# Utilizzo nel codice di simulazione
clust <- cluster_image(
  img_df_thresh = img_dat$img_df_thresh,
  k_cell_types  = cfg$k_cell_types,
  random_seed   = cfg$random_seed,
  clustering_method = "dbscan_graph"  # Nuovo metodo
)
```

#### Parametri Biologicamente Controllabili
1. **DBSCAN (Fase 1)**: Identifica regioni dense
   - `eps_factor = 1.4`: Dimensione regioni dense
   - `min_samples = 4`: Numero minimo di celle per cluster

2. **Graph Clustering (Fase 2)**: Raffina strutture locali
   - `k_neighbors = 10`: Dimensione microambiente locale
   - `resolution = 1.0`: Granularità clustering Louvain

#### Vantaggi vs k-means
- **Forme irregolari**: Trova vasi sanguigni, ramificazioni, infiltrazione immune
- **Controllo biologico**: Parametri interpretabili biologicamente
- **Robustezza**: Fallback automatico a spatial_kmeans se necessario
- **Performance**: Scala a migliaia di punti con O(n log n)

### Metodi di Clustering Disponibili
- `"spatial_kmeans"`: K-means con peso spaziale (default precedente)
- `"kmeans++"`: K-means standard migliorato
- `"slic"`: Superpixel clustering
- `"dbscan_graph"`: **NUOVO** - Pipeline consecutiva per forme biologiche

## Next Development Steps

1. **✅ Advanced Clustering (COMPLETATO)**:
   - ✅ Implementata pipeline DBSCAN + Graph clustering
   - ✅ Integrazione con sistema di simulazione esistente
   - ✅ Parametri biologicamente interpretabili

2. **Validation & Optimization**:
   - Validazione biologica risultati vs k-means
   - Performance optimization per large-scale simulations
   - Sviluppo metriche di validazione automatiche

3. **Complete remaining modularization**:
   - Finish modularizing `analyze_and_compare_clusters.R`
   - Finish modularizing `generate_synthetic_tissue.R`

4. **Implement package infrastructure**:
   - Complete roxygen documentation
   - Set up proper package namespace
   - Create package installation workflows

5. **Add example data and vignettes**:
   - Create demo datasets
   - Write tutorial vignettes

6. **Explore advanced spatial modeling**:
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

## Simulation Scripts

### Simplified Simulation (2000 genes)
```r
Rscript R/run_full_simulation_simple.R
```
- Uses biologically validated parameters
- ~2000 genes, ~8000 UMI/cell mean library size
- Includes biological validation report

### Full-size Simulation (20000 genes)
```r
# Quick test first (recommended)
Rscript R/test_full_size_quick.R

# Full simulation
Rscript R/run_full_size_optimized.R
```
- Biologically realistic: 20k genes, 8k UMI/cell
- Memory optimized with chunking
- Full biological validation

### Biological Validation
```r
source("R/biological_validation_report.R")
generate_biological_validation_report(
  sim_results = readRDS("results/your_simulation.rds"),
  output_dir = "plots/biological_validation",
  simulation_name = "Your Simulation"
)
```

## Biologically Validated Parameters

Based on extensive testing and validation:

1. **Library Size**: 8000 UMI/cell (mean) with 30% CV
2. **Gene Expression Distribution**:
   - 85% low-expressed genes (mu ~ -4.5)
   - 10% medium-expressed genes (mu ~ -1.5)  
   - 5% high-expressed genes (mu ~ 0.5)
3. **Dropout**: 40-60% range for medium difficulty
4. **Dispersion**: 10.0-5.0 range for negative binomial

## Code Style Guidelines

- **Indentation**: 2 spaces
- **Function Names**: snake_case (e.g., `calculate_dispersion_params`)
- **Parameter Structure**: Group related parameters in lists
- **Documentation**: All function parameters documented with roxygen style
- **Error Handling**: Use `stop()` for errors, `warning()` for warnings, `tryCatch()` for exceptions
- **Testing**: Each module has corresponding tests
- **Plots**: All plots must have white backgrounds (not transparent) for consistency and clear visualization
  - Use `theme(panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA))` for ggplot2
  - Set `bg = "white"` in all `ggsave()` calls

## Technical Improvement Process

When implementing technical improvements to the codebase:

1. **Clearly define the improvement idea**
2. **Implement changes methodically** without shortcuts or data fabrication
3. **Update/add tests** to validate the improvements
4. **Run and fix tests** until they pass successfully
5. **Update documentation** to reflect the changes
6. **Update CLAUDE.md** with any new processes or guidelines

## Workflow Guidelines

- **R Package Installation**:
  - Do not install R packages. Suggest to the user to install them and then execute the commands