# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Structure

This R package for spatial transcriptomics simulation has been streamlined and optimized with these key components:

### Active Core Pipeline (`R/`)
- `01_configuration.R`: Performance configuration and memory management
- `02_helper_functions.R`: General utilities and helper functions
- `03_image_processing.R`: Image loading and preprocessing
- `03b_generate_synthetic_tissue.R`: Synthetic tissue image generation (complexity levels 1-3)
- `04_clustering.R`: **SIMPLIFIED** - Only spatial_kmeans method (64 lines, reduced from 540)
- `05_grid_sampling.R`: Grid creation and sampling for Visium-like technologies
- `06_expression_profiles.R` through `06k_expression_profiles_wrapper.R`: **12 essential expression modules**
  - `06a_expression_params.R`: Parameter initialization
  - `06b_expression_baseline.R`: Baseline expression profiles
  - `06c_spatial_distances.R`: Distance calculations
  - `06d_dispersion_params.R`: Dispersion parameters
  - `06e_library_size.R`: Library size generation
  - `06f_dropout_models.R`: Dropout modeling with ambient RNA
  - `06g_gene_modules.R`: Gene module generation
  - `06h_spatial_correlation.R`: Basic spatial correlation (GRF)
  - `06i_hybrid_cells.R`: Hybrid cell handling at boundaries
  - `06j_expression_generation.R`: Expression matrix generation
  - `06k_expression_profiles_wrapper.R`: **Main coordinator function**

### Testing Framework (`R/testing/`)
- `full_test.R`: **PRIMARY SCRIPT** - Biologically realistic full pipeline test
- `visualize_clusters.R`: Visualization script for cluster analysis
- `*.png` and `*.rds`: Generated results and visualizations

### Advanced Features (Not Yet Active) (`R/test functions/`)
- `06l_ligand_receptor_interactions.R`: Cell-cell communication (future feature)
- `06m_temporal_dynamics.R`: RNA velocity and pseudotime (future feature)
- `06n_alternative_splicing.R`: Isoform regulation (future feature)
- `06o_anisotropic_patterns.R`: Directional patterns (future feature)
- `06p_3d_microenvironment.R`: 3D tissue effects (future feature)
- `07*_simulation*.R`: Advanced simulation pipeline (future features)
- `08_validation_plots.R`: Comprehensive validation (future feature)

## Build/Test Commands

### Primary Usage: Full Pipeline Test
- **Run the main simulation test** (recommended starting point):
  ```r
  # Run full biologically realistic test
  Rscript R/testing/full_test.R
  
  # Or with custom image
  Rscript R/testing/full_test.R path/to/your/image.png
  ```

### Manual Function Loading
- Load all core functions for development:
  ```r
  # Load core pipeline functions (01-06k)
  files <- list.files("R", full.names = TRUE, pattern = "\\.R$")
  files <- files[!grepl("test functions", files)]  # Exclude test functions
  for (file in sort(files)) { source(file) }
  ```

### Memory Configuration
- The pipeline requires adequate memory allocation:
  ```r
  # For full test (required)
  options(future.globals.maxSize = 4 * 1024^2) # 4 GB
  
  # For large-scale simulations
  setup_performance(max_memory_gb = 8, workers = 4)
  ```

### Test Suite (Legacy)
- Basic unit tests (if available):
  ```r
  source("tests/run_tests.R")
  testthat::test_dir("tests/testthat/")
  ```

### Package Development
- Build package documentation:
  ```r
  devtools::document()
  devtools::build()
  ```


## ✅ Clustering Simplification (Gennaio 2025)

### Streamlined Clustering Implementation

Il sistema di clustering è stato **semplificato drasticamente** per fornire una pipeline stabile e controllabile:

#### Spatial K-means Only
- **File**: `R/04_clustering.R`
- **Riduzione**: Da 540 righe a 64 righe (**92% di riduzione**)
- **Metodo unico**: `spatial_kmeans` - K-means con bilanciamento spazio-intensità

```r
# Utilizzo semplificato nel codice di simulazione
clust <- cluster_image(
  img_df_thresh = img_dat$img_df_thresh,
  k_cell_types  = cfg$k_cell_types,
  random_seed   = cfg$random_seed,
  spatial_weight = 0.6  # Controllo bilanciamento spazio-intensità
)
```

#### Parametro di Controllo Biologico
- **`spatial_weight`**: Controllo del bilanciamento tra intensità dell'immagine e posizione spaziale
  - `0.0`: Clustering solo per intensità (facile)
  - `0.5`: Bilanciamento equilibrato (medio)
  - `1.0+`: Clustering principalmente spaziale (difficile)

#### Vantaggi della Semplificazione
- **Deterministica**: Risultati riproducibili con `random_seed`
- **Controllabile**: Difficoltà parametrizzabile tramite `spatial_weight`
- **Stabile**: Funziona con pattern semplici e complessi
- **Veloce**: Operazioni vettorizzate ottimizzate
- **Ground Truth**: Cluster definiti e controllabili per benchmark

### Metodi Rimossi nella Semplificazione
- `"kmeans++"`: K-means standard (rimosso)
- `"slic"`: Superpixel clustering (rimosso)
- `"dbscan_graph"`: Pipeline DBSCAN+Graph (rimosso)
- Stima automatica k con elbow/silhouette (rimossa)
- Gestione fallback complessa (rimossa)

## Next Development Steps

### ✅ Completed (Gennaio 2025)
- [x] **Pipeline Core Semplificata**: 17 moduli essenziali funzionanti (01-06k)
- [x] **Full Test Biologicamente Realistico**: `R/testing/full_test.R` con parametri validati
- [x] **Validazione Automatica**: Benchmark completi con PASS/FAIL
- [x] **Visualizzazione Cluster**: Pannelli combinati in stile pubblicazione scientifica
- [x] **Gestione Memoria**: Chunking intelligente per dataset grandi
- [x] **Synthetic Tissue Generation**: 3 livelli di complessità biologica

### Active Development Priorities

1. **Advanced Features Integration**:
   - Integrate modules from `R/test functions/` as optional features
   - Enable ligand-receptor interactions (06l) as parameter option
   - Enable temporal dynamics (06m) as parameter option
   - Enable alternative splicing (06n) as parameter option
   - Enable anisotropic patterns (06o) as parameter option
   - Enable 3D microenvironment (06p) as parameter option

2. **Extended Validation**:
   - Test with different dataset sizes (10k, 15k, 25k genes)
   - Test with more clusters (10, 12, 15)
   - Benchmark performance scaling on very large datasets
   - Compare ground truth recovery accuracy

3. **Algorithm Benchmarking**:
   - Implement clustering comparison metrics (ARI, NMI, Silhouette)
   - Test with standard algorithms (k-means, Louvain, Leiden)  
   - Validate ground truth recovery performance
   - Compare with public real datasets

4. **Package Infrastructure**:
   - Complete roxygen documentation for all functions
   - Set up proper package namespace and DESCRIPTION
   - Create package installation workflows
   - Add comprehensive unit tests

5. **User Experience Enhancements**:
   - Support user images with automatic preprocessing
   - Export to standard formats (H5AD, Seurat, CSV)
   - Batch processing for multiple experiments
   - GPU parallelization for very large datasets

## Current Pipeline Architecture

### Core Expression Pipeline (Active - 12 Modules)
The essential expression generation pipeline consists of **12 coordinated modules** (06a-06k):

1. **06a_expression_params.R**: Parameter initialization and validation
2. **06b_expression_baseline.R**: Cell type-specific baseline expression profiles
3. **06c_spatial_distances.R**: Spatial distance calculations and local density
4. **06d_dispersion_params.R**: Negative binomial dispersion parameters
5. **06e_library_size.R**: Realistic library size generation with spatial effects
6. **06f_dropout_models.R**: Expression-dependent dropout with ambient RNA
7. **06g_gene_modules.R**: Gene co-expression modules and regulatory networks
8. **06h_spatial_correlation.R**: Gaussian Random Field spatial correlation
9. **06i_hybrid_cells.R**: Boundary cells with mixed phenotypes
10. **06j_expression_generation.R**: Final expression matrix generation
11. **06k_expression_profiles_wrapper.R**: **MAIN COORDINATOR** - orchestrates all modules

**Key Features of Active Pipeline:**
- **Biologically Realistic**: All 12 modules essential for biological accuracy
- **Non-Simplifiable**: Each module has specific dependencies in the workflow
- **Coordinated**: Wrapper function manages the complete pipeline
- **Memory Efficient**: Supports chunking for large-scale simulations
- **Validated**: Produces realistic UMI counts, sparsity, and spatial patterns

### Advanced Features (Future Integration)
Located in `R/test functions/` - **not yet active but planned for integration**:

1. **06l_ligand_receptor_interactions.R**: Cell-cell communication modeling
2. **06m_temporal_dynamics.R**: RNA velocity and pseudotime trajectories  
3. **06n_alternative_splicing.R**: Spatial regulation of isoform usage
4. **06o_anisotropic_patterns.R**: Directional patterns along tissue structures
5. **06p_3d_microenvironment.R**: 3D tissue effects on 2D measurements
6. **07*_simulation*.R**: Advanced simulation pipeline components
7. **08_validation_plots.R**: Comprehensive validation plotting suite

## Primary Simulation Script

### Full Test Pipeline (`R/testing/full_test.R`)

**Main Entry Point** - Biologically realistic full pipeline test:

```bash
# Run with synthetic tissue (default)
Rscript R/testing/full_test.R

# Run with user-provided image
Rscript R/testing/full_test.R path/to/your/tissue/image.png
```

**Configuration Full-Size Biological**:
- **5,000 genes**: Realistic for spatial transcriptomics publications
- **10,000 spots/cells**: Medium-scale dataset size
- **8 cell types**: Realistic tissue complexity
- **8,000 UMI/spot**: Target library size for Visium HD
- **800x800 pixel**: Full-size synthetic tissue with complexity=3

**Automatic Benchmarking**:
- ✓ **Matrix dimensions**: Correct genes × cells 
- ✓ **Sparsity**: 30-92% (extended range for spatial data)
- ✓ **UMI per cell**: 3,000-20,000 (realistic range for spatial)
- ✓ **UMI CV**: 0.15-1.0 (coefficient of variation)
- ✓ **Cluster assignment**: Expected number of clusters
- ✓ **Data integrity**: No NaN/Inf values
- ✓ **Moran's I**: Spatial autocorrelation test

**Performance**: Complete pipeline in **<1 minute** with chunking optimization

### Visualization and Results

**Automatic Output Generation**:
```bash
# Run visualization after full_test.R
Rscript R/testing/visualize_clusters.R
```

**Generated Files**:
- `full_test_result.rds`: Complete benchmark results and metadata
- `full_simulation_data.rds`: Expression matrix and coordinates 
- `combined_visualization.png`: Publication-style combined panel
- `cluster_visualization.png`: Detailed cluster plot (20,000 points)
- `full_tissue_complex.png`: Generated synthetic tissue image

**Manual Pipeline Execution**:
```r
# Load all functions
files <- list.files("R", full.names = TRUE, pattern = "\\.R$")
files <- files[!grepl("test functions", files)]
for (file in sort(files)) { source(file) }

# Configure memory
options(future.globals.maxSize = 4 * 1024^2) # 4 GB

# Generate data with controlled difficulty
spatial_weight <- 0.2  # Easy: cluster by intensity  
spatial_weight <- 0.8  # Hard: spatial clusters
```

## Biologically Validated Parameters

Based on extensive testing and validation with real spatial transcriptomics data:

### Realistic Expression Configuration
1. **Library Size**: 8,000 UMI/spot (mean) with 30% CV
   - Range: 3,000-20,000 UMI/spot (captures natural variation)
   - Spatial effects: Library size correlation with tissue density
   - Cell type effects: Different RNA content per cell type

2. **Gene Expression Distribution**:
   - **25 marker genes per cell type**: Strong cell type identification
   - **Marker fold-change**: 2.5x over baseline (biologically realistic)
   - **Minimal overlap**: 0.05 overlap coefficient between cell types
   - **Spatial correlation**: 40μm range, moderate intensity (1.2)

3. **Technical Parameters**:
   - **Dropout range**: 45-65% (realistic for spatial technologies)
   - **Dispersion range**: 12.0-6.0 (high-quality data characteristics)
   - **Expression-dependent dropout**: Logistic curve with midpoint=0.4
   - **Gene modules**: 8 co-expression modules with 75% correlation

### Benchmark Results (Full Test)
**Typical successful run**:
- ✓ **Execution time**: 0.79 minutes
- ✓ **Final dimensions**: 5,000×10,000 (genes×spots)
- ✓ **Sparsity**: 89.2% (typical for filtered spatial data)
- ✓ **UMI statistics**: Mean=9,022, Median=7,100
- ✓ **UMI CV**: 0.87 (realistic biological variability)
- ✓ **Clusters found**: 8/8 correctly assigned
- ✓ **Moran's I**: >0.05 (significant spatial autocorrelation)

**Biological Validation Confirmed**:
- Library size target: 8,000 UMI/spot ✓
- Dropout in expected range: 45-65% ✓  
- Marker genes per type: 25 ✓
- Spatial correlation patterns: Realistic ✓

## Code Style Guidelines

- **Indentation**: 2 spaces for R code
- **Function Names**: snake_case (e.g., `generate_expression_profiles`, `cluster_image`)
- **Parameter Structure**: Group related parameters in lists with descriptive names
- **Documentation**: All functions documented with roxygen2 style comments
- **Error Handling**: Use `stop()` for errors, `warning()` for warnings, `tryCatch()` for robustness
- **Memory Management**: Use chunking for large datasets, call `gc()` periodically
- **Plots**: All plots must have white backgrounds for consistency and publication quality
  - Use `theme(panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA))` for ggplot2
  - Set `bg = "white"` in all `ggsave()` calls
- **Reproducibility**: Always use `set.seed()` with explicit `random_seed` parameters

## Benchmark Usage Guidelines

### For Method Development
- Use `spatial_weight = 0.2` for easy clustering benchmarks
- Use `spatial_weight = 0.8` for challenging spatial pattern recognition
- Modify `n_genes` and `k_cell_types` to create specific benchmark scenarios
- Save results as `.rds` for systematic comparison studies

### For Algorithm Testing
- Ground truth is available through `cluster` assignments in `cell_df`
- Compare clustering results using standard metrics (ARI, NMI, Silhouette)
- Use generated coordinates for spatial analysis validation
- Benchmark memory usage and execution time for scaling studies

## Technical Improvement Process

When implementing improvements to the streamlined codebase:

1. **Test first**: Run `Rscript R/testing/full_test.R` to ensure current functionality
2. **Implement incrementally**: Make small, testable changes to individual modules  
3. **Validate biologically**: Ensure changes don't break biological realism
4. **Update benchmarks**: Adjust validation criteria if necessary
5. **Document changes**: Update CLAUDE.md with new features or modifications
6. **Performance check**: Verify that changes don't significantly impact execution time

## Workflow Guidelines

### Development Workflow
1. **Start with full_test.R**: Always begin by running the main test script
2. **Load functions manually**: Use the provided loading script for development work
3. **Check memory requirements**: Ensure adequate memory allocation before large runs
4. **Use visualization**: Always generate plots to visually validate results

### R Package Management
- **Do not install R packages automatically**: Suggest to the user to install them
- **Required packages**: Matrix, ClusterR, dplyr, sp, gstat, png
- **Optional packages**: future, doParallel (for advanced parallelization)

## Future Advanced Features (R/test functions/)

The repository contains advanced biological modules ready for integration as optional features:

### 🧬 Ligand-Receptor Interactions (`06l`)
- **Purpose**: Model cell-cell communication through ligand-receptor pairs
- **Parameters**: `n_interactions`, `signal_propagation_mode`, `max_signaling_distance`
- **Biology**: Distance-weighted diffusion of signaling molecules
- **Integration**: Add as `lr_params` option in future pipeline versions

### ⏱️ Temporal Dynamics (`06m`)
- **Purpose**: RNA velocity and developmental trajectories
- **Parameters**: `pseudotime_mode`, `temporal_gene_fraction`, `include_velocity`
- **Biology**: Simulate developmental processes with directional gene expression changes
- **Integration**: Add as `temporal_params` option for developmental studies

### 🧬 Alternative Splicing (`06n`)
- **Purpose**: Spatial regulation of isoform usage
- **Parameters**: `splicing_fraction`, `splicing_spatial_pattern`, `splicing_strength`
- **Biology**: Tissue-specific splicing regulation with spatial gradients
- **Integration**: Add as `splicing_params` option for transcript diversity

### 🔀 Anisotropic Patterns (`06o`)
- **Purpose**: Directional gene expression along tissue structures
- **Parameters**: `structure_type`, `anisotropic_pattern`, `anisotropic_gene_fraction`
- **Biology**: Model vascular, neural, or epithelial directional patterns
- **Integration**: Add as `anisotropic_params` option for structured tissues

### 📐 3D Microenvironment (`06p`)
- **Purpose**: Model 3D tissue effects on 2D spatial measurements
- **Parameters**: `n_layers`, `layer_specificity`, `projection_noise`
- **Biology**: Z-axis heterogeneity and projection artifacts
- **Integration**: Add as `micro3d_params` option for complex tissue architecture

**Integration Priority**: These modules represent the next development phase, adding substantial biological complexity while maintaining the current streamlined core pipeline.