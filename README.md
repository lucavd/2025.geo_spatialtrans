# Spatial Transcriptomics Simulation and Analysis Framework

This repository contains a comprehensive set of tools for simulating and analyzing spatial transcriptomics data with various statistical properties. The framework enables the generation of biologically plausible spatial gene expression patterns that reflect real-world transcriptomic data, and provides tools for analyzing this data with different computational approaches.

## Project Structure

The project is organized into a series of interconnected R scripts, each handling a specific aspect of the simulation and analysis pipeline:

1. **Data Simulation**
   - `01_simulate_data.R` - Main wrapper for data simulation
   - `01_simulate_data_simple.R` - Simplified version of the simulator
   - `image_simulations.R` - Core simulation engine with the `simulate_spatial_transcriptomics` function

2. **Biological Analysis Methods**
   - `02_bio_gauss_seurat.R` - Seurat-based analysis assuming Gaussian distribution
   - `03_bio_nongauss_tweedie.R` - Analysis with non-Gaussian Tweedie distribution

3. **Geospatial Analysis Methods**
   - `04_geo_gauss_gam.R` - Gaussian GAM-based spatial analysis
   - `05_geo_nongauss_lgcp.R` - Log-Gaussian Cox Process spatial analysis
   - `06_geo_nongauss_rf.R` - Random Forest-based spatial analysis

4. **Advanced Analysis Methods**
   - `07_GNN.R` - Graph Neural Network-based analysis
   - `analyze_and_compare_clusters.R` - Evaluation framework comparing clustering approaches

5. **Results**
   - `results.R` - Aggregation and visualization of analysis results

## Core Simulation Framework

The heart of this project is the `simulate_spatial_transcriptomics` function in `image_simulations.R`, which generates realistic spatial transcriptomics data from image inputs. This simulation addresses key challenges in spatial transcriptomics:

### Theoretical Background

#### Spatial Transcriptomics

Spatial transcriptomics technologies measure gene expression while preserving the spatial information of cells within tissues. Unlike traditional bulk RNA sequencing or single-cell RNA sequencing, these methods provide crucial insights into how gene expression varies across the spatial organization of tissues.

#### Challenges in Simulation

Simulating realistic spatial transcriptomic data presents several challenges:

1. **Count Distribution**: Gene expression counts in real data follow complex distributions, typically exhibiting overdispersion (variance exceeding mean).
2. **Spatial Correlation**: Gene expression in neighboring cells tends to be correlated due to similar microenvironments and cell-cell interactions.
3. **Dropout Effects**: Technical artifacts often result in false zeroes (dropouts), which occur non-randomly across the spatial domain.
4. **Gene-specific Variation**: Different genes exhibit different statistical properties, with some showing high stability (sub-Poisson variance) and others showing high variability.
5. **Cluster Boundaries**: Cells at the boundaries between different cell types or regions often display transitional expression patterns and higher noise levels.

### Simulation Features and Implementation

The simulation framework offers the following features:

#### Customizable Difficulty Levels

The framework provides three pre-defined difficulty levels (`easy`, `medium`, `hard`) that control how challenging the resulting data is for downstream analysis:

```r
# Example: dataset with "hard" difficulty level
simulate_spatial_transcriptomics(
  image_path = "images/colon.png",
  difficulty_level = "hard",
  n_cells = 30000,
  n_genes = 100
)
```

Each difficulty level adjusts multiple parameters simultaneously:

1. **Easy**: Strong marker genes, clear boundaries, low noise, low dropout
2. **Medium**: Moderate marker genes, some boundary overlap, medium noise and dropout
3. **Hard**: Weak marker genes, significant boundary overlap, high noise and dropout

#### Statistical Models Used

The simulation employs several statistical models to create realistic data:

1. **Negative Binomial Distribution**: Models overdispersed count data, controlled by mean (μ) and dispersion (size) parameters. Captures the variance structure observed in real transcriptomic data.
2. **Binomial Distribution**: Used for modeling constitutively expressed genes with sub-Poisson variance (lower than what a Poisson distribution would predict).
3. **Gaussian Process**: Provides a principled way to model continuous spatial correlation patterns across the tissue, allowing for smooth transitions between regions.
4. **Spatially-varying Dispersion**: Models higher biological noise at cluster boundaries and transitional zones by varying the dispersion parameter based on spatial location.
5. **Spatially-varying Dropout**: Models technical artifacts that tend to be more prevalent in certain regions (e.g., tissue edges) by making dropout probability dependent on spatial position.

#### Advanced Customization

For users requiring fine-grained control, the simulation provides parameter lists for advanced customization:

```r
# Example: dataset with custom parameters
simulate_spatial_transcriptomics(
  image_path = "images/tissue.png",
  difficulty_level = "medium",  # base parameters
  n_cells = 25000,
  n_genes = 150,
  k_cell_types = 7,
  marker_params = list(
    marker_genes_per_type = 8,
    marker_expression_fold = 1.0,
    marker_overlap_fold = 0.3
  ),
  spatial_params = list(
    spatial_noise_intensity = 1.2,
    spatial_range = 25,
    random_noise_sd = 0.3
  ),
  dropout_params = list(
    dropout_range = c(0.25, 0.6),
    dispersion_range = c(1.8, 0.9)
  )
)
```

#### Hybrid Cells at Boundaries

A unique feature of this framework is the explicit modeling of hybrid cells at region boundaries:

```r
# Hybrid cells feature can be enabled/disabled
hybrid_params = list(
  use_hybrid_cells = TRUE,
  max_hybrid_pairs = 1000,
  hybrid_intensity_range = c(0.2, 0.5)
)
```

This creates cells that have expression profiles partially matching neighboring clusters, reflecting the biological reality of transitional zones in tissues.

### Hard Clustering Challenge

The `hard` difficulty level implements several specific challenges for clustering algorithms:

1. **Reduced Marker Gene Contrast**
   - Fewer marker genes per cell type (5 instead of 10)
   - Lower expression differential (+0.8 log fold vs. +2.0)
   - Partial expression overlap between adjacent clusters

2. **Increased Spatial Noise**
   - More intense spatial noise effect
   - Shorter spatial correlation range (more irregularity)
   - Additional random noise component

3. **Higher Dropout and Variability**
   - Higher dropout probability throughout the tissue (0.3-0.7)
   - Higher dispersion (more variability) in gene expression

4. **Hybrid Cell Effects**
   - Explicit modeling of cells at cluster boundaries with mixed profiles
   - Gradual transitions rather than sharp boundaries

These modifications make clustering much more challenging for methods like Seurat, better reflecting the complexities of real tissue data.

## Analysis and Cluster Comparison

The `analyze_and_compare_clusters.R` script provides a framework for evaluating clustering performance on simulated data:

```r
# Example usage
analyze_and_compare_clusters(
  rds_path = "data/simulated_image_correlation.rds",
  output_path = "results/metrics_comparison.csv",
  k_cell_types = 5,                # Number of cell types (optional, can be read from file)
  seurat_resolution = 0.3,         # Resolution for Seurat clustering
  hdbscan_minPts = 10              # minPts parameter for HDBSCAN
)
```

This function:

1. Loads simulated data (ground truth)
2. Processes the data through a standard Seurat pipeline
3. Applies both Seurat clustering and HDBSCAN algorithms
4. Computes performance metrics comparing to ground truth
5. Generates visualizations of the clustering results

The metrics calculated include True Positives (TP), False Positives (FP), False Negatives (FN), and Precision for each cluster type, allowing quantitative comparison of different clustering approaches.

## Statistical Properties of the Simulation

### 1. Overdispersion

The Negative Binomial distribution used for most genes models overdispersion in gene expression, a key property of real transcriptomic data. The probability mass function is:

$$P(X=k) = \binom{k+r-1}{k} p^r (1-p)^k$$

Where:
- $r$ is the dispersion parameter (called `size` in R)
- $p$ is related to the mean μ by $p = \frac{r}{r+\mu}$

The variance is related to the mean by:

$$\text{Var}(X) = \mu + \frac{\mu^2}{r}$$

This allows the variance to exceed the mean, unlike the Poisson distribution where variance = mean.

### 2. Sub-Poisson Variation

For stable genes, we use a Binomial distribution with parameters calculated to achieve a sub-Poisson variance:

$$P(X=k) = \binom{n}{k} p^k (1-p)^{n-k}$$

With variance:

$$\text{Var}(X) = np(1-p)$$

When p is large (e.g., 0.9), the variance is less than the mean (np), providing a realistic model for constitutively expressed genes.

### 3. Spatial Correlation via Gaussian Process

The Gaussian Process implements a continuous spatial correlation structure. The exponential covariance function is:

$$C(d) = \sigma^2 \exp(-d/\rho)$$

Where:
- $d$ is the distance between points
- $\rho$ is the range parameter
- $\sigma^2$ is the variance

This creates smooth spatial patterns that mimic the continuous nature of biological processes across tissue.

### 4. Varying Noise at Boundaries

The spatially-varying dispersion and dropout models capture important biological realities:
- Expression is more variable at tissue boundaries
- Technical artifacts (dropouts) are more common at edges
- There's a gradual transition between different cell types/regions

## Advantages Over Alternative Approaches

### 1. Compared to Simple Poisson Models

Traditional simulations often use Poisson distributions, which:
- Constrain variance to equal the mean
- Cannot model overdispersion seen in real data
- Miss the heterogeneous noise structure across tissue

### 2. Compared to Uniform Dropout Models

Many simulations use uniform random dropout, while our approach:
- Models spatially-varying dropout rates
- Creates realistic patterns of missing data
- Reflects the true technical biases in spatial transcriptomics

### 3. Compared to Discrete Spatial Correlations

Simple approaches might model correlation only within predefined regions, while our GP approach:
- Creates continuous correlation patterns
- Allows smooth transitions between regions
- Better models the true biological continuity of tissues

## Conclusion

This simulation and analysis framework provides a sophisticated toolkit for generating realistic spatial transcriptomics data that captures key statistical properties observed in real experiments. By incorporating spatially-varying dispersion, dropout, and gene-specific variation patterns, it produces data that serves as a robust benchmark for developing and testing spatial transcriptomics analysis methods.

The framework is highly parameterized to allow for customization to different tissue types, technological platforms, and biological questions, making it a versatile tool for methodological research in spatial genomics.