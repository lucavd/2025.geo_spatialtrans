# Advanced Spatial Transcriptomics Simulation Framework: Theoretical Foundations and Implementation

## 1. Introduction to Spatial Transcriptomics

### 1.1 Biological and Technical Context

Spatial transcriptomics represents a revolutionary class of technologies that measure gene expression while preserving the spatial context of cells within tissues. Unlike conventional RNA sequencing methodologies (bulk or single-cell), spatial transcriptomics technologies—such as 10x Visium, Slide-seq, MERFISH, and Spatial Transcriptomics (ST)—provide crucial information about the spatial arrangement of gene expression patterns. This spatial information is essential for understanding complex tissue architecture, cell-cell interactions, developmental processes, and disease progression mechanisms that cannot be fully elucidated through spatially-agnostic approaches.

The technological landscape of spatial transcriptomics encompasses several methodological categories:
- **In situ hybridization-based methods** (e.g., MERFISH, seqFISH+): Directly visualize and quantify transcripts within intact tissues using oligonucleotide probes
- **In situ sequencing methods** (e.g., STARmap): Sequence RNA directly within tissue sections using localized amplification
- **Spatial barcoding methods** (e.g., 10x Visium, Slide-seq): Capture and sequence mRNA using spatially barcoded capture spots or beads
- **Computational integration methods**: Combine conventional single-cell RNA-seq with spatial reference data to infer spatial gene expression patterns

Each of these technologies offers unique trade-offs between spatial resolution, gene throughput, sensitivity, and tissue compatibility. The appropriate technology depends on the specific biological questions being addressed.

### 1.2 The Need for Sophisticated Simulation Frameworks

The development of analytical methods for spatial transcriptomics data presents significant computational challenges that require rigorous benchmarking using data with known ground truth. However, experimental generation of ground truth spatial transcriptomics data is prohibitively expensive, time-consuming, and often practically infeasible due to technical limitations in simultaneously validating expression patterns for thousands of genes.

Computational simulations offer an alternative approach by generating synthetic data with known properties that can be used to benchmark analytical methods. However, existing simulation approaches have several limitations:

- **Oversimplified statistical properties**: Many simulations use simple Poisson models that fail to capture the complex variance structure of real transcriptomic data
- **Biologically unrealistic spatial patterns**: Simulations often use simplistic geometric patterns rather than mimicking the complex tissue architectures found in vivo
- **Homogeneous technical artifacts**: Most simulations apply uniform dropout rates and noise models, whereas real data exhibits spatially variable technical artifacts
- **Discrete boundaries**: Many simulations create sharp boundaries between cell types, whereas real tissues often have gradual transitions and interface regions

Our framework addresses these limitations by implementing sophisticated statistical models informed by the biological and technical characteristics of real spatial transcriptomics data, enabling more accurate evaluation of analytical methods.

## 2. Framework Architecture and Components

This repository provides a comprehensive set of tools for simulating and analyzing spatial transcriptomics data with various statistical properties. The framework is organized into a modular architecture with the following components:

### 2.1 Project Structure

```
R/
├── analyze_and_compare_clusters.R    # Cluster evaluation framework
├── image_simulations.R               # Core simulation engine
├── image_simulations_partial_images.R # Visualization utilities
└── testing/                          # Testing modules for simulation validation
```

The design philosophy of the framework emphasizes:
- **Modularity**: Each component focuses on a specific task within the pipeline
- **Flexibility**: Parameters can be customized to simulate different biological scenarios
- **Theoretical foundation**: Statistical models are selected based on empirical observations from real data
- **Benchmarking capabilities**: Ground truth knowledge enables quantitative evaluation of analytical methods

### 2.2 Key Components

#### 2.2.1 Simulation Core

The heart of the framework is the `simulate_spatial_transcriptomics()` function in `image_simulations.R`, which generates realistic spatial transcriptomics data from image inputs. This function implements sophisticated statistical models to create synthetic data with complex properties reflecting real-world biological and technical characteristics.

#### 2.2.2 Evaluation Framework

The `analyze_and_compare_clusters.R` script provides a comprehensive framework for evaluating clustering performance using:
- Multiple clustering algorithms (Seurat and HDBSCAN)
- Quantitative performance metrics
- Spatial visualization of results

## 3. Visium HD Mode: High-Resolution Simulation

The framework now includes dedicated support for simulating Visium HD data, a high-resolution spatial transcriptomics technology from 10x Genomics that provides 2μm x 2μm resolution.

### 3.1 Grid-Based Spatial Structure

Visium HD uses a regular grid of spots rather than the scattered pattern of original Visium. Our framework implements this structure through:

```r
# Grid mode with high resolution (2μm)
simulate_spatial_transcriptomics(
  image_path = "images/colon.png",
  grid_mode = TRUE,              # Enables grid-based simulation
  grid_resolution = 2,           # 2μm x 2μm spots like Visium HD
  grid_spacing = 0,              # No gap between adjacent spots
  # Other parameters...
)
```

This approach:
- Creates a regular grid of adjacent spots following the exact Visium HD geometry
- Maps image regions to grid locations deterministically rather than randomly
- Matches the tissue structure precision of high-resolution technologies
- Enables accurate modeling of spatial autocorrelation at microscopic scales

### 3.2 Advanced Spatial Correlation Models

The framework now offers two complementary methods for spatial correlation modeling:

#### 3.2.1 Gaussian Random Fields (GRF)

GRF implements a continuous spatial correlation model through Gaussian processes with exponential covariance:

```r
correlation_method = "grf",
spatial_params = list(
  spatial_noise_intensity = 1.0,  # Magnitude of spatial effect
  spatial_range = 30,            # Correlation length scale (μm)
  random_noise_sd = 0.2          # Cell-specific random variation
)
```

This model is particularly suitable for:
- Smooth gradient patterns across tissue regions
- Diffusion-like processes (morphogen gradients, secreted signals)
- Continuous biological processes that vary gradually in space

#### 3.2.2 Conditional Autoregressive Models (CAR)

The CAR model is newly implemented for grid-based simulations:

```r
correlation_method = "car",
spatial_params = list(
  spatial_noise_intensity = 1.2,  # Controls variance
  spatial_range = 20             # Controls neighborhood size
)
```

This model is ideal for:
- Capturing local dependencies between adjacent spots
- Modeling interacting cellular neighborhoods
- Simulating lattice-based processes common in grid-structured data

### 3.3 Gradient-Based Region Transitions

A key advancement in the framework is the implementation of gradient-based transitions between tissue regions:

```r
spatial_params = list(
  gradient_regions = TRUE,      # Enable gradient transitions
  gradient_width = 5            # Width of gradient zone (in grid units)
)
```

This feature:
- Calculates distance to region boundaries for each grid point
- Applies gradual phenotypic transitions at region interfaces
- Mixes expression profiles of adjacent regions based on distance
- Creates realistic cell state transitions rather than artificial sharp boundaries

The implementation automatically identifies boundary regions and creates distance maps for smooth transitions, consistent with observations in real tissue interfaces.

### 3.4 Library Size and Dropout Modeling

The framework now incorporates sophisticated modeling of library size (sequencing depth) and dropout effects:

#### 3.4.1 Spatially Varying Library Size

```r
library_size_params = list(
  mean_library_size = 10000,     # Mean UMI count per spot
  library_size_cv = 0.3,         # Coefficient of variation
  spatial_effect_on_library = 0.5 # Spatial correlation in library size
)
```

This accounts for:
- Log-normal distribution of library sizes observed in real data
- Spatial correlation in sequencing depth due to tissue properties
- Impact of library size on expression level and zero counts

#### 3.4.2 Expression-Dependent Dropout

Dropout (zero counts) in real data depends strongly on expression level. The framework now models this relationship:

```r
dropout_params = list(
  dropout_range = c(0.1, 0.4),    # Base dropout rates (min, max)
  expression_dependent_dropout = TRUE,  # Enable expression-dependent dropout
  dropout_curve_midpoint = 0.5,   # Expression level at 50% dropout probability
  dropout_curve_steepness = 5     # Steepness of logistic dropout curve
)
```

This creates realistic patterns where:
- Lowly expressed genes have higher dropout probability
- Relationship follows logistic function, matching empirical observations
- Dropout combines spatial effects with expression-level dependency
- Dropout rate increases at tissue borders, mimicking edge artifacts

### 3.5 Enhanced Validation Metrics

The framework now calculates spatial autocorrelation metrics for validation:

```r
# Examine spatial autocorrelation in simulated data
result$spatial_autocorrelation$moran_i  # Moran's I values for selected genes
```

This provides:
- Quantitative measurement of spatial structure
- Comparability with real Visium HD data
- Validation of the spatial correlation models
- Benchmark for spatial analysis methods

### 3.6 Performance Optimization

The new implementation includes significant computational optimizations:

- **Vectorized border detection**: Efficiently identifies region boundaries using matrix operations
- **Optimized distance calculations**: Uses fast distance matrix computation with spatial indexing
- **Vectorized expression generation**: Replaces loops with matrix operations for expression calculation
- **Efficient gradient application**: Applies gradients to all border points simultaneously
- **Block processing**: Handles large datasets through block-wise calculations

These optimizations enable simulating large Visium HD datasets (100,000+ spots) with reasonable computational resources.

## 4. Theoretical Foundations of the Simulation Framework

### 4.1 Count Distribution Models

The core statistical challenge in transcriptomics simulation is accurately modeling the distribution of gene expression counts. Our framework implements multiple distributions based on empirical observations from real data:

#### 4.1.1 Negative Binomial Distribution

The primary distribution used in our simulation is the Negative Binomial (NB) distribution, which naturally models the overdispersion (variance > mean) consistently observed in transcriptomic data. The probability mass function is:

$$P(X=k) = \binom{k+r-1}{k} p^r (1-p)^k$$

Where:
- $r$ is the dispersion parameter (called `size` in R)
- $p$ is related to the mean $\mu$ by $p = \frac{r}{r+\mu}$

The variance has the fundamental relationship:

$$\text{Var}(X) = \mu + \frac{\mu^2}{r}$$

This framework provides critical flexibility, as the variance structure can be modulated through the dispersion parameter $r$:
- When $r \to \infty$, the NB approaches a Poisson distribution (variance = mean)
- When $r$ is small, the variance can be much larger than the mean, matching empirical observations

The biological justification for using the NB distribution includes:
- **Stochastic bursting**: Transcription occurs in episodic bursts rather than continuous production
- **Regulatory network variability**: Cell-to-cell differences in regulatory network states
- **Microenvironmental heterogeneity**: Local variations in the cellular microenvironment even within the same nominal cell type

#### 4.1.2 Sub-Poisson Model for Constitutive Genes

For a subset of genes (approximately 10%), our simulation implements a Binomial model with high success probability (p=0.9). This model captures the behavior of constitutively expressed "housekeeping" genes that show remarkably stable expression with lower-than-Poisson variance (variance-to-mean ratio < 1).

The probability mass function for the Binomial distribution is:

$$P(X=k) = \binom{n}{k} p^k (1-p)^{n-k}$$

The variance is given by:

$$\text{Var}(X) = np(1-p)$$

When $p$ is high (e.g., 0.9), the variance becomes lower than the mean, creating the sub-Poisson effect observed in genes under tight homeostatic regulation. The parameters are selected to maintain the desired mean while reducing variance:

```r
p <- 0.9
n_trial <- round(exp(mu_vals)/(1-p))
expression_data[, g] <- rbinom(N, n_trial, p)
```

The biological rationale for implementing sub-Poisson models includes:
- **Homeostatic feedback regulation**: Tight control of essential genes through negative feedback mechanisms
- **Redundant regulatory mechanisms**: Multiple parallel regulatory pathways ensuring stable expression
- **High-frequency transcriptional initiation**: Consistent production with reduced bursting behavior

### 4.2 Spatial Correlation Models

Real tissues exhibit complex patterns of spatial correlation in gene expression due to intercellular communication, developmental gradients, and tissue organization. Our framework implements sophisticated spatial correlation models to capture these biological realities.

#### 4.2.1 Gaussian Process Implementation

We implement a Gaussian Process (GP) model with an exponential covariance function to create spatially correlated random fields:

$$C(d) = \sigma^2 \exp(-d/\rho)$$

Where:
- $d$ is the Euclidean distance between spatial locations
- $\rho$ is the range parameter controlling the correlation length scale
- $\sigma^2$ is the variance parameter controlling the magnitude of spatial variation

The implementation uses the `gstat` package in R:

```r
gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                beta = 0, model = vgm(psill=spatial_params$spatial_noise_intensity, 
                                     range=spatial_params$spatial_range, 
                                     model="Exp"), 
                nmax=10)
```

Key parameters that can be adjusted include:
- `spatial_noise_intensity`: Controls the magnitude of spatial variation (psill)
- `spatial_range`: Controls the correlation length scale (range)

The biological justification for using a GP with exponential covariance includes:
- **Diffusion physics**: Molecular gradients in tissues often follow exponential decay patterns
- **Hierarchical organization**: Tissues show multi-scale organization with different correlation ranges
- **Continuous transitions**: Expression changes gradually across spatial domains rather than discontinuously

#### 4.2.2 Conditional Autoregressive Model

For grid-based data like Visium HD, the CAR model provides an alternative approach to spatial correlation:

$$X_i | X_{-i} \sim \mathcal{N}\left(\alpha + \rho \sum_{j \in N_i} w_{ij}(X_j - \alpha), \tau^2\right)$$

Where:
- $X_i$ is the value at location $i$
- $X_{-i}$ represents all values except at location $i$
- $N_i$ is the neighborhood of location $i$
- $w_{ij}$ are spatial weights
- $\rho$ controls the strength of spatial dependence
- $\tau^2$ is the conditional variance

The CAR implementation is particularly suited for regular grid data with well-defined neighborhood structures, making it ideal for high-resolution spatial transcriptomics simulations.

#### 4.2.3 Gradient and Hybrid Cell Implementation

Our framework implements two complementary approaches to modeling transitions between spatial domains:

1. **Gradient-based transitions**: Uses distance to boundary to create smooth transitions:

```r
gradient_weight <- (1 - cell_df$boundary_dist[i])^2
base_expr[i] <- base_expr[i] * (1 - gradient_weight) + other_expr * gradient_weight
```

2. **Hybrid cell approach**: Explicitly models cells at boundaries as mixtures:

```r
hybrid_effect <- hybrid_matrix %*% all_cluster_expr
hybrid_weight <- rowSums(hybrid_matrix)
base_expr[hybrid_cells] <- base_expr[hybrid_cells] * (1 - hybrid_weight[hybrid_cells]) + 
                         hybrid_effect[hybrid_cells]
```

The biological rationale includes:
- **Transitional cell states**: Cells at interfaces often exhibit intermediate phenotypes
- **Cell-cell communication**: Signaling between adjacent cells can induce partial phenotypic shifts
- **Plasticity gradients**: Cells may show varying degrees of commitment to particular lineages

### 4.3 Technical Artifact Models

Real spatial transcriptomics data contains various technical artifacts that can confound analysis. Our framework explicitly models these artifacts with spatial dependence to create realistic challenges for analytical methods.

#### 4.3.1 Spatially-varying and Expression-dependent Dropout

Dropout (false zeros) in spatial transcriptomics combines spatial effects with expression-level dependency:

```r
# Expression-dependent component
norm_expr <- scale01_vec(expression_data[, g])
dropout_prob_expr <- 1 / (1 + exp((norm_expr - dropout_params$dropout_curve_midpoint) * 
                                 dropout_params$dropout_curve_steepness))

# Combine with spatial component
dropout_prob <- 0.7 * dropout_prob_expr + 0.3 * base_dropout
```

This sophisticated model captures the dual nature of dropout in real data:
- **Expression dependence**: Lower expressed genes have higher dropout probability
- **Spatial effects**: Tissue edges and processing artifacts create spatial patterns in dropout
- **Combined effect**: The final model combines both factors with appropriate weighting

#### 4.3.2 Spatially-varying Dispersion

The variability in gene expression (dispersion) also exhibits spatial dependence in real tissues. Our framework implements spatially-varying dispersion using the distance-based metric:

```r
# Near boundaries: higher variability (lower dispersion parameter)
dispersion_param <- dropout_params$dispersion_range[2] + 
  cell_df$boundary_dist * (dropout_params$dispersion_range[1] - dropout_params$dispersion_range[2])
```

The biological justification includes:
- **Border instability**: Cells at tissue interfaces show higher transcriptional variability
- **Stress response heterogeneity**: Variability in stress responses at tissue edges
- **Identity ambiguity**: Cells in transitional zones show less stable gene expression patterns

### 4.4 Library Size Modeling

The framework models library size (sequencing depth) with both global and spatial components:

```r
# Log-normal distribution for global variation
library_size <- rlnorm(N, meanlog = log_mean, sdlog = log_sd)

# Spatial effect using Gaussian Process
lib_effect <- library_size_params$spatial_effect_on_library * lib_noise
library_size <- library_size * exp(lib_effect)

# Apply to expression counts
scaled_counts <- raw_counts * (library_size / mean(library_size))
```

This accounts for:
- **Log-normal global distribution**: Matching empirical observations in real data
- **Spatial correlation**: Areas with better RNA preservation or higher cell density show correlated sequencing depth
- **Scaling effect**: Library size acts as a scaling factor on observed counts

## 5. Customizable Difficulty Levels

A key feature of our framework is the ability to simulate data with varying levels of analytical challenge through pre-defined difficulty tiers.

### 5.1 Difficulty Parameterization

The framework provides three difficulty levels with comprehensive parameter adjustments:

```r
# Example: dataset with custom difficulty level
simulate_spatial_transcriptomics(
  image_path = "images/colon.png",
  difficulty_level = "hard",  # One of "easy", "medium", "hard"
  grid_mode = TRUE,           # Use Visium HD grid mode
  n_genes = 100
)
```

#### 5.1.1 Easy Difficulty

The "easy" setting creates data with clear cell type boundaries and strong marker genes:
- 10 marker genes per cell type with +2.0 log-fold expression
- No marker overlap between cell types
- Low spatial noise (0.5 intensity, 50 range)
- Minimal technical artifacts (10-30% dropout)
- Low biological variation (dispersion range 3.0-1.5)

This setting is appropriate for benchmarking basic analytical methods or educational purposes.

#### 5.1.2 Medium Difficulty

The "medium" setting introduces moderate challenges:
- 7 marker genes per cell type with +1.2 log-fold expression
- 20% marker overlap between adjacent types
- Moderate spatial noise (1.0 intensity, 30 range)
- Moderate technical artifacts (20-50% dropout)
- Substantial biological variation (dispersion range 2.0-1.0)

This setting approximates high-quality real-world datasets from technologies like 10x Visium.

#### 5.1.3 Hard Difficulty

The "hard" setting creates extremely challenging data that reflects difficult real-world scenarios:
- Only 5 marker genes per cell type with +0.8 log-fold expression
- 40% marker overlap between adjacent types
- Strong spatial noise (1.5 intensity, 15 range)
- Extreme technical artifacts (30-70% dropout)
- Very high biological variation (dispersion range 1.5-0.8)
- Additional random noise component (SD 0.4)

This setting mimics challenging datasets from tissues with subtle biological differences, significant technical noise, or technologies with high dropout rates.

## 6. Advanced Customization

For users with specific research needs, the framework provides extensive customization options beyond the pre-defined difficulty levels.

### 6.1 Full Parameter Customization

The function accepts detailed parameter lists for fine-grained control:

```r
# Example with fully customized parameters for Visium HD simulation
simulate_spatial_transcriptomics(
  # Core parameters
  image_path = "images/tissue.png",
  output_path = "data/custom_simulation.rds",
  output_plot = "results/custom_distribution.png",
  grid_mode = TRUE,
  grid_resolution = 2,
  n_genes = 150,
  k_cell_types = 7,
  correlation_method = "grf",
  
  # Marker gene parameters
  marker_params = list(
    marker_genes_per_type = 8,    # Number of marker genes per type
    marker_expression_fold = 1.0, # Log-fold change for markers
    marker_overlap_fold = 0.3     # Overlap between adjacent types
  ),
  
  # Spatial correlation parameters
  spatial_params = list(
    spatial_noise_intensity = 1.2, # Intensity of spatial effect
    spatial_range = 25,            # Correlation range
    random_noise_sd = 0.3,         # Cell-specific random noise
    gradient_regions = TRUE,       # Enable gradient transitions
    gradient_width = 8             # Width of transition regions
  ),
  
  # Modeling library size and technical artifacts
  library_size_params = list(
    mean_library_size = 8000,       # Average UMI count per spot
    library_size_cv = 0.4,          # Coefficient of variation
    spatial_effect_on_library = 0.6 # Spatial correlation in library size
  ),
  
  # Dropout and technical artifacts
  dropout_params = list(
    dropout_range = c(0.15, 0.5),    # Base dropout probabilities
    dispersion_range = c(2.0, 0.9),  # Range of NB dispersion values
    expression_dependent_dropout = TRUE,  # Enable expression-dependent dropout
    dropout_curve_midpoint = 0.4,    # Expression level at 50% dropout
    dropout_curve_steepness = 6      # Steepness of dropout curve
  ),
  
  # Cell-specific variation
  cell_specific_params = list(
    cell_specific_noise_sd = 0.25  # SD of cell-specific random effect
  )
)
```

This extensive parameterization allows researchers to:
- Test specific hypotheses about technical artifacts
- Model particular biological scenarios
- Create custom benchmarking challenges
- Simulate data resembling specific technology platforms

### 6.2 Image-Based Customization

The framework uses image inputs to define spatial domains, offering several advantages:

1. **Biological Realism**: Real tissue images can be used to create simulations that match actual tissue architecture
2. **Geometric Complexity**: Images can define complex geometric arrangements impossible to parameterize directly
3. **Multi-domain Structures**: Nested or interdigitated tissue structures can be represented
4. **Custom Challenges**: Artificial images can create specific test cases for algorithm evaluation

The image processing pipeline includes:
- Grayscale conversion for multi-channel images
- Intensity thresholding to define tissue regions
- K-means++ clustering to define spatial domains
- Projection onto regular grid for Visium HD mode

## 7. Validation and Analysis Tools

### 7.1 Spatial Autocorrelation Analysis

The framework calculates Moran's I statistic to quantify spatial autocorrelation:

$$I = \frac{n}{S_0} \frac{\sum_i \sum_j w_{ij} (x_i - \bar{x}) (x_j - \bar{x})}{\sum_i (x_i - \bar{x})^2}$$

Where:
- $n$ is the number of spots/cells
- $w_{ij}$ are spatial weights
- $S_0$ is the sum of all weights
- $x_i$ is the gene expression at location $i$
- $\bar{x}$ is the mean expression

This provides a standardized measure of spatial structure that can be:
- Compared between simulated and real datasets
- Used to validate spatial correlation models
- Applied to benchmark spatial analysis methods

### 7.2 Enhanced Visualization

The framework includes advanced visualization capabilities:

1. **Grid visualization**: Shows spot-level data with appropriate geometry
2. **Library size maps**: Visualizes spatial patterns in sequencing depth
3. **Autocorrelation plots**: Displays spatial structure metrics
4. **Multiple output formats**: Saves visualizations for further analysis

These visualizations provide both qualitative and quantitative assessment of simulation quality and realism.

## 8. Conclusion and Future Directions

This simulation and analysis framework provides a sophisticated toolkit for generating realistic spatial transcriptomics data that captures key statistical properties observed in real experiments. By incorporating spatially-varying dispersion, dropout, and gene-specific variation patterns, it produces data that serves as a robust benchmark for developing and testing spatial transcriptomics analysis methods.

### 8.1 Key Contributions

The framework makes several important contributions to spatial transcriptomics methodology:
1. **Realistic statistical properties** that match real data characteristics
2. **Grid-based simulation for Visium HD** that models high-resolution technologies
3. **Advanced spatial correlation models** using both GRF and CAR approaches
4. **Enhanced gradient-based transitions** with non-linear functions and multi-type cell mixing
5. **Gene co-expression modules** that simulate realistic transcriptional programs
6. **Cell-type specific effects** on library size and dispersion
7. **Expression-dependent dropout modeling** that reflects empirical observations
8. **Library size variation** with both spatial and cell type dependencies
9. **Optimized implementation** using vectorized operations for efficiency
10. **Comprehensive evaluation tools** with spatial autocorrelation metrics

### 8.2 Future Extensions

Having implemented several key biological improvements such as gene co-expression modules and cell-type specific effects, future development of this framework could include:
1. **Dynamic temporal components** to model developmental processes and cellular responses
2. **Ligand-receptor interaction modeling** for realistic cell-cell communication networks
3. **Multi-scale spatial patterns** that capture hierarchical tissue organization
4. **Automated parameter inference** from real Visium HD datasets
5. **Multi-omic integration** for simultaneous simulation of transcriptomic, proteomic, and epigenomic data
6. **Context-aware simulation** that incorporates histological features from input images

### 8.3 Applications

This framework is designed to support:
1. **Method development** for spatial transcriptomics analysis
2. **Benchmarking** of clustering and domain detection algorithms
3. **Educational use** for teaching spatial transcriptomics concepts
4. **Hypothesis testing** for experimental design optimization
5. **Technical artifact correction** method development

The highly parameterized design makes it adaptable to a wide range of research questions, technological platforms, and biological systems.