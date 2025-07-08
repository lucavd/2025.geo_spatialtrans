# Cluster Improvement Strategy for Spatial Transcriptomics Simulation

## Problem Statement

The current clustering mechanism uses k-means with spatial weighting, which inherently produces circular/spherical clusters. This is biologically unrealistic because real tissues exhibit:
- Elongated structures (vessels, ducts, muscle fibers)
- Layered organizations (epithelium, cortex)
- Irregular boundaries between cell types
- Non-convex regions and infiltrative patterns

## Proposed Solution: Three-Tier Pattern Generation System

### Tier 1: Geometric Primitive Generators (Simple & Fast)

These provide direct control over tissue architecture through mathematical functions.

#### 1.1 Layered/Laminar Structures

```r
# File: R/functions/04b_tissue_patterns_geometric.R

#' Generate layered tissue structure
#' @param grid_coords Data frame with x, y coordinates
#' @param n_layers Number of tissue layers
#' @param orientation "horizontal", "vertical", or angle in degrees
#' @param layer_thickness Vector of relative thickness for each layer
#' @param boundary_noise Amount of noise at layer boundaries (0-1)
#' @return Vector of cluster assignments
generate_layered_structure <- function(grid_coords, 
                                     n_layers = 3,
                                     orientation = "horizontal",
                                     layer_thickness = NULL,
                                     boundary_noise = 0.1) {
  
  # Default equal thickness
  if (is.null(layer_thickness)) {
    layer_thickness <- rep(1/n_layers, n_layers)
  }
  
  # Convert orientation to angle
  angle <- switch(orientation,
                  "horizontal" = 0,
                  "vertical" = 90,
                  as.numeric(orientation))
  
  # Rotate coordinates
  angle_rad <- angle * pi / 180
  x_rot <- grid_coords$x * cos(angle_rad) + grid_coords$y * sin(angle_rad)
  
  # Add noise to boundaries using Perlin noise
  if (boundary_noise > 0) {
    # Simple sine-based noise for now
    noise <- sin(grid_coords$y * 0.1) * boundary_noise * diff(range(x_rot))
    x_rot <- x_rot + noise
  }
  
  # Normalize to 0-1
  x_norm <- (x_rot - min(x_rot)) / (max(x_rot) - min(x_rot))
  
  # Assign to layers based on cumulative thickness
  cumulative_thickness <- cumsum(layer_thickness)
  clusters <- findInterval(x_norm, c(0, cumulative_thickness[1:(n_layers-1)], 1))
  
  return(clusters)
}

# Example usage:
# cortical_layers <- generate_layered_structure(
#   grid_coords, 
#   n_layers = 6,  # Cortical layers I-VI
#   layer_thickness = c(0.1, 0.15, 0.25, 0.2, 0.2, 0.1),
#   boundary_noise = 0.05
# )
```

#### 1.2 Branching/Vascular Structures

```r
#' Generate branching structure (e.g., blood vessels, ducts)
#' @param grid_coords Data frame with x, y coordinates
#' @param n_main_branches Number of main branches
#' @param branch_width Width of branches in spatial units
#' @param branch_probability Probability of secondary branching
#' @param seed Random seed for reproducibility
#' @return Vector of cluster assignments (0 = background, 1+ = branches)
generate_branching_structure <- function(grid_coords,
                                       n_main_branches = 3,
                                       branch_width = 5,
                                       branch_probability = 0.3,
                                       seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize cluster assignments
  clusters <- rep(0, nrow(grid_coords))
  
  # Generate main branches using random walk
  x_range <- range(grid_coords$x)
  y_range <- range(grid_coords$y)
  
  for (i in 1:n_main_branches) {
    # Start point
    start_x <- runif(1, x_range[1], x_range[2])
    start_y <- y_range[1]
    
    # Random walk parameters
    current_x <- start_x
    current_y <- start_y
    step_size <- diff(y_range) / 50
    
    # Create branch path
    while (current_y < y_range[2]) {
      # Find cells within branch width
      distances <- sqrt((grid_coords$x - current_x)^2 + 
                       (grid_coords$y - current_y)^2)
      clusters[distances < branch_width] <- i
      
      # Random walk step
      current_x <- current_x + rnorm(1, 0, step_size)
      current_y <- current_y + step_size
      
      # Secondary branching
      if (runif(1) < branch_probability) {
        # Branch off at angle
        branch_angle <- sample(c(-1, 1), 1) * runif(1, 20, 70) * pi/180
        branch_length <- runif(1, 10, 30) * step_size
        
        # Trace branch
        for (j in 1:10) {
          bx <- current_x + j * branch_length/10 * cos(branch_angle)
          by <- current_y + j * branch_length/10 * sin(branch_angle)
          distances <- sqrt((grid_coords$x - bx)^2 + (grid_coords$y - by)^2)
          clusters[distances < branch_width * 0.7] <- i
        }
      }
    }
  }
  
  return(clusters)
}
```

#### 1.3 Infiltrative Patterns

```r
#' Generate infiltrative pattern (e.g., immune infiltration, tumor invasion)
#' @param base_clusters Existing cluster assignments to infiltrate
#' @param infiltration_source Cluster ID that will infiltrate
#' @param infiltration_target Cluster IDs that can be infiltrated
#' @param infiltration_prob Base probability of infiltration
#' @param spatial_correlation Spatial correlation of infiltration (0-1)
#' @return Modified cluster assignments
generate_infiltrative_pattern <- function(base_clusters,
                                        grid_coords,
                                        infiltration_source,
                                        infiltration_target,
                                        infiltration_prob = 0.1,
                                        spatial_correlation = 0.5) {
  
  # Create smooth probability field using Gaussian blur
  # For now, use distance-based probability
  source_cells <- which(base_clusters == infiltration_source)
  target_cells <- which(base_clusters %in% infiltration_target)
  
  if (length(source_cells) == 0 || length(target_cells) == 0) {
    return(base_clusters)
  }
  
  # Calculate distance to nearest source cell
  distances_to_source <- sapply(target_cells, function(i) {
    min(sqrt((grid_coords$x[i] - grid_coords$x[source_cells])^2 +
             (grid_coords$y[i] - grid_coords$y[source_cells])^2))
  })
  
  # Convert distance to probability with spatial decay
  max_dist <- quantile(distances_to_source, 0.9)
  infiltration_probs <- infiltration_prob * exp(-distances_to_source / (max_dist * spatial_correlation))
  
  # Stochastic infiltration
  infiltrated <- runif(length(target_cells)) < infiltration_probs
  base_clusters[target_cells[infiltrated]] <- infiltration_source
  
  return(base_clusters)
}
```

### Tier 2: Topology-Aware Clustering (Data-Driven)

These methods discover non-convex clusters from data.

#### 2.1 DBSCAN Implementation

```r
# File: R/functions/04c_tissue_patterns_topology.R

#' DBSCAN clustering for irregular tissue structures
#' @param grid_coords Data frame with x, y coordinates
#' @param eps Maximum distance for neighborhood (biological: cell-cell interaction distance)
#' @param min_samples Minimum cells to form a cluster (biological: minimum niche size)
#' @return Vector of cluster assignments (-1 = noise/isolated cells)
dbscan_clustering <- function(grid_coords, eps = 10, min_samples = 5) {
  
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Please install it.")
  }
  
  # Run DBSCAN
  coords_matrix <- as.matrix(grid_coords[, c("x", "y")])
  db_result <- dbscan::dbscan(coords_matrix, eps = eps, minPts = min_samples)
  
  return(db_result$cluster)
}

# Example with biological interpretation:
# tissue_clusters <- dbscan_clustering(
#   grid_coords,
#   eps = 15,  # ~3 cell diameters for local interaction
#   min_samples = 10  # Minimum 10 cells to form stable micro-niche
# )
```

#### 2.2 Graph-Based Clustering

```r
#' Graph-based clustering using spatial neighborhoods
#' @param grid_coords Data frame with x, y coordinates  
#' @param k_neighbors Number of nearest neighbors for graph construction
#' @param resolution Resolution parameter for community detection (higher = more clusters)
#' @return Vector of cluster assignments
graph_based_clustering <- function(grid_coords, k_neighbors = 15, resolution = 1.0) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Please install it.")
  }
  
  # Build k-NN graph
  coords_matrix <- as.matrix(grid_coords[, c("x", "y")])
  distances <- as.matrix(dist(coords_matrix))
  
  # Create adjacency matrix (k nearest neighbors)
  adj_matrix <- matrix(0, nrow = nrow(coords_matrix), ncol = nrow(coords_matrix))
  for (i in 1:nrow(coords_matrix)) {
    neighbors <- order(distances[i,])[2:(k_neighbors+1)]  # Exclude self
    adj_matrix[i, neighbors] <- 1
    adj_matrix[neighbors, i] <- 1  # Symmetric
  }
  
  # Create igraph object
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # Run Louvain clustering
  communities <- igraph::cluster_louvain(g, resolution = resolution)
  
  return(igraph::membership(communities))
}
```

### Tier 3: Biomimetic Process Simulation (Highest Fidelity)

These simulate actual biological processes.

#### 3.1 Reaction-Diffusion System

```r
# File: R/functions/04d_tissue_patterns_biomimetic.R

#' Generate tissue patterns using reaction-diffusion (Turing patterns)
#' @param grid_size Grid dimensions (assumes square grid)
#' @param pattern_type "spots", "stripes", or "labyrinth"
#' @param n_iterations Number of simulation iterations
#' @param D_ratio Ratio of diffusion constants (controls pattern type)
#' @return Matrix of morphogen concentrations (threshold to get clusters)
reaction_diffusion_pattern <- function(grid_size = 100, 
                                     pattern_type = "spots",
                                     n_iterations = 5000,
                                     D_ratio = NULL) {
  
  # Gray-Scott model parameters
  params <- switch(pattern_type,
    "spots" = list(F = 0.030, k = 0.062, Du = 0.2, Dv = 0.1),
    "stripes" = list(F = 0.035, k = 0.065, Du = 0.2, Dv = 0.1),
    "labyrinth" = list(F = 0.040, k = 0.060, Du = 0.2, Dv = 0.1)
  )
  
  if (!is.null(D_ratio)) {
    params$Dv <- params$Du / D_ratio
  }
  
  # Initialize concentrations
  U <- matrix(1, grid_size, grid_size)
  V <- matrix(0, grid_size, grid_size)
  
  # Add random perturbation in center
  center <- grid_size / 2
  size <- grid_size / 10
  perturbation_area <- which(
    abs(row(U) - center) < size & abs(col(U) - center) < size,
    arr.ind = TRUE
  )
  V[perturbation_area] <- runif(nrow(perturbation_area), 0.2, 0.3)
  
  # Laplacian operator (with periodic boundary conditions)
  laplacian <- function(M) {
    L <- matrix(0, nrow(M), ncol(M))
    # Shift operations for 5-point stencil
    L <- L + cbind(M[,-1], M[,1]) + cbind(M[,ncol(M)], M[,-ncol(M)])
    L <- L + rbind(M[-1,], M[1,]) + rbind(M[nrow(M),], M[-nrow(M),])
    L <- L - 4 * M
    return(L)
  }
  
  # Time evolution
  dt <- 1.0
  for (iter in 1:n_iterations) {
    # Reaction-diffusion equations
    Lu <- laplacian(U)
    Lv <- laplacian(V)
    
    dU <- params$Du * Lu - U * V^2 + params$F * (1 - U)
    dV <- params$Dv * Lv + U * V^2 - (params$F + params$k) * V
    
    U <- U + dt * dU
    V <- V + dt * dV
    
    # Ensure stability
    U[U < 0] <- 0
    V[V < 0] <- 0
  }
  
  return(V)  # Return activator concentration
}

# Convert to clusters
#' @param morphogen_field Output from reaction_diffusion_pattern
#' @param n_clusters Number of clusters to create by thresholding
pattern_to_clusters <- function(morphogen_field, grid_coords, n_clusters = 3) {
  # Interpolate morphogen field to actual grid coordinates
  # For simplicity, assume grid_coords are on regular grid
  
  # Threshold into clusters
  thresholds <- quantile(as.vector(morphogen_field), 
                        probs = seq(0, 1, length.out = n_clusters + 1))
  
  # Map to grid coordinates (simplified - assumes regular grid)
  clusters <- cut(as.vector(morphogen_field), 
                 breaks = thresholds, 
                 labels = FALSE,
                 include.lowest = TRUE)
  
  return(clusters)
}
```

### Integration Strategy

#### Modify Main Simulation Function

```r
# In 07c_simulation_pipeline.R, add:

simulate_spatial_transcriptomics <- function(
  image_path = NULL,
  pattern_method = c("spatial_kmeans", "kmeans++", "slic", 
                     "layered", "branching", "dbscan", "graph",
                     "reaction_diffusion"),
  pattern_params = list(),
  n_clusters = 8,
  ...) {
  
  pattern_method <- match.arg(pattern_method)
  
  # Route to appropriate pattern generator
  if (!is.null(image_path)) {
    # Image-based methods
    clusters <- switch(pattern_method,
      "spatial_kmeans" = spatial_kmeans(...),
      "kmeans++" = kmeans_plusplus(...),
      "slic" = slic_clustering(...),
      stop("Image-based method not recognized")
    )
  } else {
    # Synthetic pattern methods
    clusters <- switch(pattern_method,
      "layered" = do.call(generate_layered_structure, pattern_params),
      "branching" = do.call(generate_branching_structure, pattern_params),
      "dbscan" = do.call(dbscan_clustering, pattern_params),
      "graph" = do.call(graph_based_clustering, pattern_params),
      "reaction_diffusion" = {
        morph <- do.call(reaction_diffusion_pattern, pattern_params)
        pattern_to_clusters(morph, grid_coords, n_clusters)
      },
      stop("Pattern method not recognized")
    )
  }
  
  # Continue with rest of pipeline...
}
```

## ✅ IMPLEMENTED: Pipeline Consecutiva DBSCAN + Graph Clustering

### Status Implementazione (Gennaio 2025)

**✅ COMPLETATO - Fase 1**: Pipeline consecutiva DBSCAN + Graph clustering
- **File implementato**: `R/functions/04_clustering.R`
- **Funzione principale**: `dbscan_graph_pipeline()`
- **Integrazione**: `run_full_size_optimized.R` usa `clustering_method = "dbscan_graph"`

### Strategia Consecutiva Implementata

```r
# Pipeline DBSCAN + Graph - Implementata
dbscan_graph_pipeline <- function(img_df_thresh, k_cell_types, random_seed = 123) {
  # Fase 1: DBSCAN per identificare regioni dense
  dbscan_result <- dbscan_clustering(
    img_df_thresh, 
    k_cell_types, 
    random_seed = random_seed,
    eps_factor = 1.4,        # Controllo biologico
    min_samples = 4          # Minima nicchia cellulare
  )
  
  # Fase 2: Graph clustering per raffinare
  final_result <- graph_refine_clustering(
    dbscan_result,
    img_df_thresh, 
    k_cell_types, 
    random_seed = random_seed,
    k_neighbors = 10,        # Vicinato locale
    resolution = 1.0         # Granularità Louvain
  )
  
  return(final_result)
}
```

### Parametri Biologicamente Interpretabili

1. **DBSCAN (Fase 1)**:
   - `eps_factor = 1.4`: Dimensione regioni dense (40% sopra mediana k-NN)
   - `min_samples = 4`: Numero minimo di celle per formare un cluster

2. **Graph Clustering (Fase 2)**:
   - `k_neighbors = 10`: Dimensione microambiente locale
   - `resolution = 1.0`: Granularità clustering Louvain

### Vantaggi Ottenuti

✅ **Forme irregolari**: DBSCAN trova vasi, ramificazioni, infiltrazioni  
✅ **Raffinamento topologico**: Graph clustering perfeziona strutture locali  
✅ **Controllo matematico**: Parametri biologicamente interpretabili  
✅ **Fallback robusti**: Se un metodo fallisce, usa spatial_kmeans  
✅ **Drop-in replacement**: Stessa interfaccia di k-means

### Come Usare

```r
# Nel file di simulazione principale
clust <- cluster_image(
  img_df_thresh = img_dat$img_df_thresh,
  k_cell_types  = cfg$k_cell_types,
  random_seed   = cfg$random_seed,
  clustering_method = "dbscan_graph"  # Nuovo metodo
)
```

## Implementation Priority - AGGIORNATA

1. **✅ Phase 1 (COMPLETATA)**: 
   - ✅ Implementata pipeline consecutiva DBSCAN + Graph clustering
   - ✅ Integrata nel sistema di simulazione principale
   - ✅ Parametri biologicamente interpretabili

2. **Phase 2 (Prossimi sviluppi)**:
   - Validazione biologica dei risultati vs k-means
   - Implementazione Tier 1 (layered, branching patterns) 
   - Ottimizzazione performance per simulazioni large-scale

3. **Phase 3 (Future)**:
   - Full reaction-diffusion implementation
   - Agent-based modeling framework
   - Pattern combination system

## Testing Strategy

```r
# Test file: tests/testthat/test-tissue_patterns.R

test_that("Layered structures have correct properties", {
  coords <- expand.grid(x = 1:100, y = 1:100)
  clusters <- generate_layered_structure(coords, n_layers = 3)
  
  expect_equal(length(unique(clusters)), 3)
  expect_true(all(clusters %in% 1:3))
  
  # Test that layers are contiguous
  # ...
})

test_that("DBSCAN finds non-convex clusters", {
  # Create C-shaped point cloud
  theta <- seq(0, 1.5*pi, length.out = 100)
  coords <- data.frame(
    x = c(cos(theta), 0.5*cos(theta)),
    y = c(sin(theta), 0.5*sin(theta))
  )
  
  clusters <- dbscan_clustering(coords, eps = 0.3)
  
  # Should find 1 C-shaped cluster, not 2 circular ones
  expect_equal(length(unique(clusters[clusters > 0])), 1)
})
```

## Biological Validation Criteria

For each pattern type, validate:

1. **Morphological metrics**: Shape index, solidity, convexity
2. **Spatial statistics**: Moran's I, Geary's C for spatial autocorrelation  
3. **Interface characteristics**: Roughness, fractal dimension
4. **Biological plausibility**: Cell density, cluster size distributions

## Notes for Future Implementation

- Consider caching expensive computations (e.g., reaction-diffusion)
- Add visualization functions for each pattern type
- Create pattern "presets" for common tissue types
- Document biological interpretation of all parameters
- Consider GPU acceleration for reaction-diffusion