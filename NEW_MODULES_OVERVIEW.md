# New Modules for Enhanced Biological Realism

This document provides an overview of the five new modules implemented to enhance the biological realism of the spatial transcriptomics simulation framework. Each module addresses a specific aspect of biological complexity observed in real spatial transcriptomics data.

## 1. Ligand-Receptor Interactions (06l_ligand_receptor_interactions.R)

Models cell-to-cell communication through signaling networks, affecting gene expression patterns based on the spatial arrangement of cells.

```r
# Usage example:
lr_results <- generate_lr_interactions(
  cell_df = spatial_data,
  expr_matrix = base_expression,
  dist_mat = distance_matrix,
  lr_params = list(
    use_lr_interactions = TRUE,
    n_lr_interactions = 20,
    lr_distance_decay = "exponential",
    lr_effect_strength = 0.8
  )
)
```

## 2. Temporal Dynamics (06m_temporal_dynamics.R)

Simulates developmental trajectories and temporal variation in gene expression, including RNA velocity for dynamic processes.

```r
# Usage example:
temporal_results <- generate_temporal_dynamics(
  cell_df = spatial_data,
  expr_matrix = base_expression,
  temporal_params = list(
    use_temporal_dynamics = TRUE,
    pseudotime_mode = "gradient",
    temporal_gene_fraction = 0.7,
    pattern_distribution = c(monotonic = 0.6, transient = 0.4),
    include_velocity = TRUE
  )
)
```

## 3. Alternative Splicing (06n_alternative_splicing.R)

Models spatial regulation of alternative splicing, creating variant-specific expression patterns with biological coherence.

```r
# Usage example:
splicing_results <- generate_alternative_splicing(
  cell_df = spatial_data,
  expr_matrix = base_expression,
  splicing_params = list(
    use_alternative_splicing = TRUE,
    splicing_fraction = 0.3,
    n_splicing_variants = 2,
    splicing_spatial_pattern = "gradient"
  )
)
```

## 4. Anisotropic Patterns (06o_anisotropic_patterns.R)

Creates directional gene expression along biologically relevant structures like vessels, nerves, and tissue boundaries.

```r
# Usage example:
anisotropic_results <- generate_anisotropic_patterns(
  cell_df = spatial_data,
  expr_matrix = base_expression,
  anisotropic_params = list(
    use_anisotropic_patterns = TRUE,
    n_structures = 3,
    structure_type = "branched",
    anisotropic_pattern = "gradient"
  )
)
```

## 5. 3D Microenvironment (06p_3d_microenvironment.R)

Simulates the effects of 3D tissue architecture on 2D spatial data, including cell overlap and layer-specific expression.

```r
# Usage example:
microenv_results <- generate_3d_microenvironment(
  cell_df = spatial_data,
  expr_matrix = base_expression,
  microenvironment_params = list(
    use_3d_microenvironment = TRUE,
    depth_pattern = "terrain",
    depth_range = c(0, 30),
    overlap_intensity = 0.4
  )
)
```

## Integration into Main Pipeline

These modules can be integrated into the main simulation pipeline by adding them to the appropriate stage in the expression generation process:

```r
# Example integration in simulation pipeline
simulation_result <- simulate_spatial_transcriptomics(
  input_image = "path/to/image.png",
  simulation_params = list(
    # Basic parameters...
    
    # New module parameters
    use_lr_interactions = TRUE,
    use_temporal_dynamics = TRUE,
    use_alternative_splicing = TRUE,
    use_anisotropic_patterns = TRUE,
    use_3d_microenvironment = TRUE,
    
    # Advanced parameters for each module...
  )
)
```

## Scientific Impact

These modules significantly enhance the biological realism of simulated spatial transcriptomics data by incorporating:

1. **Cell-Cell Communication**: Models how cells influence gene expression in neighboring cells
2. **Developmental Time**: Captures dynamic aspects of gene expression in spatial context
3. **RNA Processing Complexity**: Models post-transcriptional regulation in spatial domains
4. **Structural Organization**: Reflects how tissue architecture influences gene expression
5. **3D Context**: Accounts for the projection of 3D tissues onto 2D data

Together, these enhancements will allow for more realistic benchmarking of spatial transcriptomics analysis methods and better understanding of the biological factors influencing spatial gene expression patterns.