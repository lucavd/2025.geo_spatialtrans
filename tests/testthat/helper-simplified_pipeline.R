# Funzione di pipeline semplificata per i test di integrazione
run_simulation_pipeline_test <- function(config, difficulty_config, use_spatial_correlation = FALSE) {
  # Imposta il seed
  set.seed(config$random_seed)
  
  # 1. Crea dati immagine simulati
  img_df_thresh <- data.frame(
    x = sample(1:50, 2500, replace = TRUE),
    y = sample(1:50, 2500, replace = TRUE),
    intensity = runif(2500, 0, 1)
  )
  
  # Assegna cluster casuali
  img_df_thresh$cluster <- sample(1:config$k_cell_types, nrow(img_df_thresh), replace = TRUE)
  
  # 2. Simula il campionamento delle celle
  cell_df <- data.frame(
    x = sample(1:50, config$n_cells, replace = TRUE),
    y = sample(1:50, config$n_cells, replace = TRUE),
    intensity = runif(config$n_cells, 0, 1)
  )
  
  # Assegna cluster casuali
  cell_df$intensity_cluster <- factor(sample(1:config$k_cell_types, nrow(cell_df), replace = TRUE))
  levels(cell_df$intensity_cluster) <- paste0("cells_", letters[1:config$k_cell_types])
  
  # 3. Genera matrice di espressione simulata
  expr_matrix <- matrix(
    runif(config$n_cells * config$n_genes, 0, 10),
    nrow = config$n_cells,
    ncol = config$n_genes
  )
  colnames(expr_matrix) <- paste0("gene_", 1:config$n_genes)
  
  # 4. Crea struttura di risultati simulati dell'espressione
  expression_results <- list(
    expression = expr_matrix,
    library_size = runif(config$n_cells, 5000, 20000),
    dispersion_param = runif(config$n_genes, 0.1, 2),
    gene_modules = list(
      module_assignments = sample(1:5, config$n_genes, replace = TRUE),
      module_correlation = matrix(runif(25, -1, 1), 5, 5)
    )
  )
  
  # 5. Inizializza metadati dei moduli
  module_metadata <- list()
  
  # 6. Calcola la matrice di distanza
  dist_mat <- as.matrix(dist(cell_df[, c("x", "y")]))
  
  # 7. Applica i moduli biologici
  
  # 7.1. Ligand-Receptor Interactions
  if (config$lr_params$use_lr_interactions) {
    cat("Applicazione di interazioni ligando-recettore...\n")
    
    # Simula le interazioni L-R
    interaction_db <- data.frame(
      ligand = sample(paste0("gene_", 1:config$n_genes), 5),
      receptor = sample(paste0("gene_", 1:config$n_genes), 5),
      strength = runif(5, 0.5, 2)
    )
    
    # Simula gli effetti di segnalazione
    signaling_effects <- matrix(
      runif(config$n_cells * 5, 0.5, 2),
      nrow = config$n_cells,
      ncol = 5
    )
    
    # Applica effetti all'espressione (effetto moltiplicativo semplice)
    for (i in 1:nrow(interaction_db)) {
      ligand_idx <- which(colnames(expression_results$expression) == interaction_db$ligand[i])
      receptor_idx <- which(colnames(expression_results$expression) == interaction_db$receptor[i])
      
      if (length(ligand_idx) > 0 && length(receptor_idx) > 0) {
        # Aumenta l'espressione in base all'interazione
        expression_results$expression[, receptor_idx] <- 
          expression_results$expression[, receptor_idx] * signaling_effects[, i]
      }
    }
    
    # Salva metadati
    module_metadata$lr_interactions <- list(
      interaction_db = interaction_db,
      signaling_effects = signaling_effects
    )
  }
  
  # 7.2. Temporal Dynamics
  if (config$temporal_params$use_temporal_dynamics) {
    cat("Applicazione di dinamiche temporali...\n")
    
    # Simula il pseudotime
    pseudotime <- runif(config$n_cells, 0, 1)
    
    # Simula RNA velocity
    velocity <- list(
      vx = runif(config$n_cells, -1, 1),
      vy = runif(config$n_cells, -1, 1)
    )
    
    # Simula RNA unspliced
    unspliced <- matrix(
      runif(config$n_cells * config$n_genes, 0, 5),
      nrow = config$n_cells,
      ncol = config$n_genes
    )
    colnames(unspliced) <- paste0("gene_", 1:config$n_genes)
    
    # Salva metadati
    module_metadata$temporal_dynamics <- list(
      pseudotime = pseudotime,
      velocity = velocity,
      unspliced = unspliced
    )
  }
  
  # 7.3. Alternative Splicing
  if (config$splicing_params$use_alternative_splicing) {
    cat("Applicazione di splicing alternativo...\n")
    
    # Selezione geni con varianti
    n_genes_with_variants <- ceiling(config$n_genes * 0.3)
    genes_with_variants <- paste0("gene_", sample(1:config$n_genes, n_genes_with_variants))
    
    # Crea matrici di varianti
    variant_matrices <- list()
    for (gene in genes_with_variants) {
      # Crea 2 varianti con proporzioni complementari
      var1 <- runif(config$n_cells, 0, 1)
      var2 <- 1 - var1
      variant_matrices[[gene]] <- cbind(var1, var2)
    }
    
    # Salva metadati
    module_metadata$alternative_splicing <- list(
      genes_with_variants = genes_with_variants,
      variant_matrices = variant_matrices
    )
  }
  
  # 7.4. Anisotropic Patterns
  if (config$anisotropic_params$use_anisotropic_patterns) {
    cat("Applicazione di pattern anisotropici...\n")
    
    # Simula strutture spaziali
    structure_mask <- runif(config$n_cells, 0, 1)
    structure_mask[structure_mask < 0.7] <- 0
    
    # Simula matrici di distanza
    distance_matrices <- list()
    distance_matrices[[1]] <- dist_mat
    
    # Salva metadati
    module_metadata$anisotropic_patterns <- list(
      structure_mask = structure_mask,
      distance_matrices = distance_matrices
    )
  }
  
  # 7.5. 3D Microenvironment
  if (config$micro3d_params$use_3d_microenvironment) {
    cat("Applicazione di effetti microambiente 3D...\n")
    
    # Simula posizioni Z
    z_positions <- runif(config$n_cells, 0, config$micro3d_params$n_layers)
    
    # Assegnazione strati
    layer_assignments <- ceiling(z_positions)
    
    # Salva metadati
    module_metadata$microenvironment_3d <- list(
      z_positions = z_positions,
      layer_assignments = layer_assignments
    )
  }
  
  # 8. Prepara il risultato
  result <- list(
    coordinates = cell_df[, c("x", "y")],
    intensity_cluster = cell_df$intensity_cluster,
    expression = expression_results$expression,
    threshold_used = config$threshold_value,
    library_size = expression_results$library_size,
    dispersion_param = expression_results$dispersion_param,
    gene_modules = expression_results$gene_modules,
    
    # Aggiungi i metadati dei moduli biologici aggiuntivi
    module_data = module_metadata,
    
    parameters = list(
      image_path = config$image_path,
      pixel_size_um = config$pixel_size_um,
      n_cells = config$n_cells,
      n_genes = config$n_genes,
      k_cell_types = config$k_cell_types,
      difficulty_level = difficulty_config$difficulty_level,
      grid_mode = config$grid_mode,
      grid_resolution = config$grid_resolution,
      use_fixed_grid = config$use_fixed_grid,
      fixed_grid_width_mm = config$fixed_grid_width_mm,
      fixed_grid_height_mm = config$fixed_grid_height_mm,
      use_spatial_correlation = use_spatial_correlation,
      correlation_method = "grf",
      marker_params = difficulty_config$marker_params,
      spatial_params = difficulty_config$spatial_params,
      dropout_params = difficulty_config$dropout_params,
      library_size_params = list(
        mean_library_size = 10000,
        library_size_cv = 0.3,
        spatial_effect_on_library = 0.5,
        cell_type_effect = TRUE
      ),
      hybrid_params = list(
        use_hybrid_cells = TRUE,
        max_hybrid_pairs = 1000,
        hybrid_intensity_range = c(0.2, 0.5)
      ),
      cell_specific_params = difficulty_config$cell_specific_params,
      
      # Parametri dei nuovi moduli biologici
      lr_params = config$lr_params,
      temporal_params = config$temporal_params,
      splicing_params = config$splicing_params,
      anisotropic_params = config$anisotropic_params,
      micro3d_params = config$micro3d_params
    )
  )
  
  # Registra tempo di esecuzione fittizio
  result$execution_time <- 0.5
  
  return(result)
}