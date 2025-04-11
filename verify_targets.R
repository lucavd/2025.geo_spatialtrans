# Esempio di utilizzo con targets - quando pronto

library(targets)
library(tarchetypes)

# Definisci il file _targets.R
tar_script({
  
  # Carica librerie
  library(dplyr)
  library(ggplot2)
  
  # Carica i file delle funzioni
  lapply(list.files("R/functions", full.names = TRUE, pattern = "\\.R$"), source)
  
  # Pipeline principale
  list(
    # Configurazione di base
    tar_target(
      config,
      initialize_simulation_config(
        image_path = "images/colon.png",
        output_path = "results/simulation_out.rds",
        n_cells = 5000,
        n_genes = 50,
        k_cell_types = 4,
        grid_mode = TRUE,
        random_seed = 123
      )
    ),
    
    # Configurazione difficoltà
    tar_target(
      difficulty,
      configure_difficulty_level(
        difficulty_level = "medium",
        marker_params = list(marker_genes_per_type = 6),
        dropout_params = list(dropout_range = c(0.15, 0.4))
      )
    ),
    
    # Preparazione immagine
    tar_target(
      image_data,
      prepare_image(config$image_path, config$threshold_value)
    ),
    
    # Clustering dell'immagine
    tar_target(
      clustered_data,
      cluster_image(
        img_df_thresh = image_data$img_df_thresh,
        k_cell_types = config$k_cell_types,
        random_seed = config$random_seed
      )
    ),
    
    # Creazione griglia
    tar_target(
      cell_df,
      create_sampling_grid(
        img_df_thresh = clustered_data,
        img_array = image_data$img_array,
        img_width = image_data$width,
        img_height = image_data$height,
        grid_mode = config$grid_mode,
        n_cells = config$n_cells,
        k_cell_types = config$k_cell_types,
        grid_resolution = config$grid_resolution,
        grid_spacing = config$grid_spacing,
        use_fixed_grid = config$use_fixed_grid,
        fixed_grid_width_mm = config$fixed_grid_width_mm,
        fixed_grid_height_mm = config$fixed_grid_height_mm,
        pixel_size_um = config$pixel_size_um,
        threshold_value = config$threshold_value,
        random_seed = config$random_seed
      )
    ),
    
    # Calcolo distanze
    tar_target(
      spatial_distances,
      calculate_spatial_distances(cell_df, random_seed = config$random_seed)
    ),
    
    # Creazione profili base
    tar_target(
      baseline_profiles,
      generate_baseline_expression(
        n_genes = config$n_genes,
        k_cell_types = config$k_cell_types,
        marker_params = difficulty$marker_params,
        random_seed = config$random_seed
      )
    ),
    
    # Calcolo dispersione
    tar_target(
      dispersion_params,
      calculate_dispersion_params(
        cell_df,
        spatial_distances$mean_dist,
        dropout_params = difficulty$dropout_params,
        k_cell_types = config$k_cell_types,
        random_seed = config$random_seed
      )
    ),
    
    # Generazione dimensioni libreria
    tar_target(
      library_sizes,
      generate_library_sizes(
        cell_df,
        library_size_params = list(
          mean_library_size = 8000,
          library_size_cv = 0.3,
          spatial_effect_on_library = 0.5,
          cell_type_effect = TRUE
        ),
        spatial_params = difficulty$spatial_params,
        k_cell_types = config$k_cell_types,
        random_seed = config$random_seed
      )
    ),
    
    # Probabilità dropout
    tar_target(
      dropout_probs,
      calculate_dropout_probabilities(
        cell_df,
        spatial_distances$mean_dist,
        dropout_params = difficulty$dropout_params,
        spatial_params = difficulty$spatial_params
      )
    ),
    
    # Moduli genici
    tar_target(
      gene_modules,
      generate_gene_modules(
        config$n_genes,
        nrow(cell_df),
        cell_specific_params = difficulty$cell_specific_params,
        random_seed = config$random_seed
      )
    ),
    
    # Correlazione spaziale
    tar_target(
      spatial_corr,
      generate_spatial_correlation(
        cell_df,
        difficulty$spatial_params,
        correlation_method = "grf",
        random_seed = config$random_seed
      )
    ),
    
    # Cellule ibride
    tar_target(
      hybrid_cells,
      generate_hybrid_cells(
        cell_df,
        spatial_distances$dist_mat,
        k_cell_types = config$k_cell_types,
        hybrid_params = list(
          use_hybrid_cells = TRUE,
          max_hybrid_pairs = 500,
          hybrid_intensity_range = c(0.2, 0.5)
        ),
        random_seed = config$random_seed
      )
    ),
    
    # Matrice di espressione finale
    tar_target(
      expression_matrix,
      generate_expression_matrix(
        cell_df,
        baseline_profiles,
        config$n_genes,
        library_sizes,
        dispersion_params,
        dropout_probs,
        hybrid_cells,
        gene_modules$module_noise,
        spatial_corr,
        difficulty$spatial_params,
        difficulty$dropout_params,
        difficulty$cell_specific_params,
        use_spatial_correlation = TRUE,
        random_seed = config$random_seed
      )
    ),
    
    # Plot finale
    tar_target(
      final_plot,
      create_simulation_plots(
        cell_df,
        config,
        difficulty
      )
    ),
    
    # Risultato finale
    tar_target(
      final_result,
      list(
        coordinates = cell_df[, c("x", "y")],
        intensity_cluster = cell_df$intensity_cluster,
        expression = expression_matrix,
        threshold_used = config$threshold_value,
        library_size = library_sizes,
        dispersion_param = dispersion_params
      )
    ),
    
    # Salvataggio risultati
    tar_target(
      saved_result,
      save_simulation_results(final_result, config$output_path)
    )
  )
}, ask = FALSE)

# Per eseguire la pipeline:
# targets::tar_make()

# Per visualizzare la struttura della pipeline:
# targets::tar_visnetwork()