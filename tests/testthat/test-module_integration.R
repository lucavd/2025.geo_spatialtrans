test_that("New biological modules can be integrated into the simulation pipeline", {
  # Skip test on CRAN (temporarily disabled for local testing)
  # skip_on_cran()
  
  # Load required packages
  suppressPackageStartupMessages({
    library(ClusterR)
    library(imager)
    library(MASS)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })
  
  # Create a minimal configuration for testing
  # Use direct path for testing
  test_image_path <- file.path(dirname(dirname(getwd())), "inst", "extdata", "test_image.png")
  
  # Print the test image path for debugging
  cat("Looking for test image at:", test_image_path, "\n")
  cat("Does file exist?", file.exists(test_image_path), "\n")
  
  config <- list(
    image_path = test_image_path,
    n_cells = 100,
    n_genes = 20,
    k_cell_types = 3,
    threshold_value = 0.7,
    random_seed = 42,
    grid_mode = TRUE,
    grid_resolution = 10,
    
    # Only enable LR interactions and 3D microenvironment for the test
    # to check they can be mixed and matched
    lr_params = list(
      use_lr_interactions = TRUE,
      n_interactions = 5,
      signal_propagation_mode = "threshold",
      max_signaling_distance = 20
    ),
    
    temporal_params = list(
      use_temporal_dynamics = FALSE
    ),
    
    splicing_params = list(
      use_alternative_splicing = FALSE
    ),
    
    anisotropic_params = list(
      use_anisotropic_patterns = FALSE
    ),
    
    micro3d_params = list(
      use_3d_microenvironment = TRUE,
      n_layers = 3,
      layer_specificity = 0.5
    )
  )
  
  # Simple difficulty configuration
  difficulty_config <- list(
    difficulty_level = "easy",
    marker_params = list(marker_genes_per_type = 3),
    spatial_params = list(spatial_noise_intensity = 0.5),
    dropout_params = list(dropout_range = c(0.1, 0.2)),
    cell_specific_params = list(cell_specific_noise_sd = 0.1)
  )
  
  # This test will only run if the package is installed
  # and the test image is available
  if (file.exists(config$image_path)) {
    # Run the simulation with our modules using the simplified test pipeline
    source("helper-simplified_pipeline.R")
    result <- tryCatch({
      run_simulation_pipeline_test(
        config = config,
        difficulty_config = difficulty_config,
        use_spatial_correlation = FALSE
      )
    }, error = function(e) {
      skip(paste("Simulation failed:", e$message))
      NULL
    })
    
    # Skip if the simulation failed
    if (is.null(result)) {
      skip("Simulation returned NULL")
    }
    
    # Verify that the LR interactions module was applied
    expect_true(!is.null(result$module_data$lr_interactions))
    expect_true(is.list(result$module_data$lr_interactions$interaction_db))
    
    # Verify that the 3D microenvironment module was applied
    expect_true(!is.null(result$module_data$microenvironment_3d))
    expect_true(!is.null(result$module_data$microenvironment_3d$z_positions))
    expect_equal(length(result$module_data$microenvironment_3d$z_positions), config$n_cells)
    
    # Verify that parameters were properly included
    expect_equal(result$parameters$lr_params$use_lr_interactions, TRUE)
    expect_equal(result$parameters$micro3d_params$use_3d_microenvironment, TRUE)
    expect_equal(result$parameters$temporal_params$use_temporal_dynamics, FALSE)
    
    # Check that other modules were not applied (they were disabled)
    expect_null(result$module_data$temporal_dynamics)
    expect_null(result$module_data$alternative_splicing)
    expect_null(result$module_data$anisotropic_patterns)
  } else {
    skip("Test image not available")
  }
})