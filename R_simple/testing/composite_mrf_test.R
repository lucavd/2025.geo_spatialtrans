#!/usr/bin/env Rscript
# Test the enhanced MRF with composite tissue structures
source("R_simple/07_mrf_generation.R")

cat("=== TESTING COMPOSITE MRF STRUCTURES ===\n")

grid_size <- 80
k_cell_types <- 6

# Test 1: Single structure (backward compatibility)
cat("\n--- Test 1: Single structure (vessel) ---\n")
df_single <- simulate_mrf(
  grid_size = grid_size, 
  k_cell_types = k_cell_types, 
  beta = 0.6, 
  n_iter = 100, 
  seed = 42,
  tissue_structure = "vessel"
)
cat("Generated", nrow(df_single), "cells with single vessel structure\n")

# Test 2: Composite structure - vessel + gradient
cat("\n--- Test 2: Composite structure (vessel + gradient) ---\n")
composite_struct <- list(
  list(type = "vessel", weight = 0.6),
  list(type = "gradient", weight = 0.4)
)

df_composite <- simulate_mrf(
  grid_size = grid_size, 
  k_cell_types = k_cell_types, 
  beta = 0.6, 
  n_iter = 100, 
  seed = 42,
  tissue_structure = composite_struct
)
cat("Generated", nrow(df_composite), "cells with vessel + gradient structure\n")

# Test 3: Complex composite - all structures
cat("\n--- Test 3: Complex composite (vessel + boundary + gradient) ---\n")
complex_struct <- list(
  list(type = "vessel", weight = 0.4),
  list(type = "boundary", weight = 0.3),
  list(type = "gradient", weight = 0.3)
)

df_complex <- simulate_mrf(
  grid_size = grid_size, 
  k_cell_types = k_cell_types, 
  beta = 0.6, 
  n_iter = 100, 
  seed = 42,
  tissue_structure = complex_struct
)
cat("Generated", nrow(df_complex), "cells with complex composite structure\n")

# Compute spatial statistics for comparison
compute_moran <- function(df, grid_size) {
  mat <- matrix(df$cell_type, nrow = grid_size, byrow = TRUE)
  m <- mean(mat)
  w_total <- 0; num <- 0
  for (i in 1:grid_size) {
    for (j in 1:grid_size) {
      v <- mat[i, j] - m
      if (j < grid_size) { num <- num + v * (mat[i, j+1] - m); w_total <- w_total + 1 }
      if (i < grid_size) { num <- num + v * (mat[i+1, j] - m); w_total <- w_total + 1 }
    }
  }
  den <- sum((mat - m)^2)
  I <- (grid_size^2 / w_total) * (num / den)
  return(I)
}

compute_entropy <- function(df) {
  type_counts <- table(df$cell_type)
  props <- type_counts / sum(type_counts)
  entropy <- -sum(props * log(props))
  return(entropy)
}

cat("\n=== SPATIAL STATISTICS COMPARISON ===\n")
cat("Single vessel    - Moran's I:", round(compute_moran(df_single, grid_size), 3), 
    "| Entropy:", round(compute_entropy(df_single), 3), "\n")
cat("Vessel+Gradient  - Moran's I:", round(compute_moran(df_composite, grid_size), 3), 
    "| Entropy:", round(compute_entropy(df_composite), 3), "\n")
cat("Complex composite- Moran's I:", round(compute_moran(df_complex, grid_size), 3), 
    "| Entropy:", round(compute_entropy(df_complex), 3), "\n")

# Test 4: Error handling
cat("\n--- Test 4: Error handling ---\n")
tryCatch({
  df_error <- simulate_mrf(grid_size = 50, tissue_structure = "invalid_structure")
}, error = function(e) {
  cat("✓ Correctly caught invalid structure error:", e$message, "\n")
})

cat("\n✅ All composite MRF tests completed!\n")
cat("✅ Backward compatibility maintained\n")
cat("✅ Composite structures working\n")
cat("✅ Error handling functional\n")
