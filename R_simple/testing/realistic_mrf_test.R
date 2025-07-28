#!/usr/bin/env Rscript
# Test the enhanced MRF with realistic biological patterns
source("R_simple/07_mrf_generation.R")

cat("=== TESTING REALISTIC MRF PATTERNS ===\n")

# Test different tissue structures
structures <- c("vessel", "boundary", "gradient", "uniform")
grid_size <- 100
k_cell_types <- 6

for (structure in structures) {
  cat("\n--- Testing", structure, "structure ---\n")
  
  start_time <- Sys.time()
  df <- simulate_mrf(
    grid_size = grid_size, 
    k_cell_types = k_cell_types, 
    beta = 0.6, 
    n_iter = 150, 
    seed = 42,
    fast_mode = TRUE,
    tissue_structure = structure
  )
  end_time <- Sys.time()
  
  cat("Generated", nrow(df), "cells in", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), "seconds\n")
  
  # Compute spatial statistics
  mat <- matrix(df$cell_type, nrow = grid_size, byrow = TRUE)
  
  # Moran's I
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
  
  # Cell type diversity (Shannon entropy)
  type_counts <- table(df$cell_type)
  props <- type_counts / sum(type_counts)
  entropy <- -sum(props * log(props))
  
  cat("Moran's I:", round(I, 3), "| Entropy:", round(entropy, 3), 
      "| Types present:", length(type_counts), "/", k_cell_types, "\n")
}

cat("\n=== BIOLOGICAL INTERACTION MATRIX TEST ===\n")
# Test interaction matrix
interaction_mat <- create_biological_interactions(6, 42)
cat("Interaction matrix (6x6):\n")
print(round(interaction_mat, 2))

cat("\n✅ All realistic MRF tests completed!\n")
