#!/usr/bin/env Rscript
# Quick smoke test for the new MRF spatial engine
source("R_simple/07_mrf_generation.R")

df <- simulate_mrf(grid_size = 64, k_cell_types = 5, beta = 0.7, n_iter = 100, seed = 42, fast_mode = TRUE)

cat("\n[mrf_test] Generated", nrow(df), "cells\n")

# Simple Moran's I implementation (rook adjacency on lattice)
library(stats)

grid_size <- 64
mat <- matrix(df$cell_type, nrow = grid_size, byrow = TRUE)

# Compute mean
m <- mean(mat)
# Neighbor pairs count
w_total <- 0
num <- 0
for (i in 1:grid_size) {
  for (j in 1:grid_size) {
    v <- mat[i, j] - m
    # right neighbor
    if (j < grid_size) { num <- num + v * (mat[i, j+1] - m); w_total <- w_total + 1 }
    # down neighbor
    if (i < grid_size) { num <- num + v * (mat[i+1, j] - m); w_total <- w_total + 1 }
  }
}

den <- sum((mat - m)^2)
I <- (grid_size^2 / w_total) * (num / den)
cat("[mrf_test] Moran's I:", round(I, 3), "(>0 indicates positive spatial autocorrelation)\n")
