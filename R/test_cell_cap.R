#!/usr/bin/env Rscript
# Test rapido del cap a livello di cella

library(Matrix)
source("R/functions/06j_expression_generation.R")

cat("=== TEST CAP A LIVELLO DI CELLA ===\n\n")

# Crea una matrice di test con alcune celle estreme
set.seed(123)
n_genes <- 1000
n_cells <- 10

# Crea matrice con alcune celle con valori molto alti
test_matrix <- Matrix(0, nrow = n_genes, ncol = n_cells, sparse = TRUE)

# Cella 1: normale (10k UMI totali)
test_matrix[sample(n_genes, 200), 1] <- rpois(200, lambda = 50)

# Cella 2: alta ma sotto il limite (40k UMI)
test_matrix[sample(n_genes, 400), 2] <- rpois(400, lambda = 100)

# Cella 3: sopra il limite (100k UMI)
test_matrix[sample(n_genes, 500), 3] <- rpois(500, lambda = 200)

# Cella 4: molto sopra il limite (500k UMI)
test_matrix[sample(n_genes, 800), 4] <- rpois(800, lambda = 625)

# Cella 5: estrema (1M UMI)
test_matrix[sample(n_genes, 900), 5] <- rpois(900, lambda = 1111)

# Calcola totali prima
totals_before <- colSums(test_matrix)
cat("UMI totali PRIMA del cap:\n")
for (i in 1:5) {
  cat(sprintf("  Cella %d: %d UMI\n", i, totals_before[i]))
}

# Applica validazione biologica
cat("\nApplicazione validate_biological_plausibility...\n")
test_matrix_validated <- validate_biological_plausibility(test_matrix)

# Calcola totali dopo
totals_after <- colSums(test_matrix_validated)
cat("\nUMI totali DOPO il cap:\n")
for (i in 1:5) {
  cat(sprintf("  Cella %d: %d UMI (era %d)\n", i, totals_after[i], totals_before[i]))
}

# Verifica
cat("\n=== VERIFICA ===\n")
if (all(totals_after <= 50000)) {
  cat("✓ TEST PASSATO: Tutte le celle sono sotto 50k UMI\n")
} else {
  cat("✗ TEST FALLITO: Alcune celle ancora sopra 50k UMI\n")
  cat("  Celle problematiche:\n")
  prob_cells <- which(totals_after > 50000)
  for (cell in prob_cells) {
    cat(sprintf("    Cella %d: %d UMI\n", cell, totals_after[cell]))
  }
}