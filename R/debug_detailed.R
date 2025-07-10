#!/usr/bin/env Rscript
# Debug dettagliato per capire la normalizzazione

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
})

files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# Test su dataset piccolissimo per analisi dettagliata
cat("=== ANALISI DETTAGLIATA NORMALIZZAZIONE ===\n")

# Simula 10 celle, 10 geni
set.seed(42)
N <- 10
n_genes <- 10
k_cell_types <- 2

# Mock data
cell_df <- data.frame(
  x = 1:N,
  y = 1:N,
  intensity_cluster = factor(rep(1:k_cell_types, length.out = N))
)

# Mock mean expression (geni con espressioni diverse)
mean_expression_list <- list(
  # Tipo 1: alcuni geni alti, alcuni bassi
  c(-1, -2, -3, -4, -5, -6, -7, -8, -9, -10),
  # Tipo 2: pattern diverso
  c(-2, -1, -4, -3, -6, -5, -8, -7, -10, -9)
)

# Mock library sizes (target 1000 UMI)
target_library_size <- rep(1000, N)

cat("Target library sizes:\n")
print(target_library_size)

# Simula la logica di normalizzazione attuale
all_mean_expr <- matrix(0, nrow = n_genes, ncol = k_cell_types)
for (g in seq_len(n_genes)) {
  for (k in seq_len(k_cell_types)) {
    all_mean_expr[g, k] <- mean_expression_list[[k]][g]
  }
}

cl <- as.integer(cell_df$intensity_cluster)

cat("\nMatrice espressioni medie per tipo:\n")
print(all_mean_expr)

cat("\nCluster assignments:\n")
print(cl)

# Calcola normalizzazione attuale
total_exp_mu_per_cell <- sapply(1:N, function(cell_idx) {
  sum(exp(all_mean_expr[, cl[cell_idx]]))
})

cat("\nSomma exp(mu) per cella (total_exp_mu_per_cell):\n")
print(total_exp_mu_per_cell)

# Test per un gene specifico (gene 1)
g <- 1
mu_vals <- all_mean_expr[g, cl]
cat("\nPer gene", g, ":\n")
cat("mu_vals per cella:", mu_vals, "\n")
cat("exp(mu_vals):", exp(mu_vals), "\n")

# Calcola lambda con normalizzazione corrente
lambda_current <- exp(mu_vals) * target_library_size / total_exp_mu_per_cell
cat("lambda attuale:", lambda_current, "\n")

# Confronta con normalizzazione semplice (vecchia)
lambda_old <- exp(mu_vals) * target_library_size / n_genes
cat("lambda vecchia:", lambda_old, "\n")

# Calcola somma teorica se tutti i geni contribuissero con lambda_current
total_expected_per_cell <- lambda_current * n_genes
cat("UMI totali attesi se tutti i geni contribuissero così:", total_expected_per_cell, "\n")

cat("\nRapporto target/atteso:\n")
print(target_library_size / total_expected_per_cell)