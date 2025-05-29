#!/usr/bin/env Rscript
# Test rapido della nuova distribuzione genica per simulazioni full-size

cat("=== TEST NUOVA DISTRIBUZIONE GENICA ===\n\n")

# Carica le funzioni
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# Test con parametri full-size
n_genes <- 20000
k_cell_types <- 10

# Genera baseline expression
baseline <- generate_baseline_expression(
  n_genes = n_genes,
  k_cell_types = k_cell_types,
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  ),
  random_seed = 42
)

# Debug: verifica struttura del risultato
cat("Nomi in baseline:", names(baseline), "\n")
cat("Classe baseline:", class(baseline), "\n")

# La funzione restituisce direttamente la lista, non un oggetto con mean_expression_list
if (is.list(baseline) && is.null(names(baseline))) {
  # È una lista senza nomi, quindi è la mean_expression_list
  mu_values <- baseline[[1]]
} else {
  mu_values <- baseline$mean_expression_list[[1]]
}

cat("\nClasse di mu_values:", class(mu_values), "\n")
cat("Lunghezza:", length(mu_values), "\n")

exp_values <- exp(mu_values)  # Converti da log-space

cat("DISTRIBUZIONE MU (log-space) PER 20,000 GENI:\n")
cat("- Geni con mu <= -10 (praticamente zero):", sum(mu_values <= -10), 
    "(", round(sum(mu_values <= -10)/n_genes*100, 1), "%)\n")
cat("- Geni con -10 < mu <= -4 (low):", sum(mu_values > -10 & mu_values <= -4),
    "(", round(sum(mu_values > -10 & mu_values <= -4)/n_genes*100, 1), "%)\n")
cat("- Geni con -4 < mu <= -1 (medium):", sum(mu_values > -4 & mu_values <= -1),
    "(", round(sum(mu_values > -4 & mu_values <= -1)/n_genes*100, 1), "%)\n")
cat("- Geni con mu > -1 (high):", sum(mu_values > -1),
    "(", round(sum(mu_values > -1)/n_genes*100, 1), "%)\n")

cat("\nVALORI ATTESI DI ESPRESSIONE:\n")
cat("- Media dei valori mu:", round(mean(mu_values), 2), "\n")
cat("- Mediana dei valori mu:", round(median(mu_values), 2), "\n")
cat("- Media espressione (exp(mu)):", round(mean(exp_values), 3), "\n")
cat("- Mediana espressione (exp(mu)):", round(median(exp_values), 6), "\n")

# Simula una piccola matrice per verificare UMI
cat("\nSIMULAZIONE TEST (100 celle):\n")
n_test_cells <- 100
test_expression <- matrix(0, nrow = n_genes, ncol = n_test_cells)

# Simula espressione per ogni gene
library_size <- 8000

# Conta solo i geni effettivamente espressi (non-zero)
expressed_genes <- which(mu_values > -10)
n_expressed <- length(expressed_genes)
cat("Geni effettivamente espressi:", n_expressed, "\n")

for (g in 1:n_genes) {
  if (mu_values[g] > -10) {
    # Dividi la library size solo tra i geni espressi
    lambda <- exp_values[g] * library_size / n_expressed
    test_expression[g, ] <- rpois(n_test_cells, lambda)
  }
}

# Calcola UMI
umi_per_cell <- colSums(test_expression)
genes_detected <- colSums(test_expression > 0)

cat("- UMI medio per cella:", round(mean(umi_per_cell)), "\n")
cat("- UMI mediano per cella:", round(median(umi_per_cell)), "\n")
cat("- Geni rilevati medio per cella:", round(mean(genes_detected)), "\n")
cat("- Sparsità:", round(mean(test_expression == 0) * 100, 1), "%\n")

# Verifica quanti geni sono completamente non rilevati
genes_never_detected <- sum(rowSums(test_expression > 0) == 0)
cat("- Geni mai rilevati in 100 celle:", genes_never_detected,
    "(", round(genes_never_detected/n_genes*100, 1), "%)\n")

cat("\nPREVISIONE PER SIMULAZIONE COMPLETA:\n")
cat("Con questa distribuzione, ci aspettiamo:\n")
cat("- ~60% dei geni non rilevati o rarissimi\n")
cat("- UMI per cella più realistici (~6-10k)\n")
cat("- Alta sparsità ma biologicamente plausibile\n")