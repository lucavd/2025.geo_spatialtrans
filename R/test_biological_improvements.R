#!/usr/bin/env Rscript
# Test rapido per verificare i miglioramenti biologici

suppressPackageStartupMessages({
  library(ggplot2)
  library(Matrix)
  library(dplyr)
  library(sp)
  library(gstat)
})

cat("=== TEST RAPIDO MIGLIORAMENTI BIOLOGICI ===\n\n")

# Carica tutte le funzioni
cat("Caricamento funzioni...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (file in sort(files)) { source(file) }

# Parametri test ridotti
n_genes <- 20000  # Deve essere >10000 per attivare la nuova distribuzione
n_cells <- 100
k_cell_types <- 5

cat(sprintf("\nTest con %d geni e %d celle\n", n_genes, n_cells))

# 1. Test distribuzione genica migliorata
cat("\n1. Test distribuzione genica migliorata...\n")
baseline_expr <- generate_baseline_expression(
  n_genes = n_genes,
  k_cell_types = k_cell_types,
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  )
)

# Analizza distribuzione mu
all_mu <- unlist(baseline_expr)
mu_zero <- sum(all_mu <= -15)
mu_low <- sum(all_mu > -15 & all_mu <= -1)
mu_med <- sum(all_mu > -1 & all_mu <= 1.5)
mu_high <- sum(all_mu > 1.5)

cat("   Distribuzione geni:\n")
cat(sprintf("   - Non-espressi: %d (%.1f%%) [target: 40%%]\n", mu_zero, 100*mu_zero/length(all_mu)))
cat(sprintf("   - Low: %d (%.1f%%) [target: 35%%]\n", mu_low, 100*mu_low/length(all_mu)))
cat(sprintf("   - Medium: %d (%.1f%%) [target: 20%%]\n", mu_med, 100*mu_med/length(all_mu)))
cat(sprintf("   - High: %d (%.1f%%) [target: 5%%]\n", mu_high, 100*mu_high/length(all_mu)))

# 2. Test generazione espressione con limiti biologici
cat("\n2. Test generazione espressione con limiti biologici...\n")

# Crea dataframe celle di test
cell_df <- data.frame(
  original_index = 1:n_cells,
  x = runif(n_cells, 0, 100),
  y = runif(n_cells, 0, 100),
  cluster = sample(1:k_cell_types, n_cells, replace = TRUE),
  intensity_cluster = factor(sample(1:k_cell_types, n_cells, replace = TRUE))
)

# Genera espressione
expr_results <- generate_expression_profiles(
  cell_df = cell_df,
  n_genes = n_genes,
  k_cell_types = k_cell_types,
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5
  ),
  spatial_params = list(
    spatial_range = 0.2,  # Corretto nome parametro
    sill_spatial = 0.5,
    nugget_spatial = 0.1
  ),
  dropout_params = list(
    dropout_rate_min = 0.4,
    dropout_rate_max = 0.6,
    dropout_spatial_effect = 0.1,
    cell_type_effect = TRUE,
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.5,
    dropout_curve_steepness = 5
  ),
  library_size_params = list(
    mean_library_size = 8000,
    library_size_cv = 0.3,
    spatial_effect_on_library = 0.1,
    cell_type_effect = TRUE
  ),
  use_spatial_correlation = FALSE,
  random_seed = 42
)

# Analizza risultati
expr_matrix <- expr_results$expression_matrix
umi_counts <- Matrix::colSums(expr_matrix)

cat("\n   Statistiche UMI:\n")
cat(sprintf("   - Media: %.0f [target: 6000-10000]\n", mean(umi_counts)))
cat(sprintf("   - Mediana: %.0f [target: ~5000]\n", median(umi_counts)))
cat(sprintf("   - Max: %.0f [dovrebbe essere <50000]\n", max(umi_counts)))
cat(sprintf("   - CV: %.2f\n", sd(umi_counts)/mean(umi_counts)))

# Verifica valori massimi per gene
max_per_gene <- apply(expr_matrix, 1, max)
genes_over_5k <- sum(max_per_gene > 5000)
genes_over_10k <- sum(max_per_gene > 10000)

cat("\n   Controllo limiti per gene:\n")
cat(sprintf("   - Geni con max >5000 UMI: %d (%.1f%%)\n", genes_over_5k, 100*genes_over_5k/n_genes))
cat(sprintf("   - Geni con max >10000 UMI: %d (%.1f%%)\n", genes_over_10k, 100*genes_over_10k/n_genes))

# Statistiche espressione non-zero
non_zero_vals <- expr_matrix@x
cat("\n   Espressione (valori non-zero):\n")
cat(sprintf("   - Media: %.2f [target: 2-10]\n", mean(non_zero_vals)))
cat(sprintf("   - Mediana: %.2f [target: 2-5]\n", median(non_zero_vals)))
cat(sprintf("   - Max: %.0f\n", max(non_zero_vals)))

# 3. Verifica plausibilità biologica
cat("\n3. Verifica plausibilità biologica...\n")

# Conta celle in range realistici
realistic_cells <- sum(umi_counts >= 1000 & umi_counts <= 15000)
cat(sprintf("   - Celle con UMI realistici [1k-15k]: %d (%.1f%%)\n", 
            realistic_cells, 100*realistic_cells/length(umi_counts)))

# Verifica se i controlli biologici funzionano
if (max(non_zero_vals) > 10000) {
  cat("   ⚠️  ATTENZIONE: Ancora valori >10k nonostante i limiti!\n")
} else {
  cat("   ✓ Tutti i valori sotto 10k UMI per gene\n")
}

if (median(non_zero_vals) < 1) {
  cat("   ⚠️  ATTENZIONE: Mediana espressione ancora troppo bassa!\n")
} else {
  cat("   ✓ Mediana espressione in range accettabile\n")
}

# 4. Plot diagnostici
cat("\n4. Generazione plot diagnostici...\n")

# Plot distribuzione UMI
p1 <- ggplot(data.frame(umi = umi_counts), aes(x = umi)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(5000, mean(umi_counts)), 
             linetype = c("dashed", "solid"), 
             color = c("red", "blue")) +
  labs(title = "Distribuzione UMI con miglioramenti biologici",
       subtitle = sprintf("Media: %.0f (blu), Target mediana: 5000 (rosso)", mean(umi_counts)),
       x = "UMI per cella", y = "Conteggio") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# Plot valori espressione
p2 <- ggplot(data.frame(expr = non_zero_vals[1:min(10000, length(non_zero_vals))]), 
             aes(x = expr)) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
  scale_x_log10() +
  geom_vline(xintercept = median(non_zero_vals), color = "red", linetype = "dashed") +
  labs(title = "Distribuzione valori espressione (non-zero)",
       subtitle = sprintf("Mediana: %.2f", median(non_zero_vals)),
       x = "UMI (scala log)", y = "Frequenza") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("plots/test_biological_umi.png", p1, width = 8, height = 6, dpi = 150, bg = "white")
ggsave("plots/test_biological_expr.png", p2, width = 8, height = 6, dpi = 150, bg = "white")

# 5. Conclusioni
cat("\n=== CONCLUSIONI TEST ===\n")

test_passed <- TRUE

if (abs(mean(umi_counts) - 8000) / 8000 > 0.5) {
  cat("❌ Media UMI fuori range target\n")
  test_passed <- FALSE
} else {
  cat("✓ Media UMI nel range target\n")
}

if (median(umi_counts) < 2000) {
  cat("❌ Mediana UMI troppo bassa\n")
  test_passed <- FALSE
} else {
  cat("✓ Mediana UMI accettabile\n")
}

if (genes_over_10k > 0) {
  cat("❌ Presenza di valori >10k UMI per gene\n")
  test_passed <- FALSE
} else {
  cat("✓ Nessun valore >10k UMI per gene\n")
}

if (median(non_zero_vals) < 1) {
  cat("❌ Mediana espressione troppo bassa\n")
  test_passed <- FALSE
} else {
  cat("✓ Mediana espressione accettabile\n")
}

if (test_passed) {
  cat("\n✅ TEST SUPERATO! Procedere con simulazione full-size.\n")
} else {
  cat("\n⚠️  TEST FALLITO. Verificare i parametri prima di procedere.\n")
}

cat("\nPlot salvati in:\n")
cat("- plots/test_biological_umi.png\n")
cat("- plots/test_biological_expr.png\n")