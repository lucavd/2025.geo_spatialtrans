#!/usr/bin/env Rscript

# Analyze the improved simulation results in detail

library(Matrix)
library(ggplot2)

cat("=== ANALISI DETTAGLIATA RISULTATI MIGLIORATI ===\n\n")

# Load results
sim_results <- readRDS("results/visiumHD_biological.rds")

# Calculate UMI counts
umi_counts <- Matrix::colSums(sim_results$expression)

# Remove extreme outliers for better analysis
# These might be technical artifacts or edge cases
umi_filtered <- umi_counts[umi_counts > 0 & umi_counts < 100000]

cat("Statistiche UMI complete:\n")
cat(sprintf("- N. celle totali: %d\n", length(umi_counts)))
cat(sprintf("- N. celle con UMI = 0: %d (%.1f%%)\n", 
            sum(umi_counts == 0), 100*sum(umi_counts == 0)/length(umi_counts)))
cat(sprintf("- N. celle con UMI > 100k: %d (%.1f%%)\n", 
            sum(umi_counts > 100000), 100*sum(umi_counts > 100000)/length(umi_counts)))

cat("\nStatistiche UMI filtrate (0 < UMI < 100k):\n")
cat(sprintf("- N. celle analizzate: %d (%.1f%%)\n", 
            length(umi_filtered), 100*length(umi_filtered)/length(umi_counts)))
cat(sprintf("- Media: %.0f\n", mean(umi_filtered)))
cat(sprintf("- Mediana: %.0f\n", median(umi_filtered)))
cat(sprintf("- SD: %.0f\n", sd(umi_filtered)))
cat(sprintf("- Quantili: 25%%=%.0f, 75%%=%.0f, 95%%=%.0f\n", 
            quantile(umi_filtered, 0.25), quantile(umi_filtered, 0.75), quantile(umi_filtered, 0.95)))

# Check gene detection
genes_detected_per_cell <- Matrix::colSums(sim_results$expression > 0)
genes_never_detected <- sum(Matrix::rowSums(sim_results$expression > 0) == 0)

cat("\nStatistiche rilevamento geni:\n")
cat(sprintf("- Geni mai rilevati: %d (%.1f%%)\n", 
            genes_never_detected, 100*genes_never_detected/nrow(sim_results$expression)))
cat(sprintf("- Media geni/cella: %.0f\n", mean(genes_detected_per_cell)))
cat(sprintf("- Mediana geni/cella: %.0f\n", median(genes_detected_per_cell)))

# Analyze expression distribution
expression_values <- sim_results$expression@x  # Non-zero values only
cat("\nDistribuzione valori di espressione (non-zero):\n")
cat(sprintf("- Min: %.0f\n", min(expression_values)))
cat(sprintf("- Max: %.0f\n", max(expression_values)))
cat(sprintf("- Media: %.2f\n", mean(expression_values)))
cat(sprintf("- Mediana: %.2f\n", median(expression_values)))

# Target comparison
cat("\n=== CONFRONTO CON TARGET ===\n")
target_median <- 5000
target_mean_range <- c(6000, 10000)
realistic_range <- c(1000, 15000)

# Use filtered data for fair comparison
realistic_cells <- sum(umi_filtered >= realistic_range[1] & umi_filtered <= realistic_range[2])
realistic_pct <- 100 * realistic_cells / length(umi_filtered)

cat(sprintf("Target mediano UMI: %d (Raggiunto: %.0f - %s)\n", 
            target_median, median(umi_filtered), 
            ifelse(median(umi_filtered) >= target_median * 0.8, "✓", "✗")))
cat(sprintf("Target medio UMI: [%d, %d] (Raggiunto: %.0f - %s)\n", 
            target_mean_range[1], target_mean_range[2], mean(umi_filtered),
            ifelse(mean(umi_filtered) >= target_mean_range[1] && mean(umi_filtered) <= target_mean_range[2], "✓", "✗")))
cat(sprintf("Celle con UMI realistici [%d-%d]: %.1f%%\n", 
            realistic_range[1], realistic_range[2], realistic_pct))

# Plot filtered distribution
p <- ggplot(data.frame(umi = umi_filtered), aes(x = umi)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(target_median, mean(umi_filtered)), 
             linetype = c("dashed", "solid"), color = c("red", "blue"), size = 1) +
  geom_vline(xintercept = realistic_range, linetype = "dotted", color = "green") +
  scale_x_log10() +
  labs(title = "Distribuzione UMI (filtrata, scala log)",
       subtitle = sprintf("Mediana: %.0f (rosso), Media: %.0f (blu), Range realistico: %d-%d (verde)",
                         median(umi_filtered), mean(umi_filtered), realistic_range[1], realistic_range[2]),
       x = "UMI per cella (log scale)", y = "Numero di celle") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("plots/umi_distribution_filtered.png", p, width = 10, height = 6, dpi = 300, bg = "white")

# Diagnosis
cat("\n=== DIAGNOSI ===\n")
if (sum(umi_counts == 0) > length(umi_counts) * 0.1) {
  cat("⚠️  Troppe celle con 0 UMI - possibile problema con dropout eccessivo\n")
}
if (sum(umi_counts > 100000) > 10) {
  cat("⚠️  Presenza di outlier estremi - possibile instabilità numerica\n")
}
if (median(umi_filtered) < 2000) {
  cat("⚠️  Mediana ancora bassa - la distribuzione potrebbe necessitare ulteriori aggiustamenti\n")
} else {
  cat("✓ Mediana UMI in range accettabile\n")
}

# Check the gene expression distribution that was used
cat("\n=== VERIFICA DISTRIBUZIONE GENICA UTILIZZATA ===\n")
cat("Controllare i messaggi DEBUG durante la simulazione per confermare l'uso della nuova distribuzione.\n")
cat("Dovrebbe mostrare: '[DEBUG] Using improved gene distribution: 50% non-expressed, 35% low, 12% medium, 3% high'\n")