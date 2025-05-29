#!/usr/bin/env Rscript

# Simple test for improved gene expression distribution
# Focus on validating the distribution parameters only

suppressPackageStartupMessages({
  library(ggplot2)
  library(Matrix)
})

# Load required functions only
source("R/functions/06b_expression_baseline.R")
source("R/functions/06j_expression_generation.R")

# Test parameters
n_genes <- 20000
n_cells <- 100
k_cell_types <- 5

cat("\n=== Testing Improved Gene Expression Distribution ===\n")
cat(sprintf("Genes: %d, Cells: %d, Cell types: %d\n\n", n_genes, n_cells, k_cell_types))

# 1. Test baseline expression generation
cat("1. Testing baseline expression generation...\n")
baseline_expr <- generate_baseline_expression(
  n_genes = n_genes,
  k_cell_types = k_cell_types,
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  )
)

# Analyze mu distribution
all_mu <- unlist(baseline_expr)
mu_cell1 <- baseline_expr[[1]]  # First cell type

cat(sprintf("   Total mu values: %d\n", length(all_mu)))
cat(sprintf("   Mu range: [%.2f, %.2f]\n", min(all_mu), max(all_mu)))
cat(sprintf("   Mu quantiles: 25%%=%.2f, 50%%=%.2f, 75%%=%.2f\n", 
            quantile(all_mu, 0.25), quantile(all_mu, 0.5), quantile(all_mu, 0.75)))

# Count genes by expression category
mu_zero <- sum(mu_cell1 <= -15)
mu_low <- sum(mu_cell1 > -15 & mu_cell1 <= -1)
mu_med <- sum(mu_cell1 > -1 & mu_cell1 <= 1.5)
mu_high <- sum(mu_cell1 > 1.5)

cat("\n   Gene distribution by expression level (cell type 1):\n")
cat(sprintf("   - Non-expressed (mu <= -15): %d (%.1f%%)\n", mu_zero, 100*mu_zero/n_genes))
cat(sprintf("   - Low (mu > -15 & <= -1): %d (%.1f%%)\n", mu_low, 100*mu_low/n_genes))
cat(sprintf("   - Medium (mu > -1 & <= 1.5): %d (%.1f%%)\n", mu_med, 100*mu_med/n_genes))
cat(sprintf("   - High (mu > 1.5): %d (%.1f%%)\n", mu_high, 100*mu_high/n_genes))

# 2. Simulate expression for a few cells to estimate UMI
cat("\n2. Simulating expression for UMI estimation...\n")

# Simulate for each cell type
library_size <- 8000
all_umi <- numeric()

for (ct in 1:k_cell_types) {
  # Get mu values for this cell type
  mu_vals <- baseline_expr[[ct]]
  
  # Count expressed genes
  n_expressed <- sum(mu_vals > -15)
  
  # Generate expression using the same logic as generate_expression_matrix
  # Lambda calculation (corrected version)
  lambda <- exp(mu_vals) * library_size / n_expressed
  
  # Simulate for 20 cells of this type
  for (i in 1:20) {
    # Apply dispersion (negative binomial)
    dispersion <- runif(n_genes, min = 5, max = 10)
    counts <- rnbinom(n_genes, size = dispersion, mu = lambda)
    
    # Apply dropout
    dropout_prob <- runif(n_genes, min = 0.4, max = 0.6)
    dropout_mask <- runif(n_genes) > dropout_prob
    counts[dropout_mask] <- 0
    
    # Record total UMI
    all_umi <- c(all_umi, sum(counts))
  }
}

cat(sprintf("\n   Estimated UMI statistics (100 simulated cells):\n"))
cat(sprintf("   - Mean: %.0f\n", mean(all_umi)))
cat(sprintf("   - Median: %.0f\n", median(all_umi)))
cat(sprintf("   - SD: %.0f\n", sd(all_umi)))
cat(sprintf("   - Range: [%.0f, %.0f]\n", min(all_umi), max(all_umi)))

# Check targets
target_median <- 5000
target_mean <- 8000
cat("\n   Target validation:\n")
cat(sprintf("   - Target median: %d (Achieved: %.0f - %s)\n", 
            target_median, median(all_umi), 
            ifelse(median(all_umi) >= target_median * 0.8, "GOOD", "LOW")))
cat(sprintf("   - Target mean: %d (Achieved: %.0f - %s)\n", 
            target_mean, mean(all_umi),
            ifelse(abs(mean(all_umi) - target_mean) < 2000, "GOOD", "OFF")))

# 3. Visualize distribution
cat("\n3. Creating visualization...\n")

# Plot mu distribution
p1 <- ggplot(data.frame(mu = mu_cell1), aes(x = mu)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(-15, -1, 1.5), linetype = "dashed", color = "red") +
  labs(title = "Gene Expression Distribution (mu values)",
       subtitle = "50% non-expressed, 35% low, 12% medium, 3% high",
       x = "mu (log-scale expression)", y = "Count") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# Save plot
ggsave("plots/improved_distribution_mu_simple.png", p1, width = 8, height = 6, dpi = 300, bg = "white")

cat("\n=== Test Complete ===\n")

# Improvement analysis
previous_mean <- 3821
previous_median <- 393
cat(sprintf("\nImprovement from previous simulation:\n"))
cat(sprintf("- Mean: %.0f -> %.0f (%.1f%% change)\n", 
            previous_mean, mean(all_umi), 100*(mean(all_umi) - previous_mean)/previous_mean))
cat(sprintf("- Median: %.0f -> %.0f (%.1f%% change)\n", 
            previous_median, median(all_umi), 100*(median(all_umi) - previous_median)/previous_median))

if (median(all_umi) < target_median * 0.8) {
  cat("\nRecommendations:\n")
  cat("- Consider further reducing non-expressed genes to 40-45%\n")
  cat("- Or increase mu values by another 0.5\n")
} else {
  cat("\nDistribution looks good! Ready for full simulation.\n")
}