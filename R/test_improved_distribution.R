#!/usr/bin/env Rscript

# Test improved gene expression distribution
# This script validates the new 50% non-expressed distribution

suppressPackageStartupMessages({
  library(ggplot2)
  library(Matrix)
  library(dplyr)
  library(sp)
  library(gstat)
})

# Carica tutte le funzioni
cat("Loading functions...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (file in sort(files)) { source(file) }

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
cat(sprintf("   Total mu values: %d\n", length(all_mu)))
cat(sprintf("   Mu range: [%.2f, %.2f]\n", min(all_mu), max(all_mu)))
cat(sprintf("   Mu quantiles: 25%%=%.2f, 50%%=%.2f, 75%%=%.2f\n", 
            quantile(all_mu, 0.25), quantile(all_mu, 0.5), quantile(all_mu, 0.75)))

# Count genes by expression category
mu_zero <- sum(all_mu <= -15)
mu_low <- sum(all_mu > -15 & all_mu <= -1)
mu_med <- sum(all_mu > -1 & all_mu <= 1.5)
mu_high <- sum(all_mu > 1.5)

cat("\n   Gene distribution by expression level:\n")
cat(sprintf("   - Non-expressed (mu <= -15): %d (%.1f%%)\n", mu_zero, 100*mu_zero/length(all_mu)))
cat(sprintf("   - Low (mu > -15 & <= -1): %d (%.1f%%)\n", mu_low, 100*mu_low/length(all_mu)))
cat(sprintf("   - Medium (mu > -1 & <= 1.5): %d (%.1f%%)\n", mu_med, 100*mu_med/length(all_mu)))
cat(sprintf("   - High (mu > 1.5): %d (%.1f%%)\n", mu_high, 100*mu_high/length(all_mu)))

# 2. Simulate a small expression matrix to test UMI counts
cat("\n2. Testing expression generation with improved distribution...\n")

# Create simple cell dataframe
cell_df <- data.frame(
  original_index = 1:n_cells,
  x = runif(n_cells, 0, 100),
  y = runif(n_cells, 0, 100),
  cluster = sample(1:k_cell_types, n_cells, replace = TRUE)
)
print(table(cell_df$cluster))

# Generate expression with realistic library size
expr_results <- generate_expression_profiles(
  cell_df = cell_df,
  n_genes = n_genes,
  k_cell_types = k_cell_types,
  marker_params = list(
    marker_genes_per_type = 10,
    marker_expression_fold = 1.5,
    marker_overlap_fold = 0.2
  ),
  spatial_params = list(
    range_spatial = 0.2,
    sill_spatial = 0.5,
    nugget_spatial = 0.1
  ),
  dropout_params = list(
    dropout_rate_min = 0.4,
    dropout_rate_max = 0.6,
    dropout_spatial_effect = 0.1,
    cell_type_effect = TRUE
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

# Analyze results
expr_matrix <- expr_results$expression_matrix
umi_counts <- Matrix::colSums(expr_matrix)

cat("\n   UMI count statistics:\n")
cat(sprintf("   - Mean: %.0f\n", mean(umi_counts)))
cat(sprintf("   - Median: %.0f\n", median(umi_counts)))
cat(sprintf("   - Min: %.0f\n", min(umi_counts)))
cat(sprintf("   - Max: %.0f\n", max(umi_counts)))
cat(sprintf("   - SD: %.0f\n", sd(umi_counts)))
cat(sprintf("   - CV: %.2f\n", sd(umi_counts)/mean(umi_counts)))

# Count expressed genes
genes_detected <- Matrix::rowSums(expr_matrix > 0)
never_detected <- sum(genes_detected == 0)
cat(sprintf("\n   Genes never detected: %d (%.1f%%)\n", never_detected, 100*never_detected/n_genes))

# Expected UMI count calculation
n_expressed_in_baseline <- sum(all_mu[1:n_genes] > -15)  # Use first cell type's mu
expected_umi <- 8000 * n_expressed_in_baseline / n_expressed_in_baseline  # Should be ~8000
cat(sprintf("\n   Expected UMI per cell (based on distribution): ~%.0f\n", expected_umi))

# Check if we're hitting the target ranges
target_median <- 5000
target_mean_range <- c(6000, 10000)
cat("\n   Target validation:\n")
cat(sprintf("   - Target median UMI: %d (Achieved: %.0f - %s)\n", 
            target_median, median(umi_counts), 
            ifelse(median(umi_counts) >= target_median * 0.8, "GOOD", "LOW")))
cat(sprintf("   - Target mean range: [%d, %d] (Achieved: %.0f - %s)\n", 
            target_mean_range[1], target_mean_range[2], mean(umi_counts),
            ifelse(mean(umi_counts) >= target_mean_range[1] && mean(umi_counts) <= target_mean_range[2], 
                   "GOOD", "OUTSIDE RANGE")))

# 3. Create visualization
cat("\n3. Creating visualization...\n")

# Plot mu distribution
p1 <- ggplot(data.frame(mu = all_mu[1:n_genes]), aes(x = mu)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(-15, -1, 1.5), linetype = "dashed", color = "red") +
  labs(title = "Gene Expression Distribution (mu values)",
       subtitle = "50% non-expressed, 35% low, 12% medium, 3% high",
       x = "mu (log-scale expression)", y = "Count") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# Plot UMI distribution
p2 <- ggplot(data.frame(umi = umi_counts), aes(x = umi)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  geom_vline(xintercept = mean(umi_counts), linetype = "solid", color = "red", size = 1) +
  geom_vline(xintercept = median(umi_counts), linetype = "dashed", color = "blue", size = 1) +
  labs(title = "UMI Count Distribution",
       subtitle = sprintf("Mean: %.0f (red), Median: %.0f (blue)", mean(umi_counts), median(umi_counts)),
       x = "UMI counts per cell", y = "Count") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# Save plots
cat("\n4. Saving plots...\n")
ggsave("plots/improved_distribution_mu.png", p1, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("plots/improved_distribution_umi.png", p2, width = 8, height = 6, dpi = 300, bg = "white")

cat("\n=== Test Complete ===\n")
cat("\nRecommendations:\n")
if (median(umi_counts) < target_median * 0.8) {
  cat("- Median UMI is still below target. Consider further adjustments:\n")
  cat("  * Reduce non-expressed genes to 40-45%\n")
  cat("  * Increase mu values by 0.5-1.0\n")
} else {
  cat("- Distribution parameters look good!\n")
  cat("- Proceed with full simulation\n")
}

# Final summary
cat("\nSummary:\n")
cat(sprintf("- Gene categories achieved: %.0f%% / %.0f%% / %.0f%% / %.0f%%\n",
            100*mu_zero/length(all_mu), 100*mu_low/length(all_mu), 
            100*mu_med/length(all_mu), 100*mu_high/length(all_mu)))
cat(sprintf("- UMI mean: %.0f, median: %.0f\n", mean(umi_counts), median(umi_counts)))
cat(sprintf("- Improvement from previous: mean +%.0f%%, median +%.0f%%\n",
            100*(mean(umi_counts) - 3821)/3821,  # Previous mean was 3821
            100*(median(umi_counts) - 393)/393))  # Previous median was 393