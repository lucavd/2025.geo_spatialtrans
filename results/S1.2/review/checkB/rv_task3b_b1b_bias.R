# rv_task3b_b1b_bias.R — revisione avversariale S1.2 / checkB, compito 3 (diagnosi).
# Perche' la densita' seed 42 di B1b supera totale x (1 - esclusa) in 29/30 ROI? Ipotesi: (H1) lambda usa l'area dei
# poligoni SEMPLIFICATI (simplify_tol_um = 0.5), non l'area in pixel; (H2) stesso seed in tutti i ROI = numeri casuali
# comuni (stessi uniformi per l'arrotondamento). Usa extract_regions()/seed_centroids() del progetto solo per generare
# l'output da verificare. 40 seed per ROI (42, 1..9 come il progetto, piu' 1001..1030).
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task3b_b1b_bias.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R")
IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
TOT <- c(A1 = 12419, A2 = 8924, A3 = 3100, A4 = 28096, A5 = 1185, A6 = 988)
rois <- read.csv("results/R2/R2_rois_checked.csv")
SEEDS <- c(42, 1:9, 1001:1030)
res <- mclapply(seq_len(nrow(rois)), function(i) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]
  o <- readRDS(file.path(IN, sprintf("roi_%s_%s.rds", A, roi)))
  reg <- extract_regions(o$clust, pixel_size_um = o$pixel_size_um, cluster_col = o$cluster_col, stride = 1L, verbose = FALSE)
  rd <- reg$region_df; at <- reg$info$area_total_um2
  comp <- data.frame(cluster_id = unique(as.character(rd$cluster_id)), cell_type = "c", fraction = 1)
  n <- vapply(SEEDS, function(s) seed_centroids(reg, data.frame(cell_type = "c", density = TOT[[A]]), comp, random_seed = s, verbose = FALSE)$info$n_cells, 0)
  lam <- TOT[[A]] * rd$area_um2 / 1e6
  data.frame(archetype = A, roi_id = roi, n_regions = nrow(rd), area_poly_over_px = sum(rd$area_um2) / sum(rd$area_px_um2),
             E_n_poly = sum(lam), E_n_pixel = TOT[[A]] * sum(rd$area_px_um2) / 1e6, sd_round = sqrt(sum((lam - floor(lam)) * (1 - (lam - floor(lam))))),
             n_seed42 = n[1], n_mean40 = mean(n), n_sd40 = sd(n), rank_seed42_of_10 = rank(n[1:10], ties.method = "max")[1],
             z_seed42 = (n[1] - sum(lam)) / sqrt(sum((lam - floor(lam)) * (1 - (lam - floor(lam))))), area_total_um2 = at)
}, mc.cores = 6L)
out <- do.call(rbind, res)
set.seed(42); u <- runif(5000)
out$mean_first_uniforms_seed42 <- NA
write.csv(out, "results/S1.2/review/checkB/rv_task3b_b1b_bias.csv", row.names = FALSE)
cat(sprintf("mean(runif) seed 42: primi 500 %.4f, primi 1000 %.4f, primi 3000 %.4f\n", mean(u[1:500]), mean(u[1:1000]), mean(u[1:3000])))
print(out[, c("archetype", "roi_id", "n_regions", "area_poly_over_px", "E_n_poly", "E_n_pixel", "sd_round", "n_seed42", "n_mean40", "n_sd40", "rank_seed42_of_10", "z_seed42")], digits = 5)
cat("media z seed42:", mean(out$z_seed42), "; media (n_mean40 - E_n_poly)/sd*sqrt(40):", mean((out$n_mean40 - out$E_n_poly) / out$sd_round * sqrt(40)), "\n")
