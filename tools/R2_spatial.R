#!/usr/bin/env Rscript
# tools/R2_spatial.R — statistiche spaziali della sessione R2 (eseguire: Rscript --vanilla tools/R2_spatial.R)
# - g(r) (pair correlation, spatstat::pcf) dei centroidi nucleari per ROI, finestra = maschera di tessuto (owin da immagine binaria)
# - envelope CSR 95% (39 simulazioni Poisson omogeneo con la stessa intensita') -> CP-R2.4
# - test di bimodalita' (Hartigan dip test) sulle aree nucleari (log) per ROI -> B-R2.4 (A5), calcolato su tutti gli archetipi
# Input: results/R2/R2_nuclei_all.parquet (metodo cellpose_rgb, scala nativa, keep), maschere tessuto da python (results/R2/tissue_masks/*.png, 1/4)
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({library(spatstat.geom); library(spatstat.explore); library(spatstat.random); library(arrow); library(diptest); library(png)})
set.seed(20260920)
nuc <- read_parquet("results/R2/R2_nuclei_all.parquet")
nuc <- nuc[nuc$method == "cellpose_rgb" & nuc$scale == 1 & nuc$keep, ]
rois <- read.csv("results/R2/R2_rois_checked.csv")
dir.create("results/R2/spatial", showWarnings = FALSE)
r_grid <- seq(0, 60, by = 0.5)
pcf_rows <- list(); dip_rows <- list()
for (i in seq_len(nrow(rois))) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; upp <- rois$um_per_px[i]; side_um <- rois$side_px[i] * upp
  d <- nuc[nuc$archetype == A & nuc$roi_id == roi, ]
  if (nrow(d) < 50) next
  # finestra: maschera tessuto (png 1/4, bianco = tessuto) -> owin
  m <- readPNG(sprintf("results/R2/tissue_masks/%s_%s_valid_ds4.png", A, roi))
  if (length(dim(m)) == 3) m <- m[, , 1]
  m <- m[nrow(m):1, ] > 0.5                      # riga 1 = y=0 in alto -> owin vuole y crescente verso l'alto
  W <- owin(mask = m, xrange = c(0, side_um), yrange = c(0, side_um))   # righe = y crescente verso l'alto (dopo il flip), colonne = x; NIENTE t() (bug trovato dal revisore R2)
  X <- ppp(d$x_um, side_um - d$y_um, window = W)  # y invertita coerentemente
  n_lost <- nrow(d) - X$n
  if (n_lost / nrow(d) > 0.01) stop(sprintf("%s_%s: %d nuclei keep fuori dalla finestra (%.1f %%): finestra incoerente", A, roi, n_lost, 100 * n_lost / nrow(d)))
  lam <- intensity(X)
  g <- pcf(X, r = r_grid, correction = "translate", divisor = "d")
  E <- envelope(X, pcf, r = r_grid, nsim = 39, correction = "translate", divisor = "d", verbose = FALSE, savefuns = FALSE)
  pcf_rows[[length(pcf_rows) + 1]] <- data.frame(archetype = A, roi_id = roi, r_um = g$r, g_obs = g$trans, g_lo = E$lo, g_hi = E$hi, n = X$n, intensity_per_mm2 = lam * 1e6, area_mm2 = area(W) / 1e6)
  # raggio di inibizione: primo r in cui g_obs supera g_lo (esce dal deficit)
  la <- log(d$area_um2); dp <- dip.test(la)
  nn <- nndist(X)
  dip_rows[[length(dip_rows) + 1]] <- data.frame(archetype = A, roi_id = roi, n = X$n, n_keep = nrow(d), n_lost_window = n_lost, intensity_per_mm2 = lam * 1e6,
      r_inhib_um = suppressWarnings(min(g$r[g$r > 0 & g$trans >= E$lo], na.rm = TRUE)),
      g_at_3um = g$trans[which.min(abs(g$r - 3))], g_at_5um = g$trans[which.min(abs(g$r - 5))], g_at_10um = g$trans[which.min(abs(g$r - 10))], g_at_30um = g$trans[which.min(abs(g$r - 30))],
      frac_r_below_lo_lt5 = mean(g$trans[g$r > 0 & g$r <= 5] < E$lo[g$r > 0 & g$r <= 5], na.rm = TRUE),
      nn_median_um = median(nn), nn_median_csr_um = 0.5 / sqrt(lam), clark_evans_R = median(nn) / (0.5 / sqrt(lam)),
      dip_stat = unname(dp$statistic), dip_p = dp$p.value)
  cat(sprintf("%s_%s n=%d lambda=%.0f/mm2 g(3)=%.2f g(5)=%.2f g(30)=%.2f dip p=%.3g\n", A, roi, X$n, lam * 1e6, g$trans[which.min(abs(g$r - 3))], g$trans[which.min(abs(g$r - 5))], g$trans[which.min(abs(g$r - 30))], dp$p.value))
}
write.csv(do.call(rbind, pcf_rows), "results/R2/spatial/R2_pcf.csv", row.names = FALSE)
write.csv(do.call(rbind, dip_rows), "results/R2/spatial/R2_spatial_summary.csv", row.names = FALSE)
