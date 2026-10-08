# tools/R3_posthoc_P4.R — R3 POST HOC P4: sensibilita' al ritaglio. La maschera valida R2 = tessuto ∧ ¬bolla ha
# 34–221 buchi per ROI in A1, A5, A6 con area mediana ~258 µm² (= disco minimo della bubble_mask dopo dilatazione 5 µm +
# margine 3 µm): falsi positivi della maschera bolle (BL-026/BL-040) che rendono "ritagliati" molti territori.
# Qui: ritaglio sulla sola maschera di tessuto (tissue_ds4 di R2, esportata da tools/R3_export_tissue.py), segmentatore
# primario, stesse metriche. Uscita: results/R3/R3_posthoc_P4_tissue_only.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel); library(arrow); library(spatstat.geom); library(spatstat.explore) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R")
OUT <- "/mnt/micron/geo_spatialtrans/R3"; MDIR <- file.path(OUT, "posthoc", "masks_tissue")
PRIM <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
nuc <- as.data.frame(read_parquet("results/R2/R2_nuclei_all.parquet", col_select = c("archetype", "roi_id", "method", "scale", "keep", "label", "x_um", "y_um", "area_um2")))
nuc <- nuc[nuc$scale == 1 & nuc$keep, ]
res <- mclapply(seq_len(nrow(rois)), function(i) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; side_um <- rois$side_um[i]
  m <- png::readPNG(file.path(MDIR, sprintf("%s_%s_tissue_ds4.png", A, roi))); if (length(dim(m)) == 3) m <- m[, , 1]; m <- m > 0.5
  px <- side_um / ncol(m); idx <- which(m, arr.ind = TRUE)
  reg <- extract_regions(data.frame(x = idx[, 2], y = idx[, 1], cl = 1L), pixel_size_um = px, min_region_area_um2 = 0, cluster_col = "cl",
                         simplify_tol_um = 0, stride = 1L, verbose = FALSE)
  w <- list(m = m, px = px, side_um = side_um, poly = st_union(reg$region_polygons),
            frame = st_sfc(st_polygon(list(rbind(c(0, 0), c(side_um, 0), c(side_um, side_um), c(0, side_um), c(0, 0))))))
  d <- nuc[nuc$archetype == A & nuc$roi_id == roi & nuc$method == PRIM[[A]], ]; d <- d[r3_in_mask(d$x_um, d$y_um, m, px), ]
  tv <- r3_tessellate(d$x_um, d$y_um, w$poly, w$frame); li <- r3_local_intensity(d$x_um, d$y_um, w)
  s <- r3_summary(tv, li$lambda_loc); ii <- tv$interior
  data.frame(archetype = A, roi_id = roi, n = nrow(d), conserv = abs(sum(tv$area) / as.numeric(st_area(w$poly)) - 1), t(s),
             nc_median = median(d$area_um2[ii] / tv$area[ii]), median_eq_r_all = median(sqrt(tv$area / pi)))
}, mc.cores = 30L)
bad <- vapply(res, function(z) !is.data.frame(z), TRUE); if (any(bad)) stop(res[bad][[1]])
out <- do.call(rbind, res); write.csv(out, "results/R3/R3_posthoc_P4_tissue_only.csv", row.names = FALSE)
print(aggregate(cbind(frac_interior, median_eq_r, cv_loc, median_ecc_T, nc_median) ~ archetype, out, median))
