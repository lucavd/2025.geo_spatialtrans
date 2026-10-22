# tools/S1.5b_export_layers.R - S1.5b: celle del seme 1 (bracci per tipo; braccio a un tipo da S1.5 con le etichette del nullo) per
# figure e descrittive. Uscita: /mnt/micron/geo_spatialtrans/S1.5b/cells_seed1.parquet
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R")
M <- "/mnt/micron/geo_spatialtrans"
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
null_types <- function(A, n, seed) { sp <- tissue_preset_spec(PRESET[[A]], types = "multi", size_arm = "group")
  set.seed(seed + 7L); sample(sp$composition$cell_type, n, replace = TRUE, prob = sp$composition$fraction) }   # come tools/S1.5b_bench.R
out <- list()
for (k in seq_len(nrow(rois))) { A <- rois$archetype[k]; r <- rois$roi_id[k]; w <- r3_roi_window(A, r, rois)
  for (arm in c("group", "recentred", "single")) {
    f <- if (arm == "single") file.path(M, "S1.5/layers", sprintf("%s_%s_preserve_area.rds", A, r)) else file.path(M, "S1.5b/layers", sprintf("%s_%s_%s.rds", A, r, arm))
    L <- readRDS(f); cd <- L$cell_df; td <- L$territory_df
    in_frame <- geos::geos_within(geos::as_geos_geometry(L$tiles), geos::as_geos_geometry(w$frame))
    it <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
    ty <- if (arm == "single") null_types(A, nrow(cd), 20261021 + 1000 * k + 1) else as.character(cd$cell_type)
    out[[length(out) + 1]] <- data.frame(archetype = A, roi_id = r, arm = arm, x = cd$x, y = cd$y, cell_type = ty, a_n = cd$nucleus_area,
                                         a_t = cd$territory_area, interior = it) } }
write_parquet(do.call(rbind, out), file.path(M, "S1.5b/cells_seed1.parquet")); cat("ok\n")
