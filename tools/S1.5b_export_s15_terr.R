# tools/S1.5b_export_s15_terr.R - S1.5b: territori delle celle simulate in S1.5 (seme 1, braccio preserve_area) per le stime a mano degli attesi.
# Sola lettura dei layer di S1.5; regola d'interna di R3/S1.5 (invariata). Uso: Rscript --vanilla tools/S1.5b_export_s15_terr.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({library(sf); library(geos); library(arrow)})
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R")
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
out <- list()
for (k in seq_len(nrow(rois))) { A <- rois$archetype[k]; r <- rois$roi_id[k]
  L <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.5/layers/%s_%s_preserve_area.rds", A, r))
  w <- r3_roi_window(A, r, rois); td <- L$territory_df
  in_frame <- geos::geos_within(geos::as_geos_geometry(L$tiles), geos::as_geos_geometry(w$frame))
  it <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
  out[[k]] <- data.frame(archetype = A, roi_id = r, a_t = L$cell_df$territory_area, a_n = L$cell_df$nucleus_area, interior = it) }
dir.create("/mnt/micron/geo_spatialtrans/S1.5b", showWarnings = FALSE)
write_parquet(do.call(rbind, out), "/mnt/micron/geo_spatialtrans/S1.5b/s15_seed1_terr.parquet"); cat("ok\n")
