# review/code/08_at_vs_geometry.R — a_t (territory_df di tessellate_voronoi, usato per N/C, tetto di K-2, B-1) contro l'area
# della geometria del territorio (usata per il ritaglio K-1/K-2 e per d_gen), nei 30 tassellamenti reali di S1.4.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
rois <- read.csv("results/R2/R2_rois_checked.csv")
prim <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
out <- list()
for (i in seq_len(nrow(rois))) { A <- rois$archetype[i]; roi <- rois$roi_id[i]
  r <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.4/real/%s_%s_%s.rds", A, roi, prim[[A]])); cc <- r$cells
  g <- as_geos_geometry(r$terr); ag <- geos_area(g); np <- geos_num_geometries(g)
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  rel <- abs(ag - cc$a_t) / cc$a_t; bad <- rel > 1e-6
  out[[i]] <- data.frame(archetype = A, roi = roi, n = nrow(cc), n_multipart = sum(np > 1), n_bad = sum(bad), n_bad_ev = sum(bad & ev),
                         n_bad_interior = sum(bad & cc$interior), max_rel = max(rel), sum_at = sum(cc$a_t), sum_geom = sum(ag),
                         geom_lt_at = sum(ag < cc$a_t * (1 - 1e-6)), geom_gt_at = sum(ag > cc$a_t * (1 + 1e-6)))
}
d <- do.call(rbind, out); write.csv(d, "results/S1.4c/review/code/08_at_vs_geometry.csv", row.names = FALSE); print(d)
