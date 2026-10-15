# tools/S1.4_perf.R — S1.4 C-perf: derive_cell_geometry() sulla sezione S (ricetta di S1.3d: RSA, 28 096 /mm², seme 42,
# 4 000 x 4 000 µm, 449 536 celle). Uso (root del repo), un processo per variante, memoria con /usr/bin/time -v:
#   S14_STEP=prep  Rscript --vanilla tools/S1.4_perf.R            # tassella e salva S (non cronometrato)
#   S14_STEP=circle|ellipse /usr/bin/time -v Rscript --vanilla tools/S1.4_perf.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
STEP <- Sys.getenv("S14_STEP"); F <- "/mnt/micron/geo_spatialtrans/S1.4/perf_S.rds"
dir.create(dirname(F), showWarnings = FALSE, recursive = TRUE)
if (STEP == "prep") {
  source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
  s <- 4000; sqp <- st_polygon(list(rbind(c(0, 0), c(s, 0), c(s, s), c(0, s), c(0, 0))))
  reg <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = s^2), region_polygons = st_sfc(sqp))
  o <- seed_centroids(reg, data.frame(cell_type = "a4", density = 28096), data.frame(cluster_id = 1, cell_type = "a4", fraction = 1), random_seed = 42, verbose = FALSE)
  tv <- tessellate_voronoi(o$centroids, reg, verbose = FALSE); stopifnot(nrow(o$centroids) == 449536)
  saveRDS(list(cen = o$centroids, tv = tv[c("cell_territories", "territory_df")]), F); cat("prep ok", nrow(o$centroids), "\n")
} else {
  source("R/04b4_derive_cell_geometry.R")
  S <- readRDS(F); cen <- S$cen; cen$cell_type <- "a4"
  ct <- data.frame(cell_type = "a4", density = 28096, nucleus_to_eq_ratio = 0.6, nucleus_area_slope = 0.5, nucleus_area_sdlog = 0.3,
                   nucleus_aspect = if (STEP == "ellipse") 1.6 else 1, nucleus_orient = if (STEP == "ellipse") "territory" else "random")
  gc(reset = TRUE); t0 <- proc.time()[["elapsed"]]
  o <- derive_cell_geometry(cen, S$tv, ct, random_seed = 42, verbose = FALSE)
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("PERF %s n=%d elapsed_s=%.1f clamped=%.4f not_within=%d\n", STEP, o$info$n, el, o$info$frac_clamped, o$info$n_not_within))
}
