# review/code/02b_voronoi_diag.R — diagnosi della differenza di B-1 (1 cella) fra obs_territories() e il Voronoi shapely di 02.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow) })
source("R/04b1_extract_regions.R"); source("R/04b3_tessellate_voronoi.R"); source("tools/R3_voronoi_metrics.R"); source("tools/S1.4_mutants.R"); source("tools/S1.4c_shape.R")
E <- dcg_env(); args <- commandArgs(TRUE); b <- args[1]; p <- strsplit(b, "_")[[1]]
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
cc <- read_parquet(sprintf("/mnt/micron/geo_spatialtrans/S1.4c_review/code/exp/cells_%s_%s.parquet", p[1], p[2])); n <- nrow(cc)
m <- as.data.frame(read_parquet(file.path("/mnt/micron/geo_spatialtrans/S1.4c/sim_out", paste0(b, ".parquet"))))
cx <- rep(NA_real_, n); cy <- cx; cx[m$label] <- m$cx; cy[m$label] <- m$cy
ob <- obs_territories(cx, cy, r3_roi_window(p[1], p[2], rois), E)
write_parquet(data.frame(a_t = ob$a_t, d = ob$d), sprintf("/mnt/micron/geo_spatialtrans/S1.4c_review/code/02b_%s.parquet", b))
cat("ok\n")
