# Revisione avversariale S1.1 — 07b_perf_worst.R: prestazioni su mappe frammentate (non coperte da C11)
# Uso: Rscript --vanilla 07b_perf_worst.R <real_full_A6|null_full_A6|rand4_px1|rand4_px8>
.libPaths("/home/user/2025.geo_spatialtrans/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(sf))
setwd("/home/user/2025.geo_spatialtrans"); source("R/04b1_extract_regions.R")
id <- commandArgs(trailingOnly = TRUE)[1]
if (grepl("^rand4", id)) {
  set.seed(1); n <- 4000L
  clust <- data.frame(x = rep(seq_len(n), times = n), y = rep(seq_len(n), each = n), k = sample.int(4L, n * n, replace = TRUE))
  px <- if (id == "rand4_px8") 8 else 1; cc <- "k"
} else {
  o <- readRDS(file.path("/mnt/micron/geo_spatialtrans/S1.1/inputs", paste0(id, ".rds"))); clust <- o$clust; px <- o$pixel_size_um; cc <- o$cluster_col
}
gc(); t0 <- proc.time()[["elapsed"]]
r <- extract_regions(clust, px, cluster_col = cc, verbose = FALSE)
el <- proc.time()[["elapsed"]] - t0
cat(sprintf("PERF2,%s,%d,%d,%d,%.2f\n", id, nrow(clust), r$info$n_components, r$info$n_regions, el))

