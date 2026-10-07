# S1.1 — C11: prestazioni di extract_regions() su griglie piene 4000x4000
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent /usr/bin/time -v Rscript --vanilla tools/S1.1_perf.R <id>
# Stampa una riga CSV: id,n_px,n_components,n_regions,elapsed_extract_s,rss_after_load_mb
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(sf))
source("R/04b1_extract_regions.R")
id <- commandArgs(trailingOnly = TRUE)[1]
o <- readRDS(file.path("/mnt/micron/geo_spatialtrans/S1.1/inputs", paste0(id, ".rds")))
rss <- function() as.numeric(sub("kB", "", grep("^VmRSS", readLines("/proc/self/status"), value = TRUE) |> sub(pattern = "VmRSS:", replacement = ""))) / 1024
gc(); r_load <- rss()
t0 <- proc.time()[["elapsed"]]
r <- extract_regions(o$clust, o$pixel_size_um, cluster_col = o$cluster_col)   # parametri d'uso (min 100, tol 0.5)
el <- proc.time()[["elapsed"]] - t0
cat(sprintf("PERF,%s,%d,%d,%d,%.2f,%.0f\n", id, nrow(o$clust), r$info$n_components, r$info$n_regions, el, r_load))
