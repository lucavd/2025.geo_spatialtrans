# tools/S1.2_perf.R — C-perf: 4 000 × 4 000 µm a densita' A4 (28 096/mm², d regola) in processo separato
# Uso: /usr/bin/time -v Rscript --vanilla tools/S1.2_perf.R  (stampa una riga PERF,...)
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(sf))
source("R/04b2_seed_centroids.R")
s <- 4000
reg <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = s^2),
            region_polygons = st_sfc(st_polygon(list(rbind(c(0, 0), c(s, 0), c(s, s), c(0, s), c(0, 0))))))
t <- system.time(o <- seed_centroids(reg, data.frame(cell_type = "a4", density = 28096),
                                     data.frame(cluster_id = 1, cell_type = "a4", fraction = 1), random_seed = 42))[["elapsed"]]
cat(sprintf("PERF,%d,%d,%d,%.1f\n", o$info$n_target, o$info$n_cells, o$info$n_attempts, t))
