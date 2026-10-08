# tools/S1.3_perf.R — C-perf e curva di scala (C-10c): un quadrato a densita' A4 (28 096/mm²), RSA regola, seed 42.
# Uso: /usr/bin/time -v Rscript --vanilla tools/S1.3_perf.R <backend> <n>   (stampa una riga PERF,...)
# n = 449536 → lato 4 000 µm = pattern C-perf di S1.2; altrimenti lato = 4000·sqrt(n/449536).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(sf))
source("R/04b2_seed_centroids.R"); source("tools/S1.3_mutants.R"); source("tools/S1.3_variants.R")
a <- commandArgs(TRUE); backend <- a[1]; n <- as.numeric(a[2])
s <- 4000 * sqrt(n / 449536)
reg <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = s^2),
            region_polygons = st_sfc(st_polygon(list(rbind(c(0, 0), c(s, 0), c(s, s), c(0, s), c(0, 0))))))
o <- seed_centroids(reg, data.frame(cell_type = "a4", density = 28096), data.frame(cluster_id = 1, cell_type = "a4", fraction = 1),
                    random_seed = 42, verbose = FALSE)
E <- tv_variant(tv_env(), backend); stopifnot(E$.tv_engine() == backend)
t <- system.time(tv <- E$tessellate_voronoi(o$centroids, reg, verbose = FALSE))[["elapsed"]]
c2 <- abs(sum(tv$territory_df$territory_area) - s^2) / s^2
cat(sprintf("PERF,%s,%d,%.1f,%.3e,%d,%d\n", backend, nrow(o$centroids), t, c2, tv$info$n_multipart, tv$info$n_snapped))
