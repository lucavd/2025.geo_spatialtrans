
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
L <- readRDS("/mnt/micron/geo_spatialtrans/S2.6/sections/sec_A4_s34235060_orig.rds")
i <- which(L$cell_df$cell_id == 624286L); stopifnot(length(i) == 1L)
fx <- list(cell_id = 624286L, territory = sf::st_set_crs(L$cell_territories[i], NA), nucleus = sf::st_set_crs(L$cell_nuclei[i], NA),
           source = "/mnt/micron/geo_spatialtrans/S2.6/sections/sec_A4_s34235060_orig.rds")
dir.create("R/testing/fixtures", showWarnings = FALSE)
saveRDS(fx, "R/testing/fixtures/S2.7_cell624286.rds")
cat("fixture ok", geos::geos_num_coordinates(geos::as_geos_geometry(fx$territory)), "\n")
