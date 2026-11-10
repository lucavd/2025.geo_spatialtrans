# tools/S2.6_fix_cell.R - S2.6, deviazione dichiarata: nel livello della sezione del linfonodo il territorio della cellula 624286 ha uno
# spigolo di 8e-5 um su cui la triangolazione vincolata di GEOS fallisce (place_and_move_transcripts, .s23_triangles). Si semplifica
# quel solo territorio con tolleranza 1e-4 um (area +0.00033 um2, nucleo invariato); l'originale resta in <file>_orig.rds.
# Uso: Rscript --vanilla tools/S2.6_fix_cell.R <layer.rds> <cell_id> <tol_um>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE); f <- a[1]; id <- as.integer(a[2]); tol <- as.numeric(a[3])
fo <- sub("\\.rds$", "_orig.rds", f); if (!file.exists(fo)) file.copy(f, fo)
L <- readRDS(fo); i <- which(L$cell_df$cell_id == id); stopifnot(length(i) == 1L)
T <- geos::as_geos_geometry(sf::st_set_crs(L$cell_territories[i], NA)); T2 <- geos::geos_simplify_preserve_topology(T, tol)
N <- geos::as_geos_geometry(sf::st_set_crs(L$cell_nuclei[i], NA)); stopifnot(geos::geos_within(N, geos::geos_buffer(T2, 1e-9)))
stopifnot(!inherits(try(.s23_triangles(geos::geos_difference(T2, N)), silent = TRUE), "try-error"))
L$cell_territories[i] <- sf::st_as_sfc(geos::geos_write_wkb(T2))
L$metadata$s26_fixed_cells <- data.frame(cell_id = id, tol_um = tol, d_area_um2 = geos::geos_area(T2) - geos::geos_area(T))
saveRDS(L, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
cat(sprintf("cellula %d: territorio semplificato (%d -> %d vertici, area %+.6f um2)\n", id, geos::geos_num_coordinates(T), geos::geos_num_coordinates(T2), geos::geos_area(T2) - geos::geos_area(T)))
