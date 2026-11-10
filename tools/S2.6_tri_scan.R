# tools/S2.6_tri_scan.R - S2.6: cellule su cui la triangolazione vincolata di GEOS (.s23_triangles di place_and_move_transcripts) fallisce.
# Uso: Rscript --vanilla tools/S2.6_tri_scan.R <layer.rds> <out.csv>. Blocchi di 2000 cellule; nei blocchi che falliscono, cellula per cellula.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE); L <- readRDS(a[1])
Tt <- geos::as_geos_geometry(sf::st_set_crs(L$cell_territories, NA)); Nn <- geos::as_geos_geometry(sf::st_set_crs(L$cell_nuclei, NA))
Cy <- geos::geos_difference(Tt, Nn)
bad <- list()
for (part in c("nucleus", "cytoplasm")) {
  g <- if (part == "nucleus") Nn else Cy
  ch <- split(seq_along(g), ceiling(seq_along(g) / 2000))
  for (k in ch) {
    if (!inherits(try(.s23_triangles(g[k]), silent = TRUE), "try-error")) next
    for (i in k) { e <- try(.s23_triangles(g[i]), silent = TRUE)
      if (inherits(e, "try-error")) bad[[length(bad) + 1]] <- data.frame(cell_index = i, cell_id = L$cell_df$cell_id[i], part = part,
        n_vertices = geos::geos_num_coordinates(g[i]), area_um2 = geos::geos_area(g[i]), type = geos::geos_type(g[i]), valid = geos::geos_is_valid(g[i]),
        msg = substr(gsub("\n", " ", as.character(e)), 1, 120)) }
  }
}
B <- if (length(bad)) do.call(rbind, bad) else data.frame()
write.csv(B, a[2], row.names = FALSE)
cat(sprintf("%s: %d cellule, %d geometrie che falliscono\n", basename(a[1]), nrow(L$cell_df), nrow(B)))
