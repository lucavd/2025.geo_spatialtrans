# tools/F3.3_regions.R - F3.3: regioni delle aree dalla maschera valida a 1/4 (PNG di tools/F3.3_segment.py), con extract_regions()
# chiamata come i livelli di S1.6 (place_cells con min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1; finestra come
# tools/R3_voronoi_metrics.R: px = lato in um / colonne della maschera). WKB esadecimale come tools/S2.3_export_layers.R.
# Uso (root del repo): Rscript --vanilla tools/F3.3_regions.R <area,area,...|IN|holdout>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE); A <- read.csv("results/F3.3/F3.3_areas.csv", stringsAsFactors = FALSE)
A <- if (a[1] == "IN") A[A$level == "IN", ] else if (a[1] == "holdout") A[A$level != "IN", ] else A[A$area %in% strsplit(a[1], ",")[[1]], ]
hex <- function(g) vapply(sf::st_as_binary(sf::st_set_crs(g, NA)), function(b) paste(as.character(b), collapse = ""), character(1))
for (k in seq_len(nrow(A))) {
  o <- file.path(Sys.getenv("F33_ROOT", "/mnt/micron/geo_spatialtrans/F3.3"), "seg", A$area[k])
  m <- png::readPNG(file.path(o, "valid_ds4.png")); if (length(dim(m)) == 3) m <- m[, , 1]
  m <- m > 0.5; stopifnot(nrow(m) == ncol(m)); px <- A$side_um[k] / ncol(m)
  idx <- which(m, arr.ind = TRUE); clust <- data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L)
  reg <- extract_regions(clust, pixel_size_um = px, min_region_area_um2 = 0, cluster_col = "intensity_cluster", simplify_tol_um = 0, stride = 1L, verbose = FALSE)
  write.table(data.frame(region = seq_along(reg$region_polygons), wkb = hex(reg$region_polygons)), file.path(o, "regions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat(A$area[k], length(reg$region_polygons), round(sum(as.numeric(sf::st_area(reg$region_polygons)))), "\n")
}
