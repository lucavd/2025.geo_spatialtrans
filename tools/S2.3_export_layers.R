# tools/S2.3_export_layers.R - S2.3: esporta regioni e nuclei delle 30 aree dello Step 1 (S1.6, U6) in WKB esadecimale per gli script Python
# (statistiche del pannello, docs/steps/S2.3.md). Uso (root del repo): Rscript --vanilla tools/S2.3_export_layers.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf) })
a <- commandArgs(TRUE); OUT <- if (length(a) >= 1) a[1] else "/mnt/micron/geo_spatialtrans/S2.3/layers"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
hex <- function(g) vapply(sf::st_as_binary(sf::st_set_crs(g, NA)), function(b) paste(as.character(b), collapse = ""), character(1))
for (p in sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE))) {
  tag <- sub("_U6\\.rds$", "", basename(p)); L <- readRDS(p)
  write.table(data.frame(region = seq_along(L$region_polygons), wkb = hex(L$region_polygons)), file.path(OUT, paste0(tag, "_regions.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cd <- L$cell_df
  write.table(data.frame(cell_id = cd$cell_id, cell_type = cd$cell_type, wkb = hex(L$cell_nuclei)), file.path(OUT, paste0(tag, "_nuclei.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(tag, length(L$region_polygons), nrow(cd), "\n")
}
