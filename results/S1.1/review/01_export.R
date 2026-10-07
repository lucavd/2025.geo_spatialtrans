# Revisione avversariale S1.1 — 01_export.R
# Legge SOLO gli RDS di input (nessuna funzione del repo) e li esporta in parquet/json per Python.
.libPaths("/home/user/2025.geo_spatialtrans/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(arrow))
IN  <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
OUT <- "/mnt/micron/geo_spatialtrans/S1.1/review/export"
ids <- commandArgs(trailingOnly = TRUE)
for (id in ids) {
  o <- readRDS(file.path(IN, paste0(id, ".rds")))
  cl <- o$clust
  df <- data.frame(x = as.numeric(cl$x), y = as.numeric(cl$y), cluster = as.character(cl[[o$cluster_col]]))
  write_parquet(df, file.path(OUT, paste0(id, ".parquet")))
  if (!is.null(o$truth)) write.csv(o$truth, file.path(OUT, paste0(id, "_truth.csv")), row.names = FALSE)
  writeLines(jsonlite::toJSON(list(id = o$id, group = o$group, pixel_size_um = o$pixel_size_um,
             cluster_col = o$cluster_col, cluster_class = class(cl[[o$cluster_col]]),
             names = names(cl), meta = o$meta), auto_unbox = TRUE, null = "null"),
             file.path(OUT, paste0(id, ".json")))
  cat(id, nrow(df), "\n")
}

