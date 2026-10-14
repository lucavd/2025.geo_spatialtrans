# tools/S1.3d_validity.R — S1.3d (correzione D-1 dell'analisi, vedi dev_notes): validità delle uscite geometriche dei candidati.
# F1 e F3 preparano la stessa geometria (geos_set_precision per cella, g = 1e-6): celle invalide contate una volta per entrambi.
# F2: still_invalid già in res; F4: n_invalid di f4_build() già in res (c5). Uso: Rscript --vanilla tools/S1.3d_validity.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow); library(jsonlite); library(parallel) })
source("tools/S1.3d_candidates.R")
O <- commandArgs(TRUE)[1]; if (is.na(O)) O <- "/mnt/micron/geo_spatialtrans/S1.3d"; GEOM <- file.path(O, "geom")
keys <- sub("_meta\\.json$", "", list.files(GEOM, pattern = "_meta\\.json$"))
R <- mclapply(keys, function(k) { m <- fromJSON(file.path(GEOM, paste0(k, "_meta.json")))
  t <- df_to_sfc(as.data.frame(read_parquet(file.path(GEOM, paste0(k, "_cells.parquet")))), m$n)
  s <- geos_set_precision(as_geos_geometry(t), 1e-6); v <- geos_is_valid(s)
  data.frame(key = k, n = m$n, n_invalid_F1F3 = sum(!v), n_empty_F1F3 = sum(geos_is_empty(s))) }, mc.cores = 40, mc.preschedule = FALSE)
R <- do.call(rbind, R); write.csv(R, file.path(O, "checks", "validity_F1F3.csv"), row.names = FALSE)
cat("casi", nrow(R), " invalide F1/F3:", sum(R$n_invalid_F1F3), " vuote:", sum(R$n_empty_F1F3), "\n")
