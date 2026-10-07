# Revisione avversariale S1.1 — 04_export_outputs.R
# Legge gli output gia' salvati dal test (outputs/<id>.rds: r0_df = tutte le componenti, esatte) e li esporta.
.libPaths("/home/user/2025.geo_spatialtrans/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(arrow))
OUTD <- "/mnt/micron/geo_spatialtrans/S1.1/outputs"
EXP  <- "/mnt/micron/geo_spatialtrans/S1.1/review/export"
for (id in commandArgs(trailingOnly = TRUE)) {
  z <- readRDS(file.path(OUTD, paste0(id, ".rds")))
  d <- z$r0_df; d$cluster_id <- as.character(d$cluster_id)
  write_parquet(d[, c("region_id", "cluster_id", "n_px", "area_um2", "area_px_um2", "n_holes")], file.path(EXP, paste0(id, "_r0df.parquet")))
  cat(id, nrow(d), "\n")
}

