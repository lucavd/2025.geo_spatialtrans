# tools/R6_export_terr.R - export real territories (S1.4/real, R3b/S1.4/real) to parquet for R6: idx, x, y, a_t, a_n, interior, matched, ecc_n
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(arrow))
M <- "/mnt/micron/geo_spatialtrans"; dir.create(file.path(M, "R6", "terr"), showWarnings = FALSE, recursive = TRUE)
rois <- read.csv("results/R2/R2_rois_checked.csv")
src <- data.frame(archetype = paste0("A", 1:6), method = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "consensus", "cellpose_rgb", "consensus"),
                  base = file.path(M, c("S1.4", "S1.4", "S1.4", "R3b/S1.4", "S1.4", "R3b/S1.4")), stringsAsFactors = FALSE)
for (k in seq_len(nrow(src))) { A <- src$archetype[k]
  for (r in rois$roi_id[rois$archetype == A]) {
    z <- readRDS(file.path(src$base[k], "real", sprintf("%s_%s_%s.rds", A, r, src$method[k]))); cc <- as.data.frame(z$cells)
    cc <- cc[, intersect(c("idx", "x", "y", "a_t", "a_n", "interior", "matched", "ecc_n"), names(cc))]
    write_parquet(cc, file.path(M, "R6", "terr", sprintf("%s_%s.parquet", A, r))); cat(A, r, nrow(cc), "\n") } }
