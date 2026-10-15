.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; O <- file.path(D, "review_recalc", "v0viol"); dir.create(O, showWarnings = FALSE)
for (k in readLines(file.path(D, "keys_K.txt"))) { f <- file.path(D, "res", paste0(k, c("_o12.rds", "_o1.rds"))); f <- f[file.exists(f)][1]; v <- readRDS(f)$viol$o1_V0
  if (is.null(v)) next; v <- v[v$what %in% c("pair_ij", "pair_ji"), ]; p <- unique(data.frame(i = pmin(v$i, v$j), j = pmax(v$i, v$j))); write.csv(p, file.path(O, paste0(k, ".csv")), row.names = FALSE) }
cat(length(list.files(O)), "file\n")