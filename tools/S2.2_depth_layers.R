# S2.2, regola R-3 (fattore di profondità): per ogni area del pannello dello Step 1 (S1.6, predefinito dei preset)
# numero di cellule per tipo e area di tessuto coperta dai territori. Uso: Rscript --vanilla tools/S2.2_depth_layers.R
# Scrive results/S2.2/S2.2_depth_layers.csv. Pre-registrazione: docs/steps/S2.2.md (regola R-3).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
f <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE))
out <- do.call(rbind, lapply(f, function(p) {
  x <- readRDS(p); cd <- x$cell_df
  tag <- sub("_U6\\.rds$", "", basename(p)); parts <- strsplit(tag, "_")[[1]]
  tb <- table(cd$cell_type)
  data.frame(archetype = parts[1], roi_id = parts[2], cell_type = names(tb), n_cells = as.integer(tb),
             territory_area_um2 = sum(cd$territory_area), stringsAsFactors = FALSE)
}))
dir.create("results/S2.2", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "results/S2.2/S2.2_depth_layers.csv", row.names = FALSE)
cat(nrow(out), "righe da", length(f), "aree\n")
