# S2.2, G-3 (complemento): confronto del CONTENUTO dei file del pannello fra prima esecuzione e ripetizione (conteggi e cellule);
# i byte dei file differiscono perche' la ripetizione ha il campo metadata$totals_by aggiunto dopo la prima esecuzione.
# Uso: Rscript --vanilla tools/S2.2_repro_content.R -> appende a results/S2.2/S2.2_G3_repro.txt
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); suppressMessages(library(Matrix))
a <- "/mnt/micron/geo_spatialtrans/S2.2/panel"; b <- "/mnt/micron/geo_spatialtrans/S2.2/repro/panel"
fs <- list.files(a, pattern = "rds$"); d <- 0; md <- c(); dm <- 0
for (f in fs) { x <- readRDS(file.path(a, f)); y <- readRDS(file.path(b, f))
  if (!identical(x$counts, y$counts) || !identical(x$cells, y$cells)) d <- d + 1
  md <- union(md, setdiff(names(y$metadata), names(x$metadata)))
  if (!identical(x$metadata, y$metadata[names(x$metadata)])) dm <- dm + 1 }
msg <- sprintf("pannello, contenuto: %d file, %d diversi nei conteggi o nelle cellule; metadati comuni diversi in %d; campi nuovi nella ripetizione: %s",
               length(fs), d, dm, paste(md, collapse = ","))
cat(msg, "\n"); cat(msg, "\n", file = "results/S2.2/S2.2_G3_repro.txt", append = TRUE)
