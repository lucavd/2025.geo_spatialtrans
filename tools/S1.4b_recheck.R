# tools/S1.4b_recheck.R — S1.4b: i grezzi ricalcolati con il codice corretto (D-4, D-5) sono identici agli ufficiali?
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
A <- "/mnt/micron/geo_spatialtrans/S1.4b"; B <- "/mnt/micron/geo_spatialtrans/S1.4b_recheck"
res <- do.call(rbind, lapply(c("prep", "loo"), function(st) { fs <- list.files(file.path(A, st))
  do.call(rbind, lapply(fs, function(f) { x <- readRDS(file.path(A, st, f)); y <- readRDS(file.path(B, st, f))
    data.frame(stage = st, file = f, identical = identical(x, y)) })) }))
write.csv(res, "results/S1.4b/S1.4b_recheck.csv", row.names = FALSE)
cat(sprintf("RECHECK identici: %d/%d\n", sum(res$identical), nrow(res)))
