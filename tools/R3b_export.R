# tools/R3b_export.R — R3b: esporta in CSV (results/R3b/) le uscite RDS del banco: righe del nucleo simulato (P-4, invarianti),
# righe di CP-1, parametri M2 per fold, phi, regressione G-3 (parte R). Uso: Rscript --vanilla tools/R3b_export.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
O <- "/mnt/micron/geo_spatialtrans/R3b"; R <- "results/R3b"; dir.create(R, showWarnings = FALSE)
g <- do.call(rbind, lapply(sort(list.files(file.path(O, "gen"), full.names = TRUE)), readRDS))
write.csv(g, file.path(R, "R3b_gen_rows.csv"), row.names = FALSE)
cp <- lapply(sort(list.files(file.path(O, "cp1"), full.names = TRUE)), readRDS)
write.csv(do.call(rbind, lapply(cp, function(z) z$rows)), file.path(R, "R3b_cp1_rows.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(cp, function(z) data.frame(win_id = z$win_id, n = z$n, n_co_kept = z$n_co_kept, n_co_removed = z$n_co_removed,
  n_manual = z$n_manual, n_manual_in_mask = z$n_manual_in_mask))), file.path(R, "R3b_cp1_inputs.csv"), row.names = FALSE)
l <- do.call(rbind, lapply(sort(list.files(file.path(O, "S1.4", "loo"), full.names = TRUE)), function(f) { z <- readRDS(f)
  data.frame(archetype = z$archetype, method = z$method, fold = z$fold, roi_test = z$roi_test, seed = z$seed, t(z$fits$M2),
             n_ev_train = z$n_ev_train, n_ev_test = z$n_ev_test) }))
write.csv(l, file.path(R, "R3b_loo_M2.csv"), row.names = FALSE)
file.copy(file.path(O, "phi.csv"), file.path(R, "R3b_phi.csv"), overwrite = TRUE)
for (f in c("G3_R.csv", "G3_phi_sr.csv")) file.copy(file.path(O, "regress", f), file.path(R, paste0("R3b_", f)), overwrite = TRUE)
cat("export ok:", nrow(g), "righe gen,", length(cp), "finestre CP-1,", nrow(l), "fold\n")
