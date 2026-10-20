# tools/R3b_posthoc_sdarea.R — R3b, analisi POST HOC di sensibilita' (dopo tools/R3b_posthoc_area.R): maschera CO per centri e
# territori, ma modello di taglia M2 (per fold, sugli altri 4 ROI) e tetto phi (q99 sui 5 ROI) stimati sui soli nuclei di forma
# StarDist; poi nucleo ritagliato e gonfiato con gli stessi semi di P-3. Non cambia previsioni, verdetti ne' regola.
# Uscite: /mnt/micron/geo_spatialtrans/R3b/ph_sdarea/{loo,gen}/, gen_in_consensus_sdarea/; results/R3b/R3b_posthoc_sdarea_{loo,gen}.csv
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 Rscript --vanilla tools/R3b_posthoc_sdarea.R posthoc   (argomento obbligatorio: R3b_bench.R lo legge come STAGE)
source("tools/R3b_bench.R")                     # STAGE non riconosciuto: definisce solo funzioni e costanti
PH <- file.path(OUT, "ph_sdarea"); for (d in c("loo", "gen")) dir.create(file.path(PH, d), recursive = TRUE, showWarnings = FALSE)
real_sd <- function(A, roi) { z <- load_real(A, roi, "consensus")
  g <- as.data.frame(read_parquet(rp(A, roi, "consensus", "gen"))); mp <- as.data.frame(read_parquet(file.path(OUT, "masks", sprintf("%s_%s_consensus_map.parquet", A, roi))))
  stopifnot(identical(z$cells$idx, g$idx)); z$cells$matched <- mp$src[match(g$label, mp$label)] == "sd"; z }
phi_sd <- function(A) { nc <- unlist(lapply(ids_of(A), function(r) { cc <- real_sd(A, r)$cells; e <- evset(cc); cc$a_n[e] / cc$a_t[e] }))
  min(1, round(unname(quantile(nc, 0.99, type = 7)), 3)) }
ph_one <- function(A, k) {
  f <- file.path(PH, "gen", sprintf("%s_fold%d.rds", A, k)); if (file.exists(f)) return(invisible(NULL))
  ids <- ids_of(A); ai <- match(A, paste0("A", 1:6)); rs <- lapply(ids, function(r) real_sd(A, r)); fs <- fit_set(rs[-k], "full")
  seed0 <- BASE14 + 1000 * ai + 10 * k; phi <- phi_sd(A); f2 <- c(fs$fit, a_ref = fs$a_ref)
  saveRDS(data.frame(archetype = A, fold = k, roi_test = ids[k], seed = seed0, t(f2), phi = phi, n_ev_train = fs$n_ev), file.path(PH, "loo", sprintf("%s_fold%d.rds", A, k)))
  te <- load_real(A, ids[k], "consensus"); inp <- inputs_of(te); ct <- ct_m2(f2, phi)
  raster <- any(WIN$archetype == A & WIN$roi_id == ids[k])
  gdir <- file.path(OUT, "gen_in_consensus_sdarea"); dir.create(gdir, showWarnings = FALSE); rows <- list()
  for (s in seq_len(NSEED)) { sd <- seed0 + 100000 * s
    o <- E$derive_cell_geometry(inp$cen, inp$tess, ct, random_seed = sd, verbose = FALSE); cd <- o$cell_df
    rows[[s]] <- data.frame(archetype = A, fold = k, roi = ids[k], seed_idx = s, seed = sd, phi = phi,
                            trunc_all = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target))
    if (raster) write_nuc(file.path(gdir, sprintf("%s_%s_fold%d_s%02d_K2.parquet", A, ids[k], k, s)), as_geos_geometry(o$cell_nuclei)) }
  saveRDS(do.call(rbind, rows), f); invisible(NULL)
}
par_run(expand.grid(A = TIS, k = 1:5, stringsAsFactors = FALSE), ph_one, "ph_sdarea")
write.csv(do.call(rbind, lapply(sort(list.files(file.path(PH, "loo"), full.names = TRUE)), readRDS)), "results/R3b/R3b_posthoc_sdarea_loo.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(sort(list.files(file.path(PH, "gen"), full.names = TRUE)), readRDS)), "results/R3b/R3b_posthoc_sdarea_gen.csv", row.names = FALSE)
