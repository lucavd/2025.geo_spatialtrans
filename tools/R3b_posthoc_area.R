# tools/R3b_posthoc_area.R — R3b, analisi POST HOC (scritta dopo i verdetti): il 27-29 % dei nuclei della maschera CO ha la forma
# Space Ranger (nessun StarDist corrispondente) e area mediana piu' piccola dei nuclei StarDist. Quanto cambiano modello di taglia
# M2 (rho, beta, sigma) e phi se si stimano solo sui nuclei di forma StarDist (territori CO invariati)? Quanta area hanno perso i
# poligoni Space Ranger sotto i pixel StarDist? Non cambia previsioni, verdetti ne' regola.
# Uscita: results/R3b/R3b_posthoc_area.csv.  Uso: Rscript --vanilla tools/R3b_posthoc_area.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(sf); library(geos) })
source("tools/S1.4d_mutants.R"); E <- dcg_env_d()
O <- "/mnt/micron/geo_spatialtrans/R3b"
rois <- read.csv("results/R2/R2_rois_checked.csv")
out <- do.call(rbind, lapply(c("A4", "A6"), function(A) {
  ids <- rois$roi_id[rois$archetype == A]
  cc <- do.call(rbind, lapply(ids, function(r) {
    z <- readRDS(file.path(O, "S1.4", "real", sprintf("%s_%s_consensus.rds", A, r)))
    g <- as.data.frame(read_parquet(file.path(O, "R3", "real", sprintf("%s_%s_consensus_gen.parquet", A, r))))
    mp <- as.data.frame(read_parquet(file.path(O, "masks", sprintf("%s_%s_consensus_map.parquet", A, r))))
    stopifnot(identical(z$cells$idx, g$idx))
    d <- z$cells; d$src <- mp$src[match(g$label, mp$label)]; d$trim <- 1 - mp$n_px[match(g$label, mp$label)] / mp$n_px_sr[match(g$label, mp$label)]
    d$roi <- r; d$n_roi <- z$n; d }))
  a_ref <- sum(tapply(cc$a_t, cc$roi, sum)) / sum(tapply(cc$n_roi, cc$roi, `[`, 1))
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  f_all <- E$fit_nucleus_size_model(cc$a_n[ev], cc$a_t[ev], a_ref, "full")
  e2 <- ev & cc$src == "sd"; f_sd <- E$fit_nucleus_size_model(cc$a_n[e2], cc$a_t[e2], a_ref, "full")
  e3 <- ev & cc$src == "sr"
  data.frame(archetype = A, n_ev = sum(ev), frac_sr_ev = mean(cc$src[ev] == "sr"),
             rho_all = f_all[["rho"]], beta_all = f_all[["beta"]], sigma_all = f_all[["sigma"]],
             rho_sd = f_sd[["rho"]], beta_sd = f_sd[["beta"]], sigma_sd = f_sd[["sigma"]],
             phi_all = unname(quantile(cc$a_n[ev] / cc$a_t[ev], 0.99)), phi_sd = unname(quantile(cc$a_n[e2] / cc$a_t[e2], 0.99)),
             med_area_sd = median(cc$a_n[e2]), med_area_sr = median(cc$a_n[e3]), med_at_sd = median(cc$a_t[e2]), med_at_sr = median(cc$a_t[e3]),
             med_trim_sr = median(cc$trim[e3]), q90_trim_sr = unname(quantile(cc$trim[e3], 0.9)))
}))
write.csv(out, "results/R3b/R3b_posthoc_area.csv", row.names = FALSE); print(format(out, digits = 4), row.names = FALSE)
