# tools/S1.6_explore_bound.R - S1.6, esplorazione E4 (BL-161), prima della pre-registrazione: per ogni regione dei 30 strati salvati del braccio Q di
# S1.5e (seme 1), limite euristico di oggi wmax = 1.1 x massimo di w sulla griglia contro il limite esatto max_t rho_t dei tipi della regione
# (w = 1 / sum_t p_t / rho_t <= max_t rho_t; nei preset la densita' non varia, dens_sigma = 0, e i vuoti valgono 0 o 1).
# Il rapporto rho_max / wmax e' il fattore per cui cambia il numero di proposte per cellula accettata (< 1: meno proposte).
# Uscita: results/S1.6/S1.6_explore_bound.csv (per regione), stampa il rapporto pesato per area per tessuto.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({library(sf)})
fs <- list.files("/mnt/micron/geo_spatialtrans/S1.5e/layers", pattern = "_Q\\.rds$", full.names = TRUE)
out <- do.call(rbind, lapply(fs, function(f) {
  L <- readRDS(f); ct <- L$cell_types_resolved; ci <- L$metadata$seed_info$compartments; comp <- L$region_composition; rdf <- L$region_df
  do.call(rbind, lapply(ci$regions, function(r) {
    cl <- rdf$cluster_id[rdf$region_id == r$region_id]; ty <- comp$cell_type[comp$cluster_id == cl]
    data.frame(layer = basename(f), archetype = substr(basename(f), 1, 2), region_id = r$region_id, area_um2 = rdf$area_um2[rdf$region_id == r$region_id],
               wmax_old = r$wmax, w_grid_max = r$w_grid_max, w_mean = r$w_mean, rho_max = max(ct$density[match(ty, ct$cell_type)]))
  }))
}))
out$ratio <- out$rho_max / out$wmax_old
dir.create("results/S1.6", showWarnings = FALSE); write.csv(out, "results/S1.6/S1.6_explore_bound.csv", row.names = FALSE)
cat(sprintf("E4: regioni %d; w_grid_max <= rho_max in tutte: %s\n", nrow(out), all(out$w_grid_max <= out$rho_max * (1 + 1e-12))))
for (A in sort(unique(out$archetype))) { o <- out[out$archetype == A, ]
  cat(sprintf("E4 %s: rapporto rho_max / wmax pesato per area %.3f (min %.3f, max %.3f, %d regioni)\n", A, sum(o$ratio * o$area_um2) / sum(o$area_um2), min(o$ratio), max(o$ratio), nrow(o))) }
