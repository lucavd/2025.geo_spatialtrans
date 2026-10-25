# tools/S1.5d_export_r3b_cells.R - S1.5d: celle su cui e' stimato il modello di taglia dei preset del linfonodo e del cuore (R3b, variante
# «sdarea», DEC-006: celle interne, nuclei di forma StarDist della maschera CO, a_n > 0; stessa selezione di s14_metrics() in tools/R3b_refs.R),
# piu' l'area di riferimento a_ref di quella stima. Uscita: /mnt/micron/geo_spatialtrans/S1.5d/r3b_sdarea_<A>.parquet (roi, a_n, a_t, a_ref).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow) })
M <- "/mnt/micron/geo_spatialtrans"; O <- file.path(M, "R3b"); OUT <- file.path(M, "S1.5d"); dir.create(OUT, showWarnings = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")
for (A in c("A4", "A6")) {
  ids <- rois$roi_id[rois$archetype == A]
  rs <- lapply(ids, function(r) { z <- readRDS(file.path(O, "S1.4", "real", sprintf("%s_%s_consensus.rds", A, r)))
    g <- as.data.frame(read_parquet(file.path(O, "R3", "real", sprintf("%s_%s_consensus_gen.parquet", A, r))))
    mp <- as.data.frame(read_parquet(file.path(O, "masks", sprintf("%s_%s_consensus_map.parquet", A, r)))); z$cells$matched <- mp$src[match(g$label, mp$label)] == "sd"
    z$cells$roi <- r; z })
  a_ref <- sum(vapply(rs, function(r) sum(r$cells$a_t), 0)) / sum(vapply(rs, function(r) r$n, 0))
  cc <- do.call(rbind, lapply(rs, function(r) { k <- r$cells$interior & r$cells$matched & is.finite(r$cells$a_n) & r$cells$a_n > 0
    data.frame(roi = r$cells$roi[k], a_n = r$cells$a_n[k], a_t = r$cells$a_t[k]) }))
  cc$a_ref <- a_ref; write_parquet(cc, file.path(OUT, sprintf("r3b_sdarea_%s.parquet", A))); cat(A, nrow(cc), "celle, a_ref", a_ref, "densita'", 1e6 / a_ref, "\n")
}
