# tools/R3b_refs.R — R3b: riferimenti di linfonodo (A4) e cuore (A6) ricalcolati sulla maschera CO, con lo stesso codice che li ha
# prodotti (tools/R3_analysis.R: celle interne del Voronoi R3, mediane sul pool e CV locale = mediana dei 5 ROI; tools/S1.4_analysis.R:
# M2 sui 5 ROI, C8 = [q05, q95] di eq_r delle celle interne di S1.4; q99 del rapporto nucleo/territorio come B-121...B-126).
# Controllo: lo stesso codice sulle uscite originali (Space Ranger, R3 e S1.4) deve ridare il valore delle righe di BIO_REFERENCES.
# Variante «sdarea»: modello di taglia, phi e N/C sui soli nuclei di forma StarDist della CO (decisione di Luca, DEC-006).
# Uscita: results/R3b/R3b_refs.csv.  Uso: Rscript --vanilla tools/R3b_refs.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(sf); library(geos) })
source("tools/S1.4d_mutants.R"); E <- dcg_env_d()
M <- "/mnt/micron/geo_spatialtrans"; O <- file.path(M, "R3b")
rois <- read.csv("results/R2/R2_rois_checked.csv")
bio <- readLines("docs/BIO_REFERENCES.md"); bio <- bio[grepl("^\\| B-[0-9]{3} \\|", bio)]
bval <- function(id) { r <- bio[startsWith(bio, paste0("| ", id, " |"))]; trimws(strsplit(r, "\\|")[[1]][5]) }
OLD <- list(A4 = c(area = "B-073", eq_r = "B-074", cv_loc = "B-075", ecc_T = "B-076", nc = "B-077", ratio = "B-078", rho = "B-100", beta = "B-101", sigma = "B-102", c8 = "B-112"),
            A6 = c(area = "B-085", eq_r = "B-086", cv_loc = "B-087", ecc_T = "B-088", nc = "B-089", ratio = "B-090", rho = "B-106", beta = "B-107", sigma = "B-108", c8 = "B-114"))
r3_metrics <- function(A, method, base, sd_only = FALSE) {
  ids <- rois$roi_id[rois$archetype == A]
  d <- do.call(rbind, lapply(ids, function(r) { tv <- as.data.frame(read_parquet(file.path(base, "real", sprintf("%s_%s_%s_cells.parquet", A, r, method))))
    if (sd_only) { mp <- as.data.frame(read_parquet(file.path(O, "masks", sprintf("%s_%s_consensus_map.parquet", A, r)))); tv <- tv[mp$src[match(tv$label, mp$label)] == "sd", ] }
    tv <- tv[tv$interior, ]; data.frame(roi = r, area = tv$area, ecc_T = tv$ecc_T, area_nuc = tv$area_nuc) }))
  cvl <- vapply(ids, function(r) readRDS(file.path(base, "real", sprintf("%s_%s_%s_roi.rds", A, r, method)))$cv_loc, 0)
  rm <- tapply(d$area, d$roi, median); nc <- d$area_nuc / d$area
  c(n_interior = nrow(d), area = median(d$area), area_q25 = unname(quantile(d$area, .25)), area_q75 = unname(quantile(d$area, .75)),
    area_roi_min = min(rm), area_roi_max = max(rm), eq_r = median(sqrt(d$area / pi)), cv_loc = median(cvl), cv_loc_min = min(cvl), cv_loc_max = max(cvl),
    ecc_T = median(d$ecc_T), nc = median(nc), nc_q25 = unname(quantile(nc, .25)), nc_q75 = unname(quantile(nc, .75)), ratio = median(sqrt(nc)))
}
s14_metrics <- function(A, method, base, sd_only = FALSE) {
  ids <- rois$roi_id[rois$archetype == A]
  rs <- lapply(ids, function(r) { z <- readRDS(file.path(base, "real", sprintf("%s_%s_%s.rds", A, r, method)))
    if (sd_only) { g <- as.data.frame(read_parquet(file.path(O, "R3", "real", sprintf("%s_%s_consensus_gen.parquet", A, r))))
      mp <- as.data.frame(read_parquet(file.path(O, "masks", sprintf("%s_%s_consensus_map.parquet", A, r)))); z$cells$matched <- mp$src[match(g$label, mp$label)] == "sd" }
    z })
  a_ref <- sum(vapply(rs, function(r) sum(r$cells$a_t), 0)) / sum(vapply(rs, function(r) r$n, 0))
  cc <- do.call(rbind, lapply(rs, function(r) r$cells[r$cells$interior & r$cells$matched & is.finite(r$cells$a_n) & r$cells$a_n > 0, ]))
  f <- E$fit_nucleus_size_model(cc$a_n, cc$a_t, a_ref, "full")
  eq <- unlist(lapply(rs, function(r) sqrt(r$cells$a_t[r$cells$interior] / pi))); q <- quantile(eq, c(.05, .5, .95))
  c(n_ev = nrow(cc), density_mm2 = 1e6 / a_ref, rho = f[["rho"]], beta = f[["beta"]], sigma = f[["sigma"]],
    phi = unname(quantile(cc$a_n / cc$a_t, 0.99, type = 7)), c8_q05 = q[[1]], c8_q50 = q[[2]], c8_q95 = q[[3]])
}
out <- list()
for (A in c("A4", "A6")) {
  old <- c(r3_metrics(A, "spaceranger", file.path(M, "R3")), s14_metrics(A, "spaceranger", file.path(M, "S1.4")))
  new <- c(r3_metrics(A, "consensus", file.path(O, "R3")), s14_metrics(A, "consensus", file.path(O, "S1.4")))
  sda <- c(r3_metrics(A, "consensus", file.path(O, "R3"), sd_only = TRUE)[c("nc", "nc_q25", "nc_q75", "ratio")], s14_metrics(A, "consensus", file.path(O, "S1.4"), sd_only = TRUE))
  for (k in names(new)) out[[length(out) + 1]] <- data.frame(archetype = A, metric = k, bio_id = if (k %in% names(OLD[[A]])) OLD[[A]][[k]] else if (k == "c8_q05") OLD[[A]][["c8"]] else "",
    bio_value = if (k %in% names(OLD[[A]])) bval(OLD[[A]][[k]]) else if (k == "c8_q05") bval(OLD[[A]][["c8"]]) else "",
    old_recomputed = old[[k]], consensus = new[[k]], consensus_sdarea = if (k %in% names(sda)) sda[[k]] else NA_real_,
    rel_change = new[[k]] / old[[k]] - 1)
}
out <- do.call(rbind, out); write.csv(out, "results/R3b/R3b_refs.csv", row.names = FALSE)
print(format(out, digits = 4), row.names = FALSE)
