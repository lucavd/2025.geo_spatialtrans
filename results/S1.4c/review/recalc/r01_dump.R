
# results/S1.4c/review/recalc/r01_dump.R — revisione avversariale S1.4c (recalc): SOLO I/O. Scarica i grezzi in parquet/csv
# per il ricalcolo indipendente in Python (r02_recalc.py). Nessuna funzione di analisi/decisione del banco e' usata.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(arrow) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R")   # solo r3_roi_window() (finestra = input)
SCR <- "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/dump"; dir.create(SCR, recursive = TRUE, showWarnings = FALSE)
B <- "/mnt/micron/geo_spatialtrans"
ARCH <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
# obs (900)
fs <- list.files(file.path(B, "S1.4c/obs"), full.names = TRUE, pattern = "\\.rds$")
obs <- do.call(rbind, lapply(fs, readRDS)); write.csv(obs, file.path(SCR, "obs_all.csv"), row.names = FALSE)
cat("obs", nrow(obs), "\n")
# loo S1.4: yard + roi_test + seed ; gen S1.4c: b4, phi ; S1.4b loo: V-c rows
yd <- list(); b4 <- list(); ph <- list(); vc <- list()
for (A in names(ARCH)) for (k in 1:5) {
  o <- readRDS(file.path(B, "S1.4/loo", sprintf("%s_%s_fold%d.rds", A, ARCH[[A]], k)))
  yd[[length(yd) + 1]] <- data.frame(archetype = A, fold = k, roi_test = o$roi_test, seed = o$seed, o$yard)
  g <- readRDS(file.path(B, "S1.4c/gen", sprintf("%s_%s_fold%d.rds", A, ARCH[[A]], k)))
  b4[[length(b4) + 1]] <- data.frame(archetype = A, fold = k, roi_test = g$roi_test, g$b4)
  ph[[length(ph) + 1]] <- data.frame(archetype = A, fold = k, phi = g$phi, q99 = g$q99, n = g$n, n_ev = g$n_ev)
  ob <- readRDS(file.path(B, "S1.4b/loo", sprintf("%s_%s_fold%d.rds", A, ARCH[[A]], k)))
  m <- ob$m[ob$m$variant == "V-c", ]; vc[[length(vc) + 1]] <- data.frame(archetype = A, fold = k, m)
}
write.csv(do.call(rbind, yd), file.path(SCR, "loo_yard.csv"), row.names = FALSE)
write.csv(do.call(rbind, b4), file.path(SCR, "gen_b4.csv"), row.names = FALSE)
write.csv(do.call(rbind, ph), file.path(SCR, "gen_phi.csv"), row.names = FALSE)
vcd <- do.call(rbind, vc); write.csv(vcd, file.path(SCR, "s14b_Vc.csv"), row.names = FALSE)
# per ROI: celle S1.4, prep S1.4c, territori, finestra
pp <- list()
for (i in seq_len(nrow(rois))) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; met <- ARCH[[A]]
  r <- readRDS(file.path(B, "S1.4/real", sprintf("%s_%s_%s.rds", A, roi, met))); cc <- r$cells
  p <- readRDS(file.path(B, "S1.4c/prep", sprintf("%s_%s_%s.rds", A, roi, met)))
  stopifnot(length(p$ev) == nrow(cc))
  keep <- intersect(c("idx", "x", "y", "a_t", "interior", "matched", "a_n", "a_own", "d_gen"), names(cc))
  d <- cc[, keep]; d$ev_prep <- p$ev; d$shape_label <- p$shape_label
  write_parquet(d, file.path(SCR, sprintf("cells_%s_%s.parquet", A, roi)))
  if (met == "spaceranger") {
    tw <- st_as_binary(st_geometry(r$terr))
    write_parquet(data.frame(i = seq_along(tw), wkb = I(unclass(tw))), file.path(SCR, sprintf("terr_%s_%s.parquet", A, roi)))
  }
  w <- r3_roi_window(A, roi, rois)
  write_parquet(data.frame(area_um2 = w$area_um2, wkb = I(unclass(st_as_binary(w$poly)))), file.path(SCR, sprintf("win_%s_%s.parquet", A, roi)))
  pp[[i]] <- data.frame(archetype = A, roi = roi, n = p$n, n_ev = sum(p$ev), b7_real = p$b7_real, b8_real = p$b8_real,
                        win_area = w$area_um2, win_poly_area = as.numeric(st_area(w$poly)), side_um = rois$side_um[i])
}
write.csv(do.call(rbind, pp), file.path(SCR, "prep_scalars.csv"), row.names = FALSE)
cat("names cells:", names(cc), "\n"); cat("names S1.4b m:", names(vcd), "\n"); cat("DONE\n")
