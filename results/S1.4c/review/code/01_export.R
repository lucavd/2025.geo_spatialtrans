# review/code/01_export.R — esporta in parquet/csv (scratch) gli oggetti RDS che servono ai ricalcoli indipendenti in Python.
# Sola lettura dei grezzi; nessuna funzione di tools/S1.4c_analysis.R o tools/S1.4c_decision.R.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(sf); library(geos) })
source("R/04b1_extract_regions.R"); source("R/04b3_tessellate_voronoi.R"); source("tools/R3_voronoi_metrics.R")
SC <- "/mnt/micron/geo_spatialtrans/S1.4c_review/code"; dir.create(file.path(SC, "exp"), showWarnings = FALSE, recursive = TRUE)
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"; S14C <- "/mnt/micron/geo_spatialtrans/S1.4c"
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
prim <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
pr_rows <- list(); reg_rows <- list(); mi <- list()
for (i in seq_len(nrow(rois))) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; m <- prim[[A]]
  r <- readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, roi, m)))
  p <- readRDS(file.path(S14C, "prep", sprintf("%s_%s_%s.rds", A, roi, m)))
  cc <- r$cells[, c("idx", "x", "y", "a_t", "interior", "matched", "a_n", "a_own", "d_gen")]
  cc$ev_prep <- p$ev; cc$shape_label <- p$shape_label
  write_parquet(cc, file.path(SC, "exp", sprintf("cells_%s_%s.parquet", A, roi)))
  pr_rows[[i]] <- data.frame(archetype = A, roi = roi, n = p$n, b7_real = p$b7_real, b8_real = p$b8_real)
  w <- r3_roi_window(A, roi, rois)
  reg_rows[[i]] <- data.frame(archetype = A, roi = roi, area_um2 = w$area_um2, wkt = st_as_text(w$poly[[1]]))
  if (!is.null(r$match_info)) mi[[length(mi) + 1]] <- data.frame(archetype = A, roi = roi, n_sr = r$n, t(r$match_info))
}
write.csv(do.call(rbind, pr_rows), file.path(SC, "exp", "prep.csv"), row.names = FALSE)
write.csv(do.call(rbind, reg_rows), file.path(SC, "exp", "regions.csv"), row.names = FALSE)
write.csv(do.call(rbind, mi), file.path(SC, "exp", "match_info.csv"), row.names = FALSE)
lo <- list(); yd <- list(); gn <- list(); ck <- list(); b4 <- list()
for (A in names(prim)) for (k in 1:5) {
  o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, prim[[A]], k)))
  lo[[length(lo) + 1]] <- data.frame(archetype = A, fold = k, roi_test = o$roi_test, seed = o$seed, t(o$fits$M2))
  yd[[length(yd) + 1]] <- data.frame(archetype = A, fold = k, roi_test = o$roi_test, pos = seq_len(nrow(o$yard)), o$yard)
  g <- readRDS(file.path(S14C, "gen", sprintf("%s_%s_fold%d.rds", A, prim[[A]], k)))
  gn[[length(gn) + 1]] <- data.frame(archetype = A, fold = k, roi_test = g$roi_test, seed = g$seed, phi = g$phi, q99 = g$q99, n = g$n, n_ev = g$n_ev)
  z <- g$checks; z$archetype <- A; z$fold <- k; ck[[length(ck) + 1]] <- z
  z <- g$b4; z$archetype <- A; z$fold <- k; b4[[length(b4) + 1]] <- z
}
bind <- function(z) { nm <- unique(unlist(lapply(z, names))); do.call(rbind, lapply(z, function(d) { for (m in setdiff(nm, names(d))) d[[m]] <- NA; d[, nm] })) }
write.csv(do.call(rbind, lo), file.path(SC, "exp", "loo.csv"), row.names = FALSE)
write.csv(do.call(rbind, yd), file.path(SC, "exp", "yard14.csv"), row.names = FALSE)
write.csv(do.call(rbind, gn), file.path(SC, "exp", "gen.csv"), row.names = FALSE)
write.csv(bind(ck), file.path(SC, "exp", "gen_checks.csv"), row.names = FALSE)
write.csv(do.call(rbind, b4), file.path(SC, "exp", "gen_b4.csv"), row.names = FALSE)
mo <- do.call(rbind, lapply(list.files(file.path(S14C, "mut"), pattern = "_MU5\\.rds$", full.names = TRUE), readRDS))
write.csv(mo, file.path(SC, "exp", "mut_MU5_obs.csv"), row.names = FALSE)
mg <- lapply(list.files(file.path(S14C, "mut"), pattern = "_MU[12]\\.rds$", full.names = TRUE), function(f) { g <- readRDS(f); z <- g$checks; z$mutant <- g$mutant; z$archetype <- g$archetype; z })
write.csv(bind(mg), file.path(SC, "exp", "mut_MU12_checks.csv"), row.names = FALSE)
cat("EXPORT_DONE\n")
