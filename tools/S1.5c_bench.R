# tools/S1.5c_bench.R - S1.5c: banco del pannello con la semina a compartimenti (G-3, P-1...P-10, CP-1, descrittive). Definizioni: docs/steps/S1.5c.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5c_bench.R <stadio>
#   md5before | gen | check | md5after   (stadi rilanciabili: saltano le unita' fatte)
# Bracci: U (semina uniforme di S1.5b, ricentrato), C-a, C-b (compartimenti, distanza del tessuto; gruppo / ricentrato), N-a, N-b
# (compartimenti, distanza per nucleo), G (compartimenti, etichette cieche, ricentrato: controprova). Parametri dei compartimenti della piega
# lascia-fuori del ROI simulato (tools/S1.5c_bench_lib.R). Uscite stabili in /mnt/micron/geo_spatialtrans/S1.5c (un file per unita').
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5c_bench_lib.R")
options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")   # S1.6: predefiniti di S1.5c-S1.5e, per riprodurre i numeri di questo step
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); STAGE <- if (length(args)) args[1] else "none"
M <- "/mnt/micron/geo_spatialtrans"; OUT <- file.path(M, "S1.5c"); RES <- "results/S1.5c"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
wk_file <- file.path(OUT, "WORKERS"); if (!file.exists(wk_file)) writeLines("3", wk_file)
NW <- function() min(28L, as.integer(readLines(wk_file)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
ARMS <- data.frame(arm = c("U", "C-a", "C-b", "N-a", "N-b", "G"), spatial = c("uniform", rep("compartments", 5)),
                   mode = c("tissue", "tissue", "tissue", "nucleus", "nucleus", "tissue"),
                   size_arm = c("recentred", "group", "recentred", "group", "recentred", "recentred"),
                   labels = c("field", "field", "field", "field", "field", "blind"), stringsAsFactors = FALSE)
SEEDS <- 1:10; BASE <- 20261021
seed_of <- function(i, s) BASE + 1000 * i + s                      # come S1.5 e S1.5b
THR <- 59.1                                                         # banda del cervello (S1.5b, small_threshold_um2)
WIN <- read.csv("results/R2b/R2b_windows.csv"); WIN_ROI <- unique(paste(WIN$archetype, WIN$roi_id))
clust_of <- function(w) { idx <- which(w$m, arr.ind = TRUE); data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L) }
layer_md5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  on.exit(unlink(f)); unname(tools::md5sum(f)) }
run_pc <- function(i, s, arm, w) {
  A <- rois$archetype[i]; a <- ARMS[ARMS$arm == arm, ]
  oc <- if (a$spatial == "compartments") compart_for(A, rois$roi_id[i]) else NULL
  suppressWarnings(place_cells(clust_of(w), PRESET[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                               small_components = "exclude", random_seed = seed_of(i, s), types = "multi", size_arm = a$size_arm,
                               spatial = a$spatial, min_dist_mode = a$mode, labels = a$labels, override_compartments = oc,
                               keep_tiles = TRUE, strict = FALSE, verbose = FALSE))
}
ecc_of <- function(g) vapply(seq_along(g), function(k) { xy <- sf::st_coordinates(g[[k]])[, 1:2]
  m <- r3_poly_moments(xy[, 1], xy[, 2]); r3_shape(m[["mxx"]], m[["myy"]], m[["mxy"]])$ecc }, 0)
metrics_of <- function(L, w, A) {
  cd <- L$cell_df; td <- L$territory_df
  in_frame <- geos::geos_within(geos::as_geos_geometry(L$tiles), geos::as_geos_geometry(w$frame))
  it <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
  a <- cd$territory_area[it]; nc <- cd$nucleus_area / cd$territory_area; vt <- L$metadata$validation$table
  nofrag <- it & td$n_pieces_gained == 0; a2 <- cd$territory_area[nofrag]
  ty <- as.character(cd$cell_type); grp <- L$cell_types_resolved; gi <- grp$size_group[match(ty, grp$cell_type)]
  ci <- L$metadata$seed_info$compartments
  u <- data.frame(n = nrow(cd), n_interior = sum(it), median_area_nuc = median(cd$nucleus_area[it]), median_nc = median(nc[it]),
                  frac_nc095 = mean(nc[it] >= 0.95), frac_trunc = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target),
                  cv_terr = sd(a) / mean(a), cv_terr_nofrag = sd(a2) / mean(a2), median_ecc_T = median(ecc_of(L$cell_territories[which(it)])),
                  n_failed = sum(L$metadata$seed_info$n_failed), n_error = L$metadata$validation$n_error, n_warn = L$metadata$validation$n_warn,
                  statuses = paste(vt$check, vt$status, sep = ":", collapse = " "), elapsed_s = L$metadata$elapsed_s, md5 = layer_md5(L),
                  n_exceed_wmax = if (is.null(ci)) NA_integer_ else ci$n_exceed_wmax, frac_band = NA_real_, ratio_neuron_glia = NA_real_)
  if (A == "A5") {
    an <- cd$nucleus_area[it]; u$frac_band <- mean(an > THR / sqrt(2) & an < THR * sqrt(2))
    gN <- grp$size_group[grp$cell_type == "neuroni"]; gg <- gi[it]; u$ratio_neuron_glia <- median(an[gg == gN]) / median(an[gg != gN])
  }
  tt <- do.call(rbind, lapply(sort(unique(ty)), function(t) { k <- it & ty == t
    data.frame(cell_type = t, n_type = sum(ty == t), n_interior = sum(k), median_area_nuc = if (any(k)) median(cd$nucleus_area[k]) else NA_real_,
               median_nc = if (any(k)) median(nc[k]) else NA_real_, median_terr = if (any(k)) median(cd$territory_area[k]) else NA_real_) }))
  cells <- data.frame(x = cd$x, y = cd$y, cell_type = ty, interior = it, a_n = cd$nucleus_area, a_t = cd$territory_area)
  list(unit = u, type = tt, cells = cells)
}
par_run <- function(jobs, fun, label) {
  t0 <- Sys.time(); r <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch({ fun(jobs[k, ]); "ok" }, error = function(e) paste("ERR", conditionMessage(e))),
                                  mc.cores = NW(), mc.preschedule = FALSE)
  r <- unlist(r); cat(sprintf("[%s] %d unita', %d errori, %.1f min\n", label, length(r), sum(grepl("^ERR", r)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (any(grepl("^ERR", r))) print(head(cbind(jobs[grepl("^ERR", r), ], err = r[grepl("^ERR", r)]), 10))
}
md5_inputs <- function() {
  f <- c(list.files("results/R2/tissue_masks", full.names = TRUE), list.files(file.path(M, "R6/labels"), pattern = "exp_umi20_seed20261022_nuclei_ph1.parquet$", full.names = TRUE),
         list.files(file.path(M, "R6/terr"), full.names = TRUE), "results/R2b/R2b_windows.csv", "results/R2c/R2c_dmax.csv",
         "results/R2b/annotations/R2b_annotations_Luca.json", "results/S1.5/S1.5_real_refs_roi.csv", "results/S1.5b/S1.5b_real_refs_roi_type.csv",
         "results/S1.5b/S1.5b_panel_units.csv", list.files(RES, pattern = "^S1.5c_(compartments|real_|type_nuc)", full.names = TRUE), "docs/BIO_REFERENCES.md")
  data.frame(file = f, md5 = unname(tools::md5sum(f)))
}
if (STAGE == "md5before") { write.csv(md5_inputs(), file.path(OUT, "md5_before.csv"), row.names = FALSE); cat("md5 prima:", nrow(md5_inputs()), "file\n") }
if (STAGE == "md5after") {
  b <- read.csv(file.path(OUT, "md5_before.csv")); a <- md5_inputs(); write.csv(a, file.path(OUT, "md5_after.csv"), row.names = FALSE)
  same <- nrow(a) == nrow(b) && identical(a$file, b$file) && identical(a$md5, b$md5)
  cat(sprintf("%s md5 dei dati letti invariati: %d/%d\n", if (same) "PASS" else "FAIL", sum(a$md5 == b$md5[match(a$file, b$file)], na.rm = TRUE), nrow(b)))
}
if (STAGE == "gen") {
  for (arm in ARMS$arm) for (d in c("gen", "cells", "nuc")) dir.create(file.path(OUT, d, arm), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUT, "layers"), showWarnings = FALSE)
  jobs <- expand.grid(i = rois$roi_idx, s = SEEDS, arm = ARMS$arm, stringsAsFactors = FALSE)
  jobs <- jobs[order(jobs$s, -(rois$archetype[jobs$i] %in% c("A4", "A5")), jobs$arm), ]
  par_run(jobs, function(j) {
    i <- j$i; s <- j$s; arm <- j$arm; A <- rois$archetype[i]; roi <- rois$roi_id[i]
    f <- file.path(OUT, "gen", arm, sprintf("%s_%s_s%02d.rds", A, roi, s)); if (file.exists(f)) return(invisible())
    w <- r3_roi_window(A, roi, rois); L <- run_pc(i, s, arm, w); m <- metrics_of(L, w, A)
    key <- data.frame(archetype = A, roi_id = roi, roi_idx = i, seed = s, arm = arm, random_seed = seed_of(i, s))
    write_parquet(m$cells, file.path(OUT, "cells", arm, sprintf("%s_%s_s%02d.parquet", A, roi, s)))
    if (paste(A, roi) %in% WIN_ROI) {
      g <- geos::as_geos_geometry(L$cell_nuclei)
      write_parquet(data.frame(idx = seq_along(g), wkb = I(geos::geos_write_wkb(g))), file.path(OUT, "nuc", arm, sprintf("%s_%s_s%02d.parquet", A, roi, s)))
    }
    if (s == 1) saveRDS(L, file.path(OUT, "layers", sprintf("%s_%s_%s.rds", A, roi, arm)))
    saveRDS(list(unit = cbind(key, m$unit), type = cbind(key, m$type)), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  }, "gen")
  z <- lapply(list.files(file.path(OUT, "gen"), recursive = TRUE, pattern = "rds$", full.names = TRUE), readRDS)
  write.csv(do.call(rbind, lapply(z, `[[`, "unit")), file.path(RES, "S1.5c_panel_units.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(z, `[[`, "type")), file.path(RES, "S1.5c_panel_types.csv"), row.names = FALSE)
  cat("unita':", length(z), "\n")
}
if (STAGE == "check") {
  U <- read.csv(file.path(RES, "S1.5c_panel_units.csv"), stringsAsFactors = FALSE)
  e1 <- sum(grepl("C2:ERROR|C3:ERROR|C5:ERROR|C9:ERROR", U$statuses)); nf <- sum(U$n_failed)
  S <- read.csv("results/S1.5b/S1.5b_panel_units.csv", stringsAsFactors = FALSE); S <- S[S$arm == "recentred", ]
  uu <- U[U$arm == "U", ]; k <- match(paste(uu$archetype, uu$roi_id, uu$seed), paste(S$archetype, S$roi_id, S$seed))
  compat <- sum(uu$md5 == S$md5[k], na.rm = TRUE)
  pick <- do.call(rbind, lapply(split(rois, rois$archetype), function(d) d[1, ]))
  arms6 <- c("C-a", "C-b", "N-a", "N-b", "G", "U")
  rep_ok <- vapply(seq_len(nrow(pick)), function(r) { i <- pick$roi_idx[r]; arm <- arms6[r]
    w <- r3_roi_window(pick$archetype[r], pick$roi_id[r], rois); L <- run_pc(i, 3, arm, w)
    identical(layer_md5(L), U$md5[U$roi_idx == i & U$seed == 3 & U$arm == arm]) }, TRUE)
  out <- data.frame(check = c("units", "no_error_C2_C3_C5_C9", "n_failed_total", "U_equals_S1.5b_recentred_md5", "rerun_identical"),
                    value = c(nrow(U), nrow(U) - e1, nf, compat, sum(rep_ok)), expected = c(1800, 1800, 0, 300, 6))
  out$pass <- out$value == out$expected; write.csv(out, file.path(RES, "S1.5c_G3_check.csv"), row.names = FALSE); print(out)
  cat(sprintf("%s G-3 (md5 dei dati: stadio md5after)\n", if (all(out$pass)) "PASS" else "FAIL"))
}
