# tools/S1.6_bench.R - S1.6: banco del consolidamento (G-3, P-1, CP-2, D-1...D-3). Definizioni: docs/steps/S1.6.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.6_bench.R <stadio> [bracci]
#   md5before | gen <bracci> | check | md5after   (stadi rilanciabili: saltano le unita' fatte)
# Bracci: U6 (semina dei preset, deve essere U di S1.5e bit per bit), Q6 (compartimenti con i vuoti RV della piega come Q di S1.5e, con i
# predefiniti nuovi: senza quote e con il limite esatto del diradamento; nessun argomento e nessuna opzione di sessione). Per il controllo
# di riproducibilita' (check) anche RVl e Ql: RV e Q di S1.5e con le opzioni esplicite di S1.5e (quote, limite sulla griglia).
# Funzioni di misura copiate da tools/S1.5e_bench.R (identiche), salvo la frazione dei vuoti pesata per area (BL-169, RC-6) e due colonne
# nuove (regioni con il limite sulla griglia, quote attive). Uscite stabili in /mnt/micron/geo_spatialtrans/S1.6. Semi come S1.5c-S1.5e.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5c_bench_lib.R"); source("tools/S1.5d_bench_lib.R"); source("tools/S1.5e_bench_lib.R")
stopifnot(is.null(getOption("geo.spatialtrans.quotas")), is.null(getOption("geo.spatialtrans.thinning_bound")))   # predefiniti nuovi del pacchetto
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); STAGE <- if (length(args)) args[1] else "none"
M <- "/mnt/micron/geo_spatialtrans"; OUT <- file.path(M, "S1.6"); RES <- "results/S1.6"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
wk_file <- file.path(OUT, "WORKERS"); if (!file.exists(wk_file)) writeLines("20", wk_file)
NW <- function() min(28L, as.integer(readLines(wk_file)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
ARMS <- data.frame(arm = c("U6", "Q6", "RVl", "Ql"), spatial = c("uniform", rep("compartments", 3)), het = c(NA, "RV", "RV", "RV"),
                   size_arm = c("recentred", rep("group", 3)), quotas = c(NA, NA, TRUE, FALSE), bound = c(NA, NA, "grid", "grid"),
                   ref = c("U", "Q", "RV", "Q"), stringsAsFactors = FALSE)   # quotas/bound NA: predefiniti del pacchetto (S1.6)
het_of <- function(A, roi, fam) het_for(A, roi, "RV")                  # vuoti RV di S1.5d per piega, come RV e Q di S1.5e
SEEDS <- 1:10; BASE <- 20261021
seed_of <- function(i, s) BASE + 1000 * i + s
THR <- 59.1
WIN <- read.csv("results/R2b/R2b_windows.csv"); WIN_ROI <- unique(paste(WIN$archetype, WIN$roi_id))
clust_of <- function(w) { idx <- which(w$m, arr.ind = TRUE); data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L) }
layer_md5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  on.exit(unlink(f)); unname(tools::md5sum(f)) }
run_pc <- function(i, s, arm, w) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARMS[ARMS$arm == arm, ]
  if (a$spatial == "uniform")
    return(suppressWarnings(place_cells(clust_of(w), PRESET[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                               small_components = "exclude", random_seed = seed_of(i, s), types = "multi", size_arm = a$size_arm,
                               spatial = "uniform", size_resid = "gaussian", keep_tiles = TRUE, strict = FALSE, verbose = FALSE)))
  extra <- list(); if (!is.na(a$quotas)) extra$quotas <- a$quotas; if (!is.na(a$bound)) extra$thinning_bound <- a$bound
  suppressWarnings(do.call(place_cells, c(list(clust_of(w), PRESET[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                               small_components = "exclude", random_seed = seed_of(i, s), types = "multi", size_arm = a$size_arm,
                               spatial = "compartments", min_dist_mode = "nucleus", override_compartments = compart_for_d(A, roi),
                               override_heterogeneity = het_of(A, roi, a$het), size_resid = "gaussian", resid_mode = "decile",
                               keep_tiles = TRUE, strict = FALSE, verbose = FALSE), extra)))
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
  ci <- L$metadata$seed_info$compartments; rr <- if (is.null(ci)) NULL else ci$regions
  nrem <- if (is.null(rr)) 0L else sum(vapply(rr, function(z) if (is.null(z$n_removed)) 0 else z$n_removed, 0))
  rarea <- if (is.null(rr)) NULL else L$region_df$area_um2[match(vapply(rr, function(z) z$region_id, 0), L$region_df$region_id)]
  vgf <- if (is.null(rr)) 0 else stats::weighted.mean(vapply(rr, function(z) z$void_grid_frac, 0), rarea)   # S1.6: pesata per area (BL-169, RC-6)
  n_grid_bound <- if (is.null(rr)) NA_integer_ else sum(vapply(rr, function(z) !identical(z$bound, "exact"), TRUE))
  capH <- if (is.null(rr) || is.null(rr[[1]]$cap)) NA_real_ else rr[[1]]$cap$H
  u <- data.frame(n = nrow(cd), n_interior = sum(it), median_area_nuc = median(cd$nucleus_area[it]), median_nc = median(nc[it]),
                  frac_nc095 = mean(nc[it] >= 0.95), frac_trunc = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target),
                  cv_terr = sd(a) / mean(a), cv_terr_nofrag = sd(a2) / mean(a2), median_ecc_T = median(ecc_of(L$cell_territories[which(it)])),
                  n_failed = sum(L$metadata$seed_info$n_failed), n_error = L$metadata$validation$n_error, n_warn = L$metadata$validation$n_warn,
                  statuses = paste(vt$check, vt$status, sep = ":", collapse = " "), elapsed_s = L$metadata$elapsed_s, md5 = layer_md5(L),
                  n_exceed_wmax = if (is.null(ci)) NA_integer_ else ci$n_exceed_wmax, n_attempts = L$metadata$seed_info$n_attempts,
                  n_quota_reject = if (is.null(ci)) NA_integer_ else ci$n_quota_reject, n_thin_accept = if (is.null(ci)) NA_integer_ else ci$n_thin_accept,
                  n_void_disabled = if (is.null(ci)) NA_integer_ else ci$n_void_disabled, n_removed_void = nrem, void_grid_frac = vgf, cap_H = capH,
                  frac_band = NA_real_, ratio_neuron_glia = NA_real_, n_regions_grid_bound = n_grid_bound,
                  quotas_on = if (is.null(ci)) NA else isTRUE(ci$quotas))
  hz <- L$metadata$heterogeneity; gz <- function(k, d) if (is.null(hz) || is.null(hz[[k]]) || !length(hz[[k]])) d else hz[[k]][1]
  u$het_void_frac <- gz("void_frac", 0); u$het_void_ell_um <- gz("void_ell_um", NA_real_); u$het_void_aspect <- gz("void_aspect", 1)
  u$het_void_mode <- gz("void_mode", "move"); u$het_dens_sigma <- gz("dens_sigma", 0); u$het_dens_ell_um <- gz("dens_ell_um", NA_real_); u$het_dens_cap <- gz("dens_cap", FALSE)
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
         "results/R2b/annotations/R2b_annotations_Luca.json", "results/R2c/annotations/R2c_answers_Luca.json", "results/R2c/R2c_focal.csv",
         "results/S1.5/S1.5_real_refs_roi.csv", "results/S1.5b/S1.5b_real_refs_roi_type.csv", "results/S1.5c/S1.5c_real_intervals.csv",
         "results/S1.5d/S1.5d_panel_units.csv", list.files("results/S1.5d", pattern = "^S1.5d_(compartments_A|het_RV\\.csv|real_)", full.names = TRUE),
         list.files("results/S1.6", pattern = "^S1.6_real_", full.names = TRUE), "results/S1.5e/S1.5e_panel_units.csv",
         "inst/extdata/size_resid.csv")
  data.frame(file = f, md5 = unname(tools::md5sum(f)))
}
if (STAGE == "md5before") { write.csv(md5_inputs(), file.path(OUT, "md5_before.csv"), row.names = FALSE); cat("md5 prima:", nrow(md5_inputs()), "file\n") }
if (STAGE == "md5after") {
  b <- read.csv(file.path(OUT, "md5_before.csv")); a <- md5_inputs(); write.csv(a, file.path(OUT, "md5_after.csv"), row.names = FALSE)
  k <- match(b$file, a$file); same <- sum(a$md5[k] == b$md5, na.rm = TRUE)
  cat(sprintf("%s md5 dei dati letti invariati: %d/%d\n", if (same == nrow(b)) "PASS" else "FAIL", same, nrow(b)))
  if (same < nrow(b)) print(b$file[is.na(k) | a$md5[k] != b$md5])
}
if (STAGE == "gen") {
  for (arm in ARMS$arm) for (d in c("gen", "cells", "nuc")) dir.create(file.path(OUT, d, arm), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUT, "layers"), showWarnings = FALSE)
  GEN_ARMS <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else stop("gen: indicare i bracci (U6,Q6)")
  stopifnot(all(GEN_ARMS %in% ARMS$arm))
  jobs <- expand.grid(i = rois$roi_idx, s = SEEDS, arm = GEN_ARMS, stringsAsFactors = FALSE)
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
  dir.create(RES, showWarnings = FALSE)
  write.csv(do.call(rbind, lapply(z, `[[`, "unit")), file.path(RES, "S1.6_panel_units.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(z, `[[`, "type")), file.path(RES, "S1.6_panel_types.csv"), row.names = FALSE)
  cat("unita':", length(z), "\n")
}
if (STAGE == "check") {
  U <- read.csv(file.path(RES, "S1.6_panel_units.csv"), stringsAsFactors = FALSE)
  S <- read.csv("results/S1.5e/S1.5e_panel_units.csv", stringsAsFactors = FALSE)
  e1 <- sum(grepl("C2:ERROR|C3:ERROR|C5:ERROR|C9:ERROR", U$statuses))
  cmp <- function(a, ref) { uu <- U[U$arm == a, ]; ss <- S[S$arm == ref, ]; k <- match(paste(uu$archetype, uu$roi_id, uu$seed), paste(ss$archetype, ss$roi_id, ss$seed))
    sum(uu$md5 == ss$md5[k], na.rm = TRUE) }
  q6 <- U[U$arm == "Q6", ]
  pick <- do.call(rbind, lapply(split(rois, rois$archetype), function(d) d[1, ]))
  rows8 <- pick[rep(seq_len(nrow(pick)), length.out = 8), ]
  rr <- function(arms, cmpf) vapply(seq_along(arms), function(r) { i <- rows8$roi_idx[r]; arm <- arms[r]
    w <- r3_roi_window(rows8$archetype[r], rows8$roi_id[r], rois); cmpf(layer_md5(run_pc(i, 3, arm, w)), i, arm) }, TRUE)
  rep_ok <- rr(rep(c("U6", "Q6"), 4), function(m, i, arm) identical(m, U$md5[U$roi_idx == i & U$seed == 3 & U$arm == arm]))
  leg_ok <- rr(rep(c("RVl", "Ql"), 4), function(m, i, arm) identical(m, S$md5[S$roi_idx == i & S$seed == 3 & S$arm == ARMS$ref[ARMS$arm == arm]]))
  c7w <- sum(grepl("C7:WARN", U$statuses[U$arm == "U6"])); c7w_old <- sum(grepl("C7:WARN", S$statuses[S$arm == "U"]))
  out <- data.frame(check = c("units", "no_error_C2_C3_C5_C9", "U6_equals_S1.5e_U_md5", "Q6_exceed_wmax", "Q6_regions_grid_bound", "Q6_quotas_on",
                              "rerun_identical", "legacy_RV_Q_of_S1.5e_md5", "desc_Q6_n_failed", "desc_C7_warn_U6", "desc_C7_warn_U_S1.5e"),
                    value = c(nrow(U), nrow(U) - e1, cmp("U6", "U"), sum(q6$n_exceed_wmax), sum(q6$n_regions_grid_bound), sum(q6$quotas_on),
                              sum(rep_ok), sum(leg_ok), sum(q6$n_failed), c7w, c7w_old),
                    expected = c(600, nrow(U), 300, 0, 0, 0, 8, 8, NA, NA, NA))
  out$pass <- ifelse(is.na(out$expected), NA, out$value == out$expected)
  write.csv(out, file.path(RES, "S1.6_G3_check.csv"), row.names = FALSE); print(out)
  cat(sprintf("%s G-3 (md5 dei dati: stadio md5after)\n", if (all(out$pass, na.rm = TRUE)) "PASS" else "FAIL"))
}
