# tools/S1.4f_bench.R - S1.4f: banco del nucleo orientato in parte (taratura delle leggi V0, V1, V2; bracci K2, ET, EP). Attese: docs/steps/S1.4f.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 Rscript --vanilla tools/S1.4f_bench.R <stadio>
#   md5before | calB | genB | genA | check | md5after   (stadi rilanciabili: saltano le unita' fatte)
# calB: semi 1-5 del banco B (150 unita'), place_cells() come S1.4e, nuclei ridisegnati per 3 leggi x 12 valori di lambda (S14F_GRID).
# genB: le 300 unita'; K2 ed ET ridisegnati (identita' con S1.4e) e V0, V1, V2 con il lambda della piega lascia-fuori-un-ROI (S1.4f_lambda.csv).
# genA: le 24 finestre di R2b (centri come S1.4e) per K2, ET, EP (legge vincente della regola O), e le stesse senza i doppioni al bordo (BL-174).
# Uscite stabili in /mnt/micron/geo_spatialtrans/S1.4f (un file per unita'). Worker dal file WORKERS (tetto 28).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
SRC <- Sys.getenv("S14F_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5c_bench_lib.R"); source("tools/S1.5d_bench_lib.R"); source("tools/S1.4e_bench_lib.R")
source(file.path(SRC, "tools/S1.4f_bench_lib.R"))
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); STAGE <- if (length(args)) args[1] else "none"
M <- "/mnt/micron/geo_spatialtrans"; OUT <- file.path(M, "S1.4f"); E14 <- file.path(M, "S1.4e"); RES <- "results/S1.4f"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
wk_file <- file.path(OUT, "WORKERS"); if (!file.exists(wk_file)) writeLines("20", wk_file)
NW <- function() min(28L, as.integer(readLines(wk_file)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
ids_of <- function(A) rois$roi_id[rois$archetype == A]
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
SEEDS <- 1:10; SEEDS_CAL <- 1:5; BASE <- 20261021; seed_of <- function(i, s) BASE + 1000 * i + s            # come S1.5d, S1.5e, S1.4e
WIN <- read.csv("results/R2b/R2b_windows.csv", stringsAsFactors = FALSE)
LAWS <- names(S14F_LAWS)
clust_of <- function(w) { idx <- which(w$m, arr.ind = TRUE); data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L) }
layer_md5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  on.exit(unlink(f)); unname(tools::md5sum(f)) }
wkb_list <- function(nuc) lapply(as.list(geos::geos_write_wkb(geos::as_geos_geometry(nuc))), as.vector)
same_as_file <- function(nuc, f) { g0 <- read_parquet(f); identical(wkb_list(nuc), lapply(as.list(g0$wkb), as.vector)) }
write_nuc <- function(f, nuc) { g <- geos::as_geos_geometry(nuc)
  write_parquet(data.frame(idx = seq_along(g), wkb = I(geos::geos_write_wkb(g))), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f) }
run_pc_RV <- function(i, s, w) {       # = run_pc_RV() di tools/S1.4e_bench.R
  A <- rois$archetype[i]; roi <- rois$roi_id[i]
  suppressWarnings(place_cells(clust_of(w), PRESET[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                               small_components = "exclude", random_seed = seed_of(i, s), types = "multi", size_arm = "group",
                               spatial = "compartments", min_dist_mode = "nucleus", override_compartments = compart_for_d(A, roi),
                               override_heterogeneity = het_for(A, roi, "RV"), quotas = TRUE, size_resid = "gaussian", resid_mode = "decile",
                               keep_tiles = TRUE, strict = FALSE, verbose = FALSE))
}
interior_of <- function(L, w) {
  td <- L$territory_df; in_frame <- geos::geos_within(geos::as_geos_geometry(L$tiles), geos::as_geos_geometry(w$frame))
  !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
}
par_run <- function(jobs, fun, label) {
  t0 <- Sys.time(); r <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch({ fun(jobs[k, , drop = FALSE]); "ok" }, error = function(e) paste("ERR", conditionMessage(e))),
                                  mc.cores = NW(), mc.preschedule = FALSE)
  r <- unlist(r); cat(sprintf("[%s] %d unita', %d errori, %.1f min\n", label, length(r), sum(grepl("^ERR", r)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (any(grepl("^ERR", r))) { print(head(cbind(jobs[grepl("^ERR", r), , drop = FALSE], err = r[grepl("^ERR", r)]), 10)); quit(status = 1) }
}
lambda_of <- function(LT, law, A, roi) { v <- LT$lambda[LT$law == law & LT$archetype == A & LT$roi_id == roi]; stopifnot(length(v) == 1); if (is.na(v)) Inf else v }
md5_inputs <- function() {
  f <- c(list.files("results/R2/tissue_masks", full.names = TRUE), "results/R2b/R2b_windows.csv", "results/R2b/R2b_points.csv", "results/R2c/R2c_dmax.csv",
         "results/R2b/annotations/R2b_annotations_Luca.json", "results/R2c/annotations/R2c_answers_Luca.json", "results/R2c/R2c_items_Luca.csv",
         "results/R2c/R2c_focal.csv", "results/R2c/R2c_pairs.csv", "results/R2/R2_rois_checked.csv",
         list.files("results/R2/nuclei", pattern = "_(cellpose_rgb|stardist_he)_native\\.parquet$", full.names = TRUE),
         "results/S1.5d/S1.5d_panel_units.csv", list.files("results/S1.5d", pattern = "^S1.5d_(compartments_A|het_RV\\.csv)", full.names = TRUE),
         "inst/extdata/size_resid.csv", "docs/BIO_REFERENCES.md", list.files(file.path(M, "S1.4/real"), pattern = "cellpose_rgb\\.rds$", full.names = TRUE),
         list.files(file.path(M, "S1.4/loo"), pattern = "cellpose_rgb_fold", full.names = TRUE), list.files(file.path(M, "R3b/S1.4/real"), pattern = "consensus\\.rds$", full.names = TRUE),
         list.files(file.path(M, "R3b/S1.4/loo"), pattern = "consensus_fold", full.names = TRUE), file.path(M, "R3b/phi.csv"),
         list.files("results/S1.4e", pattern = "\\.csv$", full.names = TRUE), "results/S1.4e/explore/S1.4e_E5_nuclei.csv.gz",
         list.files(file.path(E14, "B/nuc"), recursive = TRUE, full.names = TRUE), list.files(file.path(E14, "A/nuc"), recursive = TRUE, full.names = TRUE),
         list.files(file.path(E14, "B/raster"), recursive = TRUE, full.names = TRUE), list.files(file.path(E14, "B/pairs"), recursive = TRUE, full.names = TRUE),
         list.files(file.path(E14, "A/pairs"), recursive = TRUE, full.names = TRUE))
  data.frame(file = f, md5 = unname(tools::md5sum(f)))
}
if (STAGE == "md5before") { x <- md5_inputs(); write.csv(x, file.path(OUT, "md5_before.csv"), row.names = FALSE); cat("md5 prima:", nrow(x), "file\n") }
if (STAGE == "md5after") {
  b <- read.csv(file.path(OUT, "md5_before.csv")); a <- md5_inputs(); write.csv(a, file.path(OUT, "md5_after.csv"), row.names = FALSE)
  k <- match(b$file, a$file); same <- sum(a$md5[k] == b$md5, na.rm = TRUE)
  cat(sprintf("%s md5 dei dati letti invariati: %d/%d\n", if (same == nrow(b)) "PASS" else "FAIL", same, nrow(b)))
  if (same < nrow(b)) print(b$file[is.na(k) | a$md5[k] != b$md5])
}
# ---- taratura ---------------------------------------------------------------------------------------------------------------
if (STAGE == "calB") {
  for (d in c("cal/units", "cal/cells", file.path("cal/nuc", outer(LAWS, sprintf("k%02d", seq_along(S14F_GRID)), paste, sep = "/"))))
    dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
  jobs <- expand.grid(i = rois$roi_idx, s = SEEDS_CAL); jobs <- jobs[order(jobs$s, -(rois$archetype[jobs$i] %in% c("A4", "A1"))), ]
  par_run(jobs, function(j) {
    i <- j$i; s <- j$s; A <- rois$archetype[i]; roi <- rois$roi_id[i]; key <- sprintf("%s_%s_s%02d", A, roi, s)
    f <- file.path(OUT, "cal", "units", paste0(key, ".rds")); if (file.exists(f)) return(invisible())
    w <- r3_roi_window(A, roi, rois); L <- run_pc_RV(i, s, w); it <- interior_of(L, w); cd <- L$cell_df
    write_parquet(data.frame(cell_id = cd$cell_id, x = cd$x, y = cd$y, interior = it), file.path(OUT, "cal", "cells", paste0(key, ".parquet")))
    rows <- list(); md5 <- layer_md5(L)
    for (law in LAWS) for (k in seq_along(S14F_GRID)) {
      o <- redraw_nuclei(L, law_catalog(L$cell_types_resolved, A, law, S14F_GRID[k]))
      write_nuc(file.path(OUT, "cal", "nuc", law, sprintf("k%02d", k), paste0(key, ".parquet")), o$cell_nuclei)
      rows[[length(rows) + 1]] <- cbind(data.frame(archetype = A, roi_id = roi, roi_idx = i, seed = s, law = law, k = k, lambda = S14F_GRID[k], md5_layer = md5),
                                        nuc_invariants(o, L$cell_territories))
    }
    saveRDS(do.call(rbind, rows), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  }, "calB")
  z <- do.call(rbind, lapply(list.files(file.path(OUT, "cal", "units"), pattern = "rds$", full.names = TRUE), readRDS))
  write.csv(z, file.path(RES, "S1.4f_cal_units.csv"), row.names = FALSE); cat("righe:", nrow(z), "\n")
}
# ---- banco B ---------------------------------------------------------------------------------------------------------------
if (STAGE == "genB") {
  LT <- read.csv(file.path(RES, "S1.4f_lambda.csv"), stringsAsFactors = FALSE)
  for (d in c("B/units", "B/cells", paste0("B/nuc/", c("K2", "ET", LAWS)), paste0("B/g4/", LAWS))) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
  jobs <- expand.grid(i = rois$roi_idx, s = SEEDS); jobs <- jobs[order(jobs$s, -(rois$archetype[jobs$i] %in% c("A4", "A1"))), ]
  par_run(jobs, function(j) {
    i <- j$i; s <- j$s; A <- rois$archetype[i]; roi <- rois$roi_id[i]; key <- sprintf("%s_%s_s%02d", A, roi, s)
    f <- file.path(OUT, "B", "units", paste0(key, ".rds")); if (file.exists(f)) return(invisible())
    t0 <- proc.time()[[3]]; w <- r3_roi_window(A, roi, rois); L <- run_pc_RV(i, s, w); t_pc <- proc.time()[[3]] - t0
    it <- interior_of(L, w); cd <- L$cell_df
    write_parquet(data.frame(cell_id = cd$cell_id, x = cd$x, y = cd$y, cell_type = as.character(cd$cell_type), interior = it,
                             territory_area = cd$territory_area), file.path(OUT, "B", "cells", paste0(key, ".parquet")))
    rows <- list(); md5 <- layer_md5(L)
    for (arm in c("K2", "ET", LAWS)) {
      ct <- if (arm %in% c("K2", "ET")) arm_catalog(L$cell_types_resolved, arm, A) else law_catalog(L$cell_types_resolved, A, arm, lambda_of(LT, arm, A, roi))
      t1 <- proc.time()[[3]]; o <- redraw_nuclei(L, ct); t_d <- proc.time()[[3]] - t1
      same <- if (arm %in% c("K2", "ET")) same_as_file(o$cell_nuclei, file.path(E14, "B", "nuc", arm, paste0(key, ".parquet"))) else NA
      write_nuc(file.path(OUT, "B", "nuc", arm, paste0(key, ".parquet")), o$cell_nuclei)
      if (s == 1 && arm %in% LAWS) {                                   # G-4: dati per la seconda implementazione
        tp <- .dcg_as_polys(L$cell_territories); g <- geos::as_geos_geometry(tp); od <- o$cell_df; mo <- .dcg_territory_moments(tp)
        w_pkg <- .dcg_orient_weight(rep(S14F_LAWS[[arm]], nrow(od)), mo$aspect, od$nucleus_area_target, od$center_dist)
        write_parquet(data.frame(cell_id = od$cell_id, x = od$x, y = od$y, nucleus_orient = od$nucleus_orient, center_dist = od$center_dist,
                                 terr_axis_pkg = mo$theta, terr_aspect_pkg = mo$aspect, w_pkg = w_pkg,
                                 nucleus_area_target = od$nucleus_area_target, lambda = lambda_of(LT, arm, A, roi), law = S14F_LAWS[[arm]],
                                 aspect = ct$nucleus_aspect[match(od$cell_type, ct$cell_type)], terr_wkb = I(geos::geos_write_wkb(g))),
                      file.path(OUT, "B", "g4", arm, paste0(key, ".parquet")))
      }
      rows[[arm]] <- cbind(data.frame(archetype = A, roi_id = roi, roi_idx = i, seed = s, arm = arm, lambda = if (arm %in% LAWS) lambda_of(LT, arm, A, roi) else NA,
                                      md5_layer = md5, same_as_S14e = same, t_place_s = t_pc, t_redraw_s = t_d, n_interior = sum(it)),
                           nuc_invariants(o, L$cell_territories))
    }
    saveRDS(do.call(rbind, rows), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  }, "genB")
  z <- do.call(rbind, lapply(list.files(file.path(OUT, "B", "units"), pattern = "rds$", full.names = TRUE), readRDS))
  write.csv(z, file.path(RES, "S1.4f_B_units.csv"), row.names = FALSE); cat("righe:", nrow(z), "\n")
}
# ---- banco A ---------------------------------------------------------------------------------------------------------------
SRC_A <- function(A) if (A %in% c("A4", "A6")) list(method = "consensus", base = file.path(M, "R3b", "S1.4")) else list(method = "cellpose_rgb", base = file.path(M, "S1.4"))
PHI_CP <- c(A1 = 0.757, A2 = 1, A3 = 0.947, A5 = 0.474)                                     # come S1.4e
phi_A <- function(A) { if (A %in% names(PHI_CP)) return(PHI_CP[[A]])
  p <- read.csv(file.path(M, "R3b", "phi.csv")); min(1, round(p$phi[p$archetype == A & p$method == "consensus"], 3)) }
ct_m2 <- function(f2, phi) data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]],
  nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]], nucleus_nc_max = phi)
one_region <- function(w) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2), region_polygons = st_sfc(w$poly[[1]]))
if (STAGE == "genA") {
  LT <- read.csv(file.path(RES, "S1.4f_lambda.csv"), stringsAsFactors = FALSE); WINNER <- read.csv(file.path(RES, "S1.4f_ruleO.csv"))$winner[1]
  stopifnot(WINNER %in% LAWS); cat("legge di EP (regola O):", WINNER, "\n")
  for (b in c("A", "A_dd")) for (d in c("map", "units", paste0("nuc/", c("K2", "ET", "EP")))) dir.create(file.path(OUT, b, d), recursive = TRUE, showWarnings = FALSE)
  ann <- jsonlite::fromJSON("results/R2b/annotations/R2b_annotations_Luca.json", simplifyVector = FALSE)$windows
  pts_all <- read.csv("results/R2b/R2b_points.csv", stringsAsFactors = FALSE); pts_all <- pts_all[pts_all$rater == "Luca", ]
  jobs <- expand.grid(win_id = WIN$win_id, bench = c("A", "A_dd"), stringsAsFactors = FALSE)
  par_run(jobs, function(j) {
    wv <- WIN[WIN$win_id == j$win_id, ]; A <- wv$archetype; roi <- wv$roi_id; sa <- SRC_A(A); bench <- j$bench
    f <- file.path(OUT, bench, "units", paste0(wv$win_id, ".rds")); if (file.exists(f)) return(invisible())
    k <- match(roi, ids_of(A)); o14 <- readRDS(file.path(sa$base, "loo", sprintf("%s_%s_fold%d.rds", A, sa$method, k))); stopifnot(o14$roi_test == roi)
    w <- r3_roi_window(A, roi, rois); te <- readRDS(file.path(sa$base, "real", sprintf("%s_%s_%s.rds", A, roi, sa$method)))
    pts <- pts_all[pts_all$win_id == wv$win_id, ]
    cc <- centers_A(wv, te$cells[, c("x", "y")], pts, w, ann); n0 <- nrow(cc)
    if (bench == "A_dd") cc <- s14f_dedup(cc)
    n <- nrow(cc); cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c")
    tv <- tessellate_voronoi(cen, one_region(w), verbose = FALSE); terr <- .dcg_as_polys(tv$cell_territories)
    tess <- list(cell_territories = terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = tv$territory_df$territory_area))
    write_parquet(cbind(data.frame(cell_id = seq_len(n)), cc), file.path(OUT, bench, "map", paste0(wv$win_id, ".parquet")))
    ct0 <- ct_m2(o14$fits$M2, phi_A(A)); rows <- list()
    for (s_ in SEEDS) {
      sd <- o14$seed + 100000 * s_
      for (arm in c("K2", "ET", "EP")) {
        ct <- if (arm == "EP") law_catalog(ct0, A, WINNER, lambda_of(LT, WINNER, A, roi)) else arm_catalog(ct0, arm, A)
        o <- derive_cell_geometry(cen, tess, ct, random_seed = sd, verbose = FALSE)
        fn <- sprintf("%s_%s_%s_s%02d.parquet", A, roi, wv$win_id, s_); write_nuc(file.path(OUT, bench, "nuc", arm, fn), o$cell_nuclei)
        same <- if (bench == "A" && arm != "EP") same_as_file(o$cell_nuclei, file.path(E14, "A", "nuc", arm, fn)) else NA
        rows[[length(rows) + 1]] <- cbind(data.frame(bench = bench, win_id = wv$win_id, archetype = A, roi_id = roi, seed_idx = s_, seed = sd, arm = arm, n = n,
                                                     n_removed_dd = n0 - n, n_manual = sum(cc$src == "manual"), same_as_S14e = same), nuc_invariants(o, terr))
      }
    }
    saveRDS(do.call(rbind, rows), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  }, "genA")
  z <- do.call(rbind, lapply(c(list.files(file.path(OUT, "A", "units"), pattern = "rds$", full.names = TRUE),
                               list.files(file.path(OUT, "A_dd", "units"), pattern = "rds$", full.names = TRUE)), readRDS))
  write.csv(z, file.path(RES, "S1.4f_A_units.csv"), row.names = FALSE); cat("righe:", nrow(z), "\n")
}
# ---- check: cancello G-3 (parte R) ------------------------------------------------------------------------------------------
if (STAGE == "check") {
  C <- read.csv(file.path(RES, "S1.4f_cal_units.csv"), stringsAsFactors = FALSE); B <- read.csv(file.path(RES, "S1.4f_B_units.csv"), stringsAsFactors = FALSE)
  Aa <- read.csv(file.path(RES, "S1.4f_A_units.csv"), stringsAsFactors = FALSE)
  S <- read.csv("results/S1.5d/S1.5d_panel_units.csv", stringsAsFactors = FALSE); S <- S[S$arm == "RV", ]
  bk <- B[B$arm == "K2", ]; k <- match(paste(bk$archetype, bk$roi_id, bk$seed), paste(S$archetype, S$roi_id, S$seed))
  ck_ <- C[!duplicated(paste(C$archetype, C$roi_id, C$seed)), ]; kc <- match(paste(ck_$archetype, ck_$roi_id, ck_$seed), paste(S$archetype, S$roi_id, S$seed))
  inv <- function(d) sum(d$n_out + d$n_overlap + d$n_area_bad + d$n_invalid + d$n_empty)
  out <- data.frame(check = c("B_units", "B_RV_equals_S1.5d_RV_md5", "cal_RV_equals_S1.5d_RV_md5", "B_K2_equals_S1.4e", "B_ET_equals_S1.4e",
                              "A_K2_equals_S1.4e", "A_ET_equals_S1.4e", "cal_rows", "A_rows", "invariant_violations_cal", "invariant_violations_B", "invariant_violations_A"),
                    value = c(nrow(bk), sum(bk$md5_layer == S$md5[k], na.rm = TRUE), sum(ck_$md5_layer == S$md5[kc], na.rm = TRUE),
                              sum(B$same_as_S14e[B$arm == "K2"], na.rm = TRUE), sum(B$same_as_S14e[B$arm == "ET"], na.rm = TRUE),
                              sum(Aa$same_as_S14e[Aa$bench == "A" & Aa$arm == "K2"], na.rm = TRUE), sum(Aa$same_as_S14e[Aa$bench == "A" & Aa$arm == "ET"], na.rm = TRUE),
                              nrow(C), nrow(Aa), inv(C), inv(B), inv(Aa)),
                    expected = c(300, 300, 150, 300, 300, 240, 240, 150 * 3 * length(S14F_GRID), 2 * 24 * 10 * 3, 0, 0, 0))
  out$pass <- out$value == out$expected
  write.csv(out, file.path(RES, "S1.4f_G3_check.csv"), row.names = FALSE); print(out)
  cat(sprintf("%s G-3 (md5 dei dati: stadio md5after)\n", if (all(out$pass)) "PASS" else "FAIL"))
}
