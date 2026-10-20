# tools/R3b_bench.R — R3b: catena della maschera unica per linfonodo (A4) e cuore (A6). Generatori e territori (R3), nucleo
# reale per cellula (S1.4), modello di taglia M2 per fold, tetto phi, nucleo ritagliato e gonfiato (S1.4d), controprova CP-1
# (punti manuali di Luca come centri). Pre-registrazione: docs/steps/R3b.md (cancelli 48573f5, ok di Luca ea1e22a).
# Le funzioni sono copie di tools/R3_run.R (real_one), tools/S1.4_bench.R (real_one, fit_set) e tools/S1.4d_bench.R (gen_one,
# invariants) con i percorsi parametrizzati; il cancello G-3 (stadio regress) verifica che con la maschera StarDist riproducano
# i file di R3, S1.4 e S1.4d.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 Rscript --vanilla tools/R3b_bench.R <stage>
#   stage = r3 | real | loo | phi | gen | cp1 | regress ; uscite in /mnt/micron/geo_spatialtrans/R3b ; unita' gia' scritte saltate.
#   Worker: file /mnt/micron/geo_spatialtrans/R3b/WORKERS (letto a ogni blocco; default 6; server condiviso con MMM).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow); library(spatstat.geom)
                                 library(spatstat.explore) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
source("tools/R3_voronoi_metrics.R"); source("tools/S1.4d_mutants.R")
E <- dcg_env_d()
STAGE <- commandArgs(TRUE)[1]
M <- "/mnt/micron/geo_spatialtrans"; OUT <- file.path(M, "R3b")
R3 <- file.path(M, "R3"); S13 <- file.path(M, "S1.3"); S14 <- file.path(M, "S1.4"); S14D <- file.path(M, "S1.4d")
for (d in c("R3/real", "S1.4/real", "S1.4/loo", "gen", "cp1", "regress")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
BASE14 <- 20261015; NSEED <- 10L; TIS <- c("A4", "A6")
PHI_SR <- c(A4 = 0.816, A6 = 0.180)          # phi di S1.4d (B-124, B-126), usato solo per la regressione G-3 di StarDist
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
WIN <- read.csv("results/R2b/R2b_windows.csv"); WIN <- WIN[WIN$archetype %in% TIS, ]
ids_of <- function(A) rois$roi_id[rois$archetype == A]
one_region <- function(w) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2),
                               region_polygons = st_sfc(w$poly[[1]]))

#' Esecuzione a blocchi con numero di worker letto da OUT/WORKERS (come tools/S1.4d_bench.R).
par_run <- function(jobs, fun, label) {
  idx <- seq_len(nrow(jobs)); errs <- character(0)
  while (length(idx)) {
    w <- suppressWarnings(as.integer(readLines(file.path(OUT, "WORKERS"), warn = FALSE)[1])); if (is.na(w) || w < 1) w <- 6L
    blk <- idx[seq_len(min(length(idx), w))]; idx <- setdiff(idx, blk)
    r <- mclapply(blk, function(i) { t0 <- proc.time()[[3]]
      x <- tryCatch(do.call(fun, as.list(jobs[i, , drop = FALSE])), error = function(e) paste("ERRORE:", conditionMessage(e)))
      cat(sprintf("[%s] %s %.0fs %s\n", label, paste(unlist(jobs[i, ]), collapse = " "), proc.time()[[3]] - t0, if (is.character(x)) x else "ok")); x },
      mc.cores = w, mc.preschedule = FALSE)
    errs <- c(errs, unlist(Filter(is.character, r)))
  }
  cat(sprintf("[%s] %d unita', %d errori\n", label, nrow(jobs), length(errs))); if (length(errs)) { print(errs); quit(status = 1) }
}

# ---- r3: generatori reali e territori di Voronoi (copia di tools/R3_run.R::real_one) ----------------------------------
load_nuc <- function(A, roi, method) {
  cols <- c("archetype", "roi_id", "method", "scale", "keep", "label", "x_um", "y_um", "area_um2", "eccentricity", "orientation", "um_per_px")
  f <- if (method == "consensus") file.path(OUT, "nuclei", sprintf("%s_%s_consensus_native.parquet", A, roi)) else
    sprintf("results/R2/nuclei/%s_%s_%s_native.parquet", A, roi, method)
  d <- as.data.frame(read_parquet(f, col_select = cols)); d[d$scale == 1 & d$keep, ]
}
r3_one <- function(A, roi, method) {
  RO <- file.path(OUT, "R3", "real")
  f_cells <- file.path(RO, sprintf("%s_%s_%s_cells.parquet", A, roi, method)); if (file.exists(f_cells)) return(invisible(NULL))
  w <- r3_roi_window(A, roi, rois)
  d <- load_nuc(A, roi, method); n_keep <- nrow(d)
  d <- d[r3_in_mask(d$x_um, d$y_um, w$m, w$px), ]
  gen <- data.frame(idx = seq_len(nrow(d)), label = d$label, x = d$x_um, y = d$y_um, area_nuc = d$area_um2,
                    ecc_parquet = d$eccentricity, orient_parquet = d$orientation, um_per_px = d$um_per_px)
  write_parquet(gen, file.path(RO, sprintf("%s_%s_%s_gen.parquet", A, roi, method)))
  tv <- r3_tessellate(gen$x, gen$y, w$poly, w$frame)
  li <- r3_local_intensity(gen$x, gen$y, w)
  tv$lambda_loc <- li$lambda_loc; tv$label <- gen$label; tv$area_nuc <- gen$area_nuc
  X <- li$X; nn <- nndist(X); bd <- bdist.points(X)
  ce <- tryCatch(clarkevans(X, correction = "cdf"), error = function(e) NA_real_)
  dip <- tryCatch(diptest::dip.test(tv$area[tv$interior])$p.value, error = function(e) NA_real_)
  write_parquet(tv, f_cells)
  roi_row <- data.frame(archetype = A, roi_id = roi, method = method, n_keep = n_keep, n = nrow(gen),
                        area_poly_um2 = w$area_um2, area_px_um2 = w$n_valid_px * w$px^2, sum_area = sum(tv$area),
                        gen_in_own = mean(tv$gen_in_own), sigma_loc_um = li$sigma, lambda_mm2 = li$lambda * 1e6,
                        nn_median_minus = median(nn[nn < bd]), nn_n_minus = sum(nn < bd), clark_evans_cdf = as.numeric(ce),
                        dip_p_area = dip, t(r3_summary(tv, tv$lambda_loc)))
  saveRDS(roi_row, file.path(RO, sprintf("%s_%s_%s_roi.rds", A, roi, method)))
  invisible(NULL)
}

# ---- real: territori del pacchetto + nucleo reale per cellula (copia di tools/S1.4_bench.R::real_one) ----------------
rp <- function(A, roi, method, k) file.path(OUT, "R3", "real", sprintf("%s_%s_%s_%s.parquet", A, roi, method, k))
real_one <- function(A, roi, method) {
  f <- file.path(OUT, "S1.4", "real", sprintf("%s_%s_%s.rds", A, roi, method)); if (file.exists(f)) return(invisible(NULL))
  w <- r3_roi_window(A, roi, rois)
  gen <- as.data.frame(read_parquet(rp(A, roi, method, "gen"))); nucp <- as.data.frame(read_parquet(rp(A, roi, method, "nuc")))
  stopifnot(identical(nucp$idx, gen$idx)); n <- nrow(gen)
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = gen$x, y = gen$y, cell_type = "c")
  tv <- tessellate_voronoi(cen, one_region(w), keep_tiles = TRUE, verbose = FALSE)
  td <- tv$territory_df; terr <- E$.dcg_as_polys(tv$cell_territories)
  in_frame <- lengths(st_within(tv$tiles, w$frame)) > 0
  interior <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
  f13 <- file.path(S13, "real", sprintf("%s_%s_%s.rds", A, roi, method))
  c0 <- if (file.exists(f13)) { ref13 <- readRDS(f13)$geos$cells
    c(rel_area_vs_S13 = max(abs(td$territory_area - ref13$area) / ref13$area), interior_identical_S13 = identical(interior, ref13$interior))
  } else c(rel_area_vs_S13 = NA_real_, interior_identical_S13 = NA)
  a_own <- gen$area_nuc
  matched <- rep(TRUE, n); a_n <- a_own; ecc_n <- nucp$ecc_N; th_n <- nucp$theta_N; match_info <- NULL   # maschera unica
  gt <- as_geos_geometry(terr)
  p <- polylabelr::poi(terr, precision = 1e-3); px <- vapply(p, function(q) q$x, 0); py <- vapply(p, function(q) q$y, 0)
  bd <- geos_boundary(gt)
  cells <- data.frame(idx = seq_len(n), x = gen$x, y = gen$y, a_t = td$territory_area, tile_area = td$tile_area,
                      interior = interior, matched = matched, a_n = a_n, a_own = a_own, ecc_n = ecc_n, theta_n = th_n,
                      theta_t = E$.dcg_territory_axis(terr), poi_x = px, poi_y = py,
                      d_poi = geos_distance(geos_make_point(px, py), bd), d_gen = geos_distance(geos_make_point(gen$x, gen$y), bd))
  saveRDS(list(archetype = A, roi_id = roi, method = method, n = n, area_window = w$area_um2, c0 = c0,
               match_info = match_info, cells = cells, terr = terr), f)
  invisible(NULL)
}
load_real <- function(A, roi, method, base = file.path(OUT, "S1.4")) readRDS(file.path(base, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
evset <- function(cc) cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
fit_set <- function(rs, model = "full") {
  a_ref <- sum(vapply(rs, function(r) sum(r$cells$a_t), 0)) / sum(vapply(rs, function(r) r$n, 0))
  cc <- do.call(rbind, lapply(rs, function(r) r$cells[evset(r$cells), ]))
  list(fit = E$fit_nucleus_size_model(cc$a_n, cc$a_t, a_ref, model), a_ref = a_ref, n_ev = nrow(cc))
}

# ---- loo: modello di taglia M2 stimato sugli altri 4 ROI (come S1.4, solo M2) ----------------------------------------
loo_one <- function(A, method, k) {
  f <- file.path(OUT, "S1.4", "loo", sprintf("%s_%s_fold%d.rds", A, method, k)); if (file.exists(f)) return(invisible(NULL))
  ids <- ids_of(A); ai <- match(A, paste0("A", 1:6))
  rs <- lapply(ids, function(r) load_real(A, r, method)); fs <- fit_set(rs[-k], "full")
  saveRDS(list(archetype = A, method = method, fold = k, roi_test = ids[k], seed = BASE14 + 1000 * ai + 10 * k,
               fits = list(M2 = c(fs$fit, a_ref = fs$a_ref)), n_ev_train = fs$n_ev, n_ev_test = sum(evset(rs[[k]]$cells))), f)
  invisible(NULL)
}

# ---- phi: q99 del rapporto area del nucleo / area del territorio delle celle interne (come B-121...B-126) -------------
phi_of <- function(A, method, base = file.path(OUT, "S1.4")) {
  nc <- unlist(lapply(ids_of(A), function(r) { cc <- load_real(A, r, method, base)$cells; e <- evset(cc); cc$a_n[e] / cc$a_t[e] }))
  c(phi = unname(quantile(nc, 0.99, type = 7)), n = length(nc))
}
PHI <- function(A, method) { p <- read.csv(file.path(OUT, "phi.csv")); min(1, round(p$phi[p$archetype == A & p$method == method], 3)) }

# ---- gen: nucleo ritagliato e gonfiato (copia di tools/S1.4d_bench.R::gen_one) ---------------------------------------
ct_m2 <- function(f2, phi = NULL) { ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]],
  nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]]); if (!is.null(phi)) ct$nucleus_nc_max <- phi; ct }
inputs_of <- function(te) { n <- te$n; cc <- te$cells
  list(cen = data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c"),
       tess = list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t))) }
wkb_md5 <- function(o) { f <- tempfile(); saveRDS(list(o$cell_df, geos_write_wkb(as_geos_geometry(o$cell_nuclei))), f, compress = FALSE)
  m <- unname(tools::md5sum(f)); unlink(f); m }
write_nuc <- function(f, g) { if (file.exists(f)) return(invisible(NULL))
  write_parquet(data.frame(idx = seq_along(g), wkb = I(geos_write_wkb(g))), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f) }
invariants <- function(o, terr) {
  g <- as_geos_geometry(o$cell_nuclei); gt <- as_geos_geometry(E$.dcg_as_polys(terr)); cd <- o$cell_df; n <- nrow(cd)
  w2 <- st_within(o$cell_nuclei, terr); c3_sf <- sum(!vapply(seq_len(n), function(i) i %in% w2[[i]], TRUE))
  tr <- geos_strtree(g); q <- geos_strtree_query(tr, g); pi_ <- rep(seq_along(q), lengths(q)); pj <- unlist(q); k <- pi_ < pj
  ov <- if (any(k)) geos_area(geos_intersection(g[pi_[k]], g[pj[k]])) else 0
  data.frame(n = n, c3_geos = sum(!geos_within(g, gt)), c3_sf = c3_sf, no_pairs = sum(ov > 1e-9),
             n_invalid = sum(!geos_is_valid(g)), n_empty = sum(geos_is_empty(g)), n_capped = sum(cd$nucleus_capped))
}
gen_one <- function(A, method, k, phimode) {
  f <- file.path(OUT, "gen", sprintf("%s_%s_fold%d_%s.rds", A, method, k, phimode)); if (file.exists(f)) return(invisible(NULL))
  o14 <- readRDS(file.path(OUT, "S1.4", "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  te <- load_real(A, o14$roi_test, method); cc <- te$cells; ev <- evset(cc); inp <- inputs_of(te)
  phi <- if (phimode == "sr") PHI_SR[[A]] else PHI(A, method)
  ct <- ct_m2(o14$fits$M2, phi)
  raster <- any(WIN$archetype == A & WIN$roi_id == o14$roi_test)
  gdir <- file.path(OUT, sprintf("gen_in_%s_%s", method, phimode)); dir.create(gdir, showWarnings = FALSE)
  rows <- list()
  for (s in seq_len(NSEED)) {
    sd <- o14$seed + 100000 * s
    o <- E$derive_cell_geometry(inp$cen, inp$tess, ct, random_seed = sd, verbose = FALSE); cd <- o$cell_df
    rows[[s]] <- cbind(data.frame(archetype = A, method = method, fold = k, roi = o14$roi_test, phimode = phimode, phi = phi, seed_idx = s, seed = sd,
                                  trunc_all = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target),
                                  nc95_ev = mean((cd$nucleus_area / cd$territory_area)[ev] >= 0.95), n_ev = sum(ev),
                                  md5 = if (s == 1L) wkb_md5(o) else NA_character_, elapsed = o$info$elapsed_s, stringsAsFactors = FALSE),
                       invariants(o, te$terr))
    if (raster) write_nuc(file.path(gdir, sprintf("%s_%s_fold%d_s%02d_K2.parquet", A, o14$roi_test, k, s)), as_geos_geometry(o$cell_nuclei))
  }
  saveRDS(do.call(rbind, rows), f); invisible(NULL)
}

# ---- cp1: centri = punti di Luca nella finestra (fuori dalle esclusioni), centri CO altrove ---------------------------
excl_um <- function(wv) {
  ann <- jsonlite::fromJSON("results/R2b/annotations/R2b_annotations_Luca.json", simplifyVector = FALSE)$windows
  aj <- ann[[which(vapply(ann, function(z) z$win_id, "") == wv$win_id)]]; m0 <- aj$margin_px; upp <- aj$um_per_px
  if (!length(aj$exclusions)) return(NULL)
  st_union(st_sfc(lapply(aj$exclusions, function(e) {
    p <- do.call(rbind, lapply(e$pts, unlist)); p <- cbind(wv$x0_um + (p[, 1] - m0) * upp, wv$y0_um + (p[, 2] - m0) * upp)
    st_make_valid(st_polygon(list(rbind(p, p[1, ])))) })))
}
cp1_one <- function(win_id) {
  wv <- WIN[WIN$win_id == win_id, ]; A <- wv$archetype; roi <- wv$roi_id
  f <- file.path(OUT, "cp1", sprintf("%s.rds", win_id)); if (file.exists(f)) return(invisible(NULL))
  k <- match(roi, ids_of(A)); o14 <- readRDS(file.path(OUT, "S1.4", "loo", sprintf("%s_consensus_fold%d.rds", A, k)))
  stopifnot(o14$roi_test == roi)
  w <- r3_roi_window(A, roi, rois); te <- load_real(A, roi, "consensus")
  s <- wv$side_um; inwin <- function(x, y) x >= wv$x0_um & x < wv$x0_um + s & y >= wv$y0_um & y < wv$y0_um + s
  ex <- excl_um(wv)
  inex <- function(x, y) if (is.null(ex)) rep(FALSE, length(x)) else lengths(st_intersects(st_as_sf(data.frame(x = x, y = y), coords = c("x", "y")), ex)) > 0
  co <- te$cells[, c("x", "y")]; keep_co <- !inwin(co$x, co$y) | inex(co$x, co$y)
  pts <- read.csv("results/R2b/R2b_points.csv"); pts <- pts[pts$rater == "Luca" & pts$win_id == win_id, ]
  ok <- r3_in_mask(pts$x_um, pts$y_um, w$m, w$px)
  gx <- c(co$x[keep_co], pts$x_um[ok]); gy <- c(co$y[keep_co], pts$y_um[ok]); src <- c(rep("CO", sum(keep_co)), rep("manual", sum(ok)))
  n <- length(gx); cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = gx, y = gy, cell_type = "c")
  tv <- tessellate_voronoi(cen, one_region(w), verbose = FALSE); terr <- E$.dcg_as_polys(tv$cell_territories)
  tess <- list(cell_territories = terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = tv$territory_df$territory_area))
  ct <- ct_m2(o14$fits$M2, PHI(A, "consensus"))
  gdir <- file.path(OUT, "gen_in_cp1"); dir.create(gdir, showWarnings = FALSE)
  rows <- list()
  for (s_ in seq_len(NSEED)) {
    sd <- o14$seed + 100000 * s_
    o <- E$derive_cell_geometry(cen, tess, ct, random_seed = sd, verbose = FALSE); cd <- o$cell_df
    rows[[s_]] <- cbind(data.frame(win_id = win_id, seed_idx = s_, seed = sd, trunc_all = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target)),
                        invariants(o, terr))
    write_nuc(file.path(gdir, sprintf("%s_%s_fold%d_s%02d_K2.parquet", A, roi, k, s_)), as_geos_geometry(o$cell_nuclei))
  }
  saveRDS(list(win_id = win_id, n = n, n_co_kept = sum(keep_co), n_co_removed = sum(!keep_co), n_manual = nrow(pts),
               n_manual_in_mask = sum(ok), src = src, x = gx, y = gy, rows = do.call(rbind, rows)), f)
  invisible(NULL)
}

# ---- regress: cancello G-3 (parte R) -----------------------------------------------------------------------------
regress <- function() {
  out <- list()
  for (A in TIS) for (roi in ids_of(A)) {
    g1 <- as.data.frame(read_parquet(rp(A, roi, "stardist_he", "gen"))); g0 <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_stardist_he_gen.parquet", A, roi))))
    c1 <- as.data.frame(read_parquet(rp(A, roi, "stardist_he", "cells"))); c0 <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_stardist_he_cells.parquet", A, roi))))
    n1 <- as.data.frame(read_parquet(rp(A, roi, "stardist_he", "nuc"))); n0 <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_stardist_he_nuc.parquet", A, roi))))
    r1 <- load_real(A, roi, "stardist_he"); r0 <- load_real(A, roi, "stardist_he", S14)
    out[[length(out) + 1]] <- data.frame(check = "R3+S1.4 real", archetype = A, roi = roi,
      gen = isTRUE(all.equal(g1, g0, tolerance = 0)), cells = isTRUE(all.equal(c1, c0, tolerance = 0)), nuc = isTRUE(all.equal(n1, n0, tolerance = 0)),
      s14_cells = isTRUE(all.equal(r1$cells, r0$cells, tolerance = 0)),
      s14_terr = identical(geos_write_wkb(as_geos_geometry(r1$terr)), geos_write_wkb(as_geos_geometry(r0$terr))))
  }
  for (A in TIS) for (k in 1:5) {
    l1 <- readRDS(file.path(OUT, "S1.4", "loo", sprintf("%s_stardist_he_fold%d.rds", A, k))); l0 <- readRDS(file.path(S14, "loo", sprintf("%s_stardist_he_fold%d.rds", A, k)))
    m1 <- readRDS(file.path(OUT, "gen", sprintf("%s_stardist_he_fold%d_sr.rds", A, k))); m0 <- readRDS(file.path(S14D, "gen", sprintf("%s_stardist_he_fold%d_tissue.rds", A, k)))
    out[[length(out) + 1]] <- data.frame(check = "S1.4 loo M2 + S1.4d gen", archetype = A, roi = l1$roi_test,
      gen = identical(l1$roi_test, l0$roi_test) && identical(l1$seed, l0$seed), cells = identical(unname(l1$fits$M2), unname(l0$fits$M2)),
      nuc = identical(m1$md5[1], m0$md5[1]), s14_cells = isTRUE(all.equal(m1$trunc_all, m0$trunc_all, tolerance = 0)), s14_terr = NA)
  }
  ph <- rbind(data.frame(archetype = TIS, method = "spaceranger_S1.4", t(sapply(TIS, function(A) phi_of(A, "spaceranger", S14)))))
  write.csv(do.call(rbind, out), file.path(OUT, "regress", "G3_R.csv"), row.names = FALSE)
  write.csv(ph, file.path(OUT, "regress", "G3_phi_sr.csv"), row.names = FALSE)
  print(do.call(rbind, out)); print(ph)
}

jobs_tm <- function(methods) do.call(rbind, lapply(TIS, function(A) expand.grid(A = A, roi = ids_of(A), method = methods, stringsAsFactors = FALSE)))
if (STAGE == "r3") par_run(jobs_tm(c("consensus", "stardist_he")), r3_one, "r3")
if (STAGE == "real") par_run(jobs_tm(c("consensus", "stardist_he")), real_one, "real")
if (STAGE == "loo") par_run(expand.grid(A = TIS, method = c("consensus", "stardist_he"), k = 1:5, stringsAsFactors = FALSE), loo_one, "loo")
if (STAGE == "phi") {
  p <- do.call(rbind, lapply(TIS, function(A) do.call(rbind, lapply(c("consensus", "stardist_he"), function(m) data.frame(archetype = A, method = m, t(phi_of(A, m)))))))
  write.csv(p, file.path(OUT, "phi.csv"), row.names = FALSE); print(p)
}
if (STAGE == "gen") {
  j <- rbind(expand.grid(A = TIS, method = "consensus", k = 1:5, phimode = "own", stringsAsFactors = FALSE),
             expand.grid(A = TIS, method = "stardist_he", k = 1:5, phimode = c("own", "sr"), stringsAsFactors = FALSE))
  par_run(j, gen_one, "gen")
}
if (STAGE == "cp1") par_run(data.frame(win_id = WIN$win_id, stringsAsFactors = FALSE), cp1_one, "cp1")
if (STAGE == "regress") regress()
