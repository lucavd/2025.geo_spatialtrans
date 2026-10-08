# tools/S1.3_roi.R — S1.3: banco di prova sui 30 ROI R2 (stadi real, null, cp3).
# Pre-registrazione: results/S1.3/S1.3_preregistration.md (commit e2f32fc). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.3_roi.R <stage> [outdir]
#   stage = real | null | cp3 ; outdir default /mnt/micron/geo_spatialtrans/S1.3 ; stadi rilanciabili (saltano i file fatti)
# Definizioni operative (dal testo pre-registrato, banco di prova ROI):
#  - finestra e maschera da r3_roi_window(); semina su w$reg (componenti) come R3; tassellazione con UNA regione = w$poly;
#  - cellula interna come R3: tile dentro la regione, oppure tile dentro il quadrato del ROI e area/tile > 0.999;
#    il tile e' quello del pacchetto (keep_tiles = TRUE); forma (ecc_T, theta_T) e lati dal tile, come R3;
#  - real: generatori = tabelle *_gen.parquet di R3 (40 = 30 primari + 10 StarDist), confronto con *_cells.parquet;
#    entrambi i motori (C-5 su ciascuno; C-10a = geos vs deldir);
#  - null: G1 CSR e G2 RSA(d*) con codice e seed di tools/R3_run.R; G3 RSA regola di default (seed +300);
#    motore geos; replica 1 anche con deldir (C-10a, K-3); tabelle per cellula delle repliche 1 e 2 (B-5);
#  - cp3: centroidi reali primari, cs = 1 con e senza contenimento (mutante M4); lacune a cs = 1/3, 2/3, 1.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel); library(arrow); library(spatstat.geom); library(spatstat.explore) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
source("tools/R3_voronoi_metrics.R"); source("tools/S1.3_mutants.R"); source("tools/S1.3_variants.R"); source("tools/S1.3_arbiter.R")
ENG <- list(geos = tv_variant(tv_env(), "geos"), deldir = tv_variant(tv_env(), "deldir"))
args <- commandArgs(TRUE); STAGE <- args[1]
OUT <- if (length(args) > 1) args[2] else "/mnt/micron/geo_spatialtrans/S1.3"
R3 <- "/mnt/micron/geo_spatialtrans/R3"
for (d in c("real", "null", "cp3")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
NCORES <- as.integer(Sys.getenv("S13_CORES", "40")); NREP <- 20L; BASE_SEED <- 20261008; TOL <- 1e-9
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   sens = c(NA, NA, NA, "stardist_he", NA, "stardist_he"),
                   dstar = c(2.5, 3.5, 3.5, 4, 7, 2.5), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))

one_region <- function(w) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2),
                               region_polygons = st_sfc(w$poly[[1]]))

#' Tassellazione del pacchetto + metriche nel formato di r3_tessellate (stesse colonne usate da r3_summary).
s13_tessellate <- function(x, y, w, backend = "geos", cs = 0) {
  cen <- data.frame(cell_id = seq_along(x), region_id = 1L, x = x, y = y)
  tv <- ENG[[backend]]$tessellate_voronoi(cen, one_region(w), corner_smoothing = cs, keep_tiles = TRUE, verbose = FALSE)
  stopifnot(tv$info$engine == backend)
  td <- tv$territory_df
  in_frame <- lengths(st_within(tv$tiles, w$frame)) > 0
  interior <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
  mom <- t(vapply(tv$tiles, function(p) { m <- unclass(p)[[1]]; m <- m[-nrow(m), , drop = FALSE]; r3_poly_moments(m[, 1], m[, 2]) }, numeric(6)))
  nsides <- vapply(tv$tiles, function(p) nrow(unclass(p)[[1]]) - 1L, 1L)
  sh <- r3_shape(mom[, "mxx"], mom[, "myy"], mom[, "mxy"])
  list(tv = tv, df = data.frame(idx = seq_along(x), x = x, y = y, area_tile = td$tile_area, area = td$territory_area,
                                interior = interior, nsides = nsides, ecc_T = sh$ecc, theta_T = sh$theta,
                                clipped = td$clipped, n_lost = td$n_pieces_lost, n_gained = td$n_pieces_gained))
}
relmax <- function(a, b) if (!length(a)) 0 else max(abs(a - b) / pmax(abs(b), 1e-300))

# controlli C-1/C-2/C-3 su una tassellazione a regione unica
c123 <- function(tv, w, x, y) {
  A <- w$area_um2; t <- tv$cell_territories
  u <- st_union(t)
  pts <- st_sfc(lapply(seq_along(x), function(i) st_point(c(x[i], y[i]))))
  own <- st_intersects(pts, t)
  c(c1 = abs(A - sum(tv$territory_df$territory_area)) / A,
    c2_overlap = (sum(tv$territory_df$territory_area) - as.numeric(st_area(u))) / A,
    c2_symdiff = sum(as.numeric(st_area(st_sym_difference(u, st_sfc(w$poly[[1]]))))) / A,
    c3_poly = mean(as.character(st_geometry_type(t)) == "POLYGON"), c3_valid = mean(st_is_valid(t)),
    c3_gen_in_own = mean(vapply(seq_along(own), function(i) i %in% own[[i]], TRUE)),
    n_multipart = tv$info$n_multipart, n_repaired = tv$info$n_repaired, n_snapped = tv$info$n_snapped, n_fragments = tv$info$n_fragments,
    n_reassigned = tv$info$n_fragments_reassigned, elapsed_s = tv$info$elapsed_s)
}

# ---- stadio real (C-5, C-10a) ----------------------------------------------------------------------
real_one <- function(A, roi, method) {
  f <- file.path(OUT, "real", sprintf("%s_%s_%s.rds", A, roi, method))
  if (file.exists(f)) return(NULL)
  w <- r3_roi_window(A, roi, rois)
  gen <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", A, roi, method))))
  ref <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_cells.parquet", A, roi, method))))
  stopifnot(nrow(gen) == nrow(ref), identical(ref$x, gen$x))
  li <- r3_local_intensity(gen$x, gen$y, w)
  out <- list(archetype = A, roi_id = roi, method = method, n = nrow(gen))
  for (b in c("geos", "deldir")) {
    s <- s13_tessellate(gen$x, gen$y, w, backend = b)
    d <- s$df; mono <- d$n_lost == 0 & d$n_gained == 0
    sumr <- r3_summary(d, li$lambda_loc); sumref <- r3_summary(ref, li$lambda_loc)
    out[[b]] <- list(c123 = c123(s$tv, w, gen$x, gen$y),
                     c5 = c(interior_identical = identical(d$interior, ref$interior),
                            n_interior_diff = sum(d$interior != ref$interior),
                            rel_area_interior = relmax(d$area[ref$interior], ref$area[ref$interior]),
                            rel_area_mono = relmax(d$area[mono], ref$area[mono]),
                            n_frag_cells = sum(!mono), rel_sum_area = abs(sum(d$area) - sum(ref$area)) / sum(ref$area),
                            rel_summary = relmax(sumr, sumref)),
                     cells = d[, c("idx", "area", "interior", "clipped", "n_lost", "n_gained", "nsides", "ecc_T")], tiles = s$tv$tiles,
                     summary = sumr)
  }
  inside <- lengths(st_within(out$geos$tiles, w$frame)) > 0
  out$arb <- arbitrate(out$geos$tiles, out$deldir$tiles, gen$x, gen$y, inside)
  out$geos$tiles <- NULL; out$deldir$tiles <- NULL
  g <- out$geos$cells; dl <- out$deldir$cells
  out$c10a <- c(rel_area = relmax(g$area, dl$area), interior_identical = identical(g$interior, dl$interior),
                frag_identical = identical(g$n_lost, dl$n_lost) && identical(g$n_gained, dl$n_gained),
                nsides_interior_identical = identical(g$nsides[g$interior], dl$nsides[dl$interior]),
                t_geos = out$geos$c123[["elapsed_s"]], t_deldir = out$deldir$c123[["elapsed_s"]])
  out$deldir$cells <- NULL
  saveRDS(out, f); invisible(NULL)
}
if (STAGE == "real") {
  jobs <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) {
    a <- ARCH[ARCH$archetype == rois$archetype[i], ]
    data.frame(archetype = rois$archetype[i], roi_id = rois$roi_id[i], method = c(a$primary, if (!is.na(a$sens)) a$sens), stringsAsFactors = FALSE)
  }))
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(real_one(j$archetype, j$roi_id, j$method), error = function(e) conditionMessage(e))
    cat(sprintf("[real] %s %s %s %.0fs %s\n", j$archetype, j$roi_id, j$method, proc.time()[["elapsed"]] - t0, if (is.character(r)) r else "ok")); r
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  errs <- vapply(res, is.character, TRUE); cat("real: errori", sum(errs), "su", length(res), "\n"); if (any(errs)) quit(status = 1)
}

# ---- stadio null (C-6, check B, C-10a) -------------------------------------------------------------
gen_points <- function(i, model, rp, w, n_obs, a) {
  if (model == "CSR") {
    seed <- BASE_SEED + 1000 * i + 100 + rp
    set.seed(seed); idx <- which(w$m, arr.ind = TRUE); s <- sample.int(nrow(idx), n_obs, replace = TRUE)
    return(list(x = (idx[s, 2] - 1 + runif(n_obs)) * w$px, y = (idx[s, 1] - 1 + runif(n_obs)) * w$px, n_failed = 0L, seed = seed))
  }
  comp <- data.frame(cluster_id = unique(as.character(w$reg$region_df$cluster_id)), cell_type = "c", fraction = 1)
  if (model == "RSA") {
    seed <- BASE_SEED + 1000 * i + 200 + rp
    ct <- data.frame(cell_type = "c", density = n_obs / w$area_um2 * 1e6, min_dist_um = a$dstar)
  } else {                                   # G3: regola di default d = 2/3 eq_r (nessun min_dist_um)
    seed <- BASE_SEED + 1000 * i + 300 + rp
    ct <- data.frame(cell_type = "c", density = n_obs / w$area_um2 * 1e6)
  }
  o <- suppressWarnings(seed_centroids(w$reg, ct, comp, random_seed = seed, verbose = FALSE))
  list(x = o$centroids$x, y = o$centroids$y, n_failed = o$info$n_failed, seed = seed)
}
null_one <- function(i, model, rp) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]
  f <- file.path(OUT, "null", sprintf("%s_%s_%s_%02d.rds", A, roi, model, rp))
  if (file.exists(f)) return(invisible(NULL))
  w <- r3_roi_window(A, roi, rois)
  n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  p <- gen_points(i, model, rp, w, n_obs, a)
  s <- s13_tessellate(p$x, p$y, w, backend = "geos")
  li <- r3_local_intensity(p$x, p$y, w)
  sm <- r3_summary(s$df, li$lambda_loc)
  out <- list(archetype = A, roi_id = roi, model = model, rep = rp, seed = p$seed, n_failed = p$n_failed, summary = sm,
              nsides = tabulate(pmin(s$df$nsides[s$df$interior], 15), 15), c123 = c123(s$tv, w, p$x, p$y))
  if (rp == 1) {
    s2 <- s13_tessellate(p$x, p$y, w, backend = "deldir")
    out$c10a <- c(rel_area = relmax(s$df$area, s2$df$area), interior_identical = identical(s$df$interior, s2$df$interior),
                  frag_identical = identical(s$df$n_lost, s2$df$n_lost) && identical(s$df$n_gained, s2$df$n_gained),
                  nsides_interior_identical = identical(s$df$nsides[s$df$interior], s2$df$nsides[s2$df$interior]),
                  t_geos = s$tv$info$elapsed_s, t_deldir = s2$tv$info$elapsed_s)
    out$c123_deldir <- c123(s2$tv, w, p$x, p$y)
    out$arb <- arbitrate(s$tv$tiles, s2$tv$tiles, p$x, p$y, lengths(st_within(s$tv$tiles, w$frame)) > 0)
  }
  if (rp <= 2) out$cells <- s$df[, c("x", "y", "area", "interior", "nsides", "ecc_T", "theta_T")]   # B-5: rumore = replica 1 vs 2
  saveRDS(out, f); invisible(NULL)
}
if (STAGE == "null") {
  jobs <- expand.grid(i = seq_len(nrow(rois)), model = c("CSR", "RSA", "RSArule"), rp = seq_len(NREP), stringsAsFactors = FALSE)
  jobs <- jobs[order(jobs$rp != 1, -jobs$i), ]          # repliche 1 (con deldir) per prime, ROI densi in testa
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    r <- tryCatch(null_one(j$i, j$model, j$rp), error = function(e) conditionMessage(e))
    if (is.character(r)) cat(sprintf("[null] ERRORE %d %s %d: %s\n", j$i, j$model, j$rp, r)); r
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  errs <- vapply(res, is.character, TRUE); cat("null: errori", sum(errs), "su", length(res), "\n"); if (any(errs)) quit(status = 1)
}

# ---- stadio cp3 (smussatura sui ROI reali) ----------------------------------------------------------
cp3_one <- function(i) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]
  f <- file.path(OUT, "cp3", sprintf("%s_%s.rds", A, roi)); if (file.exists(f)) return(invisible(NULL))
  w <- r3_roi_window(A, roi, rois)
  gen <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", A, roi, a$primary))))
  cen <- data.frame(cell_id = seq_len(nrow(gen)), region_id = 1L, x = gen$x, y = gen$y); reg <- one_region(w)
  t0 <- tessellate_voronoi(cen, reg, verbose = FALSE)
  ta <- t0$territory_df$territory_area
  convex <- abs(as.numeric(st_area(st_convex_hull(t0$cell_territories))) - ta) <= 1e-9 * ta
  overl <- function(tv) {
    pr <- st_relate(tv$cell_territories, tv$cell_territories, pattern = "2********")
    p <- do.call(rbind, lapply(seq_along(pr), function(k) { j <- pr[[k]]; j <- j[j > k]; if (length(j)) cbind(k, j) }))
    if (is.null(p)) return(list(n = 0L, n_bad = 0L))
    ar <- vapply(seq_len(nrow(p)), function(r) as.numeric(st_area(st_sfc(st_intersection(tv$cell_territories[[p[r, 1]]], tv$cell_territories[[p[r, 2]]])))), 0)
    p <- p[ar > 1e-6, , drop = FALSE]
    inv <- (t0$territory_df$clipped & !convex)
    list(n = nrow(p), n_bad = sum(!(inv[p[, 1]] | inv[p[, 2]])))
  }
  E4 <- tv_mutate(tv_env(), "M4")
  s_ok <- tessellate_voronoi(cen, reg, corner_smoothing = 1, verbose = FALSE)
  s_m4 <- E4$tessellate_voronoi(cen, reg, corner_smoothing = 1, verbose = FALSE)
  gaps <- vapply(c(1/3, 2/3, 1), function(cs) sum(tessellate_voronoi(cen, reg, corner_smoothing = cs, verbose = FALSE)$region_check$gap_area) / w$area_um2, 0)
  out <- list(archetype = A, roi_id = roi, n = nrow(gen), n_clipped_nonconvex = sum(t0$territory_df$clipped & !convex),
              ov_contained = overl(s_ok), ov_m4 = overl(s_m4), gap_frac = setNames(gaps, c("cs1/3", "cs2/3", "cs1")))
  saveRDS(out, f); invisible(NULL)
}
if (STAGE == "cp3") {
  res <- mclapply(seq_len(nrow(rois)), function(i) {
    r <- tryCatch(cp3_one(i), error = function(e) conditionMessage(e)); cat(sprintf("[cp3] %d %s\n", i, if (is.character(r)) r else "ok")); r
  }, mc.cores = min(NCORES, 30L), mc.preschedule = FALSE)
  errs <- vapply(res, is.character, TRUE); cat("cp3: errori", sum(errs), "su", length(res), "\n"); if (any(errs)) quit(status = 1)
}
