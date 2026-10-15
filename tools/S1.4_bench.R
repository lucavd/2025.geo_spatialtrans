# tools/S1.4_bench.R — S1.4: banco di derive_cell_geometry() sui 30 ROI R2 (stadi real, loo, rec, null).
# Pre-registrazione: results/S1.4/S1.4_preregistration.md (commit c72c01a); definizioni operative in
# results/S1.4/S1.4_dev_notes.md (fissate prima degli stadi ufficiali). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S14_MUTANT=none|M1] Rscript --vanilla tools/S1.4_bench.R <stage> [outdir]
#   stage = real | loo | rec | null ; outdir default /mnt/micron/geo_spatialtrans/S1.4 ; stadi rilanciabili.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
source("tools/R3_voronoi_metrics.R"); source("tools/S1.4_mutants.R")
MUT <- Sys.getenv("S14_MUTANT", "none")
E <- dcg_mutate(dcg_env(), MUT)
args <- commandArgs(TRUE); STAGE <- args[1]
OUT <- if (length(args) > 1) args[2] else "/mnt/micron/geo_spatialtrans/S1.4"
for (d in c("real", "loo", "rec", "null")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
R3 <- "/mnt/micron/geo_spatialtrans/R3"; S13 <- "/mnt/micron/geo_spatialtrans/S1.3"
NCORES <- as.integer(Sys.getenv("S14_CORES", "40")); BASE13 <- 20261008; BASE14 <- 20261015
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   sens = c(NA, NA, NA, "stardist_he", NA, "stardist_he"),
                   dstar = c(2.5, 3.5, 3.5, 4, 7, 2.5), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
one_region <- function(w) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2),
                               region_polygons = st_sfc(w$poly[[1]]))
rp <- function(A, roi, method, k) file.path(R3, "real", sprintf("%s_%s_%s_%s.parquet", A, roi, method, k))
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b)))
  max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
spear <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
c3_sf <- function(nuc, terr) { w <- st_within(nuc, terr); sum(vapply(seq_along(w), function(i) i %in% w[[i]], TRUE)) }

# ---- stadio real: territori del pacchetto sui generatori reali + nucleo reale per cellula ------------------------------
real_one <- function(A, roi, method) {
  f <- file.path(OUT, "real", sprintf("%s_%s_%s.rds", A, roi, method)); if (file.exists(f)) return(invisible(NULL))
  w <- r3_roi_window(A, roi, rois)
  gen <- as.data.frame(read_parquet(rp(A, roi, method, "gen"))); nucp <- as.data.frame(read_parquet(rp(A, roi, method, "nuc")))
  stopifnot(identical(nucp$idx, gen$idx)); n <- nrow(gen)
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = gen$x, y = gen$y, cell_type = "c")
  tv <- tessellate_voronoi(cen, one_region(w), keep_tiles = TRUE, verbose = FALSE)
  td <- tv$territory_df; terr <- E$.dcg_as_polys(tv$cell_territories)
  in_frame <- lengths(st_within(tv$tiles, w$frame)) > 0
  interior <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
  ref13 <- readRDS(file.path(S13, "real", sprintf("%s_%s_%s.rds", A, roi, method)))$geos$cells
  c0 <- c(rel_area_vs_S13 = max(abs(td$territory_area - ref13$area) / ref13$area), interior_identical_S13 = identical(interior, ref13$interior))
  a_own <- gen$area_nuc
  if (method == "spaceranger") {             # BL-060: area e forma dal nucleo StarDist il cui centroide cade nel territorio
    sg <- as.data.frame(read_parquet(rp(A, roi, "stardist_he", "gen"))); sn <- as.data.frame(read_parquet(rp(A, roi, "stardist_he", "nuc")))
    hit <- st_intersects(st_as_sf(sg, coords = c("x", "y")), terr)
    h1 <- vapply(hit, function(z) if (length(z) == 1) z else NA_integer_, 1L)
    cnt <- tabulate(h1[!is.na(h1)], n)
    matched <- cnt == 1
    j <- match(seq_len(n), h1); j[!matched] <- NA
    a_n <- sg$area_nuc[j]; ecc_n <- sn$ecc_N[j]; th_n <- sn$theta_N[j]
    match_info <- c(n_sd = nrow(sg), n_sd_in = sum(!is.na(h1)), n_terr_0 = sum(cnt == 0), n_terr_1 = sum(cnt == 1), n_terr_2p = sum(cnt >= 2))
  } else {
    matched <- rep(TRUE, n); a_n <- a_own; ecc_n <- nucp$ecc_N; th_n <- nucp$theta_N; match_info <- NULL
  }
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
load_real <- function(A, roi, method) readRDS(file.path(OUT, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
tess_of <- function(r) list(cell_territories = r$terr, territory_df = data.frame(cell_id = seq_len(r$n), territory_area = r$cells$a_t))
cen_of <- function(r) data.frame(cell_id = seq_len(r$n), region_id = 1L, x = r$cells$x, y = r$cells$y, cell_type = "c")
ctab <- function(fp, a_ref, ...) { d <- data.frame(cell_type = "c", density = 1e6 / a_ref, nucleus_to_eq_ratio = fp[["rho"]],
  nucleus_area_slope = fp[["beta"]], nucleus_area_sdlog = fp[["sigma"]]); a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]; d }
fit_set <- function(rs, model = "full", permute_seed = NULL) {
  a_ref <- sum(vapply(rs, function(r) sum(r$cells$a_t), 0)) / sum(vapply(rs, function(r) r$n, 0))
  cc <- do.call(rbind, lapply(rs, function(r) r$cells[r$cells$interior & r$cells$matched & is.finite(r$cells$a_n) & r$cells$a_n > 0, ]))
  at <- cc$a_t; if (!is.null(permute_seed)) { set.seed(permute_seed); at <- sample(at) }
  list(fit = E$fit_nucleus_size_model(cc$a_n, at, a_ref, model), a_ref = a_ref, cc = cc)
}

# ---- stadio loo: validazione leave-one-ROI-out (B-1..B-5, CP-1) ------------------------------------------------------
MODELS <- c(M0 = "design", M1 = "beta1", M2 = "full", M3 = "beta0")
loo_one <- function(A, method, k) {
  f <- file.path(OUT, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)); if (file.exists(f)) return(invisible(NULL))
  ids <- rois$roi_id[rois$archetype == A]; ai <- match(A, ARCH$archetype)
  rs <- lapply(ids, function(r) load_real(A, r, method)); te <- rs[[k]]; tr <- rs[-k]
  ev <- te$cells$interior & te$cells$matched & is.finite(te$cells$a_n) & te$cells$a_n > 0
  real_nc <- te$cells$a_n[ev] / te$cells$a_t[ev]; real_an <- te$cells$a_n[ev]; real_sp <- spear(real_an, te$cells$a_t[ev])
  yard <- do.call(rbind, lapply(seq_along(rs)[-k], function(j) { cj <- rs[[j]]$cells
    e2 <- cj$interior & cj$matched & is.finite(cj$a_n) & cj$a_n > 0
    data.frame(j = ids[j], ks_nc = ks_d(real_nc, cj$a_n[e2] / cj$a_t[e2]), ks_an = ks_d(real_an, cj$a_n[e2]),
               d_sp = abs(real_sp - spear(cj$a_n[e2], cj$a_t[e2]))) }))
  seed <- BASE14 + 1000 * ai + 10 * k
  cen <- cen_of(te); tess <- tess_of(te)
  run <- function(ct, sd, ...) E$derive_cell_geometry(cen, tess, ct, random_seed = sd, verbose = FALSE, ...)
  ev_sim <- function(o) { a <- o$cell_df$nucleus_area[ev]
    c(ks_nc = ks_d(a / te$cells$a_t[ev], real_nc), ks_an = ks_d(a, real_an), d_sp = abs(spear(a, te$cells$a_t[ev]) - real_sp),
      clamp_all = mean(o$cell_df$nucleus_clamped), clamp_ev = mean(o$cell_df$nucleus_clamped[ev]),
      c3_viol = o$info$n_not_within, min_area = min(o$cell_df$nucleus_area)) }
  fits <- list(); sims <- list()
  for (m in names(MODELS)) {
    fs <- fit_set(tr, MODELS[[m]]); fits[[m]] <- c(fs$fit, a_ref = fs$a_ref)
    o <- run(ctab(fs$fit, fs$a_ref), seed)
    sims[[m]] <- ev_sim(o)
    if (m == "M2") { o2 <- o; fs2 <- fs }
  }
  # C3 indipendente (sf) e area dei poligoni sul modello M2
  c3 <- c(n = te$n, within_sf = c3_sf(o2$cell_nuclei, te$terr),
          rel_area = max(abs(as.numeric(st_area(o2$cell_nuclei)) - o2$cell_df$nucleus_area) / o2$cell_df$nucleus_area))
  # sensibilita' al seme (M2, 10 semi)
  seeds <- t(vapply(1:10, function(s) ev_sim(run(ctab(fs2$fit, fs2$a_ref), seed + 100000 * s))[c("ks_nc", "ks_an", "d_sp", "clamp_all")], numeric(4)))
  # B-5 descrittive: ellisse (aspetto dalla mediana di ecc_N del training), orientazione territorio / casuale; centro = generatore
  e_med <- median(fs2$cc$ecc_n, na.rm = TRUE); asp <- 1 / sqrt(1 - e_med^2)
  dth <- function(a, b) r3_axis_diff(a, b)
  oe_t <- run(ctab(fs2$fit, fs2$a_ref, nucleus_aspect = asp, nucleus_orient = "territory"), seed)
  oe_r <- run(ctab(fs2$fit, fs2$a_ref, nucleus_aspect = asp, nucleus_orient = "random"), seed)
  og <- run(ctab(fs2$fit, fs2$a_ref), seed, nucleus_center = "generator")
  eqr <- sqrt(te$cells$a_t / pi)
  b5 <- c(e_med_train = e_med, aspect = asp,
          real_med_dtheta = median(dth(te$cells$theta_n[ev], te$cells$theta_t[ev]), na.rm = TRUE),
          sim_terr_med_dtheta = median(dth(oe_t$cell_df$nucleus_orient[ev], te$cells$theta_t[ev])),
          sim_rand_med_dtheta = median(dth(oe_r$cell_df$nucleus_orient[ev], te$cells$theta_t[ev])),
          clamp_ell_terr = mean(oe_t$cell_df$nucleus_clamped), clamp_ell_rand = mean(oe_r$cell_df$nucleus_clamped),
          clamp_gen = mean(og$cell_df$nucleus_clamped),
          ks_an_ell_terr = ks_d(oe_t$cell_df$nucleus_area[ev], real_an), ks_an_gen = ks_d(og$cell_df$nucleus_area[ev], real_an),
          poi_offset_med = median(sqrt((te$cells$poi_x - te$cells$x)^2 + (te$cells$poi_y - te$cells$y)^2)[ev] / eqr[ev]),
          poi_offset_q90 = unname(quantile(sqrt((te$cells$poi_x - te$cells$x)^2 + (te$cells$poi_y - te$cells$y)^2)[ev] / eqr[ev], 0.9)))
  saveRDS(list(archetype = A, method = method, fold = k, roi_test = ids[k], seed = seed, n_ev = sum(ev), n_cells = te$n,
               fits = fits, sims = sims, yard = yard, real_sp = real_sp, c3 = c3, seeds = seeds, b5 = b5,
               qs = list(real_nc = quantile(real_nc, c(.05, .25, .5, .75, .95)), sim_nc = quantile(o2$cell_df$nucleus_area[ev] / te$cells$a_t[ev], c(.05, .25, .5, .75, .95)))), f)
  invisible(NULL)
}

# ---- stadio rec: recupero dei parametri (C-rec) ------------------------------------------------------------------------
REC <- expand.grid(beta = c(0, 0.5, 1), sigma = c(0.2, 0.4))
rec_one <- function(A) {
  f <- file.path(OUT, "rec", sprintf("%s.rds", A)); if (file.exists(f)) return(invisible(NULL))
  a <- ARCH[ARCH$archetype == A, ]; ids <- rois$roi_id[rois$archetype == A]; ai <- match(A, ARCH$archetype)
  rs <- lapply(ids, function(r) load_real(A, r, a$primary))
  a_ref <- sum(vapply(rs, function(r) sum(r$cells$a_t), 0)) / sum(vapply(rs, function(r) r$n, 0))
  out <- do.call(rbind, lapply(seq_len(nrow(REC)), function(g) {
    fp <- c(rho = 0.3, beta = REC$beta[g], sigma = REC$sigma[g])
    sim <- do.call(rbind, lapply(seq_along(rs), function(i) { r <- rs[[i]]
      o <- E$derive_cell_geometry(cen_of(r), tess_of(r), ctab(fp, a_ref), random_seed = BASE14 + 50000 + 1000 * ai + 10 * g + i, verbose = FALSE)
      data.frame(a_n = o$cell_df$nucleus_area, a_t = r$cells$a_t, interior = r$cells$interior, clamped = o$cell_df$nucleus_clamped) }))
    s <- sim[sim$interior, ]
    ft <- E$fit_nucleus_size_model(s$a_n, s$a_t, a_ref, "full")
    data.frame(archetype = A, beta = fp[["beta"]], sigma = fp[["sigma"]], rho = 0.3, n = nrow(s),
               rho_hat = ft[["rho"]], beta_hat = ft[["beta"]], sigma_hat = ft[["sigma"]], clamp = mean(s$clamped))
  }))
  saveRDS(out, f); invisible(NULL)
}

# ---- stadio null: C3, C5, C-area, C7 su tassellazioni nulle (CSR e RSA(d*), replica 1, semi di S1.3) ----------------
null_one <- function(i, model) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]
  f <- file.path(OUT, "null", sprintf("%s_%s_%s.rds", A, roi, model)); if (file.exists(f)) return(invisible(NULL))
  w <- r3_roi_window(A, roi, rois)
  n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  if (model == "CSR") {
    seed <- BASE13 + 1000 * i + 100 + 1
    set.seed(seed); idx <- which(w$m, arr.ind = TRUE); s <- sample.int(nrow(idx), n_obs, replace = TRUE)
    x <- (idx[s, 2] - 1 + runif(n_obs)) * w$px; y <- (idx[s, 1] - 1 + runif(n_obs)) * w$px
  } else {
    seed <- BASE13 + 1000 * i + 200 + 1
    comp <- data.frame(cluster_id = unique(as.character(w$reg$region_df$cluster_id)), cell_type = "c", fraction = 1)
    o <- suppressWarnings(seed_centroids(w$reg, data.frame(cell_type = "c", density = n_obs / w$area_um2 * 1e6, min_dist_um = a$dstar),
                                         comp, random_seed = seed, verbose = FALSE))
    x <- o$centroids$x; y <- o$centroids$y
  }
  n <- length(x); cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = x, y = y, cell_type = "c")
  tv <- tessellate_voronoi(cen, one_region(w), verbose = FALSE)
  dens <- n / w$area_um2 * 1e6
  rs <- lapply(rois$roi_id[rois$archetype == A], function(r) load_real(A, r, a$primary)); fs <- fit_set(rs, "full")
  res <- list()
  for (v in c("circle", "ellipse")) {
    ct <- ctab(fs$fit, 1e6 / dens, nucleus_aspect = if (v == "circle") 1 else 1.6, nucleus_orient = if (v == "circle") "random" else "territory")
    o <- tryCatch(E$derive_cell_geometry(cen, tv, ct, random_seed = seed, verbose = FALSE, strict = FALSE), error = function(e) conditionMessage(e))
    if (is.character(o)) { res[[v]] <- c(error = 1); next }
    cd <- o$cell_df; un <- !cd$nucleus_clamped
    ga <- as.numeric(st_area(o$cell_nuclei))
    res[[v]] <- c(error = 0, n = n, within_sf = c3_sf(o$cell_nuclei, tv$cell_territories), n_area_le0 = sum(cd$nucleus_area <= 0),
                  rel_area_poly = max(abs(ga - cd$nucleus_area) / cd$nucleus_area),
                  rel_area_unclamped = if (any(un)) max(abs(cd$nucleus_area[un] - cd$nucleus_area_target[un]) / cd$nucleus_area_target[un]) else 0,
                  min_cyto = min(cd$cytoplasm_area), clamp = mean(cd$nucleus_clamped), elapsed = o$info$elapsed_s)
  }
  eq <- sqrt(tv$territory_df$territory_area / pi); expct <- sqrt(1e6 / (pi * dens))
  res$c7 <- c(p90_eq_r = unname(quantile(eq, 0.9)), expected = expct, ratio = unname(quantile(eq, 0.9)) / expct)
  res$meta <- c(seed = seed, n = n, density = dens)
  saveRDS(res, f); invisible(NULL)
}

# ---- esecuzione ----------------------------------------------------------------------------------------------------
par_run <- function(jobs, fun, label, cores = NCORES) {
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(do.call(fun, as.list(j)), error = function(e) conditionMessage(e))
    cat(sprintf("[%s] %s %.0fs %s\n", label, paste(unlist(j), collapse = " "), proc.time()[["elapsed"]] - t0, if (is.character(r)) r else "ok")); r
  }, mc.cores = cores, mc.preschedule = FALSE)
  errs <- vapply(res, is.character, TRUE); cat(label, ": errori", sum(errs), "su", length(res), "\n"); if (any(errs)) quit(status = 1)
}
if (STAGE == "real") {
  jobs <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) { a <- ARCH[ARCH$archetype == rois$archetype[i], ]
    data.frame(A = rois$archetype[i], roi = rois$roi_id[i], method = c(a$primary, if (!is.na(a$sens)) a$sens), stringsAsFactors = FALSE) }))
  par_run(jobs, real_one, "real")
}
if (STAGE == "loo") {
  jobs <- do.call(rbind, lapply(seq_len(nrow(ARCH)), function(i) { a <- ARCH[i, ]
    expand.grid(A = a$archetype, method = c(a$primary, if (!is.na(a$sens)) a$sens), k = 1:5, stringsAsFactors = FALSE) }))
  par_run(jobs, loo_one, "loo")
}
if (STAGE == "rec") par_run(data.frame(A = ARCH$archetype, stringsAsFactors = FALSE), rec_one, "rec")
if (STAGE == "null") par_run(expand.grid(i = seq_len(nrow(rois)), model = c("CSR", "RSA"), stringsAsFactors = FALSE), null_one, "null")
