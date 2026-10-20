# tools/S1.4d_bench.R — S1.4d: banco del modello "clip" di derive_cell_geometry() (cancelli G-3, G-4; P-2, P-3, P-5).
# Pre-registrazione: docs/steps/S1.4d.md (commit 5118bce; ok di Luca eae40f6). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 Rscript --vanilla tools/S1.4d_bench.R <stage> [outdir]
#   stage = gen | legacy | null | perf | repro ; outdir default /mnt/micron/geo_spatialtrans/S1.4d ; rilanciabili
#   (le unita' gia' scritte si saltano). Worker: file <outdir>/WORKERS (letto a ogni blocco; default 6).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
source("tools/R3_voronoi_metrics.R"); source("tools/S1.4d_mutants.R")
E <- dcg_env_d()                                          # pacchetto corrente (modello "clip" e "scale")
args <- commandArgs(TRUE); STAGE <- args[1]
OUT <- if (length(args) > 1 && !is.na(args[2])) args[2] else "/mnt/micron/geo_spatialtrans/S1.4d"
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"; S14B <- "/mnt/micron/geo_spatialtrans/S1.4b"; S14C <- "/mnt/micron/geo_spatialtrans/S1.4c"
R3 <- "/mnt/micron/geo_spatialtrans/R3"
for (d in c("gen", "gen_in_tissue", "gen_in_sd", "legacy", "null", "perf", "repro")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
NSEED <- 10L; GRID <- 1e-6; BASE13 <- 20261008
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   sens = c(NA, NA, NA, "stardist_he", NA, "stardist_he"),
                   dstar = c(2.5, 3.5, 3.5, 4, 7, 2.5), stringsAsFactors = FALSE)
# phi del tessuto: q99 del rapporto nucleo/territorio reale, B-121...B-126 (tre decimali come nella riga), tetto 1
PHI <- c(A1 = 0.757, A2 = 1, A3 = 0.947, A4 = 0.816, A5 = 0.474, A6 = 0.180)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
WIN_ROI <- unique(read.csv("results/R2b/R2b_windows.csv")[, c("archetype", "roi_id")])   # ROI con finestre di R2c: raster
load_real <- function(A, roi, method) readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
evset <- function(cc) cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b)))
  max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
spear <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
ct_m2 <- function(f2, phi = NULL) { ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]],
  nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]]); if (!is.null(phi)) ct$nucleus_nc_max <- phi; ct }
inputs_of <- function(te) { n <- te$n; cc <- te$cells
  list(cen = data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c"),
       tess = list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t))) }
wkb_md5 <- function(o) { f <- tempfile(); saveRDS(list(o$cell_df, geos_write_wkb(as_geos_geometry(o$cell_nuclei))), f, compress = FALSE)
  m <- unname(tools::md5sum(f)); unlink(f); m }
write_nuc <- function(f, g) { if (file.exists(f)) return(invisible(NULL))
  write_parquet(data.frame(idx = seq_along(g), wkb = I(geos_write_wkb(g))), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f) }

#' Invarianti di G-3 su un'uscita del modello "clip": C3 (geos; sf via WKB, indipendente), sovrapposizioni, area, validita'.
invariants <- function(o, terr) {
  g <- as_geos_geometry(o$cell_nuclei); gt <- as_geos_geometry(E$.dcg_as_polys(terr)); cd <- o$cell_df; n <- nrow(cd)
  w2 <- st_within(o$cell_nuclei, terr); c3_sf <- sum(!vapply(seq_len(n), function(i) i %in% w2[[i]], TRUE))
  tr <- geos_strtree(g); q <- geos_strtree_query(tr, g); pi_ <- rep(seq_along(q), lengths(q)); pj <- unlist(q); k <- pi_ < pj
  ov <- if (any(k)) geos_area(geos_intersection(g[pi_[k]], g[pj[k]])) else 0
  rel <- abs(cd$nucleus_area - cd$nucleus_area_cap) / cd$nucleus_area_cap
  data.frame(n = n, c3_geos = sum(!geos_within(g, gt)), c3_sf = c3_sf, no_pairs = sum(ov > 1e-9), no_max = max(c(0, ov)),
             area_rel_max_uncapped = max(c(0, rel[!cd$nucleus_capped])), area_rel_max_capped = max(c(0, rel[cd$nucleus_capped])),
             n_invalid = sum(!geos_is_valid(g)), n_empty = sum(geos_is_empty(g)), n_multipart = sum(geos_num_geometries(g) > 1L),
             n_clipped = sum(cd$nucleus_clipped), n_capped = sum(cd$nucleus_capped), max_iter = max(cd$nucleus_n_iter),
             mean_iter_clipped = if (any(cd$nucleus_clipped)) mean(cd$nucleus_n_iter[cd$nucleus_clipped]) else 0)
}

#' Esecuzione a blocchi con numero di worker letto da OUT/WORKERS prima di ogni blocco (limite termico condiviso con MMM).
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

one_region <- function(w) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2),
                               region_polygons = st_sfc(w$poly[[1]]))

# ---- gen: modello "clip" sui tassellamenti reali; phi del tessuto ("tissue") o del fold di S1.4c ("fold") -------------
gen_one <- function(A, method, k, phimode) {
  f <- file.path(OUT, "gen", sprintf("%s_%s_fold%d_%s.rds", A, method, k, phimode)); if (file.exists(f)) return(invisible(NULL))
  primary <- method == ARCH$primary[ARCH$archetype == A]
  o14 <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  te <- load_real(A, o14$roi_test, method); cc <- te$cells; ev <- evset(cc); inp <- inputs_of(te)
  pc <- if (primary) readRDS(file.path(S14C, "gen", sprintf("%s_%s_fold%d.rds", A, method, k))) else NULL
  phi <- if (phimode == "fold") pc$phi else PHI[[A]]
  ct <- ct_m2(o14$fits$M2, phi)
  raster <- phimode == "tissue" && any(WIN_ROI$archetype == A & WIN_ROI$roi_id == o14$roi_test)
  gdir <- if (primary) "gen_in_tissue" else "gen_in_sd"
  rows <- list()
  for (s in seq_len(NSEED)) {
    sd <- o14$seed + 100000 * s
    o <- E$derive_cell_geometry(inp$cen, inp$tess, ct, random_seed = sd, verbose = FALSE); cd <- o$cell_df
    r <- cbind(data.frame(archetype = A, method = method, fold = k, roi = o14$roi_test, phimode = phimode, phi = phi, seed_idx = s, seed = sd,
                          trunc_all = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target),
                          nc95_ev = mean((cd$nucleus_area / cd$territory_area)[ev] >= 0.95), n_ev = sum(ev),
                          md5 = if (s == 1L) wkb_md5(o) else NA_character_, elapsed = o$info$elapsed_s, stringsAsFactors = FALSE),
               invariants(o, te$terr))
    if (phimode == "fold") {                                 # G-4(a): contro i nuclei K-2 del prototipo S1.4c, stesso seme
      pf <- file.path(S14C, "gen_in", sprintf("%s_%s_fold%d_s%02d_K2.parquet", A, o14$roi_test, k, s))
      gp <- geos_read_wkb(unclass(read_parquet(pf)$wkb)); g <- as_geos_geometry(o$cell_nuclei); cap <- cd$nucleus_area_cap
      dA <- abs(cd$nucleus_area - geos_area(gp)); sym <- geos_area(geos_sym_difference(g, gp))
      bound <- dA + sqrt(2) * GRID * pmax(geos_length(g), geos_length(gp))
      ncp <- pc$checks$n_capped[pc$checks$variant == "K-2" & pc$checks$seed_idx == s]
      r <- cbind(r, data.frame(eq_n_proto = length(gp), eq_dA_rel_max = max(dA / cap), eq_sym_rel_max = max(sym / cap),
                               eq_sym_rel_gt1e6 = sum(sym / cap > 1e-6), eq_beyond_bound = sum(sym > bound + 1e-12),
                               eq_ncap_proto = as.integer(ncp), eq_ncap_same = isTRUE(as.integer(ncp) == sum(cd$nucleus_capped))))
    }
    if (raster) write_nuc(file.path(OUT, gdir, sprintf("%s_%s_fold%d_s%02d_K2.parquet", A, o14$roi_test, k, s)), as_geos_geometry(o$cell_nuclei))
    rows[[s]] <- r
  }
  saveRDS(do.call(rbind, rows), f); invisible(NULL)
}

# ---- legacy: modo "scale" = pacchetto al commit db4afea (G-4b) e = righe V-c di S1.4b (C-0 di S1.4c) -----------------
legacy_one <- function(A, k) {
  method <- ARCH$primary[ARCH$archetype == A]
  f <- file.path(OUT, "legacy", sprintf("%s_%s_fold%d.rds", A, method, k)); if (file.exists(f)) return(invisible(NULL))
  o14 <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  ob14 <- readRDS(file.path(S14B, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  te <- load_real(A, o14$roi_test, method); cc <- te$cells; ev <- evset(cc); inp <- inputs_of(te); ct <- ct_m2(o14$fits$M2)
  real_nc <- cc$a_n[ev] / cc$a_t[ev]; real_an <- cc$a_n[ev]; real_sp <- spear(real_an, cc$a_t[ev])
  metr <- function(a_fin, trunc, loss, rdist) c(ks_nc = ks_d(a_fin[ev] / cc$a_t[ev], real_nc), ks_an = ks_d(a_fin[ev], real_an),
    d_sp = abs(spear(a_fin[ev], cc$a_t[ev]) - real_sp), trunc_all = mean(trunc), trunc_ev = mean(trunc[ev]),
    loss_all = mean(loss), b7 = mean(sqrt(a_fin[ev] / pi) > rdist[ev]))
  rows <- list()
  for (s in 0:NSEED) {
    sd <- o14$seed + 100000 * s
    vn <- E$derive_cell_geometry(inp$cen, inp$tess, ct, random_seed = sd, nucleus_model = "scale", verbose = FALSE)
    vo <- EOLD$derive_cell_geometry(inp$cen, inp$tess, ct, random_seed = sd, verbose = FALSE)
    same <- identical(vn$cell_df, vo$cell_df) &&
      identical(geos_write_wkb(as_geos_geometry(vn$cell_nuclei)), geos_write_wkb(as_geos_geometry(vo$cell_nuclei)))
    cd <- vn$cell_df
    m_c <- metr(cd$nucleus_area, cd$nucleus_clamped, 1 - cd$nucleus_area / cd$nucleus_area_target, cd$center_dist)
    ref <- ob14$m[ob14$m$seed_idx == s & ob14$m$variant == "V-c", names(m_c)]
    rows[[s + 1]] <- data.frame(archetype = A, fold = k, seed_idx = s, identical_db4afea = same,
                                c0_identical_S14b = nrow(ref) == 1 && identical(unname(unlist(ref)), unname(m_c)))
  }
  saveRDS(do.call(rbind, rows), f); invisible(NULL)
}

# ---- null: tassellamenti CSR e RSA(d*) rigenerati con il codice e i semi di S1.4 (null_one), parametri M2 del fold 1 ----
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
  s14 <- readRDS(file.path(S14, "null", sprintf("%s_%s_%s.rds", A, roi, model)))$meta
  o14 <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold1.rds", A, a$primary)))
  o <- E$derive_cell_geometry(cen, tv, ct_m2(o14$fits$M2, PHI[[A]]), random_seed = seed, verbose = FALSE)
  saveRDS(cbind(data.frame(archetype = A, roi = roi, model = model, seed = seed, same_as_S14 = s14[["n"]] == n && s14[["seed"]] == seed),
                invariants(o, tv$cell_territories)), f); invisible(NULL)
}

# ---- perf: sezione S (449 536 celle), parametri di C-perf di S1.4/S1.4c, phi = 1; riscaldamento + 3 ripetizioni (P-5) --
perf_run <- function() {
  f <- file.path(OUT, "perf", "perf_S.rds"); if (file.exists(f)) return(invisible(NULL))
  S <- readRDS(file.path(S14, "perf_S.rds")); cen <- S$cen; cen$cell_type <- "a4"
  ct <- data.frame(cell_type = "a4", density = 28096, nucleus_to_eq_ratio = 0.6, nucleus_area_slope = 0.5, nucleus_area_sdlog = 0.3, nucleus_nc_max = 1)
  env <- function() c(load = readLines("/proc/loadavg"), tctl = paste(system("sensors 2>/dev/null | grep -i tctl", intern = TRUE), collapse = " "))
  e0 <- env(); o <- E$derive_cell_geometry(cen, S$tv, ct, random_seed = 42, verbose = FALSE)
  times <- vapply(1:3, function(r) { gc(); system.time(E$derive_cell_geometry(cen, S$tv, ct, random_seed = 42, verbose = FALSE))[["elapsed"]] }, 0)
  e1 <- env(); inv <- invariants(o, S$tv$cell_territories)
  saveRDS(list(times = times, warmup_s = o$info$elapsed_s, inv = inv, info = o$info, md5 = wkb_md5(o), env_before = e0, env_after = e1), f)
  cat(sprintf("PERF n=%d times=%s warmup=%.1f max_iter=%d n_clipped=%d n_capped=%d\n", nrow(cen), paste(round(times, 1), collapse = " "),
              o$info$elapsed_s, o$info$max_iter, o$info$n_clipped, o$info$n_capped))
}

# ---- repro: 6 unita' (fold 1 di ogni tessuto, phi del tessuto, seme 1) rieseguite in un processo nuovo --------------
repro_one <- function(A) {
  method <- ARCH$primary[ARCH$archetype == A]; o14 <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold1.rds", A, method)))
  te <- load_real(A, o14$roi_test, method); inp <- inputs_of(te)
  o <- E$derive_cell_geometry(inp$cen, inp$tess, ct_m2(o14$fits$M2, PHI[[A]]), random_seed = o14$seed + 100000, verbose = FALSE)
  g <- readRDS(file.path(OUT, "gen", sprintf("%s_%s_fold1_tissue.rds", A, method)))
  saveRDS(data.frame(archetype = A, md5_gen = g$md5[g$seed_idx == 1], md5_repro = wkb_md5(o)), file.path(OUT, "repro", paste0(A, ".rds")))
  invisible(NULL)
}

if (STAGE == "gen") {
  j1 <- expand.grid(A = ARCH$archetype, k = 1:5, phimode = c("tissue", "fold"), stringsAsFactors = FALSE)
  j1$method <- ARCH$primary[match(j1$A, ARCH$archetype)]
  j2 <- expand.grid(A = c("A4", "A6"), k = 1:5, phimode = "tissue", stringsAsFactors = FALSE); j2$method <- "stardist_he"
  jobs <- rbind(j1, j2)[, c("A", "method", "k", "phimode")]; jobs <- jobs[order(jobs$phimode != "tissue"), ]
  par_run(jobs, gen_one, "gen")
} else if (STAGE == "legacy") {
  EOLD <- dcg_env_d(file.path(OUT, "04b4_db4afea.R"))
  par_run(expand.grid(A = ARCH$archetype, k = 1:5, stringsAsFactors = FALSE), legacy_one, "legacy")
} else if (STAGE == "null") {
  par_run(expand.grid(i = seq_len(nrow(rois)), model = c("CSR", "RSA"), stringsAsFactors = FALSE), null_one, "null")
} else if (STAGE == "perf") {
  perf_run()
} else if (STAGE == "repro") {
  par_run(data.frame(A = ARCH$archetype, stringsAsFactors = FALSE), repro_one, "repro")
} else if (STAGE == "dev") {
  OUT <- file.path(OUT, "dev"); for (d in c("gen", "gen_in_tissue", "gen_in_sd")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
  NSEED <- 1L; t0 <- proc.time()[[3]]; gen_one(args[3], args[4], as.integer(args[5]), args[6])
  x <- readRDS(file.path(OUT, "gen", sprintf("%s_%s_fold%s_%s.rds", args[3], args[4], args[5], args[6])))
  print(x[, intersect(c("n", "c3_geos", "c3_sf", "no_pairs", "n_invalid", "n_empty", "n_clipped", "max_iter", "mean_iter_clipped", "elapsed"), names(x))])
  cat("dev wall", round(proc.time()[[3]] - t0, 1), "s\n")
} else stop("stadio sconosciuto: ", STAGE)
