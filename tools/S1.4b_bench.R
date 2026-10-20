# tools/S1.4b_bench.R — S1.4b (BL-098): banco delle opzioni V-c (Voronoi + clamp), V-b (Voronoi, nuclei disgiunti al
# generatore), V-a (Laguerre, stessi nuclei di V-b) sui 30 ROI di S1.4 (segmentatore primario), parametri M2 e semi di S1.4.
# Pre-registrazione: results/S1.4b/S1.4b_preregistration.md (commit ccf1c0e); definizioni operative in
# results/S1.4b/S1.4b_dev_notes.md (fissate prima degli stadi ufficiali). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S14B_MUTANT=none|MU1..MU4] Rscript --vanilla tools/S1.4b_bench.R <stage> [outdir]
#   stage = prep | loo | mut | perf ; outdir default /mnt/micron/geo_spatialtrans/S1.4b ; stadi rilanciabili (file per unita').
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R"); source("tools/S1.4_mutants.R"); source("tools/S1.4b_mutants.R")
E <- dcg_env()                                          # derive_cell_geometry() di S1.4, modo "scale" (da S1.4d il predefinito e' "clip")
args <- commandArgs(TRUE); STAGE <- args[1]
OUT <- if (length(args) > 1) args[2] else "/mnt/micron/geo_spatialtrans/S1.4b"
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"
NCORES <- as.integer(Sys.getenv("S14B_CORES", "30")); GRID <- 1e-6; GAMMA <- 0.99; NSEED <- as.integer(Sys.getenv("S14B_NSEED", "10"))
for (d in c("prep", "loo", "mut", "arb_in", "arb_in_mut", "cp1_in", "cp1", "perf")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b)))
  max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
spear <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
load_real <- function(A, roi, method) readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
evset <- function(cc) cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
region_of <- function(A, roi) { w <- r3_roi_window(A, roi, rois); geos_set_precision(as_geos_geometry(w$poly), GRID) }
k32 <- LG_k32 <- sqrt(2 * pi / (32 * sin(2 * pi / 32)))
# raggio "reale" per tutte le celle: area del nucleo appaiato (a_n) se c'e', altrimenti del nucleo-generatore (a_own)
real_area_all <- function(cc) ifelse(cc$matched & is.finite(cc$a_n) & cc$a_n > 0, cc$a_n, cc$a_own)

# ---- C3 / C-NO / C-L3 di una configurazione V-a ------------------------------------------------------------------------
va_checks <- function(L, x, y, Rc, pd) {
  nuc <- L$nucleus_polys(x, y, Rc)
  ok <- !geos_is_empty(pd$cells)
  win <- rep(FALSE, length(x)); win[ok] <- geos_within(nuc[ok], pd$cells[ok])
  # sf indipendente; conversione via WKB: st_as_sfc(<geos con griglia 1e-6>) eredita precision = 1e-6, che sf legge come
  # fattore di scala (arrotonda a 1e6 µm, area 0) — trappola trovata in sviluppo (dev notes D-2)
  c3_sf <- { w <- st_within(st_as_sfc(geos_write_wkb(nuc)), st_as_sfc(geos_write_wkb(pd$cells))); sum(!vapply(seq_along(nuc), function(i) i %in% w[[i]], TRUE)) }
  tr <- geos_strtree(nuc); q <- geos_strtree_query(tr, nuc)
  pi_ <- rep(seq_along(q), lengths(q)); pj <- unlist(q); k <- pi_ < pj
  cno <- if (any(k)) sum(geos_intersects(nuc[pi_[k]], nuc[pj[k]])) else 0L
  trc <- geos_strtree(pd$cells); qc <- geos_strtree_query(trc, pd$cells)
  ci <- rep(seq_along(qc), lengths(qc)); cj <- unlist(qc); kk <- ci < cj
  ov <- if (any(kk)) geos_area(geos_intersection_prec(pd$cells[ci[kk]], pd$cells[cj[kk]], GRID)) else 0
  A0 <- geos_area(pd$region)
  c(c3_viol = sum(!win), c3_sf_viol = c3_sf, cno_pairs = cno, cl3_rel = abs(sum(pd$area) - A0) / A0,
    cl3_max_overlap = max(c(0, ov)), n_hidden = pd$n_hidden, n_multi = pd$n_multi, n_clipped = pd$n_clipped)
}
arb_export <- function(f, set_id, x, y, w, pd) {
  write_parquet(data.frame(set_id = set_id, x = x, y = y, w = w, wkb = I(geos_write_wkb(pd$cells))), f)
  writeBin(geos_write_wkb(pd$region)[[1]], sub("\\.parquet$", "_region.wkb", f))
}

# ---- stadio prep: riferimenti reali per ROI (B-7, B-8, B-9 reale) + input di CP-1 -------------------------------------
prep_one <- function(A, roi, method) {
  f <- file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, roi, method)); if (file.exists(f)) return(invisible(NULL))
  L <- lag_env(); r <- load_real(A, roi, method); cc <- r$cells; ev <- evset(cc)
  write_parquet(data.frame(idx = as.integer(cc$idx), interior = cc$interior), file.path(OUT, "cp1_in", sprintf("%s_%s_%s.parquet", A, roi, method)))
  reg <- region_of(A, roi)
  ra <- real_area_all(cc); rr <- sqrt(ra / pi)
  b7 <- mean(sqrt(cc$a_n[ev] / pi) > cc$d_gen[ev])
  ov1 <- L$near_pairs(cc$x, cc$y, rr, gamma = 1); ovr <- L$near_pairs(cc$x, cc$y, k32 * rr, gamma = GAMMA)
  hit <- function(p) { h <- rep(FALSE, nrow(cc)); h[c(p[, "i"], p[, "j"])] <- TRUE; mean(h[ev]) }
  pdr <- L$power_diagram(cc$x, cc$y, (k32 * rr)^2, reg)
  saveRDS(list(archetype = A, roi = roi, method = method, n = nrow(cc), n_ev = sum(ev), b7_real = b7,
               b8_circle = hit(ov1), b8_rule = hit(ovr), at_lag_real = pdr$area, lag_real_hidden = pdr$n_hidden,
               reg_area = geos_area(reg)), f)
  invisible(NULL)
}

# ---- stadio loo: per fold, semi s = 0 (C-0) e 1..10 (verdetti) ------------------------------------------------------
loo_one <- function(A, method, k, mut = "none") {
  odir <- if (mut == "none") "loo" else "mut"
  f <- file.path(OUT, odir, sprintf("%s_%s_fold%d%s.rds", A, method, k, if (mut == "none") "" else paste0("_", mut)))
  if (file.exists(f)) return(invisible(NULL))
  L <- lag_mutate(lag_env(), mut)
  o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  te <- load_real(A, o$roi_test, method); cc <- te$cells; ev <- evset(cc); n <- te$n
  reg <- region_of(A, o$roi_test); bnd <- geos_boundary(reg)
  pr <- readRDS(file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, o$roi_test, method)))
  real_nc <- cc$a_n[ev] / cc$a_t[ev]; real_an <- cc$a_n[ev]; real_sp <- spear(real_an, cc$a_t[ev])
  stopifnot(isTRUE(all.equal(real_sp, o$real_sp)))
  f2 <- o$fits$M2
  ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]],
                   nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]])
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c")
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t))
  metr <- function(a_fin, trunc, loss, rdist) c(ks_nc = ks_d(a_fin[ev] / cc$a_t[ev], real_nc), ks_an = ks_d(a_fin[ev], real_an),
    d_sp = abs(spear(a_fin[ev], cc$a_t[ev]) - real_sp), trunc_all = mean(trunc), trunc_ev = mean(trunc[ev]),
    loss_all = mean(loss), b7 = mean(sqrt(a_fin[ev] / pi) > rdist[ev]))
  rows <- list(); vchk <- list(); seeds <- if (mut == "none") 0:NSEED else 1L
  ccv <- c(o$sims$M2[c("ks_nc", "ks_an", "d_sp", "clamp_all")])
  for (s in seeds) {
    sd <- o$seed + 100000 * s
    vc <- E$derive_cell_geometry(nucleus_model = "scale", cen, tess, ct, random_seed = sd, verbose = FALSE)
    a_tgt <- vc$cell_df$nucleus_area_target
    m_c <- metr(vc$cell_df$nucleus_area, vc$cell_df$nucleus_clamped, 1 - vc$cell_df$nucleus_area / a_tgt, vc$cell_df$center_dist)
    vb <- L$vab_nuclei(cc$x, cc$y, a_tgt, bnd, GAMMA); nd <- vb$df
    m_b <- metr(nd$a_n, nd$s < 1, 1 - nd$s^2, cc$d_gen)
    pd <- L$power_diagram(cc$x, cc$y, nd$Rc^2, reg)
    chk <- va_checks(L, cc$x, cc$y, nd$Rc, pd)
    nuc_b <- L$nucleus_polys(cc$x, cc$y, nd$Rc)
    c3_b <- sum(!geos_within(nuc_b, as_geos_geometry(te$terr)))
    at_l <- pd$area
    b9 <- c(dat_med = median(abs(at_l[ev] - cc$a_t[ev]) / cc$a_t[ev]), cv_vor = sd(cc$a_t[ev]) / mean(cc$a_t[ev]),
            cv_lag = sd(at_l[ev]) / mean(at_l[ev]), cv_lag_real = sd(pr$at_lag_real[ev]) / mean(pr$at_lag_real[ev]),
            ks_at_lag_real = ks_d(at_l[ev], pr$at_lag_real[ev]), ks_at_vor = ks_d(at_l[ev], cc$a_t[ev]),
            dat_real_med = median(abs(pr$at_lag_real[ev] - cc$a_t[ev]) / cc$a_t[ev]),
            trunc_pair = mean(nd$s_pair < 1), trunc_reg = mean(nd$s_reg < 1), n_pairs = vb$n_pairs)
    rows[[length(rows) + 1]] <- rbind(
      data.frame(seed_idx = s, variant = "V-c", t(m_c)), data.frame(seed_idx = s, variant = "V-b", t(m_b)),
      data.frame(seed_idx = s, variant = "V-a", t(m_b)))
    vchk[[length(vchk) + 1]] <- data.frame(seed_idx = s, t(chk), c3_b_viol = c3_b, t(b9),
      c0_match = if (s == 0) identical(unname(m_c[c("ks_nc", "ks_an", "d_sp")]), unname(ccv[c("ks_nc", "ks_an", "d_sp")])) &&
                   identical(unname(mean(vc$cell_df$nucleus_clamped)), unname(ccv[["clamp_all"]]))
                 else if (mut == "none") identical(unname(m_c[c("ks_nc", "ks_an", "d_sp", "trunc_all")]), unname(o$seeds[s, c("ks_nc", "ks_an", "d_sp", "clamp_all")])) else NA)
    if (s == 1) arb_export(file.path(OUT, if (mut == "none") "arb_in" else "arb_in_mut", sprintf("%s_%s_fold%d%s.parquet", A, method, k, if (mut == "none") "" else paste0("_", mut))),
                           1000 * match(A, ARCH$archetype) + k, cc$x, cc$y, nd$Rc^2, pd)
  }
  yard_b7 <- abs(pr$b7_real - vapply(setdiff(rois$roi_id[rois$archetype == A], o$roi_test), function(j)
    readRDS(file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, j, method)))$b7_real, 0))
  saveRDS(list(archetype = A, method = method, fold = k, roi_test = o$roi_test, seed = o$seed, mutant = mut,
               n = n, n_ev = sum(ev), m = do.call(rbind, rows), checks = do.call(rbind, vchk),
               yard = o$yard, yard_b7 = yard_b7, b7_real = pr$b7_real, b8 = c(circle = pr$b8_circle, rule = pr$b8_rule)), f)
  invisible(NULL)
}

par_run <- function(jobs, fun, label) {
  r <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(do.call(fun, as.list(jobs[i, ])), error = function(e) conditionMessage(e)),
                mc.cores = NCORES, mc.preschedule = FALSE)
  err <- vapply(r, function(z) is.character(z), TRUE)
  cat(sprintf("[%s] %d unita', %d errori\n", label, nrow(jobs), sum(err)))
  if (any(err)) { print(cbind(jobs[err, , drop = FALSE], msg = unlist(r[err]))); quit(status = 1) }
}

if (STAGE == "prep") {
  jobs <- merge(rois[, c("archetype", "roi_id")], ARCH, by = "archetype"); names(jobs) <- c("A", "roi", "method")
  par_run(jobs, prep_one, "prep")
} else if (STAGE == "loo") {
  jobs <- expand.grid(A = ARCH$archetype, k = 1:5, stringsAsFactors = FALSE); jobs$method <- ARCH$primary[match(jobs$A, ARCH$archetype)]
  par_run(jobs[, c("A", "method", "k")], loo_one, "loo")
} else if (STAGE == "mut") {
  jobs <- expand.grid(A = ARCH$archetype, mut = c("MU1", "MU2", "MU3", "MU4"), stringsAsFactors = FALSE)
  jobs$method <- ARCH$primary[match(jobs$A, ARCH$archetype)]; jobs$k <- 1L
  par_run(jobs[, c("A", "method", "k", "mut")], loo_one, "mut")
} else if (STAGE == "perf") {
  # C-perf: Laguerre su S (449 536 celle; S1.4 perf_S.rds); aree estratte con i parametri di S1.4 C-perf (non cronometrate)
  L <- lag_env(); S <- readRDS(file.path(S14, "perf_S.rds")); cen <- S$cen; cen$cell_type <- "a4"
  ct <- data.frame(cell_type = "a4", density = 28096, nucleus_to_eq_ratio = 0.6, nucleus_area_slope = 0.5, nucleus_area_sdlog = 0.3)
  o <- E$derive_cell_geometry(nucleus_model = "scale", cen, S$tv, ct, random_seed = 42, verbose = FALSE)
  reg <- geos_set_precision(as_geos_geometry(st_polygon(list(rbind(c(0, 0), c(4000, 0), c(4000, 4000), c(0, 4000), c(0, 0))))), GRID)
  gc(reset = TRUE); t0 <- proc.time()[["elapsed"]]
  vb <- L$vab_nuclei(cen$x, cen$y, o$cell_df$nucleus_area_target, geos_boundary(reg), GAMMA)
  t1 <- proc.time()[["elapsed"]]
  pd <- L$power_diagram(cen$x, cen$y, vb$df$Rc^2, reg)
  t2 <- proc.time()[["elapsed"]]
  nuc <- L$nucleus_polys(cen$x, cen$y, vb$df$Rc); nw <- sum(!geos_within(nuc, pd$cells))
  cat(sprintf("PERF laguerre n=%d nuclei_s=%.1f diagram_s=%.1f total_s=%.1f hidden=%d not_within=%d trunc=%.4f cl3_rel=%.2e\n",
              nrow(cen), t1 - t0, t2 - t1, t2 - t0, pd$n_hidden, nw, mean(vb$df$s < 1), abs(sum(pd$area) - geos_area(reg)) / geos_area(reg)))
} else if (STAGE == "dev") {
  # sviluppo: un fold (args 3..5 = A, method, k), con i prep dei ROI dell'archetipo
  A <- args[3]; method <- args[4]; k <- as.integer(args[5])
  for (r in rois$roi_id[rois$archetype == A]) { t0 <- proc.time()[[3]]; prep_one(A, r, method); cat("prep", r, round(proc.time()[[3]] - t0, 1), "s\n") }
  t0 <- proc.time()[[3]]; loo_one(A, method, k); cat("loo", round(proc.time()[[3]] - t0, 1), "s\n")
  x <- readRDS(file.path(OUT, "loo", sprintf("%s_%s_fold%d.rds", A, method, k))); print(x$m); print(x$checks)
} else stop("stadio sconosciuto: ", STAGE)
