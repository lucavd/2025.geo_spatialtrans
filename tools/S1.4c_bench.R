# tools/S1.4c_bench.R — S1.4c (BL-098 riformulata): forma del nucleo, scala (K-c) contro ritaglio nel territorio (K-1,
# K-2), sui 30 tassellamenti primari di S1.4, parametri M2 e semi di S1.4/S1.4b. Pre-registrazione:
# results/S1.4c/S1.4c_preregistration.md (commit c41cf8b); definizioni operative: results/S1.4c/S1.4c_dev_notes.md.
# Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.4c_bench.R <stage> [outdir]
#   stage = prep | gen | obs | mut | perf | dev ; outdir default /mnt/micron/geo_spatialtrans/S1.4c ; rilanciabili.
#   Ordine: prep -> gen -> (python raster real, sim) -> obs. Il raster e' tools/S1.4c_raster.py.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
source("R/04b1_extract_regions.R"); source("R/04b3_tessellate_voronoi.R"); source("tools/R3_voronoi_metrics.R")
source("tools/S1.4_mutants.R"); source("tools/S1.4c_shape.R")
E <- dcg_env()                                          # derive_cell_geometry() del pacchetto, modo "scale" (K-c; da S1.4d il predefinito e' "clip")
L <- new.env(); sys.source("tools/S1.4b_laguerre.R", envir = L)   # solo near_pairs() per B-8
args <- commandArgs(TRUE); STAGE <- args[1]
OUT <- if (length(args) > 1 && !is.na(args[2])) args[2] else "/mnt/micron/geo_spatialtrans/S1.4c"
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"; S14B <- "/mnt/micron/geo_spatialtrans/S1.4b"; R3 <- "/mnt/micron/geo_spatialtrans/R3"
NCORES <- as.integer(Sys.getenv("S14C_CORES", "30")); NSEED <- as.integer(Sys.getenv("S14C_NSEED", "10"))
for (d in c("prep", "gen", "gen_in", "obs", "mut", "perf")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   shape = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "stardist_he", "cellpose_rgb", "stardist_he"),
                   stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b)))
  max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
spear <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
load_real <- function(A, roi, method) readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
evset <- function(cc) cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
rp <- function(A, roi, method, k) file.path(R3, "real", sprintf("%s_%s_%s_%s.parquet", A, roi, method, k))
TRUNC_TOL <- 1e-6                                        # area finale < (1 - 1e-6) * target -> troncato (B-4)

# ---- stadio prep: insieme di valutazione, etichette della maschera di forma, C-OP sui reali --------------------------
prep_one <- function(A, roi, method) {
  f <- file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, roi, method)); if (file.exists(f)) return(invisible(NULL))
  r <- load_real(A, roi, method); cc <- r$cells; ev <- evset(cc); n <- r$n
  if (method == "spaceranger") {            # BL-060 come S1.4 real_one: nucleo StarDist il cui centroide cade nel territorio
    sg <- as.data.frame(read_parquet(rp(A, roi, "stardist_he", "gen")))
    hit <- st_intersects(st_as_sf(sg, coords = c("x", "y")), r$terr)
    h1 <- vapply(hit, function(z) if (length(z) == 1) z else NA_integer_, 1L)
    cnt <- tabulate(h1[!is.na(h1)], n); matched <- cnt == 1
    j <- match(seq_len(n), h1); j[!matched] <- NA
    shape_label <- as.integer(sg$label[j])
    c_match <- c(matched_identical = identical(matched, cc$matched),
                 area_identical = identical(sg$area_nuc[j][matched], cc$a_n[matched]))
  } else {
    gen <- as.data.frame(read_parquet(rp(A, roi, method, "gen")))
    shape_label <- as.integer(gen$label)
    c_match <- c(matched_identical = identical(gen$x, cc$x) && identical(gen$y, cc$y),
                 area_identical = identical(gen$area_nuc, cc$a_n))
  }
  w <- r3_roi_window(A, roi, rois)
  ob <- obs_territories(cc$x, cc$y, w, E)                 # C-OP: l'operatore dei simulati sui generatori reali = S1.4
  cop <- c(c_match, at_rel_max = max(abs(ob$a_t - cc$a_t) / cc$a_t), d_abs_max = max(abs(ob$d - cc$d_gen)),
           n_excluded = attr(ob, "n_excluded"))
  pb <- readRDS(file.path(S14B, "prep", sprintf("%s_%s_%s.rds", A, roi, method)))
  saveRDS(list(archetype = A, roi = roi, method = method, n = n, ev = ev, shape_label = shape_label, cop = cop,
               b7_real = pb$b7_real, b8_real = pb$b8_circle, nc_ev = (cc$a_n / cc$a_t)[ev]), f)
  invisible(NULL)
}

# ---- stadio gen: per fold, s = 0 (C-0) e s = 1..10: K-c, K-1, K-2; check C3, C-NO, C-area; B-4; poligoni per il raster --
nuc_checks <- function(nuc, Tg) {
  emp <- geos_is_empty(nuc)
  win <- geos_within(nuc, Tg); win[emp] <- TRUE
  w2 <- st_within(st_as_sfc(geos_write_wkb(nuc[!emp])), st_as_sfc(geos_write_wkb(Tg[!emp])))
  c3_sf <- sum(!vapply(seq_along(w2), function(i) i %in% w2[[i]], TRUE))
  tr <- geos_strtree(nuc); q <- geos_strtree_query(tr, nuc)
  pi_ <- rep(seq_along(q), lengths(q)); pj <- unlist(q); k <- pi_ < pj
  ov <- if (any(k)) geos_area(geos_intersection_prec(nuc[pi_[k]], nuc[pj[k]], S14C_GRID)) else 0
  c(n_empty = sum(emp), c3_viol = sum(!win), c3_sf_viol = c3_sf, cno_max_area = max(c(0, ov)), cno_pairs = sum(ov > 1e-9))
}
write_nuc <- function(f, nuc) {
  if (file.exists(f)) return(invisible(NULL))
  write_parquet(data.frame(idx = seq_along(nuc), wkb = I(geos_write_wkb(nuc))), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
}
gen_one <- function(A, method, k, mut = "none") {
  odir <- if (mut == "none") "gen" else "mut"
  f <- file.path(OUT, odir, sprintf("%s_%s_fold%d%s.rds", A, method, k, if (mut == "none") "" else paste0("_", mut)))
  if (file.exists(f)) return(invisible(NULL))
  o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  ob14 <- readRDS(file.path(S14B, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  te <- load_real(A, o$roi_test, method); cc <- te$cells; ev <- evset(cc); n <- te$n
  train <- setdiff(rois$roi_id[rois$archetype == A], o$roi_test)
  q99 <- unname(quantile(unlist(lapply(train, function(j) readRDS(file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, j, method)))$nc_ev)), 0.99, type = 7))
  phi <- min(1, q99)
  real_nc <- cc$a_n[ev] / cc$a_t[ev]; real_an <- cc$a_n[ev]; real_sp <- spear(real_an, cc$a_t[ev])
  f2 <- o$fits$M2
  ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]],
                   nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]])
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c")
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t))
  Tg <- geos_set_precision(as_geos_geometry(te$terr), S14C_GRID)
  metr <- function(a_fin, trunc, loss, rdist) c(ks_nc = ks_d(a_fin[ev] / cc$a_t[ev], real_nc), ks_an = ks_d(a_fin[ev], real_an),
    d_sp = abs(spear(a_fin[ev], cc$a_t[ev]) - real_sp), trunc_all = mean(trunc), trunc_ev = mean(trunc[ev]),
    loss_all = mean(loss), b7 = mean(sqrt(a_fin[ev] / pi) > rdist[ev]))
  b4 <- function(V, s, area, tgt) data.frame(seed_idx = s, variant = V, trunc_all = mean(area < (1 - TRUNC_TOL) * tgt),
    trunc_ev = mean((area < (1 - TRUNC_TOL) * tgt)[ev]), loss_all = mean(pmax(0, 1 - area / tgt)))
  rows <- list(); chk <- list(); seeds <- if (mut == "none") 0:NSEED else 1L
  for (s in seeds) {
    sd <- o$seed + 100000 * s
    vc <- E$derive_cell_geometry(nucleus_model = "scale", cen, tess, ct, random_seed = sd, verbose = FALSE)
    a_tgt <- vc$cell_df$nucleus_area_target
    m_c <- metr(vc$cell_df$nucleus_area, vc$cell_df$nucleus_clamped, 1 - vc$cell_df$nucleus_area / a_tgt, vc$cell_df$center_dist)
    ref <- ob14$m[ob14$m$seed_idx == s & ob14$m$variant == "V-c", names(m_c)]
    c0 <- if (mut == "none") identical(unname(unlist(ref)), unname(m_c)) else NA
    if (s == 0) { chk[[length(chk) + 1]] <- data.frame(seed_idx = 0L, variant = "K-c", c0 = c0); next }
    k1 <- kclip_nuclei(cc$x, cc$y, a_tgt, Tg, "K1", mutant = mut)
    k2 <- kclip_nuclei(cc$x, cc$y, a_tgt, Tg, "K2", a_t = cc$a_t, phi = phi, mutant = mut)
    nk <- !k2$df$capped
    for (V in c("K-1", "K-2")) {
      kk <- if (V == "K-1") k1 else k2
      ch <- nuc_checks(kk$nuc, Tg)
      car <- if (V == "K-2") max(c(0, abs(kk$df$area - kk$df$a_star)[nk] / kk$df$a_star[nk])) else NA
      chk[[length(chk) + 1]] <- data.frame(seed_idx = s, variant = V, c0 = NA, t(ch), c_area_rel_max = car,
        n_capped = sum(kk$df$capped), max_iter = max(kk$df$n_iter), phi = phi)
    }
    chk[[length(chk) + 1]] <- data.frame(seed_idx = s, variant = "K-c", c0 = c0)
    rows[[length(rows) + 1]] <- rbind(b4("K-c", s, vc$cell_df$nucleus_area, a_tgt), b4("K-1", s, k1$df$area, a_tgt), b4("K-2", s, k2$df$area, a_tgt))
    if (mut == "none") {
      base <- file.path(OUT, "gen_in", sprintf("%s_%s_fold%d_s%02d", A, o$roi_test, k, s))
      write_nuc(paste0(base, "_Kc.parquet"), as_geos_geometry(vc$cell_nuclei))
      write_nuc(paste0(base, "_K1.parquet"), k1$nuc); write_nuc(paste0(base, "_K2.parquet"), k2$nuc)
    }
  }
  bind <- function(z) { nm <- unique(unlist(lapply(z, names)))
    do.call(rbind, lapply(z, function(d) { for (m in setdiff(nm, names(d))) d[[m]] <- NA; d[, nm] })) }
  saveRDS(list(archetype = A, method = method, fold = k, roi_test = o$roi_test, seed = o$seed, mutant = mut, n = n,
               n_ev = sum(ev), phi = phi, q99 = q99, b4 = do.call(rbind, rows), checks = bind(chk)), f)
  invisible(NULL)
}

# ---- stadio obs: misure dei simulati con l'operatore unico, distanze dai reali del ROI di test ----------------------
obs_one <- function(file, mut = "none") {
  b <- sub("\\.parquet$", "", basename(file)); p <- strsplit(b, "_")[[1]]
  A <- p[1]; roi <- p[2]; k <- as.integer(sub("fold", "", p[3])); s <- as.integer(sub("s", "", p[4])); V <- p[5]
  f <- file.path(OUT, if (mut == "none") "obs" else "mut", paste0(b, if (mut == "none") "" else paste0("_", mut), ".rds"))
  if (file.exists(f)) return(invisible(NULL))
  method <- ARCH$primary[ARCH$archetype == A]; seg <- ARCH$shape[ARCH$archetype == A]
  te <- load_real(A, roi, method); cc <- te$cells; n <- te$n
  pr <- readRDS(file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, roi, method))); ev <- pr$ev
  m <- as.data.frame(read_parquet(file.path(OUT, "sim_out", basename(file))))
  v <- function(col) { z <- rep(NA_real_, n); z[m$label] <- m[[col]]; z }
  cx <- v("cx"); cy <- v("cy"); area <- v("area"); sol <- v("sol"); nb <- v("n_bnd"); c1 <- v("n_c1"); c2 <- v("n_c2"); fo <- v("frac_out"); npx <- v("n_px")
  ob <- obs_territories(cx, cy, r3_roi_window(A, roi, rois), E, mutant = mut, gx = cc$x, gy = cc$y)
  e <- ev & is.finite(cx) & is.finite(ob$a_t)
  rs <- as.data.frame(read_parquet(file.path(OUT, "real_shape", sprintf("%s_%s_%s_shape.parquet", A, roi, seg))))
  ri <- match(pr$shape_label[ev], rs$label); stopifnot(!anyNA(ri))
  real <- list(nc = cc$a_n[ev] / cc$a_t[ev], an = cc$a_n[ev], s1 = mean(rs$n_c2[ri] > 0), s2 = rs$n_c2[ri] / rs$n_bnd[ri],
               s3 = rs$sol[ri], s1d1 = mean(rs$n_c1[ri] > 0), s2d1 = rs$n_c1[ri] / rs$n_bnd[ri],
               m1 = sum(rs$frac_out[ri] * rs$n_px[ri]) / sum(rs$n_px[ri]), m2 = mean(rs$frac_out[ri] > 0.05))
  sim <- c(b7 = mean(sqrt(area[e] / pi) > ob$d[e]), b8 = b8_frac(cx, cy, area, ev, L), s1 = mean(c2[e] > 0),
           s1d1 = mean(c1[e] > 0), m1 = sum(fo[e] * npx[e]) / sum(npx[e]), m2 = mean(fo[e] > 0.05),
           sol_med = median(sol[e]), s2_med = median(c2[e] / nb[e]))
  dist <- c(B1 = ks_d(area[e] / ob$a_t[e], real$nc), B2 = ks_d(area[e], real$an), B7 = abs(sim[["b7"]] - pr$b7_real),
            B8 = abs(sim[["b8"]] - pr$b8_real), S1 = abs(sim[["s1"]] - real$s1), S2 = ks_d(c2[e] / nb[e], real$s2),
            S3 = ks_d(sol[e], real$s3), S1d1 = abs(sim[["s1d1"]] - real$s1d1), S2d1 = ks_d(c1[e] / nb[e], real$s2d1),
            D1m1 = abs(sim[["m1"]] - real$m1), D1m2 = abs(sim[["m2"]] - real$m2))
  saveRDS(data.frame(archetype = A, roi = roi, fold = k, seed_idx = s, variant = c(Kc = "K-c", K1 = "K-1", K2 = "K-2")[[V]],
                     mutant = mut, n_ev = sum(ev), n_eval = sum(e), n_px0 = sum(ev & !is.finite(cx)),
                     n_excluded = attr(ob, "n_excluded"), n_conflict = if (nrow(m)) m$n_conflict[1] else NA,
                     t(sim), t(dist)), f)
  invisible(NULL)
}

par_run <- function(jobs, fun, label) {
  r <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(do.call(fun, as.list(jobs[i, , drop = FALSE])), error = function(e) conditionMessage(e)),
                mc.cores = NCORES, mc.preschedule = FALSE)
  err <- vapply(r, function(z) is.character(z), TRUE)
  cat(sprintf("[%s] %d unita', %d errori\n", label, nrow(jobs), sum(err)))
  if (any(err)) { print(cbind(jobs[err, , drop = FALSE], msg = unlist(r[err]))); quit(status = 1) }
}
fold_jobs <- function() { j <- expand.grid(A = ARCH$archetype, k = 1:5, stringsAsFactors = FALSE)
  j$method <- ARCH$primary[match(j$A, ARCH$archetype)]; j[, c("A", "method", "k")] }

if (STAGE == "prep") {
  jobs <- merge(rois[, c("archetype", "roi_id")], ARCH[, c("archetype", "primary")], by = "archetype"); names(jobs) <- c("A", "roi", "method")
  par_run(jobs, prep_one, "prep")
} else if (STAGE == "gen") {
  par_run(fold_jobs(), gen_one, "gen")
} else if (STAGE == "obs") {
  fs <- sort(list.files(file.path(OUT, "sim_out"), pattern = "_s(0[1-9]|10)_K[c12]\\.parquet$"))
  par_run(data.frame(file = fs, stringsAsFactors = FALSE), obs_one, "obs")
} else if (STAGE == "mut") {
  # MU1, MU2 sul fold 1 di ogni archetipo (s = 1) nello stadio gen; MU5 nello stadio obs sul fold 1, s = 1, K-c
  j <- fold_jobs(); j <- j[j$k == 1, ]
  jm <- rbind(cbind(j, mut = "MU1"), cbind(j, mut = "MU2")); jm$mut <- as.character(jm$mut)
  par_run(jm, gen_one, "mut-gen")
  fs <- sort(list.files(file.path(OUT, "sim_out"), pattern = "_fold1_s01_Kc\\.parquet$"))
  par_run(data.frame(file = fs, mut = "MU5", stringsAsFactors = FALSE), obs_one, "mut-obs")
} else if (STAGE == "perf") {
  # C-perf: K-2 su S (449 536 celle; S1.4 perf_S.rds), aree estratte con i parametri di C-perf di S1.4 (non cronometrate)
  S <- readRDS(file.path(S14, "perf_S.rds")); cen <- S$cen; cen$cell_type <- "a4"
  ct <- data.frame(cell_type = "a4", density = 28096, nucleus_to_eq_ratio = 0.6, nucleus_area_slope = 0.5, nucleus_area_sdlog = 0.3)
  o <- E$derive_cell_geometry(nucleus_model = "scale", cen, S$tv, ct, random_seed = 42, verbose = FALSE)
  a_t <- S$tv$territory_df$territory_area[match(cen$cell_id, S$tv$territory_df$cell_id)]
  gc(reset = TRUE); t0 <- proc.time()[["elapsed"]]
  Tg <- geos_set_precision(as_geos_geometry(E$.dcg_as_polys(S$tv$cell_territories)), S14C_GRID)
  k2 <- kclip_nuclei(cen$x, cen$y, o$cell_df$nucleus_area_target, Tg, "K2", a_t = a_t, phi = 1)
  t1 <- proc.time()[["elapsed"]]
  nk <- !k2$df$capped
  cat(sprintf("PERF K2 n=%d total_s=%.1f capped=%d max_iter=%d c_area_rel_max=%.2e not_within=%d\n", nrow(cen), t1 - t0,
              sum(k2$df$capped), max(k2$df$n_iter), max(abs(k2$df$area - k2$df$a_star)[nk] / k2$df$a_star[nk]), sum(!geos_within(k2$nuc, Tg))))
} else if (STAGE == "dev") {
  A <- args[3]; method <- args[4]; k <- as.integer(args[5])
  for (r in rois$roi_id[rois$archetype == A]) { t0 <- proc.time()[[3]]; prep_one(A, r, method); cat("prep", r, round(proc.time()[[3]] - t0, 1), "s\n") }
  t0 <- proc.time()[[3]]; gen_one(A, method, k); cat("gen", round(proc.time()[[3]] - t0, 1), "s\n")
  x <- readRDS(file.path(OUT, "gen", sprintf("%s_%s_fold%d.rds", A, method, k))); print(x$checks); cat("phi", x$phi, "\n")
} else stop("stadio sconosciuto: ", STAGE)
