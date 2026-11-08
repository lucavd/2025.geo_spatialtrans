
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow) })
OUTD <- "/tmp/r7_S13_S14"
KR <- c(A1 = 1.0372736213334257, A2 = 1.0052831660148407, A3 = 1.1832782314446315, A4 = 1.211660975509211, A5 = 1.018814361785246, A6 = 1.000750783010843)
stA <- function(y, k) { mid <- (min(y) + max(y)) / 2; mid + k * (y - mid) }
KR5 <- KR; KR5[["A5"]] <- 1.361
st <- function(y, k) { mid <- (min(y) + max(y)) / 2; i1 <- which.max(y); i0 <- which.min(y); z <- y; z[i1] <- mid + k * (y[i1] - mid); z[i0] <- mid + k * (y[i0] - mid); z }
cls <- function(v, med, mx) if (v <= med) "PASS" else if (v <= mx) "WARN" else "FAIL"
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
# ---- S1.4 -----
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"
for (f in list.files(file.path(S14, "loo"), pattern = "fold", full.names = TRUE)) { o <- readRDS(f); k <- KR[[o$archetype]]
  for (m in names(o$sims)) { s <- o$sims[[m]]
    for (met in c("ks_nc", "ks_an", "d_sp")) { y <- o$yard[[met]]; y2 <- st(y, k); d <- s[[met]]
      add(step = "S1.4", archetype = o$archetype, method = o$method, fold = o$fold, roi = o$roi_test, variant = m, metric = met, d = d,
          ymin = min(y), ymed = median(y), ymax = max(y), ymed_new = median(y2), ymax_new = max(y2), kappa = k,
          v_old = cls(d, median(y), max(y)), v_new = cls(d, median(y2), max(y2)), v_alt = cls(d, median(stA(y, k)), max(stA(y, k))), v_a5 = cls(d, median(st(y, KR5[[o$archetype]])), max(st(y, KR5[[o$archetype]])))) } } }
# ---- S1.4b -----
for (f in list.files("/mnt/micron/geo_spatialtrans/S1.4b/loo", pattern = "fold", full.names = TRUE)) { o <- readRDS(f); k <- KR[[o$archetype]]
  m <- o$m[o$m$seed_idx >= 1, ]
  for (v in c("V-c", "V-a")) { x <- m[m$variant == v, ]
    d <- c(ks_nc = median(x$ks_nc), ks_an = median(x$ks_an), d_sp = median(x$d_sp), b7 = median(abs(x$b7 - o$b7_real)))
    y <- list(ks_nc = o$yard$ks_nc, ks_an = o$yard$ks_an, d_sp = o$yard$d_sp, b7 = o$yard_b7)
    for (met in names(d)) { yy <- y[[met]]; y2 <- st(yy, k)
      add(step = "S1.4b", archetype = o$archetype, method = o$method, fold = o$fold, roi = o$roi_test, variant = ifelse(v == "V-a", "V-a/V-b", v), metric = met, d = d[[met]],
          ymin = min(yy), ymed = median(yy), ymax = max(yy), ymed_new = median(y2), ymax_new = max(y2), kappa = k,
          v_old = cls(d[[met]], median(yy), max(yy)), v_new = cls(d[[met]], median(y2), max(y2)), v_alt = cls(d[[met]], median(stA(yy, k)), max(stA(yy, k))), v_a5 = cls(d[[met]], median(st(yy, KR5[[o$archetype]])), max(st(yy, KR5[[o$archetype]])))) } } }
write.csv(do.call(rbind, rows), file.path(OUTD, "folds_S14_S14b.csv"), row.names = FALSE)
cat("S1.4/S1.4b folds:", length(rows), "\n")
# ---- S1.4c (replica di tools/S1.4c_analysis.R, sola lettura) -----
source("tools/S1.4c_decision.R")
OUT <- "/mnt/micron/geo_spatialtrans/S1.4c"
ARCH <- data.frame(archetype = paste0("A", 1:6),
  primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
  shape = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "stardist_he", "cellpose_rgb", "stardist_he"), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b)))
  max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
real <- list()
for (i in seq_len(nrow(rois))) { A <- rois$archetype[i]; roi <- rois$roi_id[i]; met <- ARCH$primary[ARCH$archetype == A]; seg <- ARCH$shape[ARCH$archetype == A]
  p <- readRDS(file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, roi, met)))
  rs <- as.data.frame(read_parquet(file.path(OUT, "real_shape", sprintf("%s_%s_%s_shape.parquet", A, roi, seg))))
  ri <- match(p$shape_label[p$ev], rs$label); stopifnot(!anyNA(ri))
  real[[paste(A, roi)]] <- list(s1 = mean(rs$n_c2[ri] > 0), s2 = rs$n_c2[ri] / rs$n_bnd[ri], s3 = rs$sol[ri],
    s1d1 = mean(rs$n_c1[ri] > 0), s2d1 = rs$n_c1[ri] / rs$n_bnd[ri], m1 = sum(rs$frac_out[ri] * rs$n_px[ri]) / sum(rs$n_px[ri]),
    m2 = mean(rs$frac_out[ri] > 0.05), b7 = p$b7_real, b8 = p$b8_real) }
yards <- list(); yardsN <- list(); yardsA <- list(); yards5 <- list(); yrows <- list()
for (A in ARCH$archetype) { yards[[A]] <- list(); yardsN[[A]] <- list(); yardsA[[A]] <- list(); yards5[[A]] <- list(); met <- ARCH$primary[ARCH$archetype == A]; k <- KR[[A]]
  for (kk in 1:5) { o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, met, kk))); kr <- o$roi_test
    js <- setdiff(rois$roi_id[rois$archetype == A], kr); R <- real[[paste(A, kr)]]; Rj <- lapply(js, function(j) real[[paste(A, j)]])
    stopifnot(identical(as.character(o$yard$j), js))
    yy <- list(B1 = o$yard$ks_nc, B2 = o$yard$ks_an,
      B7 = vapply(Rj, function(z) abs(R$b7 - z$b7), 0), B8 = vapply(Rj, function(z) abs(R$b8 - z$b8), 0),
      S1 = vapply(Rj, function(z) abs(R$s1 - z$s1), 0), S2 = vapply(Rj, function(z) ks_d(R$s2, z$s2), 0),
      S3 = vapply(Rj, function(z) ks_d(R$s3, z$s3), 0), S1d1 = vapply(Rj, function(z) abs(R$s1d1 - z$s1d1), 0),
      S2d1 = vapply(Rj, function(z) ks_d(R$s2d1, z$s2d1), 0), D1m1 = vapply(Rj, function(z) abs(R$m1 - z$m1), 0),
      D1m2 = vapply(Rj, function(z) abs(R$m2 - z$m2), 0))
    yards[[A]][[kk]] <- yy; yardsN[[A]][[kk]] <- lapply(yy, st, k = k); yardsA[[A]][[kk]] <- lapply(yy, stA, k = k); yards5[[A]][[kk]] <- lapply(yy, st, k = KR5[[A]])
    for (m in names(yy)) yrows[[length(yrows) + 1]] <- data.frame(archetype = A, fold = kk, metric = m, ymin = min(yy[[m]]), ymed = median(yy[[m]]), ymax = max(yy[[m]]),
       ymed_new = median(yardsN[[A]][[kk]][[m]]), ymax_new = max(yardsN[[A]][[kk]][[m]]), kappa = k) } }
write.csv(do.call(rbind, yrows), file.path(OUTD, "S14c_yards.csv"), row.names = FALSE)
obs <- read.csv("results/S1.4c/S1.4c_obs_raw.csv", stringsAsFactors = FALSE)
stopifnot(nrow(obs) == 900)
# CP-2 con metro vecchio e nuovo
cp2v <- function(A, roi, seg) { z <- as.data.frame(read_parquet(file.path(OUT, "cp2", sprintf("%s_%s_%s_cp2.parquet", A, roi, seg)))); z <- z[z$ev, ]
  list(S1 = mean(z$n_c2 > 0), S2 = z$n_c2 / z$n_bnd, S3 = z$sol, n = nrow(z)) }
dS <- function(m, a, b) if (m == "S1") abs(a$S1 - b$S1) else ks_d(a[[m]], b[[m]])
cp2 <- list(); dem_old <- character(0); dem_new <- character(0)
for (m in c("S1", "S2", "S3")) { ex <- 0; exn <- 0
  for (A in ARCH$archetype) { P <- ARCH$shape[ARCH$archetype == A]; O <- setdiff(c("cellpose_rgb", "stardist_he"), P); rr <- rois$roi_id[rois$archetype == A]
    vP <- lapply(rr, function(r) cp2v(A, r, P)); vO <- lapply(rr, function(r) cp2v(A, r, O))
    dseg <- vapply(seq_along(rr), function(i) dS(m, vP[[i]], vO[[i]]), 0)
    yd <- lapply(seq_along(rr), function(i) vapply(setdiff(seq_along(rr), i), function(j) dS(m, vP[[i]], vP[[j]]), 0))
    ydn <- lapply(yd, st, k = KR[[A]])
    ea <- cp2_exceeds(dseg, lapply(yd, stA, k = KR[[A]])); e5 <- cp2_exceeds(dseg, lapply(yd, st, k = KR5[[A]])); e <- cp2_exceeds(dseg, yd); en <- cp2_exceeds(dseg, ydn); ex <- ex + e; exn <- exn + en
    cp2[[length(cp2) + 1]] <- data.frame(metric = m, archetype = A, n_roi_above_old = sum(dseg > vapply(yd, max, 0)), n_roi_above_new = sum(dseg > vapply(ydn, max, 0)),
      exceeds_old = e, exceeds_new = en, exceeds_alt = ea, exceeds_a5 = e5, kappa = KR[[A]]) }
  if (ex >= 3) dem_old <- c(dem_old, m); if (exn >= 3) dem_new <- c(dem_new, m) }
write.csv(do.call(rbind, cp2), file.path(OUTD, "S14c_cp2.csv"), row.names = FALSE)
panel_old <- setdiff(S14C_PANEL, dem_old); panel_new <- setdiff(S14C_PANEL, dem_new)
MET <- c(S14C_PANEL, "S1d1", "S2d1", "D1m1", "D1m2")
vt_o <- verdicts(obs, yards, metrics = MET); vt_n <- verdicts(obs, yardsN, metrics = MET); vt_a <- verdicts(obs, yardsA, metrics = MET); vt_5 <- verdicts(obs, yards5, metrics = MET)
vt <- merge(vt_o, vt_n, by = c("archetype", "variant", "metric"), suffixes = c("_old", "_new"))
names(vt_a)[4:7] <- paste0(names(vt_a)[4:7], "_alt"); vt <- merge(vt, vt_a, by = c("archetype", "variant", "metric")); v5 <- vt_5; names(v5)[4:7] <- paste0(names(v5)[4:7], "_a5"); vt <- merge(vt, v5, by = c("archetype", "variant", "metric"))
vt$in_panel_old <- vt$metric %in% panel_old; vt$in_panel_new <- vt$metric %in% panel_new
write.csv(vt, file.path(OUTD, "S14c_verdicts_oldnew.csv"), row.names = FALSE)
elig <- read.csv("results/S1.4c/S1.4c_eligibility.csv"); eligible <- setNames(elig$eligible, elig$variant)
res <- list()
for (lab in c("old", "new", "alt", "a5")) { Y <- switch(lab, old = yards, new = yardsN, alt = yardsA, a5 = yards5); pan <- if (lab == "old") panel_old else panel_new
  vv <- switch(lab, old = vt_o, new = vt_n, alt = vt_a, a5 = vt_5)
  sc <- score_of(vv, pan); bs <- boot_scores(obs, Y, pan, B = 200L)
  pb <- c("K-1" = mean(bs[, "K-1"] > bs[, "K-c"]), "K-2" = mean(bs[, "K-2"] > bs[, "K-c"]))
  dc <- decide(sc, pb, eligible); tb <- dc$table; tb$choice <- dc$choice; tb$which <- lab; tb$panel <- paste(pan, collapse = ",")
  tb$score_Kc <- as.numeric(sc[["K-c"]]); res[[lab]] <- tb
  if (lab == "old") write.csv(as.data.frame(bs), file.path(OUTD, "S14c_boot_old.csv"), row.names = FALSE) }
write.csv(do.call(rbind, res), file.path(OUTD, "S14c_decision_oldnew.csv"), row.names = FALSE)
pc <- function(vt, V, As, ms = S14C_PANEL) sum(vt$verdict[vt$variant == V & vt$archetype %in% As & vt$metric %in% ms] == "PASS")
cp1 <- do.call(rbind, lapply(list(old = vt_o, new = vt_n, alt = vt_a, a5 = vt_5), function(v) { g14 <- pc(v, "K-2", paste0("A", 1:4)) - pc(v, "K-c", paste0("A", 1:4)); g56 <- pc(v, "K-2", c("A5", "A6")) - pc(v, "K-c", c("A5", "A6")); data.frame(g14 = g14, g56 = g56, refuted = g56 >= g14 / 2) }))
cp1$which <- rownames(cp1); write.csv(cp1, file.path(OUTD, "S14c_cp1_oldnew.csv"), row.names = FALSE)
cat("DONE\n")
