# tools/S1.4c_analysis.R — S1.4c: aggregazione, check C, verdetti B, bootstrap, controprove CP-0/1/2, previsioni, decisione.
# Scritto e committato prima degli stadi ufficiali (BL-090 d); funzioni pure in tools/S1.4c_decision.R (testate).
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.4c_analysis.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(geometry) })
source("tools/S1.4c_decision.R")
args <- commandArgs(TRUE); OUT <- if (length(args)) args[1] else "/mnt/micron/geo_spatialtrans/S1.4c"
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"; RES <- "results/S1.4c"; dir.create(RES, showWarnings = FALSE)
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   shape = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "stardist_he", "cellpose_rgb", "stardist_he"), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b)))
  max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
VAR <- c("K-c", "K-1", "K-2"); NB <- 200L
chk <- list(); addc <- function(id, what, value, thr, verdict) chk[[length(chk) + 1]] <<- data.frame(check = id, what = what, value = as.character(value), threshold = thr, verdict = verdict)

# ---- riferimenti reali per ROI -------------------------------------------------------------------------------------
prep <- list(); real <- list()
for (i in seq_len(nrow(rois))) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; met <- ARCH$primary[ARCH$archetype == A]; seg <- ARCH$shape[ARCH$archetype == A]
  p <- readRDS(file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, roi, met))); prep[[paste(A, roi)]] <- p
  rs <- as.data.frame(read_parquet(file.path(OUT, "real_shape", sprintf("%s_%s_%s_shape.parquet", A, roi, seg))))
  ri <- match(p$shape_label[p$ev], rs$label); stopifnot(!anyNA(ri))
  real[[paste(A, roi)]] <- list(s1 = mean(rs$n_c2[ri] > 0), s2 = rs$n_c2[ri] / rs$n_bnd[ri], s3 = rs$sol[ri],
    s1d1 = mean(rs$n_c1[ri] > 0), s2d1 = rs$n_c1[ri] / rs$n_bnd[ri], m1 = sum(rs$frac_out[ri] * rs$n_px[ri]) / sum(rs$n_px[ri]),
    m2 = mean(rs$frac_out[ri] > 0.05), b7 = p$b7_real, b8 = p$b8_real, rs = rs, ev_lab = p$shape_label[p$ev])
}
refs <- do.call(rbind, lapply(names(real), function(k) { z <- real[[k]]; s <- strsplit(k, " ")[[1]]
  data.frame(archetype = s[1], roi = s[2], n_ev = length(z$s3), s1 = z$s1, s1d1 = z$s1d1, s2_median = median(z$s2), s3_median = median(z$s3),
             b7 = z$b7, b8 = z$b8, m1 = z$m1, m2 = z$m2) }))
write.csv(refs, file.path(RES, "S1.4c_real_refs.csv"), row.names = FALSE)

# ---- metro inter-ROI per fold ----------------------------------------------------------------------------------------
yards <- list(); fold_roi <- list()
for (A in ARCH$archetype) { yards[[A]] <- list(); met <- ARCH$primary[ARCH$archetype == A]
  for (k in 1:5) {
    o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, met, k))); kr <- o$roi_test; fold_roi[[paste(A, k)]] <- kr
    js <- setdiff(rois$roi_id[rois$archetype == A], kr); R <- real[[paste(A, kr)]]; Rj <- lapply(js, function(j) real[[paste(A, j)]])
    stopifnot(identical(as.character(o$yard$j), js))
    yards[[A]][[k]] <- list(B1 = o$yard$ks_nc, B2 = o$yard$ks_an,
      B7 = vapply(Rj, function(z) abs(R$b7 - z$b7), 0), B8 = vapply(Rj, function(z) abs(R$b8 - z$b8), 0),
      S1 = vapply(Rj, function(z) abs(R$s1 - z$s1), 0), S2 = vapply(Rj, function(z) ks_d(R$s2, z$s2), 0),
      S3 = vapply(Rj, function(z) ks_d(R$s3, z$s3), 0), S1d1 = vapply(Rj, function(z) abs(R$s1d1 - z$s1d1), 0),
      S2d1 = vapply(Rj, function(z) ks_d(R$s2d1, z$s2d1), 0), D1m1 = vapply(Rj, function(z) abs(R$m1 - z$m1), 0),
      D1m2 = vapply(Rj, function(z) abs(R$m2 - z$m2), 0))
  } }

# ---- osservazioni dei simulati -------------------------------------------------------------------------------------
obs <- do.call(rbind, lapply(list.files(file.path(OUT, "obs"), full.names = TRUE, pattern = "\\.rds$"), readRDS))
write.csv(obs, file.path(RES, "S1.4c_obs_raw.csv"), row.names = FALSE)
stopifnot(nrow(obs) == 30 * 10 * 3)
for (i in seq_len(nrow(obs))) stopifnot(fold_roi[[paste(obs$archetype[i], obs$fold[i])]] == obs$roi[i])
fv <- do.call(rbind, lapply(split(obs, list(obs$archetype, obs$variant, obs$fold), drop = TRUE), function(z) {
  A <- z$archetype[1]; k <- z$fold[1]
  do.call(rbind, lapply(c(S14C_PANEL, "S1d1", "S2d1", "D1m1", "D1m2"), function(m) { y <- yards[[A]][[k]][[m]]; d <- median(z[[m]])
    data.frame(archetype = A, variant = z$variant[1], fold = k, roi = z$roi[1], metric = m, dist_median = d,
               yard_median = median(y), yard_max = max(y), verdict = fold_verdict(d, y)) })) }))
write.csv(fv, file.path(RES, "S1.4c_fold_verdicts.csv"), row.names = FALSE)

# ---- check C ---------------------------------------------------------------------------------------------------------
gen <- lapply(list.files(file.path(OUT, "gen"), full.names = TRUE, pattern = "\\.rds$"), readRDS)
gc_ <- do.call(rbind, lapply(gen, function(g) cbind(archetype = g$archetype, fold = g$fold, g$checks)))
c0 <- gc_$c0[gc_$variant == "K-c"]
addc("C-0", "K-c = V-c di S1.4b (fold x semi 0..10)", sprintf("%d/%d", sum(c0 %in% TRUE), length(c0)), "tutti identici", if (all(c0 %in% TRUE) && length(c0) == 330) "PASS" else "FAIL")
for (V in c("K-1", "K-2")) { z <- gc_[gc_$variant == V, ]
  addc("C3", paste(V, "nucleo dentro il territorio (geos | sf)"), sprintf("%d | %d", sum(z$c3_viol), sum(z$c3_sf_viol)), "0 | 0", if (sum(z$c3_viol) + sum(z$c3_sf_viol) == 0) "PASS" else "FAIL")
  addc("C-NO", paste(V, "coppie sovrapposte (> 1e-9 µm²); max area"), sprintf("%d; %.2e", sum(z$cno_pairs), max(z$cno_max_area)), "0", if (sum(z$cno_pairs) == 0) "PASS" else "FAIL")
  addc("C-empty", paste(V, "nuclei vuoti"), sum(z$n_empty), "0", if (sum(z$n_empty) == 0) "PASS" else "FAIL") }
z <- gc_[gc_$variant == "K-2", ]
addc("C-area", "K-2 max |A - A*| / A* senza tetto", sprintf("%.2e", max(z$c_area_rel_max)), "<= 1e-6", if (max(z$c_area_rel_max) <= 1e-6) "PASS" else "FAIL")
cop <- do.call(rbind, lapply(prep, function(p) as.data.frame(t(p$cop))))
addc("C-OP", "reali: appaiamento e aree identici a S1.4 (30 ROI)", sprintf("%d/30, %d/30", sum(cop$matched_identical == 1), sum(cop$area_identical == 1)), "30/30", if (all(cop$matched_identical == 1) && all(cop$area_identical == 1)) "PASS" else "FAIL")
addc("C-OP", "reali: territori osservati = S1.4 (max rel A_T; max |d|)", sprintf("%.1e; %.1e; esclusi %d", max(cop$at_rel_max), max(cop$d_abs_max), sum(cop$n_excluded)), "<= 1e-9", if (max(cop$at_rel_max) <= 1e-9 && max(cop$d_abs_max) <= 1e-9 && sum(cop$n_excluded) == 0) "PASS" else "FAIL")
rc <- do.call(rbind, lapply(list.files(file.path(OUT, "real_shape"), pattern = "_cop\\.csv$", full.names = TRUE), read.csv))
addc("C-OP", "maschere reali -> operatore: centroidi e aree = R3 (30 ROI)", sprintf("max %.1e", max(c(rc$dx, rc$dy, rc$darea))), "0", if (nrow(rc) == 30 && max(c(rc$dx, rc$dy, rc$darea)) == 0 && all(rc$n_gen == rc$n_lab)) "PASS" else "FAIL")
addc("C-OP", "simulati: ev senza pixel | centroidi esclusi | conflitti di pixel", sprintf("%d | %d | %d", sum(obs$n_px0), sum(obs$n_excluded), sum(obs$n_conflict, na.rm = TRUE)), "riportato", "INFO")
lg <- function(f) { x <- readLines(file.path(OUT, f)); x[grepl("^RESULT", x)] }
tst <- c(R_none = lg("test_R_none.log"), D = lg("test_decision.log"), P_none = lg("test_py_none.log"))
addc("C-R", "test (R K-1/K-2 e operatore; decisione; raster Python)", paste(tst, collapse = " | "), "tutti PASS",
     if (all(grepl("RESULT (\\d+)/\\1 PASS", tst))) "PASS" else "FAIL")
mu_t <- c(MU1 = lg("test_R_MU1.log"), MU2 = lg("test_R_MU2.log"), MU3 = lg("test_py_MU3.log"), MU4 = lg("test_py_MU4.log"), MU5 = lg("test_R_MU5.log"))
mu_det <- !grepl("RESULT (\\d+)/\\1 PASS", mu_t)
mg <- do.call(rbind, lapply(list.files(file.path(OUT, "mut"), pattern = "_MU[12]\\.rds$", full.names = TRUE), function(f) { g <- readRDS(f)
  cbind(mutant = g$mutant, archetype = g$archetype, g$checks[g$checks$variant != "K-c", ]) }))
bench_det <- c(MU1 = any(mg$c_area_rel_max[mg$mutant == "MU1" & mg$variant == "K-2"] > 1e-6),
               MU2 = any(mg$c3_viol[mg$mutant == "MU2"] > 0 | mg$c3_sf_viol[mg$mutant == "MU2"] > 0))
mo <- do.call(rbind, lapply(list.files(file.path(OUT, "mut"), pattern = "_MU5\\.rds$", full.names = TRUE), readRDS))
addc("C-mut", "mutanti rilevati (test; banco MU1, MU2)", sprintf("%s; banco MU1 %s, MU2 %s", paste(names(mu_t), ifelse(mu_det, "rilevato", "NON rilevato"), collapse = ", "),
     bench_det[["MU1"]], bench_det[["MU2"]]), "5/5", if (all(mu_det)) "PASS" else "FAIL")
pf <- readLines(file.path(OUT, "perf.log")); pl <- pf[grepl("^PERF", pf)]
rss <- as.numeric(sub(".*: ", "", pf[grepl("Maximum resident", pf)])) / 1024^2
tot <- as.numeric(sub(".*total_s=([0-9.]+).*", "\\1", pl))
addc("C-perf", "K-2 su S (449 536 celle)", sprintf("%.1f s, %.2f GB; %s", tot, rss, sub("^PERF ", "", pl)), "<= 300 s, <= 8 GB (solo WARN)", if (tot <= 300 && rss <= 8) "PASS" else "WARN")
checks <- do.call(rbind, chk); write.csv(checks, file.path(RES, "S1.4c_checks.csv"), row.names = FALSE)
shared_ok <- all(checks$verdict[checks$check %in% c("C-OP", "C-R", "C-mut")] %in% c("PASS", "INFO"))
# eleggibilita' (regola 1): tutti i check C PASS, C-perf solo WARN (escluso). Correzione dopo il run (dev notes D-6): la prima
# versione selezionava le righe con grepl("^K-2", what), che include C-perf ("K-2 su S ...") e rendeva K-2 non eleggibile.
var_ok <- function(V) { z <- checks[checks$check != "C-perf" & startsWith(checks$what, paste0(V, " ")), ]; nrow(z) > 0 && all(z$verdict == "PASS") }
eligible <- c("K-c" = shared_ok && checks$verdict[checks$check == "C-0"] == "PASS", "K-1" = shared_ok && var_ok("K-1"), "K-2" = shared_ok && var_ok("K-2"))
write.csv(data.frame(variant = names(eligible), eligible = eligible, shared_ok = shared_ok), file.path(RES, "S1.4c_eligibility.csv"), row.names = FALSE)

# ---- CP-2: robustezza dei riferimenti di forma al segmentatore ---------------------------------------------------------
cp2v <- function(A, roi, seg) { z <- as.data.frame(read_parquet(file.path(OUT, "cp2", sprintf("%s_%s_%s_cp2.parquet", A, roi, seg)))); z <- z[z$ev, ]
  list(S1 = mean(z$n_c2 > 0), S2 = z$n_c2 / z$n_bnd, S3 = z$sol, n = nrow(z)) }
dS <- function(m, a, b) if (m == "S1") abs(a$S1 - b$S1) else ks_d(a[[m]], b[[m]])
cp2 <- list(); demoted <- character(0)
for (m in c("S1", "S2", "S3")) { ex <- 0
  for (A in ARCH$archetype) { P <- ARCH$shape[ARCH$archetype == A]; O <- setdiff(c("cellpose_rgb", "stardist_he"), P); rr <- rois$roi_id[rois$archetype == A]
    vP <- lapply(rr, function(r) cp2v(A, r, P)); vO <- lapply(rr, function(r) cp2v(A, r, O))
    dseg <- vapply(seq_along(rr), function(i) dS(m, vP[[i]], vO[[i]]), 0)
    yd <- lapply(seq_along(rr), function(i) vapply(setdiff(seq_along(rr), i), function(j) dS(m, vP[[i]], vP[[j]]), 0))
    e <- cp2_exceeds(dseg, yd); ex <- ex + e
    cp2[[length(cp2) + 1]] <- data.frame(metric = m, archetype = A, primary = P, other = O, d_seg = paste(signif(dseg, 3), collapse = ";"),
      yard_max = paste(signif(vapply(yd, max, 0), 3), collapse = ";"), n_roi_above = sum(dseg > vapply(yd, max, 0)), exceeds = e,
      primary_value = if (m == "S1") paste(signif(vapply(vP, function(v) v$S1, 0), 3), collapse = ";") else paste(signif(vapply(vP, function(v) median(v[[m]]), 0), 3), collapse = ";"),
      other_value = if (m == "S1") paste(signif(vapply(vO, function(v) v$S1, 0), 3), collapse = ";") else paste(signif(vapply(vO, function(v) median(v[[m]]), 0), 3), collapse = ";")) }
  if (ex >= 3) demoted <- c(demoted, m) }
cp2 <- do.call(rbind, cp2); write.csv(cp2, file.path(RES, "S1.4c_cp2.csv"), row.names = FALSE)
panel <- setdiff(S14C_PANEL, demoted)

# ---- verdetti, punteggi, bootstrap, decisione ----------------------------------------------------------------------
vt <- verdicts(obs, yards, metrics = c(S14C_PANEL, "S1d1", "S2d1", "D1m1", "D1m2"))
vt$in_panel <- vt$metric %in% panel
write.csv(vt, file.path(RES, "S1.4c_verdicts.csv"), row.names = FALSE)
score <- score_of(vt, panel)
bs <- boot_scores(obs, yards, panel, B = NB)
write.csv(as.data.frame(bs), file.path(RES, "S1.4c_bootstrap.csv"), row.names = FALSE)
pboot <- c("K-1" = mean(bs[, "K-1"] > bs[, "K-c"]), "K-2" = mean(bs[, "K-2"] > bs[, "K-c"]))
dc <- decide(score, pboot, eligible)
b4 <- do.call(rbind, lapply(gen, function(g) { z <- g$b4; z$archetype <- g$archetype; z$fold <- g$fold; z }))
b4a <- aggregate(trunc_all ~ archetype + variant + fold, b4, median); b4a <- aggregate(trunc_all ~ archetype + variant, b4a, mean)
loss <- aggregate(loss_all ~ archetype + variant, aggregate(loss_all ~ archetype + variant + fold, b4, median), mean)
b4a <- merge(b4a, loss); b4a$verdict <- ifelse(b4a$trunc_all <= 0.05, "PASS", "FAIL")
phi <- do.call(rbind, lapply(gen, function(g) data.frame(archetype = g$archetype, fold = g$fold, q99 = g$q99, phi = g$phi)))
write.csv(b4a, file.path(RES, "S1.4c_b4.csv"), row.names = FALSE); write.csv(phi, file.path(RES, "S1.4c_phi.csv"), row.names = FALSE)

# ---- CP-0: la sovrapposizione dei cerchi di pari area predice il contatto reale? -----------------------------------
cp0 <- list(); cp0roi <- list()
for (A in ARCH$archetype) { sc_ <- c(); yy <- c(); seg <- ARCH$shape[ARCH$archetype == A]
  for (roi in rois$roi_id[rois$archetype == A]) { R <- real[[paste(A, roi)]]; rs <- R$rs
    pr <- as.data.frame(read_parquet(file.path(OUT, "real_shape", sprintf("%s_%s_%s_pairs.parquet", A, roi, seg))))
    tri <- delaunayn(cbind(rs$cx, rs$cy), options = "Qt Qbb Qc Qz")
    e <- unique(t(apply(rbind(tri[, 1:2], tri[, 2:3], tri[, c(1, 3)]), 1, sort)))
    li <- rs$label[e[, 1]]; lj <- rs$label[e[, 2]]; keep <- li %in% R$ev_lab | lj %in% R$ev_lab; e <- e[keep, , drop = FALSE]
    d <- sqrt((rs$cx[e[, 1]] - rs$cx[e[, 2]])^2 + (rs$cy[e[, 1]] - rs$cy[e[, 2]])^2)
    s <- (sqrt(rs$area[e[, 1]] / pi) + sqrt(rs$area[e[, 2]] / pi)) / d
    lo <- pmin(rs$label[e[, 1]], rs$label[e[, 2]]); hi <- pmax(rs$label[e[, 1]], rs$label[e[, 2]])
    y <- paste(lo, hi) %in% paste(pr$label_i, pr$label_j)
    cp0roi[[length(cp0roi) + 1]] <- data.frame(archetype = A, roi = roi, n_pairs = length(y), frac_contact = mean(y), auc = auc(s, y))
    sc_ <- c(sc_, s); yy <- c(yy, y) }
  cp0[[A]] <- data.frame(archetype = A, n_pairs = length(yy), frac_contact = mean(yy), auc = auc(sc_, yy), auc_ok = auc(sc_, yy) >= 0.75) }
cp0 <- do.call(rbind, cp0); write.csv(cp0, file.path(RES, "S1.4c_cp0.csv"), row.names = FALSE)
write.csv(do.call(rbind, cp0roi), file.path(RES, "S1.4c_cp0_roi.csv"), row.names = FALSE)
premise_refuted <- sum(!cp0$auc_ok[cp0$archetype %in% paste0("A", 1:4)]) >= 2

# ---- CP-1: guadagno di K-2 su K-c specifico dei tessuti densi? -------------------------------------------------------
pc <- function(V, As, ms = S14C_PANEL) sum(vt$verdict[vt$variant == V & vt$archetype %in% As & vt$metric %in% ms] == "PASS")
g14 <- pc("K-2", paste0("A", 1:4)) - pc("K-c", paste0("A", 1:4)); g56 <- pc("K-2", c("A5", "A6")) - pc("K-c", c("A5", "A6"))
cp1 <- data.frame(gain_A1_A4 = g14, gain_per_pair_A1_A4 = g14 / 2, gain_A5_A6 = g56, refuted = g56 >= g14 / 2,
                  gain_A1_A4_panel = pc("K-2", paste0("A", 1:4), panel) - pc("K-c", paste0("A", 1:4), panel),
                  gain_A5_A6_panel = pc("K-2", c("A5", "A6"), panel) - pc("K-c", c("A5", "A6"), panel))
write.csv(cp1, file.path(RES, "S1.4c_cp1.csv"), row.names = FALSE)

# ---- previsioni ----------------------------------------------------------------------------------------------------
s1m <- tapply(refs$s1, refs$archetype, mean)
solm <- do.call(rbind, lapply(ARCH$archetype, function(A) do.call(rbind, lapply(c("cellpose_rgb", "stardist_he"), function(sg)
  data.frame(archetype = A, seg = sg, sol_median = median(unlist(lapply(rois$roi_id[rois$archetype == A], function(r) cp2v(A, r, sg)$S3))))))))
write.csv(solm, file.path(RES, "S1.4c_solidity_by_segmenter.csv"), row.names = FALSE)
vv <- function(V, A, m) vt$verdict[vt$variant == V & vt$archetype == A & vt$metric == m]
pred <- data.frame(id = paste0("P", 1:6), text = c("S-1 reale >= 0.30 in A1-A4 e <= 0.15 in A5, A6", "B-4 di K-2 <= 0.05 in >= 4/6",
  "K-1: B-2 FAIL in A1-A4", "punteggio(K-2) >= punteggio(K-c) + 3", "solidita' reale mediana >= 0.95, entrambi i segmentatori, 6/6",
  "CP-0 AUC >= 0.75 in A1-A4"),
  observed = c(paste(sprintf("%s %.3f", names(s1m), s1m), collapse = "; "),
    paste(sprintf("%s %.3f", b4a$archetype[b4a$variant == "K-2"], b4a$trunc_all[b4a$variant == "K-2"]), collapse = "; "),
    paste(sprintf("%s %s", paste0("A", 1:4), vapply(paste0("A", 1:4), function(A) vv("K-1", A, "B2"), "")), collapse = "; "),
    sprintf("K-2 %d, K-c %d", score[["K-2"]], score[["K-c"]]),
    paste(sprintf("%s/%s %.3f", solm$archetype, substr(solm$seg, 1, 2), solm$sol_median), collapse = "; "),
    paste(sprintf("%s %.3f", cp0$archetype, cp0$auc), collapse = "; ")),
  confirmed = c(all(s1m[paste0("A", 1:4)] >= 0.30) && all(s1m[c("A5", "A6")] <= 0.15), sum(b4a$trunc_all[b4a$variant == "K-2"] <= 0.05) >= 4,
    all(vapply(paste0("A", 1:4), function(A) vv("K-1", A, "B2"), "") == "FAIL"), score[["K-2"]] >= score[["K-c"]] + 3,
    all(solm$sol_median >= 0.95), all(cp0$auc[cp0$archetype %in% paste0("A", 1:4)] >= 0.75)))
write.csv(pred, file.path(RES, "S1.4c_predictions.csv"), row.names = FALSE)
decision <- data.frame(dc$table, choice = dc$choice, panel = paste(panel, collapse = ","), demoted = paste(demoted, collapse = ","),
                       premise_refuted = premise_refuted, cp1_refuted = cp1$refuted)
write.csv(decision, file.path(RES, "S1.4c_decision.csv"), row.names = FALSE)
print(checks[, c("check", "what", "verdict")], right = FALSE)
cat("pannello:", panel, "| punteggi:", paste(names(score), score), "| p_boot:", paste(names(pboot), pboot), "| scelta:", dc$choice,
    "| premessa smentita:", premise_refuted, "| CP-1 smentita:", cp1$refuted, "\n")
