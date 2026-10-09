# tools/S1.3_analysis.R — S1.3: verdetti (check C, C-5, C-6, C-10, check B, CP-1..3, D-1..3) e scelta del motore.
# Pre-registrazione: results/S1.3/S1.3_preregistration.md (e2f32fc) + addendum C-10. Committato prima del run (BL-059).
# Uso: Rscript --vanilla tools/S1.3_analysis.R [outdir]   → results/S1.3/S1.3_*.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(sf))
OUT <- commandArgs(TRUE)[1]; if (is.na(OUT)) OUT <- "/mnt/micron/geo_spatialtrans/S1.3"
RES <- "results/S1.3"; TOL <- 1e-9
ARCH <- data.frame(archetype = paste0("A", 1:6), primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"))
V <- list(); vr <- function(section, id, archetype, model, metric, value, threshold, outcome, prediction = NA, confirmed = NA)
  V[[length(V) + 1]] <<- data.frame(section, id, archetype, model, metric, value = as.character(value), threshold, outcome,
                                    prediction = as.character(prediction), confirmed = as.character(confirmed))
pf <- function(ok) ifelse(ok, "PASS", "FAIL")
fmt <- function(v) formatC(v, format = "e", digits = 2)

# ---- 1. check C (test) e mutanti ------------------------------------------------------------------
tst <- lapply(c("geos", "deldir"), function(b) read.csv(file.path(RES, sprintf("S1.3_test_%s_none.csv", b))))
names(tst) <- c("geos", "deldir")
for (b in names(tst)) {
  t <- tst[[b]]; t <- t[t$status != "INFO", ]
  for (ck in unique(t$check)) { k <- t$check == ck
    vr("C", ck, "-", b, "asserzioni PASS", sprintf("%d/%d", sum(t$status[k] == "PASS"), sum(k)), "100 %", pf(all(t$status[k] == "PASS"))) }
}
target <- list(M1 = "C-3", M2 = "C-1", M3 = "C-1", M4 = c("C-2", "C-7"), M5 = c("C-1", "C-4"), M6 = "C-3")
mut <- do.call(rbind, lapply(names(target), function(m) {
  f <- file.path(RES, sprintf("S1.3_test_geos_%s.csv", m))
  if (!file.exists(f)) return(data.frame(mutant = m, target = paste(target[[m]], collapse = "/"), n_fail_target = NA, detected = FALSE))
  t <- read.csv(f); nf <- sum(t$status == "FAIL" & t$check %in% target[[m]])
  data.frame(mutant = m, target = paste(target[[m]], collapse = "/"), n_fail_target = nf, detected = nf > 0)
}))
write.csv(mut, file.path(RES, "S1.3_mutants.csv"), row.names = FALSE)
vr("C", "mutanti", "-", "geos", "rilevati", sprintf("%d/6", sum(mut$detected)), "6/6", pf(all(mut$detected)))

# ---- 2. ROI reali: C-5, C-1..3, C-10a -------------------------------------------------------------
rf <- list.files(file.path(OUT, "real"), pattern = "\\.rds$", full.names = TRUE)
real <- lapply(rf, readRDS)
rr <- do.call(rbind, lapply(real, function(o) do.call(rbind, lapply(c("geos", "deldir"), function(b) {
  e <- o[[b]]; data.frame(archetype = o$archetype, roi_id = o$roi_id, method = o$method, engine = b, n = o$n, t(e$c123), t(e$c5))
}))))
rr$c5_pass <- rr$interior_identical == 1 & rr$rel_area_interior <= TOL & rr$rel_area_mono <= TOL & rr$rel_sum_area <= TOL & rr$rel_summary <= TOL
# CORREZIONE dichiarata (allineamento al testo pre-registrato, deviazione 3): nel banco ROI la regione e' MULTIPOLYGON e i territori
# multiparte (orfani su isole senza generatori) sono ammessi; la prima versione dello script li contava come FAIL.
rr$c123_pass <- rr$c1 <= TOL & rr$c2_overlap <= TOL & rr$c2_symdiff <= TOL & rr$c3_valid == 1 & rr$c3_gen_in_own == 1 & rr$n_repaired == 0
write.csv(rr, file.path(RES, "S1.3_roi_real.csv"), row.names = FALSE)
for (b in c("geos", "deldir")) { k <- rr$engine == b
  vr("C", "C-5", "A1-A6", b, "ROI identici a R3 (40)", sprintf("%d/%d; max rel area interne %s", sum(rr$c5_pass[k]), sum(k), fmt(max(rr$rel_area_interior[k]))), "40/40", pf(all(rr$c5_pass[k])))
  vr("C", "C-1..3 ROI", "A1-A6", b, "ROI con C-1/C-2/C-3 PASS", sprintf("%d/%d; max C-1 %s", sum(rr$c123_pass[k]), sum(k), fmt(max(rr$c1[k]))), "40/40", pf(all(rr$c123_pass[k]))) }
c10r <- do.call(rbind, lapply(real, function(o) data.frame(set = "real", archetype = o$archetype, roi_id = o$roi_id, model = o$method, t(o$c10a))))

# ---- 3. nulli: C-6, C-1..3, C-10a, check B ---------------------------------------------------------
nf <- list.files(file.path(OUT, "null"), pattern = "\\.rds$", full.names = TRUE)
nul <- lapply(nf, readRDS)
ns <- do.call(rbind, lapply(nul, function(o) data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep,
                                                       n_failed = o$n_failed, t(o$summary), t(o$c123))))
write.csv(ns, file.path(RES, "S1.3_null_summary.csv"), row.names = FALSE)
r3n <- read.csv("results/R3/R3_null_summary.csv")
cols <- c("n", "n_interior", "frac_interior", "mean_area", "median_area", "cv", "median_eq_r", "median_ecc_T", "mean_nsides", "var_nsides", "var_norm_area", "cv_loc")
m6 <- merge(ns[ns$model %in% c("CSR", "RSA"), c("archetype", "roi_id", "model", "rep", "n_failed", cols)],
            r3n[, c("archetype", "roi_id", "model", "rep", "n_failed", cols)], by = c("archetype", "roi_id", "model", "rep"), suffixes = c("", ".r3"))
m6$rel <- apply(sapply(cols, function(cc) abs(m6[[cc]] - m6[[paste0(cc, ".r3")]]) / pmax(abs(m6[[paste0(cc, ".r3")]]), 1e-300)), 1, max)
m6$worst_col <- cols[apply(sapply(cols, function(cc) abs(m6[[cc]] - m6[[paste0(cc, ".r3")]]) / pmax(abs(m6[[paste0(cc, ".r3")]]), 1e-300)), 1, which.max)]
write.csv(m6[, c("archetype", "roi_id", "model", "rep", "rel", "worst_col")], file.path(RES, "S1.3_c6.csv"), row.names = FALSE)
vr("C", "C-6", "A1-A6", "geos", "repliche identiche a R3", sprintf("%d/%d (attese 1200); max rel %s", sum(m6$rel <= TOL & m6$n_failed == m6$n_failed.r3), nrow(m6), fmt(max(m6$rel))),
   "1200/1200", pf(nrow(m6) == 1200 && all(m6$rel <= TOL)))
c123n <- ns$c1 <= TOL & ns$c2_overlap <= TOL & ns$c2_symdiff <= TOL & ns$c3_valid == 1 & ns$c3_gen_in_own == 1 & ns$n_repaired == 0   # multiparte ammessi (dev. 3)
vr("C", "C-1..3 nulli", "A1-A6", "geos", "tassellazioni PASS", sprintf("%d/%d", sum(c123n), nrow(ns)), "1800/1800", pf(nrow(ns) == 1800 && all(c123n)))
dn <- do.call(rbind, lapply(nul, function(o) if (!is.null(o$c123_deldir)) data.frame(t(o$c123_deldir))))
c123d <- dn$c1 <= TOL & dn$c2_overlap <= TOL & dn$c2_symdiff <= TOL & dn$c3_valid == 1 & dn$c3_gen_in_own == 1 & dn$n_repaired == 0   # multiparte ammessi (dev. 3)
vr("C", "C-1..3 nulli", "A1-A6", "deldir", "tassellazioni PASS (replica 1)", sprintf("%d/%d", sum(c123d), nrow(dn)), "90/90", pf(nrow(dn) == 90 && all(c123d)))
c10n <- do.call(rbind, lapply(nul, function(o) if (!is.null(o$c10a)) data.frame(set = "null", archetype = o$archetype, roi_id = o$roi_id, model = o$model, t(o$c10a))))
c10s <- read.csv(file.path(RES, "S1.3_c10a_syn.csv"))
c10 <- rbind(c10r, c10n, c10s[, names(c10r)])
tg <- tst$geos[tst$geos$check == "C-1" & tst$geos$metric == "max rel |area-Σterr|", ]
write.csv(c10, file.path(RES, "S1.3_c10a.csv"), row.names = FALSE)
arb <- rbind(do.call(rbind, lapply(real, function(o) if (nrow(o$arb)) data.frame(set = "real", archetype = o$archetype, roi_id = o$roi_id, model = o$method, o$arb))),
             do.call(rbind, lapply(nul, function(o) if (!is.null(o$arb) && nrow(o$arb)) data.frame(set = "null", archetype = o$archetype, roi_id = o$roi_id, model = o$model, o$arb))))
if (is.null(arb)) arb <- data.frame(set = character(0), geos_ok = logical(0), deldir_ok = logical(0), bounded = logical(0))
write.csv(arb, file.path(RES, "S1.3_c10a_arbitration.csv"), row.names = FALSE)
n_disc <- nrow(arb) + sum(c10s$n_discord); ok_g <- sum(arb$geos_ok) + sum(c10s$n_geos_ok); ok_d <- sum(arb$deldir_ok) + sum(c10s$n_deldir_ok)
vr("C10", "arbitro", "A1-A6", "semipiani", "tile discordanti: corretti geos / deldir", sprintf("%d: %d / %d (limitati %d)", n_disc, ok_g, ok_d, sum(arb$bounded)), "-", "INFO")
c10_ok <- c10$rel_area <= TOL & c10$interior_identical == 1 & c10$frag_identical == 1 & c10$nsides_interior_identical == 1
vr("C10", "C-10a", "A1-A6", "geos vs deldir", sprintf("insiemi concordi (40 reali + 90 nulli + %d sintetici/avversari)", nrow(c10s)),
   sprintf("%d/%d; max rel area %s", sum(c10_ok), nrow(c10), fmt(max(c10$rel_area))), "100 %", pf(all(c10_ok)), "PASS", all(c10_ok))

# ---- check B -----------------------------------------------------------------------------------------
rs <- read.csv("results/R3/R3_roi_summary.csv"); rs <- rs[rs$method == ARCH$primary[match(rs$archetype, ARCH$archetype)], ]
key <- c("median_eq_r", "cv", "cv_loc", "median_ecc_T")
agg <- aggregate(ns[, key], ns[, c("archetype", "roi_id", "model")], median)
mb <- merge(agg, rs[, c("archetype", "roi_id", key)], by = c("archetype", "roi_id"), suffixes = c("_sim", "_real"))
mb$B1 <- mb$median_eq_r_sim / mb$median_eq_r_real - 1; mb$B2 <- mb$cv_sim / mb$cv_real - 1
mb$B3 <- mb$cv_loc_sim / mb$cv_loc_real - 1; mb$B4 <- mb$median_ecc_T_sim - mb$median_ecc_T_real
mb$p1 <- abs(mb$B1) <= 0.05; mb$p2 <- abs(mb$B2) <= 0.10; mb$p3 <- abs(mb$B3) <= 0.10; mb$p4 <- abs(mb$B4) < 0.05
write.csv(mb, file.path(RES, "S1.3_checkB_roi.csv"), row.names = FALSE)
pred <- list(RSA = list(B1 = c("FAIL", "FAIL", "FAIL", "FAIL", "FAIL", "PASS"), B2 = c(rep("FAIL", 5), "PASS"), B3 = c(rep("FAIL", 5), "PASS"), B4 = rep("PASS", 6)),
             RSArule = list(B1 = c("FAIL", "FAIL", NA, NA, NA, NA), B2 = c(NA, NA, NA, NA, NA, "FAIL"), B3 = c(NA, NA, NA, NA, NA, "FAIL"), B4 = c(NA, NA, NA, NA, "FAIL", "FAIL")))
lab <- c(B1 = "eq_r mediano ±5 %", B2 = "CV globale ±10 %", B3 = "CV locale ±10 %", B4 = "eccentricita' |Δ|<0.05")
for (md in c("CSR", "RSA", "RSArule")) for (A in paste0("A", 1:6)) for (b in names(lab)) {
  d <- mb[mb$archetype == A & mb$model == md, ]; np <- sum(d[[sub("B", "p", b)]]); out <- pf(np >= 4)
  pr <- if (md == "CSR") NA else pred[[md]][[b]][as.integer(sub("A", "", A))]
  vr("B", b, A, c(CSR = "G1", RSA = "G2", RSArule = "G3")[md], lab[b], sprintf("%d/5 ROI; %+.3f..%+.3f", np, min(d[[b]]), max(d[[b]])), ">= 4/5", out, pr,
     if (is.na(pr)) NA else pr == out)
}
# sensibilita' (dichiarata, non pre-registrata): riferimento reale ricalcolato dal pacchetto (stessa regola dei frammenti D-S1.3.2)
rp <- do.call(rbind, lapply(real, function(o) if (o$method == ARCH$primary[ARCH$archetype == o$archetype]) data.frame(archetype = o$archetype, roi_id = o$roi_id, t(o$geos$summary[key]))))
ms <- merge(agg, rp, by = c("archetype", "roi_id"), suffixes = c("_sim", "_real"))
ms$B1 <- ms$median_eq_r_sim / ms$median_eq_r_real - 1; ms$B2 <- ms$cv_sim / ms$cv_real - 1
ms$B3 <- ms$cv_loc_sim / ms$cv_loc_real - 1; ms$B4 <- ms$median_ecc_T_sim - ms$median_ecc_T_real
ms$p1 <- abs(ms$B1) <= 0.05; ms$p2 <- abs(ms$B2) <= 0.10; ms$p3 <- abs(ms$B3) <= 0.10; ms$p4 <- abs(ms$B4) < 0.05
write.csv(ms, file.path(RES, "S1.3_checkB_roi_pkgref.csv"), row.names = FALSE)
for (md in c("RSA", "RSArule")) for (A in paste0("A", 1:6)) for (b in names(lab)) {
  d <- ms[ms$archetype == A & ms$model == md, ]; d0 <- mb[mb$archetype == A & mb$model == md, ]
  o1 <- pf(sum(d[[sub("B", "p", b)]]) >= 4); o0 <- pf(sum(d0[[sub("B", "p", b)]]) >= 4)
  vr("Bsens", b, A, c(RSA = "G2", RSArule = "G3")[md], paste(lab[b], "(rif. pacchetto)"), sprintf("%d/5 ROI", sum(d[[sub("B", "p", b)]])), ">= 4/5", o1, o0, o1 == o0)
}
# previsioni d'ordine per G3 (pre-registrate): CV (B-2, B-3) G3 < G2 in A3, A5, A6; A4 fra G2 e CSR
for (A in c("A3", "A5", "A6", "A4")) for (cc in c("cv", "cv_loc")) {
  g2 <- mb[mb$archetype == A & mb$model == "RSA", ]; g3 <- mb[mb$archetype == A & mb$model == "RSArule", ]; g1 <- mb[mb$archetype == A & mb$model == "CSR", ]
  g2 <- g2[order(g2$roi_id), ]; g3 <- g3[order(g3$roi_id), ]; g1 <- g1[order(g1$roi_id), ]
  v2 <- g2[[paste0(cc, "_sim")]]; v3 <- g3[[paste0(cc, "_sim")]]; v1 <- g1[[paste0(cc, "_sim")]]
  if (A == "A4") { ok <- sum(v3 > pmin(v1, v2) & v3 < pmax(v1, v2)); pr <- "G3 fra G2 e CSR" } else { ok <- sum(v3 < v2); pr <- "G3 < G2" }
  vr("B", paste0("ordine ", cc), A, "G3", pr, sprintf("%d/5 ROI", ok), ">= 4/5", pf(ok >= 4), "PASS", ok >= 4)
}
# B-5 KS D (descrittivo) e D-2 (lati)
load_cells <- function(o) { d <- o$cells; d <- d[d$interior, ]; d$a <- d$area / mean(d$area); d }
ks <- list(); ksd <- function(a, b) suppressWarnings(ks.test(a, b)$statistic[[1]])
for (A in paste0("A", 1:6)) {
  re <- do.call(rbind, lapply(real, function(o) if (o$archetype == A && o$method == ARCH$primary[ARCH$archetype == A]) {
    d <- o$geos$cells; d <- d[d$interior, ]; data.frame(a = d$area / mean(d$area), nsides = d$nsides) }))
  for (md in c("CSR", "RSA", "RSArule")) {
    s1 <- do.call(rbind, lapply(nul, function(o) if (o$archetype == A && o$model == md && o$rep == 1) load_cells(o)))
    s2 <- do.call(rbind, lapply(nul, function(o) if (o$archetype == A && o$model == md && o$rep == 2) load_cells(o)))
    ks[[length(ks) + 1]] <- data.frame(archetype = A, model = md, D_real_sim = ksd(re$a, s1$a), D_noise = ksd(s1$a, s2$a),
                                       mean_nsides_real = mean(re$nsides), mean_nsides_sim = mean(s1$nsides),
                                       var_nsides_real = var(re$nsides), var_nsides_sim = var(s1$nsides))
  }
}
ks <- do.call(rbind, ks); write.csv(ks, file.path(RES, "S1.3_checkB_ks_nsides.csv"), row.names = FALSE)
for (i in seq_len(nrow(ks))) vr("B", "B-5", ks$archetype[i], c(CSR = "G1", RSA = "G2", RSArule = "G3")[ks$model[i]], "KS D (rumore)",
                                sprintf("%.3f (%.3f)", ks$D_real_sim[i], ks$D_noise[i]), "descrittivo", "INFO")

# ---- 4. controprove e descrittivi ------------------------------------------------------------------
cp1 <- read.csv(file.path(RES, "S1.3_cp1_regions.csv")); k <- cp1$n_cells >= 3 & cp1$n_cells <= 1000
f_des <- mean(cp1$uncovered_design[k] > 0.01); n_pkg <- sum(abs(cp1$uncovered_pkg[cp1$n_cells > 0]) > TOL)
vr("CP", "CP-1", "sintetici", "testo del design", "frazione regioni (3<=n<=1000) scoperte > 1 %", sprintf("%.3f (n=%d); mediana %.4f", f_des, sum(k), median(cp1$uncovered_design[k])), ">= 0.5", pf(f_des >= 0.5), "PASS", f_des >= 0.5)
vr("CP", "CP-1", "sintetici", "D-S1.3.1", "regioni scoperte > 1e-9", n_pkg, "0", pf(n_pkg == 0), "PASS", n_pkg == 0)
cp2 <- read.csv(file.path(RES, "S1.3_cp2_anisotropy.csv")); a2 <- aggregate(median_ecc_T ~ k, cp2, median); a2$delta <- a2$median_ecc_T - a2$median_ecc_T[a2$k == 1]
mono <- all(diff(a2$median_ecc_T) > 0); kmin <- suppressWarnings(min(a2$k[a2$delta >= 0.05]))
vr("CP", "CP-2", "A6 (densita')", "RSA affine", "eccentricita' monotona in k", paste(sprintf("%.3f", a2$median_ecc_T), collapse = " / "), "monotona", pf(mono), "PASS", mono)
vr("CP", "CP-2", "A6 (densita')", "RSA affine", "k minimo con Δ >= 0.05", kmin, "<= 1.5", pf(is.finite(kmin) && kmin <= 1.5), "PASS", is.finite(kmin) && kmin <= 1.5)
if (file.exists(file.path(RES, "S1.3_cp2b_posthoc_regular.csv"))) {     # POST HOC dichiarato
  c2b <- read.csv(file.path(RES, "S1.3_cp2b_posthoc_regular.csv")); b2 <- aggregate(median_ecc_T ~ k, c2b, median); b2$delta <- b2$median_ecc_T - b2$median_ecc_T[b2$k == 1]
  kmin2 <- suppressWarnings(min(b2$k[b2$delta >= 0.05]))
  vr("CP", "CP-2b (post hoc)", "A6 (densita')", "RSA regolare affine", "eccentricita' mediana per k", paste(sprintf("%.3f", b2$median_ecc_T), collapse = " / "), "monotona",
     pf(all(diff(b2$median_ecc_T) > 0)), "PASS", all(diff(b2$median_ecc_T) > 0))
  vr("CP", "CP-2b (post hoc)", "A6 (densita')", "RSA regolare affine", "k minimo con Δ >= 0.05", kmin2, "<= 1.5", pf(is.finite(kmin2) && kmin2 <= 1.5), "PASS", is.finite(kmin2) && kmin2 <= 1.5)
}
cf <- list.files(file.path(OUT, "cp3"), pattern = "\\.rds$", full.names = TRUE); cp3 <- lapply(cf, readRDS)
c3 <- do.call(rbind, lapply(cp3, function(o) data.frame(archetype = o$archetype, roi_id = o$roi_id, n = o$n, n_clipped_nonconvex = o$n_clipped_nonconvex,
  ov_contained = o$ov_contained$n, ov_m4 = o$ov_m4$n, ov_m4_not_boundary = o$ov_m4$n_bad, t(o$gap_frac))))
write.csv(c3, file.path(RES, "S1.3_cp3.csv"), row.names = FALSE)
vr("CP", "CP-3", "A1-A6", "Chaikin senza contenimento", "coppie sovrapposte (ROI con > 0)", sprintf("%d (%d/30)", sum(c3$ov_m4), sum(c3$ov_m4 > 0)), "> 0", pf(sum(c3$ov_m4) > 0), "PASS", sum(c3$ov_m4) > 0)
vr("CP", "CP-3", "A1-A6", "Chaikin senza contenimento", "coppie senza territorio ritagliato non convesso", sum(c3$ov_m4_not_boundary), "0", pf(sum(c3$ov_m4_not_boundary) == 0), "PASS", sum(c3$ov_m4_not_boundary) == 0)
vr("CP", "CP-3", "A1-A6", "D-S1.3.3", "coppie sovrapposte con contenimento", sum(c3$ov_contained), "0", pf(sum(c3$ov_contained) == 0), "PASS", sum(c3$ov_contained) == 0)
d1 <- read.csv(file.path(RES, "S1.3_d1_boundary.csv"))
vr("D", "D-1", "A1 (densita')", "RSA regola", "area fascia regione 2 / regione 1 (mediana 20 seed)", sprintf("%.3f", median(d1$ratio_band)), ">= 1.1 (previsione)", "INFO", ">= 1.1", median(d1$ratio_band) >= 1.1)
vr("D", "D-3", "sintetici", "-", "frazione ritagliate: Spearman con log area regione", sprintf("%.2f", suppressWarnings(cor(cp1$clipped_frac[cp1$n_cells > 0], log(cp1$area_um2[cp1$n_cells > 0]), method = "spearman"))), "descrittivo", "INFO")

# ---- 5. prestazioni e scelta del motore -------------------------------------------------------------
pf_ <- read.csv(file.path(RES, "S1.3_perf.csv"))
slope <- sapply(c("geos", "deldir"), function(b) { d <- pf_[pf_$backend == b & pf_$n <= 160000 & is.finite(pf_$elapsed_s), ]; unname(coef(lm(log(elapsed_s) ~ log(n), d))[2]) })
cperf <- sapply(c("geos", "deldir"), function(b) { d <- pf_[pf_$backend == b & pf_$n == 449536, ]; nrow(d) == 1 && d$status == "ok" && is.finite(d$elapsed_s) && d$elapsed_s <= 150 && d$max_rss_mb <= 4096 })
for (b in c("geos", "deldir")) { d <- pf_[pf_$backend == b & pf_$n == 449536, ]
  vr("C", "C-perf", "-", b, "449 536 punti: s, MB", sprintf("%s s, %s MB (%s)", d$elapsed_s, d$max_rss_mb, d$status), "<= 150 s, <= 4096 MB", pf(cperf[[b]]))
  vr("C10", "C-10c", "-", b, "esponente di scala (1e4-1.6e5)", sprintf("%.2f", slope[[b]]), "descrittivo", "INFO", if (b == "geos") "≈ 1" else "≈ 2") }
rob <- sapply(c("geos", "deldir"), function(b) {
  t <- tst[[b]][tst[[b]]$check == "ROB", ]; v <- do.call(rbind, strsplit(t$value, "/")); s_test <- sum(as.integer(v))
  s_real <- sum(rr$n_snapped[rr$engine == b] + rr$n_multipart[rr$engine == b] + rr$n_repaired[rr$engine == b])
  s_null <- if (b == "geos") sum(ns$n_snapped + ns$n_multipart + ns$n_repaired) else sum(dn$n_snapped + dn$n_multipart + dn$n_repaired)
  c(test = s_test, real = s_real, null = s_null, total = s_test + s_real + s_null) })
K <- sapply(c("geos", "deldir"), function(b) {
  t <- tst[[b]]; t <- t[t$status != "INFO", ]
  k1 <- all(t$status == "PASS"); k2 <- all(rr$c5_pass[rr$engine == b] & rr$c123_pass[rr$engine == b])
  k3 <- if (b == "geos") all(c123n) else all(c123d)
  k4 <- (if (b == "geos") ok_g else ok_d) == n_disc        # addendum: escluso il motore che sbaglia un caso verificabile
  c(K1 = k1, K2 = k2, K3 = k3, K_arb = k4) })
elig <- colnames(K)[apply(K, 2, all)]
choice <- if (!length(elig)) NA_character_ else if (length(elig) == 1) elig else {
  if (sum(cperf[elig]) == 1) elig[cperf[elig]] else {
    tot <- rob["total", elig]; if (length(unique(tot)) > 1) elig[which.min(tot)] else {
      tt <- sapply(elig, function(b) pf_$elapsed_s[pf_$backend == b & pf_$n == 449536]); elig[which.min(tt)] } } }
eng <- data.frame(engine = c("geos", "deldir"), K1 = K["K1", ], K2 = K["K2", ], K3 = K["K3", ], K_arb = K["K_arb", ], eligible = c("geos", "deldir") %in% elig,
                  cperf = cperf, slope = slope, robustness_events = rob["total", ], rob_test = rob["test", ], rob_real = rob["real", ], rob_null = rob["null", ],
                  chosen = c("geos", "deldir") %in% choice)
write.csv(eng, file.path(RES, "S1.3_engine_choice.csv"), row.names = FALSE)
vr("C10", "scelta", "-", paste(elig, collapse = "+"), "motore del pacchetto (regola dell'addendum)", choice, "-", if (is.na(choice)) "FAIL" else "INFO", "geos", identical(choice, "geos"))

# POST HOC dichiarato: copertura per punti sui casi con C-2 (unione) fallito
if (file.exists(file.path(RES, "S1.3_c2_coverage.csv"))) {
  cv <- read.csv(file.path(RES, "S1.3_c2_coverage.csv"))
  for (e in c("geos", "deldir")) { d <- cv[cv$engine == e, ]
    vr("C", "C-2 copertura (post hoc)", "A1-A6 + I4", e, "punti coperti 0/1/≥2; coppie sovrapposte; area",
       sprintf("%d/%d/%d; %d coppie; %.2e um2 (%d casi)", sum(d$cov0), sum(d$cov1), sum(d$cov2), sum(d$n_pairs), sum(d$ov_area), nrow(d)), "0 / tutti / 0", pf(sum(d$cov0) == 0 && sum(d$cov2) == 0)) }
}
VV <- do.call(rbind, V); write.csv(VV, file.path(RES, "S1.3_verdicts.csv"), row.names = FALSE)
print(eng, row.names = FALSE)
print(VV[VV$section != "B" | VV$id != "B-5", c("section", "id", "archetype", "model", "value", "outcome", "prediction", "confirmed")], row.names = FALSE)
