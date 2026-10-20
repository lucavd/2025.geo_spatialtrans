# tools/S1.4_analysis.R — S1.4: verdetti (check C, B, controprove, previsioni) dagli stadi del banco.
# Regole: results/S1.4/S1.4_preregistration.md + definizioni operative in results/S1.4/S1.4_dev_notes.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.4_analysis.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(digest) })
source("R/04b4_derive_cell_geometry.R")
args <- commandArgs(TRUE); OUT <- if (length(args)) args[1] else "/mnt/micron/geo_spatialtrans/S1.4"
RES <- "results/S1.4"; dir.create(RES, showWarnings = FALSE)
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   sens = c(NA, NA, NA, "stardist_he", NA, "stardist_he"), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
V <- list()
vr <- function(section, id, archetype, metric, value, threshold, outcome, note = "") {
  V[[length(V) + 1]] <<- data.frame(section, id, archetype, metric, value = as.character(value), threshold, outcome, note)
  cat(sprintf("%-4s %-8s %-6s %-40s %-28s %-18s %s\n", outcome, id, archetype, substr(metric, 1, 40), substr(value, 1, 28), threshold, note))
}
pf <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
fmt <- function(v, d = 3) formatC(v, format = "g", digits = d)
fe <- function(v) formatC(v, format = "e", digits = 2)
cls <- function(v, med, mx) if (v <= med) "PASS" else if (v <= mx) "WARN" else "FAIL"
arch_v <- function(o) if (sum(o == "PASS") >= 4) "PASS" else if (sum(o == "FAIL") >= 2) "FAIL" else "WARN"
load_real <- function(A, roi, method) readRDS(file.path(OUT, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
evset <- function(cc) cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0

# ---- tabella dei fold -------------------------------------------------------------------------------------------
lf <- list.files(file.path(OUT, "loo"), pattern = "fold", full.names = TRUE)
L <- lapply(lf, readRDS)
folds <- do.call(rbind, lapply(L, function(o) do.call(rbind, lapply(names(o$sims), function(m) {
  s <- o$sims[[m]]; f <- o$fits[[m]]
  data.frame(archetype = o$archetype, method = o$method, fold = o$fold, roi_test = o$roi_test, model = m, n_ev = o$n_ev,
             rho = f[["rho"]], beta = f[["beta"]], sigma = f[["sigma"]], a_ref = f[["a_ref"]],
             ks_nc = s[["ks_nc"]], ks_an = s[["ks_an"]], d_sp = s[["d_sp"]], clamp_all = s[["clamp_all"]], c3_viol = s[["c3_viol"]], min_area = s[["min_area"]],
             y_ks_nc_med = median(o$yard$ks_nc), y_ks_nc_max = max(o$yard$ks_nc), y_ks_an_med = median(o$yard$ks_an), y_ks_an_max = max(o$yard$ks_an),
             y_d_sp_med = median(o$yard$d_sp), y_d_sp_max = max(o$yard$d_sp)) }))))
folds$S <- folds$ks_nc + folds$ks_an
folds$v_b1 <- mapply(cls, folds$ks_nc, folds$y_ks_nc_med, folds$y_ks_nc_max)
folds$v_b2 <- mapply(cls, folds$ks_an, folds$y_ks_an_med, folds$y_ks_an_max)
folds$v_b3 <- mapply(cls, folds$d_sp, folds$y_d_sp_med, folds$y_d_sp_max)
folds <- folds[order(folds$archetype, folds$method, folds$model, folds$fold), ]
write.csv(folds, file.path(RES, "S1.4_loo_folds.csv"), row.names = FALSE)
prim <- folds[folds$method == ARCH$primary[match(folds$archetype, ARCH$archetype)], ]
stopifnot(nrow(prim) == 6 * 5 * 4)

# ---- check C --------------------------------------------------------------------------------------------------------
tt <- read.csv(file.path(RES, "S1.4_test_none.csv"))
vr("C", "C-1", "-", "test unitari (test_S1.4.R)", sprintf("%d/%d", sum(tt$status == "PASS"), nrow(tt)), "tutti", pf(all(tt$status == "PASS")))
mut <- vapply(c(sprintf("M%d", 1:5), sprintf("X%d", 1:9)), function(m) { t <- read.csv(file.path(RES, sprintf("S1.4_test_%s.csv", m))); sum(t$status == "FAIL") }, 0)
vr("C", "C-mut", "-", "mutanti rilevati (M1-M5 pre-reg.; X1-X9 revisione)", sprintf("%d/14 (%s)", sum(mut > 0), paste(mut, collapse = ",")), "5/5 pre-reg.", pf(all(mut > 0)))
c3r <- do.call(rbind, lapply(L, function(o) data.frame(archetype = o$archetype, method = o$method, n = o$c3[["n"]], w = o$c3[["within_sf"]], ra = o$c3[["rel_area"]])))
viol <- sum(folds$c3_viol)
vr("C", "C3", "A1-A6", "nuclei dentro il territorio, ROI reali (sf, M2)", sprintf("%d/%d (40 set x fold)", sum(c3r$w), sum(c3r$n)), "0 violazioni", pf(all(c3r$w == c3r$n)))
vr("C", "C3", "A1-A6", "violazioni geos interne, M0-M3 reali", viol, "0", pf(viol == 0))
vr("C", "C5", "A1-A6", "min area nucleo, M0-M3 reali", fe(min(folds$min_area)), "> 0", pf(min(folds$min_area) > 0))
vr("C", "C-area", "A1-A6", "max rel |area poligono - nucleus_area| (M2 reali)", fe(max(c3r$ra)), "<= 1e-9", pf(max(c3r$ra) <= 1e-9))
nl <- list.files(file.path(OUT, "null"), full.names = TRUE)
N <- do.call(rbind, lapply(nl, function(f) { r <- readRDS(f); b <- sub(".rds", "", basename(f))
  do.call(rbind, lapply(c("circle", "ellipse"), function(v) { z <- r[[v]]
    data.frame(set = b, variant = v, error = z[["error"]], n = z[["n"]], within = z[["within_sf"]], le0 = z[["n_area_le0"]], rap = z[["rel_area_poly"]],
               rau = z[["rel_area_unclamped"]], mincyto = z[["min_cyto"]], clamp = z[["clamp"]], c7 = r$c7[["ratio"]]) })) }))
write.csv(N, file.path(RES, "S1.4_null.csv"), row.names = FALSE)
vr("C", "C3", "nulli", "nuclei dentro (CSR+RSA, cerchio+ellisse)", sprintf("%d/%d; errori %d; set %d", sum(N$within), sum(N$n), sum(N$error), nrow(N) / 2), "0 violazioni", pf(sum(N$error) == 0 && all(N$within == N$n) && nrow(N) == 120))
vr("C", "C5", "nulli", "nuclei con area <= 0", sum(N$le0), "0", pf(sum(N$le0) == 0))
vr("C", "C-area", "nulli", "max rel area poligono / target (non clampate)", sprintf("%s / %s", fe(max(N$rap)), fe(max(N$rau))), "<= 1e-9", pf(max(N$rap) <= 1e-9 && max(N$rau) <= 1e-9))
vr("C", "C-area", "nulli", "min citoplasma", fmt(min(N$mincyto)), ">= 0", pf(min(N$mincyto) >= 0))
c7 <- N[N$variant == "circle", ]
vr("C", "C7", "nulli", "p90(eq_r)/atteso, min-max su 60", sprintf("%.2f-%.2f", min(c7$c7), max(c7$c7)), "[0.5, 2]", if (all(c7$c7 >= 0.5 & c7$c7 <= 2)) "PASS" else "WARN")
rc <- do.call(rbind, lapply(list.files(file.path(OUT, "rec"), full.names = TRUE), readRDS))
rc$d_beta <- abs(rc$beta_hat - rc$beta); rc$d_sig <- abs(rc$sigma_hat - rc$sigma) / rc$sigma; rc$d_rho <- abs(rc$rho_hat - rc$rho) / rc$rho
write.csv(rc, file.path(RES, "S1.4_rec.csv"), row.names = FALSE)
okr <- rc$d_beta <= 0.05 & rc$d_sig <= 0.05 & rc$d_rho <= 0.02
vr("C", "C-rec", "A1-A6", "recupero (rho, beta, sigma) entro soglia", sprintf("%d/%d; max |dB| %.3f, dS %.3f, dR %.3f; clamp max %.3f", sum(okr), nrow(rc),
   max(rc$d_beta), max(rc$d_sig), max(rc$d_rho), max(rc$clamp)), "36/36", pf(all(okr) && nrow(rc) == 36))
# C-rip sui dati reali: due run identici (A1 r1, M2 pooled)
r1 <- load_real("A1", "r1", "cellpose_rgb")
tess <- list(cell_territories = r1$terr, territory_df = data.frame(cell_id = seq_len(r1$n), territory_area = r1$cells$a_t))
cen <- data.frame(cell_id = seq_len(r1$n), region_id = 1L, x = r1$cells$x, y = r1$cells$y, cell_type = "c")
ctr <- data.frame(cell_type = "c", density = 1e4, nucleus_to_eq_ratio = 0.45, nucleus_area_slope = 0.5, nucleus_area_sdlog = 0.3)
h <- vapply(1:2, function(i) digest(serialize(derive_cell_geometry(nucleus_model = "scale", cen, tess, ctr, random_seed = 5, verbose = FALSE)[1:2], NULL)), "")
vr("C", "C-rip", "A1", "due run sullo stesso ROI reale (md5)", substr(h[1], 1, 8), "identici", pf(h[1] == h[2]))
# C-perf
pp <- function(v) { x <- readLines(file.path(OUT, sprintf("perf_%s.log", v)))
  el <- as.numeric(sub(".*elapsed_s=([0-9.]+).*", "\\1", grep("^PERF", x, value = TRUE)))
  kb <- as.numeric(sub(".*: ", "", grep("Maximum resident", x, value = TRUE))); nw <- as.integer(sub(".*not_within=([0-9]+).*", "\\1", grep("^PERF", x, value = TRUE)))
  c(el = el, gb = kb / 1024^2, nw = nw) }
pc <- pp("circle"); pe <- pp("ellipse")
vr("C", "C-perf", "S", "449 536 celle, cerchio (tempo, RSS)", sprintf("%.1f s, %.2f GB, fuori %d", pc[["el"]], pc[["gb"]], pc[["nw"]]), "<= 300 s, <= 8 GB", pf(pc[["el"]] <= 300 && pc[["gb"]] <= 8 && pc[["nw"]] == 0))
vr("C", "C-perf", "S", "449 536 celle, ellisse 1.6 (riportata)", sprintf("%.1f s, %.2f GB, fuori %d", pe[["el"]], pe[["gb"]], pe[["nw"]]), "riportata", "INFO")
ml <- readLines(file.path(OUT, "mutM1.log"))
vr("C", "C-mut", "reali", "M1 (clamp rimosso) sullo stadio loo: errori C3", sprintf("%d/%d set x fold", length(grep("non contenuti", ml)), length(grep("^\\[loo\\]", ml))), "riportato", "INFO")

# ---- check B (M2, segmentatore primario) -----------------------------------------------------------------------------
m2 <- prim[prim$model == "M2", ]
for (A in ARCH$archetype) {
  s <- m2[m2$archetype == A, ]
  for (b in c("b1", "b2", "b3")) {
    o <- s[[paste0("v_", b)]]; met <- c(b1 = "KS N/C", b2 = "KS A_N", b3 = "|dSpearman(A_N,A_T)|")[[b]]
    val <- s[[c(b1 = "ks_nc", b2 = "ks_an", b3 = "d_sp")[[b]]]]
    vr("B", toupper(sub("b", "B-", b)), A, sprintf("%s (fold: %s)", met, paste(substr(o, 1, 1), collapse = "")), paste(fmt(val, 2), collapse = " "), "metro inter-ROI", arch_v(o))
  }
  vr("B", "B-4", A, "frazione clampata (media 5 fold, M2 cerchio)", sprintf("%.4f", mean(s$clamp_all)), "<= 0.05", pf(mean(s$clamp_all) <= 0.05))
}
# B-5 descrittive
b5 <- do.call(rbind, lapply(L, function(o) data.frame(archetype = o$archetype, method = o$method, fold = o$fold, t(o$b5))))
write.csv(b5, file.path(RES, "S1.4_b5.csv"), row.names = FALSE)
for (A in ARCH$archetype) { z <- b5[b5$archetype == A & b5$method == ARCH$primary[ARCH$archetype == A], ]
  vr("B", "B-5", A, "dTheta mediano reale | sim terr | sim rand (°)", sprintf("%.0f | %.0f | %.0f", 180 / pi * median(z$real_med_dtheta), 180 / pi * median(z$sim_terr_med_dtheta), 180 / pi * median(z$sim_rand_med_dtheta)), "descrittiva", "INFO")
  vr("B", "B-5", A, "clamp: cerchio POI | gen | ellisse terr | rand", sprintf("%.3f | %.3f | %.3f | %.3f", mean(m2$clamp_all[m2$archetype == A]), mean(z$clamp_gen), mean(z$clamp_ell_terr), mean(z$clamp_ell_rand)), "descrittiva", "INFO")
  vr("B", "B-5", A, "scarto POI-nucleo reale / eq_r (mediana, q90)", sprintf("%.2f, %.2f", median(z$poi_offset_med), median(z$poi_offset_q90)), "descrittiva", "INFO") }
# B-6: intervalli C8 per archetipo e potenza
c8 <- do.call(rbind, lapply(ARCH$archetype, function(A) { ids <- rois$roi_id[rois$archetype == A]; p <- ARCH$primary[ARCH$archetype == A]
  eq <- lapply(ids, function(r) { cc <- load_real(A, r, p)$cells; sqrt(cc$a_t[cc$interior] / pi) })
  q <- quantile(unlist(eq), c(.05, .5, .95))
  pw <- t(vapply(seq_along(ids), function(k) { qk <- quantile(unlist(eq[-k]), c(.05, .95)); md <- median(eq[[k]])
    c(spec = md >= qk[[1]] && md <= qk[[2]], x2 = 2 * md > qk[[2]], x05 = 0.5 * md < qk[[1]]) }, logical(3)))
  data.frame(archetype = A, q05 = q[[1]], q50 = q[[2]], q95 = q[[3]], spec = sum(pw[, "spec"]), det_x2 = sum(pw[, "x2"]), det_x05 = sum(pw[, "x05"])) }))
write.csv(c8, file.path(RES, "S1.4_c8.csv"), row.names = FALSE)
for (i in seq_len(nrow(c8))) vr("B", "B-6", c8$archetype[i], "C8 [q05, q95] eq_r; specif. | x2 | x0.5 (su 5)",
   sprintf("[%.2f, %.2f]; %d | %d | %d", c8$q05[i], c8$q95[i], c8$spec[i], c8$det_x2[i], c8$det_x05[i]), "riportata", "INFO")


# ---- sensibilita' al seme (10 semi M2, dichiarata «solo riportata»; aggiunta al report dopo la revisione, RC-8) -------
ss <- do.call(rbind, lapply(ARCH$archetype, function(A) { p <- ARCH$primary[ARCH$archetype == A]
  oo <- L[vapply(L, function(o) o$archetype == A && o$method == p, TRUE)]; oo <- oo[order(vapply(oo, function(o) o$fold, 1))]
  do.call(rbind, lapply(1:10, function(s) { v <- t(vapply(oo, function(o) c(
      cls(o$seeds[s, "ks_nc"], median(o$yard$ks_nc), max(o$yard$ks_nc)), cls(o$seeds[s, "ks_an"], median(o$yard$ks_an), max(o$yard$ks_an)),
      cls(o$seeds[s, "d_sp"], median(o$yard$d_sp), max(o$yard$d_sp))), character(3)))
    data.frame(archetype = A, seed = s, B1 = arch_v(v[, 1]), B2 = arch_v(v[, 2]), B3 = arch_v(v[, 3]), clamp = mean(vapply(oo, function(o) o$seeds[s, "clamp_all"], 0))) })) }))
write.csv(ss, file.path(RES, "S1.4_seed_sensitivity.csv"), row.names = FALSE)
for (A in ARCH$archetype) { z <- ss[ss$archetype == A, ]; tb <- function(x) paste(sprintf("%s %d", names(table(x)), table(x)), collapse = ", ")
  vr("B", "B-seed", A, "verdetti M2 su 10 semi: B-1 | B-2 | B-3", sprintf("%s | %s | %s", tb(z$B1), tb(z$B2), tb(z$B3)), "riportata", "INFO") }

# ---- parametri calibrati (5 ROI) e CP-2 placebo -------------------------------------------------------------------
cal <- do.call(rbind, lapply(seq_len(nrow(ARCH)), function(i) { A <- ARCH$archetype[i]
  do.call(rbind, lapply(na.omit(c(ARCH$primary[i], ARCH$sens[i])), function(p) {
    rs <- lapply(rois$roi_id[rois$archetype == A], function(r) load_real(A, r, p))
    a_ref <- sum(vapply(rs, function(r) sum(r$cells$a_t), 0)) / sum(vapply(rs, function(r) r$n, 0))
    cc <- do.call(rbind, lapply(rs, function(r) r$cells[evset(r$cells), ]))
    f <- fit_nucleus_size_model(cc$a_n, cc$a_t, a_ref, "full")
    x <- log(cc$a_t) - log(a_ref); se <- f[["sigma"]] / (sd(x) * sqrt(nrow(cc)))
    set.seed(20261015 + i); fp <- fit_nucleus_size_model(cc$a_n, sample(cc$a_t), a_ref, "full")
    own <- do.call(rbind, lapply(rs, function(r) r$cells[r$cells$interior & is.finite(r$cells$a_own) & r$cells$a_own > 0, ]))
    data.frame(archetype = A, method = p, n = nrow(cc), a_ref = a_ref, density_mm2 = 1e6 / a_ref, rho = f[["rho"]], beta = f[["beta"]], sigma = f[["sigma"]],
               se_beta = se, beta_perm = fp[["beta"]], nc_median = median(cc$a_n / cc$a_t), ratio_median = median(sqrt(cc$a_n / cc$a_t)),
               nc_median_own = median(own$a_own / own$a_t), ratio_median_own = median(sqrt(own$a_own / own$a_t)), n_own = nrow(own),
               ecc_n_median = median(cc$ecc_n, na.rm = TRUE)) })) }))
write.csv(cal, file.path(RES, "S1.4_calibration.csv"), row.names = FALSE)
cp <- cal[cal$method == ARCH$primary[match(cal$archetype, ARCH$archetype)], ]
for (i in seq_len(nrow(cp))) vr("CP", "CP-2", cp$archetype[i], "beta con A_T permutata (SE OLS)", sprintf("%.4f (SE %.4f)", cp$beta_perm[i], cp$se_beta[i]), "|b| <= 0.02", pf(abs(cp$beta_perm[i]) <= 0.02))
# ---- CP-1 modelli annidati ----------------------------------------------------------------------------------------
cp1 <- do.call(rbind, lapply(ARCH$archetype, function(A) { s <- prim[prim$archetype == A, ]; g <- function(m) s$S[s$model == m][order(s$fold[s$model == m])]
  data.frame(archetype = A, m2_beats_m0 = sum(g("M2") < g("M0")), m1_not_worse = sum(g("M1") <= g("M2")), m3_not_worse = sum(g("M3") <= g("M2")),
             S_M0 = mean(g("M0")), S_M1 = mean(g("M1")), S_M2 = mean(g("M2")), S_M3 = mean(g("M3"))) }))
write.csv(cp1, file.path(RES, "S1.4_cp1.csv"), row.names = FALSE)
for (i in seq_len(nrow(cp1))) vr("CP", "CP-1", cp1$archetype[i], "M2 batte M0 (fold); M1 | M3 non peggiori (fold)",
   sprintf("%d/5; %d | %d; S %.2f %.2f %.2f %.2f", cp1$m2_beats_m0[i], cp1$m1_not_worse[i], cp1$m3_not_worse[i], cp1$S_M0[i], cp1$S_M1[i], cp1$S_M2[i], cp1$S_M3[i]), ">= 4/5",
   pf(cp1$m2_beats_m0[i] >= 4), note = paste(c(if (cp1$m1_not_worse[i] >= 2) "M1 basta", if (cp1$m3_not_worse[i] >= 2) "M3 basta"), collapse = "; "))
vr("CP", "CP-1", "A1-A6", "smentita? (M2 batte M0 >= 4/5 in >= 5/6)", sprintf("%d/6 archetipi", sum(cp1$m2_beats_m0 >= 4)), ">= 5/6", if (sum(cp1$m2_beats_m0 >= 4) >= 5) "PASS" else "FAIL")
# ---- CP-3 sensibilita' al segmentatore -----------------------------------------------------------------------------
for (A in c("A4", "A6")) for (p in c("spaceranger", "stardist_he")) { s <- folds[folds$archetype == A & folds$method == p & folds$model == "M2", ]; ci <- cal[cal$archetype == A & cal$method == p, ]
  vr("CP", "CP-3", A, sprintf("%s: rho beta sigma; B-1 B-2", p), sprintf("%.3f %.3f %.3f; %s %s", ci$rho, ci$beta, ci$sigma, arch_v(s$v_b1), arch_v(s$v_b2)), "riportata", "INFO") }
# ---- previsioni -------------------------------------------------------------------------------------------------
bv <- do.call(rbind, V); bB <- bv[bv$section == "B", ]
m0 <- prim[prim$model == "M0", ]; m0v <- tapply(m0$v_b1, m0$archetype, arch_v)
vr("P", "P1", "A1-A6", "beta (5 ROI) in (0, 1)", sprintf("%d/6: %s", sum(cp$beta > 0 & cp$beta < 1), paste(sprintf("%.2f", cp$beta), collapse = " ")), ">= 5/6", pf(sum(cp$beta > 0 & cp$beta < 1) >= 5))
vr("P", "P2", "A1-A6", "M0 (design) FAIL B-1", sprintf("%d/6 (%s)", sum(m0v == "FAIL"), paste(m0v, collapse = " ")), "6/6", pf(all(m0v == "FAIL")))
b12 <- vapply(ARCH$archetype, function(A) all(bB$outcome[bB$archetype == A & bB$id %in% c("B-1", "B-2")] == "PASS"), TRUE)
vr("P", "P3", "A1-A6", "M2 PASS B-1 e B-2", sprintf("%d/6 (%s)", sum(b12), paste(names(b12)[b12], collapse = " ")), ">= 4/6", pf(sum(b12) >= 4))
cl <- tapply(m2$clamp_all, m2$archetype, mean)
vr("P", "P4", "A1-A6", "clamp > 5 % solo in A4", paste(sprintf("%s %.3f", names(cl), cl), collapse = "; "), "A4 > 0.05, altri <= 0.05", pf(cl[["A4"]] > 0.05 && all(cl[names(cl) != "A4"] <= 0.05)))
vr("P", "P5", "A1-A6", "C8 rileva x2 (>= 4/5 ROI)", sprintf("%d/6", sum(c8$det_x2 >= 4)), ">= 4/6", pf(sum(c8$det_x2 >= 4) >= 4))
out <- do.call(rbind, V); write.csv(out, file.path(RES, "S1.4_verdicts.csv"), row.names = FALSE)
cat("\nverdetti:", nrow(out), "| PASS", sum(out$outcome == "PASS"), "WARN", sum(out$outcome == "WARN"), "FAIL", sum(out$outcome == "FAIL"), "INFO", sum(out$outcome == "INFO"), "\n")
