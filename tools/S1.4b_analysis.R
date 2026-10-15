# tools/S1.4b_analysis.R — S1.4b: verdetti (check C, check B, CP-1, previsioni) e regola di decisione, dai grezzi del banco.
# Definizioni: results/S1.4b/S1.4b_preregistration.md, results/S1.4b/S1.4b_dev_notes.md. Uso: Rscript --vanilla tools/S1.4b_analysis.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
OUT <- "/mnt/micron/geo_spatialtrans/S1.4b"; RES <- "results/S1.4b"
ARCH <- data.frame(archetype = paste0("A", 1:6), primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"))
cls <- function(v, med, mx) if (v <= med) "PASS" else if (v <= mx) "WARN" else "FAIL"
arch_v <- function(o) if (sum(o == "PASS") >= 4) "PASS" else if (sum(o == "FAIL") >= 2) "FAIL" else "WARN"
V <- list(); vr <- function(section, id, archetype, metric, value, threshold, outcome, note = "")
  V[[length(V) + 1]] <<- data.frame(section, id, archetype, metric, value = as.character(value), threshold, outcome, note)
loo <- lapply(seq_len(nrow(ARCH)), function(i) lapply(1:5, function(k) readRDS(file.path(OUT, "loo", sprintf("%s_%s_fold%d.rds", ARCH$archetype[i], ARCH$primary[i], k)))))
names(loo) <- ARCH$archetype
all_f <- unlist(loo, recursive = FALSE)
chk <- do.call(rbind, lapply(all_f, function(o) cbind(archetype = o$archetype, fold = o$fold, o$checks)))
met <- do.call(rbind, lapply(all_f, function(o) cbind(archetype = o$archetype, fold = o$fold, roi = o$roi_test, o$m)))
write.csv(chk, file.path(RES, "S1.4b_checks_raw.csv"), row.names = FALSE); write.csv(met, file.path(RES, "S1.4b_metrics_raw.csv"), row.names = FALSE)

# ---- Check C ----------------------------------------------------------------------------------------------------------
c0 <- chk$c0_match
vr("C", "C-0", "A1-A6", "V-c = S1.4 (s = 0: sims$M2; s = 1..10: seeds)", sprintf("%d/%d", sum(c0), length(c0)), "identico", if (all(c0)) "PASS" else "FAIL")
tst <- read.csv(file.path(RES, "S1.4b_test_none.csv"))
vr("C", "C-L1", "-", "test analitici (test_S1.4b.R)", sprintf("%d/%d", sum(tst$pass), nrow(tst)), "tutti", if (all(tst$pass)) "PASS" else "FAIL")
tm <- vapply(c("MU1", "MU2", "MU3", "MU4"), function(m) { f <- file.path(RES, sprintf("S1.4b_test_%s.csv", m)); if (file.exists(f)) sum(!read.csv(f)$pass) else NA_integer_ }, 1L)
vr("C", "C-L1 mut", "-", "test rieseguiti con MU1..MU4: asserzioni FAIL", paste(sprintf("%s %d", names(tm), tm), collapse = "; "), "> 0 ciascuno", if (all(tm > 0, na.rm = TRUE) && !anyNA(tm)) "PASS" else "FAIL")
rdp <- function(f) { x <- read.csv(f); for (k in c("pass_", "area_pass")) x[[k]] <- as.logical(x[[k]]); x }
arb <- rdp(file.path(RES, "S1.4b_arbiter.csv"))
vr("C", "C-L2", "A1-A6", "arbitro: insiemi PASS (assegnati male | scoperti | doppi; |z|>4 oss/att)",
   sprintf("%d/%d (%d | %d | %d; %d / %.2f)", sum(arb$pass_), nrow(arb), sum(arb$misassigned), sum(arb$uncovered), sum(arb$double), sum(arb$n_absz_gt4), sum(arb$exp_absz_gt4)),
   "30/30", if (all(arb$pass_)) "PASS" else "FAIL")
vr("C", "C-L3", "A1-A6", "partizione: max rel |sum aree - regione|; max sovrapposizione (µm²)",
   sprintf("%.2e; %.2e", max(chk$cl3_rel), max(chk$cl3_max_overlap)), "<= 1e-9; <= 1e-9", if (max(chk$cl3_rel) <= 1e-9 && max(chk$cl3_max_overlap) <= 1e-9) "PASS" else "FAIL")
vr("C", "C3", "A1-A6", "V-a: nuclei fuori dalla cella (geos | sf), 30 fold x 11 semi", sprintf("%d | %d su %d nuclei", sum(chk$c3_viol), sum(chk$c3_sf_viol), sum(sapply(all_f, function(o) o$n)) * 11L),
   "0", if (sum(chk$c3_viol) == 0 && sum(chk$c3_sf_viol) == 0) "PASS" else "FAIL")
vr("C", "C-NO", "A1-A6", "V-a/V-b: coppie di nuclei che si intersecano", sum(chk$cno_pairs), "0", if (sum(chk$cno_pairs) == 0) "PASS" else "FAIL")
vr("C", "C3 V-b", "A1-A6", "V-b: nuclei fuori dal territorio di Voronoi (frazione, s = 1..10)",
   sprintf("%.4f", sum(chk$c3_b_viol[chk$seed_idx > 0]) / (10 * sum(sapply(all_f, function(o) o$n)))), "riportato", "INFO")
vr("C", "C-hidden", "A1-A6", "V-a: celle vuote | MULTIPOLYGON | ritagliate (totale su 30 x 11)", sprintf("%d | %d | %d", sum(chk$n_hidden), sum(chk$n_multi), sum(chk$n_clipped)), "riportato", "INFO")
# C-mut
mut <- do.call(rbind, lapply(list.files(file.path(OUT, "mut"), full.names = TRUE), function(f) { o <- readRDS(f); cbind(archetype = o$archetype, mutant = o$mutant, o$checks) }))
am <- rdp(file.path(RES, "S1.4b_arbiter_mut.csv")); am$mutant <- sub(".*_(MU[0-9])$", "\\1", am$set)
a5 <- rdp(file.path(RES, "S1.4b_arbiter_MU5.csv"))
det <- c(MU1 = any(mut$c3_viol[mut$mutant == "MU1"] > 0) || any(!am$pass_[am$mutant == "MU1"]),
         MU2 = any(mut$c3_viol[mut$mutant == "MU2"] > 0) || any(!am$pass_[am$mutant == "MU2"]),
         MU3 = any(mut$cno_pairs[mut$mutant == "MU3"] > 0),
         MU4 = any(mut$cl3_rel[mut$mutant == "MU4"] > 1e-9 | mut$cl3_max_overlap[mut$mutant == "MU4"] > 1e-9) || any(!am$pass_[am$mutant == "MU4"]),
         MU5 = any(!a5$pass_))
detail <- sprintf("MU1 C3 %d, arb %d/6 FAIL; MU2 C3 %d, arb %d/6; MU3 C-NO %d; MU4 C-L3 %d/6, arb %d/6; MU5 arb %d/30",
  sum(mut$c3_viol[mut$mutant == "MU1"]), sum(!am$pass_[am$mutant == "MU1"]), sum(mut$c3_viol[mut$mutant == "MU2"]), sum(!am$pass_[am$mutant == "MU2"]),
  sum(mut$cno_pairs[mut$mutant == "MU3"]), sum(mut$cl3_rel[mut$mutant == "MU4"] > 1e-9 | mut$cl3_max_overlap[mut$mutant == "MU4"] > 1e-9), sum(!am$pass_[am$mutant == "MU4"]), sum(!a5$pass_))
vr("C", "C-mut", "A1-A6", "mutanti rilevati", sprintf("%d/5 (%s)", sum(det), detail), "5/5", if (all(det)) "PASS" else "FAIL")
pl <- readLines(file.path(OUT, "perf.log"))
pf <- grep("^PERF", pl, value = TRUE); rss <- as.numeric(sub(".*: ", "", grep("Maximum resident", pl, value = TRUE))) / 1024^2
tot <- as.numeric(sub(".*total_s=([0-9.]+).*", "\\1", pf))
vr("C", "C-perf", "S", "Laguerre 449 536 celle (nuclei + diagramma; RSS del processo)", sprintf("%.1f s, %.2f GB; %s", tot, rss, sub("^PERF laguerre ", "", pf)),
   "<= 300 s, <= 8 GB (solo WARN)", if (tot <= 300 && rss <= 8) "PASS" else "WARN")

# ---- Check B ----------------------------------------------------------------------------------------------------------
Bm <- c(`B-1` = "ks_nc", `B-2` = "ks_an", `B-3` = "d_sp", `B-7` = "b7")
fold_tab <- do.call(rbind, lapply(all_f, function(o) { m <- o$m[o$m$seed_idx >= 1, ]
  do.call(rbind, lapply(c("V-c", "V-a"), function(v) { x <- m[m$variant == v, ]
    d <- c(ks_nc = median(x$ks_nc), ks_an = median(x$ks_an), d_sp = median(x$d_sp), b7 = median(abs(x$b7 - o$b7_real)),
           trunc_all = median(x$trunc_all), loss_all = median(x$loss_all), b7_sim = median(x$b7))
    y <- list(ks_nc = o$yard$ks_nc, ks_an = o$yard$ks_an, d_sp = o$yard$d_sp, b7 = o$yard_b7)
    data.frame(archetype = o$archetype, fold = o$fold, roi = o$roi_test, variant = v, t(d), b7_real = o$b7_real,
               v_b1 = cls(d[["ks_nc"]], median(y$ks_nc), max(y$ks_nc)), v_b2 = cls(d[["ks_an"]], median(y$ks_an), max(y$ks_an)),
               v_b3 = cls(d[["d_sp"]], median(y$d_sp), max(y$d_sp)), v_b7 = cls(d[["b7"]], median(y$b7), max(y$b7)),
               y_b7_med = median(y$b7), y_b7_max = max(y$b7)) })) }))
fold_tab$variant[fold_tab$variant == "V-a"] <- "V-a/V-b"
write.csv(fold_tab, file.path(RES, "S1.4b_folds.csv"), row.names = FALSE)
arch_tab <- do.call(rbind, lapply(split(fold_tab, list(fold_tab$archetype, fold_tab$variant)), function(z)
  data.frame(archetype = z$archetype[1], variant = z$variant[1], B1 = arch_v(z$v_b1), B2 = arch_v(z$v_b2), B3 = arch_v(z$v_b3), B7 = arch_v(z$v_b7),
             B4 = mean(z$trunc_all), B4_v = if (mean(z$trunc_all) <= 0.05) "PASS" else "FAIL", loss = mean(z$loss_all),
             b7_sim = mean(z$b7_sim), b7_real = mean(z$b7_real),
             f_b1 = paste(substr(z$v_b1, 1, 1), collapse = ""), f_b2 = paste(substr(z$v_b2, 1, 1), collapse = ""),
             f_b3 = paste(substr(z$v_b3, 1, 1), collapse = ""), f_b7 = paste(substr(z$v_b7, 1, 1), collapse = ""))))
arch_tab <- arch_tab[order(arch_tab$variant, arch_tab$archetype), ]
write.csv(arch_tab, file.path(RES, "S1.4b_archetypes.csv"), row.names = FALSE)
for (i in seq_len(nrow(arch_tab))) { a <- arch_tab[i, ]; z <- fold_tab[fold_tab$archetype == a$archetype & fold_tab$variant == a$variant, ]
  vr("B", "B-1", a$archetype, sprintf("%s: KS N/C (fold %s)", a$variant, a$f_b1), paste(signif(z$ks_nc, 2), collapse = " "), "metro inter-ROI", a$B1)
  vr("B", "B-2", a$archetype, sprintf("%s: KS A_N (fold %s)", a$variant, a$f_b2), paste(signif(z$ks_an, 2), collapse = " "), "metro inter-ROI", a$B2)
  vr("B", "B-3", a$archetype, sprintf("%s: |dSpearman| (fold %s)", a$variant, a$f_b3), paste(signif(z$d_sp, 2), collapse = " "), "metro inter-ROI", a$B3)
  vr("B", "B-7", a$archetype, sprintf("%s: frazione fuori dal Voronoi, sim | reale (fold %s)", a$variant, a$f_b7), sprintf("%.3f | %.3f", a$b7_sim, a$b7_real), "metro inter-ROI", a$B7)
  vr("B", "B-4", a$archetype, sprintf("%s: frazione troncata (perdita media di area)", a$variant), sprintf("%.4f (%.4f)", a$B4, a$loss), "<= 0.05", a$B4_v) }
# B-8, B-9 descrittive
b9 <- do.call(rbind, lapply(all_f, function(o) { x <- o$checks[o$checks$seed_idx >= 1, ]
  data.frame(archetype = o$archetype, fold = o$fold, dat_med = median(x$dat_med), dat_real_med = median(x$dat_real_med), cv_vor = median(x$cv_vor),
             cv_lag = median(x$cv_lag), cv_lag_real = median(x$cv_lag_real), ks_at_lag_real = median(x$ks_at_lag_real),
             trunc_pair = median(x$trunc_pair), trunc_reg = median(x$trunc_reg), b8_circle = o$b8[["circle"]], b8_rule = o$b8[["rule"]]) }))
write.csv(b9, file.path(RES, "S1.4b_b8_b9.csv"), row.names = FALSE)
for (A in ARCH$archetype) { z <- b9[b9$archetype == A, ]
  vr("B", "B-8", A, "cerchi reali sovrapposti: r_i + r_j > d | regola 0.99 (media ROI)", sprintf("%.3f | %.3f", mean(z$b8_circle), mean(z$b8_rule)), "descrittiva", "INFO")
  vr("B", "B-9", A, "|dA_T|/A_T mediana: Laguerre sim | Laguerre reale; CV Vor | Lag sim | Lag reale",
     sprintf("%.3f | %.3f; %.2f | %.2f | %.2f", median(z$dat_med), median(z$dat_real_med), median(z$cv_vor), median(z$cv_lag), median(z$cv_lag_real)), "descrittiva", "INFO")
  vr("B", "B-4 det", A, "V-a/V-b: troncati per coppia | per bordo della regione", sprintf("%.4f | %.4f", mean(z$trunc_pair), mean(z$trunc_reg)), "descrittiva", "INFO") }

# ---- CP-1 -------------------------------------------------------------------------------------------------------------
cp <- do.call(rbind, lapply(list.files(file.path(OUT, "cp1"), pattern = "csv$", full.names = TRUE), read.csv))
write.csv(cp, file.path(RES, "S1.4b_cp1_raw.csv"), row.names = FALSE)
r3c <- cp[cp$cond == "check_R3", ]
vr("C", "C-cp1", "A1-A6", "CP-1 con lambda = 0 = frac_out di R3 (max |diff|; frazione identica)", sprintf("%.2e; %.4f", max(r3c$m1), min(r3c$m2)), "0; 1", if (max(r3c$m1) == 0 && min(r3c$m2) == 1) "PASS" else "FAIL")
cp <- cp[cp$cond != "check_R3", ]
roi_cp <- do.call(rbind, lapply(split(cp, paste(cp$archetype, cp$roi)), function(z) { pm <- z[grepl("^perm", z$cond), ]
  g <- function(cn, col) z[[col]][z$cond == cn]
  data.frame(archetype = z$archetype[1], roi = z$roi[1], n = g("vor", "n_nuc"), m1_vor = g("vor", "m1"), m1_lag = g("lag", "m1"), m1_perm_min = min(pm$m1), m1_perm_med = median(pm$m1),
             m2_vor = g("vor", "m2"), m2_lag = g("lag", "m2"), m2_perm_min = min(pm$m2), m2_perm_med = median(pm$m2),
             m1_l05 = g("lam0.5", "m1"), m1_l15 = g("lam1.5", "m1"), m1_l2 = g("lam2", "m1"),
             m2_l05 = g("lam0.5", "m2"), m2_l15 = g("lam1.5", "m2"), m2_l2 = g("lam2", "m2"),
             beats = g("lag", "m1") < min(pm$m1) && g("lag", "m2") < min(pm$m2)) }))
write.csv(roi_cp, file.path(RES, "S1.4b_cp1_roi.csv"), row.names = FALSE)
cpa <- do.call(rbind, lapply(split(roi_cp, roi_cp$archetype), function(z) data.frame(archetype = z$archetype[1], n_beats = sum(z$beats),
  beats = sum(z$beats) >= 4, m1_vor = mean(z$m1_vor), m1_lag = mean(z$m1_lag), m1_perm = mean(z$m1_perm_med),
  m2_vor = mean(z$m2_vor), m2_lag = mean(z$m2_lag), m2_perm = mean(z$m2_perm_med), red_m2 = 1 - mean(z$m2_lag) / mean(z$m2_vor), red_m1 = 1 - mean(z$m1_lag) / mean(z$m1_vor))))
write.csv(cpa, file.path(RES, "S1.4b_cp1_archetype.csv"), row.names = FALSE)
for (i in seq_len(nrow(cpa))) { a <- cpa[i, ]
  vr("CP", "CP-1", a$archetype, "Laguerre batte 20/20 permutazioni (ROI); m1 vor | lag | perm; m2 vor | lag | perm",
     sprintf("%d/5; %.4f | %.4f | %.4f; %.3f | %.3f | %.3f", a$n_beats, a$m1_vor, a$m1_lag, a$m1_perm, a$m2_vor, a$m2_lag, a$m2_perm), ">= 4/5 ROI", if (a$beats) "PASS" else "FAIL") }
refuted <- sum(cpa$beats) < 5
vr("CP", "CP-1", "A1-A6", "(a) smentita? archetipi in cui Laguerre batte le permutazioni", sprintf("%d/6", sum(cpa$beats)), ">= 5/6 per non smentire", if (refuted) "FAIL" else "PASS",
   if (refuted) "(a) smentita" else "(a) non smentita")

# ---- Previsioni -------------------------------------------------------------------------------------------------------
g <- function(v, col) setNames(arch_tab[[col]][arch_tab$variant == v], arch_tab$archetype[arch_tab$variant == v])
p1 <- all(cpa$red_m2[cpa$archetype %in% c("A1", "A2", "A3", "A4")] >= 0.5) && sum(cpa$beats) == 6
vr("P", "P1", "A1-A6", "riduzione m2 Laguerre vs Voronoi in A1-A4 (>= 0.5); CP-1 6/6",
   sprintf("%s; %d/6", paste(sprintf("%s %.2f", cpa$archetype, cpa$red_m2), collapse = " "), sum(cpa$beats)), ">= 0.5 in A1-A4 e 6/6", if (p1) "PASS" else "FAIL")
b4c <- g("V-c", "B4"); b4a <- g("V-a/V-b", "B4"); d4 <- c("A1", "A2", "A3", "A4")
vr("P", "P2", "A1-A4", "B-4 V-a/b <= 0.5 x V-c", paste(sprintf("%s %.3f/%.3f", d4, b4a[d4], b4c[d4]), collapse = "; "), "tutti", if (all(b4a[d4] <= 0.5 * b4c[d4])) "PASS" else "FAIL")
b1a <- g("V-a/V-b", "B1"); b2a <- g("V-a/V-b", "B2")
vr("P", "P3", "A1-A6", "V-a/b: archetipi con B-1 e B-2 PASS (separati: B-1 | B-2)", sprintf("%d/6 (%d | %d)", sum(b1a == "PASS" & b2a == "PASS"), sum(b1a == "PASS"), sum(b2a == "PASS")), ">= 4/6",
   if (sum(b1a == "PASS" & b2a == "PASS") >= 4) "PASS" else "FAIL")
b7a <- g("V-a/V-b", "B7"); b7c <- g("V-c", "B7")
vr("P", "P4", "A1-A6", "B-7: V-a/b PASS o WARN; V-c FAIL in A1-A5", sprintf("V-a/b %d/6; V-c FAIL %d/5", sum(b7a != "FAIL"), sum(b7c[paste0("A", 1:5)] == "FAIL")), ">= 5/6; 5/5",
   if (sum(b7a != "FAIL") >= 5 && all(b7c[paste0("A", 1:5)] == "FAIL")) "PASS" else "FAIL")
dm <- tapply(b9$dat_med, b9$archetype, median)
vr("P", "P5", "A1-A4", "B-9 |dA_T|/A_T mediana <= 0.10", paste(sprintf("%s %.3f", d4, dm[d4]), collapse = " "), "tutti", if (all(dm[d4] <= 0.10)) "PASS" else "FAIL")
vr("P", "P6", "A5", "B-1 A5 non PASS per tutte le opzioni", sprintf("V-c %s; V-a/b %s", g("V-c", "B1")[["A5"]], b1a[["A5"]]), "non PASS", if (g("V-c", "B1")[["A5"]] != "PASS" && b1a[["A5"]] != "PASS") "PASS" else "FAIL")

# ---- Regola di decisione -----------------------------------------------------------------------------------------------
Cv <- do.call(rbind, V); cc <- Cv[Cv$section == "C", ]
elig_a <- all(cc$outcome[cc$id %in% c("C-L1", "C-L1 mut", "C-L2", "C-L3", "C3", "C-NO", "C-mut")] == "PASS")
elig_b <- cc$outcome[cc$id == "C-NO"] == "PASS"; elig_c <- cc$outcome[cc$id == "C-0"] == "PASS"
score <- function(v) sum(unlist(arch_tab[arch_tab$variant == v, c("B1", "B2", "B3", "B7")]) == "PASS")
s_ab <- score("V-a/V-b"); s_c <- score("V-c")
step2 <- if (!(elig_a || elig_b)) "c" else if (!elig_c) "ab" else if (s_ab > s_c) "ab" else "c"
choice <- if (step2 == "c") "(c) contenimento con troncamento dichiarato" else if (elig_a && !refuted) "(a) Laguerre" else if (elig_b) "(b) nucleo non contenuto" else "(a) Laguerre"
vr("D", "Regola", "-", "eleggibili a | b | c; PASS {B-1,B-2,B-3,B-7}: {a,b} vs c; CP-1",
   sprintf("%s | %s | %s; %d vs %d; %s", elig_a, elig_b, elig_c, s_ab, s_c, if (refuted) "a smentita" else "a non smentita"), "regola pre-registrata", choice)
out <- do.call(rbind, V); write.csv(out, file.path(RES, "S1.4b_verdicts.csv"), row.names = FALSE)
cat("Decisione per regola:", choice, "\n"); print(table(out$section, out$outcome))
