# tools/S1.3b_analysis.R — S1.3b (BL-069): verdetti e decisione. Pre-registrazione results/S1.3b/S1.3b_preregistration.md (d96fea1).
# Committato prima del run ufficiale (BL-059). Uso: Rscript --vanilla tools/S1.3b_analysis.R [outdir] → results/S1.3b/S1.3b_*.csv
# Definizioni operative dichiarate prima del run (2026-10-10):
#  (a) C-2 per V0/V0g: tau_ij = 1e-6 * (P_i + P_j) (stessa scala di V1/V2); con g = 0 ogni residuo di 1e-13 µm² conterebbe.
#  (b) C-0: c1 e c2 sono già relativi all'area della regione; «rel. ≤ 1e-12» = |Δ| ≤ 1e-12.
#  (c) CP-2, criterio di FAIL di V0: |ov| o symd > 1e-9·A + √2·δ·ΣP (la soglia di S1.3 più lo spostamento ammesso dalla perturbazione).
#  (d) B-2: verdetti di S1.3 ricalcolati con lo script di S1.3 sulle aree della variante (tools/S1.3b_b2.R); confronto su tutte le
#      righe B (outcome, value, confirmed); le «62 previsioni» = righe B con confirmed non NA.
#  Nota (sviluppo, prima del run): B-1 a 1e-6 relativo è incompatibile con tau_i = g·P_i per celle piccole (errore di
#  pre-registrazione); la soglia resta quella pre-registrata, l'esito si riporta com'è. B-1 non entra nella regola di eleggibilità.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
OUTB <- commandArgs(TRUE)[1]; if (is.na(OUTB)) OUTB <- "/mnt/micron/geo_spatialtrans/S1.3b"
RES <- "results/S1.3b"; dir.create(RES, showWarnings = FALSE)
V <- list(); vr <- function(section, id, variant, metric, value, threshold, outcome, prediction = NA, confirmed = NA)
  V[[length(V) + 1]] <<- data.frame(section, id, variant, metric, value = as.character(value), threshold, outcome,
                                    prediction = as.character(prediction), confirmed = as.character(confirmed))
`%||%` <- function(a, b) if (is.null(a)) b else a
pf <- function(ok) ifelse(ok, "PASS", "FAIL"); fe <- function(x) formatC(x, format = "e", digits = 2)
VS <- c("V1", "V2", "V3")

# ---- nulli --------------------------------------------------------------------------------------------------------
nf <- list.files(file.path(OUTB, "null"), pattern = "\\.rds$", full.names = TRUE); stopifnot(length(nf) == 1800)
nul <- lapply(nf, readRDS)
key <- vapply(nul, function(o) paste(o$archetype, o$roi_id, o$model, o$rep), "")
tab <- do.call(rbind, lapply(nul, function(o) do.call(rbind, lapply(names(o$res), function(v) { m <- o$res[[v]]
  pr <- m$pairs; g2 <- if (v %in% c("V0", "V0g")) 1e-6 else m$gchk
  data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep, n = o$n, A = o$A, sumP = o$sumP, Preg = o$Preg,
             reg_holes = o$reg[["n_holes"]], reg_vertices = o$reg[["n_vertices"]], variant = v,
             ov = m$ov, symd = m$symd, tau = m$tau, n_pairs = m$n_pairs, pairs_over = sum(pr$area > g2 * pr$psum),
             pair_max = if (nrow(pr)) max(pr$area) else 0,
             cov_anom = if (is.null(m$cov)) NA else m$cov[["n_anom"]], cov_far = if (is.null(m$cov)) NA else m$cov[["n_anom_far"]],
             cov_dmax = if (is.null(m$cov)) NA else m$cov[["dmax"]],
             n_invalid = m$n_invalid %||% NA, n_empty = m$n_empty %||% NA, n_out = m$n,
             types_same = if (is.null(m$types)) NA else identical(as.vector(m$types[names(o$res$V0$types)]), as.vector(o$res$V0$types)) && length(m$types) == length(o$res$V0$types),
             c5_over = m$c5[["n_over"]], c5_ratio = m$c5[["max_ratio"]], b1_max_rel = m$b1[["max_rel"]], b1_over = m$b1[["n_over"]],
             flips = m$interior_flips,
             s_median_area = m$summary[["median_area"]], s_cv = m$summary[["cv"]], s_median_eq_r = m$summary[["median_eq_r"]], s_frac_interior = m$summary[["frac_interior"]],
             v3_pairs_diff = if (is.null(m$info)) NA else m$info[["n_pairs_cross_diff"]], v3_tj = if (is.null(m$info)) NA else m$info[["n_tjunction"]],
             v3_dropped = if (is.null(m$info)) NA else m$info[["n_rings_dropped"]], v3_dmax = if (is.null(m$info)) NA else m$info[["d_max"]],
             t_total = sum(m$time), t_union = m$time[["union"]], t_build = m$time[["build"]],
             c0_c1 = abs(o$c0[["c1"]] - o$old_c123[["c1"]]), c0_c2 = abs(o$c0[["c2_overlap"]] - o$old_c123[["c2_overlap"]]),
             c0_n = o$n == o$old_summary[["n"]], c0_c2_old = o$old_c123[["c2_overlap"]], c0_sd_old = o$old_c123[["c2_symdiff"]])
}))))
write.csv(tab, file.path(RES, "S1.3b_null_table.csv"), row.names = FALSE)
t0 <- tab[tab$variant == "V0", ]
# C-0
c0ok <- t0$c0_n & t0$c0_c1 <= 1e-12 & t0$c0_c2 <= 1e-12
fail_new <- t0$ov / t0$A > 1e-9 | t0$symd / t0$A > 1e-9; fail_old <- t0$c0_c2_old > 1e-9 | t0$c0_sd_old > 1e-9
vr("C", "C-0", "V0", "nulli riprodotti (n, c1, c2 ≤ 1e-12)", sprintf("%d/1800", sum(c0ok)), "1800/1800", pf(all(c0ok)), "1800/1800", all(c0ok))
vr("C", "C-0", "V0", "nulli con C-2 > 1e-9 (stessi di S1.3)", sprintf("%d (S1.3: %d; identici: %s)", sum(fail_new), sum(fail_old), identical(fail_new, fail_old)), "15, stessi", pf(sum(fail_new) == 15 && identical(fail_new, fail_old)), "15", sum(fail_new) == 15 && identical(fail_new, fail_old))
# C-1
for (v in c("V0", "V0g")) { d <- tab[tab$variant == v, ]
  vr("C", "C-1", v, "nulli con |ov| o symd > 1e-9·A", sum(d$ov / d$A > 1e-9 | d$symd / d$A > 1e-9), "descrittivo", "INFO") }
c1 <- sapply(VS, function(v) { d <- tab[tab$variant == v, ]; sum(abs(d$ov) > d$tau | d$symd > d$tau) })
for (v in VS) { d <- tab[tab$variant == v, ]
  vr("C", "C-1", v, "nulli oltre tau; max |ov|/tau, max symd/tau", sprintf("%d; %s; %s", c1[[v]], fe(max(abs(d$ov) / d$tau)), fe(max(d$symd / d$tau))), "0", pf(c1[[v]] == 0), "0", c1[[v]] == 0) }
# C-2
for (v in c("V0", "V0g")) { d <- tab[tab$variant == v, ]; k <- which(d$pairs_over > 0)
  vr("C", "C-2", v, "coppie oltre 1e-6·(Pi+Pj) (nulli)", sprintf("%d (%s)", sum(d$pairs_over), paste(paste(d$archetype[k], d$roi_id[k], d$model[k], d$rep[k]), collapse = "; ")),
     "1 (A3 r3 CSR 3)", if (v == "V0") pf(sum(d$pairs_over) == 1 && identical(paste(d$archetype[k], d$roi_id[k], d$model[k], d$rep[k]), "A3 r3 CSR 3")) else "INFO",
     if (v == "V0") "1" else NA, if (v == "V0") sum(d$pairs_over) == 1 else NA) }
c2 <- sapply(VS, function(v) sum(tab$pairs_over[tab$variant == v]))
for (v in VS) vr("C", "C-2", v, "coppie oltre tau_ij", c2[[v]], "0", pf(c2[[v]] == 0), "0", c2[[v]] == 0)
# C-3
c3 <- sapply(c("V0", VS), function(v) { d <- tab[tab$variant == v, ]; c(far = sum(d$cov_far), anom = sum(d$cov_anom), dmax = max(d$cov_dmax)) })
vr("C", "C-3", "V0", "anomalie di copertura (lontane dal bordo / totale)", sprintf("%d / %d", c3["far", "V0"], c3["anom", "V0"]), "descrittivo", "INFO")
for (v in VS) vr("C", "C-3", v, "anomalie lontane > 1e-5 µm / totale (max distanza)", sprintf("%d / %d (%s)", c3["far", v], c3["anom", v], fe(c3["dmax", v])), "0 / <= 50",
                 pf(c3["far", v] == 0 && c3["anom", v] <= 50), "0 / <= 50", c3["far", v] == 0 && c3["anom", v] <= 50)
# C-4
c4 <- sapply(VS, function(v) { d <- tab[tab$variant == v, ]; c(inv = sum(d$n_invalid), empty = sum(d$n_empty), n_bad = sum(d$n_out != d$n), types = sum(!d$types_same)) })
for (v in VS) { ok <- c4["empty", v] == 0 && c4["n_bad", v] == 0 && c4["types", v] == 0
  out <- if (v == "V1") { if (!ok) "FAIL" else if (c4["inv", v] > 0) "WARN" else "PASS" } else pf(ok && c4["inv", v] == 0)
  vr("C", "C-4", v, "invalide; vuote; n diverso; tipi cambiati", paste(c4[, v], collapse = "; "), if (v == "V1") "0 (invalide: WARN)" else "0; 0; 0; 0", out,
     if (v == "V1") "conteggio" else "0", if (v == "V1") NA else ok && c4["inv", v] == 0) }
# C-5
c5 <- sapply(VS, function(v) sum(tab$c5_over[tab$variant == v]))
for (v in VS) vr("C", "C-5", v, "celle con |Δa| > g·P_i (max rapporto)", sprintf("%d (%s)", c5[[v]], fe(max(tab$c5_ratio[tab$variant == v]))), "0", pf(c5[[v]] == 0), "0", c5[[v]] == 0)
# C-unit
ut <- read.csv(file.path(RES, "S1.3b_test.csv"))
cu <- sapply(VS, function(v) { k <- ut$variant %in% c(v, "-") & ut$status != "INFO"; c(pass = sum(ut$status[k] == "PASS"), n = sum(k)) })
for (v in VS) vr("C", "C-unit", v, "asserzioni PASS", sprintf("%d/%d", cu["pass", v], cu["n", v]), "100 %", pf(cu["pass", v] == cu["n", v]), "100 %", cu["pass", v] == cu["n", v])
# C-6
gf <- list.files(file.path(OUTB, "grid"), pattern = "\\.rds$", full.names = TRUE); gr <- lapply(gf, readRDS); stopifnot(length(gr) == 30)
gt <- do.call(rbind, lapply(gr, function(o) do.call(rbind, lapply(c("V0", VS), function(v) data.frame(set = o$set, archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep,
  n = o$n, n_bins = o$n_bins, variant = v, t(o[[v]]))))))
write.csv(gt, file.path(RES, "S1.3b_grid.csv"), row.names = FALSE)
for (v in c("V0", VS)) { d <- gt[gt$variant == v, ]; nv <- sum(d$cell_viol) + sum(d$bin_viol)
  vr("C", "C-6", v, "violazioni celle + bin (30 nulli); max |Δ| cella, bin", sprintf("%d + %d; %s, %s", sum(d$cell_viol), sum(d$bin_viol), fe(max(d$cell_max)), fe(max(d$bin_max))),
     if (v == "V0") "nessuna previsione" else "0", if (v == "V0") "INFO" else pf(nv == 0), if (v == "V0") NA else "0", if (v == "V0") NA else nv == 0) }
c6 <- sapply(VS, function(v) { d <- gt[gt$variant == v, ]; sum(d$cell_viol) + sum(d$bin_viol) })
# mutanti
mu <- readRDS(file.path(OUTB, "mutants", "mutants.rds")); write.csv(mu, file.path(RES, "S1.3b_mutants.csv"), row.names = FALSE)
md <- tapply(mu$n_fail, mu$mutant, function(x) any(x > 0))
vr("C", "mutanti", "-", "rilevati (M1 C-5, M2 C-1/C-2, M3 C-1/C-2)", sprintf("%d/3 (%s)", sum(md), paste(names(md), md, collapse = ", ")), "3/3", pf(sum(md) == 3), "3/3", sum(md) == 3)

# ---- controprove --------------------------------------------------------------------------------------------------
cp <- readRDS(file.path(OUTB, "cp1", "cp1.rds")); k1 <- cp$case1; k2 <- cp$case2
h1 <- k1$d0_info[["n_pairs_cross_diff"]] + k1$d0_info[["n_tjunction"]]
vr("CP", "D-0", "766/1987", "vertici condivisi diversi; giunzioni a T; identici; max d", sprintf("%d; %d; %d; %s", k1$d0_info[["n_pairs_cross_diff"]], k1$d0_info[["n_tjunction"]],
   k1$d0_info[["n_pairs_cross_identical"]], fe(k1$d0_info[["d_max"]])), ">= 1 (diversi o T)", pf(h1 >= 1), ">= 1", h1 >= 1)
vr("CP", "D-0", "A5 r5 RSA 15", "vertici condivisi diversi; giunzioni a T; identici; sul bordo regione", sprintf("%d; %d; %d; %d/%d", k2$d0_info[["n_pairs_cross_diff"]], k2$d0_info[["n_tjunction"]],
   k2$d0_info[["n_pairs_cross_identical"]], sum(k2$diff_on_region_boundary), length(k2$diff_on_region_boundary)), "descrittivo", "INFO")
vr("CP", "CP-1a", "766/1987", "area intersezione dopo annodamento (V0)", sprintf("%s (V0 %.4f)", fe(k1$a_cp1a), k1$a_v0), "<= 1e-6", pf(k1$a_cp1a <= 1e-6), "<= 1e-6", k1$a_cp1a <= 1e-6)
npl <- sum(abs(k1$placebo$area - k1$a_v0) <= 0.01)
vr("CP", "CP-1b", "766/1987", "placebo con errore intatto (|area − V0| <= 0.01)", sprintf("%d/20 (spost. %d, ins. %d+%d)", npl, k1$placebo$moved[1], k1$placebo$ins1[1], k1$placebo$ins2[1]), ">= 19/20", pf(npl >= 19), ">= 19/20", npl >= 19)
ok1c <- k2$cp1c_rel_ov <= 1e-9 && k2$cp1c_rel_symd <= 1e-9
vr("CP", "CP-1c", "A5 r5 RSA 15", "unione dopo V3: |ov|/A, symd/A (V0 deficit µm²)", sprintf("%s, %s (%.1f)", fe(k2$cp1c_rel_ov), fe(k2$cp1c_rel_symd), k2$union_deficit_v0), "<= 1e-9", pf(ok1c), "<= 1e-9", ok1c)
cf <- list.files(file.path(OUTB, "cp2"), pattern = "\\.rds$", full.names = TRUE); c2r <- lapply(cf, readRDS); stopifnot(length(c2r) == 60)
ct <- do.call(rbind, lapply(c2r, function(o) { if (is.null(o$rows$error)) o$rows$error <- NA_character_
  cbind(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep, A = o$A, sumP = o$sumP, o$rows) }))
ct$fail <- ifelse(ct$variant == "V0", abs(ct$ov) > 1e-9 * ct$A + sqrt(2) * ct$delta * ct$sumP | ct$symd > 1e-9 * ct$A + sqrt(2) * ct$delta * ct$sumP,
                  abs(ct$ov) > ct$tau | ct$symd > ct$tau)
if (is.null(ct$error)) ct$error <- NA_character_
ct$fail[is.na(ct$fail) | !is.na(ct$error)] <- TRUE        # POST HOC dichiarato (2026-10-10, prima di vedere la variante): eccezione GEOS = FAIL
write.csv(ct, file.path(RES, "S1.3b_cp2.csv"), row.names = FALSE)
for (v in c("V0", VS)) vr("CP", "CP-2 eccezioni", v, "eccezioni GEOS (per δ = 1e-12 / 1e-10 / 1e-8)",
  paste(tapply(!is.na(ct$error[ct$variant == v]), ct$delta[ct$variant == v], sum), collapse = " / "), "descrittivo (contano come FAIL)", "INFO")
r0 <- tapply(ct$fail[ct$variant == "V0"], ct$delta[ct$variant == "V0"], mean)
vr("CP", "CP-2", "V0", "frequenza FAIL per δ = 1e-12 / 1e-10 / 1e-8 (basale 0.0083)", paste(sprintf("%.3f", r0), collapse = " / "), "> 0.0083 ad almeno un δ", pf(any(r0 > 15 / 1800)), "> basale", any(r0 > 15 / 1800))
cp2f <- sapply(VS, function(v) sum(ct$fail[ct$variant == v]))
for (v in VS) vr("CP", "CP-2", v, "FAIL su 60 x 3 δ (per δ)", sprintf("%d (%s)", cp2f[[v]], paste(tapply(ct$fail[ct$variant == v], ct$delta[ct$variant == v], sum), collapse = "/")), "0", pf(cp2f[[v]] == 0), "0", cp2f[[v]] == 0)
vr("CP", "CP-2", "-", "input perturbati invalidi (per δ)", paste(tapply(ct$n_invalid_input[ct$variant == "V0"], ct$delta[ct$variant == "V0"], sum), collapse = "/"), "descrittivo", "INFO")
# CP-3 descrittivo
t0$fail <- fail_new; v3 <- tab[tab$variant == "V3", ]
cp3 <- data.frame(feature = c("n celle", "buchi della regione", "vertici della regione", "V3 vertici condivisi diversi", "V3 giunzioni a T"),
  median_fail = c(median(t0$n[t0$fail]), median(t0$reg_holes[t0$fail]), median(t0$reg_vertices[t0$fail]), median(v3$v3_pairs_diff[fail_new]), median(v3$v3_tj[fail_new])),
  median_pass = c(median(t0$n[!t0$fail]), median(t0$reg_holes[!t0$fail]), median(t0$reg_vertices[!t0$fail]), median(v3$v3_pairs_diff[!fail_new]), median(v3$v3_tj[!fail_new])),
  p_wilcox = c(wilcox.test(t0$n ~ t0$fail)$p.value, wilcox.test(t0$reg_holes ~ t0$fail)$p.value, wilcox.test(t0$reg_vertices ~ t0$fail)$p.value,
               wilcox.test(v3$v3_pairs_diff ~ fail_new)$p.value, wilcox.test(v3$v3_tj ~ fail_new)$p.value))
write.csv(cp3, file.path(RES, "S1.3b_cp3.csv"), row.names = FALSE)
for (i in seq_len(nrow(cp3))) vr("CP", "CP-3", "V0", paste("FAIL vs PASS:", cp3$feature[i]), sprintf("%g vs %g (p %.2g)", cp3$median_fail[i], cp3$median_pass[i], cp3$p_wilcox[i]), "descrittivo", "INFO")

# ---- check B ------------------------------------------------------------------------------------------------------
for (v in VS) { d <- tab[tab$variant == v, ]
  vr("B", "B-1", v, "nulli con celle (a >= 1 µm²) oltre 1e-6 relativo; max", sprintf("%d (%d celle); %s", sum(d$b1_over > 0), sum(d$b1_over), fe(max(d$b1_max_rel))), "0", pf(sum(d$b1_over) == 0), "0", sum(d$b1_over) == 0) }
for (v in VS) { d <- tab[tab$variant == v, ]; d0 <- t0[match(paste(d$archetype, d$roi_id, d$model, d$rep), paste(t0$archetype, t0$roi_id, t0$model, t0$rep)), ]
  rel <- c(median_area = max(abs(d$s_median_area / d0$s_median_area - 1)), cv = max(abs(d$s_cv / d0$s_cv - 1)),
           median_eq_r = max(abs(d$s_median_eq_r / d0$s_median_eq_r - 1)), frac_interior = max(abs(d$s_frac_interior / d0$s_frac_interior - 1)))
  vr("B", "B-2", v, "max variazione relativa: mediana area, CV, mediana eq_r, frazione interne", paste(fe(rel), collapse = ", "), "<= 1e-6", pf(all(rel <= 1e-6)), "<= 1e-6", all(rel <= 1e-6)) }
ref <- read.csv("results/S1.3/S1.3_verdicts.csv")
b2 <- lapply(c("V0", VS), function(v) read.csv(file.path(OUTB, "b2", sprintf("verdicts_%s.csv", v)))); names(b2) <- c("V0", VS)
same_all <- identical(b2$V0[, c("section", "id", "value", "outcome", "confirmed")], ref[, c("section", "id", "value", "outcome", "confirmed")])
vr("B", "B-2 controllo", "V0", "verdetti S1.3 riprodotti (tutte le righe)", sprintf("%s (%d righe)", same_all, nrow(ref)), "identici", pf(same_all), "identici", same_all)
kB <- ref$section == "B" & !is.na(ref$confirmed)
for (v in VS) { x <- b2[[v]]; kb <- ref$section == "B"
  s_out <- sum(x$outcome[kB] == ref$outcome[kB] & x$confirmed[kB] == ref$confirmed[kB]); s_val <- sum(x$value[kb] == ref$value[kb])
  vr("B", "B-2", v, "previsioni B invariate (esito); righe B con valore identico", sprintf("%d/%d; %d/%d", s_out, sum(kB), s_val, sum(kb)), sprintf("%d/%d", sum(kB), sum(kB)),
     pf(s_out == sum(kB)), sprintf("%d/%d", sum(kB), sum(kB)), s_out == sum(kB)) }

# ---- decisione ----------------------------------------------------------------------------------------------------
el <- sapply(VS, function(v) c(C1 = c1[[v]] == 0, C2 = c2[[v]] == 0, C3 = c3["far", v] == 0 && c3["anom", v] <= 50,
                                C4 = c4["inv", v] == 0 && c4["empty", v] == 0 && c4["n_bad", v] == 0 && c4["types", v] == 0,
                                C5 = c5[[v]] == 0, C6 = c6[[v]] == 0, Cunit = cu["pass", v] == cu["n", v], CP2 = cp2f[[v]] == 0))
tt <- sapply(VS, function(v) sum(tab$t_build[tab$variant == v] + tab$t_union[tab$variant == v]))
dec <- data.frame(variant = VS, t(el), eligible = apply(el, 2, all), deps = c("sf", "codice nostro", "geos + libgeos (GEOS 3.11.1)"), time_build_union_s = tt)
choice <- if (any(dec$eligible)) dec$variant[dec$eligible][1] else "P0"
dec$chosen <- dec$variant == choice
write.csv(dec, file.path(RES, "S1.3b_decision.csv"), row.names = FALSE)
vr("DEC", "regola", paste(dec$variant[dec$eligible], collapse = "+"), "variante scelta (ordine V1, V3, V2; nessuna → P0)", choice, "-", "INFO")
VV <- do.call(rbind, V); write.csv(VV, file.path(RES, "S1.3b_verdicts.csv"), row.names = FALSE)
print(dec, row.names = FALSE); print(VV[, c("section", "id", "variant", "value", "outcome")], row.names = FALSE)
