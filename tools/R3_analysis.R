# tools/R3_analysis.R — R3: check B, controprove e descrittivi dalle uscite di tools/R3_run.R e
# tools/R3_containment.py. Pre-registrazione: results/R3/R3_preregistration.md (commit 3fac4e7).
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/R3_analysis.R
# Uscite in results/R3/: R3_cells_summary.csv (archetipo x metodo), R3_roi_summary.csv, R3_null_summary.csv,
# R3_verdicts.csv (asserzione, valore, soglia, esito, previsione, previsione confermata), R3_cp2_roi.csv,
# R3_cp4.csv, R3_descriptive.csv, R3_bio_references_rows.csv; per cellula (pool) su micron: R3/cells_all.parquet.
# Definizioni operative aggiunte qui (piu' specifiche del testo, non divergenti):
#  - permutazione CP-2b: 2 000 permutazioni, seed = 20261008 + 1000*indice_ROI + 300;
#  - CV per ROI = sd/media delle aree interne del ROI; valore d'archetipo = mediana dei 5 ROI;
#  - CP-4: CV del pool grezzo delle aree interne delle 4 finestre (testo pre-registrato).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow) })
source("tools/R3_voronoi_metrics.R")
OUT <- "/mnt/micron/geo_spatialtrans/R3"; RES <- "results/R3"; BASE_SEED <- 20261008
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   design_ratio = c(0.45, 0.45, 0.45, 0.55, NA, NA), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]; rois$roi_idx <- seq_len(nrow(rois))

# ---- tabella per cellula -------------------------------------------------------------------------
rr <- do.call(rbind, lapply(list.files(file.path(OUT, "real"), "_roi\\.rds$", full.names = TRUE), readRDS))
cells <- do.call(rbind, lapply(seq_len(nrow(rr)), function(k) {
  A <- rr$archetype[k]; roi <- rr$roi_id[k]; me <- rr$method[k]
  tv <- as.data.frame(read_parquet(file.path(OUT, "real", sprintf("%s_%s_%s_cells.parquet", A, roi, me))))
  nu <- as.data.frame(read_parquet(file.path(OUT, "real", sprintf("%s_%s_%s_nuc.parquet", A, roi, me))))
  nu <- nu[match(tv$idx, nu$idx), c("ecc_N", "theta_N", "frac_out", "hed_label", "hed_area_um2")]
  data.frame(archetype = A, roi_id = roi, method = me, tv[, c("idx", "x", "y", "area", "area_tile", "interior", "nsides",
             "ecc_T", "theta_T", "lambda_loc", "area_nuc")], nu)
}))
cells$eq_r <- sqrt(cells$area / pi); cells$nc <- cells$area_nuc / cells$area; cells$ratio <- sqrt(cells$nc)
write_parquet(cells, file.path(OUT, "cells_all.parquet"))
ci <- cells[cells$interior, ]
is_prim <- function(d) d$method == ARCH$primary[match(d$archetype, ARCH$archetype)]
cp <- ci[is_prim(ci), ]

q <- function(v, p) as.numeric(quantile(v, p, na.rm = TRUE))
arch_sum <- do.call(rbind, lapply(split(ci, list(ci$archetype, ci$method), drop = TRUE), function(d) {
  rm <- tapply(d$area, d$roi_id, median); re <- tapply(d$eq_r, d$roi_id, median)
  data.frame(archetype = d$archetype[1], method = d$method[1], n_interior = nrow(d),
             area_median = median(d$area), area_q25 = q(d$area, .25), area_q75 = q(d$area, .75),
             area_q05 = q(d$area, .05), area_q95 = q(d$area, .95), area_roi_med_min = min(rm), area_roi_med_max = max(rm),
             eq_r_median = median(d$eq_r), eq_r_q25 = q(d$eq_r, .25), eq_r_q75 = q(d$eq_r, .75),
             eq_r_roi_med_min = min(re), eq_r_roi_med_max = max(re),
             nc_median = median(d$nc), nc_q25 = q(d$nc, .25), nc_q75 = q(d$nc, .75),
             ratio_median = median(d$ratio), ratio_q25 = q(d$ratio, .25), ratio_q75 = q(d$ratio, .75),
             ecc_T_median = median(d$ecc_T), ecc_T_q25 = q(d$ecc_T, .25), ecc_T_q75 = q(d$ecc_T, .75),
             nsides_mean = mean(d$nsides), nsides_var = var(d$nsides),
             spearman_nuc_terr = suppressWarnings(cor(d$area_nuc, d$area, method = "spearman")),
             frac_cut = mean(d$frac_out > 0.05, na.rm = TRUE), frac_out_median = median(d$frac_out, na.rm = TRUE),
             hed_cover = mean(d$hed_label > 0), hed_ratio_median = median(d$hed_area_um2 / d$area, na.rm = TRUE))
}))
arch_sum <- arch_sum[order(arch_sum$archetype, arch_sum$method), ]
write.csv(arch_sum, file.path(RES, "R3_cells_summary.csv"), row.names = FALSE)
write.csv(rr, file.path(RES, "R3_roi_summary.csv"), row.names = FALSE)

# ---- nulli ----------------------------------------------------------------------------------------
nf <- list.files(file.path(OUT, "null"), "\\.rds$", full.names = TRUE)
nl <- lapply(nf, readRDS)
ns <- do.call(rbind, lapply(nl, function(o) data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep,
                                                        n_failed = o$n_failed, conserv = abs(o$sum_area / o$area_poly - 1), t(o$summary))))
write.csv(ns, file.path(RES, "R3_null_summary.csv"), row.names = FALSE)
nsides_null <- do.call(rbind, lapply(nl, function(o) data.frame(archetype = o$archetype, model = o$model, k = 1:15, count = o$nsides)))

V <- list()
verdict <- function(id, A, stat, value, thr, pass, pred, note = "") {
  V[[length(V) + 1]] <<- data.frame(id = id, archetype = A, statistic = stat, value = signif(value, 5), threshold = thr,
                                     assertion = if (is.na(pass)) "registrato" else if (pass) "PASS" else "FAIL",
                                     predicted = pred,
                                     prediction_confirmed = if (is.na(pass) || pred == "registrare") NA else
                                       (pred == "PASS" && pass) || (pred == "FAIL" && !pass), note = note)
}
P <- function(A) arch_sum[arch_sum$archetype == A & arch_sum$method == ARCH$primary[ARCH$archetype == A], ]

# ---- check B ----------------------------------------------------------------------------------------
for (A in ARCH$archetype) verdict("B-R3.1", A, "eq_r mediano (µm)", P(A)$eq_r_median, "[5, 25]",
                                  P(A)$eq_r_median >= 5 && P(A)$eq_r_median <= 25, if (A == "A4") "FAIL" else "PASS")
for (A in ARCH$archetype) {
  dr <- ARCH$design_ratio[ARCH$archetype == A]; v <- P(A)$ratio_median
  if (is.na(dr)) verdict("B-R3.2", A, "rapporto dei raggi mediano", v, "nessun valore di catalogo", NA, "registrare")
  else verdict("B-R3.2", A, "rapporto dei raggi mediano", v, sprintf("%.2f ± 20%%", dr), abs(v / dr - 1) <= 0.20,
               c(A1 = "FAIL", A2 = "PASS", A3 = "FAIL", A4 = "PASS")[[A]])
}
verdict("B-R3.3a", "A4", "N/C mediano", P("A4")$nc_median, "[0.8, 0.9]", P("A4")$nc_median >= 0.8 && P("A4")$nc_median <= 0.9, "FAIL")
verdict("B-R3.3b", "A3", "N/C mediano A3 (vs A1, A2)", P("A3")$nc_median,
        sprintf("< %.3f e < %.3f", P("A1")$nc_median, P("A2")$nc_median),
        P("A3")$nc_median < P("A1")$nc_median && P("A3")$nc_median < P("A2")$nc_median, "PASS")
d5 <- rr[rr$archetype == "A5" & rr$method == "cellpose_rgb", ]
verdict("B-R3.3c", "A5", "ROI con dip p < 0.05 (aree territori)", sum(d5$dip_p_area < 0.05), ">= 4/5", sum(d5$dip_p_area < 0.05) >= 4, "FAIL",
        paste(sprintf("%s p=%.3g", d5$roi_id, d5$dip_p_area), collapse = "; "))
nsup <- 0
for (A in ARCH$archetype) { s <- P(A)$spearman_nuc_terr; nsup <- nsup + (s >= 0.3)
  verdict("B-R3.4", A, "Spearman(area nucleo, area territorio)", s, ">= 0.3 (regola sostenuta)", s >= 0.3, "registrare") }
verdict("B-R3.4", "tutti", "archetipi con regola sostenuta", nsup, "<= 2/6 (previsione)", nsup <= 2, "PASS")

# ---- CP-1 -------------------------------------------------------------------------------------------
cp1 <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; me <- ARCH$primary[ARCH$archetype == A]
  real <- rr[rr$archetype == A & rr$roi_id == roi & rr$method == me, ]
  cs <- ns[ns$archetype == A & ns$roi_id == roi & ns$model == "CSR", ]; rs <- ns[ns$archetype == A & ns$roi_id == roi & ns$model == "RSA", ]
  pos <- function(v, ref) if (v < min(ref)) "sotto" else if (v > max(ref)) "sopra" else "dentro"
  data.frame(archetype = A, roi_id = roi, cv_real = real$cv, cv_loc_real = real$cv_loc,
             cv_csr_min = min(cs$cv), cv_csr_max = max(cs$cv), cv_loc_csr_min = min(cs$cv_loc), cv_loc_csr_max = max(cs$cv_loc),
             cv_loc_rsa_median = median(rs$cv_loc), cv_rsa_median = median(rs$cv),
             pos_cv_csr = pos(real$cv, cs$cv), pos_cvloc_csr = pos(real$cv_loc, cs$cv_loc),
             rel_err_rsa_cvloc = abs(median(rs$cv_loc) - real$cv_loc) / real$cv_loc,
             ecc_T_real = real$median_ecc_T, ecc_T_rsa_median = median(rs$median_ecc_T), ecc_T_csr_median = median(cs$median_ecc_T),
             n_rep_csr = nrow(cs), n_rep_rsa = nrow(rs), rsa_n_failed = sum(rs$n_failed))
}))
write.csv(cp1, file.path(RES, "R3_cp1_roi.csv"), row.names = FALSE)
for (A in ARCH$archetype) {
  d <- cp1[cp1$archetype == A, ]
  if (A %in% c("A4", "A5")) verdict("CP-1a", A, "ROI con CV_loc reale sotto CSR", sum(d$pos_cvloc_csr == "sotto"), ">= 4/5", sum(d$pos_cvloc_csr == "sotto") >= 4, "PASS")
  else if (A == "A6") verdict("CP-1a", A, "ROI con CV_loc reale sopra CSR", sum(d$pos_cvloc_csr == "sopra"), ">= 4/5", sum(d$pos_cvloc_csr == "sopra") >= 4, "PASS")
  else verdict("CP-1a", A, "posizione CV_loc reale vs CSR (sotto/dentro/sopra)", sum(d$pos_cvloc_csr == "sopra"), "registrare", NA, "registrare",
               paste(d$pos_cvloc_csr, collapse = ","))
  nrep <- sum(d$rel_err_rsa_cvloc <= 0.10)
  if (A %in% c("A4", "A5", "A6")) verdict("CP-1b", A, "ROI con RSA(d*) entro ±10% del CV_loc reale", nrep, ">= 4/5", nrep >= 4, "PASS")
  else verdict("CP-1b", A, "ROI con RSA(d*) fuori da ±10% del CV_loc reale", 5 - nrep, ">= 3/5", (5 - nrep) >= 3, "PASS")
}

# ---- CP-2 -------------------------------------------------------------------------------------------
cp2 <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; me <- ARCH$primary[ARCH$archetype == A]
  d <- cp[cp$archetype == A & cp$roi_id == roi & cp$ecc_N >= 0.8 & cp$ecc_T >= 0.5 & is.finite(cp$theta_N), ]
  dth <- r3_axis_diff(d$theta_T, d$theta_N) * 180 / pi; obs <- median(dth)
  set.seed(BASE_SEED + 1000 * rois$roi_idx[i] + 300)
  perm <- replicate(2000, median(r3_axis_diff(d$theta_T, sample(d$theta_N)) * 180 / pi))
  data.frame(archetype = A, roi_id = roi, n_eligible = nrow(d), median_dtheta = obs, perm_p = mean(perm <= obs),
             perm_median = median(perm), frac_dtheta_lt30 = mean(dth < 30))
}))
cp2 <- merge(cp2, cp1[, c("archetype", "roi_id", "ecc_T_real", "ecc_T_rsa_median", "ecc_T_csr_median")])
cp2$d_ecc_rsa <- cp2$ecc_T_real - cp2$ecc_T_rsa_median
write.csv(cp2, file.path(RES, "R3_cp2_roi.csv"), row.names = FALSE)
d6 <- cp2[cp2$archetype == "A6", ]
verdict("CP-2a", "A6", "ROI con mediana e_T reale − RSA >= 0.05", sum(d6$d_ecc_rsa >= 0.05), ">= 4/5", sum(d6$d_ecc_rsa >= 0.05) >= 4, "PASS")
for (A in c("A6", "A1")) { d <- cp2[cp2$archetype == A, ]; k <- sum(d$median_dtheta <= 35 & d$perm_p < 0.01)
  verdict("CP-2b", A, "ROI con mediana Δθ <= 35° e p perm < 0.01", k, ">= 4/5", k >= 4, "PASS",
          if (A == "A1") "bassa confidenza (pre-registrato)" else "") }
d4 <- cp2[cp2$archetype == "A4", ]
verdict("CP-2c", "A4", "ROI con |mediana e_T reale − RSA| < 0.05", sum(abs(d4$d_ecc_rsa) < 0.05), ">= 4/5", sum(abs(d4$d_ecc_rsa) < 0.05) >= 4, "PASS")
for (A in setdiff(ARCH$archetype, c("A1", "A6"))) { d <- cp2[cp2$archetype == A, ]
  verdict("CP-2b", A, "ROI con mediana Δθ <= 35° e p perm < 0.01", sum(d$median_dtheta <= 35 & d$perm_p < 0.01), "registrare", NA, "registrare") }

# ---- CP-3 -------------------------------------------------------------------------------------------
for (A in ARCH$archetype) { v <- P(A)$frac_cut
  if (A == "A1") verdict("CP-3", A, "frazione nuclei tagliati (>5% fuori)", v, ">= 0.10", v >= 0.10, "PASS")
  else if (A == "A5") verdict("CP-3", A, "frazione nuclei tagliati (>5% fuori)", v, "<= 0.05", v <= 0.05, "PASS")
  else verdict("CP-3", A, "frazione nuclei tagliati (>5% fuori)", v, "registrare", NA, "registrare") }

# ---- CP-4 -------------------------------------------------------------------------------------------
cal <- read.csv(file.path(RES, "R3_calibration_windows.csv"))
cc <- do.call(rbind, lapply(list.files(file.path(OUT, "calib"), "\\.parquet$", full.names = TRUE), function(f) as.data.frame(read_parquet(f))))
cc <- cc[cc$interior, ]; cc$archetype <- substr(cc$win_id, 1, 2); cc$kind <- ifelse(cc$source == "manual_Luca", "manual", "seg")
cp4 <- do.call(rbind, lapply(ARCH$archetype, function(A) {
  m <- cc$area[cc$archetype == A & cc$kind == "manual"]; s <- cc$area[cc$archetype == A & cc$kind == "seg"]
  data.frame(archetype = A, n_int_manual = length(m), n_int_seg = length(s), cv_manual = sd(m) / mean(m), cv_seg = sd(s) / mean(s),
             eq_r_med_manual = median(sqrt(m / pi)), eq_r_med_seg = median(sqrt(s / pi)),
             ratio_eqr_seg_over_manual = median(sqrt(s / pi)) / median(sqrt(m / pi)),
             n_manual = sum(cal$n[cal$archetype == A & cal$source == "manual_Luca"]),
             n_seg = sum(cal$n[cal$archetype == A & cal$source != "manual_Luca"]),
             conserv_max = max(abs(cal$sum_area[cal$archetype == A] / cal$area_region[cal$archetype == A] - 1)))
}))
cp4$sqrt_n_ratio <- sqrt(cp4$n_manual / cp4$n_seg)
write.csv(cp4, file.path(RES, "R3_cp4.csv"), row.names = FALSE)
ncv <- sum(cp4$cv_seg > cp4$cv_manual)
for (A in ARCH$archetype) { d <- cp4[cp4$archetype == A, ]
  verdict("CP-4", A, "CV_seg − CV_manuale", d$cv_seg - d$cv_manual, "> 0", d$cv_seg > d$cv_manual, "registrare") }
verdict("CP-4", "tutti", "archetipi con CV_seg > CV_manuale", ncv, ">= 5/6", ncv >= 5, "PASS")

Vd <- do.call(rbind, V); write.csv(Vd, file.path(RES, "R3_verdicts.csv"), row.names = FALSE)

# ---- descrittivi ------------------------------------------------------------------------------------
desc <- do.call(rbind, lapply(ARCH$archetype, function(A) {
  d <- rr[rr$archetype == A & rr$method == ARCH$primary[ARCH$archetype == A], ]
  nsr <- nsides_null[nsides_null$archetype == A, ]
  mv <- function(mod) { z <- nsr[nsr$model == mod, ]; k <- rep(z$k, z$count); c(mean(k), var(k)) }
  data.frame(archetype = A, nn_median_minus_roi_median = median(d$nn_median_minus), clark_evans_cdf_roi_median = median(d$clark_evans_cdf),
             clark_evans_min = min(d$clark_evans_cdf), clark_evans_max = max(d$clark_evans_cdf),
             frac_interior_median = median(d$frac_interior), mean_area_x_lambda = median(d$mean_area * d$lambda_mm2 / 1e6),
             nsides_mean_real = P(A)$nsides_mean, nsides_var_real = P(A)$nsides_var,
             nsides_mean_csr = mv("CSR")[1], nsides_var_csr = mv("CSR")[2], nsides_mean_rsa = mv("RSA")[1], nsides_var_rsa = mv("RSA")[2],
             hed_cover = P(A)$hed_cover, hed_ratio_median = P(A)$hed_ratio_median)
}))
write.csv(desc, file.path(RES, "R3_descriptive.csv"), row.names = FALSE)
sens <- rbind(
  data.frame(case = "A3 senza r3", archetype = "A3", { d <- cp[cp$archetype == "A3" & cp$roi_id != "r3", ]
    data.frame(n_interior = nrow(d), area_median = median(d$area), eq_r_median = median(d$eq_r), nc_median = median(d$nc), ratio_median = median(d$ratio), ecc_T_median = median(d$ecc_T)) }),
  data.frame(case = "A4 follicolo (f1, f2)", archetype = "A4", { d <- cp[cp$archetype == "A4" & cp$roi_id %in% c("f1", "f2"), ]
    data.frame(n_interior = nrow(d), area_median = median(d$area), eq_r_median = median(d$eq_r), nc_median = median(d$nc), ratio_median = median(d$ratio), ecc_T_median = median(d$ecc_T)) }),
  data.frame(case = "A4 paracorticale (p1–p3)", archetype = "A4", { d <- cp[cp$archetype == "A4" & cp$roi_id %in% c("p1", "p2", "p3"), ]
    data.frame(n_interior = nrow(d), area_median = median(d$area), eq_r_median = median(d$eq_r), nc_median = median(d$nc), ratio_median = median(d$ratio), ecc_T_median = median(d$ecc_T)) }))
write.csv(sens, file.path(RES, "R3_sensitivity.csv"), row.names = FALSE)

# ---- righe candidate per BIO_REFERENCES (intervallo [corretto, segmentatore]) -----------------------
bio <- do.call(rbind, lapply(ARCH$archetype, function(A) {
  p <- P(A); k <- cp4$ratio_eqr_seg_over_manual[cp4$archetype == A]
  data.frame(archetype = A, method = p$method, n_interior = p$n_interior,
             area_median_seg = p$area_median, area_median_corr = p$area_median / k^2, area_q25 = p$area_q25, area_q75 = p$area_q75,
             area_roi_med_min = p$area_roi_med_min, area_roi_med_max = p$area_roi_med_max,
             eq_r_median_seg = p$eq_r_median, eq_r_median_corr = p$eq_r_median / k, eq_r_q25 = p$eq_r_q25, eq_r_q75 = p$eq_r_q75,
             cv_roi_median = median(rr$cv[rr$archetype == A & rr$method == p$method]),
             cv_loc_roi_median = median(rr$cv_loc[rr$archetype == A & rr$method == p$method]),
             ecc_T_median = p$ecc_T_median, nc_median = p$nc_median, ratio_median = p$ratio_median,
             frac_cut = p$frac_cut, k_eqr = k)
}))
write.csv(bio, file.path(RES, "R3_bio_references_rows.csv"), row.names = FALSE)
cat("verdetti:", nrow(Vd), " PASS", sum(Vd$assertion == "PASS"), " FAIL", sum(Vd$assertion == "FAIL"),
    " previsioni confermate", sum(Vd$prediction_confirmed, na.rm = TRUE), "/", sum(!is.na(Vd$prediction_confirmed)), "\n")
