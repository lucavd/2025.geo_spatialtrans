# results/S1.3d/review/recalc/recalc_compare_indep.R — revisione «ricalcolo» S1.3d: confronto fra il ricalcolo indipendente (recalc_independent.R)
# e i confronti registrati nei res (letti da recalc_from_res.R -> review_recalc/recalc_tables.rds). Uso: Rscript --vanilla recalc_compare_indep.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; HEAVY <- file.path(D, "review_recalc"); REV <- "results/S1.3d/review/recalc"
TB <- readRDS(file.path(HEAVY, "recalc_tables.rds")); O1 <- TB$O1; O2 <- TB$O2
fs <- list.files(file.path(HEAVY, "indep"), full.names = TRUE)
I1 <- do.call(rbind, lapply(fs[grepl("_o1(F1)?\\.rds$", fs)], readRDS)); I2 <- do.call(rbind, lapply(fs[grepl("_o2(F1|V0)\\.rds$", fs)], readRDS))
KN <- readLines(file.path(D, "keys_N.txt")); KH <- readLines(file.path(D, "keys_H.txt")); KO2 <- readLines(file.path(D, "keys_O2.txt"))
st <- function(k) ifelse(k %in% KN, "N", ifelse(k %in% KH, "H", ifelse(k == "S_cperf", "S", "?")))
I1$set <- st(I1$key); I2$set <- st(I2$key)
cat("indipendente O1: righe", nrow(I1), "(", paste(names(table(paste(I1$set, I1$cand))), table(paste(I1$set, I1$cand)), collapse = "; "), ")\n")
cat("indipendente O2/O3: righe", nrow(I2), "(", paste(names(table(paste(I2$set, I2$cand, I2$op))), table(paste(I2$set, I2$cand, I2$op)), collapse = "; "), ")\n")
cat("coppie: elenco geom/_pairs vs st_intersects ricalcolato, chiavi con differenze:", sum(I1$pairs_mismatch > 0, na.rm = TRUE), "su", sum(!is.na(I1$pairs_mismatch)), "\n")
cat("verità celle: n_tc != n:", sum(!I1$tc_cells_ok), "\n")
m <- merge(I1, O1[, c("key", "cand", "union_err", "symd_err", "uv", "sv", "pv_ij", "pv_ji", "pv_any", "cell_viol", "exc", "silent", "fail", "fail_eff", "n_inv_out")], by = c("key", "cand"), suffixes = c(".i", ".r"))
cat("\nO1 abbinate", nrow(m), "\n")
cmpf <- function(nm, a, b) cat(sprintf("  %-28s diversi: %d\n", nm, sum(a != b, na.rm = TRUE) + sum(xor(is.na(a), is.na(b)))))
cmpf("union_viol", m$uv.i, m$uv.r); cmpf("symd_viol", m$sv.i, m$sv.r); cmpf("coppie ij", m$n_ij, m$pv_ij); cmpf("coppie ji", m$n_ji, m$pv_ji); cmpf("coppie any", m$n_any, m$pv_any)
cmpf("celle", m$n_cell, m$cell_viol); cmpf("eccezioni", m$n_exc > 0, m$exc > 0); cmpf("silent", m$silent.i, m$silent.r); cmpf("fail", m$fail.i, m$fail.r)
cat(sprintf("  union_err: max |Δ| %.3g ; identici %d/%d ; symd_err max |Δ| %.3g\n", max(abs(m$union_err.i - m$union_err.r), na.rm = TRUE), sum(m$union_err.i == m$union_err.r, na.rm = TRUE), nrow(m), max(abs(m$symd_err.i - m$symd_err.r), na.rm = TRUE)))
dd <- m[(m$silent.i != m$silent.r) | (m$n_any != m$pv_any) | (m$uv.i != m$uv.r), c("key", "cand", "set", "uv.i", "uv.r", "n_any", "pv_any", "n_cell", "cell_viol", "silent.i", "silent.r", "union_err.i", "union_err.r", "tauU")]
if (nrow(dd)) { cat("  righe discordi:\n"); print(head(dd, 30)) }
for (s in c("N", "H", "S")) for (cd in c("V0", "F1")) { d <- I1[I1$set == s & I1$cand == cd, ]; if (!nrow(d)) next
  cat(sprintf("  [indip] %s %s: casi %d, silenziosi %d, fail %d, unione %d, coppie %d, celle %d, eccezioni %d, max errore/tau coppie %.3g\n", s, cd, nrow(d), sum(d$silent), sum(d$fail), sum(d$uv), sum(d$n_any > 0), sum(d$n_cell > 0), sum(d$n_exc > 0), max(d$pmax, na.rm = TRUE))) }
f1 <- I1[I1$cand == "F1", ]; cat("  [indip] F1 celle invalide dopo geos_set_precision: N", sum(f1$n_invalid_cells[f1$set == "N"]), " H", sum(f1$n_invalid_cells[f1$set == "H"]), " S", sum(f1$n_invalid_cells[f1$set == "S"]), "\n")
# O2/O3
m2 <- merge(I2, O2[, c("key", "cand", "op", "n_truth", "n_cand", "n_viol", "n_missing", "n_spurious", "max_ratio", "fail")], by = c("key", "cand", "op"), suffixes = c(".i", ".r"))
cat("\nO2/O3 abbinate", nrow(m2), "su", nrow(I2), "\n")
cmpf("n_truth", m2$n_truth.i, m2$n_truth.r); cmpf("n_cand", m2$n_cand.i, m2$n_cand.r); cmpf("n_viol", m2$n_viol.i, m2$n_viol.r); cmpf("n_missing", m2$n_missing.i, m2$n_missing.r); cmpf("n_spurious", m2$n_spurious.i, m2$n_spurious.r)
cat(sprintf("  max_ratio: max |Δ| %.3g (rel %.3g)\n", max(abs(m2$max_ratio.i - m2$max_ratio.r)), max(abs(m2$max_ratio.i - m2$max_ratio.r) / m2$max_ratio.r)))
d2 <- m2[m2$n_viol.i != m2$n_viol.r | m2$n_cand.i != m2$n_cand.r, ]; if (nrow(d2)) print(head(d2, 20))
for (s in c("N", "H", "S")) for (cd in c("F1", "V0")) for (op in c("O2s2", "O2s8", "O3")) { d <- I2[I2$set == s & I2$cand == cd & I2$op == op, ]; if (!nrow(d)) next
  cat(sprintf("  [indip] %s %s %s: casi %d, con violazioni %d (pezzi violati %d; mancanti %d, spuri %d), max errore/tau %.3f, eccezioni %d; identità O3 (gruppi mancanti %d, extra %d, max rel %.2g)\n", s, cd, op, nrow(d), sum(d$n_viol > 0), sum(d$n_viol),
              sum(d$n_missing), sum(d$n_spurious), max(d$max_ratio), sum(d$n_exc > 0), sum(d$o3_ident_missing), sum(d$o3_ident_extra), max(d$o3_ident_maxrel))) }
h <- I1[I1$set == "H" & I1$cand == "F1", ]; h2 <- I2[I2$set == "H" & I2$cand == "F1", ]
cat(sprintf("\nCP-3 indipendente (F1 su H): O1 %d/%d ; O2 %d/%d ; O3 %d/%d\n", sum(h$fail), nrow(h), sum(h2$n_viol[h2$op != "O3"] > 0 | h2$n_exc[h2$op != "O3"] > 0), sum(h2$op != "O3"),
            sum(h2$n_viol[h2$op == "O3"] > 0 | h2$n_exc[h2$op == "O3"] > 0), sum(h2$op == "O3")))
write.csv(m, file.path(HEAVY, "recalc_indep_o1_vs_res.csv"), row.names = FALSE); write.csv(m2, file.path(HEAVY, "recalc_indep_o2_vs_res.csv"), row.names = FALSE)
cat("FINE\n")
