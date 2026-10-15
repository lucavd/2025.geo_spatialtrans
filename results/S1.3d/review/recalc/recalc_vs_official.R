# confronto della tabella ricalcolata con i CSV ufficiali (letti solo qui)
R <- read.csv("results/S1.3d/review/recalc/recalc_eligibility.csv"); O <- read.csv("results/S1.3d/S1.3d_eligibility.csv")
m <- merge(R, O, by = c("op", "cand"), suffixes = c(".rec", ".off")); cat("righe abbinate", nrow(m), "di", nrow(O), "\n")
R$n_inv_tot <- R$inv_cells_main + R$inv_na_main + R$inv_cells_cp1 + R$inv_na_cp1; m$n_inv_tot <- R$n_inv_tot[match(paste(m$op, m$cand), paste(R$op, R$cand))]
eqn <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
chk <- c(n_cases = sum(!eqn(m$n_rows, m$n_cases)), has_S = sum(!eqn(m$has_S.rec, m$has_S.off)), n_silent = sum(!eqn(m$n_silent.rec, m$n_silent.off)), n_exc = sum(!eqn(m$n_exc, m$n_exceptions)),
         n_cov = sum(!eqn(m$n_cov_fail.rec, m$n_cov_fail.off)), n_inv = sum(!eqn(m$n_inv_tot, m$n_still_invalid)), n_cp1 = sum(!eqn(m$n_cp1_rows, m$n_cp1)), n_cp1_fail = sum(!eqn(m$n_cp1_fail.rec, m$n_cp1_fail.off)),
         eligible = sum(!eqn(m$eligible.rec, m$eligible.off)), guaranteed = sum(!eqn(m$guaranteed.rec, m$guaranteed.off)), t_S = sum(!eqn(round(m$t_S.rec, 2), round(m$t_S.off, 2))), rss_S = sum(abs(m$rss_S.rec - m$rss_S.off) > 1e-9, na.rm = TRUE))
print(chk)
Dr <- read.csv("results/S1.3d/review/recalc/recalc_decision.csv"); Do <- read.csv("results/S1.3d/S1.3d_decision.csv")
cat("decisione ricalcolata:", paste(Dr$op, Dr$chosen, collapse = "; "), " | ufficiale:", paste(Do$op, Do$chosen_final, collapse = "; "), " | uguali:", identical(Dr$chosen, Do$chosen_final), "\n")
V <- read.csv("results/S1.3d/S1.3d_verdicts.csv"); print(V[V$id %in% c("C-0", "C-5", "P-V0a", "P-V0cov", "P-F4cov", "P-F2inv", "CP-1", "CP-3") , c("id", "cand", "value", "outcome")], row.names = FALSE)
