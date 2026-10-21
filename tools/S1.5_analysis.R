# tools/S1.5_analysis.R - S1.5: verdetti di P-1...P-8 e CP-1, regola sul modello di taglia, descrittive.
# Regole: docs/steps/S1.5.md (protocollo v2: stima +- 2 errori standard contro l'intervallo di rilevanza; frazioni tagliate a [0, 1];
# errore standard = deviazione standard fra unita' / sqrt(unita')). Uso: Rscript --vanilla tools/S1.5_analysis.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
RES <- "results/S1.5"; OUT <- "/mnt/micron/geo_spatialtrans/S1.5"
TIS <- c(A1 = "intestino", A2 = "tumore del colon", A3 = "stroma", A4 = "linfonodo", A5 = "cervello", A6 = "cuore")
u <- read.csv(file.path(RES, "S1.5_panel_units.csv")); rr <- read.csv(file.path(RES, "S1.5_real_refs_roi.csv"))
verdict <- function(est, se, lo, hi, frac = FALSE) {
  a <- est - 2 * se; b <- est + 2 * se; if (frac) { a <- max(0, a); b <- min(1, b) }
  if (!is.finite(est)) return(c(lo2 = NA, hi2 = NA, outcome = "non valutabile"))
  out <- if (a >= lo && b <= hi) "confermata" else if (b < lo || a > hi) "smentita" else "indeterminata"
  c(lo2 = a, hi2 = b, outcome = out)
}
se_of <- function(x) sd(x) / sqrt(length(x))
roi_means <- function(arm) aggregate(. ~ archetype + roi_id, data = u[u$arm == arm, c("archetype", "roi_id", "median_area_nuc", "median_nc", "frac_nc095",
  "frac_trunc", "cv_terr", "median_ecc_T", "median_eq_r", "frac_pale", "desc_border_ratio", "n_failed", "n_warn", "n_error")], FUN = mean)
rows <- list()
add <- function(id, A, arm, qty, est, se, lo, hi, src, frac = FALSE, n_units = NA) {
  v <- verdict(est, se, lo, hi, frac)
  rows[[length(rows) + 1]] <<- data.frame(id = id, archetype = A, tissue = TIS[[A]], arm = arm, quantity = qty, estimate = est, se = se,
                                          lo2 = as.numeric(v[["lo2"]]), hi2 = as.numeric(v[["hi2"]]), int_lo = lo, int_hi = hi, interval_source = src,
                                          outcome = v[["outcome"]], n_units = n_units)
}
for (arm in c("preserve_area", "none")) {
  m <- roi_means(arm); m <- merge(m, rr, by = c("archetype", "roi_id"), suffixes = c("", "_real"))
  for (A in names(TIS)) {
    s <- m[m$archetype == A, ]; r <- rr[rr$archetype == A, ]
    add("P-1", A, arm, "area mediana del nucleo (um2), celle interne", mean(s$median_area_nuc), se_of(s$median_area_nuc),
        min(r$median_area_nuc), max(r$median_area_nuc), "variabilita' fra ROI reali", n_units = nrow(s))
    add("P-2", A, arm, "rapporto mediano nucleo/territorio, celle interne", mean(s$median_nc), se_of(s$median_nc),
        min(r$median_nc), max(r$median_nc), "variabilita' fra ROI reali", frac = TRUE, n_units = nrow(s))
    add("P-3", A, arm, "frazione di nuclei troncati dal tetto", mean(s$frac_trunc), se_of(s$frac_trunc), 0, 0.05, "criterio di Luca (ROADMAP §4, S1.4)", frac = TRUE, n_units = nrow(s))
    add("P-4", A, arm, "frazione di cellule quasi senza citoplasma, celle interne", mean(s$frac_nc095), se_of(s$frac_nc095), 0, max(r$frac_nc095),
        "variabilita' fra ROI reali (massimo)", frac = TRUE, n_units = nrow(s))
    rel <- s$cv_terr / s$cv_terr_real - 1
    add("P-6", A, arm, "CV dell'area dei territori, simulato / reale - 1", mean(rel), se_of(rel), -0.10, 0.10, "soglia S1.3 B-2 (Luca, 2026-10-08)", n_units = nrow(s))
    de <- s$median_ecc_T - s$median_ecc_T_real
    add("P-7", A, arm, "eccentricita' mediana del territorio, simulato - reale", mean(de), se_of(de), -0.05, 0.05, "soglia S1.3 B-4 (Luca, 2026-10-08)", n_units = nrow(s))
  }
}
# P-5: finestre di R2c, media a peso uguale delle 4 finestre per seme, poi semi come unita'
B127 <- data.frame(archetype = names(TIS), lo = c(0.366, 0.217, 0.083, 0.350, 0.000, 0.000), hi = c(0.650, 0.567, 0.283, 0.750, 0.217, 0.070),
                   src = paste0("B-", 127:132, " (manuale, due stadi)"))
wf <- file.path(RES, "S1.5_windows.csv")
if (file.exists(wf)) {
  w <- read.csv(wf); ps <- aggregate(s1 ~ arm + archetype + seed, data = w, FUN = mean)
  for (arm in c("preserve_area", "none")) for (A in names(TIS)) {
    x <- ps$s1[ps$arm == arm & ps$archetype == A]; b <- B127[B127$archetype == A, ]
    add("P-5", A, arm, "frazione di nuclei a contatto, finestre di R2c a peso uguale", mean(x), se_of(x), b$lo, b$hi, b$src, frac = TRUE, n_units = length(x))
  }
}
V <- do.call(rbind, rows)
# regola della scelta 2: piu' tessuti con P-1 confermata; poi P-2; poi preserve_area
cnt <- function(arm, id) sum(V$arm == arm & V$id == id & V$outcome == "confermata")
c1 <- c(preserve_area = cnt("preserve_area", "P-1"), none = cnt("none", "P-1")); c2 <- c(preserve_area = cnt("preserve_area", "P-2"), none = cnt("none", "P-2"))
win <- if (c1[1] != c1[2]) names(which.max(c1)) else if (c2[1] != c2[2]) names(which.max(c2)) else "preserve_area"
rule <- data.frame(arm = names(c1), P1_confirmed = as.integer(c1), P2_confirmed = as.integer(c2), winner = names(c1) == win)
write.csv(rule, file.path(RES, "S1.5_size_rule.csv"), row.names = FALSE); print(rule)
V$winner_arm <- V$arm == win
# P-8: tempi delle 3 ripetizioni (riscaldamento escluso) e memoria di picco del processo
pf <- list.files(file.path(OUT, "perf"), pattern = "\\.csv$", full.names = TRUE)
if (length(pf)) {
  P <- do.call(rbind, lapply(pf, read.csv)); P <- P[!P$warmup, ]
  arch <- c(epithelial = "A1", tumor_microenv = "A2", stromal_rich = "A3", lymphoid = "A4", brain_cortex = "A5", columnar_muscle = "A6")
  for (pr in unique(P$preset)) { x <- P$elapsed_s[P$preset == pr]; A <- arch[[pr]]
    rows <- list(); add("P-8", A, "preserve_area", "tempo di place_cells (s), 4000 x 4000 px", mean(x), se_of(x), 0, 300, "design v1.1, casi limite", n_units = length(x)); V <- rbind(V, cbind(do.call(rbind, rows), winner_arm = TRUE))
    tf <- file.path(OUT, "perf", paste0(pr, ".time")); if (file.exists(tf)) { kb <- as.numeric(sub(".*: ", "", grep("Maximum resident", readLines(tf), value = TRUE)))
      rows <- list(); add("P-8", A, "preserve_area", "memoria di picco del processo (GB)", kb / 1024^2, 0, 0, 8, "design v1.1, casi limite", n_units = 1); V <- rbind(V, cbind(do.call(rbind, rows), winner_arm = TRUE)) } }
}
write.csv(V, file.path(RES, "S1.5_verdicts.csv"), row.names = FALSE)
print(V[V$winner_arm, c("id", "archetype", "estimate", "lo2", "hi2", "int_lo", "int_hi", "outcome")], digits = 4)
# CP-1
cf <- file.path(RES, "S1.5_CP1_scale.csv")
if (file.exists(cf)) {
  K <- read.csv(cf); out <- list()
  for (v in unique(K$variant)) for (s in sort(unique(K$scale_true))) for (A in names(TIS)) {
    k <- K[K$variant == v & K$scale_true == s & K$archetype == A, ]; ok <- k[!k$refused %in% TRUE & is.finite(k$ratio), ]
    est <- if (nrow(ok) >= 3) mean(ok$ratio) else NA; se <- if (nrow(ok) >= 3) se_of(ok$ratio) else NA
    vv <- if (nrow(ok) >= 3) verdict(est, se, 0.8, 1.2) else c(lo2 = NA, hi2 = NA, outcome = "rifiuto (>= 3 ROI su 5)")
    out[[length(out) + 1]] <- data.frame(variant = v, scale_true = s, archetype = A, tissue = TIS[[A]], n_roi = nrow(k), n_refused = sum(k$refused %in% TRUE),
      n_number_outside = sum(!k$refused %in% TRUE & is.finite(k$ratio) & (k$ratio < 0.8 | k$ratio > 1.2)), ratio_mean = est, se = se,
      lo2 = as.numeric(vv[["lo2"]]), hi2 = as.numeric(vv[["hi2"]]), outcome = vv[["outcome"]])
  }
  CP <- do.call(rbind, out); write.csv(CP, file.path(RES, "S1.5_CP1_verdicts.csv"), row.names = FALSE); print(CP, digits = 3)
}
# descrittive
D <- aggregate(cbind(desc_border_ratio, n_failed, n_warn, frac_pale) ~ archetype, data = u[u$arm == win, ], FUN = mean)
write.csv(D, file.path(RES, "S1.5_descriptive.csv"), row.names = FALSE); print(D, digits = 3)
