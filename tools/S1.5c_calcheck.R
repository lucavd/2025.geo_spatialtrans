# tools/S1.5c_calcheck.R - S1.5c, cancello G-4: (a) le etichette del pacchetto (stesse formule di R/04b2_seed_centroids.R: campo a caratteristiche di Fourier,
# softmax, spostamenti a_t; generatore di R) estratte sulle posizioni dei nuclei veri con i parametri di tutte le ROI riproducono le curve di calibrazione di
# tools/S1.5c_calibration.py (media di 10 repliche; |differenza dei logaritmi| <= 0.10 nei bin con almeno 200 coppie reali dello stesso tipo);
# (b) le unita' del banco (seme 1, bracci a compartimenti) hanno nei metadati i parametri della piega lascia-fuori del proprio ROI (30/30 per braccio).
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5c_calcheck.R [a|b|ab]   Uscita: results/S1.5c/S1.5c_G4_calcheck.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow); library(spatstat.geom) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/S1.5c_bench_lib.R")
options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")   # S1.6: predefiniti di S1.5c-S1.5e, per riprodurre i numeri di questo step
what <- commandArgs(TRUE); what <- if (length(what)) what[1] else "ab"
RES <- "results/S1.5c"; M <- "/mnt/micron/geo_spatialtrans"
UNTYPED <- c("non risolto", "non classificabile"); EPI <- c("epiteliale normale" = "epiteliale", "epiteliale tumorale" = "epiteliale")
BINS <- c(0, 10, 20, 40, 80, 160); NFOC <- 2000L; NREP <- 10L
rows <- list()
if (grepl("a", what)) for (A in paste0("A", 1:6)) {
  lab <- as.data.frame(read_parquet(sprintf("%s/R6/labels/%s_exp_umi20_seed20261022_nuclei_ph1.parquet", M, A)))
  lab$type <- ifelse(lab$cell_type %in% names(EPI), EPI[lab$cell_type], lab$cell_type); lab <- lab[!lab$type %in% UNTYPED, ]
  types <- sort(unique(lab$type)); f <- as.numeric(table(factor(lab$type, levels = types))) / nrow(lab)
  cp <- compart_for(A, NA, "pooled"); s <- cp$s[match(types, cp$cell_type)]; s[is.na(s)] <- 0; ell <- cp$ell_um[match(types, cp$cell_type)]
  main <- which(types %in% cp$cell_type); K <- length(types)
  num_r <- den_r <- num_s <- den_s <- matrix(0, K, length(BINS) - 1)
  set.seed(20261023)
  for (r in sort(unique(lab$roi_id))) {
    e <- lab[lab$roi_id == r, ]; L <- match(e$type, types); n <- nrow(e)
    foc <- if (n <= NFOC) seq_len(n) else sort(sample.int(n, NFOC))
    W <- owin(range(e$x_um) + c(-1, 1), range(e$y_um) + c(-1, 1))
    X <- ppp(e$x_um[foc], e$y_um[foc], window = W); Y <- ppp(e$x_um, e$y_um, window = W)
    cpp <- crosspairs(X, Y, rmax = BINS[length(BINS)], what = "ijd"); fi <- foc[cpp$i]; j <- cpp$j; k <- fi != j; fi <- fi[k]; j <- j[k]
    b <- findInterval(cpp$d[k], BINS, rightmost.closed = TRUE); b[b >= length(BINS)] <- length(BINS) - 1
    cnt <- function(Lv) { idx <- (b - 1) * K + Lv[fi]; same <- Lv[fi] == Lv[j]
      list(num = matrix(tabulate(idx[same], K * (length(BINS) - 1)), K), den = matrix(tabulate(idx, K * (length(BINS) - 1)), K)) }
    cr <- cnt(L); num_r <- num_r + cr$num; den_r <- den_r + cr$den
    for (rep in seq_len(NREP)) {
      fields <- lapply(seq_len(K), function(t) if (s[t] > 0) .sc_rff_draw(ell[t], 256L) else NULL)
      eta <- .sc_eta(fields, s, log(f), e$x_um, e$y_um); a <- .sc_solve_a(eta, f, rep(1, K))     # densita' uguali: media semplice sui punti (come la calibrazione)
      P <- .sc_prob(eta, a, rep(1, K))$p; Ls <- .sc_draw_type(P, runif(n))
      cs <- cnt(Ls); num_s <- num_s + cs$num; den_s <- den_s + cs$den
    }
  }
  fs <- f; Rr <- num_r / pmax(den_r, 1) / fs; Rs <- num_s / pmax(den_s, 1) / fs
  py <- read.csv(sprintf("%s/S1.5c_compartments_curves_%s.csv", RES, A)); py <- py[py$fold == "pooled", ]
  for (t in main) for (kb in seq_len(length(BINS) - 1)) {
    bn <- sprintf("%d-%d", BINS[kb], BINS[kb + 1]); p <- py[py$cell_type == types[t] & py$bin == bn, ]
    rows[[length(rows) + 1]] <- data.frame(part = "a", archetype = A, cell_type = types[t], bin = bn, real_pairs = num_r[t, kb],
      R_real_R = Rr[t, kb], R_real_py = p$real, R_fit_R = Rs[t, kb], R_fit_py = p$fitted, dlog = abs(log(Rs[t, kb]) - log(p$fitted)),
      checked = num_r[t, kb] >= 200, pass = num_r[t, kb] < 200 || abs(log(Rs[t, kb]) - log(p$fitted)) <= 0.10)
  }
  cat(A, "fatto\n")
}
if (grepl("b", what)) {
  rois <- read.csv("results/R2/R2_rois_checked.csv")
  for (arm in c("C-a", "C-b", "N-a", "N-b", "G")) for (k in seq_len(nrow(rois))) {
    A <- rois$archetype[k]; r <- rois$roi_id[k]; L <- readRDS(sprintf("%s/S1.5c/layers/%s_%s_%s.rds", M, A, r, arm))
    got <- L$metadata$compartments; exp <- compart_for(A, r); m <- match(exp$cell_type, got$cell_type)
    ok <- all(!is.na(m)) && all(got$s[m] == exp$s) && all(got$ell_um[m] == exp$ell_um)
    rows[[length(rows) + 1]] <- data.frame(part = "b", archetype = A, cell_type = arm, bin = r, real_pairs = NA, R_real_R = NA, R_real_py = NA,
                                           R_fit_R = NA, R_fit_py = NA, dlog = NA, checked = TRUE, pass = ok)
  }
}
out <- do.call(rbind, rows); f <- file.path(RES, sprintf("S1.5c_G4_calcheck%s.csv", if (what == "ab") "" else paste0("_", what)))
write.csv(out, f, row.names = FALSE)
for (p in unique(out$part)) { o <- out[out$part == p & out$checked, ]; cat(sprintf("G-4 parte %s: %d/%d %s\n", p, sum(o$pass), nrow(o), if (all(o$pass)) "PASS" else "FAIL")) }
