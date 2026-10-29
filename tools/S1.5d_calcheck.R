# tools/S1.5d_calcheck.R - S1.5d, cancello G-4 (docs/steps/S1.5d.md).
# (a) Etichette del pacchetto (formule di R/04b2_seed_centroids.R, generatore di R, densita' uguali come la calibrazione) sulle posizioni dei nuclei veri con
#     i parametri di tutte le ROI, 20 repliche, contro le 20 repliche di tools/S1.5d_calibration_types.py: z = differenza dei logaritmi delle medie / errore
#     standard della differenza (metodo delta dalle repliche), |z| <= 4 nei (ROI, tipo, bin) con almeno 200 coppie vere dello stesso tipo. Focali: celle
#     interne del tipo (al piu' 1 500), campione di R (non quello di numpy: rumore in piu', stessa attesa).
# (b) Layer del banco (seme 1): compartimenti e eterogeneita' della piega lascia-fuori della propria ROI (30/30 per braccio).
# (c) Tabella dei residui: (c1) come scritto nella pre-registrazione, quantili «di tools/S1.5d_size_emp.py» (residui della retta OLS sui nuclei di R6, decimi
#     e quantili di numpy) <= 1e-9 (errore di pre-registrazione D-4: quantita' diverse); (c2) seconda implementazione in R dalle stesse celle e dagli stessi
#     parametri del preset (quantili di tipo 7 = numpy lineare) <= 1e-9.
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5d_calcheck.R [a|b|c|abc]   Uscita: results/S1.5d/S1.5d_G4_calcheck[_x].csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow); library(spatstat.geom) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/S1.5d_bench_lib.R")
options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")   # S1.6: predefiniti di S1.5c-S1.5e, per riprodurre i numeri di questo step
what <- commandArgs(TRUE); what <- if (length(what)) what[1] else "abc"
RES <- "results/S1.5d"; M <- "/mnt/micron/geo_spatialtrans"
UNTYPED <- c("non risolto", "non classificabile"); EPI <- c("epiteliale normale" = "epiteliale", "epiteliale tumorale" = "epiteliale")
BINS <- c(0, 10, 20, 40, 80, 160); NB <- length(BINS) - 1; NFOC <- 1500L; NREP <- 20L
MAIN <- list(A1 = c("epiteliale", "stromale e muscolo liscio", "immunitario"), A2 = c("epiteliale", "fibroblasti"),
             A3 = c("epiteliale", "fibroblasti", "endoteliale", "immunitario"), A4 = c("linfociti B", "linfociti T"),
             A5 = c("neuroni", "oligodendrociti", "astrociti"), A6 = c("cardiomiociti", "fibroblasti", "immunitari"))
PT <- read.csv(file.path(RES, "S1.5d_preset_types.csv"), stringsAsFactors = FALSE)
load_terr <- function(A) {                                            # stesso insieme di tools/S1.5b_calibration.py load(): territori R6 + etichette PH-1
  lab <- as.data.frame(read_parquet(sprintf("%s/R6/labels/%s_exp_umi20_seed20261022_nuclei_ph1.parquet", M, A)))
  lab$type <- ifelse(lab$cell_type %in% names(EPI), EPI[lab$cell_type], lab$cell_type)
  do.call(rbind, lapply(sort(unique(lab$roi_id)), function(r) { T <- as.data.frame(read_parquet(sprintf("%s/R6/terr/%s_%s.parquet", M, A, r)))
    L <- lab[lab$roi_id == r, ]; k <- match(paste(round(T$x, 6), round(T$y, 6)), paste(round(L$x_um, 6), round(L$y_um, 6)))
    stopifnot(!anyNA(k)); data.frame(roi_id = r, x = T$x, y = T$y, interior = T$interior, type = L$type[k], a_t = T$a_t, a_n = T$a_n) }))
}
rows <- list()
if (grepl("a", what)) for (A in paste0("A", 1:6)) {
  d <- load_terr(A); pt <- PT[PT$archetype == A, ]; types <- sort(pt$cell_type); d <- d[d$type %in% types, ]
  f <- pt$fraction[match(types, pt$cell_type)]; f <- f / sum(f); K <- length(types); main <- match(MAIN[[A]], types)
  cp <- compart_for_d(A, NA, "pooled"); s <- cp$s[match(types, cp$cell_type)]; s[is.na(s)] <- 0; ell <- cp$ell_um[match(types, cp$cell_type)]
  py <- read.csv(sprintf("%s/S1.5d_compartments_curves_%s.csv", RES, A), stringsAsFactors = FALSE)
  set.seed(20261025)
  for (r in sort(unique(d$roi_id))) {
    e <- d[d$roi_id == r, ]; n <- nrow(e); Lr <- match(e$type, types)
    W <- owin(range(e$x) + c(-1, 1), range(e$y) + c(-1, 1)); Y <- ppp(e$x, e$y, window = W)
    ratio <- function(Lv) { foc <- unlist(lapply(main, function(t) { k <- which(e$interior & Lv == t); if (length(k) < 20) integer(0) else if (length(k) <= NFOC) k else sort(sample(k, NFOC)) }))
      X <- ppp(e$x[foc], e$y[foc], window = W); cpp <- crosspairs(X, Y, rmax = 160, what = "ijd"); fi <- foc[cpp$i]; j <- cpp$j; k <- fi != j
      fi <- fi[k]; j <- j[k]; b <- findInterval(cpp$d[k], BINS, rightmost.closed = TRUE); b[b > NB] <- NB
      idx <- (b - 1) * K + Lv[fi]; same <- Lv[fi] == Lv[j]
      num <- matrix(tabulate(idx[same], K * NB), K); den <- matrix(tabulate(idx, K * NB), K); fr <- tabulate(Lv, K) / n
      list(R = num / pmax(den, 1) / fr, num = num) }
    Rs <- array(NA_real_, c(NREP, K, NB))
    for (rep in seq_len(NREP)) { fields <- lapply(seq_len(K), function(t) if (s[t] > 0) .sc_rff_draw(ell[t], 256L) else NULL)
      eta <- .sc_eta(fields, s, log(f), e$x, e$y); a <- .sc_solve_a(eta, f, rep(1, K)); P <- .sc_prob(eta, a, rep(1, K))$p
      Rs[rep, , ] <- ratio(.sc_draw_type(P, runif(n)))$R }
    real <- ratio(Lr)
    for (t in main) for (kb in seq_len(NB)) { bn <- sprintf("%d-%d", BINS[kb], BINS[kb + 1])
      p <- py[py$roi_id == r & py$cell_type == types[t] & py$bin == bn, ]; xr <- Rs[, t, kb]; xp <- p$sim
      se <- sqrt(var(xr) / mean(xr)^2 / NREP + var(xp) / mean(xp)^2 / length(xp)); z <- (log(mean(xr)) - log(mean(xp))) / se
      chk <- p$n_pairs_real[1] >= 200
      rows[[length(rows) + 1]] <- data.frame(part = "a", archetype = A, roi_id = r, cell_type = types[t], bin = bn, n_pairs_real = p$n_pairs_real[1],
        R_real_py = p$real[1], R_sim_R = mean(xr), R_sim_py = mean(xp), z = z, value = z, tol = 4, checked = chk, pass = !chk || (is.finite(z) && abs(z) <= 4)) }
  }
  cat(A, "fatto\n")
}
if (grepl("b", what)) {
  rois <- read.csv("results/R2/R2_rois_checked.csv")
  for (arm in c("R", "RV", "RD", "RVD", "VF", "W3", "W3p", "Q")) for (k in seq_len(nrow(rois))) {
    A <- rois$archetype[k]; r <- rois$roi_id[k]; f <- sprintf("%s/S1.5d/layers/%s_%s_%s.rds", M, A, r, arm); if (!file.exists(f)) next
    L <- readRDS(f); got <- L$metadata$compartments                      # attesi letti dai CSV, senza le funzioni del banco (revisione RC-06)
    cc <- read.csv(sprintf("%s/S1.5d_compartments_%s.csv", RES, A), stringsAsFactors = FALSE); exp <- cc[cc$fold == paste0("loo_", r), ]
    m <- match(exp$cell_type, got$cell_type)
    ok_c <- nrow(exp) > 0 && all(!is.na(m)) && all(got$s[m] == exp$s) && all(got$ell_um[m] == exp$ell_um)
    fam <- c(RV = "RV", RD = "RD", RVD = "RVD", VF = "RV")[arm]; if (arm %in% c("W3", "W3p", "Q")) fam <- read.csv(file.path(RES, "S1.5d_ruleA.csv"))$winner[1]
    h <- L$metadata$heterogeneity
    if (is.null(h)) h <- data.frame(void_frac = 0, void_ell_um = NA_real_, dens_sigma = 0, dens_ell_um = NA_real_)   # NULL = vuoti e densita' spenti (nota D-11)
    if (arm == "R") ok_h <- h$void_frac == 0 && h$dens_sigma == 0 else {
      hh <- read.csv(sprintf("%s/S1.5d_het_%s.csv", RES, fam), stringsAsFactors = FALSE); e <- hh[hh$archetype == A & hh$fold == paste0("loo_", r), ]
      if (arm == "VF") e$void_ell_um <- 2
      ok_h <- nrow(e) == 1 && h$void_frac == e$void_frac && h$dens_sigma == e$dens_sigma &&
        (e$void_frac == 0 || h$void_ell_um == e$void_ell_um) && (e$dens_sigma == 0 || h$dens_ell_um == e$dens_ell_um) }
    rows[[length(rows) + 1]] <- data.frame(part = "b", archetype = A, roi_id = r, cell_type = arm, bin = "", n_pairs_real = NA, R_real_py = NA, R_sim_R = NA, R_sim_py = NA,
      z = NA, value = as.numeric(ok_c && ok_h), tol = 1, checked = TRUE, pass = ok_c && ok_h)
  }
}
if (grepl("c", what)) {
  tab <- read.csv("inst/extdata/size_resid.csv", stringsAsFactors = FALSE); P <- round(seq(0, 1, 0.01), 2)
  for (key in unique(tab$key)) { pre <- sub(":.*", "", key); g <- as.integer(sub(".*:", "", key)); ty <- PT[PT$preset == pre & PT$size_group == g, ][1, ]
    cc <- as.data.frame(read_parquet(sprintf("%s/S1.5d/size_resid_cells/%s_%d.parquet", M, pre, g)))
    x <- log(cc$a_t); rr <- log(cc$a_n) - (2 * log(ty$rho_g) + (1 - ty$beta_g) * log(1e6 / ty$density_model) + ty$beta_g * x)
    ed <- quantile(x, seq(0.1, 0.9, 0.1), type = 7, names = FALSE); dc <- findInterval(x, ed, left.open = TRUE) + 1L
    q2 <- c(unlist(lapply(1:10, function(k) quantile(rr[dc == k], P, type = 7, names = FALSE))), quantile(rr, P, type = 7, names = FALSE))
    tk <- tab[tab$key == key, ]; tk <- tk[order(tk$decile == 0, tk$decile, tk$p), ]; d2 <- max(abs(q2 - tk$r))
    rows[[length(rows) + 1]] <- data.frame(part = "c2", archetype = pre, roi_id = "", cell_type = key, bin = "", n_pairs_real = NA, R_real_py = NA, R_sim_R = NA, R_sim_py = NA,
      z = NA, value = d2, tol = 1e-9, checked = TRUE, pass = d2 <= 1e-9)
    A <- ty$archetype; d <- load_terr(A); d <- d[d$interior & d$type %in% PT$cell_type[PT$preset == pre & PT$size_group == g], ]   # c1: come scritto (nuclei di R6, retta OLS)
    x1 <- log(d$a_t); y1 <- log(d$a_n); b <- coef(lm(y1 ~ x1)); r1 <- y1 - b[1] - b[2] * x1
    e1 <- quantile(x1, seq(0.1, 0.9, 0.1), type = 7, names = FALSE); d1c <- findInterval(x1, e1, left.open = TRUE) + 1L
    q1 <- c(unlist(lapply(1:10, function(k) quantile(r1[d1c == k], P, type = 7, names = FALSE))), quantile(r1, P, type = 7, names = FALSE)); d1 <- max(abs(q1 - tk$r))
    rows[[length(rows) + 1]] <- data.frame(part = "c1", archetype = pre, roi_id = "", cell_type = key, bin = "", n_pairs_real = NA, R_real_py = NA, R_sim_R = NA, R_sim_py = NA,
      z = NA, value = d1, tol = 1e-9, checked = TRUE, pass = d1 <= 1e-9)
  }
}
out <- do.call(rbind, rows); fo <- file.path(RES, sprintf("S1.5d_G4_calcheck%s.csv", if (what == "abc") "" else paste0("_", what))); write.csv(out, fo, row.names = FALSE)
for (p in unique(out$part)) { o <- out[out$part == p & out$checked, ]; cat(sprintf("G-4 parte %s: %d/%d %s\n", p, sum(o$pass), nrow(o), if (all(o$pass)) "PASS" else "FAIL")) }
