# R/testing/test_S1.5c.R - S1.5c: test della semina a compartimenti (cancello G-1; mutanti G-2). Pre-registrazione: docs/steps/S1.5c.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S15C_RSRC=<cartella con R/ e tools/>] Rscript --vanilla R/testing/test_S1.5c.R
# I mutanti di G-2 sono copie modificate del codice (run_S1.5c.sh mutants), non interruttori nel codice di produzione (BL-155).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
SRC <- Sys.getenv("S15C_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source(file.path(SRC, "tools/S1.5c_bench_lib.R"))
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "t_err"))
isok <- function(o) !inherits(o, "t_err"); msg <- function(o) if (isok(o)) "" else o$err
PRE <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
square <- function(side, x0 = 0, y0 = 0) {
  p <- sf::st_sfc(sf::st_polygon(list(rbind(c(x0, y0), c(x0 + side, y0), c(x0 + side, y0 + side), c(x0, y0 + side), c(x0, y0)))))
  list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = side^2), region_polygons = p)
}
comp1 <- function(ct, f) data.frame(cluster_id = 1L, cell_type = ct, fraction = f)
nn_all <- function(x, y) { d <- as.matrix(dist(cbind(x, y))); diag(d) <- Inf; d }

# ---- T1: campo a caratteristiche di Fourier: media 0, varianza 1, correlazione exp(-r^2 / (2 l^2)) ----
set.seed(1); ell <- 20; N <- 1000
Z <- t(vapply(seq_len(N), function(k) { f <- .sc_rff_draw(ell, 256L); .sc_rff_eval(f, c(0, ell, 3 * ell), c(0, 0, 0)) }, numeric(3)))
r1 <- cor(Z[, 1], Z[, 2]); r3 <- cor(Z[, 1], Z[, 3]); se <- (1 - exp(-0.5)^2) / sqrt(N)
ck("T1a media 0 e varianza 1 del campo (1 000 campi)", abs(mean(Z[, 1])) < 3 / sqrt(N) && abs(var(Z[, 1]) - 1) < 3 * sqrt(2 / N), sprintf("media %.3f var %.3f", mean(Z[, 1]), var(Z[, 1])))
ck("T1b correlazione a distanza l = exp(-1/2)", abs(r1 - exp(-0.5)) < 3 * se, sprintf("%.3f contro %.3f (3 SE %.3f)", r1, exp(-0.5), 3 * se))
ck("T1c correlazione a distanza 3 l = exp(-9/2)", abs(r3 - exp(-4.5)) < 3 / sqrt(N), sprintf("%.3f", r3))

# ---- T2: spostamenti a_t: frazioni attese pesate con la densita' = f ----
set.seed(2); eta <- cbind(log(0.6) + 3 * rnorm(5000), log(0.3) + 2 * rnorm(5000), log(0.1) + rnorm(5000)); f <- c(0.6, 0.3, 0.1); rho <- c(9000, 3000, 1500)
a <- .sc_solve_a(eta, f, rho); pr <- .sc_prob(eta, a, rho); fr <- colSums(pr$w * pr$p) / sum(pr$w)
ck("T2 frazioni attese pesate con lambda = f (1e-6)", max(abs(fr - f)) < 1e-6, sprintf("%.2e", max(abs(fr - f))))

# ---- seme dei test con i preset su un tessuto sintetico ----
sim_preset <- function(A, side, seed, spatial = "compartments", mode = "tissue", labels = "field") {
  sp <- tissue_preset_spec(PRE[[A]], types = "multi", size_arm = "group", spatial = spatial, min_dist_mode = mode)
  seed_centroids(square(side), sp$cell_types, comp1(sp$composition$cell_type, sp$composition$fraction), random_seed = seed,
                 spatial = spatial, compartments = sp$compartments, labels = labels, verbose = FALSE)
}

# ---- T3: distanza minima rispettata da ogni coppia, nei due modi ----
viol <- c(); for (A in c("A2", "A5")) for (mode in c("tissue", "nucleus")) {
  sp <- tissue_preset_spec(PRE[[A]], types = "multi", size_arm = "group", spatial = "compartments", min_dist_mode = mode)
  sc <- sim_preset(A, if (A == "A5") 500 else 250, 11, mode = mode); cdf <- sc$centroids
  d <- nn_all(cdf$x, cdf$y); dm <- sp$cell_types$min_dist_um[match(cdf$cell_type, sp$cell_types$cell_type)]
  thr <- outer(dm, dm, "+") / 2; viol <- c(viol, sum(d < thr - 1e-9))
}
ck("T3 distanza minima della coppia (d_i + d_j)/2 rispettata (A2, A5; tessuto e per nucleo)", all(viol == 0), viol)
# T3b: la distanza della coppia e' la media, non il massimo (due tipi, d = 2 e 6, miscela densa)
ct2 <- data.frame(cell_type = c("p", "q"), density = c(8000, 8000), min_dist_um = c(2, 6))
sc <- seed_centroids(square(200), ct2, comp1(c("p", "q"), c(0.5, 0.5)), random_seed = 3, spatial = "compartments",
                     compartments = data.frame(cell_type = c("p", "q"), ell_um = 20, s = 0), verbose = FALSE)
d <- nn_all(sc$centroids$x, sc$centroids$y); cross <- outer(sc$centroids$cell_type == "p", sc$centroids$cell_type == "q") |
  outer(sc$centroids$cell_type == "q", sc$centroids$cell_type == "p")
ck("T3b coppie miste piu' vicine di 5.4 um (media 4, massimo 6)", min(d[cross]) >= 4 - 1e-9 && min(d[cross]) < 5.4, sprintf("min %.3f", min(d[cross])))

# ---- T4: numero di cellule per regione uguale al modo uniforme (stesso seme) ----
eq_n <- vapply(names(PRE), function(A) { a <- sim_preset(A, 300, 21); b <- sim_preset(A, 300, 21, spatial = "uniform")
  identical(a$region_df$n_target, b$region_df$n_target) && nrow(a$centroids) == nrow(b$centroids) }, TRUE)
ck("T4 numero di cellule per regione uguale al modo uniforme (6 preset)", all(eq_n), sum(eq_n))

# ---- T5: kappa della distanza per nucleo: sum f_t d_t = d* ----
kk <- vapply(names(PRE), function(A) { sp <- tissue_preset_spec(PRE[[A]], types = "multi", size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus")
  f <- sp$composition$fraction[match(sp$cell_types$cell_type, sp$composition$cell_type)]
  abs(sum(f * sp$cell_types$min_dist_um) - .TP_TABLE$min_dist_um[.TP_TABLE$preset == PRE[[A]]]) }, 0)
ck("T5 sum f_t d_t = d* (1e-12, 6 preset)", all(kk < 1e-12), sprintf("%.1e", max(kk)))

# ---- T6: composizione media su 200 semi entro 3 errori Monte Carlo da f ----
sp <- tissue_preset_spec("tumor_microenv", types = "multi", size_arm = "group", spatial = "compartments")
fr <- vapply(1:200, function(s) mean(sim_preset("A2", 150, 1000 + s)$centroids$cell_type == "epiteliale"), 0)
f0 <- sp$composition$fraction[sp$composition$cell_type == "epiteliale"]
ck("T6 frazione dell'epitelio (tumore del colon) entro 3 SE Monte Carlo da f (200 semi)", abs(mean(fr) - f0) < 3 * sd(fr) / sqrt(200), sprintf("%.4f contro %.4f (SE %.4f)", mean(fr), f0, sd(fr) / sqrt(200)))

qr <- vapply(1:20, function(s) { ci <- sim_preset("A2", 300, 3000 + s)$info$compartments; ci$n_quota_reject / ci$n_thin_accept }, 0)
ck("T6b proposte rifiutate per quota esaurita < 5 % in media (tumore del colon, 20 semi): gli spostamenti a_t rendono le quote quasi mai vincolanti", mean(qr) < 0.05, sprintf("media %.4f massimo %.4f", mean(qr), max(qr)))

# ---- T7: con s = 0 la semina a compartimenti e' uguale in distribuzione al modo uniforme ----
ctz <- tissue_preset_spec("epithelial", types = "multi", size_arm = "group")
nnd <- function(sp_mode, s) unlist(lapply(1:30, function(k) { cp <- data.frame(cell_type = ctz$cell_types$cell_type, ell_um = 20, s = 0)
  sc <- seed_centroids(square(120), ctz$cell_types, comp1(ctz$composition$cell_type, ctz$composition$fraction), random_seed = s + k,
                       spatial = sp_mode, compartments = cp, type_order = "blind", verbose = FALSE)
  apply(nn_all(sc$centroids$x, sc$centroids$y), 1, min) }))
ks <- suppressWarnings(ks.test(nnd("compartments", 500), nnd("uniform", 900)))
ck("T7 s = 0: distanza dal vicino uguale in distribuzione al modo uniforme (KS, p > 0.001)", ks$p.value > 0.001, sprintf("p = %.3g", ks$p.value))

# ---- T8: controllo positivo del raggruppamento ----
ctp <- data.frame(cell_type = c("p", "q"), density = c(5000, 5000), min_dist_um = c(0, 0))
same_frac <- function(s) { sc <- seed_centroids(square(400), ctp, comp1(c("p", "q"), c(0.5, 0.5)), random_seed = 8, spatial = "compartments",
  compartments = data.frame(cell_type = c("p", "q"), ell_um = 30, s = s), verbose = FALSE)
  d <- nn_all(sc$centroids$x, sc$centroids$y); j <- apply(d, 1, which.min); mean(sc$centroids$cell_type[j] == sc$centroids$cell_type) }
s4 <- same_frac(4); s0 <- same_frac(0)
ck("T8 s = 4: vicino dello stesso tipo > 0.7; s = 0: 0.45-0.55", s4 > 0.7 && s0 > 0.45 && s0 < 0.55, sprintf("%.3f / %.3f", s4, s0))

# ---- T9: controllo positivo della densita' del tipo nelle zone pure (entro il 10 % di rho_t) ----
ctd <- data.frame(cell_type = c("p", "q"), density = c(4000, 1000), min_dist_um = c(0, 0))
sc <- seed_centroids(square(1500), ctd, comp1(c("p", "q"), c(0.5, 0.5)), random_seed = 9, spatial = "compartments",
                     compartments = data.frame(cell_type = c("p", "q"), ell_um = c(150, 150), s = c(6, 6)), verbose = FALSE)
X <- cbind(sc$centroids$x, sc$centroids$y); ty <- as.character(sc$centroids$cell_type); R9 <- 40
D9 <- as.matrix(dist(X)); diag(D9) <- Inf; nb <- D9 < R9
pure <- vapply(seq_len(nrow(X)), function(i) all(ty[nb[i, ]] == ty[i]) && any(nb[i, ]), TRUE)
inner <- X[, 1] > R9 & X[, 1] < 1500 - R9 & X[, 2] > R9 & X[, 2] < 1500 - R9
dens <- rowSums(nb) / (pi * R9^2) * 1e6                     # vicini entro R (escluso se stesso) / area: stimatore di Palm senza distorsione
dp <- mean(dens[pure & inner & ty == "p"]); dq <- mean(dens[pure & inner & ty == "q"])
ck("T9 densita' nelle zone pure entro il 10 % di rho_t (4 000 e 1 000 /mm2)", abs(dp / 4000 - 1) < 0.1 && abs(dq / 1000 - 1) < 0.1 && sc$info$compartments$n_exceed_wmax == 0,
   sprintf("%.0f / %.0f; superamenti di lambda_max %d", dp, dq, sc$info$compartments$n_exceed_wmax))

# ---- T10: etichette cieche: stesse posizioni e stessi conteggi per tipo ----
a <- sim_preset("A3", 300, 31); b <- sim_preset("A3", 300, 31, labels = "blind")
ck("T10 labels = \"blind\": posizioni identiche, conteggi per tipo identici, etichette diverse",
   identical(a$centroids$x, b$centroids$x) && identical(a$centroids$y, b$centroids$y) && identical(table(a$centroids$cell_type), table(b$centroids$cell_type)) &&
   any(a$centroids$cell_type != b$centroids$cell_type))
e <- safe(sim_preset("A5", 300, 1, mode = "nucleus", labels = "blind"))
ck("T10b etichette cieche con distanze per tipo diverse: errore", !isok(e), msg(e))

# ---- T11: place_cells() sui 6 preset a compartimenti, due modi: nessun ERROR ----
clust_sq <- function(side) { g <- expand.grid(x = 1:side, y = 1:side); data.frame(x = g$x, y = g$y, intensity_cluster = 1L) }
cl <- clust_sq(250); bad <- character(0); exc <- 0
for (A in names(PRE)) for (mode in c("tissue", "nucleus")) {
  L <- safe(suppressWarnings(place_cells(cl, PRE[[A]], pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
             small_components = "exclude", random_seed = 5, types = "multi", size_arm = "group", spatial = "compartments",
             min_dist_mode = mode, strict = FALSE, verbose = FALSE)))
  if (!isok(L)) { bad <- c(bad, paste(A, mode, msg(L))); next }
  if (L$metadata$validation$n_error > 0 || sum(L$metadata$seed_info$n_failed) > 0) bad <- c(bad, paste(A, mode))
  exc <- exc + L$metadata$seed_info$compartments$n_exceed_wmax
}
ck("T11 place_cells() a compartimenti sui 6 preset x 2 distanze: 0 ERROR, 0 non piazzate", !length(bad), bad)
ck("T11b nessun superamento di lambda_max", exc == 0, exc)

# ---- T12: valori del preset = righe B-209...B-231 di BIO_REFERENCES = calibrazione su tutte le ROI ----
bio <- readLines("docs/BIO_REFERENCES.md"); bio <- bio[grepl("^\\| B-2(09|1[0-9]|2[0-9]|3[01]) \\|", bio)]
cells <- function(l) trimws(strsplit(l, "|", fixed = TRUE)[[1]][-1])
kv <- function(s, k) { m <- regmatches(s, regexpr(paste0("(^|; )", k, "=[^;]+"), s)); if (!length(m)) NA_real_ else suppressWarnings(as.numeric(sub(".*=", "", m))) }
B <- do.call(rbind, lapply(bio, function(l) { x <- cells(l); data.frame(id = x[1], arch = x[2], metric = x[3], value = x[4], dist = x[5], stringsAsFactors = FALSE) }))
B$type <- sub(".*\\*\\*(.+)\\*\\*.*", "\\1", B$metric); B$s <- as.numeric(B$value); B$ell <- vapply(B$dist, kv, 0, k = "ell"); B$a_t <- vapply(B$dist, kv, 0, k = "a_t")
ok12 <- nrow(B) == 23; mism <- character(0)
for (k in seq_len(nrow(B))) { r <- .TP_COMPART[.TP_COMPART$src == B$id[k], ]
  if (nrow(r) != 1 || r$cell_type != B$type[k] || abs(r$s - B$s[k]) > 1e-9 || !identical(is.na(r$ell_um), is.na(B$ell[k])) ||
      (!is.na(r$ell_um) && abs(r$ell_um - B$ell[k]) > 1e-9) || abs(r$nuc_median_area_um2 / B$a_t[k] - 1) > 1e-5) mism <- c(mism, B$id[k]) }
cal <- do.call(rbind, lapply(names(PRE), function(A) { d <- read.csv(sprintf("results/S1.5c/S1.5c_compartments_%s.csv", A)); d[d$fold == "pooled", ] }))
for (k in seq_len(nrow(cal))) { r <- .TP_COMPART[.TP_COMPART$preset == PRE[[cal$archetype[k]]] & .TP_COMPART$cell_type == cal$cell_type[k], ]
  if (nrow(r) != 1 || r$s != cal$s[k] || r$ell_um != cal$ell_um[k]) mism <- c(mism, paste(cal$archetype[k], cal$cell_type[k])) }
minor_ok <- all(.TP_COMPART$s[!paste(.TP_COMPART$preset, .TP_COMPART$cell_type) %in% paste(PRE[cal$archetype], cal$cell_type)] == 0)
ck("T12 preset = BIO_REFERENCES (23 righe) = calibrazione su tutte le ROI; tipi minori s = 0", ok12 && !length(mism) && minor_ok, c(nrow(B), mism))

# ---- T13: riproducibilita' con lo stesso seme ----
a <- sim_preset("A1", 200, 77); b <- sim_preset("A1", 200, 77)
ck("T13 stesso seme, stessi centri e tipi", identical(a$centroids, b$centroids))

# ---- T14: parametri della piega nel banco: lascia-fuori del ROI simulato ----
rois <- read.csv("results/R2/R2_rois_checked.csv"); okf <- 0; diff_pooled <- 0
for (k in seq_len(nrow(rois))) { A <- rois$archetype[k]; r <- rois$roi_id[k]
  cp <- compart_for(A, r); d <- read.csv(sprintf("results/S1.5c/S1.5c_compartments_%s.csv", A)); e <- d[d$fold == paste0("loo_", r), ]
  m <- match(e$cell_type, cp$cell_type); okf <- okf + (all(!is.na(m)) && all(cp$s[m] == e$s) && all(cp$ell_um[m] == e$ell_um))
  pl <- d[d$fold == "pooled", ]; diff_pooled <- diff_pooled + any(pl$s != e$s | pl$ell_um != e$ell_um) }
ck("T14 il banco legge i parametri lascia-fuori del ROI (30/30; diversi da tutte le ROI in almeno 20)", okf == 30 && diff_pooled >= 20, sprintf("%d / %d", okf, diff_pooled))

# ---- T15: errori attesi ----
e1 <- safe(tissue_preset_spec("lymphoid", types = "multi", spatial = "uniform", min_dist_mode = "nucleus"))
e2 <- safe(tissue_preset_spec("lymphoid", types = "single", spatial = "compartments"))
e3 <- safe(seed_centroids(square(100), ctp, comp1(c("p", "q"), c(0.5, 0.5)), spatial = "compartments", compartments = data.frame(cell_type = "p", ell_um = -1, s = 1), verbose = FALSE))
ck("T15 errori: distanza per nucleo senza compartimenti; compartimenti a un tipo; ell negativa", !isok(e1) && !isok(e2) && !isok(e3))
ck("T16 predefinito della semina dei preset = uniforme (prima della regola 3)", identical(tissue_preset_spec("epithelial")$spatial, "uniform"))

n <- length(res); np <- sum(unlist(res)); cat(sprintf("RESULT %d/%d PASS\n", np, n))
