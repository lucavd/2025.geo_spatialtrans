# R/testing/test_S1.5d.R - S1.5d: test dell'eterogeneita' della semina e dei residui empirici del modello di taglia (cancello G-1; mutanti G-2).
# Pre-registrazione: docs/steps/S1.5d.md. Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S15D_RSRC=<cartella con R/ e tools/>]
#   Rscript --vanilla R/testing/test_S1.5d.R      I mutanti sono copie modificate del codice (run_S1.5d.sh mutants).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow) })
SRC <- Sys.getenv("S15D_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source(file.path(SRC, "tools/S1.5d_bench_lib.R")); source("tools/R3_voronoi_metrics.R")
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "t_err"))
isok <- function(o) !inherits(o, "t_err"); msg <- function(o) if (isok(o)) "" else o$err
PRE <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
square <- function(side) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = side^2),
  region_polygons = sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(side, 0), c(side, side), c(0, side), c(0, 0))))))
spec_c <- function(A) tissue_preset_spec(PRE[[A]], size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus")
comp_of <- function(sp) data.frame(cluster_id = 1L, cell_type = sp$composition$cell_type, fraction = sp$composition$fraction)
seed_h <- function(A, side, seed, het, quotas = TRUE) { sp <- spec_c(A)
  seed_centroids(square(side), sp$cell_types, comp_of(sp), random_seed = seed, spatial = "compartments", compartments = sp$compartments,
                 heterogeneity = het, quotas = quotas, verbose = FALSE) }
H <- function(pv = 0, lv = NA, sd = 0, ld = NA) data.frame(void_frac = pv, void_ell_um = lv, dens_sigma = sd, dens_ell_um = ld)
hrel3 <- function(x, y, side) { g <- expand.grid(x = seq(0.5, side - 0.5, 1), y = seq(0.5, side - 0.5, 1))
  d <- FNN::get.knnx(cbind(x, y), g, k = 1)$nn.dist[, 1]; nn <- FNN::get.knn(cbind(x, y), k = 1)$nn.dist[, 1]; mean(d / median(nn) > 3) }

# ---- T1: nessun centro dove il campo dei vuoti supera la soglia ----
v1 <- seed_h("A2", 400, 5, H(0.15, 20)); ri <- v1$info$compartments$regions[[1]]; tau <- qnorm(1 - 0.15)
Vc <- .sc_rff_eval(ri$field_void, v1$centroids$x, v1$centroids$y)
ck("T1 nessun centro con V > tau (tumore del colon, vuoti 15 %, 20 um)", length(Vc) > 0 && max(Vc) <= tau, sprintf("max V %.3f, tau %.3f, n %d", max(Vc), tau, length(Vc)))
# ---- T2: frazione dei vuoti sulla griglia = void_frac (200 campi, 3 errori Monte Carlo) ----
set.seed(2); g <- expand.grid(x = seq(1, 999, 4), y = seq(1, 999, 4))
fr <- vapply(1:200, function(k) mean(.sc_rff_eval(.sc_rff_draw(30, 256L), g$x, g$y) > tau), 0)
ck("T2 frazione dei vuoti = 0.15 entro 3 errori Monte Carlo (200 campi)", abs(mean(fr) - 0.15) < 3 * sd(fr) / sqrt(200), sprintf("%.4f (3 ES %.4f)", mean(fr), 3 * sd(fr) / sqrt(200)))
# ---- T3: controllo positivo della densita': rapporto fra zone alte e basse di W = rapporto atteso entro il 10 % ----
ct3 <- data.frame(cell_type = "p", density = 800, min_dist_um = 0); cp3 <- data.frame(cluster_id = 1L, cell_type = "p", fraction = 1)
s3 <- seed_centroids(square(3000), ct3, cp3, random_seed = 3, spatial = "compartments", compartments = data.frame(cell_type = "p", ell_um = 20, s = 0),
                     heterogeneity = H(0, NA, 0.6, 40), verbose = FALSE)
fw <- s3$info$compartments$regions[[1]]$field_dens; gg <- expand.grid(x = seq(3, 2997, 6), y = seq(3, 2997, 6))
Wg <- .sc_rff_eval(fw, gg$x, gg$y); Wc <- .sc_rff_eval(fw, s3$centroids$x, s3$centroids$y); q <- quantile(Wg, c(1 / 3, 2 / 3))
dens_ratio <- (sum(Wc > q[2]) / mean(Wg > q[2])) / (sum(Wc < q[1]) / mean(Wg < q[1]))
exp_ratio <- mean(exp(0.6 * Wg[Wg > q[2]])) / mean(exp(0.6 * Wg[Wg < q[1]]))
ck("T3 densita' zone alte / basse di W = atteso entro il 10 %", abs(dens_ratio / exp_ratio - 1) < 0.10, sprintf("%.3f contro %.3f", dens_ratio, exp_ratio))
# ---- T4: numero di cellule per regione uguale al caso senza eterogeneita' ----
eq4 <- vapply(names(PRE), function(A) { a <- seed_h(A, 300, 21, NULL); b <- seed_h(A, 300, 21, H(0.1, 16, 0.3, 40))
  identical(a$region_df$n_target, b$region_df$n_target) && nrow(a$centroids) == nrow(b$centroids) && b$info$n_failed == 0 }, TRUE)
ck("T4 numero di cellule per regione invariato con vuoti e densita' (6 preset)", all(eq4), names(eq4)[!eq4])
# ---- T5: distanza minima rispettata con l'eterogeneita' ----
viol <- vapply(c("A2", "A4", "A5"), function(A) { sp <- spec_c(A); s <- seed_h(A, if (A == "A5") 500 else 250, 7, H(0.15, 16, 0.6, 10))
  d <- as.matrix(dist(s$centroids[, c("x", "y")])); diag(d) <- Inf; dm <- sp$cell_types$min_dist_um[match(s$centroids$cell_type, sp$cell_types$cell_type)]
  sum(d < outer(dm, dm, "+") / 2 - 1e-9) }, 0)
ck("T5 distanza minima della coppia rispettata (A2, A4, A5, vuoti e densita')", all(viol == 0), viol)
# ---- T6: frazioni attese pesate con lambda e con l'eterogeneita' = f (1e-6) ----
set.seed(6); eta <- cbind(log(0.6) + 3 * rnorm(5000), log(0.3) + 2 * rnorm(5000), log(0.1) + rnorm(5000)); f <- c(0.6, 0.3, 0.1); rho <- c(9000, 3000, 1500)
h6 <- exp(0.6 * rnorm(5000)) * (runif(5000) > 0.2); k6 <- h6 > 0
a <- .sc_solve_a(eta[k6, ], f, rho, h = h6[k6]); pr <- .sc_prob(eta[k6, ], a, rho); fr6 <- colSums(h6[k6] * pr$w * pr$p) / sum(h6[k6] * pr$w)
ck("T6 frazioni attese pesate = f con i pesi dell'eterogeneita' (1e-6)", max(abs(fr6 - f)) < 1e-6, sprintf("%.2e", max(abs(fr6 - f))))
# ---- T7: eterogeneita' NULL, a zero e del preset non calibrato identiche; layer identici a S1.5c (3 preset, seme 1) ----
same7 <- vapply(c("A1", "A3", "A6"), function(A) { a <- seed_h(A, 300, 9, NULL); b <- seed_h(A, 300, 9, H(0, NA, 0, NA)); identical(a$centroids, b$centroids) }, TRUE)
ck("T7a heterogeneity NULL = tutto a zero (bit per bit, 3 preset; con T7b = S1.5c)", all(same7), names(same7)[!same7])   # nota D-12: tolto «= preset non calibrato»
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
U5 <- read.csv("results/S1.5c/S1.5c_panel_units.csv", stringsAsFactors = FALSE)
clust_of <- function(w) { idx <- which(w$m, arr.ind = TRUE); data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L) }
layer_md5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  on.exit(unlink(f)); unname(tools::md5sum(f)) }
source("tools/S1.5c_bench_lib.R")
md7 <- vapply(list(c("A3", "r1", "N-a"), c("A6", "r1", "N-a"), c("A5", "r2", "U")), function(z) { A <- z[1]; roi <- z[2]; arm <- z[3]
  i <- which(rois$archetype == A & rois$roi_id == roi); w <- r3_roi_window(A, roi, rois)
  L <- if (arm == "U") suppressWarnings(place_cells(clust_of(w), PRE[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
            small_components = "exclude", random_seed = 20261021 + 1000 * i + 1, types = "multi", size_arm = "recentred", spatial = "uniform",
            keep_tiles = TRUE, strict = FALSE, verbose = FALSE)) else
       suppressWarnings(place_cells(clust_of(w), PRE[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
            small_components = "exclude", random_seed = 20261021 + 1000 * i + 1, types = "multi", size_arm = "group", spatial = "compartments",
            min_dist_mode = "nucleus", override_compartments = compart_for(A, roi), override_heterogeneity = H(0, NA, 0, NA), keep_tiles = TRUE, strict = FALSE, verbose = FALSE))
  identical(layer_md5(L), U5$md5[U5$archetype == A & U5$roi_id == roi & U5$seed == 1 & U5$arm == arm]) }, TRUE)
ck("T7b place_cells() predefinito = layer di S1.5c (md5; stroma e cuore a compartimenti, cervello uniforme; seme 1)", all(md7), md7)
# ---- T8: lo spazio vuoto cresce con void_frac (5 valori, 20 semi) ----
h8 <- vapply(c(0, 0.05, 0.10, 0.15, 0.20), function(pv) mean(vapply(1:20, function(sd) { s <- seed_h("A1", 200, 100 + sd, if (pv > 0) H(pv, 16) else NULL)
  hrel3(s$centroids$x, s$centroids$y, 200) }, 0)), 0)
ck("T8 H_rel_3 cresce con void_frac (intestino, 0...0.20, 20 semi)", all(diff(h8) > 0), sprintf("%.4f", h8))
# ---- T9: superamenti di lambda_max <= 1e-3 delle proposte (densita' fine: l_d = 10 um su 1 mm^2) ----
s9 <- seed_h("A2", 1000, 9, H(0.1, 16, 0.6, 10)); ex <- s9$info$compartments$n_exceed_wmax; np <- s9$info$n_attempts
ck("T9 superamenti di lambda_max <= 1e-3 delle proposte (l_d = 10 um, 1 mm^2)", ex / np <= 1e-3, sprintf("%d / %d", ex, np))
# ---- T10: quotas = FALSE senza rifiuti per quota ----
q1 <- seed_h("A4", 400, 10, NULL, quotas = TRUE); q0 <- seed_h("A4", 400, 10, NULL, quotas = FALSE)
ck("T10 quotas = FALSE: 0 rifiuti per quota (con quote > 0)", q0$info$compartments$n_quota_reject == 0 && q1$info$compartments$n_quota_reject > 0,
   sprintf("%d contro %d", q0$info$compartments$n_quota_reject, q1$info$compartments$n_quota_reject))
# ---- T11: place_cells() sui 6 preset con eterogeneita' e residui empirici, senza ERROR; l'eterogeneita' arriva alla semina ----
e11 <- vapply(names(PRE), function(A) { cl <- expand.grid(x = 1:260, y = 1:260); cl$intensity_cluster <- 1L
  L <- safe(suppressWarnings(place_cells(cl, PRE[[A]], pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L, random_seed = 11,
            size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus", override_heterogeneity = H(0.15, 16, 0.3, 40), size_resid = "empirical", verbose = FALSE)))
  isok(L) && L$metadata$validation$n_error == 0 && !is.null(L$metadata$heterogeneity) && L$metadata$size_resid == "empirical" }, TRUE)
ck("T11 place_cells() con eterogeneita' e residui empirici sui 6 preset: 0 ERROR, parametri registrati", all(e11), names(e11)[!e11])
cl <- expand.grid(x = 1:300, y = 1:300); cl$intensity_cluster <- 1L
pc11 <- function(h) suppressWarnings(place_cells(cl, "epithelial", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L, random_seed = 12,
  size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus", override_heterogeneity = h, validate = FALSE, verbose = FALSE))
a11 <- pc11(H(0, NA, 0, NA)); b11 <- pc11(H(0.2, 16))   # nota D-12: confronto con vuoti a zero, non con il preset
ck("T11b override_heterogeneity cambia la semina di place_cells() (H_rel_3 piu' alto)", hrel3(b11$cell_df$x, b11$cell_df$y, 300) > 3 * hrel3(a11$cell_df$x, a11$cell_df$y, 300),
   sprintf("%.4f contro %.4f", hrel3(b11$cell_df$x, b11$cell_df$y, 300), hrel3(a11$cell_df$x, a11$cell_df$y, 300)))
# ---- T12: valori dei preset uguali alle righe di BIO_REFERENCES (eterogeneita' e tabella dei residui) ----
bio <- readLines("docs/BIO_REFERENCES.md"); bio <- bio[grepl("^\\| B-[0-9]{3} \\|", bio)]
bval <- function(id) { r <- bio[startsWith(bio, paste0("| ", id, " |"))]; if (!length(r)) return(NA_character_); trimws(strsplit(r, "\\|")[[1]][5]) }
het_ok <- all(vapply(seq_len(nrow(.TP_HET)), function(k) { s <- .TP_HET$src[k]; if (!grepl("^B-[0-9]{3}$", s)) return(FALSE)
  v <- bval(s); h <- .TP_HET[k, ]; identical(v, sprintf("%g / %g / %g / %g", h$void_frac, h$void_ell_um, h$dens_sigma, h$dens_ell_um)) }, TRUE))
md5r <- unname(tools::md5sum(file.path("inst", "extdata", .TP_SIZE_RESID_FILE))); rb <- bio[grepl("size_resid.csv", bio, fixed = TRUE)]
ck("T12 eterogeneita' dei 6 preset e md5 della tabella dei residui uguali alle righe di BIO_REFERENCES", het_ok && length(rb) >= 1 && all(grepl(md5r, rb, fixed = TRUE)),
   sprintf("eterogeneita' %s, righe della tabella %d", het_ok, length(rb)))
# ---- T13: residui del pacchetto (via .dcg_target_area) contro i residui veri, per gruppo e decimo: Kolmogorov <= 0.02 ----
tab <- .tp_resid_table(); U <- (seq_len(200) - 0.5) / 200; worst <- c(); exact <- c()
for (A in names(PRE)) { sp <- tissue_preset_spec(PRE[[A]], size_arm = "group", size_resid = "empirical"); ct <- sp$cell_types; p <- .TP_TABLE[.TP_TABLE$preset == PRE[[A]], ]
  for (g in sort(unique(ct$size_group))) { key <- paste0(PRE[[A]], ":", g); ti <- which(ct$size_group == g)[1]; ty <- .TP_TYPES[.TP_TYPES$preset == PRE[[A]] & .TP_TYPES$size_group == g, ][1, ]
    cc <- as.data.frame(read_parquet(sprintf("/mnt/micron/geo_spatialtrans/S1.5d/size_resid_cells/%s_%d.parquet", PRE[[A]], g)))
    rr <- log(cc$a_n) - (2 * log(ty$rho_g) + (1 - ty$beta_g) * log(1e6 / p$density_model) + ty$beta_g * log(cc$a_t))
    tk <- tab[tab$key == key & tab$decile > 0, ]; hi <- tapply(tk$x_hi, tk$decile, `[`, 1); dc <- findInterval(log(cc$a_t), hi[-10], left.open = TRUE) + 1L
    for (d in 1:10) { k <- which(dc == d); if (length(k) > 500) k <- k[round(seq(1, length(k), length.out = 500))]
      at <- rep(cc$a_t[k], each = length(U)); z <- qnorm(rep(U, length(k)))
      a_tgt <- .dcg_target_area(at, ct, rep(ti, length(at)), z, tab, "decile")
      la <- log(ct$nucleus_to_eq_ratio[ti]^2 * 1e6 / ct$density[ti]) + ct$nucleus_area_slope[ti] * (log(at) - log(1e6 / ct$density[ti]))
      eps <- log(a_tgt) - la; rd <- rr[dc == d]
      qd <- tk$r[tk$decile == d][order(tk$p[tk$decile == d])]; eps_ref <- rep(approx(seq(0, 1, 0.01), qd, U)$y, length(k)); exact[paste(key, d)] <- max(abs(eps - eps_ref))
      xs <- sort(unique(c(rd, eps))); D <- max(abs(ecdf(rd)(xs) - ecdf(eps)(xs))); worst[paste(key, d)] <- D } } }
ck("T13 residui del pacchetto = residui veri per gruppo e decimo (Kolmogorov <= 0.02, 10 gruppi x 10 decimi)", max(worst) <= 0.02, sprintf("max %.4f (%s)", max(worst), names(which.max(worst))))
# T13b (aggiunto prima della corsa ufficiale, nota D-6): residuo del pacchetto = interpolazione della tabella del decimo della cella, cella per cella
ck("T13b residui del pacchetto = quantili interpolati del decimo giusto (1e-9, 10 gruppi x 10 decimi)", max(exact) <= 1e-9, sprintf("max %.2e", max(exact)))
# ---- T14: senza tabella identico bit per bit; con la tabella i bersagli cambiano ----
sq14 <- square(200); sp14 <- spec_c("A2"); sc14 <- seed_centroids(sq14, sp14$cell_types, comp_of(sp14), random_seed = 14, spatial = "compartments", compartments = sp14$compartments, verbose = FALSE)
tv14 <- tessellate_voronoi(sc14$centroids, sq14, verbose = FALSE)
d0 <- derive_cell_geometry(sc14$centroids, tv14, sp14$cell_types, random_seed = 15, verbose = FALSE)
d1 <- derive_cell_geometry(sc14$centroids, tv14, sp14$cell_types, random_seed = 15, verbose = FALSE, resid_table = tab)
ct_e <- tissue_preset_spec("tumor_microenv", size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus", size_resid = "empirical")$cell_types
d2 <- derive_cell_geometry(sc14$centroids, tv14, ct_e, random_seed = 15, verbose = FALSE, resid_table = tab)
ck("T14a tabella data ma catalogo senza nucleus_area_resid: identico bit per bit", identical(d0$cell_df, d1$cell_df) && identical(sf::st_as_binary(d0$cell_nuclei), sf::st_as_binary(d1$cell_nuclei)))
ck("T14b con i residui empirici i bersagli cambiano in almeno il 99 % delle cellule", mean(d2$cell_df$nucleus_area_target != d0$cell_df$nucleus_area_target) >= 0.99,
   sprintf("%.4f", mean(d2$cell_df$nucleus_area_target != d0$cell_df$nucleus_area_target)))
# ---- T15: celle fuori dall'intervallo vero nel decimo estremo; T16: modo pooled = decimo 0 ----
tk <- tab[tab$key == "tumor_microenv:1", ]; q1 <- tk$r[tk$decile == 1]; q10 <- tk$r[tk$decile == 10]; q0 <- tk$r[tk$decile == 0]
uu <- c(0.005, 0.5, 0.995); lo <- .dcg_resid_draw(tab, rep("tumor_microenv:1", 3), rep(-50, 3), uu, "decile"); hi2 <- .dcg_resid_draw(tab, rep("tumor_microenv:1", 3), rep(50, 3), uu, "decile")
po <- .dcg_resid_draw(tab, rep("tumor_microenv:1", 3), rep(4, 3), uu, "pooled")
ip <- function(q, u) approx(seq(0, 1, 0.01), q, u)$y
ck("T15 log A_T fuori dall'intervallo vero: decimo 1 sotto, decimo 10 sopra", max(abs(lo - ip(q1, uu))) < 1e-12 && max(abs(hi2 - ip(q10, uu))) < 1e-12)
ck("T16 modo pooled: quantili del decimo 0", max(abs(po - ip(q0, uu))) < 1e-12)
# ---- T17: correlazione dei campi dei vuoti e della densita' a distanza l = exp(-1/2) (campi registrati dalla semina) ----
cr <- t(vapply(1:300, function(k) { s <- seed_centroids(square(60), ct3, cp3, random_seed = 1000 + k, spatial = "compartments",
     compartments = data.frame(cell_type = "p", ell_um = 20, s = 0), heterogeneity = H(0.1, 15, 0.3, 25), verbose = FALSE)$info$compartments$regions[[1]]
  c(.sc_rff_eval(s$field_void, c(0, 15), c(0, 0)), .sc_rff_eval(s$field_dens, c(0, 25), c(0, 0))) }, numeric(4)))
se17 <- (1 - exp(-1)) / sqrt(300)
ck("T17 correlazione a distanza l del campo dei vuoti e della densita' = exp(-1/2) (300 regioni)", abs(cor(cr[, 1], cr[, 2]) - exp(-0.5)) < 3 * se17 && abs(cor(cr[, 3], cr[, 4]) - exp(-0.5)) < 3 * se17,
   sprintf("%.3f, %.3f (3 ES %.3f)", cor(cr[, 1], cr[, 2]), cor(cr[, 3], cr[, 4]), 3 * se17))
# ---- T18: il banco legge i parametri della piega lascia-fuori (compartimenti ed eterogeneita') ----
f18 <- vapply(c("A1", "A4"), function(A) isTRUE(safe({ d <- read.csv(sprintf("results/S1.5d/S1.5d_compartments_%s.csv", A)); e <- compart_for_d(A, if (A == "A4") "f1" else "r1")
  h <- read.csv("results/S1.5d/S1.5d_het_RVD.csv"); hh <- het_for(A, if (A == "A4") "f1" else "r1", "RVD")
  key <- paste0("loo_", if (A == "A4") "f1" else "r1"); all(e$s == d$s[d$fold == key]) && isTRUE(all.equal(hh$void_frac, h$void_frac[h$archetype == A & h$fold == key])) })), TRUE)
ck("T18 compart_for_d() e het_for() leggono la piega lascia-fuori della ROI", all(f18), names(f18)[!f18])

n <- length(res); np <- sum(unlist(res)); cat(sprintf("RESULT %d/%d PASS\n", np, n))
