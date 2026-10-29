# R/testing/test_S1.5e.R - S1.5e: test dei vuoti che tolgono cellule, dei vuoti anisotropi e del tetto della densita' (cancello G-1; mutanti G-2).
# Pre-registrazione: docs/steps/S1.5e.md. Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S15E_RSRC=<cartella con R/ e tools/>]
#   Rscript --vanilla R/testing/test_S1.5e.R      I mutanti sono copie modificate del codice (run_S1.5e.sh mutants).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow) })
SRC <- Sys.getenv("S15E_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5d_bench_lib.R"); source(file.path(SRC, "tools/S1.5e_bench_lib.R"))
options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")   # S1.6: predefiniti di S1.5c-S1.5e, per riprodurre i numeri di questo step
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "t_err"))
isok <- function(o) !inherits(o, "t_err")
PRE <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
square <- function(side) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = side^2),
  region_polygons = sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(side, 0), c(side, side), c(0, side), c(0, 0))))))
spec_c <- function(A) tissue_preset_spec(PRE[[A]], size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus")
comp_of <- function(sp) data.frame(cluster_id = 1L, cell_type = sp$composition$cell_type, fraction = sp$composition$fraction)
seed_h <- function(A, side, seed, het, quotas = TRUE) { sp <- spec_c(A)
  seed_centroids(square(side), sp$cell_types, comp_of(sp), random_seed = seed, spatial = "compartments", compartments = sp$compartments,
                 heterogeneity = het, quotas = quotas, verbose = FALSE) }
H <- function(pv = 0, lv = NA, sd = 0, ld = NA, mode = NA, k = NA, cap = NA) {
  h <- data.frame(void_frac = pv, void_ell_um = lv, dens_sigma = sd, dens_ell_um = ld)
  if (!is.na(mode)) h$void_mode <- mode; if (!is.na(k)) h$void_aspect <- k; if (!is.na(cap)) h$dens_cap <- cap; h }
elong_of <- function(v, px = 1) {        # allungamento mediano pesato per area delle componenti (>= 50 um^2) di una maschera logica
  lab <- imager::label(imager::as.cimg(v * 1), high_connectivity = TRUE); L <- matrix(as.numeric(lab), nrow(v), ncol(v)) + 1; L[!v] <- 0   # etichette da 0: +1 prima di togliere lo sfondo
  ids <- setdiff(unique(as.vector(L)), 0); out <- c(); wts <- c()
  for (i in ids) { k <- which(L == i, arr.ind = TRUE); if (nrow(k) * px^2 < 50) next
    ev <- eigen(stats::cov(k * px), symmetric = TRUE, only.values = TRUE)$values; out <- c(out, sqrt(max(ev) / max(min(ev), 1e-12))); wts <- c(wts, nrow(k)) }
  if (!length(out)) return(NA_real_); o <- order(out); cw <- cumsum(wts[o]) / sum(wts); out[o][which(cw >= 0.5)[1]] }

# ---- T1: con "remove" e campo anisotropo nessun centro dove V supera tau ----
s1 <- seed_h("A2", 400, 5, H(0.15, 20, mode = "remove", k = 4)); ri <- s1$info$compartments$regions[[1]]; tau <- qnorm(1 - 0.15)
Vc <- .sc_rff_eval(ri$field_void, s1$centroids$x, s1$centroids$y)
ck("T1 nessun centro con V > tau (togliere, kappa 4, tumore del colon)", length(Vc) > 0 && max(Vc) <= tau, sprintf("max V %.3f, tau %.3f, n %d", max(Vc), tau, length(Vc)))
# ---- T2: numero di cellule con "remove" = piano x (1 - v): per seme (scarto <= numero di tipi) e in media su 200 semi (3 errori Monte Carlo) ----
r2 <- t(vapply(1:200, function(sd) { s <- seed_h("A2", 200, 1000 + sd, H(0.15, 16, mode = "remove")); r <- s$info$compartments$regions[[1]]
  c(n = nrow(s$centroids), plan = r$n_plan, v = r$void_grid_frac, failed = s$info$n_failed) }, numeric(4)))
dev2 <- abs(r2[, "n"] - r2[, "plan"] * (1 - r2[, "v"])); rat2 <- r2[, "n"] / (r2[, "plan"] * (1 - r2[, "v"]))
ck("T2a per seme |n - piano (1 - v)| <= 2 (200 semi, 2 tipi)", all(dev2 <= 2) && all(r2[, "failed"] == 0), sprintf("max %.2f", max(dev2)))
ck("T2b in media n / (piano (1 - v)) = 1 entro 3 errori Monte Carlo", abs(mean(rat2) - 1) < 3 * sd(rat2) / sqrt(200) + 1e-12, sprintf("%.5f (3 ES %.5f)", mean(rat2), 3 * sd(rat2) / sqrt(200)))
# ---- T3: conteggi per tipo con "remove" = n_t (1 - v) arrotondato (scarto <= 1 per tipo), 6 preset ----
t3 <- vapply(names(PRE), function(A) { s <- seed_h(A, if (A %in% c("A5", "A6")) 600 else 300, 33, H(0.2, 16, mode = "remove")); r <- s$info$compartments$regions[[1]]
  sp <- spec_c(A); n0 <- seed_h(A, if (A %in% c("A5", "A6")) 600 else 300, 33, NULL)$centroids$cell_type
  e <- table(factor(n0, levels = sp$cell_types$cell_type)) * (1 - r$void_grid_frac); o <- table(factor(s$centroids$cell_type, levels = sp$cell_types$cell_type))
  max(abs(as.numeric(o) - as.numeric(e))) }, 0)
ck("T3 conteggi per tipo con togliere = n_t (1 - v) entro 1 (6 preset)", all(t3 <= 1), sprintf("%s %.2f", names(t3), t3))
# ---- T4: colonne assenti = valori neutri, bit per bit (3 preset); "move" esplicito uguale ----
same4 <- vapply(c("A1", "A3", "A6"), function(A) { a <- seed_h(A, 300, 9, H(0.15, 16, 0.3, 40)); b <- seed_h(A, 300, 9, H(0.15, 16, 0.3, 40, mode = "move", k = 1, cap = FALSE))
  identical(a$centroids, b$centroids) && identical(a$region_df, b$region_df) }, TRUE)
ck("T4 colonne nuove ai valori neutri identiche a S1.5d (bit per bit, 3 preset)", all(same4), names(same4)[!same4])
# ---- T5: covarianza del campo anisotropo lungo i due assi (200 campi, kappa 4, ell 10: lunghezze 20 e 5 um) ----
set.seed(55); co <- t(vapply(1:200, function(k) { f <- .sc_rff_draw_aniso(10, 4, 256L); th <- f$theta
  P <- cbind(runif(400, 0, 500), runif(400, 0, 500)); u1 <- c(cos(th), sin(th)); u2 <- c(-sin(th), cos(th))
  v0 <- .sc_rff_eval(f, P[, 1], P[, 2]); a <- .sc_rff_eval(f, P[, 1] + 20 * u1[1], P[, 2] + 20 * u1[2]); b <- .sc_rff_eval(f, P[, 1] + 5 * u2[1], P[, 2] + 5 * u2[2])
  c(mean(v0 * a), mean(v0 * b)) }, numeric(2)))
ck("T5 correlazione alla lunghezza dell'asse lungo e corto = exp(-1/2) entro 0.05 (200 campi)", all(abs(colMeans(co) - exp(-0.5)) < 0.05), sprintf("%.3f %.3f (atteso %.3f)", colMeans(co)[1], colMeans(co)[2], exp(-0.5)))
# ---- T6: orientamento uniforme in (0, pi), una estrazione per regione (Kolmogorov-Smirnov, 500 regioni) ----
th6 <- vapply(1:500, function(sd) seed_h("A3", 150, 6000 + sd, H(0.1, 16, mode = "remove", k = 4))$info$compartments$regions[[1]]$field_void$theta, 0)
p6 <- suppressWarnings(stats::ks.test(th6, "punif", 0, pi)$p.value)
ck("T6 orientamento del campo uniforme in (0, pi) (KS, p > 0.01)", p6 > 0.01 && length(unique(th6)) > 400, sprintf("p %.3f, distinti %d", p6, length(unique(th6))))
# ---- T7: allungamento dei vuoti del campo crescente con kappa (1, 2, 4, 8; 10 semi; 400 um, frazione 0.15, ell 16) ----
gx <- seq(0.5, 399.5, 1); G7 <- expand.grid(x = gx, y = gx)
e7 <- vapply(c(1, 2, 4, 8), function(k) mean(vapply(1:10, function(sd) { s <- seed_h("A2", 400, 700 + sd, H(0.15, 16, mode = "remove", k = k))
  fv <- s$info$compartments$regions[[1]]$field_void; v <- matrix(.sc_rff_eval(fv, G7$x, G7$y) > qnorm(0.85), 400, 400); elong_of(v) }, 0), na.rm = TRUE), 0)
ck("T7 allungamento mediano dei vuoti crescente con kappa (1, 2, 4, 8)", all(diff(e7) > 0), sprintf("%.2f", e7))
# ---- T8: tetto della densita': fattore <= H sulla griglia; caso saturo di S1.5d (linfonodo p2 seme 4, densita' 0.9 / 80 um) senza cellule perse ----
s8 <- seed_h("A4", 400, 8, H(0.05, 16, 0.9, 40, mode = "remove", cap = TRUE)); r8 <- s8$info$compartments$regions[[1]]
gg <- expand.grid(x = seq(1, 399, 2), y = seq(1, 399, 2)); hw <- .sc_het_weight(.sc_check_het(H(0.05, 16, 0.9, 40, mode = "remove", cap = TRUE)), r8$field_void, r8$field_dens, gg$x, gg$y, cap = r8$cap)
ck("T8a fattore della densita' con il tetto <= H in tutti i punti (linfonodo)", !is.null(r8$cap) && max(hw) <= r8$cap$H * (1 + 1e-12) && max(exp(0.9 * .sc_rff_eval(r8$field_dens, gg$x, gg$y))) > r8$cap$H,
   sprintf("max %.3f, H %.3f", max(hw), r8$cap$H))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]; rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
i8 <- which(rois$archetype == "A4" & rois$roi_id == "p2"); w8 <- r3_roi_window("A4", "p2", rois); sp8 <- spec_c("A4")
cp8 <- data.frame(cluster_id = w8$reg$region_df$cluster_id[1], cell_type = sp8$composition$cell_type, fraction = sp8$composition$fraction)
rd8 <- het_for("A4", "p2", "RD")
sat <- function(cap) { h <- rd8; if (cap) h$dens_cap <- TRUE
  suppressWarnings(seed_centroids(w8$reg, sp8$cell_types, cp8, random_seed = stage_seeds(20261021 + 1000 * i8 + 4)[["seed"]], spatial = "compartments",
                                  compartments = compart_for_d("A4", "p2"), heterogeneity = h, verbose = FALSE))$info$n_failed }
f8 <- c(senza = sat(FALSE), con = sat(TRUE))
ck("T8b linfonodo p2 seme 4 con la densita' di S1.5d: cellule perse senza tetto > 0, con tetto 0", f8[["senza"]] > 0 && f8[["con"]] == 0, sprintf("senza %d, con %d (parametri %s)", f8[["senza"]], f8[["con"]], paste(unlist(rd8), collapse = "/")))
# ---- T8c (dopo la revisione, RC-5): il tetto riscala il fattore a media 1 sui punti non vuoti; H ricalcolato da zero; con "remove" n_target = cellule + non piazzate ----
k8 <- r8$void_grid_frac; sp4 <- spec_c("A4"); f4 <- sp4$composition$fraction[match(sp4$cell_types$cell_type, sp4$composition$cell_type)]
dbar <- sum(f4 * sp4$cell_types$min_dist_um); ggk <- expand.grid(x = seq(1, 399, 2), y = seq(1, 399, 2))
inV <- .sc_rff_eval(r8$field_void, ggk$x, ggk$y) <= qnorm(1 - 0.05)
lam_bar <- nrow(s8$centroids) / (400^2 * (1 - r8$void_grid_frac)); H_ref <- 0.8 * 0.547 / (pi * dbar^2 / 4) / lam_bar
m8 <- mean(pmin(r8$cap$s * exp(0.9 * .sc_rff_eval(r8$field_dens, ggk$x[inV], ggk$y[inV])), r8$cap$H))
ck("T8c tetto: media del fattore fuori dai vuoti = 1 (0.02, griglia diversa), H = 0.8 * saturazione / densita' media fuori dai vuoti (1e-6)",
   abs(m8 - 1) < 0.02 && abs(r8$cap$H / H_ref - 1) < 1e-6 && s8$info$n_failed == 0, sprintf("media %.4f, H %.6f contro %.6f", m8, r8$cap$H, H_ref))
nt <- vapply(c("A1", "A2", "A5"), function(A) { s <- seed_h(A, if (A == "A5") 600 else 300, 81, H(0.2, 16, mode = "remove")); r <- s$region_df
  isTRUE(all(r$n_target == r$n_cells + r$n_failed)) && isTRUE(all(r$n_removed_void > 0)) }, TRUE)
ck("T8d con togliere: n_target della regione = cellule piazzate + non piazzate (3 preset)", all(nt), names(nt)[!nt])
# ---- T9: place_cells() sui 6 preset con le colonne nuove, senza ERROR; le colonne arrivano alla semina ----
e9 <- vapply(names(PRE), function(A) { cl <- expand.grid(x = 1:260, y = 1:260); cl$intensity_cluster <- 1L
  L <- safe(suppressWarnings(place_cells(cl, PRE[[A]], pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L, random_seed = 11,
            size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus", override_heterogeneity = H(0.15, 16, 0.3, 40, mode = "remove", k = 2, cap = TRUE), verbose = FALSE)))
  isok(L) && L$metadata$validation$n_error == 0 && identical(L$metadata$heterogeneity$void_mode, "remove") && L$metadata$heterogeneity$void_aspect == 2 }, TRUE)
ck("T9 place_cells() con togliere, kappa e tetto sui 6 preset: 0 ERROR, parametri registrati", all(e9), names(e9)[!e9])
cl <- expand.grid(x = 1:300, y = 1:300); cl$intensity_cluster <- 1L
pc <- function(h) suppressWarnings(place_cells(cl, "epithelial", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L, random_seed = 12,
  size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus", override_heterogeneity = h, validate = FALSE, verbose = FALSE))
n9 <- c(move = nrow(pc(H(0.2, 16))$cell_df), remove = nrow(pc(H(0.2, 16, mode = "remove"))$cell_df))
ck("T9b place_cells() inoltra void_mode: togliere riduce le cellule di circa il 20 %", abs(n9[["remove"]] / n9[["move"]] - 0.8) < 0.05, sprintf("%d contro %d", n9[["remove"]], n9[["move"]]))
# ---- T10: valori dei preset (tabella heterogeneity) uguali alle righe di BIO_REFERENCES ----
bio <- readLines("docs/BIO_REFERENCES.md"); bio <- bio[grepl("^\\| B-[0-9]{3} \\|", bio)]
bval <- function(id) { r <- bio[startsWith(bio, paste0("| ", id, " |"))]; if (!length(r)) return(NA_character_); trimws(strsplit(r, "\\|")[[1]][5]) }
fmt <- function(h) { s <- sprintf("%g / %g / %g / %g", h$void_frac, h$void_ell_um, h$dens_sigma, h$dens_ell_um)
  if ("void_mode" %in% names(h)) s <- sprintf("%s / %s / %g / %s", s, h$void_mode, h$void_aspect, h$dens_cap); s }
h10 <- vapply(seq_len(nrow(.TP_HET)), function(k) identical(bval(.TP_HET$src[k]), fmt(.TP_HET[k, ])), TRUE)
ck("T10 eterogeneita' dei 6 preset = righe di BIO_REFERENCES", all(h10), .TP_HET$src[!h10])
# ---- T11: il banco legge la piega lascia-fuori della ROI per X, XA, XAD ----
t11 <- vapply(c("X", "XA", "XAD", "XADs"), function(fam) { f <- sprintf("results/S1.5e/S1.5e_het_%s.csv", fam); if (!file.exists(f)) return(FALSE)
  d <- read.csv(f, stringsAsFactors = FALSE); ok <- TRUE
  for (A in unique(d$archetype)) for (roi in unique(rois$roi_id[rois$archetype == A])) { h <- het_for_e(A, roi, fam); e <- d[d$archetype == A & d$fold == paste0("loo_", roi), ]
    ok <- ok && nrow(e) == 1 && isTRUE(all.equal(c(h$void_frac, h$void_aspect, h$void_ell_um, h$dens_sigma), c(e$void_frac, e$void_aspect, e$void_ell_um, e$dens_sigma))) &&
      identical(is.na(h$dens_ell_um), is.na(e$dens_ell_um)) && (is.na(h$dens_ell_um) || h$dens_ell_um == e$dens_ell_um) && identical(h$dens_cap, as.logical(e$dens_cap)) &&
      identical(h$void_mode, e$void_mode) }   # dopo la revisione (RC-10): tutti i parametri
  ok }, TRUE)
ck("T11 het_for_e() legge la piega lascia-fuori della ROI (X, XA, XAD, XADs), tutti i parametri", all(t11), names(t11)[!t11])
cat(sprintf("RESULT %d/%d PASS\n", sum(unlist(res)), length(res)))
