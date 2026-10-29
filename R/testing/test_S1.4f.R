# R/testing/test_S1.4f.R - S1.4f: test dell'orientamento parziale nel pacchetto e del banco (cancello G-1; mutanti G-2).
# Pre-registrazione: docs/steps/S1.4f.md. Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S14F_RSRC=<cartella con R/ e tools/>]
#   Rscript --vanilla R/testing/test_S1.4f.R      I mutanti sono copie modificate del codice (run_S1.4f.sh mutants).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow) })
SRC <- Sys.getenv("S14F_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5c_bench_lib.R"); source("tools/S1.5d_bench_lib.R"); source("tools/S1.4e_bench_lib.R")
source(file.path(SRC, "tools/S1.4f_bench_lib.R"))
options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")   # S1.6: predefiniti di S1.5c-S1.5e, per riprodurre i numeri di questo step
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
wkb <- function(n) lapply(as.list(geos::geos_write_wkb(geos::as_geos_geometry(n))), as.vector)
axis_of <- function(g) { xy <- sf::st_coordinates(g)[, 1:2]; m <- r3_poly_moments(xy[, 1], xy[, 2]); 0.5 * atan2(2 * m[["mxy"]], m[["mxx"]] - m[["myy"]]) }
axd <- function(a) (a + pi / 2) %% pi - pi / 2                      # differenza assiale in [-pi/2, pi/2)
rect <- function(cx, cy, w, h, al) { p <- rbind(c(-w, -h), c(w, -h), c(w, h), c(-w, h), c(-w, -h)) / 2
  R <- matrix(c(cos(al), sin(al), -sin(al), cos(al)), 2); q <- p %*% t(R); sf::st_polygon(list(cbind(q[, 1] + cx, q[, 2] + cy))) }
ctp <- function(area, asp, orient, conc = NULL, law = NULL, phi = 1) { x <- data.frame(cell_type = "c", density = 1e3, nucleus_to_eq_ratio = NA_real_,
  nucleus_radius = sqrt(area / pi), nucleus_area_sdlog = 0, nucleus_aspect = asp, nucleus_orient = orient, nucleus_nc_max = phi, stringsAsFactors = FALSE)
  if (!is.null(conc)) x$nucleus_orient_conc <- conc; if (!is.null(law)) x$nucleus_orient_law <- law; x }
tess_of <- function(polys) { g <- sf::st_sfc(polys); list(cell_territories = g, territory_df = data.frame(cell_id = seq_along(g), territory_area = as.numeric(sf::st_area(g)))) }
cen_of <- function(x, y) data.frame(cell_id = seq_along(x), region_id = 1L, x = x, y = y, cell_type = "c")
pwc <- function(phi, rho) 0.5 + atan((1 + rho) / (1 - rho) * tan(phi / 2)) / pi          # funzione di ripartizione della Cauchy avvolta su (-pi, pi)

# ---- T1: senza colonne nuove e con lambda = Inf i nuclei non cambiano (livello di place_cells, tumore del colon) ----
sq <- expand.grid(x = 1:200, y = 1:200); sq$intensity_cluster <- 1L
L <- suppressWarnings(place_cells(sq, "tumor_microenv", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L, small_components = "exclude",
  random_seed = 11, types = "multi", size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus", override_compartments = compart_for_d("A2", "r1"),
  override_heterogeneity = het_for("A2", "r1", "RV"), quotas = TRUE, size_resid = "gaussian", resid_mode = "decile", keep_tiles = TRUE, strict = FALSE, verbose = FALSE))
o_k2 <- redraw_nuclei(L, L$cell_types_resolved)
ck("T1 K2 ridisegnato = place_cells() con il pacchetto modificato", identical(wkb(o_k2$cell_nuclei), wkb(L$cell_nuclei)), sprintf("%d nuclei", nrow(L$cell_df)))
et0 <- arm_catalog(L$cell_types_resolved, "ET", "A2"); et1 <- et0; et1$nucleus_orient_conc <- Inf; et1$nucleus_orient_law <- "constant"
o_e0 <- redraw_nuclei(L, et0); o_e1 <- redraw_nuclei(L, et1)
ck("T1 ET con lambda = Inf esplicito = ET senza colonne nuove (bit per bit)", identical(wkb(o_e0$cell_nuclei), wkb(o_e1$cell_nuclei)) &&
   identical(o_e0$cell_df$nucleus_orient, o_e1$cell_df$nucleus_orient))
f14e <- "/mnt/micron/geo_spatialtrans/S1.4e/B/nuc/ET/A3_r1_s01.parquet"
if (file.exists(f14e)) ck("T1 file di S1.4e presente per il confronto di G-3", TRUE, f14e)
# ---- T2: lambda = 0: orientamento uniforme rispetto al territorio (2 000 rettangoli allungati lungo x) ----
gx <- rep(0:49, 40) * 40; gy <- rep(0:39, each = 50) * 25
t2 <- tess_of(lapply(seq_along(gx), function(k) rect(gx[k], gy[k], 36, 20, 0)))
o2 <- derive_cell_geometry(cen_of(gx, gy), t2, ctp(20, 1.6, "territory", conc = 0, law = "constant"), random_seed = 3, verbose = FALSE)
p2 <- suppressWarnings(stats::ks.test((axd(o2$cell_df$nucleus_orient) + pi / 2) / pi, "punif")$p.value)
ck("T2 lambda = 0: orientamenti uniformi (KS)", p2 > 0.01, sprintf("p = %.3f", p2))
# ---- T3: quantile della Cauchy avvolta e estrazioni ----
uu <- c(1e-6, 0.01, 0.2, 0.5, 0.77, 0.99, 1 - 1e-6); err <- 0
for (rho in c(0.1, 0.5, 0.9, 0.99)) err <- max(err, abs(pwc(.dcg_wc_quantile(uu, rho), rho) - uu))
ck("T3 F(Q(u)) = u entro 1e-12 (rho 0.1, 0.5, 0.9, 0.99)", err < 1e-12, sprintf("max %.1e", err))
set.seed(7); ud <- runif(20000); pk <- vapply(c(0.1, 0.5, 0.9, 0.99), function(rho) suppressWarnings(ks.test(.dcg_wc_quantile(ud, rho), function(q) pwc(q, rho))$p.value), 0)
ck("T3 20 000 estrazioni contro la funzione di ripartizione (KS, p > 0.01; equivale a F(Q(u)) = u, il controllo indipendente e' G-4)", all(pk > 0.01), sprintf("p %s", paste(round(pk, 3), collapse = " ")))
ck("T3 estremi: rho = 0 uniforme, rho = 1 zero", max(abs(.dcg_wc_quantile(uu, 0) - 2 * pi * (uu - 0.5))) < 1e-12 && all(.dcg_wc_quantile(uu, 1) == 0))
o3 <- derive_cell_geometry(cen_of(gx, gy), t2, ctp(20, 1.6, "territory", conc = -log(1 - 0.6), law = "constant"), random_seed = 4, verbose = FALSE)
th3 <- vapply(seq_along(gx), function(k) axis_of(o3$cell_nuclei[[k]]), 0)
p3 <- suppressWarnings(ks.test(2 * axd(th3), function(q) pwc(q, 0.6))$p.value)
ck("T3 nuclei di derive_cell_geometry(): 2 x (asse del poligono - asse del territorio) ~ Cauchy avvolta con rho = 0.6 (KS)", p3 > 0.01 && max(o3$cell_df$nucleus_scale) < 1 + 1e-6,
   sprintf("p = %.3f, scala massima %.8f", p3, max(o3$cell_df$nucleus_scale)))
# ---- T4: pesi delle leggi e angolo estratto, ricostruiti a mano ----
asp_r <- c(1.2, 1.5, 2, 3); ok_m <- TRUE; dev <- 0
for (a in asp_r) { mo <- .dcg_territory_moments(.dcg_as_polys(sf::st_sfc(rect(3, -2, 10 * a, 10, 0.7))))
  dev <- max(dev, abs(mo$aspect - a), abs(axd(mo$theta - 0.7))) }
ck("T4 allungamento del territorio = rapporto dei lati del rettangolo, asse = rotazione (entro 1e-9)", dev < 1e-9, sprintf("max %.1e", dev))
wq <- .dcg_orient_weight(c("constant", "elongation", "compression"), c(1.8, 1.8, 1.8), c(12, 12, 12), c(2, 2, 2))
ck("T4 pesi: costante 1, allungamento a - 1, compressione sqrt(A/pi) / d", max(abs(wq - c(1, 0.8, sqrt(12 / pi) / 2))) < 1e-12, paste(round(wq, 6), collapse = " "))
n4 <- 60; x4 <- (seq_len(n4) - 1) * 50; set.seed(9); wd <- runif(n4, 10, 30); ht <- runif(n4, 6, 12); al4 <- runif(n4, -pi / 2, pi / 2)
t4 <- tess_of(lapply(seq_len(n4), function(k) rect(x4[k], 0, wd[k], ht[k], al4[k]))); c4 <- cen_of(x4, rep(0, n4))
for (law in c("constant", "elongation", "compression")) {
  o4 <- derive_cell_geometry(c4, t4, ctp(25, 1.6, "territory", conc = 0.8, law = law), random_seed = 21, verbose = FALSE)
  set.seed(21); z <- rnorm(n4); u <- runif(n4)
  w <- switch(law, constant = rep(1, n4), elongation = wd / ht - 1, compression = sqrt(o4$cell_df$nucleus_area_target / pi) / (ht / 2))
  want <- al4 + 0.5 * .dcg_wc_quantile(u, 1 - exp(-0.8 * w)); e4 <- max(abs(axd(o4$cell_df$nucleus_orient - want)))
  ck(sprintf("T4 legge %s: angolo = asse del territorio + Q(u, rho_i) / 2 ricostruito a mano (entro 1e-9)", law), e4 < 1e-9, sprintf("max %.1e", e4))
}
# ---- T5: errori espliciti ----
er <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
cc5 <- cen_of(0, 0); t5 <- tess_of(list(rect(0, 0, 20, 10, 0)))
ck("T5 lambda negativo, legge sconosciuta, colonne nuove con orientamento a caso: errore",
   er(derive_cell_geometry(cc5, t5, ctp(10, 1.5, "territory", conc = -1), verbose = FALSE)) &&
   er(derive_cell_geometry(cc5, t5, ctp(10, 1.5, "territory", conc = 1, law = "nematic"), verbose = FALSE)) &&
   er(derive_cell_geometry(cc5, t5, ctp(10, 1.5, "random", conc = 1), verbose = FALSE)) &&
   er(derive_cell_geometry(cc5, t5, ctp(10, 1.5, "random", law = "compression"), verbose = FALSE)) &&
   !er(derive_cell_geometry(cc5, t5, ctp(10, 1.5, "random"), verbose = FALSE)))
# ---- T6: ellissi orientate in parte nel modello clip su territori non convessi e con buco ----
Lsh <- sf::st_polygon(list(rbind(c(0, 0), c(14, 0), c(14, 5), c(5, 5), c(5, 14), c(0, 14), c(0, 0))))
hol <- sf::st_polygon(list(rbind(c(30, 0), c(44, 0), c(44, 12), c(30, 12), c(30, 0)), rbind(c(36, 5), c(38, 5), c(38, 7), c(36, 7), c(36, 5))))
t6 <- tess_of(list(Lsh, hol, rect(70, 6, 24, 5, 0.3))); c6 <- cen_of(c(2.5, 33, 70), c(2.5, 6, 6))
for (law in c("constant", "elongation", "compression")) for (sd in 1:3) {
  o6 <- derive_cell_geometry(c6, t6, ctp(60, 1.7, "territory", conc = 1, law = law, phi = 0.8), random_seed = sd, verbose = FALSE)
  inv <- nuc_invariants(o6, t6$cell_territories)
  ck(sprintf("T6 legge %s seme %d: dentro, nessuna sovrapposizione, area = A*, valida", law, sd),
     inv$n_out == 0 && inv$n_overlap == 0 && inv$n_area_bad == 0 && inv$n_invalid == 0 && inv$n_empty == 0, sprintf("max errore relativo %.2e", inv$max_area_rel))
}
# ---- T8: BL-174, centri doppi al bordo ----
cc8 <- data.frame(x = c(0, 0.9, 1.5, 10, 10.2, 30), y = 0, src = c("manual", "CO", "CO", "manual", "CO", "CO"), luca_idx = c(1, NA, NA, 2, NA, NA))
d8 <- s14f_dedup(cc8)
ck("T8 tolti i centri del segmentatore entro 1 um da un punto di Luca, tenuti quelli a 1.5 um", identical(d8$x, c(0, 1.5, 10, 30)), paste(d8$x, collapse = " "))
# ---- T10: catalogo dei bracci orientati in parte ----
ctz <- data.frame(cell_type = c("a", "b"), nucleus_aspect = 1, nucleus_orient = "random", stringsAsFactors = FALSE); lc <- law_catalog(ctz, "A3", "V2", 0.75)
ck("T10 catalogo EP: allungamento vero, territorio, lambda e legge", all(abs(lc$nucleus_aspect - s14e_aspect("A3")) < 1e-12) && all(lc$nucleus_orient == "territory") &&
   all(lc$nucleus_orient_conc == 0.75) && all(lc$nucleus_orient_law == "compression"))
# ---- T11: la compressione usa l'area bersaglio prima dei tetti (aggiunto dopo la revisione, mutante X11) ----
t11 <- tess_of(list(rect(0, 0, 6, 4, 0.3), rect(20, 0, 7, 4, -0.5))); c11 <- cen_of(c(0, 20), c(0, 0))
o11 <- derive_cell_geometry(c11, t11, ctp(25, 1.5, "territory", conc = 1.5, law = "compression", phi = 0.6), random_seed = 5, verbose = FALSE)
set.seed(5); z <- rnorm(2); u <- runif(2); w11 <- sqrt(o11$cell_df$nucleus_area_target / pi) / o11$cell_df$center_dist
e11 <- max(abs(axd(o11$cell_df$nucleus_orient - (c(0.3, -0.5) + 0.5 * .dcg_wc_quantile(u, 1 - exp(-1.5 * w11))))))
ck("T11 compressione con l'area bersaglio prima del tetto (tetto attivo: A* < A bersaglio)", e11 < 1e-9 && all(o11$cell_df$nucleus_area_cap < o11$cell_df$nucleus_area_target),
   sprintf("max %.1e", e11))
# ---- T1b: K2 ridisegnato = place_cells() su altri due preset (aggiunto dopo la revisione) ----
for (pr in list(list("epithelial", list(spatial = "uniform", size_resid = "gaussian")), list("brain_cortex", list(spatial = "uniform", size_resid = "empirical")))) {
  Lb <- suppressWarnings(do.call(place_cells, c(list(sq, pr[[1]], pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L, small_components = "exclude",
    random_seed = 11, types = "multi", size_arm = "group", keep_tiles = TRUE, strict = FALSE, verbose = FALSE), pr[[2]])))
  ob <- redraw_nuclei(Lb, Lb$cell_types_resolved)
  ck(paste("T1b K2 ridisegnato = place_cells():", pr[[1]]), identical(wkb(ob$cell_nuclei), wkb(Lb$cell_nuclei)), sprintf("%d nuclei", nrow(Lb$cell_df)))
}
cat(sprintf("RESULT %d/%d PASS\n", sum(unlist(res)), length(res)))
