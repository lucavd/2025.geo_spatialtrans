# S1.4 — check C (C-1) di derive_cell_geometry() su casi analitici, con mutanti (C-mut)
# Pre-registrazione: results/S1.4/S1.4_preregistration.md (approvata da Luca 2026-10-15, commit c72c01a)
# Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S14_MUTANT=none|M1..M5] Rscript --vanilla R/testing/test_S1.4.R
# Stampa PASS/FAIL per asserzione; scrive results/S1.4/S1.4_test_<mutante>.csv.
# Definizioni operative (fissate prima di eseguire):
#  - territori costruiti a mano (territory_df con l'area GEOS esatta del poligono), piu' un caso d'integrazione con
#    tessellate_voronoi() (200 centroidi uniformi, seed 1, in un quadrato 100 x 100 µm);
#  - tolleranza delle aree 1e-9 relativa; posizioni del POI entro 2 * poi_precision_um (1e-3 µm) dal valore analitico;
#  - C3 verificato con sf::st_within (indipendente dal predicato geos interno); i casi con nuclei fuori dal territorio
#    sotto mutante si leggono con strict = FALSE.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
source("tools/S1.4_mutants.R")
MUT <- Sys.getenv("S14_MUTANT", "none")
E <- dcg_mutate(dcg_env(), MUT)
DCG <- function(...) E$derive_cell_geometry(..., nucleus_model = "scale", verbose = FALSE)   # S1.4d: il predefinito e' "clip"
RES <- "results/S1.4"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)
TOUT <- Sys.getenv("GEO_TEST_OUT", RES); dir.create(TOUT, showWarnings = FALSE, recursive = TRUE)   # S1.6 (BL-165, revisione RC-05): tabelle fuori dai file versionati se GEO_TEST_OUT e' dato
TOL <- 1e-9; PREC <- 1e-3
rows <- list()
rec <- function(check, input, metric, value, threshold, status) {
  rows[[length(rows) + 1]] <<- data.frame(check = check, input = input, metric = metric, value = as.character(value),
                                          threshold = threshold, status = status)
  cat(sprintf("%-5s %-6s %-24s %-30s %s (soglia %s)\n", status, check, input, metric, value, threshold))
}
pf <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
fmt <- function(v) formatC(v, format = "e", digits = 2)
rel <- function(a, b) max(abs(a - b) / abs(b))
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)), class = "dcg_err"))
P <- function(m) st_polygon(list(rbind(m, m[1, ])))
mk <- function(polys, gx = NULL, gy = NULL, type = "t") {
  g <- st_sfc(polys); n <- length(g)
  if (is.null(gx)) { pp <- st_coordinates(st_point_on_surface(g)); gx <- pp[, 1]; gy <- pp[, 2] }
  list(cen = data.frame(cell_id = seq_len(n), region_id = 1L, x = gx, y = gy, cell_type = type),
       tess = list(cell_territories = g, territory_df = data.frame(cell_id = seq_len(n), territory_area = as.numeric(st_area(g)))))
}
CT <- function(...) { a <- list(...); d <- data.frame(cell_type = "t", density = 1e4, nucleus_to_eq_ratio = 0.4)
  for (k in names(a)) d[[k]] <- a[[k]]; d }
nuc_area <- function(o) as.numeric(st_area(o$cell_nuclei))
within_sf <- function(o, tess) vapply(seq_along(o$cell_nuclei), function(i)
  lengths(st_within(o$cell_nuclei[i], tess$cell_territories[i])) == 1, TRUE)

sq   <- P(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10)))
hexa <- P(t(sapply(0:5, function(k) 5 * c(cos(k * pi / 3), sin(k * pi / 3)))))
rect <- P(rbind(c(0, 0), c(10, 0), c(10, 1), c(0, 1)))
Lsh  <- P(rbind(c(0, 0), c(10, 0), c(10, 2), c(2, 2), c(2, 10), c(0, 10)))
hole <- st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)), rbind(c(3, 3), c(7, 3), c(7, 7), c(3, 7), c(3, 3))))
tri  <- P(rbind(c(0, 0), c(20, 0), c(10, 1)))
rot  <- function(p, a) { m <- unclass(p)[[1]]; R <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2); st_polygon(list(m %*% t(R))) }
k32 <- sqrt(2 * pi / (32 * sin(2 * pi / 32)))

# T1 quadrato: POI, d, area = target, contenimento, citoplasma
x <- mk(list(sq)); o <- safe(DCG(x$cen, x$tess, CT()))
if (inherits(o, "dcg_err")) rec("C-1", "T1 quadrato", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  rec("C-1", "T1 quadrato", "POI", sprintf("%.4f,%.4f", cd$nucleus_cx, cd$nucleus_cy), "5,5 ± 2e-3", pf(abs(cd$nucleus_cx - 5) <= 2 * PREC && abs(cd$nucleus_cy - 5) <= 2 * PREC))
  rec("C-1", "T1 quadrato", "area target", fmt(rel(cd$nucleus_area_target, 16)), "16, rel 1e-9", pf(rel(cd$nucleus_area_target, 16) <= TOL))
  rec("C-1", "T1 quadrato", "area poligono = target", fmt(rel(nuc_area(o), 16)), "rel 1e-9", pf(rel(nuc_area(o), 16) <= TOL && !cd$nucleus_clamped))
  rec("C-1", "T1 quadrato", "citoplasma", fmt(rel(cd$cytoplasm_area, 84)), "84, rel 1e-9", pf(rel(cd$cytoplasm_area, 84) <= TOL))
  rec("C-1", "T1 quadrato", "C3 st_within", sum(within_sf(o, x$tess)), "1/1", pf(all(within_sf(o, x$tess))))
}
# T2 esagono: POI al centro, d = apotema
x <- mk(list(hexa)); o <- safe(DCG(x$cen, x$tess, CT()))
if (inherits(o, "dcg_err")) rec("C-1", "T2 esagono", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  rec("C-1", "T2 esagono", "POI e d", sprintf("%.4f,%.4f d=%.4f", cd$nucleus_cx, cd$nucleus_cy, cd$center_dist), "0,0 d=4.3301",
      pf(abs(cd$nucleus_cx) <= 2 * PREC && abs(cd$nucleus_cy) <= 2 * PREC && abs(cd$center_dist - 5 * sqrt(3) / 2) <= 2 * PREC))
  rec("C-1", "T2 esagono", "area poligono = target", fmt(rel(nuc_area(o), cd$nucleus_area_target)), "rel 1e-9", pf(rel(nuc_area(o), cd$nucleus_area_target) <= TOL))
}
# T3 rettangolo 1:10, cerchio clampato: area = pi (0.85 d / k)^2
x <- mk(list(rect)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.5)))
if (inherits(o, "dcg_err")) rec("C-1", "T3 rettangolo", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df; aexp <- pi * (0.85 * cd$center_dist / k32)^2
  rec("C-1", "T3 rettangolo", "d", sprintf("%.5f", cd$center_dist), "0.5 ± 2e-3", pf(abs(cd$center_dist - 0.5) <= 2 * PREC))
  rec("C-1", "T3 rettangolo", "clamp cerchio", sprintf("clamped=%s area rel %s", cd$nucleus_clamped, fmt(rel(nuc_area(o), aexp))), "TRUE, rel 1e-9",
      pf(cd$nucleus_clamped && rel(nuc_area(o), aexp) <= TOL && rel(cd$nucleus_area, aexp) <= TOL))
  rec("C-1", "T3 rettangolo", "C3 st_within", sum(within_sf(o, x$tess)), "1/1", pf(all(within_sf(o, x$tess))))
}
# T4 rettangolo, ellisse 10:1 lungo il territorio: nessun clamp, area 2.5, orientazione 0
o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.5, nucleus_aspect = 10, nucleus_orient = "territory")))
if (inherits(o, "dcg_err")) rec("C-1", "T4 ellisse lungo", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  rec("C-1", "T4 ellisse lungo", "non clampata, area 2.5", sprintf("clamped=%s rel %s", cd$nucleus_clamped, fmt(rel(nuc_area(o), 2.5))), "FALSE, rel 1e-9",
      pf(!cd$nucleus_clamped && rel(nuc_area(o), 2.5) <= TOL))
  rec("C-1", "T4 ellisse lungo", "orientazione", fmt(abs(cd$nucleus_orient)), "0 ± 1e-9", pf(abs(cd$nucleus_orient) <= TOL))
}
# T5 stessa ellisse ortogonale (angolo fisso pi/2): clampata e contenuta
o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.5, nucleus_aspect = 10, nucleus_orient = as.character(pi / 2)), strict = FALSE))
if (inherits(o, "dcg_err")) rec("C-1", "T5 ellisse ortog.", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  rec("C-1", "T5 ellisse ortog.", "clampata e contenuta", sprintf("clamped=%s s=%.4f within=%s", cd$nucleus_clamped, cd$nucleus_scale, all(within_sf(o, x$tess))),
      "TRUE, within", pf(cd$nucleus_clamped && all(within_sf(o, x$tess))))
  rec("C-1", "T5 ellisse ortog.", "area = s^2 target", fmt(rel(nuc_area(o), cd$nucleus_scale^2 * 2.5)), "rel 1e-9", pf(rel(nuc_area(o), cd$nucleus_scale^2 * 2.5) <= TOL))
}
# T6 L: POI analitico c = 2 sqrt2 / (1 + sqrt2)
x <- mk(list(Lsh)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.2), strict = FALSE)); cL <- 2 * sqrt(2) / (1 + sqrt(2))
if (inherits(o, "dcg_err")) rec("C-1", "T6 L", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  rec("C-1", "T6 L", "POI", sprintf("%.4f,%.4f d=%.4f", cd$nucleus_cx, cd$nucleus_cy, cd$center_dist), sprintf("%.4f ± 2e-3", cL),
      pf(abs(cd$nucleus_cx - cL) <= 2 * PREC && abs(cd$nucleus_cy - cL) <= 2 * PREC && abs(cd$center_dist - cL) <= 2 * PREC))
  rec("C-1", "T6 L", "C3 st_within", sum(within_sf(o, x$tess)), "1/1", pf(all(within_sf(o, x$tess))))
}
# T7 quadrato con buco: POI analitico c = (10 + 7 sqrt2) / (1 + sqrt2)
x <- mk(list(hole)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.2), strict = FALSE)); cH <- (10 + 7 * sqrt(2)) / (1 + sqrt(2))
if (inherits(o, "dcg_err")) rec("C-1", "T7 buco", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  okp <- min(sqrt((c(cd$nucleus_cx) - c(cH, 10 - cH, cH, 10 - cH))^2 + (cd$nucleus_cy - c(cH, cH, 10 - cH, 10 - cH))^2)) <= 3 * PREC
  rec("C-1", "T7 buco", "POI (uno dei 4 angoli)", sprintf("%.4f,%.4f d=%.4f", cd$nucleus_cx, cd$nucleus_cy, cd$center_dist), sprintf("%.4f d=%.4f", cH, 10 - cH),
      pf(okp && abs(cd$center_dist - (10 - cH)) <= 2 * PREC))
  rec("C-1", "T7 buco", "C3 st_within", sum(within_sf(o, x$tess)), "1/1", pf(all(within_sf(o, x$tess))))
}
# T8 triangolo sottile: POI = incentro, d = raggio inscritto
x <- mk(list(tri)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.3), strict = FALSE))
s2 <- sqrt(10^2 + 1); rin <- 2 * 10 / (20 + 2 * s2); ix <- (20 * 10 + s2 * 0 + s2 * 20) / (20 + 2 * s2)
if (inherits(o, "dcg_err")) rec("C-1", "T8 triangolo", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  rec("C-1", "T8 triangolo", "POI = incentro", sprintf("%.4f,%.4f d=%.5f", cd$nucleus_cx, cd$nucleus_cy, cd$center_dist), sprintf("%.4f,%.4f d=%.5f", ix, rin, rin),
      pf(abs(cd$center_dist - rin) <= 2 * PREC && abs(cd$nucleus_cy - rin) <= 2 * PREC))
  rec("C-1", "T8 triangolo", "clamp e C3", sprintf("clamped=%s within=%s", cd$nucleus_clamped, all(within_sf(o, x$tess))), "TRUE, within",
      pf(cd$nucleus_clamped && all(within_sf(o, x$tess))))
}
# T9 centro = generatore
x <- mk(list(sq), gx = 2, gy = 3); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.1), nucleus_center = "generator"))
if (inherits(o, "dcg_err")) rec("C-1", "T9 generatore", "esecuzione", o$msg, "nessun errore", "FAIL") else
  rec("C-1", "T9 generatore", "centro e d", sprintf("%g,%g d=%g", o$cell_df$nucleus_cx, o$cell_df$nucleus_cy, o$cell_df$center_dist), "2,3 d=2",
      pf(o$cell_df$nucleus_cx == 2 && o$cell_df$nucleus_cy == 3 && abs(o$cell_df$center_dist - 2) <= TOL))
# T10 raggio assoluto: area = pi r^2
x <- mk(list(sq)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = NA, nucleus_radius = 1, nucleus_area_slope = 0.3)))
if (inherits(o, "dcg_err")) rec("C-1", "T10 assoluto", "esecuzione", o$msg, "nessun errore", "FAIL") else
  rec("C-1", "T10 assoluto", "area = pi", fmt(rel(nuc_area(o), pi)), "rel 1e-9", pf(rel(nuc_area(o), pi) <= TOL && rel(o$cell_df$nucleus_area, pi) <= TOL))
# T11 beta = 0: area indipendente dal territorio; T12 beta = 0.5
sq20 <- P(rbind(c(20, 0), c(40, 0), c(40, 20), c(20, 20)))
x <- mk(list(sq, sq20)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.3, nucleus_area_slope = 0)))
if (inherits(o, "dcg_err")) rec("C-1", "T11 beta0", "esecuzione", o$msg, "nessun errore", "FAIL") else
  rec("C-1", "T11 beta0", "A_N = 9 in entrambi", paste(signif(o$cell_df$nucleus_area, 10), collapse = ","), "9,9 rel 1e-9", pf(rel(o$cell_df$nucleus_area, c(9, 9)) <= TOL && rel(nuc_area(o), c(9, 9)) <= TOL))
o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.3, nucleus_area_slope = 0.5)))
if (inherits(o, "dcg_err")) rec("C-1", "T12 beta0.5", "esecuzione", o$msg, "nessun errore", "FAIL") else
  rec("C-1", "T12 beta0.5", "A_N = 9, 18", paste(signif(o$cell_df$nucleus_area, 10), collapse = ","), "9,18 rel 1e-9", pf(rel(o$cell_df$nucleus_area, c(9, 18)) <= TOL))
# T13 orientazione "territory" su un rettangolo ruotato di 30°
x <- mk(list(rot(P(rbind(c(0, 0), c(20, 0), c(20, 4), c(0, 4))), pi / 6)))
o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.3, nucleus_aspect = 2, nucleus_orient = "territory")))
if (inherits(o, "dcg_err")) rec("C-1", "T13 ruotato", "esecuzione", o$msg, "nessun errore", "FAIL") else
  rec("C-1", "T13 ruotato", "orientazione = pi/6", fmt(abs(o$cell_df$nucleus_orient - pi / 6)), "1e-9", pf(abs(o$cell_df$nucleus_orient - pi / 6) <= TOL))
# T14 sigma: 2000 quadrati identici, sd(log A_N) e media; seme e ordine
many <- lapply(0:1999, function(i) { dx <- 11 * (i %% 50); dy <- 11 * (i %/% 50); P(rbind(c(dx, dy), c(dx + 10, dy), c(dx + 10, dy + 10), c(dx, dy + 10))) })
x <- mk(many); ctS <- CT(nucleus_to_eq_ratio = 0.3, nucleus_area_sdlog = 0.3)
o1 <- safe(DCG(x$cen, x$tess, ctS, random_seed = 7)); o2 <- safe(DCG(x$cen, x$tess, ctS, random_seed = 7)); o3 <- safe(DCG(x$cen, x$tess, ctS, random_seed = 8))
if (inherits(o1, "dcg_err")) rec("C-1", "T14 sigma", "esecuzione", o1$msg, "nessun errore", "FAIL") else {
  la <- log(o1$cell_df$nucleus_area_target)
  rec("C-1", "T14 sigma", "sd(log A_N)", sprintf("%.4f", sd(la)), "0.3 ± 0.03", pf(abs(sd(la) - 0.3) <= 0.03))
  rec("C-1", "T14 sigma", "media log A_N", sprintf("%.4f", mean(la)), sprintf("%.4f ± 0.03", log(9)), pf(abs(mean(la) - log(9)) <= 0.03))
  rec("C-rip", "T14 sigma", "stesso seme identico", identical(serialize(o1[1:2], NULL), serialize(o2[1:2], NULL)), "TRUE",
      pf(identical(serialize(o1[1:2], NULL), serialize(o2[1:2], NULL))))
  rec("C-rip", "T14 sigma", "seme diverso: POI =, aree !=", sprintf("POI %s, aree %s", identical(o1$cell_df$nucleus_cx, o3$cell_df$nucleus_cx),
      identical(o1$cell_df$nucleus_area, o3$cell_df$nucleus_area)), "TRUE, FALSE",
      pf(identical(o1$cell_df$nucleus_cx, o3$cell_df$nucleus_cx) && !identical(o1$cell_df$nucleus_area, o3$cell_df$nucleus_area)))
  pr <- sample(nrow(x$cen)); o4 <- DCG(x$cen[pr, ], x$tess, ctS, random_seed = 7)
  rec("C-rip", "T14 sigma", "righe permutate", identical(o4$cell_df, o1$cell_df), "TRUE", pf(identical(o4$cell_df, o1$cell_df)))
}
# T15 stima: dati senza rumore -> parametri esatti; con rumore -> entro 3 SE
at <- exp(seq(log(20), log(400), length.out = 500)); an <- 0.3^2 * 100 * (at / 100)^0.6
fp <- E$fit_nucleus_size_model(an, at, 100)
rec("C-1", "T15 stima", "esatta (rho, beta, sigma)", sprintf("%.12f %.12f %.1e", fp[["rho"]], fp[["beta"]], fp[["sigma"]]), "0.3 0.6 0 (1e-9)",
    pf(abs(fp[["rho"]] - 0.3) <= TOL && abs(fp[["beta"]] - 0.6) <= TOL && fp[["sigma"]] <= 1e-9))
fd <- E$fit_nucleus_size_model(an, at, 100, "design"); f0 <- E$fit_nucleus_size_model(an, at, 100, "beta0")
rec("C-1", "T15 stima", "vincoli dei modelli", sprintf("design b=%g s=%g; beta0 b=%g", fd[["beta"]], fd[["sigma"]], f0[["beta"]]), "1 0; 0",
    pf(fd[["beta"]] == 1 && fd[["sigma"]] == 0 && f0[["beta"]] == 0))
# T16 integrazione con tessellate_voronoi(): 200 centroidi uniformi in 100 x 100
set.seed(1); n16 <- 200
reg <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = 1e4),
            region_polygons = st_sfc(P(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100)))))
cen <- data.frame(cell_id = seq_len(n16), region_id = 1L, x = runif(n16, 0, 100), y = runif(n16, 0, 100), cell_type = "t")
tv <- tessellate_voronoi(cen, reg, verbose = FALSE)
for (asp in c(1, 2)) {
  o <- safe(DCG(cen, tv, CT(density = 2e4, nucleus_to_eq_ratio = 0.6, nucleus_area_sdlog = 0.2, nucleus_aspect = asp), random_seed = 3, strict = FALSE))
  lab <- sprintf("T16 tess asp=%g", asp)
  if (inherits(o, "dcg_err")) { rec("C-1", lab, "esecuzione", o$msg, "nessun errore", "FAIL"); next }
  cd <- o$cell_df; w <- within_sf(o, tv); un <- !cd$nucleus_clamped
  rec("C3", lab, "st_within", sprintf("%d/%d", sum(w), n16), "200/200", pf(all(w)))
  rec("C5", lab, "area > 0", sprintf("%d/%d", sum(cd$nucleus_area > 0), n16), "200/200", pf(all(cd$nucleus_area > 0)))
  rec("C-area", lab, "area poligono = s^2 target", fmt(rel(nuc_area(o), cd$nucleus_area)), "rel 1e-9", pf(rel(nuc_area(o), cd$nucleus_area) <= TOL))
  rec("C-area", lab, "non clampate: = target", fmt(rel(cd$nucleus_area[un], cd$nucleus_area_target[un])), "rel 1e-9",
      pf(any(un) && rel(cd$nucleus_area[un], cd$nucleus_area_target[un]) <= TOL))
  rec("C-area", lab, "citoplasma >= 0, identita'", sprintf("min %.3f; rel %s", min(cd$cytoplasm_area), fmt(rel(cd$cytoplasm_area + cd$nucleus_area, cd$territory_area))),
      ">= 0; 1e-9", pf(min(cd$cytoplasm_area) >= 0 && rel(cd$cytoplasm_area + cd$nucleus_area, cd$territory_area) <= TOL))
  rec("C-1", lab, "frazione clampata", sprintf("%.3f", mean(cd$nucleus_clamped)), "riportata", "PASS")
}
# T18 sfc misto (POLYGON + MULTIPOLYGON con orfano): POI nella parte principale, C3, orientazione (D-1)
mpo <- st_multipolygon(list(list(rbind(c(30, 0), c(31, 0), c(31, 1), c(30, 1), c(30, 0))), list(rbind(c(12, 0), c(28, 0), c(28, 4), c(12, 4), c(12, 0)))))
x <- mk(list(sq, mpo), gx = c(5, 20), gy = c(5, 2))
stopifnot(inherits(x$tess$cell_territories, "sfc_GEOMETRY"))
o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.3, nucleus_aspect = 1.5, nucleus_orient = "territory")))
if (inherits(o, "dcg_err")) rec("C-1", "T18 sfc misto", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df
  rec("C-1", "T18 sfc misto", "POI parte principale", sprintf("%.3f,%.3f d=%.3f", cd$nucleus_cx[2], cd$nucleus_cy[2], cd$center_dist[2]), "x in [14,26], y=2, d=2",
      pf(cd$nucleus_cx[2] >= 14 && cd$nucleus_cx[2] <= 26 && abs(cd$nucleus_cy[2] - 2) <= 2 * PREC && abs(cd$center_dist[2] - 2) <= 2 * PREC))
  rec("C-1", "T18 sfc misto", "C3 st_within", sum(within_sf(o, x$tess)), "2/2", pf(all(within_sf(o, x$tess))))
  # atteso analitico: momenti combinati di rettangolo 16x4 (centro 20,2) e quadrato 1x1 (centro 30.5,0.5)
  Ai <- c(64, 1); cxi <- c(20, 30.5); cyi <- c(2, 0.5); ox <- c(16^2, 1) / 12; oy <- c(4^2, 1) / 12
  cxa <- sum(Ai * cxi) / 65; cya <- sum(Ai * cyi) / 65
  mxx <- sum(Ai * (ox + (cxi - cxa)^2)) / 65; myy <- sum(Ai * (oy + (cyi - cya)^2)) / 65; mxy <- sum(Ai * (cxi - cxa) * (cyi - cya)) / 65
  th_exp <- 0.5 * atan2(2 * mxy, mxx - myy)
  rec("C-1", "T18 sfc misto", "orientazione multiparte", sprintf("%.6f vs %.6f", cd$nucleus_orient[2], th_exp), "analitico ± 1e-9", pf(abs(cd$nucleus_orient[2] - th_exp) <= TOL))
}
# T19 clamp dell'ellisse: s massimo analitico (rettangolo 1x10, aspetto 20 lungo x) — X1 bisezione, X2 margine 1/0.85
x <- mk(list(rect)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.9, nucleus_aspect = 20, nucleus_orient = "0"), strict = FALSE))
if (inherits(o, "dcg_err")) rec("C-1", "T19 ellisse max", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  cd <- o$cell_df; A19 <- 0.81 * 10; a19 <- sqrt(A19 * 20 / pi); b19 <- sqrt(A19 / (20 * pi))
  s_exp <- min(1, 0.85 * min(cd$nucleus_cx, 10 - cd$nucleus_cx) / (k32 * a19), 0.85 * min(cd$nucleus_cy, 1 - cd$nucleus_cy) / (k32 * b19))
  rec("C-1", "T19 ellisse max", "s = massimo analitico", sprintf("%.10f vs %.10f", cd$nucleus_scale, s_exp), "rel 1e-9", pf(abs(cd$nucleus_scale / s_exp - 1) <= TOL))
  rec("C3", "T19 ellisse max", "C3 st_within", sum(within_sf(o, x$tess)), "1/1", pf(all(within_sf(o, x$tess))))
}
# T20 flusso casuale dichiarato: set.seed; z <- rnorm(n); u <- runif(n) — X4 ordine, X8 intervallo dell'angolo
x <- mk(many[1:50]); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.3, nucleus_area_sdlog = 0.3, nucleus_aspect = 1.2), random_seed = 11))
if (inherits(o, "dcg_err")) rec("C-1", "T20 flusso", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  set.seed(11); z20 <- rnorm(50); u20 <- runif(50)
  rec("C-rip", "T20 flusso", "aree target = rnorm del seme", fmt(rel(o$cell_df$nucleus_area_target, 9 * exp(0.3 * z20))), "rel 1e-9",
      pf(rel(o$cell_df$nucleus_area_target, 9 * exp(0.3 * z20)) <= TOL))
  rec("C-rip", "T20 flusso", "angoli = pi u - pi/2", fmt(max(abs(o$cell_df$nucleus_orient - (pi * u20 - pi / 2)))), "1e-12",
      pf(max(abs(o$cell_df$nucleus_orient - (pi * u20 - pi / 2))) <= 1e-12))
}
# T21 asse del territorio con un buco decentrato: momenti analitici ad area con segno — X5
hole2 <- st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)), rbind(c(1, 1), c(5, 1), c(5, 3), c(1, 3), c(1, 1))))
x <- mk(list(hole2)); o <- safe(DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0.2, nucleus_aspect = 1.5, nucleus_orient = "territory")))
if (inherits(o, "dcg_err")) rec("C-1", "T21 asse con buco", "esecuzione", o$msg, "nessun errore", "FAIL") else {
  Ai <- c(100, -8); cxi <- c(5, 3); cyi <- c(5, 2); ox <- c(100, 16) / 12; oy <- c(100, 4) / 12; At <- sum(Ai)
  cxa <- sum(Ai * cxi) / At; cya <- sum(Ai * cyi) / At
  mxx <- sum(Ai * (ox + (cxi - cxa)^2)) / At; myy <- sum(Ai * (oy + (cyi - cya)^2)) / At; mxy <- sum(Ai * (cxi - cxa) * (cyi - cya)) / At
  th21 <- 0.5 * atan2(2 * mxy, mxx - myy)
  rec("C-1", "T21 asse con buco", "orientazione analitica", sprintf("%.6f vs %.6f", o$cell_df$nucleus_orient, th21), "1e-9", pf(abs(o$cell_df$nucleus_orient - th21) <= TOL))
}
# T22 centro = generatore sul bordo -> errore (CODE-05: area 0 invisibile a C3)
x <- mk(list(sq), gx = 0, gy = 5)
m22 <- tryCatch({ DCG(x$cen, x$tess, CT(), nucleus_center = "generator"); "nessun errore" }, error = function(e) conditionMessage(e))
rec("C-1", "T22 generatore sul bordo", "messaggio", substr(m22, 1, 60), "sul bordo", pf(grepl("sul bordo", m22, fixed = TRUE)))
# T17 errori con messaggio atteso
x <- mk(list(sq))
err_cases <- list(
  list("tipo mancante", function() DCG(transform(x$cen, cell_type = "zz"), x$tess, CT()), "tipi senza catalogo"),
  list("rho = 0", function() DCG(x$cen, x$tess, CT(nucleus_to_eq_ratio = 0)), "servono nucleus_radius"),
  list("aspect < 1", function() DCG(x$cen, x$tess, CT(nucleus_aspect = 0.5)), "nucleus_aspect >= 1"),
  list("orient non valido", function() DCG(x$cen, x$tess, CT(nucleus_orient = "su")), "nucleus_orient non valido"),
  list("clamp_factor = 1", function() DCG(x$cen, x$tess, CT(), clamp_factor = 1), "clamp_factor"),
  list("n_vertices = 6", function() DCG(x$cen, x$tess, CT(), n_vertices = 6), "n_vertices"),
  list("cell_id non corrispondenti", function() DCG(transform(x$cen, cell_id = 5L), x$tess, CT()), "non corrispondono"),
  list("shape incoerente", function() DCG(x$cen, x$tess, CT(nucleus_shape = "circle", nucleus_aspect = 2)), "incoerente"),
  list("sigma < 0", function() DCG(x$cen, x$tess, CT(nucleus_area_sdlog = -1)), "nucleus_area_sdlog"),
  list("territorio linea", function() { y <- x; y$tess$cell_territories <- st_sfc(st_linestring(rbind(c(0, 0), c(1, 1)))); DCG(y$cen, y$tess, CT()) }, "territori non poligonali"))
for (ec in err_cases) {
  m <- tryCatch({ ec[[2]](); "nessun errore" }, error = function(e) conditionMessage(e))
  rec("C-1", ec[[1]], "messaggio", substr(m, 1, 60), ec[[3]], pf(grepl(ec[[3]], m, fixed = TRUE)))
}

tab <- do.call(rbind, rows); tab$mutant <- MUT
f <- file.path(TOUT, sprintf("S1.4_test_%s.csv", MUT)); write.csv(tab, f, row.names = FALSE)
cat(sprintf("\nTOTALE %s: %d PASS, %d FAIL su %d -> %s\n", MUT, sum(tab$status == "PASS"), sum(tab$status == "FAIL"), nrow(tab), f))
