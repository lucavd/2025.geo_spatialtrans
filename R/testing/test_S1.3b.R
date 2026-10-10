# R/testing/test_S1.3b.R — S1.3b (BL-069): C-unit, casi analitici per le varianti di precisione.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S1.3b.R
# Stampa PASS/FAIL per asserzione e scrive results/S1.3b/S1.3b_test.csv. Le righe V0 sono INFO (nessuna attesa).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
source("tools/S1.3b_overlay.R")
R <- list(); chk <- function(case, check, variant, ok, value, info = FALSE) {
  st <- if (info) "INFO" else if (isTRUE(ok)) "PASS" else "FAIL"
  R[[length(R) + 1L]] <<- data.frame(case, check, variant, status = st, value = as.character(value))
  cat(sprintf("%-4s %-22s %-14s %-4s %s\n", st, case, check, variant, value)) }
ulp <- function(x, k) { y <- x; for (s in seq_len(k)) y <- nextafter_up(y); y }
nextafter_up <- function(x) { e <- 2^(floor(log2(abs(x))) - 52); x + e }
sq <- function(x0, y0, x1, y1) st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
VARS <- c("V0", "V0g", "V1", "V2", "V3")

run_case <- function(name, terr, reg, A_true) {
  P <- s3b_perim(terr); Preg <- s3b_perim(reg)
  for (v in VARS) {
    V <- s3b_variant(terr, reg, v); m <- s3b_measure(V, P, Preg, pts = NULL, full = v != "V0g")
    info <- v %in% c("V0", "V0g"); tau <- G_CHECK[[v]] * (Preg + sum(P))
    chk(name, "C-1 unione", v, abs(m$ov) <= tau && m$symd <= tau, sprintf("ov %.3e symd %.3e tau %.1e", m$ov, m$symd, tau), info)
    pmax_ <- if (m$n_pairs) max(m$pairs$area) else 0
    chk(name, "C-2 coppie", v, m$n_pairs_over == 0, sprintf("coppie %d, max area %.3e", m$n_pairs, pmax_), info)
    chk(name, "area totale", v, abs(m$sum_a - A_true) <= max(tau, 1e-12 * A_true), sprintf("%.15g", m$sum_a), info)
    if (v != "V0g") chk(name, "C-4 valide", v, m$n_invalid == 0 && m$n_empty == 0, sprintf("invalide %d vuote %d", m$n_invalid, m$n_empty), info)
  }
}
# U1: due quadrati con lato comune; il vertice (1001, 1) del secondo spostato di k ulp in x
for (k in c(1L, 4L, 16L)) {
  a <- sq(1000, 0, 1001, 1)
  b <- st_polygon(list(rbind(c(1001, 0), c(1002, 0), c(1002, 1), c(ulp(1001, k), 1), c(1001, 0))))
  run_case(sprintf("U1 %d ulp", k), st_sfc(a, b), st_sfc(sq(1000, 0, 1002, 1)), 2)
}
# U2: giunzione a T — A senza vertice in (1002, 1); B1, B2 lo condividono con coordinate diverse di qualche ulp
a <- sq(1000, 0, 1002, 2)
b1 <- st_polygon(list(rbind(c(1002, 0), c(1004, 0), c(1004, 1), c(ulp(1002, 3), 1), c(1002, 0))))
b2 <- st_polygon(list(rbind(c(ulp(1002, 7), 1), c(1004, 1), c(1004, 2), c(1002, 2), c(ulp(1002, 7), 1))))
run_case("U2 giunzione T", st_sfc(a, b1, b2), st_sfc(sq(1000, 0, 1004, 2)), 8)
nd <- s3b_node(list(sq(1000, 0, 1004, 2), a, b1, b2), 1e-7)
chk("U2 giunzione T", "V3 inserimenti", "V3", nd$info[["n_tjunction"]] >= 1, sprintf("T %d, spostati %d", nd$info[["n_tjunction"]], nd$info[["n_moved"]]))
# U3: V1 mette le coordinate sulla griglia
set.seed(1); pr <- st_sfc(st_polygon(list(rbind(c(0, 0), c(10 + runif(1), 0), c(10, 10 + runif(1)), c(0, 0)) + 1000)))
xy <- st_coordinates(s3b_round(pr, 1e-6))[, 1:2]
chk("U3 griglia", "coordinate su 1e-6", "V1", all(abs(xy * 1e6 - round(xy * 1e6)) < 1e-3), sprintf("max scarto %.2e griglie", max(abs(xy * 1e6 - round(xy * 1e6)))))
xy2 <- st_coordinates(s3b_sf(list(kind = "geos"), geos_set_precision(as_geos_geometry(pr), 1e-6)))[, 1:2]
chk("U3 griglia", "coordinate su 1e-6", "V2", all(abs(xy2 * 1e6 - round(xy2 * 1e6)) < 1e-3), sprintf("max scarto %.2e griglie", max(abs(xy2 * 1e6 - round(xy2 * 1e6)))))
# U4: rappresentante = prima geometria
g1 <- sq(0, 0, 1, 1); g2 <- st_polygon(list(rbind(c(1 + 5e-8, 0), c(2, 0), c(2, 1), c(1, 1 + 3e-8), c(1 + 5e-8, 0))))
nd <- s3b_node(list(g1, g2), 1e-7); c2 <- unclass(nd$geoms[[2]])[[1]]
chk("U4 rappresentante", "vertici di g2 = g1", "V3", any(c2[, 1] == 1 & c2[, 2] == 0) && any(c2[, 1] == 1 & c2[, 2] == 1) && identical(unclass(nd$geoms[[1]]), unclass(g1)),
    sprintf("spostati %d", nd$info[["n_moved"]]))
nd0 <- s3b_node(list(g1, g2), 1e-9)
chk("U4 rappresentante", "tol 1e-9: nessuna fusione", "V3", nd0$info[["n_moved"]] == 0, sprintf("spostati %d", nd0$info[["n_moved"]]))
# U5: perimetri
holed <- st_polygon(list(rbind(c(0, 0), c(4, 0), c(4, 4), c(0, 4), c(0, 0)), rbind(c(1, 1), c(1, 2), c(2, 2), c(2, 1), c(1, 1))))
pp <- s3b_perim(st_sfc(sq(0, 0, 1, 1), holed))
chk("U5 perimetro", "quadrato 4, con buco 20", "-", isTRUE(all.equal(pp, c(4, 20), tolerance = 1e-14)), paste(pp, collapse = ", "))
# U6: conversione sf <-> geos senza perdita
set.seed(2); q <- st_sfc(st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0)) * pi + 1000 + runif(1))))
back <- st_as_sfc(as_wkb(as_geos_geometry(q)))
chk("U6 conversione", "coordinate identiche", "-", identical(st_coordinates(back), st_coordinates(q)), "bit per bit")
chk("U6 conversione", "area GEOS 3.11 = sf", "-", abs(geos_area(as_geos_geometry(q)) - as.numeric(st_area(q))) <= 1e-12 * as.numeric(st_area(q)),
    sprintf("%.3e", geos_area(as_geos_geometry(q)) - as.numeric(st_area(q))))
# U7: perturbazione: chiusura degli anelli mantenuta, spostamento <= delta
pt <- s3b_perturb(st_sfc(holed + 1000), 1e-8, 3); m0 <- st_coordinates(st_sfc(holed + 1000)); m1 <- st_coordinates(pt)
chk("U7 perturbazione", "|spostamento| <= delta, valido", "-", max(abs(m1[, 1:2] - m0[, 1:2])) <= 1e-8 && all(st_is_valid(pt)), sprintf("max %.2e", max(abs(m1[, 1:2] - m0[, 1:2]))))
dir.create("results/S1.3b", showWarnings = FALSE)
out <- do.call(rbind, R); write.csv(out, "results/S1.3b/S1.3b_test.csv", row.names = FALSE)
cat(sprintf("\nPASS %d  FAIL %d  INFO %d\n", sum(out$status == "PASS"), sum(out$status == "FAIL"), sum(out$status == "INFO")))
