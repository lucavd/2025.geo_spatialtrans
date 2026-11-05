# R/testing/test_S2.1.R - S2.1: batteria di test di build_hd_grid() e della verita' geometrica (cancello G-1).
# Pre-registrazione: docs/steps/S2.1.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S21_RSRC=<cartella con R/>] Rscript --vanilla R/testing/test_S2.1.R
# I mutanti sono copie modificate del codice (tools/S2.1_mutants.py) passate con S21_RSRC.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(Matrix) })
SRC <- Sys.getenv("S21_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "t_err"))
g <- 1e-6
poly <- function(m) sf::st_polygon(list(rbind(m, m[1, ])))
sq_poly <- function(x0, y0, s) poly(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s)))
# --- oracolo indipendente: ritaglio di Sutherland-Hodgman con un rettangolo, area con la formula di Gauss ---
clip_rect <- function(P, xmin, xmax, ymin, ymax) {
  clip <- function(P, inside, inter) { if (!nrow(P)) return(P); out <- NULL; n <- nrow(P)
    for (k in seq_len(n)) { a <- P[k, ]; b <- P[if (k == n) 1 else k + 1, ]
      ia <- inside(a); ib <- inside(b)
      if (ia && ib) out <- rbind(out, b) else if (ia && !ib) out <- rbind(out, inter(a, b)) else if (!ia && ib) out <- rbind(out, inter(a, b), b) }
    if (is.null(out)) matrix(numeric(0), 0, 2) else out }
  ix <- function(xc) function(a, b) { t <- (xc - a[1]) / (b[1] - a[1]); c(xc, a[2] + t * (b[2] - a[2])) }
  iy <- function(yc) function(a, b) { t <- (yc - a[2]) / (b[2] - a[2]); c(a[1] + t * (b[1] - a[1]), yc) }
  P <- clip(P, function(p) p[1] >= xmin, ix(xmin)); P <- clip(P, function(p) p[1] <= xmax, ix(xmax))
  P <- clip(P, function(p) p[2] >= ymin, iy(ymin)); P <- clip(P, function(p) p[2] <= ymax, iy(ymax)); P }
shoelace <- function(P) if (nrow(P) < 3) 0 else abs(sum(P[, 1] * c(P[-1, 2], P[1, 2]) - c(P[-1, 1], P[1, 1]) * P[, 2])) / 2
oracle <- function(grid, P) {
  out <- NULL
  for (i in 0:(grid$n_row - 1)) for (j in 0:(grid$n_col - 1)) {
    x0 <- grid$x0 + 2 * j; y0 <- grid$y0 + 2 * i; a <- shoelace(clip_rect(P, x0, x0 + 2, y0, y0 + 2))
    if (a > 0) out <- rbind(out, data.frame(bin_id = i * grid$n_col + j + 1, a = a)) }
  out }
G16 <- list(x0 = 0, y0 = 0, bin_um = 2, n_row = 16L, n_col = 16L, block = 4L, margin_um = 0)
# --- livello cellulare sintetico: Voronoi di punti casuali in un quadrato, nuclei = dischi ritagliati nel territorio ---
synth_layer <- function(side = 120, n = 150, seed = 1, shift = c(0, 0), region_from = 0) {
  set.seed(seed); pts <- cbind(runif(n, region_from, side), runif(n, 0, side))
  R <- sf::st_sfc(sq_poly(region_from, 0, side - region_from))
  P <- sf::st_sfc(lapply(seq_len(n), function(k) sf::st_point(pts[k, ])))
  vor <- sf::st_collection_extract(sf::st_voronoi(sf::st_multipoint(pts), envelope = R[[1]]), "POLYGON")
  vor <- sf::st_sfc(vor)[unlist(sf::st_intersects(P, sf::st_sfc(vor)))]
  Tt <- overlay_prec("intersection", vor, R, area = FALSE)
  Nn <- overlay_prec("intersection", sf::st_buffer(P, 2.2, nQuadSegs = 8), Tt, area = FALSE)
  Tt <- Tt + shift; Nn <- Nn + shift; R <- R + shift
  list(cell_df = data.frame(cell_id = seq_len(n), x = pts[, 1] + shift[1], y = pts[, 2] + shift[2],
                            territory_area = as.numeric(sf::st_area(Tt)), nucleus_area = as.numeric(sf::st_area(Nn))),
       region_polygons = R, cell_territories = Tt, cell_nuclei = Nn, metadata = list(pixel_size_um = 1.1, random_seed = seed))
}
tau_pieces <- function(n, P) g * P + n * g^2
# ================= T-A casi analitici =================
oA <- hd_overlap(G16, sf::st_sfc(sq_poly(2, 2, 4)), sf::st_sfc(sq_poly(2, 2, 4)), tile_um = Inf)
ck("T-A1 quadrato di 4 um allineato: 4 quadratini, aree 4 (entro tau), nucleo uguale",
   nrow(oA) == 4 && all(abs(oA$area_territory_um2 - 4) <= tau_pieces(1, 8) + 1e-12) && all(abs(oA$area_nucleus_um2 - 4) <= 1e-5) &&
     setequal(oA$bin_id, c(1 * 16 + 1 + 1, 1 * 16 + 2 + 1, 2 * 16 + 1 + 1, 2 * 16 + 2 + 1)))
oB <- hd_overlap(G16, sf::st_sfc(sq_poly(2.5, 2.5, 4)), NULL, tile_um = Inf)
w <- c(1.5, 2, 0.5); expB <- data.frame(bin_id = c(outer(1:3, 1:3, function(i, j) i * 16 + j + 1)), a = c(outer(w, w)))
mB <- merge(oB, expB, by = "bin_id")
ck("T-A2 quadrato spostato di 0.5 um: 9 quadratini con aree prodotto di 1.5, 2, 0.5",
   nrow(oB) == 9 && nrow(mB) == 9 && max(abs(mB$area_territory_um2 - mB$a)) <= tau_pieces(1, 8), sprintf("max scarto %.2e", max(abs(mB$area_territory_um2 - mB$a))))
Ptri <- rbind(c(1.3, 0.7), c(27.1, 5.9), c(9.4, 22.6))
th <- seq(0, 2 * pi, length.out = 1001)[-1001]; Pgon <- cbind(16 + 11.3 * cos(th), 15 + 11.3 * sin(th))
for (nm in c("triangolo", "poligono di 1000 lati")) {
  P <- if (nm == "triangolo") Ptri else Pgon
  o <- hd_overlap(G16, sf::st_sfc(poly(P)), NULL, tile_um = Inf); ora <- oracle(G16, P); m <- merge(o, ora, by = "bin_id", all = TRUE)
  m[is.na(m)] <- 0; per <- as.numeric(sf::st_length(sf::st_cast(sf::st_sfc(poly(P)), "LINESTRING")))
  tol <- g * (8 * nrow(m) + per) + nrow(m) * g^2
  ck(sprintf("T-A3 %s contro il ritaglio indipendente (Sutherland-Hodgman): stesse coppie, scarti entro tau", nm),
     all(abs(m$area_territory_um2 - m$a) <= g * (8 + per) + g^2) && abs(sum(m$area_territory_um2) - shoelace(P)) <= tol && sum(m$a > 1e-9 & m$area_territory_um2 == 0) == 0,
     sprintf("coppie %d, max scarto %.2e", nrow(m), max(abs(m$area_territory_um2 - m$a))))
}
# ================= T-B convenzione della griglia =================
L <- synth_layer(); H <- build_hd_grid(L)
bbL <- sf::st_bbox(c(L$cell_territories, L$region_polygons))
ck("T-B1 origine = angolo del riquadro meno 280 um; righe e colonne multiple di 4; la griglia copre riquadro + margine",
   abs(H$grid$x0 - (bbL[["xmin"]] - 280)) < 1e-12 && abs(H$grid$y0 - (bbL[["ymin"]] - 280)) < 1e-12 && H$grid$n_row %% 4 == 0 && H$grid$n_col %% 4 == 0 &&
     H$grid$x0 + 2 * H$grid$n_col >= bbL[["xmax"]] + 280 - 1e-9 && H$grid$y0 + 2 * H$grid$n_row >= bbL[["ymax"]] + 280 - 1e-9,
   sprintf("%d x %d", H$grid$n_row, H$grid$n_col))
b <- H$bins
ck("T-B2 bin_id = i * n_col + j + 1; centro (x0 + 2 j + 1, y0 + 2 i + 1); i = riga (y), j = colonna (x)",
   identical(as.integer(b$bin_id), seq_len(nrow(b))) && all(b$bin_id == b$i * H$grid$n_col + b$j + 1) &&
     max(abs(b$x - (H$grid$x0 + 2 * b$j + 1))) < 1e-9 && max(abs(b$y - (H$grid$y0 + 2 * b$i + 1))) < 1e-9 &&
     all(diff(b$x[b$i == 0]) > 0) && all(diff(b$y[b$j == 0]) > 0) && nrow(b) == H$grid$n_row * H$grid$n_col)
p <- c(L$cell_df$x[5], L$cell_df$y[5]); ip <- floor((p[2] - H$grid$y0) / 2); jp <- floor((p[1] - H$grid$x0) / 2)
idp <- ip * H$grid$n_col + jp + 1
ck("T-B3 il quadratino che contiene il centro di una cellula sta nella verita' geometrica di quella cellula",
   any(H$truth$overlap$bin_id == idp & H$truth$overlap$cell_id == 5) && abs(b$x[idp] - p[1]) <= 1 && abs(b$y[idp] - p[2]) <= 1)
ck("T-B6 indici interi (bin_id e cell_id della verita' geometrica, bin_id e bin8_id della tabella): un indice double come 2e5 diventa \"2e+05\" come testo",
   is.integer(H$truth$overlap$bin_id) && is.integer(H$truth$overlap$cell_id) && is.integer(H$bins$bin_id) && is.integer(H$bins$bin8_id) &&
     identical(as.character(200000L), "200000"))
L0 <- L; L0$metadata$pixel_size_um <- NULL
ck("T-B4 errore senza pixel_size_um (BL-005); errore con margine 279 o 4 um; la scala passa nei metadati",
   inherits(safe(build_hd_grid(L0)), "t_err") && inherits(safe(build_hd_grid(L, margin_um = 279)), "t_err") &&
     inherits(safe(build_hd_grid(L, margin_um = 4)), "t_err") && identical(H$metadata$pixel_size_um, 1.1))
# livello esatto: regione [0, 40]^2 in 4 territori quadrati, piu' una regione senza cellule [41, 51] x [0, 10]: quadratini [40, 42] a meta'
Tq <- sf::st_sfc(sq_poly(0, 0, 20), sq_poly(20, 0, 20), sq_poly(0, 20, 20), sq_poly(20, 20, 20))
Lh <- list(cell_df = data.frame(cell_id = 1:4, x = c(10, 30, 10, 30), y = c(10, 10, 30, 30), territory_area = 400, nucleus_area = 400),
           region_polygons = sf::st_sfc(sq_poly(0, 0, 40), sq_poly(41, 0, 10)), cell_territories = Tq, cell_nuclei = Tq,
           metadata = list(pixel_size_um = 0.5))
Hh <- build_hd_grid(Lh); half <- Hh$bins[Hh$bins$x == 41 & Hh$bins$y > 0 & Hh$bins$y < 10, ]
ck("T-B5 in_tissue = frazione di tessuto >= 0.5 (quadratini a meta' esatta: in tessuto); tessuto totale con due regioni",
   Hh$grid$x0 == -280 && Hh$grid$x0 + 2 * Hh$grid$n_col >= 51 + 280 && Hh$grid$n_col %% 4 == 0 && nrow(half) == 5 && all(abs(half$tissue_frac - 0.5) < 1e-6) && all(half$in_tissue) &&
     abs(sum(Hh$bins$tissue_area_um2) - 1700) < 1e-6 && sum(Hh$bins$in_tissue) == 400 + 20 + 10,
   sprintf("%d quadratini a meta', in tessuto %d", nrow(half), sum(Hh$bins$in_tissue)))
# ================= T-C somma a 8 um =================
set.seed(3); nb <- H$grid$n_row * H$grid$n_col
X <- Matrix::rsparsematrix(5, nb, density = 0.01, rand.x = function(n) as.numeric(rpois(n, 3) + 1))
S <- sum_to_8um(X, H$grid); b8 <- (b$i %/% 4) * (H$grid$n_col %/% 4) + b$j %/% 4 + 1
Sb <- vapply(seq_len(ncol(S)), function(q) Matrix::rowSums(X[, b8 == q, drop = FALSE]), numeric(5))
ck("T-C1 somma a 8 um di una matrice intera = somma per forza bruta dei 16 figli (identita' esatta); blocco = (i %/% 4, j %/% 4)",
   identical(unname(as.matrix(S)), unname(Sb)) && all(table(b8) == 16) && all(b8 == b$bin8_id) && sum(S) == sum(X))
o8 <- sum_to_8um(H$truth$overlap[, c("bin_id", "cell_id", "area_territory_um2", "area_nucleus_um2")], H$grid, by = "cell_id")
ck("T-C2 somma a 8 um della verita' geometrica: totali per cellula invariati",
   max(abs(tapply(o8$area_territory_um2, o8$cell_id, sum) - tapply(H$truth$overlap$area_territory_um2, H$truth$overlap$cell_id, sum))) < 1e-9)
# ================= T-D invarianze =================
ovk <- function(H) H$truth$overlap[, c("bin_id", "cell_id", "area_territory_um2", "area_nucleus_um2")]
H100 <- build_hd_grid(L, tile_um = 100); Hinf <- build_hd_grid(L, tile_um = Inf)
same <- function(a, b, tol) nrow(a) == nrow(b) && all(a$bin_id == b$bin_id) && all(a$cell_id == b$cell_id) &&
  max(abs(a$area_territory_um2 - b$area_territory_um2)) <= tol && max(abs(a$area_nucleus_um2 - b$area_nucleus_um2)) <= tol
ck("T-D1 tessere di 100 um, 500 um e una sola tessera: stessa verita' geometrica e stessa area di tessuto (entro 16 g)",
   same(ovk(H100), ovk(H), 16 * g) && same(ovk(Hinf), ovk(H), 16 * g) && max(abs(H100$bins$tissue_area_um2 - Hinf$bins$tissue_area_um2)) <= 16 * g)
L1 <- synth_layer(shift = c(1, 1)); H1 <- build_hd_grid(L1); L4 <- synth_layer(shift = c(1e4, 1e4)); H4 <- build_hd_grid(L4)
ck("T-D2 cellule e griglia traslate insieme di 1 um (mezzo quadratino) e di 1e4 um (BL-104): stesse coppie e stesse aree",
   same(ovk(H1), ovk(H), 16 * g) && same(ovk(H4), ovk(H), 16 * g) && max(abs(H4$bins$tissue_area_um2 - H$bins$tissue_area_um2)) <= 16 * g)
ovs <- hd_overlap(H$grid, L1$cell_territories, L1$cell_nuclei); nks <- tabulate(ovs$cell_id, nrow(L$cell_df))
perS <- as.numeric(sf::st_length(sf::st_cast(L1$cell_territories, "MULTILINESTRING")))
ck("T-D3 solo le cellule spostate di mezzo quadratino, griglia ferma: totali per cellula invariati (entro tau), coppie diverse (la prova ha potere)",
   all(abs(tapply(ovs$area_territory_um2, factor(ovs$cell_id, levels = L$cell_df$cell_id), sum) - L1$cell_df$territory_area) <= g * (8 * nks + perS) + nks * g^2) &&
     nrow(ovs) != nrow(ovk(H)), sprintf("coppie %d contro %d", nrow(ovs), nrow(ovk(H))))
H2 <- build_hd_grid(L)
ck("T-D4 due esecuzioni bit per bit uguali", identical(H2[c("grid", "bins", "bins_8um", "truth")], H[c("grid", "bins", "bins_8um", "truth")]))
# ================= T-E invarianti (come G-2 sul pannello) =================
ov <- H$truth$overlap; perT <- as.numeric(sf::st_length(sf::st_cast(L$cell_territories, "MULTILINESTRING"))); nk <- tabulate(ov$cell_id, nrow(L$cell_df))
lev <- factor(ov$cell_id, levels = L$cell_df$cell_id)
sT <- tapply(ov$area_territory_um2, lev, sum); sN <- tapply(ov$area_nucleus_um2, lev, sum)
ck("T-E1 ogni cellula ha quadratini (BL-010) e la somma sui quadratini = area del territorio e del nucleo (entro g (8 n + P) + n g^2)",
   all(nk > 0) && all(abs(sT - L$cell_df$territory_area) <= g * (8 * nk + perT) + nk * g^2) && all(abs(sN - L$cell_df$nucleus_area) <= g * (8 * nk + perT) + nk * g^2))
sb <- tapply(ov$area_territory_um2, ov$bin_id, sum); ib <- as.integer(names(sb)); nbk <- tabulate(ov$bin_id, nrow(b))[ib]
ck("T-E2 per quadratino: somma dei territori = area di tessuto, <= 4 um^2; nucleo <= territorio per coppia; tessuto totale = area della regione",
   all(abs(sb - b$tissue_area_um2[ib]) <= 16 * g * (nbk + 1)) && all(sb <= 4 + 16 * g * nbk) && all(ov$area_nucleus_um2 <= ov$area_territory_um2 + 16 * g) &&
     abs(sum(b$tissue_area_um2) - 120^2) <= 16 * g * sum(b$tissue_area_um2 > 0 & b$tissue_area_um2 < 4) + 1e-9 &&
     sum(b$tissue_area_um2[-ib]) == 0)
G8 <- modifyList(H$grid, list(bin_um = 8, n_row = H$grid$n_row %/% 4L, n_col = H$grid$n_col %/% 4L, block = 1L))
d8 <- hd_overlap(G8, L$cell_territories, L$cell_nuclei)
o8b <- o8[o8$area_territory_um2 > 0, ]; m8 <- merge(o8b, d8, by.x = c("bin8_id", "cell_id"), by.y = c("bin_id", "cell_id"))
ck("T-E3 somma a 8 um della verita' geometrica = verita' calcolata direttamente sui quadrati di 8 um",
   nrow(m8) == nrow(o8b) && nrow(m8) == nrow(d8) && max(abs(m8$area_territory_um2.x - m8$area_territory_um2.y)) <= 16 * 16 * g &&
     max(abs(m8$area_nucleus_um2.x - m8$area_nucleus_um2.y)) <= 16 * 16 * g, sprintf("coppie %d / %d", nrow(m8), nrow(d8)))
n_fail <- sum(!unlist(res)); cat(sprintf("\nS2.1 test: %d PASS, %d FAIL\n", sum(unlist(res)), n_fail))
quit(status = if (n_fail) 1L else 0L)
