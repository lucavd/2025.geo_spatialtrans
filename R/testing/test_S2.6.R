# R/testing/test_S2.6.R - S2.6: batteria degli strumenti della prova di prestazione (tools/S2.6_pieces.R, tools/S2.6_measures.R,
# tools/S2.6_section_layer.R). Cancello G-1 di docs/steps/S2.6.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.6.R
# Stampa PASS/FAIL per asserzione; esce con 1 se una fallisce. Variabile S26_TOOLS: cartella degli strumenti (mutanti).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
TOOLS <- Sys.getenv("S26_TOOLS", "tools")
source(file.path(TOOLS, "S2.6_pieces.R"))
nf <- 0L; np <- 0L
ok <- function(name, cond) { cond <- isTRUE(cond); cat(sprintf("%s  %s\n", if (cond) "PASS" else "FAIL", name)); if (cond) np <<- np + 1L else nf <<- nf + 1L }
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
mk_layer <- function(A, m, px, seed) {
  idx <- which(m, arr.ind = TRUE); cl <- data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L)
  suppressWarnings(place_cells(cl, PRESET[[A]], pixel_size_um = px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                    small_components = "exclude", random_seed = seed, types = "multi", size_arm = "recentred", spatial = "uniform",
                    size_resid = "gaussian", strict = FALSE, verbose = FALSE))
}
TMP <- tempfile("s26test"); dir.create(TMP)
run <- function(args) { invisible(NULL); st <- system2("Rscript", c("--vanilla", file.path(TOOLS, "S2.6_pieces.R"), args), stdout = TRUE, stderr = TRUE); invisible(attr(st, "status")) }
rd <- function(d) list(u2 = readBin(file.path(d, "u2.bin"), "integer", 1e8, size = 4L), c8 = readRDS(file.path(d, "c8.rds")),
                       own8 = readBin(file.path(d, "own8.bin"), "integer", 1e8, size = 4L), dup8 = readBin(file.path(d, "dup8.bin"), "integer", 1e8, size = 4L),
                       j = jsonlite::fromJSON(file.path(d, "run.json"), simplifyVector = FALSE))
# --- livelli di prova: quadrato del cuore (400 um), quadrato del tumore (240 um), due dischi del linfonodo (r = 80 um, maschera non piena)
LA6 <- mk_layer("A6", matrix(TRUE, 200, 200), 2, 11L); fA6 <- file.path(TMP, "A6.rds"); saveRDS(LA6, fA6)
LA2 <- mk_layer("A2", matrix(TRUE, 120, 120), 2, 12L); fA2 <- file.path(TMP, "A2.rds"); saveRDS(LA2, fA2)
xy <- expand.grid(i = 1:200, j = 1:200); md <- matrix(((xy$i - 50.5)^2 + (xy$j - 50.5)^2) <= 40^2 | ((xy$i - 150.5)^2 + (xy$j - 150.5)^2) <= 40^2, 200, 200)   # due dischi: in alcuni pezzi il tessuto non arriva ai lati del riquadro
LD <- mk_layer("A4", md, 2, 13L); fD <- file.path(TMP, "D.rds"); saveRDS(LD, fD)
s1 <- s26_seeds(20261201L, 1L)[1]
# T1: un pezzo unico senza bordo = simulate_visium_hd() bit per bit (cuore, tumore)
for (A in c("A6", "A2")) {
  f <- if (A == "A6") fA6 else fA2
  run(c("one", f, s1, TMP, paste0("one_", A), "0")); run(c("pieces", f, "100000", "0", "20261201", TMP, paste0("single_", A), "0"))
  o <- rd(file.path(TMP, paste0("one_", A))); p <- rd(file.path(TMP, paste0("single_", A)))
  ok(sprintf("T1 %s: un pezzo unico = un pezzo solo (UMI per quadratino identici, %d quadratini)", A, length(o$u2)), identical(o$u2, p$u2))
  ok(sprintf("T1 %s: matrice a 8 um identica", A), identical(as.matrix(o$c8 != 0), as.matrix(p$c8 != 0)) && max(abs(o$c8 - p$c8)) == 0)
  ok(sprintf("T1 %s: origine per blocco identica", A), identical(o$own8, p$own8 + p$dup8))
}
# T2...T5: pezzi con bordo sul cuore e sul disco del linfonodo
run(c("pieces", fA6, "200", "100", "20261201", TMP, "A6_T200H100", "0", "2"))
run(c("pieces", fD, "120", "56", "20261202", TMP, "D_T120H56", "1", "2"))
run(c("pieces", fD, "120", "0", "20261202", TMP, "D_T120H0", "0", "2"))
for (tg in c("A6_T200H100", "D_T120H56", "D_T120H0")) {
  p <- rd(file.path(TMP, tg)); g <- p$j$grid; nb <- g$n_row * g$n_col
  pz <- list.files(file.path(TMP, tg, "pieces"), full.names = TRUE); Z <- lapply(pz, readRDS); Z <- Z[!vapply(Z, function(z) isTRUE(z$empty), TRUE)]
  cov <- tabulate(unlist(lapply(Z, `[[`, "gbin")), nb)
  ok(sprintf("T2 %s: nessun quadratino coperto due volte; quadratini scoperti solo nei pezzi senza cellule (%d)", tg, nb),
     max(cov) == 1L && sum(cov) + sum(vapply(p$j$pieces, function(z) if (isTRUE(z$empty)) 1 else 0, 0)) >= 1 &&
     (all(cov == 1L) || any(vapply(p$j$pieces, function(z) isTRUE(z$empty), TRUE))))
  if (grepl("^D", tg)) ok(sprintf("T2 %s: ogni quadratino con il centro nell'area utile di un pezzo con cellule e' coperto", tg), {
    gg <- g; jj <- (seq_len(nb) - 1L) %% gg$n_col; ii <- (seq_len(nb) - 1L) %/% gg$n_col; xc <- gg$x0 + (jj + 0.5) * 2; yc <- gg$y0 + (ii + 0.5) * 2
    need <- rep(FALSE, nb)
    for (z in p$j$pieces) if (!isTRUE(z$empty)) { U <- unlist(z$U); need <- need | (xc >= U[1] & xc < U[3] & yc >= U[2] & yc < U[4]) }
    all(cov[need] == 1L) })
  ok(sprintf("T3 %s: UMI cuciti = somma dei pezzi", tg), sum(p$u2) == sum(vapply(Z, function(z) sum(z$X), 0)))
  n8 <- (g$n_row %/% 4L) * (g$n_col %/% 4L); id <- seq_len(nb)
  b8 <- as.integer((((id - 1L) %/% g$n_col) %/% 4L) * (g$n_col %/% 4L) + ((id - 1L) %% g$n_col) %/% 4L + 1L)
  ok(sprintf("T3 %s: matrice a 8 um = somma dei 16 quadratini", tg), all(abs(Matrix::colSums(p$c8) - tabulate(rep(b8, p$u2), n8)) == 0))
  ok(sprintf("T3 %s: origine propria + doppia <= UMI per blocco", tg), all(p$own8 + p$dup8 <= tabulate(rep(b8, p$u2), n8)))
  seeds <- vapply(p$j$pieces, function(z) as.numeric(z$seed), 0)
  ok(sprintf("T4 %s: semi dei pezzi tutti diversi (%d pezzi)", tg, length(seeds)), length(unique(seeds)) == length(seeds))
  H <- as.numeric(sub(".*H", "", tg))
  if (H > 0) ok(sprintf("T5 %s: trascritti da cellule del bordo presenti", tg), sum(p$dup8) > 0) else ok(sprintf("T5 %s: senza bordo nessuna cellula doppia", tg), sum(p$dup8) == 0)
}
# T6: ogni cellula ha un solo pezzo proprietario; cellule del pezzo = centro nel riquadro; territori tagliati dentro il riquadro
L <- LD; G <- s26_global_grid(L); gbb <- G$bbox; Tu <- 120; H <- 56
xs <- seq(gbb[["xmin"]], gbb[["xmax"]], by = Tu); if (tail(xs, 1) < gbb[["xmax"]]) xs <- c(xs, gbb[["xmax"]]); ys <- xs
own <- integer(nrow(L$cell_df)); allin <- TRUE; areas_ok <- TRUE; nclip <- 0L; n_short <- 0L
for (ix in seq_len(length(xs) - 1)) for (iy in seq_len(length(ys) - 1)) {
  U <- c(xs[ix], ys[iy], xs[ix + 1], ys[iy + 1]); box <- c(max(gbb[[1]], U[1] - H), max(gbb[[2]], U[2] - H), min(gbb[[3]], U[3] + H), min(gbb[[4]], U[4] + H))
  P <- s26_piece_layer(L, box, gbb); if (!nrow(P$cell_df)) next
  oi <- match(P$cell_df$cell_id, L$cell_df$cell_id)
  ow <- L$cell_df$x[oi] >= U[1] & (L$cell_df$x[oi] < U[3] | U[3] >= gbb[["xmax"]]) & L$cell_df$y[oi] >= U[2] & (L$cell_df$y[oi] < U[4] | U[4] >= gbb[["ymax"]])
  own[oi[ow]] <- own[oi[ow]] + 1L
  gB <- geos::as_geos_geometry(.s26_box_sfc(box))
  allin <- allin && all(geos::geos_within(geos::as_geos_geometry(P$cell_territories), geos::geos_buffer(gB, 1e-9)))
  a0 <- geos::geos_area(geos::as_geos_geometry(sf::st_set_crs(L$cell_territories[oi], NA))); a1 <- geos::geos_area(geos::as_geos_geometry(P$cell_territories))
  inside0 <- geos::geos_within(geos::as_geos_geometry(sf::st_set_crs(L$cell_territories[oi], NA)), gB)
  areas_ok <- areas_ok && all(abs(a1[inside0] - a0[inside0]) < 1e-9) && all(a1 <= a0 + 1e-9); nclip <- nclip + sum(!inside0)
  bb <- sf::st_bbox(c(sf::st_set_crs(P$cell_territories, NA), P$region_polygons))
  bt <- sf::st_bbox(sf::st_set_crs(P$cell_territories, NA)); if (any(abs(as.numeric(bt) - box) > 1)) n_short <- n_short + 1L
  if (!isTRUE(all.equal(as.numeric(bb), box, tolerance = 0))) allin <- FALSE
}
ok(sprintf("T6: ogni cellula ha un solo pezzo proprietario (%d cellule)", length(own)), all(own == 1L))
ok("T6: territori tagliati dentro il riquadro e riquadro del pezzo = riquadro di territori, regioni e ancore", allin)
ok(sprintf("T6: aree dei territori non tagliati invariate (%d tagliati)", nclip), areas_ok && nclip > 0)
ok(sprintf("T6: il caso di prova ha pezzi in cui il tessuto non arriva ai lati del riquadro (%d)", n_short), n_short > 0)
# T7: misure su un'esecuzione sintetica (quadrato di 2000 um, deficit del 10 % nella prima colonna di blocchi a destra di x = 1000)
SD <- file.path(TMP, "synth"); dir.create(SD)
g <- hd_grid_geometry(c(xmin = 0, ymin = 0, xmax = 2000, ymax = 2000), 280, 2); nb <- g$n_row * g$n_col
id <- seq_len(nb); jj <- (id - 1L) %% g$n_col; ii <- (id - 1L) %/% g$n_col; xc <- g$x0 + (jj + 0.5) * 2; yc <- g$y0 + (ii + 0.5) * 2
u2 <- rep(10L, nb); u2[xc > 1000 & xc < 1008 & abs(yc - 1000) > 16] <- 8L; u2[xc > 992 & xc < 1000 & abs(yc - 1000) > 16] <- 9L
n8 <- (g$n_row %/% 4L) * (g$n_col %/% 4L); b8 <- as.integer(((ii %/% 4L) * (g$n_col %/% 4L)) + jj %/% 4L + 1L)
c8 <- Matrix::Matrix(matrix(tabulate(rep(b8, u2), n8), 1), sparse = TRUE)
save_run(SD, u2, c8, integer(n8), integer(n8), list(tag = "synth", mode = "one", grid = g[c("x0", "y0", "n_row", "n_col")]))
invisible(system2("Rscript", c("--vanilla", file.path(TOOLS, "S2.6_measures.R"), SD, "1000"), stdout = FALSE, stderr = FALSE))
m <- read.csv(file.path(SD, "measures_T1000.csv")); pr <- read.csv(file.path(SD, "pairs_T1000.csv"))
ok("T7: fascia 0-8 um con il deficit atteso (996 blocchi: 246 a 144 a sinistra, 246 a 128 a destra, gli altri a 160)", abs(m$mean_umi[m$band == "0-8"] - (160 - 246 * 16 / 996 - 246 * 32 / 996)) < 1e-9)
ok("T7: fasce oltre 8 um senza deficit (160)", all(abs(m$mean_umi[m$band %in% c("8-16", "16-40", "40-100", "100-250")] - 160) < 1e-9))
ok("T7: coppie attraverso le giunture: due linee di 2000 um meno l'incrocio (circa 2 x 246)", pr$n_pairs >= 480 && pr$n_pairs <= 500)
unlink(TMP, recursive = TRUE)
cat(sprintf("\n%d PASS, %d FAIL\n", np, nf)); quit(status = if (nf == 0L) 0L else 1L)
