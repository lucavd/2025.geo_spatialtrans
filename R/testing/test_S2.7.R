# R/testing/test_S2.7.R - S2.7: batteria della simulazione a pezzi (R/05b6_simulate_visium_hd_pieces.R, ripiego della triangolazione in
# R/05b3, flussi per blocco in R/05b2). Cancello G-1 di docs/steps/S2.7.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.7.R [gruppi]
# Stampa PASS/FAIL per asserzione; esce con 1 se una fallisce. Variabile S27_RDIR: cartella dei sorgenti (mutanti); argomento facoltativo:
# gruppi da eseguire (per esempio "a,b"); predefinito tutti.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
RDIR <- Sys.getenv("S27_RDIR", "R")
for (f in sort(list.files(RDIR, pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
GR <- commandArgs(TRUE); GR <- if (length(GR)) strsplit(GR[1], ",")[[1]] else c("a", "b", "c", "d", "e", "g", "h")
nf <- 0L; np <- 0L
ok <- function(name, cond) { cond <- isTRUE(cond); cat(sprintf("%s  %s\n", if (cond) "PASS" else "FAIL", name)); if (cond) np <<- np + 1L else nf <<- nf + 1L }
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
mk_layer <- function(A, m, px, seed) {
  idx <- which(m, arr.ind = TRUE); cl <- data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L)
  suppressWarnings(place_cells(cl, PRESET[[A]], pixel_size_um = px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                    small_components = "exclude", random_seed = seed, types = "multi", size_arm = "recentred", spatial = "uniform",
                    size_resid = "gaussian", strict = FALSE, verbose = FALSE))
}
same_out <- function(a, b) {
  nm <- c("overlap", "overlap_8um", "origin", "origin_8um", "background", "background_8um", "arrived_2um", "cell_counts", "cells", "genes")
  identical(a$counts_2um, b$counts_2um) && identical(a$counts_8um, b$counts_8um) && identical(a$bins, b$bins) && identical(a$bins_8um, b$bins_8um) &&
    identical(a$truth[nm], b$truth[nm]) && identical(a$truth$layers, b$truth$layers)
}
SIDE <- c(A1 = 150L, A2 = 150L, A3 = 150L, A4 = 120L, A5 = 200L, A6 = 200L)        # pixel di 2 um
TS <- Sys.getenv("S27_TISSUES", ""); if (nzchar(TS)) SIDE <- SIDE[strsplit(TS, ",")[[1]]]
# --- (a) uguaglianza bit per bit fra divisioni e numero di processi; (f) controlli; (g) struttura
if ("a" %in% GR) {
  for (A in names(SIDE)) {
    Lr <- mk_layer(A, matrix(TRUE, SIDE[[A]], SIDE[[A]]), 2, 100L + match(A, names(SIDE)))
    o1 <- simulate_visium_hd_pieces(Lr, seed = 7L, piece_um = Inf, frame_open = TRUE)
    o2 <- simulate_visium_hd_pieces(Lr, seed = 7L, piece_um = 56, frame_open = TRUE, workers = 3L)
    o3 <- simulate_visium_hd_pieces(Lr, seed = 7L, piece_um = 168, frame_open = TRUE, workers = 1L)
    ok(sprintf("a %s: pezzo solo = pezzi di 56 um con 3 processi (%d cellule, %d pezzi, %.0f trascritti scambiati)", A, nrow(Lr$cell_df),
               o2$metadata$n_pieces, o2$metadata$n_sent_between_pieces), same_out(o1, o2))
    ok(sprintf("a %s: pezzo solo = pezzi di 168 um con 1 processo", A), same_out(o1, o3))
    ok(sprintf("a %s: trascritti scambiati fra pezzi > 0 con 56 um", A), o2$metadata$n_sent_between_pieces > 0)
    ok(sprintf("f %s: controlli di conservazione e copertura (3 esecuzioni)", A), all(unlist(o1$metadata$checks)) && all(unlist(o2$metadata$checks)) && all(unlist(o3$metadata$checks)))
    ok(sprintf("f %s: conservazione per cellula a mano (generati = catturati + scartati + persi)", A),
       all(o2$truth$cells$n_transcripts_generated == o2$truth$cells$n_captured + o2$truth$cells$n_thinned + o2$truth$cells$n_lost))
    ok(sprintf("f %s: profili della verita' = trascritti generati per cellula", A), all(Matrix::colSums(o1$truth$cell_counts) == o1$truth$cells$n_transcripts_generated))
    o4 <- simulate_visium_hd_pieces(Lr, seed = 8L, piece_um = Inf, frame_open = TRUE)
    ok(sprintf("a %s: un altro seme da' un altro risultato", A), !identical(o1$counts_2um, o4$counts_2um))
    if (A == "A6") {
      s5 <- simulate_visium_hd(Lr, seed = 7L, frame_open = TRUE)
      ok("g: con un pezzo solo le uscite hanno gli stessi elementi di simulate_visium_hd()", identical(names(o1)[1:5], names(s5)[1:5]) &&
         all(names(s5$truth) %in% names(o1$truth)) && identical(dim(o1$counts_2um), dim(s5$counts_2um)) && identical(dim(o1$counts_8um), dim(s5$counts_8um)) &&
         identical(rownames(o1$counts_2um), rownames(s5$counts_2um)) && all(names(s5$bins) %in% names(o1$bins)) && identical(o1$bins$bin_id, s5$bins$bin_id) &&
         identical(names(o1$truth$cells), names(s5$truth$cells)))
      ok("g: verita' geometrica uguale a build_hd_grid()", identical(o1$truth$overlap, s5$truth$overlap) && identical(o1$bins$tissue_area_um2, s5$bins$tissue_area_um2))
      o5 <- simulate_visium_hd_pieces(Lr, seed = 7L, piece_um = 112, frame_open = FALSE, workers = 2L)
      o6 <- simulate_visium_hd_pieces(Lr, seed = 7L, piece_um = Inf, frame_open = FALSE)
      ok("a A6: frame_open = FALSE, pezzo solo = pezzi di 112 um", same_out(o5, o6))
      lyr <- list(stripes = FALSE, background = FALSE, field = FALSE, second_tail = FALSE)
      ok("a A6: strati spenti, pezzo solo = pezzi di 56 um", same_out(simulate_visium_hd_pieces(Lr, seed = 7L, piece_um = Inf, layers = lyr, frame_open = TRUE),
                                                                   simulate_visium_hd_pieces(Lr, seed = 7L, piece_um = 56, layers = lyr, frame_open = TRUE, workers = 2L)))
    }
  }
  # due dischi del linfonodo: tessuto che non arriva ai lati, vetro fra i dischi, riquadro della riflessione dato
  xy <- expand.grid(i = 1:200, j = 1:200); md <- matrix(((xy$i - 50.5)^2 + (xy$j - 50.5)^2) <= 40^2 | ((xy$i - 150.5)^2 + (xy$j - 150.5)^2) <= 40^2, 200, 200)
  LD <- mk_layer("A4", md, 2, 13L)
  for (fo in c(FALSE, TRUE)) {
    fr <- if (fo) c(xmin = 0, ymin = 0, xmax = 400, ymax = 400) else NULL
    d1 <- simulate_visium_hd_pieces(LD, seed = 9L, piece_um = Inf, frame_open = fo, frame = fr)
    d2 <- simulate_visium_hd_pieces(LD, seed = 9L, piece_um = 56, frame_open = fo, frame = fr, workers = 3L)
    ok(sprintf("a dischi del linfonodo, frame_open = %s: pezzo solo = pezzi di 56 um (%d cellule)", fo, nrow(LD$cell_df)), same_out(d1, d2) && all(unlist(d2$metadata$checks)))
  }
  # seconda coda ripiegata sul riquadro dell'area: arriva sul vetro lontano dal tessuto (oltre la portata della prima coda)
  gg <- d1$metadata$grid; fr <- c(xmin = gg$x0 + 8, ymin = gg$y0 + 8, xmax = gg$x0 + gg$n_col * 2 - 8, ymax = gg$y0 + gg$n_row * 2 - 8)
  t_on <- simulate_visium_hd_pieces(LD, seed = 9L, piece_um = 112, frame_open = TRUE, frame = fr, workers = 2L)
  t_off <- simulate_visium_hd_pieces(LD, seed = 9L, piece_um = 112, frame_open = TRUE, frame = fr, workers = 2L, layers = list(second_tail = FALSE))
  ctr <- sf::st_as_sf(data.frame(x = t_on$bins$x, y = t_on$bins$y), coords = c("x", "y"))
  dist <- geos::geos_distance(geos::as_geos_geometry(sf::st_geometry(ctr)), geos::geos_unary_union(geos::geos_make_collection(geos::as_geos_geometry(sf::st_set_crs(LD$region_polygons, NA)))))
  far <- dist > 450
  u_on <- sum(Matrix::colSums(t_on$counts_2um)[far]); u_off <- sum(Matrix::colSums(t_off$counts_2um)[far])
  ok(sprintf("a: con frame_open la seconda coda arriva sul vetro a oltre 450 um dal tessuto (UMI %d con la coda, %d senza; %d quadratini)", u_on, u_off, sum(far)),
     sum(far) > 1000 && u_on > 10 * (u_off + 1))
  ok("a: lato dei pezzi non multiplo di 56 um rifiutato", inherits(try(simulate_visium_hd_pieces(LD, seed = 1L, piece_um = 100), silent = TRUE), "try-error"))
}
# --- (b) distanza dal bordo con l'alone = distanza sull'area intera dove e' < H; stessa scelta del modo ovunque
if ("b" %in% GR) {
  xy <- expand.grid(i = 1:240, j = 1:240); md <- matrix(((xy$i - 60.5)^2 + (xy$j - 60.5)^2) <= 45^2 | ((xy$i - 180.5)^2 + (xy$j - 170.5)^2) <= 50^2, 240, 240)
  R <- mk_layer("A4", md, 2, 21L)$region_polygons
  bb <- sf::st_bbox(R); g <- hd_grid_geometry(bb, 280, 2); step <- 0.5
  for (fo in c(FALSE, TRUE)) {
    env <- list(gx0 = g$x0, gy0 = g$y0, step = step, frame_open = fo, frame = c(xmin = bb[["xmin"]], ymin = bb[["ymin"]], xmax = bb[["xmax"]], ymax = bb[["ymax"]]),
                dom = c(q0 = as.integer(floor((bb[["xmin"]] - g$x0) / step)) - 1L, q1 = as.integer(ceiling((bb[["xmax"]] - g$x0) / step)) + 1L,
                        p0 = as.integer(floor((bb[["ymin"]] - g$y0) / step)) - 1L, p1 = as.integer(ceiling((bb[["ymax"]] - g$y0) / step)) + 1L))
    GG <- geos::as_geos_geometry(sf::st_set_crs(R, NA))
    P0 <- env$dom[["p0"]] - 400L; Q0 <- env$dom[["q0"]] - 400L; nr <- env$dom[["p1"]] - env$dom[["p0"]] + 800L; nc <- env$dom[["q1"]] - env$dom[["q0"]] + 800L
    Dg <- .p27_edt_window(GG, P0, nr, Q0, nc, env)                         # area intera
    Dw <- .p27_edt_window(GG, P0 + 500L, 300L, Q0 + 600L, 260L, env)        # finestra interna
    sub <- Dg[500L + seq_len(300L), 600L + seq_len(260L)]
    H <- 40                                                                # um: alone della finestra = distanza dal suo lato
    pg <- outer(pmin(seq_len(300L) - 1L, 300L - seq_len(300L)), pmin(seq_len(260L) - 1L, 260L - seq_len(260L)), pmin) * step
    inner <- sub < pg
    ok(sprintf("b frame_open = %s: dove la distanza vera e' minore della distanza dal lato della finestra, finestra = area intera (%d sottopixel)", fo, sum(inner)),
       sum(inner) > 1000 && identical(Dw[inner], sub[inner]))
    ok(sprintf("b frame_open = %s: altrove la finestra non sottostima la distanza", fo), all(Dw[!inner] >= sub[!inner] | Dw[!inner] >= pg[!inner]))
    ok(sprintf("b frame_open = %s: a pg ≥ 25 D la probabilita' della zona di bordo e' sotto il minimo del generatore", fo), exp(-25) < 2.3e-10 / 2)
    # il generatore L'Ecuyer-CMRG non estrae mai meno di 2^-32 circa: verifica su 1e7 estrazioni
    RNGkind("L'Ecuyer-CMRG"); set.seed(3); u <- runif(1e7); RNGkind("Mersenne-Twister")
    ok("b: minimo di 1e7 estrazioni uniformi > exp(-25)", min(u) > exp(-25))
  }
}
# --- (c) campo di variabilita' fine a finestre = campo sulla finestra grande
if ("c" %in% GR) {
  g <- hd_grid_geometry(c(xmin = 0, ymin = 0, xmax = 300, ymax = 300), 280, 2); L <- .p27_layout(g, Inf); tp <- 112L
  st <- .p27_streams(99L, L$n_tiles)
  for (A in c("A1", "A5", "A6")) {
    pp <- transcript_params(A)
    F0 <- .p27_field_window(500L, 700L, 480L, 650L, pp$field_s2, pp$field_L, 0.5, st, L, tp)
    F1 <- .p27_field_window(620L, 211L, 555L, 300L, pp$field_s2, pp$field_L, 0.5, st, L, tp)
    ok(sprintf("c %s: campo su una finestra = campo sulla finestra grande negli stessi sottopixel", A), identical(F1, F0[120L + seq_len(211L), 75L + seq_len(300L)]))
    ok(sprintf("c %s: media 1 e varianza vicina a s2 (%.3f contro %.3f)", A, stats::var(as.vector(F0)), pp$field_s2), abs(mean(F0) - 1) < 0.05 && abs(stats::var(as.vector(F0)) / pp$field_s2 - 1) < 0.25)
  }
}
# --- (d) ripiego della triangolazione (BL-241)
if ("d" %in% GR) {
  fx <- readRDS("R/testing/fixtures/S2.7_cell624286.rds")
  Tt <- geos::as_geos_geometry(fx$territory); Nn <- geos::as_geos_geometry(fx$nucleus); Cy <- geos::geos_difference(Tt, Nn)
  ok("d: precondizione, GEOS fallisce sul citoplasma originale della cellula 624286", inherits(try(geos::geos_constrained_delaunay_triangles(Cy), silent = TRUE), "try-error"))
  tr <- .s23_triangles(Cy)
  ok(sprintf("d: con il ripiego si triangola (%d triangoli), ripiego segnalato", nrow(tr)), nrow(tr) > 0 && identical(attr(tr, "fallback"), 1L))
  ok(sprintf("d: area dei triangoli entro 1e-3 um^2 dall'area del citoplasma (scarto %.2e)", abs(sum(tr$area) - geos::geos_area(Cy))), abs(sum(tr$area) - geos::geos_area(Cy)) < 1e-3)
  Lr <- mk_layer("A2", matrix(TRUE, 60, 60), 2, 5L)
  G0 <- geos::geos_difference(geos::as_geos_geometry(sf::st_set_crs(Lr$cell_territories, NA)), geos::as_geos_geometry(sf::st_set_crs(Lr$cell_nuclei, NA)))
  t0 <- .s23_triangles(G0); t1 <- .s23_triangles(c(G0, Cy))
  ok("d: le altre geometrie hanno gli stessi triangoli di prima", identical(as.matrix(t1[t1$cell <= length(G0), 1:8]), as.matrix(t0[, 1:8])) && identical(attr(t1, "fallback"), length(G0) + 1L))
  ok("d: senza geometrie degeneri nessun ripiego", length(attr(t0, "fallback")) == 0L)
}
# --- (e) tetto della sovra-generazione contro Monte Carlo
if ("e" %in% GR) {
  set.seed(20261222)
  for (A in c("A1", "A2", "A4", "A5", "A6")) {
    cp <- capture_params(A); Ms <- stripe_ceiling(cp)
    n <- 4e6; a <- rnorm(n, 0, cp$sd_alpha); b <- rnorm(n, 0, cp$sd_beta)
    odr <- runif(n) < 0.5; odc <- runif(n) < 0.5
    lr <- ifelse(odr, log(cp$parity_row), 0); lc <- ifelse(odc, log(cp$parity_col), 0)
    if (!is.null(cp$parity_field)) { lr <- lr + ifelse(odr, rnorm(n, 0, cp$parity_field$sd_row), 0); lc <- lc + ifelse(odc, rnorm(n, 0, cp$parity_field$sd_col), 0) }
    cc <- exp(a + b + lr + lc); mexp <- exp((cp$sd_alpha^2 + cp$sd_beta^2) / 2) * mean(exp(lr)) * mean(exp(lc)) / mean(exp(lr + lc)) * mean(cc)
    k <- sum(cc / mean(cc) > Ms)
    ok(sprintf("e %s: tetto %.4f, superamenti Monte Carlo %d su 4e6 (atteso 40 +- 25)", A, Ms, k), abs(k - 40) <= 25)
  }
  ok("e: tetto monotono in q", stripe_ceiling(capture_params("A1"), 1e-6) > stripe_ceiling(capture_params("A1"), 1e-5))
}
# --- (h) anomalia A-1: chiave unica oltre 2^53 nella funzione di S2.3 (corretta) e chiavi esatte nella funzione nuova
if ("h" %in% GR) {
  key <- 8e12; G <- 18132
  ok("h: precondizione, con key * (G + 1) + gene due geni vicini collidono a 8e12", (key * (G + 1) + 1) == (key * (G + 1) + 2))
  bi <- c(5L, 5L, 5L, 6L); ci <- c(9L, 9L, 9L, 9L); ge <- c(1L, 2L, 2L, 2L)
  ok("h: gruppi esatti della funzione nuova", identical(.p27_runs(bi, ci, ge), c(1L, 3L, 4L)))
  src <- paste(readLines(file.path(RDIR, "05b3_place_and_move_transcripts.R")), collapse = "\n")
  ok("h: place_and_move_transcripts() non usa piu' la chiave unica", !grepl("k3 <- key * (nrow(X) + 1) + gi", src, fixed = TRUE))
}
cat(sprintf("\n%d PASS, %d FAIL\n", np, nf))
if (nf > 0L) quit(status = 1L)
