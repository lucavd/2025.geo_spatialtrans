# R/05b3_place_and_move_transcripts.R - Step 2, passi 3 e 4: dove nascono i trascritti dentro la cellula e come si spostano.
# Documento di progetto dello Step 2 (docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md), sezioni 5, 6 e 9; sessione S2.3
# (docs/steps/S2.3.md). Regola della zona di bordo (DEC-039), campo di variabilita' fine (DEC-040).

# Tabella Z -> F del campo (stessa di tools/S2_3_common.py): passo 1e-3 in [-8.5, 8.5], interpolazione lineare.
.s23_ztab <- round(seq(-8500, 8500) * 1e-3, 3)

#' Parametri della posizione e dello spostamento dei trascritti per tessuto
#'
#' Valori della sezione 9 del documento di progetto dello Step 2, con tutte le cifre dei file che li hanno prodotti:
#' spostamento con la regola della zona di bordo (`S2.0/run1_wfix/params_flat.json`, regola C2; il cuore usa la mediana delle altre
#' sezioni, DEC-037; B-379...B-383), rapporto nucleo/citoplasma r (mediana delle 5 aree, `results/S2.0/wfix/S2.0_inner_summary.csv`,
#' B-384...B-394), fondo sul vetro lontano b (B-342...B-366, usato dallo strato di cattura), campo di variabilita' fine s2, L (S2.3, regola R-1).
#' @param archetype `"A1"` ... `"A6"`.
#' @return lista `r`, `sigma_in`, `sigma_e`, `D`, `w`, `ell` (um), `b` (UMI per quadratino), `field_s2`, `field_L` (um), `archetype`.
#' @export
transcript_params <- function(archetype) {
  tab <- list(
    A1 = list(r = 0.7459677413764999, sigma_e = 6.225411090229985, D = 24.74688933131723, w = 0.08719108841155758, ell = 32.85247983364391, b = 0.16145222298215922),
    A2 = list(r = 0.7546137105976738, sigma_e = 4.723658754203147, D = 8.936788049047061, w = 0.039734698098957445, ell = 24.730357429038083, b = 0.19714853642503777),
    A3 = list(r = 0.8207102541658515, sigma_e = 4.723658754203147, D = 8.936788049047061, w = 0.039734698098957445, ell = 24.730357429038083, b = 0.19714853642503777),
    A4 = list(r = 0.5458166288112142, sigma_e = 4.485594647990861, D = 22.04369695886815, w = 0.10882106662490754, ell = 41.50177210774394, b = 0.04195920540378664),
    A5 = list(r = 1.4834764850727304, sigma_e = 4.994182159099694, D = 6.860327792916313, w = 0.11230122544048957, ell = 33.0702155086906, b = 0.3086890909119761),
    A6 = list(r = 2.605117512918027, sigma_e = 4.85892045665142, D = 15.490242503957607, w = 0.09800607751823256, ell = 32.96134767116726, b = 0.1171531387445022))
  if (!archetype %in% names(tab)) stop("transcript_params(): archetipo sconosciuto: ", archetype, call. = FALSE)
  p <- tab[[archetype]]
  # campo di variabilita' fine: regola R-1 di S2.3, taratura sul nucleo (results/S2.3/S2.3_R1_choice.csv)
  fld <- list(A1 = c(0.5, 0.6375), A2 = c(0.5, 0.5), A3 = c(0.3535533905932738, 0.85), A4 = c(0.29730177875068026, 1.275),
              A5 = c(0.07432544468767006, 1.275), A6 = c(0.8408964152537145, 0.885))
  p$sigma_in <- 0.7; p$field_s2 <- fld[[archetype]][1]; p$field_L <- fld[[archetype]][2]; p$archetype <- archetype
  p
}

.s23_check_params <- function(p, field) {
  need <- c("r", "sigma_in", "sigma_e", "D", "w", "ell")
  miss <- setdiff(need, names(p)); if (length(miss)) stop("place_and_move_transcripts(): parametri mancanti: ", paste(miss, collapse = ", "), call. = FALSE)
  ok1 <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  for (k in c("r", "sigma_in", "sigma_e", "D", "ell")) if (!ok1(p[[k]]) || p[[k]] <= 0) stop(sprintf("place_and_move_transcripts(): `%s` deve essere un numero > 0.", k), call. = FALSE)
  if (!ok1(p$w) || p$w < 0 || p$w > 1) stop("place_and_move_transcripts(): `w` deve stare in [0, 1].", call. = FALSE)
  # seconda coda (S2.4, BL-221): facoltativa; assente o w2 = 0 = comportamento di S2.3
  if (!is.null(p$w2)) {
    if (!ok1(p$w2) || p$w2 < 0 || p$w2 >= 1) stop("place_and_move_transcripts(): `w2` deve stare in [0, 1).", call. = FALSE)
    if (p$w2 > 0 && (!ok1(p$ell2) || p$ell2 <= 0)) stop("place_and_move_transcripts(): `ell2` deve essere un numero > 0.", call. = FALSE)
  }
  if (isTRUE(field)) {
    if (!ok1(p$field_s2) || p$field_s2 < 0) stop("place_and_move_transcripts(): campo acceso ma `field_s2` non valido (campo non tarato?).", call. = FALSE)
    if (!ok1(p$field_L) || p$field_L <= 0) stop("place_and_move_transcripts(): campo acceso ma `field_L` non valido.", call. = FALSE)
  }
  invisible(TRUE)
}

# Riflessione ripetuta nell'intervallo [a, b] (S2.4): x' = a + W - |((x - a) mod 2W) - W|, W = b - a.
.s24_fold <- function(x, a, b) { W <- b - a; t <- (x - a) %% (2 * W); a + W - abs(t - W) }

# Triangoli (triangolazione vincolata di GEOS) di un vettore di poligoni: una riga per triangolo, con l'indice del poligono.
.s23_triangles <- function(g) {
  tri <- geos::geos_constrained_delaunay_triangles(g)
  k <- wk::wk_coords(tri)
  if (!nrow(k)) return(data.frame(cell = integer(0), ax = numeric(0), ay = numeric(0), bx = numeric(0), by = numeric(0), cx = numeric(0), cy = numeric(0), area = numeric(0)))
  stopifnot(nrow(k) %% 4L == 0L)
  i1 <- seq(1L, nrow(k), by = 4L)
  t <- data.frame(cell = k$feature_id[i1], ax = k$x[i1], ay = k$y[i1], bx = k$x[i1 + 1L], by = k$y[i1 + 1L], cx = k$x[i1 + 2L], cy = k$y[i1 + 2L])
  t$area <- abs((t$bx - t$ax) * (t$cy - t$ay) - (t$cx - t$ax) * (t$by - t$ay)) / 2
  t[t$area > 0, , drop = FALSE]
}

# Maschera "centro del sottopixel nei poligoni" (regola pari-dispari, scansione per righe): matrice nr x nc, riga = asse y.
.s23_raster_center_in <- function(g, x0, y0, step, nr, nc) {
  k <- wk::wk_coords(geos::geos_boundary(g))
  M <- matrix(FALSE, nr, nc)
  if (!nrow(k)) return(M)
  # lati consecutivi dello stesso anello (ring_id unico per anello)
  n <- nrow(k); same <- k$ring_id[-1] == k$ring_id[-n] & k$feature_id[-1] == k$feature_id[-n] & k$part_id[-1] == k$part_id[-n]
  x1 <- k$x[-n][same]; y1 <- k$y[-n][same]; x2 <- k$x[-1][same]; y2 <- k$y[-1][same]
  yc0 <- y0 + 0.5 * step
  lo <- pmin(y1, y2); hi <- pmax(y1, y2)
  r_lo <- pmax(0L, as.integer(ceiling((lo - yc0) / step))); r_hi <- pmin(nr - 1L, as.integer(ceiling((hi - yc0) / step)) - 1L)
  keep <- r_hi >= r_lo & y1 != y2
  x1 <- x1[keep]; y1 <- y1[keep]; x2 <- x2[keep]; y2 <- y2[keep]; r_lo <- r_lo[keep]; r_hi <- r_hi[keep]
  len <- r_hi - r_lo + 1L
  e <- rep(seq_along(len), len); rr <- r_lo[e] + sequence(len) - 1L
  yc <- yc0 + rr * step
  # attraversa se (y1 <= yc) != (y2 <= yc): con lo <= yc < hi vale sempre
  xs <- x1[e] + (yc - y1[e]) * (x2[e] - x1[e]) / (y2[e] - y1[e])
  o <- order(rr, xs); rr <- rr[o]; xs <- xs[o]
  if (length(rr) %% 2L) stop(".s23_raster_center_in(): numero dispari di attraversamenti (poligono non chiuso?).", call. = FALSE)
  a <- seq(1L, length(rr), by = 2L)
  stopifnot(all(rr[a] == rr[a + 1L]))
  c_lo <- pmax(0L, as.integer(ceiling((xs[a] - x0) / step - 0.5))); c_hi <- pmin(nc - 1L, as.integer(ceiling((xs[a + 1L] - x0) / step - 0.5)) - 1L)
  ok <- c_hi >= c_lo; r_ <- rr[a][ok]; c_lo <- c_lo[ok]; c_hi <- c_hi[ok]
  nb <- nr * (nc + 1L)
  D <- tabulate((c_lo) * nr + r_ + 1L, nb) - tabulate((c_hi + 1L) * nr + r_ + 1L, nb)
  D <- matrix(D, nr, nc + 1L)
  S <- t(apply(D, 1, cumsum))
  S[, seq_len(nc), drop = FALSE] > 0L
}

# Riflessione ai bordi come numpy 'reflect' (d c b | a b c d | c b a), p pixel per lato.
.s23_reflect_pad <- function(M, p) {
  ri <- function(n) { i <- seq(-p, n - 1L + p); m <- 2L * (n - 1L); if (m == 0L) return(rep(1L, length(i))); i <- abs(i) %% m; i <- ifelse(i > n - 1L, m - i, i); i + 1L }
  M[ri(nrow(M)), ri(ncol(M)), drop = FALSE]
}

# Distanza euclidea esatta (um) dal centro di ogni sottopixel di tessuto al centro del sottopixel di non tessuto piu' vicino.
.s23_edt <- function(M, step) {
  if (all(M)) return(matrix(Inf, nrow(M), ncol(M)))
  d <- as.matrix(imager::distance_transform(imager::as.cimg(1 * !M), value = 1, metric = 2))
  stopifnot(identical(dim(d), dim(M)))
  d * step
}

#' Campo di variabilita' fine (design sezione 5, DEC-040)
#'
#' Z = rumore bianco gaussiano convoluto con una gaussiana di deviazione standard L (pesi troncati a 4 deviazioni standard), riportato
#' a varianza 1; F = quantile della gamma di media 1 e varianza s2 in Phi(Z) (interpolazione lineare sulla tabella `.s23_ztab`).
#' Usa il generatore casuale corrente.
#' @param nr,nc dimensioni della sottogriglia.
#' @param s2 varianza puntuale di F (0 = campo spento, F = 1).
#' @param L lunghezza (um) della gaussiana.
#' @param step lato del sottopixel (um).
#' @return matrice nr x nc.
#' @export
s23_gamma_field <- function(nr, nc, s2, L, step = 0.5) {
  if (s2 <= 0) return(matrix(1, nr, nc))
  Z <- .s23_gaussian_z(nr, nc, L, step)
  tabF <- .s23_gamma_table(s2)
  Z <- pmin(pmax(Z, .s23_ztab[1]), .s23_ztab[length(.s23_ztab)])
  matrix(stats::approx(.s23_ztab, tabF, xout = as.vector(Z))$y, nr, nc)
}

# Tabella Z -> F: per Z >= 0 con le code superiori (precise vicino a 1, finite fino a Z = 8.5), come tools/S2_3_common.py.
.s23_gamma_table <- function(s2) {
  lo <- stats::qgamma(stats::pnorm(.s23_ztab), shape = 1 / s2, scale = s2)
  hi <- stats::qgamma(stats::pnorm(.s23_ztab, lower.tail = FALSE), shape = 1 / s2, scale = s2, lower.tail = FALSE)
  ifelse(.s23_ztab >= 0, hi, lo)
}

# Pesi della gaussiana di deviazione standard L (pixel L / step), troncata a 4 deviazioni standard (come tools/S2_3_common.py).
.s23_taps <- function(L, step = 0.5) { sd <- L / step; h <- as.integer(ceiling(4 * sd)); exp(-0.5 * (seq(-h, h) / sd)^2) }

# Z gaussiano a varianza 1: rumore bianco su un riquadro allargato, convoluzione separabile, riportato a varianza 1.
.s23_gaussian_z <- function(nr, nc, L, step = 0.5) {
  g <- .s23_taps(L, step); h <- (length(g) - 1L) %/% 2L
  W <- matrix(stats::rnorm((nr + 2L * h) * (nc + 2L * h)), nr + 2L * h, nc + 2L * h)
  Tm <- matrix(0, nr, nc + 2L * h)
  for (a in seq_along(g)) Tm <- Tm + g[a] * W[a - 1L + seq_len(nr), , drop = FALSE]
  Z <- matrix(0, nr, nc)
  for (a in seq_along(g)) Z <- Z + g[a] * Tm[, a - 1L + seq_len(nc), drop = FALSE]
  Z / sum(g^2)
}

#' Posizione e spostamento dei trascritti (Step 2, passi 3 e 4)
#'
#' Ogni trascritto di una cellula nasce nel suo territorio con densita' r volte piu' alta nel nucleo che nel citoplasma, moltiplicata per
#' il campo F (estrazione esatta: triangolazione vincolata, punto uniforme nel triangolo, accettazione con F / max F della cellula);
#' poi si sposta (regola della zona di bordo, DEC-039): con probabilita' w coda esponenziale (raggio Gamma(2, ell), angolo uniforme),
#' altrimenti gaussiana isotropa con sigma_e per asse con probabilita' exp(-d / D), con sigma_in con probabilita' 1 - exp(-d / D), dove d e'
#' la distanza del punto di nascita dal bordo del tessuto (trasformata di distanza della maschera del tessuto sulla sottogriglia di 0.5 um,
#' come il prototipo di S2.0); infine si conta nel quadratino di 2 um in cui atterra. I trascritti che escono dalla griglia sono persi.
#' Il gene non entra nella posizione. Cellule a blocchi di `chunk_size`, un flusso casuale L'Ecuyer-CMRG per blocco e uno per il campo:
#' il risultato non dipende da `workers`; lo stato del generatore dell'utente e' ripristinato all'uscita.
#' @param cell_layer uscita di `place_cells()` (territori, nuclei, regioni, `cell_df`, `metadata$pixel_size_um`).
#' @param hd uscita di `build_hd_grid()` sullo stesso livello.
#' @param cell_profiles uscita di `simulate_cell_profiles()` (colonne di `counts` nell'ordine di `cell_df$cell_id`).
#' @param params lista come `transcript_params()`.
#' @param seed intero.
#' @param field `TRUE` per lo strato di variabilita' fine (servono `params$field_s2`, `params$field_L`), `FALSE` per spegnerlo.
#' @param frame_open `TRUE`: il riquadro delle regioni non conta come bordo (maschera riflessa oltre il riquadro per `pad_um`), per aree
#'   ritagliate da tessuto che continua.
#' @param field_matrix campo gia' realizzato sulla sottogriglia (per i test); sostituisce quello generato.
#' @param step lato del sottopixel (um) di campo e maschera.
#' @param pad_um riflessione della maschera con `frame_open = TRUE`.
#' @param chunk_size cellule per blocco (fa parte della definizione del risultato).
#' @param workers processi (`parallel::mclapply`).
#' @param return_points `TRUE` aggiunge `points` (nascita, distanza, modo, arrivo per trascritto) per i test.
#' @return lista `transcripts` (data.frame `bin_id`, `cell_id`, `gene`, `n`), `origin` (`bin_id`, `cell_id`, `n`, `n_born_nucleus`),
#'   `lost` (`cell_id`, `n_lost`, `n_lost_born_nucleus`), `genes` (identificativi), `field_grid` (origine, passo, dimensioni), `metadata`; con `return_points`, `points`.
#' @export
place_and_move_transcripts <- function(cell_layer, hd, cell_profiles, params, seed, field = TRUE, frame_open = FALSE, field_matrix = NULL,
                                       step = 0.5, pad_um = 100, chunk_size = 2000L, workers = 1L, return_points = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  px <- cell_layer$metadata$pixel_size_um
  if (is.null(px) || !is.numeric(px) || length(px) != 1L || !is.finite(px) || px <= 0)
    stop("place_and_move_transcripts(): `cell_layer$metadata$pixel_size_um` manca o non e' valido (C-A): nessun default.", call. = FALSE)
  if (missing(seed) || length(seed) != 1L || !is.finite(seed)) stop("place_and_move_transcripts(): serve un `seed`.", call. = FALSE)
  .s23_check_params(params, field && is.null(field_matrix))
  cd <- cell_layer$cell_df; n <- nrow(cd)
  X <- cell_profiles$counts
  if (is.null(X) || ncol(X) != n || !identical(colnames(X), as.character(cd$cell_id)))
    stop("place_and_move_transcripts(): `cell_profiles$counts` non allineato a `cell_df$cell_id`.", call. = FALSE)
  Tt <- geos::as_geos_geometry(sf::st_set_crs(cell_layer$cell_territories, NA)); Nn <- geos::as_geos_geometry(sf::st_set_crs(cell_layer$cell_nuclei, NA))
  if (length(Tt) != n || length(Nn) != n) stop("place_and_move_transcripts(): territori, nuclei e cell_df di lunghezza diversa.", call. = FALSE)
  if (any(geos::geos_is_empty(Tt))) stop("place_and_move_transcripts(): territori vuoti.", call. = FALSE)
  G <- geos::as_geos_geometry(sf::st_set_crs(cell_layer$region_polygons, NA))
  r <- params$r
  # (1) geometria: triangoli di nucleo e citoplasma, aree esatte
  Cy <- geos::geos_difference(Tt, Nn)
  A_n <- geos::geos_area(Nn); A_n[geos::geos_is_empty(Nn)] <- 0; A_c <- geos::geos_area(Cy); A_c[geos::geos_is_empty(Cy)] <- 0
  trN <- .s23_triangles(Nn); trC <- .s23_triangles(Cy)
  pn <- ifelse(A_n > 0, r * A_n / (r * A_n + A_c), 0)
  cumN <- cumsum(trN$area); cumC <- cumsum(trC$area)
  firstN <- match(seq_len(n), trN$cell); lastN <- length(trN$cell) - match(seq_len(n), rev(trN$cell)) + 1L
  firstC <- match(seq_len(n), trC$cell); lastC <- length(trC$cell) - match(seq_len(n), rev(trC$cell)) + 1L
  startN <- ifelse(is.na(firstN), NA, c(0, cumN)[firstN]); startC <- ifelse(is.na(firstC), NA, c(0, cumC)[firstC])
  totN <- ifelse(is.na(firstN), 0, cumN[pmax(1L, lastN)] - startN); totC <- ifelse(is.na(firstC), 0, cumC[pmax(1L, lastC)] - startC)
  pn[is.na(firstC) | totC <= 0] <- 1; pn[is.na(firstN) | totN <= 0] <- 0
  Nc0 <- as.integer(round(Matrix::colSums(X)))
  dead <- which(Nc0 > 0L & totN <= 0 & totC <= 0)
  if (length(dead)) stop(sprintf("place_and_move_transcripts(): %d cellule con trascritti e area nulla (cell_id %s).", length(dead), paste(utils::head(cd$cell_id[dead], 5), collapse = ", ")), call. = FALSE)
  # (2) sottogriglia allineata alla griglia dei quadratini, sul riquadro delle regioni e dei territori
  bb <- sf::st_bbox(c(sf::st_set_crs(cell_layer$cell_territories, NA), sf::st_set_crs(cell_layer$region_polygons, NA)))
  gx0 <- hd$grid$x0; gy0 <- hd$grid$y0
  q0 <- as.integer(floor((bb[["xmin"]] - gx0) / step)) - 1L; q1 <- as.integer(ceiling((bb[["xmax"]] - gx0) / step)) + 1L
  p0 <- as.integer(floor((bb[["ymin"]] - gy0) / step)) - 1L; p1 <- as.integer(ceiling((bb[["ymax"]] - gy0) / step)) + 1L
  fx0 <- gx0 + q0 * step; fy0 <- gy0 + p0 * step; nr <- p1 - p0; nc <- q1 - q0
  Mt <- .s23_raster_center_in(G, fx0, fy0, step, nr, nc)
  if (frame_open) {
    # il riquadro delle regioni non e' bordo: i sottopixel con il centro fuori dal riquadro valgono tessuto, e la maschera si allarga di
    # pad_um di tessuto (revisione G-5, DSC-1: prima i sottopixel tagliati dal lato del riquadro restavano vetro)
    xc <- fx0 + (seq_len(nc) - 0.5) * step; yc <- fy0 + (seq_len(nr) - 0.5) * step
    out_box <- outer(yc < bb[["ymin"]] | yc > bb[["ymax"]], xc < bb[["xmin"]] | xc > bb[["xmax"]], `|`)
    Mo <- Mt | out_box
    pp <- as.integer(round(pad_um / step))
    Mp <- matrix(TRUE, nr + 2L * pp, nc + 2L * pp); Mp[pp + seq_len(nr), pp + seq_len(nc)] <- Mo
    Dm <- .s23_edt(Mp, step)[pp + seq_len(nr), pp + seq_len(nc), drop = FALSE]
  } else {
    Dm <- .s23_edt(Mt, step)
  }
  # flussi casuali: uno per il campo, uno per blocco
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  s <- get(".Random.seed", envir = globalenv())
  Fm <- if (!is.null(field_matrix)) {
    if (!identical(dim(field_matrix), c(nr, nc))) stop(sprintf("place_and_move_transcripts(): `field_matrix` deve essere %d x %d.", nr, nc), call. = FALSE)
    field_matrix
  } else if (isTRUE(field)) { assign(".Random.seed", s, envir = globalenv()); s23_gamma_field(nr, nc, params$field_s2, params$field_L, step) } else NULL
  # massimo del campo nel riquadro di ogni cellula
  Fmax <- rep(1, n)
  if (!is.null(Fm)) {
    allt <- rbind(trN[, c("cell", "ax", "ay", "bx", "by", "cx", "cy")], trC[, c("cell", "ax", "ay", "bx", "by", "cx", "cy")])
    xmn <- tapply(pmin(allt$ax, allt$bx, allt$cx), allt$cell, min); xmx <- tapply(pmax(allt$ax, allt$bx, allt$cx), allt$cell, max)
    ymn <- tapply(pmin(allt$ay, allt$by, allt$cy), allt$cell, min); ymx <- tapply(pmax(allt$ay, allt$by, allt$cy), allt$cell, max)
    ids <- as.integer(names(xmn))
    a0 <- pmax(1L, as.integer(floor((ymn - fy0) / step)) + 1L); a1 <- pmin(nr, as.integer(floor((ymx - fy0) / step)) + 1L)
    b0 <- pmax(1L, as.integer(floor((xmn - fx0) / step)) + 1L); b1 <- pmin(nc, as.integer(floor((xmx - fx0) / step)) + 1L)
    for (k in seq_along(ids)) Fmax[ids[k]] <- max(Fm[a0[k]:a1[k], b0[k]:b1[k]])
  }
  Nc <- Nc0
  w2 <- if (is.null(params$w2)) 0 else params$w2
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  streams <- vector("list", length(chunks)); s <- parallel::nextRNGStream(s)
  for (k in seq_along(chunks)) { streams[[k]] <- s; s <- parallel::nextRNGStream(s) }
  gn <- hd$grid$n_col; gnr <- hd$grid$n_row; bin_um <- hd$grid$bin_um
  Xp <- X@p; Xi <- X@i; Xx <- as.integer(round(X@x))
  NX <- as.list(trN[, c("ax", "ay", "bx", "by", "cx", "cy")]); CX <- as.list(trC[, c("ax", "ay", "bx", "by", "cx", "cy")])
  one_chunk <- function(k) {
    assign(".Random.seed", streams[[k]], envir = globalenv())
    idx <- chunks[[k]]; need <- Nc[idx]; names(need) <- idx
    acc_cell <- list(); acc_x <- list(); acc_y <- list(); acc_nuc <- list(); it <- 0L
    while (any(need > 0L)) {
      it <- it + 1L; if (it > 200L) stop("place_and_move_transcripts(): accettazione non converge.", call. = FALSE)
      cc <- idx[need > 0L]; m <- as.integer(ceiling(need[need > 0L] * Fmax[cc] * 1.2)) + 2L
      cell <- rep(cc, m); nm <- length(cell)
      isn <- stats::runif(nm) < pn[cell]
      u <- stats::runif(nm)
      ti <- integer(nm)
      if (any(isn)) { cn <- cell[isn]; tg <- startN[cn] + u[isn] * totN[cn]; j <- findInterval(tg, cumN, left.open = TRUE) + 1L; ti[isn] <- pmin(pmax(j, firstN[cn]), lastN[cn]) }
      if (any(!isn)) { cn <- cell[!isn]; tg <- startC[cn] + u[!isn] * totC[cn]; j <- findInterval(tg, cumC, left.open = TRUE) + 1L; ti[!isn] <- pmin(pmax(j, firstC[cn]), lastC[cn]) }
      r1 <- stats::runif(nm); r2 <- stats::runif(nm); fl <- r1 + r2 > 1; r1[fl] <- 1 - r1[fl]; r2[fl] <- 1 - r2[fl]
      x <- numeric(nm); y <- numeric(nm)
      if (any(isn)) { j <- ti[isn]; x[isn] <- NX$ax[j] + r1[isn] * (NX$bx[j] - NX$ax[j]) + r2[isn] * (NX$cx[j] - NX$ax[j]); y[isn] <- NX$ay[j] + r1[isn] * (NX$by[j] - NX$ay[j]) + r2[isn] * (NX$cy[j] - NX$ay[j]) }
      if (any(!isn)) { j <- ti[!isn]; x[!isn] <- CX$ax[j] + r1[!isn] * (CX$bx[j] - CX$ax[j]) + r2[!isn] * (CX$cx[j] - CX$ax[j]); y[!isn] <- CX$ay[j] + r1[!isn] * (CX$by[j] - CX$ay[j]) + r2[!isn] * (CX$cy[j] - CX$ay[j]) }
      ua <- stats::runif(nm)
      ok <- if (is.null(Fm)) rep(TRUE, nm) else {
        pa <- pmin(nr, pmax(1L, as.integer(floor((y - fy0) / step)) + 1L)); qa <- pmin(nc, pmax(1L, as.integer(floor((x - fx0) / step)) + 1L))
        ua < Fm[cbind(pa, qa)] / Fmax[cell] }
      cell <- cell[ok]; x <- x[ok]; y <- y[ok]; isn <- isn[ok]
      rk <- sequence(rle(cell)$lengths)
      take <- rk <= need[as.character(cell)]
      acc_cell[[it]] <- cell[take]; acc_x[[it]] <- x[take]; acc_y[[it]] <- y[take]; acc_nuc[[it]] <- isn[take]
      got <- tabulate(match(cell[take], idx), length(idx)); need <- need - got
    }
    cell <- unlist(acc_cell); x <- unlist(acc_x); y <- unlist(acc_y); isn <- unlist(acc_nuc)
    o <- order(cell); cell <- cell[o]; x <- x[o]; y <- y[o]; isn <- isn[o]
    # (3) distanza dal bordo al punto di nascita e spostamento
    pa <- pmin(nr, pmax(1L, as.integer(floor((y - fy0) / step)) + 1L)); qa <- pmin(nc, pmax(1L, as.integer(floor((x - fx0) / step)) + 1L))
    d <- Dm[cbind(pa, qa)]
    nt <- length(cell); u1 <- stats::runif(nt); u2 <- stats::runif(nt)
    tail <- u1 < params$w; wide <- !tail & (u2 < exp(-d / params$D))
    mode <- ifelse(tail, 2L, ifelse(wide, 1L, 0L))
    dx <- numeric(nt); dy <- numeric(nt)
    nT <- sum(tail); if (nT) { rho <- stats::rgamma(nT, shape = 2, scale = params$ell); th <- stats::runif(nT, 0, 2 * pi); dx[tail] <- rho * cos(th); dy[tail] <- rho * sin(th) }
    nG <- sum(!tail); if (nG) { sg <- ifelse(wide[!tail], params$sigma_e, params$sigma_in); dx[!tail] <- stats::rnorm(nG) * sg; dy[!tail] <- stats::rnorm(nG) * sg }
    x1 <- x + dx; y1 <- y + dy
    # seconda coda (S2.4, BL-221, DEC-051 e regola R-3 di S2.4): con probabilita' w2 raggio Gamma(2, ell2) e angolo uniforme (modo 3), estratta
    # dopo le estrazioni di S2.3 (con w2 assente o 0 nessuna estrazione in piu': risultato identico a S2.3). Con frame_open l'arrivo si ripiega
    # nel riquadro per riflessione sui lati (scelta di Luca, 2026-11-06); altrimenti chi esce dalla griglia e' perso.
    if (w2 > 0) {
      sec <- stats::runif(nt) < w2; nS <- sum(sec)
      if (nS) {
        rho2 <- stats::rgamma(nS, shape = 2, scale = params$ell2); th2 <- stats::runif(nS, 0, 2 * pi)
        x1[sec] <- x[sec] + rho2 * cos(th2); y1[sec] <- y[sec] + rho2 * sin(th2)
        if (frame_open) { x1[sec] <- .s24_fold(x1[sec], bb[["xmin"]], bb[["xmax"]]); y1[sec] <- .s24_fold(y1[sec], bb[["ymin"]], bb[["ymax"]]) }
        mode[sec] <- 3L
      }
    }
    # (4) conteggio
    ii <- floor((y1 - gy0) / bin_um); jj <- floor((x1 - gx0) / bin_um)
    inside <- ii >= 0 & ii < gnr & jj >= 0 & jj < gn
    bin <- ifelse(inside, ii * gn + jj + 1, NA_real_)
    # geni nell'ordine dei trascritti della cellula (le posizioni sono estrazioni indipendenti e identiche)
    gene <- unlist(lapply(idx, function(cc) { a <- Xp[cc] + 1L; b <- Xp[cc + 1L]; if (b < a) integer(0) else rep(Xi[a:b] + 1L, Xx[a:b]) }))
    stopifnot(length(gene) == nt)
    pts <- if (return_points) data.frame(cell = cell, x0 = x, y0 = y, nucleus = isn, d = d, mode = mode, x1 = x1, y1 = y1, bin_id = bin, gene = gene) else NULL
    lost <- tabulate(match(cell[!inside], idx), length(idx)); lost_nuc <- tabulate(match(cell[!inside & isn], idx), length(idx))
    lost_gene <- tabulate(gene[!inside], nrow(X))
    ci <- cell[inside]; bi <- bin[inside]; gi <- gene[inside]; ni <- isn[inside]
    key <- (bi - 1) * n + (ci - 1)
    ko <- order(key, gi); key <- key[ko]; gi <- gi[ko]; ni <- ni[ko]
    k3 <- key * (nrow(X) + 1) + gi
    rl <- rle(k3); e3 <- cumsum(rl$lengths)
    tr <- data.frame(bin_id = as.integer(key[e3] %/% n + 1), cell = as.integer(key[e3] %% n + 1), gene = as.integer(gi[e3]), n = as.integer(rl$lengths))
    rk2 <- rle(key); e2 <- cumsum(rk2$lengths)
    nb <- diff(c(0L, cumsum(ni)[e2]))
    og <- data.frame(bin_id = as.integer(key[e2] %/% n + 1), cell = as.integer(key[e2] %% n + 1), n = as.integer(rk2$lengths), n_born_nucleus = as.integer(nb))
    list(tr = tr, og = og, lost = lost, lost_nuc = lost_nuc, lost_gene = lost_gene, pts = pts, n_gen = nt, it = it)
  }
  out <- if (workers > 1L) parallel::mclapply(seq_along(chunks), one_chunk, mc.cores = workers, mc.preschedule = FALSE) else lapply(seq_along(chunks), one_chunk)
  bad <- vapply(out, function(o) inherits(o, "try-error"), logical(1))
  if (any(bad)) stop("place_and_move_transcripts(): errore in un blocco: ", as.character(out[[which(bad)[1]]]), call. = FALSE)
  tr <- do.call(rbind, lapply(out, `[[`, "tr")); og <- do.call(rbind, lapply(out, `[[`, "og"))
  tr <- tr[order(tr$bin_id, tr$cell, tr$gene), ]; og <- og[order(og$bin_id, og$cell), ]
  lost <- unlist(lapply(out, `[[`, "lost")); lost_nuc <- unlist(lapply(out, `[[`, "lost_nuc"))
  tr$cell_id <- cd$cell_id[tr$cell]; og$cell_id <- cd$cell_id[og$cell]
  res <- list(transcripts = data.frame(bin_id = tr$bin_id, cell_id = tr$cell_id, gene = tr$gene, n = tr$n),
              origin = data.frame(bin_id = og$bin_id, cell_id = og$cell_id, n = og$n, n_born_nucleus = og$n_born_nucleus),
              lost = data.frame(cell_id = cd$cell_id, n_lost = as.integer(lost), n_lost_born_nucleus = as.integer(lost_nuc)),
              genes = rownames(X), lost_by_gene = Reduce(`+`, lapply(out, `[[`, "lost_gene")),
              field_grid = list(x0 = fx0, y0 = fy0, step = step, n_row = nr, n_col = nc),
              cells = data.frame(cell_id = cd$cell_id, area_nucleus_um2 = A_n, area_cytoplasm_um2 = A_c, p_nucleus_no_field = pn, n_transcripts = Nc),
              metadata = list(step = "S2.3", seed = seed, params = params, field = !is.null(Fm), frame_open = frame_open, pad_um = pad_um,
                              chunk_size = chunk_size, n_transcripts = sum(Nc), n_lost = sum(lost),
                              package_version = "geo.spatialtrans step2-visium-grid (S2.3)", elapsed_s = proc.time()[["elapsed"]] - t0))
  rownames(res$transcripts) <- NULL; rownames(res$origin) <- NULL
  if (return_points) { pts <- do.call(rbind, lapply(out, `[[`, "pts")); pts$cell_id <- cd$cell_id[pts$cell]; res$points <- pts; res$field_matrix <- Fm; res$distance_matrix <- Dm; res$tissue_matrix <- Mt }
  res
}
