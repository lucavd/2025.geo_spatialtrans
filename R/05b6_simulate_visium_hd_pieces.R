# R/05b6_simulate_visium_hd_pieces.R - Step 2, S2.7: simulazione a pezzi dentro il pacchetto.
# File dello step: docs/steps/S2.7.md; architettura DEC-074; tetto della sovra-generazione DEC-075; requisiti DEC-073.
# Riusa le funzioni di R/05b1...R/05b5: griglia (hd_grid_geometry), verita' geometrica (stessa aritmetica di hd_tissue_area e hd_overlap),
# strisce (simulate_stripes), profili (simulate_cell_profiles con i flussi per tessera), triangoli (.s23_triangles con il ripiego di BL-241),
# campo (stessa aritmetica di .s23_gaussian_z e s23_gamma_field), distanza dal bordo (.s23_edt), cattura (stessa logica di capture_and_bin).
# Tutte le estrazioni casuali sono legate a tessere fisse di 56 um ancorate all'origine della griglia: il risultato non dipende dal lato dei
# pezzi ne' dal numero di processi (G-1, G-2 di S2.7). Dipendenze: Matrix, sf, geos, wk, imager, parallel.

.P27_TILE_UM <- 56

#' Tetto della sovra-generazione delle strisce (DEC-075)
#'
#' Quantile `1 - q` della legge del fattore di striscia di un quadratino, c = p_r(i mod 2) p_c(j mod 2) exp(alpha_i + beta_j) / media:
#' alpha e beta gaussiani con le deviazioni standard di `capture_params()` (parte lenta compresa, se c'e'), campo della parita' gaussiano
#' sulle righe e colonne dispari dove c'e' (linfonodo), media della legge = 1. Calcolo esatto (somma di quattro code gaussiane, radice in
#' scala logaritmica).
#' @param params lista come `capture_params()`.
#' @param q probabilita' di superare il tetto (predefinito 1e-5, DEC-075).
#' @return numero > 1.
#' @export
stripe_ceiling <- function(params, q = 1e-5) {
  .s24_check_capture_params(params)
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0 || q >= 0.25) stop("stripe_ceiling(): `q` deve stare in (0, 0.25).", call. = FALSE)
  va <- params$sd_alpha^2; vb <- params$sd_beta^2
  if (!is.null(params$slow)) {
    if (!is.null(params$slow$sd_alpha)) va <- va + params$slow$sd_alpha^2
    if (!is.null(params$slow$sd_beta)) vb <- vb + params$slow$sd_beta^2
  }
  pf <- params$parity_field
  er <- c(0, if (is.null(pf)) 0 else pf$sd_row^2); ec <- c(0, if (is.null(pf)) 0 else pf$sd_col^2)
  lr <- c(0, log(params$parity_row)); lc <- c(0, log(params$parity_col))
  lm <- (va + vb) / 2 + log(mean(exp(lr + er / 2))) + log(mean(exp(lc + ec / 2)))      # log della media della legge
  ltail <- function(lt) {
    l <- c(); for (a in 1:2) for (b in 1:2) l <- c(l, stats::pnorm(lt + lm, mean = lr[a] + lc[b], sd = sqrt(va + vb + er[a] + ec[b]), lower.tail = FALSE, log.p = TRUE))
    mx <- max(l); mx + log(sum(exp(l - mx))) - log(4)
  }
  r <- stats::uniroot(function(lt) ltail(lt) - log(q), c(0, 30), tol = 1e-13)
  exp(r$root)
}

# flussi L'Ecuyer-CMRG 1...n dal seme (stessa costruzione di .s24_with_stream: seme rimescolato con Mersenne-Twister)
.p27_streams <- function(seed, n) .s24_with_stream(seed, n, function(st) st)

# semi delle sei famiglie di flussi, estratti con Mersenne-Twister dal seme dell'utente
.p27_draw_seeds <- function(seed) {
  old_kind <- RNGkind(); had <- exists(".Random.seed", envir = globalenv())
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had) assign(".Random.seed", old, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  s <- sample.int(1e8, 6L)
  c(stripes = s[1], profiles = s[2], place = s[3], capture = s[4], field = s[5], background = s[6])
}

# reticolo delle tessere e dei pezzi sulla griglia (indici 0-based; tessera di kS quadratini, pezzo di kP tessere per lato)
.p27_layout <- function(grid, piece_um, tile_um = .P27_TILE_UM) {
  kS <- as.integer(round(tile_um / grid$bin_um))
  if (abs(kS * grid$bin_um - tile_um) > 1e-9 || kS %% grid$block != 0L) stop("simulate_visium_hd_pieces(): tessere non multiple dei blocchi di 8 um.", call. = FALSE)
  ntr <- as.integer(ceiling(grid$n_row / kS)); ntc <- as.integer(ceiling(grid$n_col / kS))
  if (is.null(piece_um) || !is.finite(piece_um)) kP <- max(ntr, ntc) else {
    if (!is.numeric(piece_um) || length(piece_um) != 1L || piece_um <= 0 || abs(piece_um / tile_um - round(piece_um / tile_um)) > 1e-9)
      stop(sprintf("simulate_visium_hd_pieces(): `piece_um` deve essere un multiplo di %g um (tessere della casualita', DEC-074).", tile_um), call. = FALSE)
    kP <- as.integer(round(piece_um / tile_um))
  }
  npr <- as.integer(ceiling(ntr / kP)); npc <- as.integer(ceiling(ntc / kP))
  list(kS = kS, ntr = ntr, ntc = ntc, kP = kP, npr = npr, npc = npc, n_tiles = ntr * ntc, n_pieces = npr * npc, tile_um = tile_um,
       piece_um = if (is.null(piece_um) || !is.finite(piece_um)) Inf else piece_um)
}
.p27_bin_tile <- function(bin_id, grid, L) { b <- bin_id - 1L; ((b %/% grid$n_col) %/% L$kS) * L$ntc + (b %% grid$n_col) %/% L$kS + 1L }
.p27_tile_piece <- function(tile, L) { t <- tile - 1L; ((t %/% L$ntc) %/% L$kP) * L$npc + (t %% L$ntc) %/% L$kP + 1L }
# righe e colonne (0-based, estremi inclusi) dei quadratini di un pezzo
.p27_piece_bins <- function(k, grid, L) {
  pr <- (k - 1L) %/% L$npc; pc <- (k - 1L) %% L$npc; s <- L$kP * L$kS
  c(i0 = pr * s, i1 = min(grid$n_row, (pr + 1L) * s) - 1L, j0 = pc * s, j1 = min(grid$n_col, (pc + 1L) * s) - 1L)
}
.p27_ids <- function(rg, grid) { ii <- rep(rg[["i0"]]:rg[["i1"]], each = rg[["j1"]] - rg[["j0"]] + 1L); jj <- rep(rg[["j0"]]:rg[["j1"]], times = rg[["i1"]] - rg[["i0"]] + 1L)
  as.integer(ii * grid$n_col + jj + 1L) }

# --- verita' geometrica di un pezzo (stessa aritmetica di hd_tissue_area() e hd_overlap(), ristretta ai quadratini del pezzo) ---
.p27_tissue_area <- function(grid, region_polygons, rg, tile_um = 500, grid_um = 1e-6) {
  ids <- .p27_ids(rg, grid); out <- numeric(length(ids)); a_full <- grid$bin_um^2
  loc <- function(id) { i <- (id - 1L) %/% grid$n_col - rg[["i0"]]; j <- (id - 1L) %% grid$n_col - rg[["j0"]]; i * (rg[["j1"]] - rg[["j0"]] + 1L) + j + 1L }
  R <- geos::as_geos_geometry(sf::st_set_crs(region_polygons, NA)); bbR <- .hd_bboxes(region_polygons)
  tl <- .hd_tiles(grid, tile_um)
  for (r in seq_along(R)) {
    rr <- .hd_bbox_bins(grid, bbR[r, , drop = FALSE])
    for (t in seq_len(nrow(tl))) {
      i0 <- max(tl$i0[t], rr[1, "i0"], rg[["i0"]]); i1 <- min(tl$i1[t] - 1L, rr[1, "i1"], rg[["i1"]])
      j0 <- max(tl$j0[t], rr[1, "j0"], rg[["j0"]]); j1 <- min(tl$j1[t] - 1L, rr[1, "j1"], rg[["j1"]])
      if (i1 < i0 || j1 < j0) next
      ii <- rep(i0:i1, each = j1 - j0 + 1L); jj <- rep(j0:j1, times = i1 - i0 + 1L)
      sq <- .hd_squares(grid, ii, jj); id <- ii * grid$n_col + jj + 1L; li <- loc(id)
      inside <- geos::geos_prepared_contains(R[r], sq); touch <- geos::geos_prepared_intersects(R[r], sq)
      out[li[inside]] <- out[li[inside]] + a_full
      bd <- which(touch & !inside)
      if (length(bd)) {
        org <- c(grid$x0 + tl$j0[t] * grid$bin_um, grid$y0 + tl$i0[t] * grid$bin_um)
        out[li[bd]] <- out[li[bd]] + overlay_prec("intersection", .hd_sfc(sq[bd]), region_polygons[r], grid_um = grid_um, origin = org)
      }
    }
  }
  data.frame(bin_id = ids, tissue_area_um2 = out)
}

.p27_overlap <- function(grid, territories, nuclei, gidx, rg, tile_um = 500, grid_um = 1e-6) {
  n <- length(territories)
  empty <- data.frame(bin_id = integer(0), cell = integer(0), area_territory_um2 = numeric(0), area_nucleus_um2 = numeric(0))
  if (!n) return(empty)
  bb <- .hd_bbox_bins(grid, .hd_bboxes(territories))
  bb[, "i0"] <- pmax(bb[, "i0"], rg[["i0"]]); bb[, "i1"] <- pmin(bb[, "i1"], rg[["i1"]])
  bb[, "j0"] <- pmax(bb[, "j0"], rg[["j0"]]); bb[, "j1"] <- pmin(bb[, "j1"], rg[["j1"]])
  ni <- bb[, "i1"] - bb[, "i0"] + 1L; nj <- bb[, "j1"] - bb[, "j0"] + 1L; ok <- ni > 0 & nj > 0
  if (!any(ok)) return(empty)
  cand <- rep(seq_len(n)[ok], (ni * nj)[ok]); off <- sequence((ni * nj)[ok]) - 1L
  ci <- bb[cand, "i0"] + off %/% nj[cand]; cj <- bb[cand, "j0"] + off %% nj[cand]
  nt <- if (is.null(tile_um) || !is.finite(tile_um)) max(grid$n_row, grid$n_col) else max(1L, as.integer(round(tile_um / grid$bin_um)))
  tile <- (ci %/% nt) * (grid$n_col %/% nt + 1L) + cj %/% nt
  res <- list()
  for (t in sort(unique(tile))) {
    k <- which(tile == t); i <- ci[k]; j <- cj[k]; cc <- cand[k]
    org <- c(grid$x0 + (min(j) %/% nt) * nt * grid$bin_um, grid$y0 + (min(i) %/% nt) * nt * grid$bin_um)
    sq <- .hd_sfc(.hd_squares(grid, i, j))
    aT <- overlay_prec("intersection", sq, territories[cc], grid_um = grid_um, origin = org)
    keep <- aT > 0; aN <- numeric(sum(keep))
    if (any(keep)) aN <- overlay_prec("intersection", sq[keep], nuclei[cc[keep]], grid_um = grid_um, origin = org)
    res[[length(res) + 1L]] <- data.frame(bin_id = as.integer(i[keep] * grid$n_col + j[keep] + 1L), cell = as.integer(gidx[cc[keep]]),
                                          area_territory_um2 = aT[keep], area_nucleus_um2 = aN)
  }
  out <- if (length(res)) do.call(rbind, res) else empty
  out <- out[order(out$bin_id, out$cell), ]; rownames(out) <- NULL
  out
}

# --- campo di variabilita' fine su una finestra (rumore bianco per tessera della sottogriglia, stessa aritmetica di .s23_gaussian_z) ---
# finestra: righe P0...P0+nr-1 e colonne Q0...Q0+nc-1 della sottogriglia globale di passo `step` (0-based dall'origine della griglia)
.p27_field_window <- function(P0, nr, Q0, nc, s2, Lc, step, streams, L, tile_px) {
  if (s2 <= 0) return(matrix(1, nr, nc))
  g <- .s23_taps(Lc, step); h <- (length(g) - 1L) %/% 2L
  r0 <- P0 - h; c0 <- Q0 - h; NR <- nr + 2L * h; NC <- nc + 2L * h
  if (r0 < 0L || c0 < 0L || (r0 + NR) > L$ntr * tile_px || (c0 + NC) > L$ntc * tile_px)
    stop("simulate_visium_hd_pieces(): finestra del campo fuori dal reticolo delle tessere.", call. = FALSE)
  W <- matrix(0, NR, NC)
  for (tp in (r0 %/% tile_px):((r0 + NR - 1L) %/% tile_px)) for (tq in (c0 %/% tile_px):((c0 + NC - 1L) %/% tile_px)) {
    assign(".Random.seed", streams[[tp * L$ntc + tq + 1L]], envir = globalenv())
    z <- matrix(stats::rnorm(tile_px * tile_px), tile_px, tile_px)
    gr <- tp * tile_px + seq_len(tile_px) - 1L; gc <- tq * tile_px + seq_len(tile_px) - 1L
    kr <- which(gr >= r0 & gr < r0 + NR); kc <- which(gc >= c0 & gc < c0 + NC)
    W[gr[kr] - r0 + 1L, gc[kc] - c0 + 1L] <- z[kr, kc]
  }
  Tm <- matrix(0, nr, NC)
  for (a in seq_along(g)) Tm <- Tm + g[a] * W[a - 1L + seq_len(nr), , drop = FALSE]
  Z <- matrix(0, nr, nc)
  for (a in seq_along(g)) Z <- Z + g[a] * Tm[, a - 1L + seq_len(nc), drop = FALSE]
  Z <- Z / sum(g^2)
  tabF <- .s23_gamma_table(s2)
  Z <- pmin(pmax(Z, .s23_ztab[1]), .s23_ztab[length(.s23_ztab)])
  matrix(stats::approx(.s23_ztab, tabF, xout = as.vector(Z))$y, nr, nc)
}

# maschera "centro del sottopixel nei poligoni" su una finestra della sottogriglia globale (come .s23_raster_center_in, con i centri
# calcolati dall'origine della griglia, cosi' ogni sottopixel ha lo stesso valore in ogni finestra)
.p27_raster_center_in <- function(g, gx0, gy0, step, P0, nr, Q0, nc) {
  M <- matrix(FALSE, nr, nc)
  if (!length(g)) return(M)
  k <- wk::wk_coords(geos::geos_boundary(g))
  if (!nrow(k)) return(M)
  n <- nrow(k); same <- k$ring_id[-1] == k$ring_id[-n] & k$feature_id[-1] == k$feature_id[-n] & k$part_id[-1] == k$part_id[-n]
  x1 <- k$x[-n][same]; y1 <- k$y[-n][same]; x2 <- k$x[-1][same]; y2 <- k$y[-1][same]
  lo <- pmin(y1, y2); hi <- pmax(y1, y2)
  r_lo <- pmax(P0, as.integer(ceiling((lo - gy0) / step - 0.5))); r_hi <- pmin(P0 + nr - 1L, as.integer(ceiling((hi - gy0) / step - 0.5)) - 1L)
  keep <- r_hi >= r_lo & y1 != y2
  x1 <- x1[keep]; y1 <- y1[keep]; x2 <- x2[keep]; y2 <- y2[keep]; r_lo <- r_lo[keep]; r_hi <- r_hi[keep]
  if (!length(x1)) return(M)
  len <- r_hi - r_lo + 1L
  e <- rep(seq_along(len), len); rr <- r_lo[e] + sequence(len) - 1L
  yc <- gy0 + (rr + 0.5) * step
  xs <- x1[e] + (yc - y1[e]) * (x2[e] - x1[e]) / (y2[e] - y1[e])
  o <- order(rr, xs); rr <- rr[o]; xs <- xs[o]
  if (length(rr) %% 2L) stop(".p27_raster_center_in(): numero dispari di attraversamenti (poligono non chiuso?).", call. = FALSE)
  a <- seq(1L, length(rr), by = 2L)
  stopifnot(all(rr[a] == rr[a + 1L]))
  c_lo <- pmax(Q0, as.integer(ceiling((xs[a] - gx0) / step - 0.5))); c_hi <- pmin(Q0 + nc - 1L, as.integer(ceiling((xs[a + 1L] - gx0) / step - 0.5)) - 1L)
  ok <- c_hi >= c_lo; r_ <- rr[a][ok] - P0; c_lo <- c_lo[ok] - Q0; c_hi <- c_hi[ok] - Q0
  nb <- nr * (nc + 1L)
  D <- tabulate(c_lo * nr + r_ + 1L, nb) - tabulate((c_hi + 1L) * nr + r_ + 1L, nb)
  D <- matrix(D, nr, nc + 1L)
  S <- t(apply(D, 1, cumsum))
  S[, seq_len(nc), drop = FALSE] > 0L
}

# distanza dal bordo del tessuto su una finestra: tessuto = centro nelle regioni; fuori dal dominio (riquadro di S2.3) e, con frame_open,
# fuori dal riquadro dell'area, i sottopixel valgono tessuto (non sono sorgenti della distanza)
.p27_edt_window <- function(G, P0, nr, Q0, nc, env) {
  M <- .p27_raster_center_in(G, env$gx0, env$gy0, env$step, P0, nr, Q0, nc)
  pg <- P0 + seq_len(nr) - 1L; qg <- Q0 + seq_len(nc) - 1L
  out_dom <- outer(pg < env$dom[["p0"]] | pg >= env$dom[["p1"]], qg < env$dom[["q0"]] | qg >= env$dom[["q1"]], `|`)
  M <- M | out_dom
  if (env$frame_open) {
    xc <- env$gx0 + (qg + 0.5) * env$step; yc <- env$gy0 + (pg + 0.5) * env$step
    fb <- env$frame
    M <- M | outer(yc < fb[["ymin"]] | yc > fb[["ymax"]], xc < fb[["xmin"]] | xc > fb[["xmax"]], `|`)
  }
  .s23_edt(M, env$step)
}

# gruppi consecutivi uguali di piu' colonne intere ordinate: indice dell'ultimo elemento di ogni gruppo
.p27_runs <- function(...) {
  v <- list(...); n <- length(v[[1]]); if (!n) return(integer(0))
  d <- rep(FALSE, n - 1L); for (x in v) d <- d | (x[-1L] != x[-n])
  c(which(d), n)
}

# --- nascita, spostamento e conteggio dei trascritti delle cellule di una tessera (stessa logica di place_and_move_transcripts()) ---
.p27_place_tile <- function(idx, gi, Nc, Amn, Amc, Fmax, sN, sC, Fm, FP0, FQ0, Dm, DP0, DQ0, env, Xp, Xi, Xx) {
  p <- env$params; step <- env$step; gx0 <- env$gx0; gy0 <- env$gy0
  m <- length(idx); pos <- seq_len(m)
  fN <- match(idx, sN$cell); lN <- nrow(sN) - match(idx, rev(sN$cell)) + 1L
  fC <- match(idx, sC$cell); lC <- nrow(sC) - match(idx, rev(sC$cell)) + 1L
  cumN <- cumsum(sN$area); cumC <- cumsum(sC$area)
  stN <- ifelse(is.na(fN), NA, c(0, cumN)[fN]); stC <- ifelse(is.na(fC), NA, c(0, cumC)[fC])
  ttN <- ifelse(is.na(fN), 0, cumN[pmax(1L, lN)] - stN); ttC <- ifelse(is.na(fC), 0, cumC[pmax(1L, lC)] - stC)
  pn <- ifelse(Amn > 0, p$r * Amn / (p$r * Amn + Amc), 0)
  pn[is.na(fC) | ttC <= 0] <- 1; pn[is.na(fN) | ttN <= 0] <- 0
  dead <- which(Nc > 0L & ttN <= 0 & ttC <= 0)
  if (length(dead)) stop(sprintf("simulate_visium_hd_pieces(): %d cellule con trascritti e area nulla.", length(dead)), call. = FALSE)
  NX <- sN[, c("ax", "ay", "bx", "by", "cx", "cy")]; CX <- sC[, c("ax", "ay", "bx", "by", "cx", "cy")]
  nrF <- if (is.null(Fm)) 0L else nrow(Fm); ncF <- if (is.null(Fm)) 0L else ncol(Fm)
  need <- Nc; acc_cell <- list(); acc_x <- list(); acc_y <- list(); acc_nuc <- list(); it <- 0L
  while (any(need > 0L)) {
    it <- it + 1L; if (it > 200L) stop("simulate_visium_hd_pieces(): accettazione non converge.", call. = FALSE)
    cc <- pos[need > 0L]; mm <- as.integer(ceiling(need[need > 0L] * Fmax[cc] * 1.2)) + 2L
    cell <- rep(cc, mm); nm <- length(cell)
    isn <- stats::runif(nm) < pn[cell]
    u <- stats::runif(nm)
    ti <- integer(nm)
    if (any(isn)) { cn <- cell[isn]; tg <- stN[cn] + u[isn] * ttN[cn]; j <- findInterval(tg, cumN, left.open = TRUE) + 1L; ti[isn] <- pmin(pmax(j, fN[cn]), lN[cn]) }
    if (any(!isn)) { cn <- cell[!isn]; tg <- stC[cn] + u[!isn] * ttC[cn]; j <- findInterval(tg, cumC, left.open = TRUE) + 1L; ti[!isn] <- pmin(pmax(j, fC[cn]), lC[cn]) }
    r1 <- stats::runif(nm); r2 <- stats::runif(nm); fl <- r1 + r2 > 1; r1[fl] <- 1 - r1[fl]; r2[fl] <- 1 - r2[fl]
    x <- numeric(nm); y <- numeric(nm)
    if (any(isn)) { j <- ti[isn]; x[isn] <- NX$ax[j] + r1[isn] * (NX$bx[j] - NX$ax[j]) + r2[isn] * (NX$cx[j] - NX$ax[j]); y[isn] <- NX$ay[j] + r1[isn] * (NX$by[j] - NX$ay[j]) + r2[isn] * (NX$cy[j] - NX$ay[j]) }
    if (any(!isn)) { j <- ti[!isn]; x[!isn] <- CX$ax[j] + r1[!isn] * (CX$bx[j] - CX$ax[j]) + r2[!isn] * (CX$cx[j] - CX$ax[j]); y[!isn] <- CX$ay[j] + r1[!isn] * (CX$by[j] - CX$ay[j]) + r2[!isn] * (CX$cy[j] - CX$ay[j]) }
    ua <- stats::runif(nm)
    ok <- if (is.null(Fm)) rep(TRUE, nm) else {
      pa <- pmin(nrF, pmax(1L, as.integer(floor((y - gy0) / step)) - FP0 + 1L)); qa <- pmin(ncF, pmax(1L, as.integer(floor((x - gx0) / step)) - FQ0 + 1L))
      ua < Fm[cbind(pa, qa)] / Fmax[cell] }
    cell <- cell[ok]; x <- x[ok]; y <- y[ok]; isn <- isn[ok]
    rk <- sequence(rle(cell)$lengths)
    take <- rk <= need[cell]
    acc_cell[[it]] <- cell[take]; acc_x[[it]] <- x[take]; acc_y[[it]] <- y[take]; acc_nuc[[it]] <- isn[take]
    need <- need - tabulate(cell[take], m)
  }
  cell <- unlist(acc_cell); x <- unlist(acc_x); y <- unlist(acc_y); isn <- unlist(acc_nuc)
  o <- order(cell); cell <- cell[o]; x <- x[o]; y <- y[o]; isn <- isn[o]
  # distanza dal bordo alla nascita (finestra con l'alone) e spostamento
  pa <- pmin(nrow(Dm), pmax(1L, as.integer(floor((y - gy0) / step)) - DP0 + 1L)); qa <- pmin(ncol(Dm), pmax(1L, as.integer(floor((x - gx0) / step)) - DQ0 + 1L))
  d <- Dm[cbind(pa, qa)]
  nt <- length(cell); u1 <- stats::runif(nt); u2 <- stats::runif(nt)
  tail <- u1 < p$w; wide <- !tail & (u2 < exp(-d / p$D))
  dx <- numeric(nt); dy <- numeric(nt)
  nT <- sum(tail); if (nT) { rho <- stats::rgamma(nT, shape = 2, scale = p$ell); th <- stats::runif(nT, 0, 2 * pi); dx[tail] <- rho * cos(th); dy[tail] <- rho * sin(th) }
  nG <- sum(!tail); if (nG) { sg <- ifelse(wide[!tail], p$sigma_e, p$sigma_in); dx[!tail] <- stats::rnorm(nG) * sg; dy[!tail] <- stats::rnorm(nG) * sg }
  x1 <- x + dx; y1 <- y + dy
  if (env$w2 > 0) {
    sec <- stats::runif(nt) < env$w2; nS <- sum(sec)
    if (nS) {
      rho2 <- stats::rgamma(nS, shape = 2, scale = p$ell2); th2 <- stats::runif(nS, 0, 2 * pi)
      x1[sec] <- x[sec] + rho2 * cos(th2); y1[sec] <- y[sec] + rho2 * sin(th2)
      if (env$frame_open) { fb <- env$frame; x1[sec] <- .s24_fold(x1[sec], fb[["xmin"]], fb[["xmax"]]); y1[sec] <- .s24_fold(y1[sec], fb[["ymin"]], fb[["ymax"]]) }
    }
  }
  g <- env$grid
  ii <- floor((y1 - g$y0) / g$bin_um); jj <- floor((x1 - g$x0) / g$bin_um)
  inside <- ii >= 0 & ii < g$n_row & jj >= 0 & jj < g$n_col
  gene <- unlist(lapply(idx, function(cc) { a <- Xp[cc] + 1L; b <- Xp[cc + 1L]; if (b < a) integer(0) else rep(Xi[a:b] + 1L, Xx[a:b]) }))
  stopifnot(length(gene) == nt)
  lost <- tabulate(cell[!inside], m); lost_nuc <- tabulate(cell[!inside & isn], m)
  lost_gene <- tabulate(gene[!inside], env$nG)
  bi <- as.integer(ii[inside] * g$n_col + jj[inside] + 1); ci <- gi[cell[inside]]; ge <- gene[inside]; ni <- isn[inside]
  ko <- order(bi, ci, ge); bi <- bi[ko]; ci <- ci[ko]; ge <- ge[ko]; ni <- ni[ko]
  e3 <- .p27_runs(bi, ci, ge); e2 <- .p27_runs(bi, ci)
  tr <- data.frame(bin_id = bi[e3], cell = ci[e3], gene = ge[e3], n = as.integer(diff(c(0L, e3))))
  og <- data.frame(bin_id = bi[e2], cell = ci[e2], n = as.integer(diff(c(0L, e2))), n_born_nucleus = as.integer(diff(c(0L, cumsum(ni)[e2]))))
  list(tr = tr, og = og, lost = lost, lost_nuc = lost_nuc, lost_gene = lost_gene, it = it)
}

.p27_peak_mb <- function() { s <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0)); v <- s[grepl("^VmHWM", s)]
  if (length(v)) as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", v)) / 1024 else NA_real_ }
.p27_rss_mb <- function() { s <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0)); v <- s[grepl("^VmRSS", s)]
  if (length(v)) as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", v)) / 1024 else NA_real_ }

# esegue FUN(k) per i pezzi `ks` (un processo figlio nuovo per pezzo con workers > 1); FUN scrive il suo file e ritorna una riga di stato
.p27_lapply <- function(ks, FUN, workers) {
  out <- if (workers > 1L && length(ks) > 1L) parallel::mclapply(ks, FUN, mc.cores = workers, mc.preschedule = FALSE) else lapply(ks, FUN)
  bad <- vapply(out, function(o) inherits(o, "try-error") || is.null(o), logical(1))
  if (any(bad)) stop("simulate_visium_hd_pieces(): errore in un pezzo: ", paste(as.character(out[[which(bad)[1]]]), collapse = " "), call. = FALSE)
  do.call(rbind, out)
}
.p27_f <- function(env, ...) file.path(env$wd, paste0(...))
.p27_tab <- function(key, w, n) { v <- numeric(n); if (length(key)) { s <- rowsum(as.numeric(w), key); v[as.integer(rownames(s))] <- s[, 1] }; v }

# fase 0: area di tessuto e verita' geometrica dei quadratini del pezzo
.p27_stage0 <- function(k, env) {
  f <- .p27_f(env, "s0_", k, ".rds")
  if (file.exists(f)) return(data.frame(piece = k, stage = "0", elapsed_s = NA_real_, peak_mb = NA_real_, reused = TRUE))
  t0 <- proc.time()[["elapsed"]]
  inp <- readRDS(.p27_f(env, "in_", k, ".rds"))
  rg <- .p27_piece_bins(k, env$grid, env$L)
  ta <- .p27_tissue_area(env$grid, env$regions, rg)
  ov <- .p27_overlap(env$grid, inp$touch$terr, inp$touch$nuc, inp$touch$gi, rg)
  saveRDS(list(tissue = ta, overlap = ov), f)
  data.frame(piece = k, stage = "0", elapsed_s = proc.time()[["elapsed"]] - t0, peak_mb = .p27_peak_mb(), reused = FALSE)
}

# fase A: profili, nascita, spostamento delle cellule del pezzo; i trascritti si scrivono per il pezzo di arrivo
.p27_stageA <- function(k, env) {
  f <- .p27_f(env, "A_", k, ".rds")
  if (file.exists(f)) return(data.frame(piece = k, stage = "A", elapsed_s = NA_real_, peak_mb = NA_real_, reused = TRUE))
  t0 <- proc.time()[["elapsed"]]
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  inp <- readRDS(.p27_f(env, "in_", k, ".rds")); own <- inp$own; m <- length(own$gi); G <- env$nG
  res <- list(gene_gen = numeric(G), lost_gene = numeric(G), cells = data.frame(cell = integer(0), subtype = character(0), n_transcripts = integer(0),
              n_lost = integer(0), n_lost_born_nucleus = integer(0)), fallback = integer(0), sent = 0, counts = NULL, it_max = 0L)
  if (m) {
    ut <- unique(own$tile); chunks <- unname(split(seq_len(m), factor(own$tile, levels = ut)))
    prof <- simulate_cell_profiles(own$cell_df, env$profiles, seed = 1L, depth_factor = env$depth, chunks = chunks, streams = env$st_prof[ut])
    X <- prof$counts
    Tt <- geos::as_geos_geometry(own$terr); Nn <- geos::as_geos_geometry(own$nuc)
    Cy <- geos::geos_difference(Tt, Nn)
    A_n <- geos::geos_area(Nn); A_n[geos::geos_is_empty(Nn)] <- 0; A_c <- geos::geos_area(Cy); A_c[geos::geos_is_empty(Cy)] <- 0
    trN <- .s23_triangles(Nn); trC <- .s23_triangles(Cy)
    res$fallback <- own$gi[unique(c(attr(trN, "fallback"), attr(trC, "fallback")))]
    step <- env$step; gx0 <- env$gx0; gy0 <- env$gy0
    bbT <- sf::st_bbox(own$terr)
    Q0 <- as.integer(floor((bbT[["xmin"]] - gx0) / step)) - 1L; Q1 <- as.integer(ceiling((bbT[["xmax"]] - gx0) / step)) + 1L
    P0 <- as.integer(floor((bbT[["ymin"]] - gy0) / step)) - 1L; P1 <- as.integer(ceiling((bbT[["ymax"]] - gy0) / step)) + 1L
    nr <- P1 - P0; nc <- Q1 - Q0
    Fm <- if (env$field) .p27_field_window(P0, nr, Q0, nc, env$params$field_s2, env$params$field_L, step, env$st_field, env$L, env$tile_px) else NULL
    Fmax <- rep(1, m)
    if (!is.null(Fm)) {
      allt <- rbind(trN[, c("cell", "ax", "ay", "bx", "by", "cx", "cy")], trC[, c("cell", "ax", "ay", "bx", "by", "cx", "cy")])
      xmn <- tapply(pmin(allt$ax, allt$bx, allt$cx), allt$cell, min); xmx <- tapply(pmax(allt$ax, allt$bx, allt$cx), allt$cell, max)
      ymn <- tapply(pmin(allt$ay, allt$by, allt$cy), allt$cell, min); ymx <- tapply(pmax(allt$ay, allt$by, allt$cy), allt$cell, max)
      ids <- as.integer(names(xmn))
      a0 <- pmax(1L, as.integer(floor((ymn - gy0) / step)) - P0 + 1L); a1 <- pmin(nr, as.integer(floor((ymx - gy0) / step)) - P0 + 1L)
      b0 <- pmax(1L, as.integer(floor((xmn - gx0) / step)) - Q0 + 1L); b1 <- pmin(nc, as.integer(floor((xmx - gx0) / step)) - Q0 + 1L)
      for (h in seq_along(ids)) Fmax[ids[h]] <- max(Fm[a0[h]:a1[h], b0[h]:b1[h]])
    }
    Hp <- as.integer(ceiling(env$halo_um / step))
    Dm <- .p27_edt_window(env$regions_geos, P0 - Hp, nr + 2L * Hp, Q0 - Hp, nc + 2L * Hp, env)
    Nc <- as.integer(round(Matrix::colSums(X)))
    Xp <- X@p; Xi <- X@i; Xx <- as.integer(round(X@x))
    outs <- vector("list", length(chunks))
    # righe dei triangoli di ogni tessera (cellule contigue, triangoli ordinati per cellula): una sola ricerca vettoriale
    lo_c <- vapply(chunks, function(v) v[1L], 1L); hi_c <- vapply(chunks, function(v) v[length(v)], 1L)
    nN0 <- findInterval(lo_c - 0.5, trN$cell) + 1L; nN1 <- findInterval(hi_c + 0.5, trN$cell)
    nC0 <- findInterval(lo_c - 0.5, trC$cell) + 1L; nC1 <- findInterval(hi_c + 0.5, trC$cell)
    for (h in seq_along(chunks)) {
      idx <- chunks[[h]]
      sN <- trN[if (nN1[h] >= nN0[h]) nN0[h]:nN1[h] else integer(0), , drop = FALSE]; sC <- trC[if (nC1[h] >= nC0[h]) nC0[h]:nC1[h] else integer(0), , drop = FALSE]
      assign(".Random.seed", env$st_place[[ut[h]]], envir = globalenv())
      outs[[h]] <- .p27_place_tile(idx, own$gi[idx], Nc[idx], A_n[idx], A_c[idx], Fmax[idx], sN, sC, Fm, P0, Q0, Dm, P0 - Hp, Q0 - Hp, env, Xp, Xi, Xx)
    }
    tr <- do.call(rbind, lapply(outs, `[[`, "tr")); og <- do.call(rbind, lapply(outs, `[[`, "og"))
    res$lost_gene <- Reduce(`+`, lapply(outs, `[[`, "lost_gene"))
    res$gene_gen <- as.numeric(Matrix::rowSums(X))
    res$cells <- data.frame(cell = own$gi, subtype = prof$cells$subtype, n_transcripts = prof$cells$n_transcripts,
                            n_lost = as.integer(unlist(lapply(outs, `[[`, "lost"))), n_lost_born_nucleus = as.integer(unlist(lapply(outs, `[[`, "lost_nuc"))))
    res$it_max <- max(vapply(outs, `[[`, integer(1), "it"))
    if (env$keep_profiles) res$counts <- X
    # trascritti per pezzo di arrivo
    dt <- .p27_tile_piece(.p27_bin_tile(tr$bin_id, env$grid, env$L), env$L); do_ <- .p27_tile_piece(.p27_bin_tile(og$bin_id, env$grid, env$L), env$L)
    for (d in sort(unique(c(dt, do_)))) saveRDS(list(tr = tr[dt == d, , drop = FALSE], og = og[do_ == d, , drop = FALSE]), .p27_f(env, "AB_", k, "_", d, ".rds"))
    res$sent <- sum(tr$n[dt != k])
    res$n_arrived_rows <- nrow(tr)
  }
  res$elapsed_s <- proc.time()[["elapsed"]] - t0; res$peak_mb <- .p27_peak_mb()
  saveRDS(res, f)
  data.frame(piece = k, stage = "A", elapsed_s = res$elapsed_s, peak_mb = res$peak_mb, reused = FALSE)
}

# fase B: cattura dei trascritti arrivati nei quadratini del pezzo (suoi e ricevuti), fondo, somma a 8 um
.p27_stageB <- function(k, env) {
  f <- .p27_f(env, "B_", k, ".rds")
  if (file.exists(f)) return(data.frame(piece = k, stage = "B", elapsed_s = NA_real_, peak_mb = NA_real_, reused = TRUE))
  t0 <- proc.time()[["elapsed"]]
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  g <- env$grid; L <- env$L; G <- env$nG
  ff <- list.files(env$wd, pattern = paste0("^AB_[0-9]+_", k, "\\.rds$"), full.names = TRUE)
  parts <- lapply(ff, readRDS)
  tr <- do.call(rbind, c(list(data.frame(bin_id = integer(0), cell = integer(0), gene = integer(0), n = integer(0))), lapply(parts, `[[`, "tr")))
  og <- do.call(rbind, c(list(data.frame(bin_id = integer(0), cell = integer(0), n = integer(0), n_born_nucleus = integer(0))), lapply(parts, `[[`, "og")))
  rm(parts)
  rg <- .p27_piece_bins(k, g, L); ids <- .p27_ids(rg, g); bt_ids <- .p27_bin_tile(ids, g, L)
  btr <- .p27_bin_tile(tr$bin_id, g, L); bog <- .p27_bin_tile(og$bin_id, g, L)
  o <- order(btr, tr$bin_id, tr$cell, tr$gene); tr <- tr[o, , drop = FALSE]; btr <- btr[o]
  o <- order(bog, og$bin_id, og$cell); og <- og[o, , drop = FALSE]; bog <- bog[o]
  kept <- tr$n; kept_o <- og$n; kept_nuc <- og$n_born_nucleus
  ct <- env$c_t; Ms <- env$Mstar
  bins_by_tile <- split(ids, factor(bt_ids, levels = sort(unique(bt_ids))))
  bg_bin <- integer(0); bg_n <- integer(0); bg_gene <- integer(0); bg_genebin <- integer(0)
  tiles <- as.integer(names(bins_by_tile))
  # righe di ogni tessera nelle terne e nell'origine ordinate per tessera: una sola ricerca vettoriale
  t0r <- findInterval(tiles - 0.5, btr) + 1L; t1r <- findInterval(tiles + 0.5, btr)
  o0r <- findInterval(tiles - 0.5, bog) + 1L; o1r <- findInterval(tiles + 0.5, bog)
  for (h in seq_along(tiles)) {
    t <- tiles[h]
    if (env$thin) {
      it <- if (t1r[h] >= t0r[h]) t0r[h]:t1r[h] else integer(0); io <- if (o1r[h] >= o0r[h]) o0r[h]:o1r[h] else integer(0)
      if (length(it)) {
        assign(".Random.seed", env$st_cap[[t]], envir = globalenv())
        kt <- stats::rbinom(length(it), tr$n[it], ct[tr$bin_id[it]] / Ms)
        bt <- tr$bin_id[it]; cl <- tr$cell[it]; nt <- length(it)
        grp <- cumsum(c(TRUE, bt[-1L] != bt[-nt] | cl[-1L] != cl[-nt]))
        first <- which(!duplicated(grp))
        if (length(first) != length(io) || any(bt[first] != og$bin_id[io]) || any(cl[first] != og$cell[io]))
          stop("simulate_visium_hd_pieces(): terne e verita' di origine non allineate.", call. = FALSE)
        kk <- as.numeric(rowsum(kt, grp, reorder = FALSE))
        kn <- stats::rhyper(length(io), og$n_born_nucleus[io], og$n[io] - og$n_born_nucleus[io], kk)
        kept[it] <- kt; kept_o[io] <- kk; kept_nuc[io] <- kn
      }
    }
    if (env$background) {
      bb <- bins_by_tile[[as.character(t)]]
      assign(".Random.seed", env$st_bg[[t]], envir = globalenv())
      K <- stats::rpois(length(bb), env$b * ct[bb])
      if (sum(K) > 0) {
        gb <- sample.int(G, sum(K), replace = TRUE, prob = env$comp)
        w <- which(K > 0); bg_bin <- c(bg_bin, bb[w]); bg_n <- c(bg_n, as.integer(K[w]))
        bg_genebin <- c(bg_genebin, rep(bb, K)); bg_gene <- c(bg_gene, gb)
      }
    }
  }
  kt <- kept > 0
  ii <- c(tr$gene[kt], bg_gene); jj <- c(tr$bin_id[kt], bg_genebin); xx <- c(as.numeric(kept[kt]), rep(1, length(bg_gene)))
  # somma per (gene, quadratino) e per (gene, blocco di 8 um), ordinate
  o <- order(jj, ii); ii <- ii[o]; jj <- jj[o]; xx <- xx[o]; e <- .p27_runs(jj, ii)
  X2 <- data.frame(gene = ii[e], bin_id = jj[e], x = diff(c(0, cumsum(xx)[e])))
  nc8 <- g$n_col %/% g$block
  b8 <- as.integer(((X2$bin_id - 1L) %/% g$n_col) %/% g$block * nc8 + ((X2$bin_id - 1L) %% g$n_col) %/% g$block + 1L)
  o <- order(b8, X2$gene); b8s <- b8[o]; gs <- X2$gene[o]; xs <- X2$x[o]; e <- .p27_runs(b8s, gs)
  X8 <- data.frame(gene = gs[e], bin8_id = b8s[e], x = diff(c(0, cumsum(xs)[e])))
  ko <- kept_o > 0
  origin <- data.frame(bin_id = og$bin_id[ko], cell = og$cell[ko], n = as.integer(kept_o[ko]), n_born_nucleus = as.integer(kept_nuc[ko]))
  origin <- origin[order(origin$bin_id, origin$cell), , drop = FALSE]
  arr_bin <- if (nrow(tr)) { a <- rowsum(tr$n, tr$bin_id); data.frame(bin_id = as.integer(rownames(a)), n = as.integer(a[, 1])) } else data.frame(bin_id = integer(0), n = integer(0))
  cellsum <- if (nrow(tr)) { a <- rowsum(cbind(arr = as.numeric(tr$n), cap = as.numeric(kept)), tr$cell); data.frame(cell = as.integer(rownames(a)), arrived = a[, "arr"], captured = a[, "cap"]) } else
    data.frame(cell = integer(0), arrived = numeric(0), captured = numeric(0))
  gsum <- data.frame(arrived = .p27_tab(tr$gene, tr$n, G), captured = .p27_tab(tr$gene, kept, G), background = tabulate(bg_gene, G))
  trunc <- if (env$thin && nrow(arr_bin)) c(sum(arr_bin$n * env$c_raw[arr_bin$bin_id]), sum(arr_bin$n * (env$c_raw[arr_bin$bin_id] - ct[arr_bin$bin_id]))) else c(0, 0)
  u2 <- integer(length(ids)); if (nrow(X2)) { s <- rowsum(X2$x, X2$bin_id); u2[match(as.integer(rownames(s)), ids)] <- as.integer(s[, 1]) }
  res <- list(X2 = X2, X8 = X8, origin = origin, background = data.frame(bin_id = bg_bin, n = bg_n)[order(bg_bin), , drop = FALSE], arrived_bin = arr_bin,
              cells = cellsum, genes = gsum, trunc = trunc, ids = ids, u2 = u2, received = sum(as.numeric(tr$n)),
              elapsed_s = proc.time()[["elapsed"]] - t0, peak_mb = .p27_peak_mb())
  saveRDS(res, f)
  data.frame(piece = k, stage = "B", elapsed_s = res$elapsed_s, peak_mb = res$peak_mb, reused = FALSE)
}

#' Simulatore di Visium HD a pezzi (Step 2, S2.7)
#'
#' Stesso modello di `simulate_visium_hd()` (S2.5), calcolato a pezzi quadrati di lato `piece_um`: ogni cellula appartiene al pezzo che
#' contiene il suo centro e si simula una volta sola; i trascritti che arrivano nei quadratini di un altro pezzo passano a quel pezzo;
#' strisce, normalizzazione e composizione del fondo sono dell'area intera; la seconda coda si ripiega sul riquadro dell'area (`frame_open`)
#' o si perde fuori dalla griglia. Tutte le estrazioni casuali sono legate a tessere fisse di 56 um (DEC-074): il risultato non dipende da
#' `piece_um` ne' da `workers`. Sovra-generazione con il tetto fisso del tessuto (`stripe_ceiling()`, DEC-075): i fattori di striscia sopra il
#' tetto sono tagliati al tetto. Differenze dichiarate rispetto a `simulate_visium_hd()`: realizzazione casuale diversa (stessa legge salvo il
#' tetto); distanza dal bordo calcolata su finestre con un alone di 25 D (stessa scelta del modo); riquadro della riflessione dato da `frame`.
#' @param cell_layer uscita di `place_cells()`, oppure il percorso di un file RDS che la contiene (letto e liberato dopo la divisione: conviene
#'   per le aree grandi, cosi' i processi dei pezzi non portano con se' il livello intero).
#' @param seed intero.
#' @param piece_um lato dei pezzi (um), multiplo di 56; `Inf` = un pezzo solo. Predefinito 504 (DEC-076: memoria di un pezzo nel bilancio di 8 GB in tutti i tessuti).
#' @param archetype,layers,profiles,depth_factor come `simulate_visium_hd()`.
#' @param frame_open `TRUE` per aree ritagliate da un tessuto che continua: i lati del riquadro `frame` non sono bordo del tessuto.
#' @param frame riquadro dell'area (vettore `xmin`, `ymin`, `xmax`, `ymax`); predefinito: riquadro di territori e regioni, come
#'   `simulate_visium_hd()`. Con `frame_open = TRUE` vale due volte: la seconda coda si ripiega sui suoi lati, e fuori dal riquadro i
#'   sottopixel contano come tessuto nella distanza dal bordo (un lato del riquadro non e' un bordo del tessuto, quindi i trascritti nati
#'   vicino a un lato non ricevono la sfocatura larga della zona di bordo). Per un'area ritagliata da una sezione vera conviene dare la
#'   finestra del ritaglio: dove il tessuto non arriva al lato della finestra, il predefinito tratta come aperto un bordo vero (BL-237).
#'   Con `frame_open = FALSE` non ha effetto.
#' @param ceiling_q probabilita' di superare il tetto della sovra-generazione (DEC-075).
#' @param out_dir cartella dei file per pezzo (rilanciabile: i file gia' scritti si riusano); predefinito: cartella temporanea, cancellata.
#' @param workers processi (un processo figlio per pezzo e per fase).
#' @param assemble `"full"` (uscite come `simulate_visium_hd()`) o `"light"` (UMI per quadratino, matrice a 8 um, tabelle; per le sezioni).
#' @param keep_profiles `TRUE` mette nelle uscite la matrice dei trascritti generati per cellula (`truth$cell_counts`).
#' @param verbose messaggi di avanzamento.
#' @return con `assemble = "full"`, lista come `simulate_visium_hd()`; con `"light"`, lista `umi_2um`, `counts_8um`, `bins`, `bins_8um`,
#'   `truth` (`cells`, `genes`, `layers`), `metadata`. In `metadata`: pezzi, tetto, quota troncata, ripieghi della triangolazione, tempi e memoria per fase e pezzo, controlli.
#'   Strisce e sovra-generazione (DEC-075): `bins$c` e' il fattore di striscia usato, tagliato al tetto M*; `bins$c_raw` quello grezzo;
#'   `truth$layers$stripes$M` = `metadata$overgeneration` = M* (sovra-generazione usata, come in `simulate_visium_hd()`), `$M_raw` e `$C`
#'   il massimo grezzo dei fattori dell'area. Il fondo e' Poisson(b c) per quadratino con il fattore tagliato `c`, come la cattura.
#'   Ripresa (`out_dir`): l'impronta `run.rds` contiene argomenti, impronta del codice (`.geo_code_hash()`), dei profili e del livello
#'   cellulare (identificativi, tipi, centri, territori, nuclei, regioni); se la cartella
#'   contiene un'esecuzione diversa la funzione si ferma e dice quali campi cambiano (BL-247).
#' @export
simulate_visium_hd_pieces <- function(cell_layer, seed, piece_um = 504, archetype = NULL,
                                      layers = list(stripes = TRUE, background = TRUE, field = TRUE, second_tail = TRUE),
                                      frame_open = FALSE, frame = NULL, profiles = NULL, depth_factor = NULL, ceiling_q = 1e-5,
                                      out_dir = NULL, workers = 1L, assemble = c("full", "light"), keep_profiles = NULL, verbose = FALSE) {
  t0 <- proc.time()[["elapsed"]]; tm <- c()
  assemble <- match.arg(assemble)
  if (is.null(keep_profiles)) keep_profiles <- assemble == "full"
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }, add = TRUE)
  from_file <- is.character(cell_layer)
  if (from_file) cell_layer <- readRDS(cell_layer)
  px <- cell_layer$metadata$pixel_size_um
  if (is.null(px) || !is.numeric(px) || length(px) != 1L || !is.finite(px) || px <= 0)
    stop("simulate_visium_hd_pieces(): `cell_layer$metadata$pixel_size_um` manca o non e' valido (C-A): nessun default.", call. = FALSE)
  if (missing(seed) || length(seed) != 1L || !is.finite(seed)) stop("simulate_visium_hd_pieces(): serve un `seed`.", call. = FALSE)
  ly <- list(stripes = TRUE, background = TRUE, field = TRUE, second_tail = TRUE)
  if (!is.list(layers) || length(setdiff(names(layers), names(ly))) || any(!vapply(layers, function(v) is.logical(v) && length(v) == 1L && !is.na(v), logical(1))))
    stop("simulate_visium_hd_pieces(): `layers` accetta solo stripes, background, field, second_tail (TRUE o FALSE).", call. = FALSE)
  ly[names(layers)] <- layers
  if (is.null(archetype)) archetype <- visium_hd_archetype(cell_layer)
  if (!archetype %in% paste0("A", 1:6)) stop("simulate_visium_hd_pieces(): archetipo sconosciuto: ", archetype, call. = FALSE)
  sd_ <- .p27_draw_seeds(seed)
  cd <- cell_layer$cell_df; n <- nrow(cd)
  Tt <- sf::st_set_crs(cell_layer$cell_territories, NA); Nn <- sf::st_set_crs(cell_layer$cell_nuclei, NA); R <- sf::st_set_crs(cell_layer$region_polygons, NA)
  if (length(Tt) != n || length(Nn) != n) stop("simulate_visium_hd_pieces(): territori, nuclei e cell_df di lunghezza diversa.", call. = FALSE)
  if (any(geos::geos_is_empty(geos::as_geos_geometry(Tt)))) stop("simulate_visium_hd_pieces(): territori vuoti.", call. = FALSE)
  bb <- sf::st_bbox(c(Tt, R))
  g <- hd_grid_geometry(bb, margin_um = 280, bin_um = 2)
  L <- .p27_layout(g, piece_um)
  fb <- if (is.null(frame)) c(xmin = bb[["xmin"]], ymin = bb[["ymin"]], xmax = bb[["xmax"]], ymax = bb[["ymax"]]) else {
    fr <- frame[c("xmin", "ymin", "xmax", "ymax")]
    if (anyNA(fr) || any(!is.finite(fr)) || fr[["xmax"]] <= fr[["xmin"]] || fr[["ymax"]] <= fr[["ymin"]] ||
        fr[["xmin"]] < g$x0 || fr[["ymin"]] < g$y0 || fr[["xmax"]] > g$x0 + g$n_col * g$bin_um || fr[["ymax"]] > g$y0 + g$n_row * g$bin_um)
      stop("simulate_visium_hd_pieces(): `frame` deve essere un riquadro xmin, ymin, xmax, ymax dentro la griglia.", call. = FALSE)
    stats::setNames(as.numeric(fr), c("xmin", "ymin", "xmax", "ymax"))
  }
  cp <- capture_params(archetype); p <- transcript_params(archetype)
  p$w2 <- if (ly$second_tail) cp$w2 else 0; p$ell2 <- cp$ell2
  .s23_check_params(p, ly$field)
  Ms <- if (ly$stripes) stripe_ceiling(cp, ceiling_q) else 1
  prof_default <- is.null(profiles)
  if (prof_default) profiles <- load_cell_profiles(archetype)
  df0 <- if (!is.null(depth_factor)) depth_factor else if (prof_default) visium_hd_depth_factor(archetype) else profiles$depth_factor
  step <- 0.5; tile_px <- as.integer(round(L$tile_um / step))
  # cellule -> tessere -> pezzi; pezzi toccati da ogni territorio (verita' geometrica)
  ti <- as.integer(floor((cd$y - g$y0) / L$tile_um)); tj <- as.integer(floor((cd$x - g$x0) / L$tile_um))
  if (any(ti < 0L | ti >= L$ntr | tj < 0L | tj >= L$ntc)) stop("simulate_visium_hd_pieces(): centri delle cellule fuori dalla griglia.", call. = FALSE)
  ctile <- ti * L$ntc + tj + 1L; cpiece <- .p27_tile_piece(ctile, L)
  tb <- .hd_bbox_bins(g, .hd_bboxes(Tt)); s <- L$kP * L$kS
  pr0 <- tb[, "i0"] %/% s; pr1 <- tb[, "i1"] %/% s; pc0 <- tb[, "j0"] %/% s; pc1 <- tb[, "j1"] %/% s
  nr_ <- pr1 - pr0 + 1L; nc_ <- pc1 - pc0 + 1L
  tc <- rep(seq_len(n), nr_ * nc_); off <- sequence(nr_ * nc_) - 1L
  tp <- (pr0[tc] + off %/% nc_[tc]) * L$npc + pc0[tc] + off %% nc_[tc] + 1L
  # cartella di lavoro e impronta dell'esecuzione
  tmp_wd <- is.null(out_dir)
  wd <- if (tmp_wd) tempfile("s27_pieces_") else out_dir
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  if (tmp_wd) on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  fp <- list(seed = seed, piece_um = L$piece_um, archetype = archetype, layers = ly, frame_open = frame_open, frame = fb, ceiling_q = ceiling_q,
             Mstar = Ms, depth = df0 * Ms, n_cells = n, grid = g[c("x0", "y0", "n_row", "n_col")], keep_profiles = keep_profiles,
             cell_xy = sum(cd$x * seq_len(n)) + sum(cd$y), step1_seed = cell_layer$metadata$random_seed,
             code = .geo_code_hash("simulate_visium_hd_pieces"), profiles = .geo_object_hash(profiles),
             layer = .geo_object_hash(list(cd[, intersect(c("cell_id", "cell_type", "x", "y"), names(cd))], unclass(sf::st_as_binary(Tt)),
                                           unclass(sf::st_as_binary(Nn)), unclass(sf::st_as_binary(R)))))
  fpf <- file.path(wd, "run.rds")
  if (file.exists(fpf)) {
    fo <- readRDS(fpf)
    dif <- union(names(fp)[!vapply(names(fp), function(k) isTRUE(all.equal(fo[[k]], fp[[k]], tolerance = 0)), logical(1))], setdiff(names(fo), names(fp)))
    if (length(dif)) stop("simulate_visium_hd_pieces(): `out_dir` contiene un'esecuzione diversa (campi: ", paste(dif, collapse = ", "),
                          "); usa un'altra cartella o cancellala.", call. = FALSE)
  } else saveRDS(fp, fpf)
  t1 <- proc.time()[["elapsed"]]
  ord <- order(ctile, seq_len(n))
  for (k in seq_len(L$n_pieces)) {
    f <- file.path(wd, paste0("in_", k, ".rds")); if (file.exists(f)) next
    gi <- ord[cpiece[ord] == k]; tg <- sort(unique(tc[tp == k]))
    saveRDS(list(own = list(gi = gi, tile = ctile[gi], cell_df = cd[gi, c("cell_id", "cell_type"), drop = FALSE], terr = Tt[gi], nuc = Nn[gi]),
                 touch = list(gi = tg, terr = Tt[tg], nuc = Nn[tg])), f)
  }
  cell_info <- cd[, c("cell_id", "cell_type")]; step1_meta <- cell_layer$metadata[c("tissue_preset_used", "random_seed", "package_version")]
  rm(Tt, Nn, tc, tp, ord); if (from_file) rm(cell_layer)
  invisible(gc(FALSE))
  tm["split"] <- proc.time()[["elapsed"]] - t1
  if (verbose) message(sprintf("simulate_visium_hd_pieces(): %s, %d cellule, %d pezzi di %s um, tetto %.4f", archetype, n, L$n_pieces, format(L$piece_um), Ms))
  # flussi delle tessere
  t1 <- proc.time()[["elapsed"]]
  env <- new.env(parent = emptyenv())
  env$wd <- wd; env$grid <- g; env$L <- L; env$step <- step; env$gx0 <- g$x0; env$gy0 <- g$y0; env$tile_px <- tile_px
  env$regions <- R; env$regions_geos <- geos::as_geos_geometry(R)
  env$dom <- c(q0 = as.integer(floor((bb[["xmin"]] - g$x0) / step)) - 1L, q1 = as.integer(ceiling((bb[["xmax"]] - g$x0) / step)) + 1L,
               p0 = as.integer(floor((bb[["ymin"]] - g$y0) / step)) - 1L, p1 = as.integer(ceiling((bb[["ymax"]] - g$y0) / step)) + 1L)
  env$frame_open <- frame_open; env$frame <- fb; env$params <- p; env$w2 <- p$w2; env$field <- ly$field
  env$halo_um <- 25 * p$D; env$profiles <- profiles; env$depth <- df0 * Ms; env$nG <- ncol(profiles$P); env$keep_profiles <- keep_profiles
  env$st_prof <- .p27_streams(sd_[["profiles"]], L$n_tiles); env$st_place <- .p27_streams(sd_[["place"]], L$n_tiles)
  env$st_field <- if (ly$field) .p27_streams(sd_[["field"]], L$n_tiles) else NULL
  env$st_cap <- .p27_streams(sd_[["capture"]], L$n_tiles); env$st_bg <- .p27_streams(sd_[["background"]], L$n_tiles)
  tm["streams"] <- proc.time()[["elapsed"]] - t1
  rss_fork <- .p27_rss_mb()
  # fase 0
  t1 <- proc.time()[["elapsed"]]
  st0 <- .p27_lapply(seq_len(L$n_pieces), function(k) .p27_stage0(k, env), workers)
  tm["stage0"] <- proc.time()[["elapsed"]] - t1
  # tessuto e strisce dell'area intera
  t1 <- proc.time()[["elapsed"]]
  nb <- g$n_row * g$n_col; tarea <- numeric(nb)
  for (k in seq_len(L$n_pieces)) { s0 <- readRDS(file.path(wd, paste0("s0_", k, ".rds"))); tarea[s0$tissue$bin_id] <- s0$tissue$tissue_area_um2 }
  in_tissue <- tarea / g$bin_um^2 >= 0.5
  St <- if (ly$stripes) simulate_stripes(list(grid = g, bins = list(in_tissue = in_tissue)), cp, seed = sd_[["stripes"]]) else NULL
  env$c_raw <- if (is.null(St)) rep(1, nb) else St$c
  env$c_t <- pmin(env$c_raw, Ms); env$Mstar <- Ms; env$thin <- !is.null(St)
  env$background <- ly$background; env$b <- cp$b
  n_over <- sum(env$c_raw > Ms)
  tm["stripes"] <- proc.time()[["elapsed"]] - t1
  # fase A
  t1 <- proc.time()[["elapsed"]]
  stA <- .p27_lapply(seq_len(L$n_pieces), function(k) .p27_stageA(k, env), workers)
  tm["stageA"] <- proc.time()[["elapsed"]] - t1
  A <- lapply(seq_len(L$n_pieces), function(k) readRDS(file.path(wd, paste0("A_", k, ".rds"))))
  comp <- Reduce(`+`, lapply(A, `[[`, "gene_gen"))
  env$comp <- comp / sum(comp)
  # fase B
  t1 <- proc.time()[["elapsed"]]
  stB <- .p27_lapply(seq_len(L$n_pieces), function(k) .p27_stageB(k, env), workers)
  tm["stageB"] <- proc.time()[["elapsed"]] - t1
  # cucitura
  t1 <- proc.time()[["elapsed"]]
  out <- .p27_assemble(env, A, cd = cell_info, assemble = assemble, tarea = tarea, in_tissue = in_tissue, St = St)
  tm["assemble"] <- proc.time()[["elapsed"]] - t1
  stats_ <- rbind(st0, stA, stB)
  lay <- list(stripes = if (is.null(St)) NULL else c(St[c("alpha", "beta", "slow_alpha", "slow_beta", "z_row", "z_col", "norm", "C")],
                                                   list(M = Ms, M_raw = St$M), St[c("seed", "params")]),
              background = list(on = ly$background, b = if (ly$background) cp$b else 0), second_tail = p[c("w2", "ell2")], overgeneration = Ms,
              diffusion = p[c("r", "sigma_in", "sigma_e", "D", "w", "ell", "w2", "ell2")], field = list(on = ly$field, s2 = p$field_s2, L = p$field_L),
              depth = list(depth_factor = df0, overgeneration = Ms, profiles_winner = profiles$meta$winner, noise = profiles$noise))
  out$truth$layers <- lay
  out$metadata <- c(list(step = "S2.7", archetype = archetype, seed = seed, component_seeds = sd_, layers = ly, frame_open = frame_open, frame = fb,
                         piece_um = L$piece_um, tile_um = L$tile_um, n_pieces = L$n_pieces, halo_um = env$halo_um, overgeneration = Ms, ceiling_q = ceiling_q,
                         n_bins_over_ceiling = n_over, pixel_size_um = px, grid = g, n_cells = n, step1 = step1_meta, piece_stats = stats_,
                         main_rss_at_fork_mb = rss_fork, main_peak_mb = .p27_peak_mb(),
                         package_version = .geo_version(), code_hash = fp$code, profiles_hash = fp$profiles, elapsed_s = c(tm, total = proc.time()[["elapsed"]] - t0)),
                    out$metadata)
  out
}

# cucitura delle uscite dei pezzi (ordine globale: indipendente dalla divisione)
.p27_assemble <- function(env, A, cd, assemble, tarea, in_tissue, St) {
  g <- env$grid; L <- env$L; wd <- env$wd; nb <- g$n_row * g$n_col; G <- env$nG; Ms <- env$Mstar; n <- nrow(cd)
  B <- lapply(seq_len(L$n_pieces), function(k) readRDS(file.path(wd, paste0("B_", k, ".rds"))))
  gene_ids <- env$profiles$genes$gene_id
  bins <- hd_bin_table(g); bins$tissue_area_um2 <- tarea; bins$tissue_frac <- tarea / g$bin_um^2; bins$in_tissue <- in_tissue
  bins$c <- env$c_t; bins$c_raw <- env$c_raw
  b8 <- sum_to_8um(data.frame(bin_id = bins$bin_id, tissue_area_um2 = bins$tissue_area_um2, n_in_tissue = as.numeric(bins$in_tissue)), g,
                   cols = c("tissue_area_um2", "n_in_tissue"))
  nc8 <- g$n_col %/% g$block
  b8$i8 <- (b8$bin8_id - 1L) %/% nc8; b8$j8 <- (b8$bin8_id - 1L) %% nc8
  b8$x <- g$x0 + (b8$j8 + 0.5) * g$bin_um * g$block; b8$y <- g$y0 + (b8$i8 + 0.5) * g$bin_um * g$block
  b8$in_tissue <- b8$n_in_tissue > 0
  b8 <- b8[, c("bin8_id", "i8", "j8", "x", "y", "tissue_area_um2", "n_in_tissue", "in_tissue")]
  # cellule e geni
  cA <- do.call(rbind, lapply(A, `[[`, "cells")); cB <- do.call(rbind, lapply(B, `[[`, "cells"))
  arr <- .p27_tab(cB$cell, cB$arrived, n); cap <- .p27_tab(cB$cell, cB$captured, n)
  ca <- cA[match(seq_len(n), cA$cell), , drop = FALSE]
  cells <- data.frame(cell_id = cd$cell_id, cell_type = cd$cell_type, subtype = ca$subtype, n_transcripts_generated = ca$n_transcripts,
                      n_generated = arr + ca$n_lost, n_lost = ca$n_lost, n_thinned = arr - cap, n_captured = cap, stringsAsFactors = FALSE)
  cells$n_expected_captured <- cells$n_transcripts_generated / Ms
  lostg <- Reduce(`+`, lapply(A, function(a) as.numeric(a$lost_gene)))
  gB <- Reduce(function(a, b) Map(`+`, a, b), lapply(B, function(b) lapply(as.list(b$genes), as.numeric)))
  genes <- data.frame(gene = gene_ids, n_generated = gB$arrived + lostg, n_lost = lostg, n_thinned = gB$arrived - gB$captured,
                      n_captured = gB$captured, n_background = gB$background)
  ids_all <- unlist(lapply(B, `[[`, "ids"))
  u2 <- integer(nb); u2[ids_all] <- unlist(lapply(B, `[[`, "u2"))
  bg <- do.call(rbind, lapply(B, `[[`, "background")); bg <- bg[order(bg$bin_id), , drop = FALSE]; rownames(bg) <- NULL
  trunc <- Reduce(`+`, lapply(B, `[[`, "trunc"))
  sent <- sum(vapply(A, function(a) as.numeric(a$sent), numeric(1)))
  arrived_rows <- sum(vapply(A, function(a) if (is.null(a$n_arrived_rows)) 0 else as.numeric(a$n_arrived_rows), numeric(1)))
  received <- sum(vapply(B, `[[`, numeric(1), "received"))
  fallback <- sort(unlist(lapply(A, `[[`, "fallback")))
  X8t <- do.call(rbind, lapply(B, `[[`, "X8"))
  X8 <- Matrix::sparseMatrix(i = X8t$gene, j = X8t$bin8_id, x = as.numeric(X8t$x), dims = c(G, nrow(b8)), dimnames = list(gene_ids, NULL))
  checks <- list(
    cells_once = identical(sort(cA$cell), seq_len(n)),
    bins_once = identical(sort(ids_all), seq_len(nb)),
    conservation_cell = isTRUE(all(cells$n_transcripts_generated == cells$n_generated)),
    conservation_total = sum(cells$n_generated) == sum(cells$n_captured) + sum(cells$n_thinned) + sum(cells$n_lost),
    sent_received = isTRUE(all.equal(sum(arr), received, tolerance = 0)),
    umi_total = sum(as.numeric(u2)) == sum(cells$n_captured) + sum(as.numeric(bg$n)),
    genes_cells = sum(genes$n_captured) == sum(cells$n_captured) && sum(genes$n_generated) == sum(cells$n_generated),
    no_na = !anyNA(cells) && !anyNA(u2))
  meta <- list(n_generated = sum(cells$n_generated), n_captured = sum(cells$n_captured), n_background = sum(as.numeric(bg$n)), n_lost = sum(cells$n_lost),
               n_sent_between_pieces = sent, truncated_share = if (trunc[1] > 0) trunc[2] / trunc[1] else 0,
               fallback_cell_ids = cd$cell_id[fallback], n_fallback = length(fallback), checks = checks)
  if (assemble == "light") {
    meta$checks$sum_8um <- sum(X8t$x) == sum(as.numeric(u2))
    return(list(umi_2um = u2, counts_8um = X8, bins = bins, bins_8um = b8, truth = list(cells = cells, genes = genes, background = bg), metadata = meta))
  }
  X2t <- do.call(rbind, lapply(B, `[[`, "X2")); rm(B)
  X2 <- Matrix::sparseMatrix(i = X2t$gene, j = X2t$bin_id, x = as.numeric(X2t$x), dims = c(G, nb), dimnames = list(gene_ids, NULL)); rm(X2t)
  X8s <- sum_to_8um(X2, g); dimnames(X8s) <- list(gene_ids, NULL)
  meta$checks$sum_8um <- Matrix::nnzero(X8s - X8) == 0L
  S0 <- lapply(seq_len(L$n_pieces), function(k) readRDS(file.path(wd, paste0("s0_", k, ".rds")))$overlap)
  ov <- do.call(rbind, S0); ov <- ov[order(ov$bin_id, ov$cell), , drop = FALSE]; rownames(ov) <- NULL
  ov <- data.frame(bin_id = ov$bin_id, cell_id = cd$cell_id[ov$cell], area_territory_um2 = ov$area_territory_um2, area_nucleus_um2 = ov$area_nucleus_um2)
  ov$frac_territory <- ov$area_territory_um2 / g$bin_um^2; ov$frac_nucleus <- ov$area_nucleus_um2 / g$bin_um^2
  ov8 <- .s25_to_8um(ov[, c("bin_id", "cell_id", "area_territory_um2", "area_nucleus_um2")], g, c("area_territory_um2", "area_nucleus_um2"))
  ov8$frac_territory <- ov8$area_territory_um2 / (g$bin_um * g$block)^2; ov8$frac_nucleus <- ov8$area_nucleus_um2 / (g$bin_um * g$block)^2
  B2 <- lapply(seq_len(L$n_pieces), function(k) readRDS(file.path(wd, paste0("B_", k, ".rds")))[c("origin", "arrived_bin")])
  og <- do.call(rbind, lapply(B2, `[[`, "origin")); og <- og[order(og$bin_id, og$cell), , drop = FALSE]; rownames(og) <- NULL
  og <- data.frame(bin_id = og$bin_id, cell_id = cd$cell_id[og$cell], n = og$n, n_born_nucleus = og$n_born_nucleus)
  og8 <- .s25_to_8um(og, g, c("n", "n_born_nucleus"))
  ab <- do.call(rbind, lapply(B2, `[[`, "arrived_bin")); arr2 <- integer(nb); arr2[ab$bin_id] <- ab$n
  bg8 <- if (nrow(bg)) { s <- sum_to_8um(bg, g, cols = "n"); s$n <- as.integer(s$n); s } else data.frame(bin8_id = integer(0), n = integer(0))
  cc <- NULL
  if (env$keep_profiles) {
    cc <- do.call(cbind, lapply(A, `[[`, "counts"))
    if (!is.null(cc)) cc <- cc[, match(as.character(cd$cell_id), colnames(cc)), drop = FALSE]
  }
  list(counts_2um = X2, counts_8um = X8s, bins = bins, bins_8um = b8,
       truth = list(overlap = ov, overlap_8um = ov8, origin = og, origin_8um = og8, background = bg, background_8um = bg8, arrived_2um = arr2,
                    cell_counts = cc, cells = cells, genes = genes),
       metadata = meta)
}
