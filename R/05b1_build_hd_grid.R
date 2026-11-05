# R/05b1_build_hd_grid.R - Step 2, S2.1: griglia di Visium HD (quadratini di 2 um, blocchi di 8 um) e verita' geometrica.
# Design dello Step 2 v0.2 (docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md), sezioni 1, 3, 8, 12, 13;
# definizioni operative e controlli in docs/steps/S2.1.md. Sostituisce R/05_grid_sampling.R (BL-005, BL-010).
# Ogni sovrapposizione passa da overlay_prec() (R/04b5, arrotondamento a 1e-6 um, D2-12), in coordinate locali della tessera (BL-104).
# Dipendenze: sf, geos, Matrix.

#' Geometria della griglia di Visium HD
#'
#' Quadrati di `bin_um` senza spazi; i = riga (asse y), j = colonna (asse x), numerate da 0 dall'origine (x0, y0) = angolo minimo di
#' `bbox` meno `margin_um`. Righe e colonne arrotondate per eccesso a un multiplo di 4, cosi' ogni blocco di 8 um ha i suoi 16 figli.
#' @param bbox vettore `xmin`, `ymin`, `xmax`, `ymax` (um).
#' @param margin_um margine di vetro (um), multiplo di 8 (design sezione 3, DEC-038 punto 3).
#' @param bin_um lato del quadratino (um).
#' @param block numero di quadratini per lato di un blocco di 8 um.
#' @return lista `x0`, `y0`, `bin_um`, `n_row`, `n_col`, `block`, `margin_um`.
#' @export
hd_grid_geometry <- function(bbox, margin_um = 280, bin_um = 2, block = 4L) {
  bb <- as.numeric(bbox[c("xmin", "ymin", "xmax", "ymax")])
  if (any(!is.finite(bb)) || bb[3] < bb[1] || bb[4] < bb[2]) stop("hd_grid_geometry(): riquadro non valido.", call. = FALSE)
  if (!is.numeric(bin_um) || length(bin_um) != 1L || !is.finite(bin_um) || bin_um <= 0)
    stop("hd_grid_geometry(): `bin_um` deve essere un numero > 0.", call. = FALSE)
  step8 <- bin_um * block
  if (!is.numeric(margin_um) || length(margin_um) != 1L || !is.finite(margin_um) || margin_um < 0 ||
      abs(margin_um / step8 - round(margin_um / step8)) > 1e-9)
    stop(sprintf("hd_grid_geometry(): `margin_um` deve essere un multiplo di %g um (somma a 8 um completa).", step8), call. = FALSE)
  x0 <- bb[1] - margin_um; y0 <- bb[2] - margin_um
  nb <- function(lo, hi) { n <- ceiling((hi + margin_um - lo) / bin_um - 1e-9); max(block, block * ceiling(n / block)) }
  list(x0 = x0, y0 = y0, bin_um = bin_um, n_row = as.integer(nb(y0, bb[4])), n_col = as.integer(nb(x0, bb[3])),
       block = as.integer(block), margin_um = margin_um)
}

#' Tabella dei quadratini di una griglia
#'
#' @param grid uscita di `hd_grid_geometry()`.
#' @return data.frame `bin_id` (= i * n_col + j + 1), `i`, `j`, `x`, `y` (centro, um), `i8`, `j8`, `bin8_id`.
#' @export
hd_bin_table <- function(grid) {
  i <- rep(seq_len(grid$n_row) - 1L, each = grid$n_col); j <- rep(seq_len(grid$n_col) - 1L, times = grid$n_row)
  i8 <- i %/% grid$block; j8 <- j %/% grid$block
  data.frame(bin_id = i * grid$n_col + j + 1L, i = i, j = j,
             x = grid$x0 + (j + 0.5) * grid$bin_um, y = grid$y0 + (i + 0.5) * grid$bin_um,
             i8 = i8, j8 = j8, bin8_id = i8 * (grid$n_col %/% grid$block) + j8 + 1L)
}

# quadrati dei quadratini (i, j) come geometrie geos (indici 0-based)
.hd_squares <- function(grid, i, j) {
  b <- grid$bin_um; x <- grid$x0 + j * b; y <- grid$y0 + i * b; n <- length(i)
  xs <- c(rbind(x, x + b, x + b, x, x)); ys <- c(rbind(y, y, y + b, y + b, y))
  geos::geos_make_polygon(xs, ys, feature_id = rep(seq_len(n), each = 5L))
}
.hd_sfc <- function(g) sf::st_as_sfc(geos::geos_write_wkb(g))

# tessere: intervalli di righe e colonne (0-based, estremo superiore escluso) di lato tile_um
.hd_tiles <- function(grid, tile_um) {
  nt <- if (is.null(tile_um) || !is.finite(tile_um)) max(grid$n_row, grid$n_col) else max(1L, as.integer(round(tile_um / grid$bin_um)))
  ri <- seq(0L, grid$n_row - 1L, by = nt); cj <- seq(0L, grid$n_col - 1L, by = nt)
  g <- expand.grid(k = seq_along(cj), h = seq_along(ri))
  data.frame(i0 = ri[g$h], i1 = pmin(ri[g$h] + nt, grid$n_row), j0 = cj[g$k], j1 = pmin(cj[g$k] + nt, grid$n_col))
}

# intervalli di righe e colonne toccati dai riquadri delle geometrie (bb: matrice n x 4 xmin ymin xmax ymax)
.hd_bbox_bins <- function(grid, bb) {
  b <- grid$bin_um
  m <- cbind(i0 = pmax(0, floor((bb[, 2] - grid$y0) / b)), i1 = pmin(grid$n_row - 1, floor((bb[, 4] - grid$y0) / b)),
             j0 = pmax(0, floor((bb[, 1] - grid$x0) / b)), j1 = pmin(grid$n_col - 1, floor((bb[, 3] - grid$x0) / b)))
  storage.mode(m) <- "integer"   # indici interi: un bin_id double come 2e5 diventa "2e+05" come testo (fattori, nomi), S2.1
  m
}
.hd_bboxes <- function(g) {
  gg <- geos::as_geos_geometry(sf::st_set_crs(g, NA))
  cbind(geos::geos_xmin(gg), geos::geos_ymin(gg), geos::geos_xmax(gg), geos::geos_ymax(gg))
}

#' Area di tessuto per quadratino
#'
#' Unione delle regioni (che non si sovrappongono): quadratini interi dentro una regione valgono `bin_um`^2, quelli fuori 0, quelli di
#' bordo l'area dell'intersezione con `overlay_prec()` (origine locale della tessera).
#' @param grid uscita di `hd_grid_geometry()`.
#' @param region_polygons `sfc` delle regioni.
#' @param tile_um lato delle tessere (um); `Inf` = una sola tessera.
#' @param grid_um passo dell'arrotondamento degli overlay.
#' @return vettore numerico lungo n_row * n_col (ordine di `bin_id`), area in um^2.
#' @export
hd_tissue_area <- function(grid, region_polygons, tile_um = 500, grid_um = 1e-6) {
  out <- numeric(grid$n_row * grid$n_col); a_full <- grid$bin_um^2
  R <- geos::as_geos_geometry(sf::st_set_crs(region_polygons, NA)); bbR <- .hd_bboxes(region_polygons)
  for (r in seq_along(R)) {
    rg <- .hd_bbox_bins(grid, bbR[r, , drop = FALSE]); if (rg[1, "i1"] < rg[1, "i0"] || rg[1, "j1"] < rg[1, "j0"]) next
    tl <- .hd_tiles(grid, tile_um)
    for (t in seq_len(nrow(tl))) {
      i0 <- max(tl$i0[t], rg[1, "i0"]); i1 <- min(tl$i1[t] - 1L, rg[1, "i1"]); j0 <- max(tl$j0[t], rg[1, "j0"]); j1 <- min(tl$j1[t] - 1L, rg[1, "j1"])
      if (i1 < i0 || j1 < j0) next
      ii <- rep(i0:i1, each = j1 - j0 + 1L); jj <- rep(j0:j1, times = i1 - i0 + 1L)
      sq <- .hd_squares(grid, ii, jj); id <- ii * grid$n_col + jj + 1L
      inside <- geos::geos_prepared_contains(R[r], sq); touch <- geos::geos_prepared_intersects(R[r], sq)
      out[id[inside]] <- out[id[inside]] + a_full
      bd <- which(touch & !inside)
      if (length(bd)) {
        org <- c(grid$x0 + tl$j0[t] * grid$bin_um, grid$y0 + tl$i0[t] * grid$bin_um)
        out[id[bd]] <- out[id[bd]] + overlay_prec("intersection", .hd_sfc(sq[bd]), region_polygons[r], grid_um = grid_um, origin = org)
      }
    }
  }
  out
}

#' Verita' geometrica: area di ogni quadratino coperta da ogni territorio e da ogni nucleo
#'
#' Coppie candidate dai riquadri dei territori, per tessere (ogni quadratino in una sola tessera; origine locale della tessera, BL-104);
#' aree con `overlay_prec()`. Restano le coppie con area del territorio > 0.
#' @param grid uscita di `hd_grid_geometry()`.
#' @param territories,nuclei `sfc` della stessa lunghezza (cellula k = elemento k); `nuclei` puo' essere NULL.
#' @param tile_um lato delle tessere (um); `Inf` = una sola tessera.
#' @param grid_um passo dell'arrotondamento degli overlay.
#' @return data.frame `bin_id`, `cell_id`, `area_territory_um2`, `area_nucleus_um2`, ordinato per `bin_id`, `cell_id`.
#' @export
hd_overlap <- function(grid, territories, nuclei = NULL, tile_um = 500, grid_um = 1e-6) {
  n <- length(territories)
  if (!is.null(nuclei) && length(nuclei) != n) stop("hd_overlap(): territori e nuclei di lunghezza diversa.", call. = FALSE)
  rg <- .hd_bbox_bins(grid, .hd_bboxes(territories))
  ni <- rg[, "i1"] - rg[, "i0"] + 1L; nj <- rg[, "j1"] - rg[, "j0"] + 1L; ok <- ni > 0 & nj > 0
  cand_cell <- rep(seq_len(n)[ok], (ni * nj)[ok])
  off <- sequence((ni * nj)[ok]) - 1L
  ci <- rg[cand_cell, "i0"] + off %/% nj[cand_cell]; cj <- rg[cand_cell, "j0"] + off %% nj[cand_cell]
  nt <- if (is.null(tile_um) || !is.finite(tile_um)) max(grid$n_row, grid$n_col) else max(1L, as.integer(round(tile_um / grid$bin_um)))
  tile <- (ci %/% nt) * (grid$n_col %/% nt + 1L) + cj %/% nt
  res <- vector("list", 0L)
  for (t in sort(unique(tile))) {
    k <- which(tile == t); i <- ci[k]; j <- cj[k]; cc <- cand_cell[k]
    org <- c(grid$x0 + (min(j) %/% nt) * nt * grid$bin_um, grid$y0 + (min(i) %/% nt) * nt * grid$bin_um)
    sq <- .hd_sfc(.hd_squares(grid, i, j))
    aT <- overlay_prec("intersection", sq, territories[cc], grid_um = grid_um, origin = org)
    keep <- aT > 0
    aN <- numeric(sum(keep))
    if (!is.null(nuclei) && any(keep)) aN <- overlay_prec("intersection", sq[keep], nuclei[cc[keep]], grid_um = grid_um, origin = org)
    res[[length(res) + 1L]] <- data.frame(bin_id = as.integer(i[keep] * grid$n_col + j[keep] + 1L), cell_id = as.integer(cc[keep]),
                                          area_territory_um2 = aT[keep], area_nucleus_um2 = aN)
  }
  out <- if (length(res)) do.call(rbind, res) else
    data.frame(bin_id = integer(0), cell_id = integer(0), area_territory_um2 = numeric(0), area_nucleus_um2 = numeric(0))
  out <- out[order(out$bin_id, out$cell_id), ]; rownames(out) <- NULL
  out
}

#' Somma a 8 um
#'
#' Somma dei 16 figli di 2 um di ogni blocco (i %/% 4, j %/% 4), come Space Ranger (design sezione 8).
#' @param x matrice (anche sparsa) con una colonna per quadratino di 2 um in ordine di `bin_id`, oppure data.frame con `bin_id`.
#' @param grid uscita di `hd_grid_geometry()`.
#' @param cols per i data.frame: colonne da sommare (obbligatorio: nessuna colonna si somma per difetto, revisione G-5 di S2.1).
#' @param by per i data.frame: colonne di raggruppamento oltre al blocco (per esempio `"cell_id"`).
#' @return matrice con una colonna per blocco (ordine di `bin8_id`), oppure data.frame con `bin8_id`, `by` e le somme di `cols`.
#' @export
sum_to_8um <- function(x, grid, cols = NULL, by = NULL) {
  nb <- grid$n_row * grid$n_col; nc8 <- grid$n_col %/% grid$block; n8 <- (grid$n_row %/% grid$block) * nc8
  b8_of <- function(id) { i <- (id - 1L) %/% grid$n_col; j <- (id - 1L) %% grid$n_col; as.integer((i %/% grid$block) * nc8 + j %/% grid$block + 1L) }
  if (is.data.frame(x)) {
    if (is.null(cols) || !all(cols %in% names(x)) || any(cols %in% c("bin_id", by)))
      stop("sum_to_8um(): per un data.frame indicare in `cols` le colonne da sommare (esistenti, diverse da bin_id e da `by`).", call. = FALSE)
    key <- data.frame(bin8_id = b8_of(x$bin_id), x[by])
    s <- stats::aggregate(x[cols], by = key, FUN = sum)
    s <- s[do.call(order, s[c("bin8_id", by)]), c("bin8_id", by, cols)]; rownames(s) <- NULL
    return(s)
  }
  if (ncol(x) != nb) stop("sum_to_8um(): la matrice deve avere una colonna per quadratino di 2 um.", call. = FALSE)
  A <- Matrix::sparseMatrix(i = seq_len(nb), j = b8_of(seq_len(nb)), x = 1, dims = c(nb, n8))
  x %*% A
}

#' Griglia di Visium HD e verita' geometrica di un livello cellulare (Step 2, passo 1)
#'
#' @param cell_layer uscita di `place_cells()`: servono `cell_territories`, `cell_nuclei`, `region_polygons`, `cell_df` e
#'   `metadata$pixel_size_um` (obbligatoria, decisione 10 del design dello Step 1 / C-A).
#' @param margin_um margine di vetro (um, multiplo di 8; 280 = DEC-038).
#' @param bin_um lato del quadratino (um; Visium HD = 2).
#' @param tile_um lato delle tessere di calcolo (um; design sezione 12).
#' @param grid_um passo dell'arrotondamento degli overlay (um).
#' @param verbose messaggi di avanzamento.
#' @return lista `grid`, `bins` (`hd_bin_table()` piu' `tissue_area_um2`, `tissue_frac`, `in_tissue` = frazione >= 0.5), `bins_8um`
#'   (blocchi con area di tessuto e `in_tissue` = almeno un figlio in tessuto, come Space Ranger), `truth$overlap` (piu' le frazioni
#'   `frac_territory`, `frac_nucleus` = area / bin_um^2), `metadata`.
#' @export
build_hd_grid <- function(cell_layer, margin_um = 280, bin_um = 2, tile_um = 500, grid_um = 1e-6, verbose = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  px <- cell_layer$metadata$pixel_size_um
  if (is.null(px) || !is.numeric(px) || length(px) != 1L || !is.finite(px) || px <= 0)
    stop("build_hd_grid(): `cell_layer$metadata$pixel_size_um` manca o non e' valido (decisione 10 / C-A, BL-005): nessun default.", call. = FALSE)
  Tt <- cell_layer$cell_territories; Nn <- cell_layer$cell_nuclei; R <- cell_layer$region_polygons
  if (length(Tt) != nrow(cell_layer$cell_df) || length(Nn) != length(Tt))
    stop("build_hd_grid(): territori, nuclei e cell_df di lunghezza diversa.", call. = FALSE)
  gT <- geos::as_geos_geometry(sf::st_set_crs(Tt, NA))
  if (any(geos::geos_is_empty(gT)))
    stop(sprintf("build_hd_grid(): %d territori vuoti (cell_id %s): ogni cellula dello Step 1 deve avere un territorio.",
                 sum(geos::geos_is_empty(gT)), paste(utils::head(cell_layer$cell_df$cell_id[geos::geos_is_empty(gT)], 5), collapse = ", ")), call. = FALSE)
  bb <- sf::st_bbox(c(sf::st_set_crs(Tt, NA), sf::st_set_crs(R, NA)))
  grid <- hd_grid_geometry(bb, margin_um = margin_um, bin_um = bin_um)
  bins <- hd_bin_table(grid)
  bins$tissue_area_um2 <- hd_tissue_area(grid, R, tile_um = tile_um, grid_um = grid_um)
  bins$tissue_frac <- bins$tissue_area_um2 / bin_um^2
  bins$in_tissue <- bins$tissue_frac >= 0.5
  if (verbose) message(sprintf("build_hd_grid(): %d x %d quadratini, tessuto %.0f um^2", grid$n_row, grid$n_col, sum(bins$tissue_area_um2)))
  ov <- hd_overlap(grid, Tt, Nn, tile_um = tile_um, grid_um = grid_um)
  # il nucleo si calcola solo dove il territorio copre il quadratino (lo Step 1 contiene il nucleo nel territorio, C3): se un nucleo sporge,
  # la parte esterna andrebbe persa in silenzio; si controlla entro la tolleranza a n pezzi (DEC-038) e si avvisa
  nk <- tabulate(ov$cell_id, length(Tt)); sN <- as.numeric(tapply(ov$area_nucleus_um2, factor(ov$cell_id, levels = seq_along(Tt)), sum)); sN[is.na(sN)] <- 0
  gN <- geos::as_geos_geometry(sf::st_set_crs(Nn, NA)); aN <- geos::geos_area(gN); PN <- geos::geos_length(gN)
  lost <- which(aN - sN > grid_um * (8 * nk + 2 * PN) + (nk + 1) * grid_um^2)
  if (length(lost)) warning(sprintf("build_hd_grid(): %d nuclei sporgono dal loro territorio (cell_id %s): la parte esterna non entra in truth$overlap.",
                                    length(lost), paste(utils::head(cell_layer$cell_df$cell_id[lost], 5), collapse = ", ")), call. = FALSE)
  ov$cell_id <- cell_layer$cell_df$cell_id[ov$cell_id]
  ov$frac_territory <- ov$area_territory_um2 / bin_um^2; ov$frac_nucleus <- ov$area_nucleus_um2 / bin_um^2
  b8 <- sum_to_8um(data.frame(bin_id = bins$bin_id, tissue_area_um2 = bins$tissue_area_um2, n_in_tissue = as.numeric(bins$in_tissue)), grid,
                   cols = c("tissue_area_um2", "n_in_tissue"))
  nc8 <- grid$n_col %/% grid$block
  b8$i8 <- (b8$bin8_id - 1L) %/% nc8; b8$j8 <- (b8$bin8_id - 1L) %% nc8
  b8$x <- grid$x0 + (b8$j8 + 0.5) * bin_um * grid$block; b8$y <- grid$y0 + (b8$i8 + 0.5) * bin_um * grid$block
  b8$in_tissue <- b8$n_in_tissue > 0
  b8 <- b8[, c("bin8_id", "i8", "j8", "x", "y", "tissue_area_um2", "n_in_tissue", "in_tissue")]
  list(grid = grid, bins = bins, bins_8um = b8, truth = list(overlap = ov),
       metadata = list(pixel_size_um = px, margin_um = margin_um, bin_um = bin_um, tile_um = tile_um, grid_um = grid_um,
                       n_cells = length(Tt), step1_random_seed = cell_layer$metadata$random_seed,
                       step1_preset = cell_layer$metadata$tissue_preset_used,
                       package_version = "geo.spatialtrans step2-visium-grid (S2.1)", elapsed_s = proc.time()[["elapsed"]] - t0))
}
