# tools/R3_voronoi_metrics.R — funzioni di R3 (Voronoi sui nuclei reali), riusabili dai check B di S1.3–S1.4.
# Pre-registrazione: results/R3/R3_preregistration.md (commit 3fac4e7). Convenzione: µm, y verso il basso
# (righe dell'immagine), come extract_regions() e R2_nuclei_all.parquet. Nessuna dipendenza oltre sf, deldir.
# Mutanti (C-R3.7) attivati da Sys.getenv("R3_MUTANT"): M1 tile permutati, M2 niente ritaglio,
# M4 flag interne ignorato, M5 rapporto delle deviazioni standard al posto delle varianze in e_T
# (M3 e' nel contenimento, tools/R3_containment.py).

r3_mutant <- function() Sys.getenv("R3_MUTANT", "none")

#' Momenti esatti di un poligono semplice (anello aperto, qualunque verso): area, centroide,
#' momenti centrali del secondo ordine (per unita' di area).
r3_poly_moments <- function(x, y) {
  x0 <- mean(x); y0 <- mean(y); x <- x - x0; y <- y - y0     # origine locale: evita la cancellazione in m - c^2
  x2 <- c(x[-1], x[1]); y2 <- c(y[-1], y[1])
  cr <- x * y2 - x2 * y
  A <- sum(cr) / 2
  cx <- sum((x + x2) * cr) / (6 * A); cy <- sum((y + y2) * cr) / (6 * A)
  Ixx <- sum((x^2 + x * x2 + x2^2) * cr) / 12
  Iyy <- sum((y^2 + y * y2 + y2^2) * cr) / 12
  Ixy <- sum((x * y2 + 2 * x * y + 2 * x2 * y2 + x2 * y) * cr) / 24
  c(area = abs(A), cx = cx + x0, cy = cy + y0, mxx = Ixx / A - cx^2, myy = Iyy / A - cy^2, mxy = Ixy / A - cx * cy)
}

#' Eccentricita' (definizione di regionprops: sqrt(1 - l2/l1), autovalori del tensore dei momenti) e
#' orientazione dell'asse maggiore (rad, da +x verso +y, y in basso; in (-pi/2, pi/2]).
r3_shape <- function(mxx, myy, mxy) {
  h <- sqrt(((mxx - myy) / 2)^2 + mxy^2); tr2 <- (mxx + myy) / 2
  l1 <- tr2 + h; l2 <- pmax(tr2 - h, 0)
  ratio <- ifelse(l1 > 0, l2 / l1, 1)
  if (r3_mutant() == "M5") ratio <- sqrt(ratio)
  list(ecc = sqrt(pmax(1 - ratio, 0)), theta = 0.5 * atan2(2 * mxy, mxx - myy))
}

#' Differenza assiale fra due orientazioni (rad) in [0, pi/2].
r3_axis_diff <- function(a, b) { d <- abs(a - b) %% pi; pmin(d, pi - d) }

#' Poligono della maschera valida di un ROI R2 (come tools/S1.2_checkB.R).
r3_roi_window <- function(A, roi, rois) {
  side_um <- rois$side_um[rois$archetype == A & rois$roi_id == roi]
  m <- png::readPNG(sprintf("results/R2/tissue_masks/%s_%s_valid_ds4.png", A, roi))
  if (length(dim(m)) == 3) m <- m[, , 1]
  m <- m > 0.5; stopifnot(nrow(m) == ncol(m))
  px <- side_um / ncol(m)
  idx <- which(m, arr.ind = TRUE)
  clust <- data.frame(x = idx[, 2], y = idx[, 1], cl = 1L)
  reg <- extract_regions(clust, pixel_size_um = px, min_region_area_um2 = 0, cluster_col = "cl",
                         simplify_tol_um = 0, stride = 1L, verbose = FALSE)
  poly <- sf::st_union(reg$region_polygons)
  frame <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(side_um, 0), c(side_um, side_um), c(0, side_um), c(0, 0)))))
  list(m = m, px = px, side_um = side_um, reg = reg, poly = poly, frame = frame,
       n_valid_px = sum(m), area_um2 = as.numeric(sf::st_area(poly)))
}

#' Centroidi dentro la maschera (pixel valido, convenzione dei lati come extract_regions()).
r3_in_mask <- function(x, y, m, px) {
  ci <- floor(x / px) + 1; ri <- floor(y / px) + 1
  ok <- ci >= 1 & ri >= 1 & ci <= ncol(m) & ri <= nrow(m)
  out <- rep(FALSE, length(x)); out[ok] <- m[cbind(ri[ok], ci[ok])]; out
}

#' Tassellazione di Voronoi ritagliata.
#' @param x,y generatori (µm); @param clip sfc (1 geometria) regione di ritaglio; @param frame sfc quadrato
#' del ROI/finestra (una cella interna deve stare dentro); @param pad allargamento del riquadro deldir.
#' @return data.frame per generatore (ordine di ingresso): area clip, area tile, interna, lati, forma.
r3_tessellate <- function(x, y, clip, frame, pad = 50) {
  n <- length(x); stopifnot(n >= 3, length(clip) == 1)
  bb <- sf::st_bbox(frame)
  rw <- c(bb[["xmin"]] - pad, bb[["xmax"]] + pad, bb[["ymin"]] - pad, bb[["ymax"]] + pad)
  dd <- deldir::deldir(x, y, rw = rw, round = FALSE)
  tl <- deldir::tile.list(dd)
  pt <- vapply(tl, function(t) t$ptNum, 1L)
  if (length(tl) != n || any(sort(pt) != seq_len(n))) stop("r3_tessellate(): deldir ha scartato o duplicato generatori")
  tl <- tl[order(pt)]
  if (r3_mutant() == "M1") tl <- tl[c(2:n, 1)]
  mom <- t(vapply(tl, function(t) r3_poly_moments(t$x, t$y), numeric(6)))
  nsides <- vapply(tl, function(t) length(t$x), 1L)
  tiles <- sf::st_sfc(lapply(tl, function(t) sf::st_polygon(list(cbind(c(t$x, t$x[1]), c(t$y, t$y[1]))))))
  within <- lengths(sf::st_within(tiles, clip)) > 0
  a_clip <- mom[, "area"]
  nb <- which(!within)
  if (length(nb) && r3_mutant() != "M2") {
    cg <- clip[[1]]
    a_clip[nb] <- vapply(nb, function(i) {
      tb <- sf::st_bbox(tiles[i])
      loc <- suppressWarnings(sf::st_crop(sf::st_sfc(cg), tb))
      if (length(loc) == 0 || sf::st_is_empty(loc)) return(0)
      g <- suppressWarnings(sf::st_intersection(tiles[i], loc))
      if (length(g) == 0) 0 else sum(as.numeric(sf::st_area(g)))
    }, 0)
  }
  in_frame <- lengths(sf::st_within(tiles, frame)) > 0
  interior <- within | (in_frame & a_clip / mom[, "area"] > 0.999)
  if (r3_mutant() == "M4") interior[] <- TRUE
  sh <- r3_shape(mom[, "mxx"], mom[, "myy"], mom[, "mxy"])
  own <- mapply(function(t, xi, yi) sp_in_poly(xi, yi, t$x, t$y), tl, x, y)
  data.frame(idx = seq_len(n), x = x, y = y, area_tile = mom[, "area"], area = a_clip, interior = interior,
             nsides = nsides, ecc_T = sh$ecc, theta_T = sh$theta, cx_T = mom[, "cx"], cy_T = mom[, "cy"],
             gen_in_own = own)
}

#' Punto in poligono (ray casting; i tile di Voronoi sono convessi e il generatore e' interno).
sp_in_poly <- function(px, py, x, y) {
  n <- length(x); j <- c(n, seq_len(n - 1)); xi <- x; yi <- y; xj <- x[j]; yj <- y[j]
  sum(((yi > py) != (yj > py)) & (px < (xj - xi) * (py - yi) / (yj - yi) + xi)) %% 2 == 1
}

#' Intensita' locale al generatore (kernel gaussiano, banda 5/sqrt(lambda), correzione di bordo).
#' Finestra spatstat = maschera (y verso l'alto: y_up = side - y).
r3_local_intensity <- function(x, y, w) {
  W <- spatstat.geom::owin(mask = w$m[nrow(w$m):1, ], xrange = c(0, w$side_um), yrange = c(0, w$side_um))
  X <- spatstat.geom::ppp(x, w$side_um - y, window = W, check = FALSE)
  lam <- X$n / spatstat.geom::area(W)
  sig <- 5 / sqrt(lam)
  d <- spatstat.explore::density.ppp(X, sigma = sig, at = "points", edge = TRUE, leaveoneout = FALSE)
  list(lambda_loc = as.numeric(d), sigma = sig, lambda = lam, X = X)
}

#' Riassunto di un pattern tassellato (solo interne): usato per reale e nulli.
r3_summary <- function(tv, lam_loc = NULL) {
  i <- tv$interior; a <- tv$area[i]
  out <- c(n = nrow(tv), n_interior = sum(i), frac_interior = mean(i), mean_area = mean(a), median_area = median(a),
           cv = sd(a) / mean(a), median_eq_r = median(sqrt(a / pi)), median_ecc_T = median(tv$ecc_T[i]),
           mean_nsides = mean(tv$nsides[i]), var_nsides = var(tv$nsides[i]), var_norm_area = var(a / mean(a)))
  if (!is.null(lam_loc)) { al <- a * lam_loc[i]; out["cv_loc"] <- sd(al) / mean(al) }
  out
}
