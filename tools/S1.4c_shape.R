# tools/S1.4c_shape.R — S1.4c: nuclei ritagliati nel territorio (K-1, K-2) e operatore di osservazione dei simulati
# (Voronoi dei centroidi osservati con tessellate_voronoi(), lo stesso operatore dei reali di S1.4 `real_one`).
# Prototipo: R/04b4* non si tocca (implementazione in S1.4d). Mutanti (argomento `mutant`): MU1 tetto ignorato,
# MU2 ritaglio contro il territorio del vicino, MU5 territori osservati dai generatori invece che dai centroidi.
S14C_GRID <- 1e-6

#' Poligono regolare a nv vertici di raggio circoscritto R (stessi angoli di .dcg_polygons: 2*pi*k/nv).
kc_ngon <- function(x, y, R, nv = 32L) {
  t <- 2 * pi * (seq_len(nv) - 1) / nv
  X <- x + outer(R, cos(t)); Y <- y + outer(R, sin(t))
  geos::geos_make_polygon(as.vector(t(X)), as.vector(t(Y)), feature_id = rep(seq_along(x), each = nv))
}
#' Fattore del raggio circoscritto: un nv-gono di raggio kp * sqrt(A / pi) ha area A.
kc_poly_factor <- function(nv) sqrt(2 * pi / (nv * sin(2 * pi / nv)))

#' Nuclei ritagliati. K1: nv-gono di area a_tgt al generatore, intersecato con il territorio (griglia F1).
#' K2: come K1, raggio aumentato per bisezione finche' area(N) = a_star = min(a_tgt, phi * a_t) (tetto).
#' tol = 1e-7: l'errore dello snap sulla griglia F1 e' ~2e-8 relativo su 50 µm² (sviluppo); soglia C-area 1e-6.
#' @param terr_g territori (geos, gia' sulla griglia S14C_GRID), ordine delle celle. @param inset rientro (µm) del ritaglio.
#' @return list(nuc = geos, df = data.frame(area, a_star, capped, R, n_iter))
kclip_nuclei <- function(x, y, a_tgt, terr_g, mode = c("K1", "K2"), a_t = NULL, phi = 1, nv = 32L,
                         tol = 1e-7, max_iter = 200L, mutant = "none", inset = 1e-5) {
  mode <- match.arg(mode); n <- length(x); kp <- kc_poly_factor(nv)
  stopifnot(length(y) == n, length(a_tgt) == n, length(terr_g) == n, all(a_tgt > 0))
  ti <- seq_len(n); if (mutant == "MU2" && n > 1) ti <- c(2:n, 1L)
  # ritaglio contro il territorio rientrato di `inset` (dev notes D-1): lo snap-rounding della griglia F1 sposta i vertici
  # nuovi fino a ~1e-6 µm fuori dal territorio (C3 36/1 000 in sviluppo); con 1e-5 il nucleo resta dentro per costruzione
  # il buffer si calcola SENZA modello di precisione: con la griglia 1e-6 ereditata, geos_buffer(-1e-5) restituisce in silenzio
  # poligoni vuoti (23/674 in sviluppo, dev notes D-2); poi si torna sulla griglia F1
  Tg <- geos::geos_set_precision(geos::geos_buffer(geos::geos_set_precision(terr_g[ti], 0), -inset,
                                                   params = geos::geos_buffer_params(quad_segs = 2L)), S14C_GRID)
  if (any(geos::geos_is_empty(Tg))) stop("kclip_nuclei(): territorio rientrato vuoto per ", sum(geos::geos_is_empty(Tg)), " celle")
  cap_on <- mode == "K2" && mutant != "MU1"
  if (cap_on) stopifnot(length(a_t) == n, phi > 0, phi <= 1)
  # tetto sull'area del territorio rientrato (differisce da a_t di <= perimetro * inset): con phi = 1 il target resta raggiungibile
  a_star <- if (cap_on) pmin(a_tgt, phi * pmin(a_t, geos::geos_area(Tg))) else a_tgt
  R <- kp * sqrt(a_star / pi)
  nuc <- geos::geos_intersection_prec(kc_ngon(x, y, R, nv), Tg, S14C_GRID)
  area <- geos::geos_area(nuc); n_iter <- integer(n)
  if (mode == "K2") {
    need <- which(area < a_star * (1 - tol))
    if (length(need)) {
      xs <- x[need]; ys <- y[need]; Tn <- Tg[need]; as_ <- a_star[need]
      ext <- geos::geos_extent(Tn)
      dmax <- sqrt(pmax((xs - ext$xmin)^2, (xs - ext$xmax)^2) + pmax((ys - ext$ymin)^2, (ys - ext$ymax)^2))
      lo <- R[need]; hi <- dmax / cos(pi / nv) * (1 + 1e-9)
      a_hi <- geos::geos_area(geos::geos_intersection_prec(kc_ngon(xs, ys, hi, nv), Tn, S14C_GRID))
      a_lo <- area[need]; it <- integer(length(need))
      for (k in seq_len(max_iter)) {
        act <- which((a_hi - as_) / as_ > tol & (hi - lo) / hi > 1e-15)
        if (!length(act)) break
        mid <- (lo[act] + hi[act]) / 2
        am <- geos::geos_area(geos::geos_intersection_prec(kc_ngon(xs[act], ys[act], mid, nv), Tn[act], S14C_GRID))
        up <- am < as_[act]
        lo[act[up]] <- mid[up]; a_lo[act[up]] <- am[up]
        hi[act[!up]] <- mid[!up]; a_hi[act[!up]] <- am[!up]
        it[act] <- it[act] + 1L
      }
      best <- ifelse(abs(a_lo - as_) < abs(a_hi - as_), lo, hi)
      R[need] <- best; n_iter[need] <- it
      nuc[need] <- geos::geos_intersection_prec(kc_ngon(xs, ys, best, nv), Tn, S14C_GRID)
      area[need] <- geos::geos_area(nuc[need])
    }
  }
  list(nuc = nuc, df = data.frame(area = area, a_star = a_star, capped = a_star < a_tgt, R = R, n_iter = n_iter))
}

#' Territori osservati: Voronoi dei centroidi (cx, cy) dei nuclei osservati, ritagliato alla regione, con
#' tessellate_voronoi() (operatore dei reali). Centroidi fuori dalla regione o mancanti: esclusi (NA, contati).
#' @param region list(area_um2, poly) come r3_roi_window(). @param gx,gy generatori (solo per il mutante MU5).
#' @return data.frame(a_t, d, interior_obs) per cella (NA dove escluso) e attributo n_excluded.
obs_territories <- function(cx, cy, region, dcg, mutant = "none", gx = NULL, gy = NULL) {
  n <- length(cx)
  if (mutant == "MU5") { cx <- gx; cy <- gy }
  ok <- is.finite(cx) & is.finite(cy)
  ok[ok] <- geos::geos_intersects(geos::geos_make_point(cx[ok], cy[ok]), geos::as_geos_geometry(region$poly))
  id <- which(ok)
  cen <- data.frame(cell_id = id, region_id = 1L, x = cx[id], y = cy[id], cell_type = "c")
  reg <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = region$area_um2),
              region_polygons = sf::st_sfc(region$poly[[1]]))
  tv <- tessellate_voronoi(cen, reg, verbose = FALSE)
  td <- tv$territory_df; stopifnot(identical(as.integer(td$cell_id), as.integer(id)))
  terr <- geos::as_geos_geometry(dcg$.dcg_as_polys(tv$cell_territories))
  a_t <- rep(NA_real_, n); d <- rep(NA_real_, n); cl <- rep(NA, n)
  a_t[id] <- td$territory_area; cl[id] <- td$clipped
  d[id] <- geos::geos_distance(geos::geos_make_point(cx[id], cy[id]), geos::geos_boundary(terr))
  out <- data.frame(a_t = a_t, d = d, clipped = cl); attr(out, "n_excluded") <- sum(!ok); out
}

#' B-8: frazione delle celle `ev` il cui cerchio di pari area (centro, area) interseca il cerchio di un'altra cella.
b8_frac <- function(x, y, area, ev, L) {
  ok <- is.finite(x) & is.finite(area) & area > 0; id <- which(ok)
  p <- L$near_pairs(x[id], y[id], sqrt(area[id] / pi), gamma = 1)
  hit <- rep(FALSE, length(x)); hit[id[c(p[, "i"], p[, "j"])]] <- TRUE
  mean(hit[ev & ok])
}
