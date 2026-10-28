# tools/S1.4e_bench_lib.R - S1.4e: bracci del disegno del nucleo (K2, ER, ET), ridisegno dei nuclei di un livello di place_cells(),
# invarianti geometriche, centri del banco A (punti di Luca nella finestra, centri del segmentatore altrove). Definizioni: docs/steps/S1.4e.md.
# Richiede: R/04b*.R, tools/R3_voronoi_metrics.R gia' caricati.

#' Eccentricita' mediana vera per tessuto (B-115...B-120, segmentatore primario) e allungamento = asse maggiore / minore = 1 / sqrt(1 - e^2).
S14E_ECC <- c(A1 = 0.80, A2 = 0.78, A3 = 0.81, A4 = 0.60, A5 = 0.56, A6 = 0.81)
s14e_aspect <- function(A) 1 / sqrt(1 - S14E_ECC[[A]]^2)
S14E_ARMS <- c("K2", "ER", "ET")

#' Catalogo dei tipi per un braccio: K2 invariato; ER ellisse orientata a caso; ET ellisse lungo l'asse maggiore del territorio.
#' Stesso allungamento per tutti i tipi del tessuto (semplificazione dichiarata nello step).
arm_catalog <- function(ct, arm, A) {
  arm <- match.arg(arm, S14E_ARMS)
  if (arm == "K2") return(ct)
  ct$nucleus_aspect <- s14e_aspect(A)
  ct$nucleus_orient <- if (arm == "ET") "territory" else "random"
  if ("nucleus_shape" %in% names(ct)) ct$nucleus_shape <- "ellipse"      # coerenza richiesta da derive_cell_geometry()
  ct
}

#' Ridisegna i nuclei di un livello di place_cells() con un catalogo dei tipi: stessi centri, territori, seme dello stadio "size",
#' modello del nucleo, tabella dei residui e modo dei residui. Con il catalogo del livello riproduce i nuclei (T1 di test_S1.4e.R).
redraw_nuclei <- function(L, ct) {
  cd <- L$cell_df
  cen <- cd[, intersect(c("cell_id", "region_id", "cluster_id", "x", "y", "cell_type"), names(cd))]
  tess <- list(cell_territories = L$cell_territories, territory_df = L$territory_df)
  rtab <- if ("nucleus_area_resid" %in% names(ct) && any(!is.na(ct$nucleus_area_resid))) .tp_resid_table() else NULL
  derive_cell_geometry(cen, tess, ct, random_seed = L$metadata$stage_seeds[["size"]], nucleus_model = L$metadata$nucleus_model,
                       verbose = FALSE, resid_table = rtab, resid_mode = L$metadata$resid_mode)
}

#' Invarianti del cancello G-3: dentro il proprio territorio, nessuna coppia sovrapposta (intersezione > 1e-6 um^2),
#' |A - A*| <= max(1e-6 A*, 1e-6 um^2) con A* = nucleus_area_cap, geometrie valide e non vuote.
nuc_invariants <- function(o, terr) {
  g <- geos::as_geos_geometry(o$cell_nuclei); gt <- geos::as_geos_geometry(terr); cd <- o$cell_df
  wi <- geos::geos_within(g, gt)
  tr <- geos::geos_strtree(g); q <- geos::geos_strtree_query(tr, g)
  pi_ <- rep(seq_along(q), lengths(q)); pj <- unlist(q); k <- pi_ < pj
  ov <- if (any(k)) geos::geos_area(geos::geos_intersection(g[pi_[k]], g[pj[k]])) else 0
  a <- geos::geos_area(g); as <- cd$nucleus_area_cap
  data.frame(n = nrow(cd), n_out = sum(!wi), n_overlap = sum(ov > 1e-6), max_overlap = max(c(0, ov)),
             n_area_bad = sum(abs(a - as) > pmax(1e-6 * as, 1e-6)), max_area_rel = max(abs(a - as) / as),
             n_invalid = sum(!geos::geos_is_valid(g)), n_empty = sum(geos::geos_is_empty(g)))
}

#' Esclusioni annotate da Luca in una finestra di R2b (um del ROI), come tools/R3b_bench.R::excl_um.
s14e_excl_um <- function(wv, ann) {
  aj <- ann[[which(vapply(ann, function(z) z$win_id, "") == wv$win_id)]]; m0 <- aj$margin_px; upp <- aj$um_per_px
  if (!length(aj$exclusions)) return(NULL)
  sf::st_union(sf::st_sfc(lapply(aj$exclusions, function(e) {
    p <- do.call(rbind, lapply(e$pts, unlist)); p <- cbind(wv$x0_um + (p[, 1] - m0) * upp, wv$y0_um + (p[, 2] - m0) * upp)
    sf::st_make_valid(sf::st_polygon(list(rbind(p, p[1, ])))) })))
}

#' Centri del banco A per una finestra: i punti di Luca (R2b) dentro la maschera del ROI; i centri del segmentatore (co: x, y) fuori dal
#' quadrato della finestra o dentro le esclusioni di Luca. Come tools/R3b_bench.R::cp1_one. Restituisce x, y, src, luca_idx.
centers_A <- function(wv, co, pts, w, ann) {
  s <- wv$side_um; inwin <- function(x, y) x >= wv$x0_um & x < wv$x0_um + s & y >= wv$y0_um & y < wv$y0_um + s
  ex <- s14e_excl_um(wv, ann)
  inex <- function(x, y) if (is.null(ex)) rep(FALSE, length(x)) else
    lengths(sf::st_intersects(sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y")), ex)) > 0
  keep_co <- !inwin(co$x, co$y) | inex(co$x, co$y)
  ok <- r3_in_mask(pts$x_um, pts$y_um, w$m, w$px)
  data.frame(x = c(co$x[keep_co], pts$x_um[ok]), y = c(co$y[keep_co], pts$y_um[ok]),
             src = c(rep("CO", sum(keep_co)), rep("manual", sum(ok))), luca_idx = c(rep(NA_integer_, sum(keep_co)), pts$idx[ok]),
             stringsAsFactors = FALSE)
}
