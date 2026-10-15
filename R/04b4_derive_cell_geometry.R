# R/04b4_derive_cell_geometry.R — Step 1, fase 4: nucleo (centro, taglia, forma, clamp) per ogni territorio.
# Design: docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md §4.4, §5.4 (v1.4, S1.4: D-S1.4.1, D-S1.4.2).
# Pre-registrazione: results/S1.4/S1.4_preregistration.md (commit c72c01a).
# Convenzioni: µm, y verso il basso (come extract_regions()); angoli in rad da +x verso +y, in (-pi/2, pi/2].
# Dipendenze: sf, geos, polylabelr.

#' Geometria del nucleo per ogni cellula
#'
#' Modello di taglia (D-S1.4.1), per tipo t con `A_ref = 1e6 / density_t`:
#' `log A_N = log(rho^2 * A_ref) + beta * (log A_T - log A_ref) + sigma * z`, z ~ N(0, 1);
#' con `nucleus_radius` assoluto: `log A_N = log(pi * r^2) + sigma * z`. Con beta = 1, sigma = 0 il modello
#' coincide con il design (`nucleus_radius = rho * eq_radius`).
#' Forma: poligono di `n_vertices` vertici su un'ellisse di semiassi a = sqrt(A * asp / pi), b = sqrt(A / (asp * pi)),
#' moltiplicati per k = sqrt(2 pi / (n sin(2 pi / n))) cosi' che l'area del poligono sia esattamente A.
#' Clamp (D-S1.4.2, in tutti i modi): fattore di scala uniforme s <= 1. Cerchio: raggio circoscritto
#' k * a * s <= clamp_factor * d (d = distanza del centro dal bordo del territorio). Ellisse: s parte dal valore del
#' cerchio circoscritto e cresce per bisezione finche' il poligono dilatato di 1/clamp_factor resta dentro il territorio.
#'
#' @param centroids data.frame con `cell_id` (unico), `region_id`, `x`, `y`, `cell_type`; `cluster_id` opzionale.
#' @param territories output di `tessellate_voronoi()`: `cell_territories` (sfc, ordine di `cell_id` crescente) e
#'   `territory_df` (`cell_id`, `territory_area`).
#' @param cell_types catalogo: `cell_type`, `density` (cell/mm²), `nucleus_to_eq_ratio` (rho) e/o `nucleus_radius`
#'   (µm, NA = usa rho); opzionali `nucleus_area_slope` (beta, default 1), `nucleus_area_sdlog` (sigma, default 0),
#'   `nucleus_aspect` (>= 1, default 1 = cerchio), `nucleus_orient` (`"random"` default, `"territory"` = asse maggiore
#'   del territorio, oppure un angolo in rad, numerico o stringa numerica).
#' @param random_seed seme; un solo flusso: `z <- rnorm(n)` e poi `u <- runif(n)`, in ordine di `cell_id`.
#' @param nucleus_center `"poi"` (polo di inaccessibilita', design) o `"generator"` (il centroide seminato).
#' @param n_vertices vertici del poligono del nucleo (>= 8).
#' @param clamp_factor frazione della distanza dal bordo concessa al nucleo (design: 0.85), in (0, 1).
#' @param poi_precision_um precisione di `polylabelr::poi()` (µm).
#' @param strict se TRUE (default) errore se un nucleo non e' dentro il proprio territorio (C3).
#' @param verbose stampa un riepilogo.
#' @return lista con `cell_df` (schema del design §4.2 piu' `nucleus_area_target`, `nucleus_scale`, `nucleus_clamped`,
#'   `center_dist`), `cell_nuclei` (sfc_POLYGON, ordine di `cell_id`), `info`.
derive_cell_geometry <- function(centroids, territories, cell_types, random_seed = 42,
                                 nucleus_center = c("poi", "generator"), n_vertices = 32L, clamp_factor = 0.85,
                                 poi_precision_um = 1e-3, strict = TRUE, verbose = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  nucleus_center <- match.arg(nucleus_center)
  cen <- .dcg_check_args(centroids, territories, n_vertices, clamp_factor, poi_precision_um)
  ct <- .dcg_check_types(cell_types, unique(cen$cell_type))
  terr <- territories$cell_territories
  n <- nrow(cen)
  a_t <- territories$territory_df$territory_area[match(cen$cell_id, territories$territory_df$cell_id)]
  if (any(!is.finite(a_t) | a_t <= 0)) stop("derive_cell_geometry(): territory_area non positiva o mancante")
  ti <- match(cen$cell_type, ct$cell_type)
  set.seed(random_seed); z <- stats::rnorm(n); u <- stats::runif(n)
  gt <- geos::as_geos_geometry(terr)
  ctr <- .dcg_center(terr, cen, nucleus_center, poi_precision_um)
  pts <- geos::geos_make_point(ctr$x, ctr$y)
  if (!all(geos::geos_intersects(pts, gt))) stop("derive_cell_geometry(): centro del nucleo fuori dal territorio")
  d <- geos::geos_distance(pts, geos::geos_boundary(gt))
  a_tgt <- .dcg_target_area(a_t, ct, ti, z)
  asp <- ct$nucleus_aspect[ti]
  theta <- .dcg_orientation(terr, ct, ti, u, asp)
  kp <- .dcg_poly_factor(n_vertices)
  a_eq <- sqrt(a_tgt * asp / pi); b_eq <- sqrt(a_tgt / (asp * pi))
  s <- .dcg_clamp(a_eq, b_eq, theta, ctr, d, gt, kp, n_vertices, clamp_factor)
  nuc_g <- .dcg_polygons(ctr$x, ctr$y, kp * a_eq * s, kp * b_eq * s, theta, n_vertices)
  within <- geos::geos_within(nuc_g, gt)
  if (strict && !all(within)) stop(sprintf("derive_cell_geometry(): %d nuclei non contenuti nel territorio (C3)", sum(!within)))
  a_n <- s^2 * a_tgt
  circ <- asp == 1
  cell_df <- data.frame(
    cell_id = cen$cell_id, region_id = cen$region_id,
    cluster_id = if ("cluster_id" %in% names(cen)) cen$cluster_id else NA_integer_,
    cell_type = cen$cell_type, x = cen$x, y = cen$y,
    territory_area = a_t, eq_radius = sqrt(a_t / pi),
    nucleus_shape = ifelse(circ, "circle", "ellipse"),
    nucleus_radius = ifelse(circ, sqrt(a_n / pi), NA_real_),
    nucleus_major_r = ifelse(circ, NA_real_, a_eq * s), nucleus_minor_r = ifelse(circ, NA_real_, b_eq * s),
    nucleus_orient = ifelse(circ, NA_real_, theta),
    nucleus_cx = ctr$x, nucleus_cy = ctr$y, nucleus_area = a_n, cytoplasm_area = a_t - a_n,
    nucleus_area_target = a_tgt, nucleus_scale = s, nucleus_clamped = s < 1, center_dist = d,
    stringsAsFactors = FALSE)
  info <- list(n = n, n_clamped = sum(s < 1), frac_clamped = mean(s < 1), n_not_within = sum(!within),
               nucleus_center = nucleus_center, random_seed = random_seed, n_vertices = as.integer(n_vertices),
               clamp_factor = clamp_factor, poly_factor = kp, poi_precision_um = poi_precision_um,
               elapsed_s = proc.time()[["elapsed"]] - t0)
  if (verbose) message(sprintf("derive_cell_geometry(): %d cellule, %d clampate (%.2f %%), centro = %s, %.1f s",
                               n, info$n_clamped, 100 * info$frac_clamped, nucleus_center, info$elapsed_s))
  list(cell_df = cell_df, cell_nuclei = sf::st_as_sfc(nuc_g), info = info)
}

#' Stima dei parametri del modello di taglia (D-S1.4.1) da aree osservate.
#'
#' OLS di `log A_N` su `log A_T - log A_ref`. `model`: `"full"` (beta, sigma stimati), `"beta1"` (beta = 1),
#' `"beta0"` (beta = 0), `"design"` (beta = 1, sigma = 0). Per `"beta1"`/`"design"` rho e' la media geometrica di
#' sqrt(A_N / A_T); sigma e' la deviazione standard dei residui (denominatore n - numero di parametri stimati).
#' @param a_n,a_t aree di nucleo e territorio (µm², > 0). @param a_ref area di riferimento (µm², = 1e6 / density).
#' @return c(rho, beta, sigma, n).
fit_nucleus_size_model <- function(a_n, a_t, a_ref, model = c("full", "beta1", "beta0", "design")) {
  model <- match.arg(model)
  ok <- is.finite(a_n) & is.finite(a_t) & a_n > 0 & a_t > 0
  y <- log(a_n[ok]); x <- log(a_t[ok]) - log(a_ref); n <- sum(ok)
  if (n < 3) stop("fit_nucleus_size_model(): meno di 3 osservazioni valide")
  if (model == "full") {
    b <- stats::cov(x, y) / stats::var(x); a <- mean(y) - b * mean(x); r <- y - a - b * x; k <- 2
  } else if (model == "beta0") {
    b <- 0; a <- mean(y); r <- y - a; k <- 1
  } else {
    b <- 1; a <- mean(y - x); r <- y - a - x; k <- 1
  }
  sig <- if (model == "design") 0 else sqrt(sum(r^2) / (n - k))
  c(rho = sqrt(exp(a) / a_ref), beta = b, sigma = sig, n = n)
}

# ---- interni ------------------------------------------------------------------------------------------------------

.dcg_check_args <- function(centroids, territories, n_vertices, clamp_factor, poi_precision_um) {
  need <- c("cell_id", "region_id", "x", "y", "cell_type")
  if (!is.data.frame(centroids) || !all(need %in% names(centroids)))
    stop("derive_cell_geometry(): centroids deve avere le colonne ", paste(need, collapse = ", "))
  if (anyDuplicated(centroids$cell_id)) stop("derive_cell_geometry(): cell_id duplicati")
  if (!is.list(territories) || is.null(territories$cell_territories) || is.null(territories$territory_df))
    stop("derive_cell_geometry(): territories deve essere l'output di tessellate_voronoi()")
  cen <- centroids[order(centroids$cell_id), , drop = FALSE]
  td <- territories$territory_df
  if (length(territories$cell_territories) != nrow(cen) || nrow(td) != nrow(cen) ||
      !identical(as.numeric(sort(td$cell_id)), as.numeric(cen$cell_id)) ||
      !identical(as.numeric(td$cell_id), as.numeric(sort(td$cell_id))))
    stop("derive_cell_geometry(): territori e centroidi non corrispondono per cell_id (ordine crescente richiesto)")
  if (any(sf::st_is_empty(territories$cell_territories))) stop("derive_cell_geometry(): territori vuoti")
  if (!is.numeric(n_vertices) || length(n_vertices) != 1 || n_vertices < 8 || n_vertices != round(n_vertices))
    stop("derive_cell_geometry(): n_vertices deve essere un intero >= 8")
  if (!is.numeric(clamp_factor) || length(clamp_factor) != 1 || clamp_factor <= 0 || clamp_factor >= 1)
    stop("derive_cell_geometry(): clamp_factor deve stare in (0, 1)")
  if (!is.numeric(poi_precision_um) || poi_precision_um <= 0) stop("derive_cell_geometry(): poi_precision_um > 0")
  cen
}

.dcg_check_types <- function(cell_types, used) {
  if (!is.data.frame(cell_types) || !all(c("cell_type", "density") %in% names(cell_types)))
    stop("derive_cell_geometry(): cell_types deve avere almeno cell_type e density")
  ct <- cell_types
  miss <- setdiff(used, ct$cell_type)
  if (length(miss)) stop("derive_cell_geometry(): tipi senza catalogo: ", paste(miss, collapse = ", "))
  dflt <- list(nucleus_to_eq_ratio = NA_real_, nucleus_radius = NA_real_, nucleus_area_slope = 1,
               nucleus_area_sdlog = 0, nucleus_aspect = 1, nucleus_orient = "random")
  for (nm in names(dflt)) if (!nm %in% names(ct)) ct[[nm]] <- dflt[[nm]]
  ct$nucleus_orient <- as.character(ct$nucleus_orient); ct$nucleus_orient[is.na(ct$nucleus_orient)] <- "random"
  if ("nucleus_shape" %in% names(ct)) {
    bad <- (ct$nucleus_shape == "circle" & ct$nucleus_aspect != 1) | (ct$nucleus_shape == "ellipse" & ct$nucleus_aspect <= 1)
    if (any(bad, na.rm = TRUE)) stop("derive_cell_geometry(): nucleus_shape incoerente con nucleus_aspect per ",
                                     paste(ct$cell_type[which(bad)], collapse = ", "))
  }
  u <- ct[ct$cell_type %in% used, , drop = FALSE]
  okr <- !is.na(u$nucleus_radius) & u$nucleus_radius > 0
  okq <- !is.na(u$nucleus_to_eq_ratio) & u$nucleus_to_eq_ratio > 0 & !is.na(u$density) & u$density > 0
  if (any(!okr & !okq)) stop("derive_cell_geometry(): servono nucleus_radius > 0 oppure nucleus_to_eq_ratio > 0 e density > 0 per ",
                             paste(u$cell_type[!okr & !okq], collapse = ", "))
  if (any(!is.finite(u$nucleus_area_slope)) || any(!is.finite(u$nucleus_area_sdlog) | u$nucleus_area_sdlog < 0))
    stop("derive_cell_geometry(): nucleus_area_slope finito e nucleus_area_sdlog >= 0")
  if (any(!is.finite(u$nucleus_aspect) | u$nucleus_aspect < 1)) stop("derive_cell_geometry(): nucleus_aspect >= 1")
  om <- u$nucleus_orient
  okm <- om %in% c("random", "territory") | !is.na(suppressWarnings(as.numeric(om)))
  if (any(!okm)) stop("derive_cell_geometry(): nucleus_orient non valido: ", paste(om[!okm], collapse = ", "))
  ct
}

.dcg_center <- function(terr, cen, mode, prec) {
  if (mode == "generator") return(list(x = cen$x, y = cen$y))
  p <- polylabelr::poi(terr, precision = prec)
  list(x = vapply(p, function(q) q$x, 0), y = vapply(p, function(q) q$y, 0))
}

.dcg_target_area <- function(a_t, ct, ti, z) {
  rho <- ct$nucleus_to_eq_ratio[ti]; r_abs <- ct$nucleus_radius[ti]
  beta <- ct$nucleus_area_slope[ti]; sig <- ct$nucleus_area_sdlog[ti]
  a_ref <- 1e6 / ct$density[ti]
  abs_mode <- !is.na(r_abs) & r_abs > 0
  la <- ifelse(abs_mode, log(pi * r_abs^2), log(rho^2 * a_ref) + beta * (log(a_t) - log(a_ref)))
  exp(la + sig * z)
}

.dcg_poly_factor <- function(nv) sqrt(2 * pi / (nv * sin(2 * pi / nv)))

.dcg_orientation <- function(terr, ct, ti, u, asp) {
  om <- ct$nucleus_orient[ti]
  th <- pi * u - pi / 2
  num <- suppressWarnings(as.numeric(om)); fx <- !is.na(num); th[fx] <- num[fx]
  tt <- om == "territory" & asp > 1
  if (any(tt)) th[tt] <- .dcg_territory_axis(terr[tt])
  th[asp == 1] <- 0
  th
}

# asse maggiore del territorio dai momenti del secondo ordine (anello esterno meno buchi, tutte le parti)
.dcg_territory_axis <- function(g) {
  co <- sf::st_coordinates(g)
  lv <- setdiff(colnames(co), c("X", "Y"))
  fid <- co[, lv[length(lv)]]
  ring <- interaction(as.data.frame(co[, lv, drop = FALSE]), drop = TRUE, lex.order = TRUE)
  nxt <- c(seq_len(nrow(co))[-1], 1L)
  same <- c(ring[-1] == ring[-length(ring)], FALSE)
  x0 <- stats::ave(co[, "X"], fid); y0 <- stats::ave(co[, "Y"], fid)
  x <- co[, "X"] - x0; y <- co[, "Y"] - y0; x2 <- x[nxt]; y2 <- y[nxt]
  cr <- ifelse(same, x * y2 - x2 * y, 0)
  ra <- stats::ave(cr, ring, FUN = sum)                 # 2 x area con segno di ogni anello
  hole <- co[, "L1"] > 1                                 # L1 = anello nel poligono (1 = esterno)
  cr <- cr * ifelse(hole, -1, 1) * sign(ra)              # esterno con area positiva, buchi negativa
  A <- rowsum(cr, fid)[, 1] / 2
  sx <- rowsum((x + x2) * cr, fid)[, 1] / 6; sy <- rowsum((y + y2) * cr, fid)[, 1] / 6
  ixx <- rowsum((x^2 + x * x2 + x2^2) * cr, fid)[, 1] / 12
  iyy <- rowsum((y^2 + y * y2 + y2^2) * cr, fid)[, 1] / 12
  ixy <- rowsum((x * y2 + 2 * x * y + 2 * x2 * y2 + x2 * y) * cr, fid)[, 1] / 24
  cx <- sx / A; cy <- sy / A
  mxx <- ixx / A - cx^2; myy <- iyy / A - cy^2; mxy <- ixy / A - cx * cy
  as.numeric(0.5 * atan2(2 * mxy, mxx - myy))
}

.dcg_polygons <- function(cx, cy, a, b, theta, nv) {
  n <- length(cx); t <- 2 * pi * (seq_len(nv) - 1) / nv
  ct <- cos(t); st <- sin(t)
  lx <- outer(a, ct); ly <- outer(b, st)
  co <- cos(theta); si <- sin(theta)
  X <- cx + lx * co - ly * si; Y <- cy + lx * si + ly * co
  geos::geos_make_polygon(as.vector(t(X)), as.vector(t(Y)), feature_id = rep(seq_len(n), each = nv))
}

.dcg_clamp <- function(a_eq, b_eq, theta, ctr, d, gt, kp, nv, cf, n_iter = 40L) {
  s <- pmin(1, cf * d / (kp * a_eq))                     # cerchio circoscritto: sempre ammissibile
  e <- which(s < 1 & a_eq > b_eq)                         # ellissi clampate: si prova ad allargare
  if (length(e)) {
    fits <- function(sc, idx) geos::geos_within(
      .dcg_polygons(ctr$x[idx], ctr$y[idx], kp * a_eq[idx] * sc / cf, kp * b_eq[idx] * sc / cf, theta[idx], nv), gt[idx])
    f1 <- fits(rep(1, length(e)), e)
    s[e[f1]] <- 1
    e <- e[!f1]
    if (length(e)) {
      lo <- s[e]; hi <- rep(1, length(e))
      for (it in seq_len(n_iter)) {
        mid <- (lo + hi) / 2; ok <- fits(mid, e)
        lo[ok] <- mid[ok]; hi[!ok] <- mid[!ok]
      }
      s[e] <- lo
    }
  }
  s
}
