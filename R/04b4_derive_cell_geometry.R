# R/04b4_derive_cell_geometry.R — Step 1, fase 4: nucleo (centro, taglia, forma) per ogni territorio.
# Design: docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md §4.4, §5.4 (v1.5, S1.4d: modello "clip" = K-2 di
# S1.4c, decisione 2026-10-16; v1.4, S1.4: D-S1.4.1, D-S1.4.2 per il modello "scale").
# Pre-registrazioni: results/S1.4/S1.4_preregistration.md (c72c01a); docs/steps/S1.4d.md (5118bce, eae40f6).
# Convenzioni: µm, y verso il basso (come extract_regions()); angoli in rad da +x verso +y, in (-pi/2, pi/2].
# Dipendenze: sf, geos, polylabelr.

#' Geometria del nucleo per ogni cellula
#'
#' Modello di taglia (D-S1.4.1), per tipo t con `A_ref = 1e6 / density_t`:
#' `log A_N = log(rho^2 * A_ref) + beta * (log A_T - log A_ref) + sigma * z`, z ~ N(0, 1);
#' con `nucleus_radius` assoluto: `log A_N = log(pi * r^2) + sigma * z`. Con beta = 1, sigma = 0 il modello
#' coincide con il design (`nucleus_radius = rho * eq_radius`). L'area bersaglio non dipende da `nucleus_model`.
#'
#' Forma: poligono di `n_vertices` vertici su un'ellisse di semiassi a = sqrt(A * asp / pi), b = sqrt(A / (asp * pi)),
#' moltiplicati per k = sqrt(2 pi / (n sin(2 pi / n))) cosi' che l'area del poligono sia esattamente A.
#'
#' `nucleus_model = "clip"` (predefinito da S1.4d; K-2 di S1.4c): centro = generatore; tetto
#' `A* = min(A_N, phi * min(A_T, area(T_in)))` con phi = `nucleus_nc_max` del tipo e T_in = territorio rientrato di
#' `clip_inset_um` (calcolato senza modello di precisione, poi sulla griglia `snap_grid_um`); nucleo = poligono(s) ∩ T_in
#' (intersezione sulla griglia), con la scala s >= 1 scelta perche' area = A* (|area - A*| / A* <= `area_rel_tol`).
#' Solo le celle ritagliate a s = 1 cercano s: primo passo s = sqrt(A* / area(1)), poi falsa posizione di Illinois
#' nell'intervallo, con ricaduta sulla bisezione; arresto anche quando lo spostamento dei vertici fra gli estremi e' sotto
#' un quarto della griglia (euristica: `area_rel_tol` e' un obiettivo, non una garanzia; S1.4d). Nuclei contenuti nel proprio territorio e disgiunti per costruzione.
#'
#' `nucleus_model = "scale"` (D-S1.4.2, il modello fino a S1.4c, invariato): fattore di scala uniforme s <= 1. Cerchio:
#' raggio circoscritto k * a * s <= clamp_factor * d (d = distanza del centro dal bordo del territorio). Ellisse: s parte
#' dal valore del cerchio circoscritto e cresce per bisezione finche' il poligono dilatato di 1/clamp_factor resta dentro.
#'
#' @param centroids data.frame con `cell_id` (unico), `region_id`, `x`, `y`, `cell_type`; `cluster_id` opzionale.
#' @param territories output di `tessellate_voronoi()`: `cell_territories` (sfc, ordine di `cell_id` crescente) e
#'   `territory_df` (`cell_id`, `territory_area`).
#' @param cell_types catalogo: `cell_type`, `density` (cell/mm²), `nucleus_to_eq_ratio` (rho) e/o `nucleus_radius`
#'   (µm, NA = usa rho); opzionali `nucleus_area_slope` (beta, default 1), `nucleus_area_sdlog` (sigma, default 0),
#'   `nucleus_aspect` (>= 1, default 1 = cerchio), `nucleus_orient` (`"random"` default, `"territory"` = asse maggiore
#'   del territorio, oppure un angolo in rad, numerico o stringa numerica), `nucleus_nc_max` (phi in (0, 1], default 1;
#'   solo "clip"). S1.4f, solo con `nucleus_orient = "territory"`: `nucleus_orient_conc` (lambda >= 0, default Inf = esattamente
#'   lungo il territorio) e `nucleus_orient_law` (`"constant"` default, `"elongation"`, `"compression"`): angolo = asse del
#'   territorio + delta, 2 delta ~ Cauchy avvolta di lunghezza risultante rho_i = 1 - exp(-lambda w_i), estratta dalla stessa u
#'   dell'orientamento a caso (quantile in forma chiusa); w_i = 1, allungamento del territorio - 1 (radice del rapporto dei
#'   momenti principali), oppure raggio equivalente dell'area bersaglio / distanza del centro dal bordo del territorio.
#' @param random_seed seme; un solo flusso: `z <- rnorm(n)` e poi `u <- runif(n)`, in ordine di `cell_id`.
#' @param nucleus_model `"clip"` (ritaglio nel territorio a pari area) o `"scale"` (scala con clamp).
#' @param nucleus_center `"generator"` o `"poi"` (polo di inaccessibilita'); NULL = `"generator"` con "clip", `"poi"`
#'   con "scale".
#' @param n_vertices vertici del poligono del nucleo (>= 8).
#' @param clamp_factor frazione della distanza dal bordo concessa al nucleo (solo "scale", design: 0.85), in (0, 1).
#' @param clip_inset_um rientro del territorio per il ritaglio (µm, > 0; S1.4c D-1, BL-118a).
#' @param snap_grid_um griglia delle intersezioni (µm, > 0; F1 di S1.3d / V2 di S1.3b).
#' @param area_rel_tol tolleranza relativa obiettivo sull'area del nucleo ritagliato (> 0, < 1e-3); misurato in S1.4d: 94.1 % delle
#'   celle reali entro 1e-7, tutte entro 7.4e-7; per nuclei sotto il µm² il rumore della griglia supera 1e-6.
#' @param poi_precision_um precisione di `polylabelr::poi()` (µm).
#' @param resid_table S1.5d: NULL (predefinito, residui gaussiani sigma * z) o tabella dei residui veri (`key`, `decile`, `x_lo`, `x_hi`,
#'   `p`, `r`; `tools/S1.5d_size_resid_table.py`); i tipi con `nucleus_area_resid` = key usano r = quantile empirico in pnorm(z) del
#'   decimo di log A_T (assoluto) della cella al posto di sigma * z; stessa z, quindi stessa estrazione casuale del modello gaussiano.
#' @param resid_mode `"decile"` (predefinito) o `"pooled"` (residui di tutti i decimi insieme: controprova CP-3 di S1.5d).
#'   Colonna opzionale del catalogo `nucleus_area_resid_shift`: spostamento di log A_T per la ricerca del decimo (revisione RC-03:
#'   nel modello ricentrato log(m_t / m_sim)); assente o NA = 0.
#' @param strict se TRUE (default) errore se un nucleo non e' dentro il proprio territorio (C3).
#' @param verbose stampa un riepilogo.
#' @return lista con `cell_df`, `cell_nuclei` (sfc, ordine di `cell_id`; sfc_POLYGON, oppure sfc_MULTIPOLYGON se
#'   qualche nucleo ritagliato e' in piu' parti), `info`. Colonne di `cell_df`: schema del design §4.2 piu'
#'   `nucleus_area_target`, `nucleus_scale`, `center_dist`; con "scale" anche `nucleus_clamped`; con "clip" anche
#'   `nucleus_area_cap` (A*), `nucleus_capped` (A* < A_N), `nucleus_clipped` (ritagliato a s = 1), `nucleus_n_iter`.
derive_cell_geometry <- function(centroids, territories, cell_types, random_seed = 42,
                                 nucleus_model = c("clip", "scale"), nucleus_center = NULL, n_vertices = 32L,
                                 clamp_factor = 0.85, clip_inset_um = 1e-5, snap_grid_um = 1e-6, area_rel_tol = 1e-7,
                                 poi_precision_um = 1e-3, strict = TRUE, verbose = TRUE, resid_table = NULL,
                                 resid_mode = c("decile", "pooled")) {
  resid_mode <- match.arg(resid_mode)
  t0 <- proc.time()[["elapsed"]]
  nucleus_model <- match.arg(nucleus_model)
  if (is.null(nucleus_center)) nucleus_center <- if (nucleus_model == "clip") "generator" else "poi"
  nucleus_center <- match.arg(nucleus_center, c("poi", "generator"))
  cen <- .dcg_check_args(centroids, territories, n_vertices, clamp_factor, poi_precision_um)
  .dcg_check_clip_args(clip_inset_um, snap_grid_um, area_rel_tol)
  ct <- .dcg_check_types(cell_types, unique(cen$cell_type))
  terr <- .dcg_as_polys(territories$cell_territories)
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
  if (any(d <= 0)) stop("derive_cell_geometry(): centro del nucleo sul bordo del territorio (d = 0) per ", sum(d <= 0), " cellule")
  a_tgt <- .dcg_target_area(a_t, ct, ti, z, resid_table, resid_mode)
  asp <- ct$nucleus_aspect[ti]
  theta <- .dcg_orientation(terr, ct, ti, u, asp, d = d, a_tgt = a_tgt)
  kp <- .dcg_poly_factor(n_vertices)
  if (nucleus_model == "clip")
    return(.dcg_run_clip(cen, ct, ti, a_t, a_tgt, asp, theta, ctr, d, gt, kp, n_vertices, clip_inset_um, snap_grid_um,
                         area_rel_tol, strict, verbose, t0,
                         list(nucleus_center = nucleus_center, random_seed = random_seed, poi_precision_um = poi_precision_um)))
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
               elapsed_s = proc.time()[["elapsed"]] - t0, nucleus_model = "scale")
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
               nucleus_area_sdlog = 0, nucleus_aspect = 1, nucleus_orient = "random", nucleus_nc_max = 1, nucleus_area_resid = NA_character_,
               nucleus_orient_conc = Inf, nucleus_orient_law = "constant")
  for (nm in names(dflt)) if (!nm %in% names(ct)) ct[[nm]] <- dflt[[nm]]
  ct$nucleus_orient <- as.character(ct$nucleus_orient); ct$nucleus_orient[is.na(ct$nucleus_orient)] <- "random"
  if (!is.numeric(ct$nucleus_orient_conc)) stop("derive_cell_geometry(): nucleus_orient_conc deve essere numerico (lambda >= 0, Inf)")
  ct$nucleus_orient_conc[is.na(ct$nucleus_orient_conc)] <- Inf
  ct$nucleus_orient_law <- as.character(ct$nucleus_orient_law); ct$nucleus_orient_law[is.na(ct$nucleus_orient_law)] <- "constant"
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
  if (any(!is.finite(u$nucleus_nc_max) | u$nucleus_nc_max <= 0 | u$nucleus_nc_max > 1))
    stop("derive_cell_geometry(): nucleus_nc_max (phi) deve stare in (0, 1]")
  om <- u$nucleus_orient
  okm <- om %in% c("random", "territory") | !is.na(suppressWarnings(as.numeric(om)))
  if (any(!okm)) stop("derive_cell_geometry(): nucleus_orient non valido: ", paste(om[!okm], collapse = ", "))
  if (any(u$nucleus_orient_conc < 0)) stop("derive_cell_geometry(): nucleus_orient_conc (lambda) deve essere >= 0")
  okl <- u$nucleus_orient_law %in% c("constant", "elongation", "compression")
  if (any(!okl)) stop("derive_cell_geometry(): nucleus_orient_law non valida: ", paste(u$nucleus_orient_law[!okl], collapse = ", "))
  part <- is.finite(u$nucleus_orient_conc) | u$nucleus_orient_law != "constant"
  if (any(part & om != "territory"))
    stop("derive_cell_geometry(): nucleus_orient_conc / nucleus_orient_law valgono solo con nucleus_orient = \"territory\" (tipi ",
         paste(u$cell_type[part & om != "territory"], collapse = ", "), ")")
  ct
}

# sfc misto (POLYGON + MULTIPOLYGON, orfani senza confine di D-S1.3.2) -> MULTIPOLYGON: polylabelr e st_coordinates
# non accettano sfc_GEOMETRY (difetto D-1 di S1.4, trovato sui ROI reali)
.dcg_as_polys <- function(g) {
  if (inherits(g, c("sfc_POLYGON", "sfc_MULTIPOLYGON"))) return(g)
  ty <- unique(as.character(sf::st_geometry_type(g)))
  if (!all(ty %in% c("POLYGON", "MULTIPOLYGON"))) stop("derive_cell_geometry(): territori non poligonali: ", paste(ty, collapse = ", "))
  sf::st_cast(g, "MULTIPOLYGON")
}

.dcg_center <- function(terr, cen, mode, prec) {
  if (mode == "generator") return(list(x = cen$x, y = cen$y))
  p <- polylabelr::poi(terr, precision = prec)
  list(x = vapply(p, function(q) q$x, 0), y = vapply(p, function(q) q$y, 0))
}

.dcg_target_area <- function(a_t, ct, ti, z, resid_table = NULL, resid_mode = "decile") {
  rho <- ct$nucleus_to_eq_ratio[ti]; r_abs <- ct$nucleus_radius[ti]
  beta <- ct$nucleus_area_slope[ti]; sig <- ct$nucleus_area_sdlog[ti]
  a_ref <- 1e6 / ct$density[ti]
  abs_mode <- !is.na(r_abs) & r_abs > 0
  la <- ifelse(abs_mode, log(pi * r_abs^2), log(rho^2 * a_ref) + beta * (log(a_t) - log(a_ref)))
  eps <- sig * z
  key <- if (!is.null(resid_table) && "nucleus_area_resid" %in% names(ct)) as.character(ct$nucleus_area_resid[ti]) else rep(NA_character_, length(ti))
  use <- !is.na(key) & !abs_mode
  sh <- if ("nucleus_area_resid_shift" %in% names(ct)) as.numeric(ct$nucleus_area_resid_shift[ti]) else rep(0, length(ti))   # revisione RC-03
  sh[is.na(sh)] <- 0
  if (any(use)) eps[use] <- .dcg_resid_draw(resid_table, key[use], log(a_t[use]) + sh[use], stats::pnorm(z[use]), resid_mode)
  exp(la + eps)
}

#' Residuo empirico (S1.5d): quantile della tabella del gruppo nel decimo di log A_T, interpolato linearmente in u.
#' @keywords internal
.dcg_resid_draw <- function(tab, key, lx, u, mode) {
  out <- numeric(length(key))
  for (k in unique(key)) {
    w <- which(key == k); tk <- tab[tab$key == k, , drop = FALSE]
    if (!nrow(tk)) stop("derive_cell_geometry(): tabella dei residui senza la chiave ", k, call. = FALSE)
    tk <- if (mode == "pooled") tk[tk$decile == 0L, , drop = FALSE] else tk[tk$decile > 0L, , drop = FALSE]
    tk <- tk[order(tk$decile, tk$p), , drop = FALSE]
    dec <- sort(unique(tk$decile)); np <- sum(tk$decile == dec[1])
    Q <- matrix(tk$r, nrow = length(dec), ncol = np, byrow = TRUE)
    if (mode == "pooled") di <- rep(1L, length(w)) else {
      hi <- tk$x_hi[match(dec, tk$decile)]                                          # limiti superiori dei decimi (log A_T assoluto)
      di <- findInterval(lx[w], hi[-length(hi)], left.open = TRUE) + 1L             # fuori dall'intervallo vero: decimo estremo
    }
    pos <- u[w] * (np - 1); j <- pmin(floor(pos), np - 2); fr <- pos - j
    out[w] <- Q[cbind(di, j + 1L)] + fr * (Q[cbind(di, j + 2L)] - Q[cbind(di, j + 1L)])
  }
  out
}

.dcg_poly_factor <- function(nv) sqrt(2 * pi / (nv * sin(2 * pi / nv)))

.dcg_orientation <- function(terr, ct, ti, u, asp, d = NULL, a_tgt = NULL) {
  om <- ct$nucleus_orient[ti]
  th <- pi * u - pi / 2
  num <- suppressWarnings(as.numeric(om)); fx <- !is.na(num); th[fx] <- num[fx]
  tt <- om == "territory" & asp > 1
  if (any(tt)) {
    it <- which(tt); lam <- ct$nucleus_orient_conc[ti][it]; pe <- is.finite(lam)
    if (!any(pe)) th[it] <- .dcg_territory_axis(terr[it]) else {          # S1.4f: orientamento parziale
      mo <- .dcg_territory_moments(terr[it]); th[it] <- mo$theta
      ip <- it[pe]; w <- .dcg_orient_weight(ct$nucleus_orient_law[ti][ip], mo$aspect[pe], a_tgt[ip], d[ip])
      rho <- 1 - exp(-lam[pe] * w)
      th[ip] <- th[ip] + 0.5 * .dcg_wc_quantile(u[ip], rho)
    }
  }
  th[asp == 1] <- 0
  th
}

#' Peso w_i della concentrazione (S1.4f): "constant" 1; "elongation" allungamento del territorio - 1; "compression" raggio
#' equivalente dell'area bersaglio / distanza del centro dal bordo del territorio.
#' @keywords internal
.dcg_orient_weight <- function(law, aspect_t, a_tgt, d) {
  w <- rep(1, length(law))
  el <- law == "elongation"; w[el] <- aspect_t[el] - 1
  co <- law == "compression"; w[co] <- sqrt(a_tgt[co] / pi) / d[co]
  w
}

#' Quantile della Cauchy avvolta su (-pi, pi) con centro 0 e lunghezza risultante media rho in [0, 1] (S1.4f):
#' Q(u) = 2 atan((1 - rho) / (1 + rho) tan(pi (u - 1/2))); rho = 0 uniforme, rho = 1 punto in 0.
#' @keywords internal
.dcg_wc_quantile <- function(u, rho) 2 * atan((1 - rho) / (1 + rho) * tan(pi * (u - 0.5)))

# asse maggiore del territorio dai momenti del secondo ordine (anello esterno meno buchi, tutte le parti)
.dcg_territory_axis <- function(g) .dcg_territory_moments(g)$theta

# angolo dell'asse maggiore e allungamento (radice del rapporto degli autovalori) dai momenti del secondo ordine
.dcg_territory_moments <- function(g) {
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
  tq <- sqrt(((mxx - myy) / 2)^2 + mxy^2); l1 <- (mxx + myy) / 2 + tq; l2 <- pmax((mxx + myy) / 2 - tq, 1e-15)
  list(theta = as.numeric(0.5 * atan2(2 * mxy, mxx - myy)), aspect = as.numeric(sqrt(l1 / l2)))
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

# ---- modello "clip" (S1.4d; K-2 di S1.4c) --------------------------------------------------------------------------

.dcg_check_clip_args <- function(inset, grid, tol) {
  ok1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  if (!ok1(inset) || inset <= 0) stop("derive_cell_geometry(): clip_inset_um deve essere > 0")
  if (!ok1(grid) || grid <= 0) stop("derive_cell_geometry(): snap_grid_um deve essere > 0")
  if (!ok1(tol) || tol <= 0 || tol >= 1e-3) stop("derive_cell_geometry(): area_rel_tol deve stare in (0, 1e-3)")
  invisible(TRUE)
}

# territorio rientrato: il buffer negativo si calcola SENZA modello di precisione (con la griglia ereditata
# geos_buffer(-eps) restituisce in silenzio poligoni vuoti, BL-118b), poi si torna sulla griglia; vuoti = errore
.dcg_inset <- function(gt, inset, grid) {
  tg <- geos::geos_set_precision(gt, grid)
  tin <- geos::geos_set_precision(geos::geos_buffer(geos::geos_set_precision(tg, 0), -inset,
                                                    params = geos::geos_buffer_params(quad_segs = 2L)), grid)
  emp <- geos::geos_is_empty(tin)
  if (any(emp)) stop(sprintf("derive_cell_geometry(): territorio rientrato vuoto per %d celle (BL-118)", sum(emp)))
  tin
}

# scala s per cella: area(poligono(s) ∩ tin) = a_star. Restituisce nuclei (geos), aree, s, iterazioni, ritagliati.
.dcg_clip_solve <- function(cx, cy, a0, b0, theta, kp, nv, tin, a_star, grid, tol, max_iter = 100L) {
  n <- length(cx); s <- rep(1, n); it <- integer(n)
  clip_at <- function(idx, sc) geos::geos_intersection_prec(
    .dcg_polygons(cx[idx], cy[idx], kp * a0[idx] * sc, kp * b0[idx] * sc, theta[idx], nv), tin[idx], grid)
  nuc <- clip_at(seq_len(n), s)
  area <- geos::geos_area(nuc)
  need <- which(area < a_star * (1 - tol))
  if (length(need)) {
    A <- a_star[need]; m <- length(need); rv <- kp * pmax(a0[need], b0[need])   # raggio del vertice piu' lontano a s = 1
    lo <- rep(1, m); alo <- area[need]
    ext <- geos::geos_extent(tin[need]); xs <- cx[need]; ys <- cy[need]
    dmax <- sqrt(pmax((xs - ext$xmin)^2, (xs - ext$xmax)^2) + pmax((ys - ext$ymin)^2, (ys - ext$ymax)^2))
    # il poligono di scala s contiene l'ellisse inscritta di semiasse minore k * b0 * s * cos(pi / nv)
    hi <- pmax(dmax / (kp * b0[need] * cos(pi / nv)) * (1 + 1e-9), 1 + 1e-9)
    ahi <- geos::geos_area(clip_at(need, hi))
    fl <- alo - A; fh <- ahi - A                     # valori per l'interpolazione (dimezzati da Illinois)
    side <- integer(m); wa <- rep(Inf, m); wb <- rep(Inf, m)
    done <- abs(ahi - A) / A <= tol                  # tetto = area del territorio rientrato: hi e' gia' la soluzione
    for (k in seq_len(max_iter)) {
      # arresto anche quando lo spostamento dei vertici fra gli estremi e' < griglia / 4 (nota D-1 di S1.4d). E' un'euristica:
      # i vertici d'intersezione scorrono lungo il lato del territorio piu' dello spostamento radiale quando il taglio e'
      # radente (revisione S1.4d, F-02), quindi area_rel_tol e' un obiettivo, non una garanzia: sui tassellamenti reali di
      # S1.4d il 94.1 % delle celle e' entro 1e-7 e tutte entro 7.4e-7; sotto il µm² il rumore della griglia supera 1e-6
      act <- which(!done & (hi - lo) * rv > grid / 4)
      if (!length(act)) break
      L <- lo[act]; H <- hi[act]; W <- H - L
      cand <- if (k == 1L) L * sqrt(A[act] / pmax(alo[act], .Machine$double.xmin))
              else (L * fh[act] - H * fl[act]) / (fh[act] - fl[act])
      bis <- !is.finite(cand) | cand <= L | cand >= H | W > 0.5 * wb[act]
      cand[bis] <- (L[bis] + H[bis]) / 2
      am <- geos::geos_area(clip_at(need[act], cand))
      it[need[act]] <- it[need[act]] + 1L
      fm <- am - A[act]; up <- fm < 0
      iu <- act[up]; id <- act[!up]
      lo[iu] <- cand[up]; alo[iu] <- am[up]; fl[iu] <- fm[up]
      fh[iu] <- ifelse(side[iu] == -1L, fh[iu] / 2, fh[iu]); side[iu] <- -1L
      hi[id] <- cand[!up]; ahi[id] <- am[!up]; fh[id] <- fm[!up]
      fl[id] <- ifelse(side[id] == 1L, fl[id] / 2, fl[id]); side[id] <- 1L
      wb[act] <- wa[act]; wa[act] <- W
      done[act] <- abs(fm) / A[act] <= tol
    }
    left <- sum(!done & (hi - lo) * rv > grid / 4)
    if (left) warning(sprintf("derive_cell_geometry(): %d celle non convergenti dopo %d iterazioni", left, max_iter))
    best <- ifelse(abs(alo - A) <= abs(ahi - A), lo, hi)
    s[need] <- best
    nuc[need] <- clip_at(need, best)
    area[need] <- geos::geos_area(nuc[need])
  }
  list(nuc = nuc, area = area, s = s, n_iter = it, clipped = seq_len(n) %in% need)
}

# geos -> sf via WKB: con la griglia 1e-6 ereditata st_as_sfc() legge la precisione come fattore di scala (nota S1.4b)
.dcg_to_sfc <- function(g) {
  out <- sf::st_as_sfc(geos::geos_write_wkb(g))
  ty <- as.character(sf::st_geometry_type(out))
  if (all(ty == "POLYGON")) sf::st_cast(out, "POLYGON") else sf::st_cast(out, "MULTIPOLYGON")
}

.dcg_run_clip <- function(cen, ct, ti, a_t, a_tgt, asp, theta, ctr, d, gt, kp, nv, inset, grid, tol, strict, verbose,
                          t0, meta) {
  n <- nrow(cen)
  tin <- .dcg_inset(gt, inset, grid)
  phi <- ct$nucleus_nc_max[ti]
  a_star <- pmin(a_tgt, phi * pmin(a_t, geos::geos_area(tin)))
  a0 <- sqrt(a_star * asp / pi); b0 <- sqrt(a_star / (asp * pi))
  r <- .dcg_clip_solve(ctr$x, ctr$y, a0, b0, theta, kp, nv, tin, a_star, grid, tol)
  emp <- geos::geos_is_empty(r$nuc)
  if (any(emp)) stop(sprintf("derive_cell_geometry(): %d nuclei vuoti", sum(emp)))
  within <- geos::geos_within(r$nuc, gt)
  if (strict && !all(within)) stop(sprintf("derive_cell_geometry(): %d nuclei non contenuti nel territorio (C3)", sum(!within)))
  multi <- geos::geos_num_geometries(r$nuc) > 1L
  circ <- asp == 1
  cell_df <- data.frame(
    cell_id = cen$cell_id, region_id = cen$region_id,
    cluster_id = if ("cluster_id" %in% names(cen)) cen$cluster_id else NA_integer_,
    cell_type = cen$cell_type, x = cen$x, y = cen$y,
    territory_area = a_t, eq_radius = sqrt(a_t / pi),
    nucleus_shape = ifelse(circ, "circle", "ellipse"),
    nucleus_radius = ifelse(circ, sqrt(r$area / pi), NA_real_),
    nucleus_major_r = ifelse(circ, NA_real_, a0 * r$s), nucleus_minor_r = ifelse(circ, NA_real_, b0 * r$s),
    nucleus_orient = ifelse(circ, NA_real_, theta),
    nucleus_cx = ctr$x, nucleus_cy = ctr$y, nucleus_area = r$area, cytoplasm_area = a_t - r$area,
    nucleus_area_target = a_tgt, nucleus_area_cap = a_star, nucleus_capped = a_star < a_tgt,
    nucleus_clipped = r$clipped, nucleus_scale = r$s, nucleus_n_iter = r$n_iter, center_dist = d,
    stringsAsFactors = FALSE)
  rel <- abs(r$area - a_star) / a_star
  info <- c(list(n = n, nucleus_model = "clip", n_clipped = sum(r$clipped), n_capped = sum(a_star < a_tgt),
                 frac_capped = mean(a_star < a_tgt), n_multipart = sum(multi), n_not_within = sum(!within),
                 max_iter = max(c(0L, r$n_iter)), total_iter = sum(r$n_iter), max_area_rel_err = max(rel),
                 n_vertices = as.integer(nv), poly_factor = kp, clip_inset_um = inset, snap_grid_um = grid,
                 area_rel_tol = tol), meta, list(elapsed_s = proc.time()[["elapsed"]] - t0))
  if (verbose) message(sprintf("derive_cell_geometry(): %d cellule, clip: %d ritagliate, %d al tetto (%.2f %%), max %d iterazioni, %.1f s",
                               n, info$n_clipped, info$n_capped, 100 * info$frac_capped, info$max_iter, info$elapsed_s))
  list(cell_df = cell_df, cell_nuclei = .dcg_to_sfc(r$nuc), info = info)
}
