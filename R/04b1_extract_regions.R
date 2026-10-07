#' Estrae le regioni connesse da una mappa di cluster in pixel (Step 1, S1.1)
#'
#' Trasforma l'output di `cluster_image()` (un pixel per riga, con il cluster
#' di appartenenza) in **regioni connesse** (4-vicinato) per cluster, restituite
#' come poligoni `sf` in µm. Il contorno segue i **lati dei pixel** (non i
#' centri): l'area del poligono non semplificato coincide esattamente con
#' n_px × (s × pixel_size_um)^2. Componenti con area < `min_region_area_um2`
#' sono escluse (spazio extracellulare macro per lo Step 2) e riportate in
#' `excluded_df`. Design: `docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md`
#' §5.1, con le deviazioni dichiarate in `results/S1.1/S1.1_preregistration.md`.
#'
#' Convenzione spaziale: il pixel (x, y) copre
#' [(x-1)*px, (x-1+s)*px] x [(y-1)*px, (y-1+s)*px], con y verso il basso
#' (righe dell'immagine), px = `pixel_size_um`, s = stride.
#'
#' @param clust data.frame con colonne intere `x`, `y` e la colonna dei cluster.
#'   I pixel assenti sono fondo (sotto soglia).
#' @param pixel_size_um µm per pixel. **Obbligatorio** (decisione 10 / C-A).
#' @param min_region_area_um2 area minima (µm²) di una regione; filtro sull'area
#'   esatta in pixel. Default 100 (design, senza fonte biologica).
#' @param cluster_col nome della colonna dei cluster (factor, intero o carattere).
#' @param simplify_tol_um tolleranza di `sf::st_simplify()` (preserveTopology).
#'   0 = solo rimozione dei vertici collineari (nessuna perdita di area).
#' @param stride passo di campionamento dei pixel; `NULL` = rilevato dai dati
#'   (MCD delle differenze fra coordinate distinte).
#' @param verbose stampa un riepilogo.
#' @return lista con `region_df`, `region_polygons` (sfc, µm, senza CRS),
#'   `excluded_df`, `info`.
#' @export
extract_regions <- function(clust,
                            pixel_size_um,
                            min_region_area_um2 = 100,
                            cluster_col = "intensity_cluster",
                            simplify_tol_um = 0.5,
                            stride = NULL,
                            verbose = TRUE) {
  t0 <- proc.time()[["elapsed"]]                        # include la validazione
  if (missing(pixel_size_um)) {
    stop("extract_regions(): `pixel_size_um` e' obbligatorio (decisione 10 / C-A), nessun default.",
         call. = FALSE)
  }
  .er_check_args(clust, pixel_size_um, min_region_area_um2, cluster_col,
                 simplify_tol_um, stride)

  grid <- .er_build_grid(clust, cluster_col, stride)
  runs <- .er_label_runs(grid)
  cc   <- .er_run_components(runs)
  comp <- cc$df

  px    <- pixel_size_um
  e_um  <- grid$stride * px
  comp$area_px_um2 <- comp$n_px * e_um^2
  keep  <- comp$area_px_um2 >= min_region_area_um2

  polys <- .er_component_polygons(runs, cc$comp_of_run, which(keep), grid)
  polys <- polys * px                                  # unita' pixel -> µm
  if (simplify_tol_um > 0) {
    polys <- sf::st_simplify(polys, preserveTopology = TRUE,
                             dTolerance = simplify_tol_um)
  }

  kc <- comp[keep, , drop = FALSE]
  region_df <- data.frame(
    region_id    = seq_len(nrow(kc)),
    cluster_id   = grid$cluster_ids[kc$code],
    n_px         = kc$n_px,
    area_um2     = if (nrow(kc)) as.numeric(sf::st_area(polys)) else numeric(0),
    area_px_um2  = kc$area_px_um2,
    perimeter_um = if (nrow(kc)) as.numeric(sf::st_length(sf::st_boundary(polys))) else numeric(0),
    n_holes      = if (nrow(kc)) lengths(polys) - 1L else integer(0)
  )
  ec <- comp[!keep, , drop = FALSE]
  excluded_df <- data.frame(
    excluded_id = seq_len(nrow(ec)),
    cluster_id  = grid$cluster_ids[ec$code],
    n_px        = ec$n_px,
    area_px_um2 = ec$area_px_um2
  )

  area_tot <- sum(comp$area_px_um2)
  info <- list(
    pixel_size_um       = px,
    stride              = grid$stride,
    effective_pixel_um  = e_um,
    min_region_area_um2 = min_region_area_um2,
    simplify_tol_um     = simplify_tol_um,
    cluster_ids         = grid$cluster_ids,
    n_pixels            = grid$n_pixels,
    n_components        = nrow(comp),
    n_regions           = nrow(region_df),
    n_excluded          = nrow(excluded_df),
    area_total_um2      = area_tot,
    area_excluded_um2   = sum(ec$area_px_um2),
    frac_area_excluded  = if (area_tot > 0) sum(ec$area_px_um2) / area_tot else NA_real_,
    origin_px           = c(x = grid$x0, y = grid$y0),
    grid_dim            = c(nx = grid$nx, ny = grid$ny),
    elapsed_s           = proc.time()[["elapsed"]] - t0
  )
  if (verbose) {
    cat(sprintf("[extract_regions] %d pixel, stride %d, %d componenti -> %d regioni, %d escluse (%.2f%% dell'area), %.1f s\n",
                info$n_pixels, info$stride, info$n_components, info$n_regions,
                info$n_excluded, 100 * info$frac_area_excluded, info$elapsed_s))
  }
  list(region_df = region_df, region_polygons = polys,
       excluded_df = excluded_df, info = info)
}

# ---- helper interni ---------------------------------------------------------

.er_check_args <- function(clust, pixel_size_um, min_region_area_um2,
                           cluster_col, simplify_tol_um, stride) {
  if (!is.numeric(pixel_size_um) || length(pixel_size_um) != 1L ||
      !is.finite(pixel_size_um) || pixel_size_um <= 0) {
    stop("extract_regions(): `pixel_size_um` deve essere un numero finito > 0 (µm/px).",
         call. = FALSE)
  }
  if (!is.data.frame(clust)) stop("extract_regions(): `clust` deve essere un data.frame.", call. = FALSE)
  miss <- setdiff(c("x", "y", cluster_col), names(clust))
  if (length(miss)) {
    stop(sprintf("extract_regions(): colonne mancanti in `clust`: %s.",
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  if (nrow(clust) == 0L) stop("extract_regions(): `clust` non ha righe.", call. = FALSE)
  for (v in c("x", "y")) {
    z <- clust[[v]]
    if (!is.numeric(z) || anyNA(z) || any(!is.finite(z)) || any(abs(z - round(z)) > 1e-8)) {
      stop(sprintf("extract_regions(): `%s` deve contenere coordinate pixel intere, senza NA.", v),
           call. = FALSE)
    }
  }
  if (anyNA(clust[[cluster_col]])) {
    stop("extract_regions(): la colonna dei cluster contiene NA.", call. = FALSE)
  }
  # chiave numerica (x, y): anyDuplicated su un data.frame incolla 16 M stringhe
  # (24 s e ~5 GB su 4000x4000, misurato in S1.1)
  xr <- round(clust$x); yr <- round(clust$y)
  if (anyDuplicated((xr - min(xr)) * (max(yr) - min(yr) + 1) + (yr - min(yr)))) {
    stop("extract_regions(): coppie (x, y) duplicate in `clust`.", call. = FALSE)
  }
  if (!is.numeric(min_region_area_um2) || length(min_region_area_um2) != 1L || min_region_area_um2 < 0) {
    stop("extract_regions(): `min_region_area_um2` deve essere un numero >= 0.", call. = FALSE)
  }
  if (!is.numeric(simplify_tol_um) || length(simplify_tol_um) != 1L || simplify_tol_um < 0) {
    stop("extract_regions(): `simplify_tol_um` deve essere un numero >= 0.", call. = FALSE)
  }
  if (!is.null(stride) && (!is.numeric(stride) || length(stride) != 1L ||
                           stride < 1 || stride != round(stride))) {
    stop("extract_regions(): `stride` deve essere NULL o un intero >= 1.", call. = FALSE)
  }
  invisible(TRUE)
}

.er_gcd <- function(v) {
  g <- 0
  for (a in v) { b <- g; while (a > 0) { t <- b %% a; b <- a; a <- t }; g <- b }
  g
}

.er_axis_stride <- function(z) {
  u <- sort(unique(z))
  if (length(u) < 2L) return(NA_integer_)
  as.integer(.er_gcd(unique(diff(u))))
}

#' Matrice delle etichette (nx x ny; x varia piu' in fretta) e metadati di griglia
.er_build_grid <- function(clust, cluster_col, stride) {
  x <- as.integer(round(clust$x)); y <- as.integer(round(clust$y))
  cl <- clust[[cluster_col]]
  if (is.factor(cl)) cl <- as.character(cl)
  num <- suppressWarnings(as.numeric(cl))
  cl_val <- if (!anyNA(num)) num else cl
  cluster_ids <- sort(unique(cl_val))
  code <- match(cl_val, cluster_ids)

  if (is.null(stride)) {
    sx <- .er_axis_stride(x); sy <- .er_axis_stride(y)
    if (is.na(sx) && is.na(sy)) s <- 1L
    else if (is.na(sx)) s <- sy
    else if (is.na(sy)) s <- sx
    else if (sx != sy) {
      stop(sprintf("extract_regions(): stride diverso sui due assi (x = %d, y = %d).", sx, sy),
           call. = FALSE)
    } else s <- sx
  } else s <- as.integer(stride)
  x0 <- min(x); y0 <- min(y)
  if (any((x - x0) %% s != 0L) || any((y - y0) %% s != 0L)) {
    stop(sprintf("extract_regions(): coordinate non allineate allo stride %d.", s), call. = FALSE)
  }
  i <- (x - x0) %/% s + 1L; j <- (y - y0) %/% s + 1L
  nx <- max(i); ny <- max(j)
  lab <- matrix(0L, nrow = nx, ncol = ny)
  lab[cbind(i, j)] <- code
  list(lab = lab, nx = nx, ny = ny, stride = s, x0 = x0, y0 = y0,
       cluster_ids = cluster_ids, n_pixels = length(x))
}

#' Run orizzontali (stessa riga, stesso cluster, pixel contigui) e coppie di run
#' adiacenti in verticale (4-vicinato)
.er_label_runs <- function(grid) {
  v  <- as.vector(grid$lab); nx <- grid$nx; n <- length(v)
  pos <- seq_len(n)
  row_start <- (pos - 1L) %% nx == 0L
  row_end   <- pos %% nx == 0L
  nz <- v != 0L
  st <- nz & (c(TRUE, v[-1L] != v[-n]) | row_start)
  en <- nz & (c(v[-1L] != v[-n], TRUE) | row_end)
  ps <- which(st); pe <- which(en)
  stopifnot(length(ps) == length(pe))
  rid <- cumsum(st); rid[!nz] <- 0L
  rm(pos, row_start, row_end, st, en)
  if (grid$ny > 1L) {
    a <- rid[seq_len(n - nx)]; b <- rid[(nx + 1L):n]
    k <- a != 0L & v[seq_len(n - nx)] == v[(nx + 1L):n]
    a <- a[k]; b <- b[k]
    key <- unique(as.numeric(a) * (length(ps) + 1) + b)
    ea <- as.integer(key %/% (length(ps) + 1)); eb <- as.integer(key %% (length(ps) + 1))
  } else {
    ea <- integer(0); eb <- integer(0)
  }
  list(j  = (ps - 1L) %/% nx + 1L,
       i0 = (ps - 1L) %% nx + 1L,
       i1 = (pe - 1L) %% nx + 1L,
       code = v[ps],
       ea = ea, eb = eb)
}

#' Union-find vettorizzato (aggancio delle radici + compressione dei puntatori)
.er_union_find <- function(n, a, b) {
  lab <- seq_len(n)
  if (!length(a)) return(lab)
  repeat {
    la <- lab[a]; lb <- lab[b]
    d <- la != lb
    if (!any(d)) break
    hi <- pmax(la[d], lb[d]); lo <- pmin(la[d], lb[d])
    o <- order(hi, lo); hi <- hi[o]; lo <- lo[o]
    f <- !duplicated(hi)
    lab[hi[f]] <- pmin(lab[hi[f]], lo[f])
    repeat { nl <- lab[lab]; if (identical(nl, lab)) break; lab <- nl }
  }
  lab
}

#' Componenti connesse sulle run; ordine: cluster, poi primo pixel in ordine di riga.
#' Restituisce list(df = data.frame(comp_id, code, n_px), comp_of_run = int).
.er_run_components <- function(runs) {
  nr <- length(runs$j)
  root <- .er_union_find(nr, runs$ea, runs$eb)
  roots <- unique(root)
  roots <- roots[order(runs$code[roots], roots)]
  comp_of_run <- match(root, roots)
  len <- runs$i1 - runs$i0 + 1L
  n_px <- as.integer(rowsum(len, comp_of_run, reorder = TRUE)[, 1])
  list(df = data.frame(comp_id = seq_along(roots), code = runs$code[roots], n_px = n_px),
       comp_of_run = comp_of_run)
}

#' Rimuove i vertici collineari da un anello assiale chiuso (coordinate intere)
.er_drop_collinear <- function(m) {
  m <- m[-nrow(m), , drop = FALSE]                    # apre l'anello
  repeat {
    k <- nrow(m); if (k <= 4L) break
    pr <- m[c(k, seq_len(k - 1L)), , drop = FALSE]
    nx <- m[c(seq(2L, k), 1L), , drop = FALSE]
    col <- (pr[, 1] == m[, 1] & m[, 1] == nx[, 1]) | (pr[, 2] == m[, 2] & m[, 2] == nx[, 2])
    if (!any(col)) break
    m <- m[!col, , drop = FALSE]
  }
  rbind(m, m[1L, ])
}

#' Poligoni delle componenti selezionate, in unita' pixel originali (interi).
#' Rettangolo di una run: x in [x0 + (i0-1)s - 1, x0 + i1*s - 1], y analogo.
#' Una run sola -> rettangolo; piu' run -> unione GEOS dei rettangoli (non 'coverage': i lati condivisi
#' hanno giunzioni a T, che CoverageUnion non gestisce; vedi smoke test S1.1).
.er_component_polygons <- function(runs, comp_of_run, sel, grid) {
  if (!length(sel)) return(sf::st_sfc())
  s <- grid$stride
  xmin <- grid$x0 + (runs$i0 - 1L) * s - 1L; xmax <- grid$x0 + runs$i1 * s - 1L
  ymin <- grid$y0 + (runs$j  - 1L) * s - 1L; ymax <- grid$y0 + runs$j  * s - 1L
  rect <- function(r) {
    sf::st_polygon(list(matrix(c(xmin[r], ymin[r], xmax[r], ymin[r], xmax[r], ymax[r],
                                 xmin[r], ymax[r], xmin[r], ymin[r]), ncol = 2L, byrow = TRUE)))
  }
  runs_by_comp <- split(seq_along(comp_of_run), comp_of_run)[as.character(sel)]
  geoms <- lapply(runs_by_comp, function(r) {
    if (length(r) == 1L) return(rect(r))
    # run consecutive in verticale con lo stesso intervallo -> un solo rettangolo
    key <- paste(runs$i0[r], runs$i1[r])
    o <- order(key, runs$j[r]); r <- r[o]; key <- key[o]
    brk <- c(TRUE, key[-1L] != key[-length(key)] | diff(runs$j[r]) != 1L)
    grp <- cumsum(brk)
    rr <- lapply(split(r, grp), function(q) {
      a <- q[1L]; b <- q[length(q)]
      sf::st_polygon(list(matrix(c(xmin[a], ymin[a], xmax[a], ymin[a], xmax[a], ymax[b],
                                   xmin[a], ymax[b], xmin[a], ymin[a]), ncol = 2L, byrow = TRUE)))
    })
    if (length(rr) == 1L) return(rr[[1L]])
    u <- sf::st_union(sf::st_sfc(rr))[[1L]]
    if (inherits(u, "POLYGON")) u <- sf::st_polygon(lapply(u, .er_drop_collinear))
    u
  })
  sf::st_sfc(unname(geoms))
}
