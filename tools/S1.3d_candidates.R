# tools/S1.3d_candidates.R — S1.3d (BL-084): candidati di overlay V0, F1..F6 con interfaccia unica per O1 e O2.
# Pre-registrazione: results/S1.3d/S1.3d_preregistration.md (ee6208b). Definizioni operative: results/S1.3d/S1.3d_dev_notes.md.
# Richiede tools/S1.3d_common.R (sf, geos, polyclip, exactextractr, terra; helper .tv_* di R/04b3).
CANDS_O1 <- c("V0", "F1", "F2", "F3", "F4", "F5")
CANDS_O2 <- c("V0", "F1", "F2", "F3", "F4", "F5", "F6")

# ---- I/O bit-esatto delle geometrie (parquet double) ------------------------------------------------------------
#' sfc POLYGON/MULTIPOLYGON -> data.frame(cell, part, ring, k, x, y) senza il vertice di chiusura.
geom_to_df <- function(sfc, ids = seq_along(sfc)) {
  out <- vector("list", length(sfc))
  for (q in seq_along(sfc)) {
    g <- sfc[[q]]; parts <- if (inherits(g, "MULTIPOLYGON")) unclass(g) else list(unclass(g)); rr <- list()
    for (p in seq_along(parts)) for (r in seq_along(parts[[p]])) {
      m <- parts[[p]][[r]]; nv <- nrow(m) - 1L
      rr[[length(rr) + 1L]] <- data.frame(cell = ids[q], part = p, ring = r, k = seq_len(nv), x = m[seq_len(nv), 1], y = m[seq_len(nv), 2])
    }
    out[[q]] <- do.call(rbind, rr)
  }
  do.call(rbind, out)
}
#' inverso di geom_to_df (POLYGON se una parte, MULTIPOLYGON altrimenti); celle senza righe -> POLYGON vuoto.
df_to_sfc <- function(df, n = max(df$cell)) {
  sp <- split(df, factor(df$cell, levels = seq_len(n)))
  st_sfc(lapply(sp, function(d) {
    if (!nrow(d)) return(st_polygon())
    parts <- lapply(split(d, d$part), function(dp) lapply(split(dp, dp$ring), function(dr) { m <- cbind(dr$x[order(dr$k)], dr$y[order(dr$k)]); rbind(m, m[1, ]) }))
    if (length(parts) == 1) st_polygon(parts[[1]]) else st_multipolygon(unname(parts))
  }))
}

# ---- anelli per polyclip (esterni CCW, buchi CW, senza chiusura) ------------------------------------------------
sarea <- function(x, y) { n <- length(x); k <- c(2:n, 1); sum(x * y[k] - x[k] * y) / 2 }
rings_of <- function(g) {
  parts <- if (inherits(g, "MULTIPOLYGON")) unclass(g) else list(unclass(g)); out <- list()
  for (p in parts) for (r in seq_along(p)) {
    m <- p[[r]]; m <- m[-nrow(m), , drop = FALSE]; x <- m[, 1]; y <- m[, 2]; a <- sarea(x, y)
    if ((r == 1 && a < 0) || (r > 1 && a > 0)) { x <- rev(x); y <- rev(y) }
    out[[length(out) + 1L]] <- list(x = x, y = y)
  }
  out
}
pc_area  <- function(R) if (!length(R)) 0 else sum(vapply(R, function(r) sarea(r$x, r$y), 0))
pc_perim <- function(R) if (!length(R)) 0 else sum(vapply(R, function(r) { n <- length(r$x); k <- c(2:n, 1); sum(sqrt((r$x[k] - r$x)^2 + (r$y[k] - r$y)^2)) }, 0))
PC <- function(A, B, op) polyclip(A, B, op = op, fillA = "nonzero", fillB = "nonzero", x0 = 0, y0 = 0, eps = G_GRID)

# ---- preparazione dei candidati ----------------------------------------------------------------------------------
#' Rappresentazione del candidato. terr/reg = sfc V0 (geometria di tessellate_voronoi()); f4 = sfc F4 (solo per F4).
cand_prepare <- function(cand, terr, reg, f4 = NULL) {
  t0 <- proc.time()[[3]]; g <- G_GRID; info <- list()
  C <- switch(cand,
    V0 = list(kind = "sf", terr = terr, reg = reg),
    F1 = list(kind = "geosprec", terr = geos_set_precision(as_geos_geometry(terr), g), reg = geos_set_precision(as_geos_geometry(reg), g)),
    F2 = { tr <- s3b_round(terr, g); rr <- s3b_round(reg, g); bad <- which(!st_is_valid(tr) | is.na(st_is_valid(tr)))
           info$invalid_idx <- bad
           if (length(bad)) { fx <- st_make_valid(tr[bad])        # parte poligonale; se il collasso non lascia poligoni: cella vuota
             tr[bad] <- st_sfc(lapply(seq_along(bad), function(q) { g <- fx[[q]]
               if (inherits(g, c("POLYGON", "MULTIPOLYGON"))) return(g)
               if (inherits(g, "GEOMETRYCOLLECTION")) { p <- Filter(function(z) inherits(z, c("POLYGON", "MULTIPOLYGON")), unclass(g))
                 if (length(p)) return(st_union(st_sfc(p))[[1]]) }
               st_polygon() })) }
           info$still_invalid <- sum(!st_is_valid(tr)); list(kind = "sf", terr = tr, reg = rr) },
    F3 = list(kind = "sf", terr = st_as_sfc(as_wkb(geos_set_precision(as_geos_geometry(terr), g))),
              reg = st_as_sfc(as_wkb(geos_set_precision(as_geos_geometry(reg), g)))),
    F4 = { stopifnot(!is.null(f4)); list(kind = "sf", terr = f4, reg = reg) },
    F5 = list(kind = "clip", rings = lapply(seq_along(terr), function(i) rings_of(terr[[i]])), reg_rings = rings_of(reg[[1]]), terr = terr),
    F6 = list(kind = "ee", terr = terr),
    stop("candidato sconosciuto ", cand))
  C$cand <- cand; C$info <- info; C$t_prepare <- proc.time()[[3]] - t0
  C
}

# ---- O1: unione, differenza simmetrica, coppie in entrambi gli ordini, aree per cella -----------------------------
#' pairs = matrice (i, j), i < j, delle coppie che si toccano (predicato st_intersects su V0). Eccezioni registrate per parte.
o1_run <- function(C, pairs) {
  err <- character(0); tm <- c(); safe <- function(tag, expr) tryCatch(expr, error = function(e) { err[[tag]] <<- conditionMessage(e); NA_real_ })
  t0 <- proc.time()[[3]]
  if (C$kind == "sf") {
    a <- as.numeric(st_area(C$terr))
    U <- safe("union", if (C$cand == "F4") st_union(C$terr, is_coverage = TRUE) else st_union(C$terr))
    ua <- if (inherits(U, "sfc")) sum(as.numeric(st_area(U))) else NA_real_
    sd <- if (inherits(U, "sfc")) safe("symd", sum(as.numeric(st_area(st_sym_difference(U, C$reg))))) else NA_real_
    tm["union"] <- proc.time()[[3]] - t0; t0 <- proc.time()[[3]]
    pa <- function(x, y) vapply(seq_along(x), function(k) tryCatch(sum(as.numeric(st_area(st_intersection(x[k], y[k])))), error = function(e) NA_real_), 0)
    pij <- pa(C$terr[pairs[, 1]], C$terr[pairs[, 2]]); pji <- pa(C$terr[pairs[, 2]], C$terr[pairs[, 1]])
  } else if (C$kind == "geosprec") {
    g <- G_GRID; a <- geos_area(C$terr)
    U <- safe("union", geos_unary_union_prec(geos_make_collection(C$terr), g))
    ua <- if (inherits(U, "geos_geometry")) geos_area(U) else NA_real_
    sd <- if (inherits(U, "geos_geometry")) safe("symd", geos_area(geos_sym_difference_prec(U, C$reg, g))) else NA_real_
    tm["union"] <- proc.time()[[3]] - t0; t0 <- proc.time()[[3]]
    pij <- safe("pairs", geos_area(geos_intersection_prec(C$terr[pairs[, 1]], C$terr[pairs[, 2]], g)))
    pji <- safe("pairs", geos_area(geos_intersection_prec(C$terr[pairs[, 2]], C$terr[pairs[, 1]], g)))
  } else if (C$kind == "clip") {
    a <- vapply(C$rings, pc_area, 0)
    R <- C$rings; U <- safe("union", PC(unlist(R[-1], recursive = FALSE), R[[1]], "union"))
    ua <- if (is.list(U)) pc_area(U) else NA_real_
    sd <- if (is.list(U)) safe("symd", pc_area(PC(U, C$reg_rings, "xor"))) else NA_real_
    tm["union"] <- proc.time()[[3]] - t0; t0 <- proc.time()[[3]]
    pa <- function(I, J) vapply(seq_along(I), function(k) tryCatch(pc_area(PC(R[[I[k]]], R[[J[k]]], "intersection")), error = function(e) NA_real_), 0)
    pij <- pa(pairs[, 1], pairs[, 2]); pji <- pa(pairs[, 2], pairs[, 1])
  } else stop("O1 non definita per ", C$cand)
  tm["pairs"] <- proc.time()[[3]] - t0
  n_pair_err <- sum(is.na(pij)) + sum(is.na(pji))
  list(cand = C$cand, area = a, union_a = ua, symd_a = sd, pij = pij, pji = pji, errors = err, n_pair_err = n_pair_err,
       time = c(prepare = C$t_prepare, tm), info = C$info)
}

# ---- O2: cella x bin -----------------------------------------------------------------------------------------------
#' Griglia allineata: origine (X0, Y0) multiplo di 8 µm, nx8 x ny8 bin da 8 µm; lato s in {2, 8}; indice bin k = ix + iy*nx (0-based).
grid_spec <- function(bb) { X0 <- floor(bb[["xmin"]] / 8) * 8; Y0 <- floor(bb[["ymin"]] / 8) * 8
  list(X0 = X0, Y0 = Y0, nx8 = ceiling((bb[["xmax"]] - X0) / 8), ny8 = ceiling((bb[["ymax"]] - Y0) / 8)) }
grid_sfc <- function(gs, s) {
  nx <- gs$nx8 * 8 / s; ny <- gs$ny8 * 8 / s; ix <- rep(0:(nx - 1), times = ny); iy <- rep(0:(ny - 1), each = nx)
  x0 <- gs$X0 + s * ix; y0 <- gs$Y0 + s * iy
  sq <- st_sfc(lapply(seq_along(ix), function(k) st_polygon(list(rbind(c(x0[k], y0[k]), c(x0[k] + s, y0[k]), c(x0[k] + s, y0[k] + s), c(x0[k], y0[k] + s), c(x0[k], y0[k]))))))
  list(sfc = sq, ix = ix, iy = iy, nx = nx, ny = ny)
}
o2_run <- function(C, gs, s, grid = NULL) {
  t0 <- proc.time()[[3]]; err <- character(0)
  if (C$kind != "ee" && is.null(grid)) grid <- grid_sfc(gs, s)
  res <- tryCatch({
    if (C$kind == "sf") {
      pc <- st_intersection(C$terr, grid$sfc); ix <- attr(pc, "idx")
      keep <- as.character(st_geometry_type(pc)) %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION")
      ar <- as.numeric(st_area(pc)); pe <- geos_length(as_geos_geometry(pc))
      data.frame(i = ix[, 1], b = ix[, 2], area = ar, perim = pe)[keep & ar > 0, ]
    } else if (C$kind == "geosprec") {
      gg <- as_geos_geometry(grid$sfc); im <- geos_intersects_matrix(C$terr, gg); ii <- rep(seq_along(im), lengths(im)); jj <- unlist(im)
      pc <- geos_intersection_prec(C$terr[ii], gg[jj], G_GRID); ar <- geos_area(pc)
      data.frame(i = ii, b = jj, area = ar, perim = geos_length(pc))[ar > 0, ]
    } else if (C$kind == "clip") {
      im <- st_intersects(C$terr, grid$sfc); ii <- rep(seq_along(im), lengths(im)); jj <- unlist(im)
      sqr <- function(b) { x0 <- gs$X0 + s * grid$ix[b]; y0 <- gs$Y0 + s * grid$iy[b]; list(list(x = c(x0, x0 + s, x0 + s, x0), y = c(y0, y0, y0 + s, y0 + s))) }
      ap <- t(vapply(seq_along(ii), function(k) { R <- PC(C$rings[[ii[k]]], sqr(jj[k]), "intersection"); c(pc_area(R), pc_perim(R)) }, numeric(2)))
      data.frame(i = ii, b = jj, area = ap[, 1], perim = ap[, 2])[ap[, 1] > 0, ]
    } else if (C$kind == "ee") {
      nx <- gs$nx8 * 8 / s; ny <- gs$ny8 * 8 / s
      r <- terra::rast(xmin = gs$X0, xmax = gs$X0 + 8 * gs$nx8, ymin = gs$Y0, ymax = gs$Y0 + 8 * gs$ny8, resolution = s, crs = "")
      terra::values(r) <- 1
      ee <- exact_extract(r, st_sf(id = seq_along(C$terr), geometry = C$terr), include_cell = TRUE, progress = FALSE)
      ii <- rep(seq_along(ee), vapply(ee, nrow, 1L)); cell <- unlist(lapply(ee, `[[`, "cell")); cf <- unlist(lapply(ee, `[[`, "coverage_fraction"))
      row <- (cell - 1) %/% nx; col <- (cell - 1) %% nx; iy <- ny - 1 - row
      data.frame(i = ii, b = col + iy * nx + 1, area = cf * s^2, perim = NA_real_)[cf > 0, ]
    } else stop("O2 non definita per ", C$cand)
  }, error = function(e) { err[["o2"]] <<- conditionMessage(e); NULL })
  nx <- gs$nx8 * 8 / s
  if (!is.null(res)) { res$ix <- (res$b - 1) %% nx; res$iy <- (res$b - 1) %/% nx; res$b <- NULL }
  list(cand = C$cand, s = s, pieces = res, errors = err, time = c(prepare = C$t_prepare, o2 = proc.time()[[3]] - t0))
}
#' O3: aggregazione 2 -> 8 µm dei pezzi a 2 µm del candidato
o3_aggregate <- function(p2) {
  if (is.null(p2)) return(NULL)
  k <- paste(p2$i, p2$ix %/% 4, p2$iy %/% 4); ag <- tapply(p2$area, k, sum)
  sp <- do.call(rbind, strsplit(names(ag), " ", fixed = TRUE))
  data.frame(i = as.integer(sp[, 1]), ix = as.integer(sp[, 2]), iy = as.integer(sp[, 3]), area = as.numeric(ag), perim = NA_real_)
}
