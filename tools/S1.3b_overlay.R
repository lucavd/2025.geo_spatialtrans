# tools/S1.3b_overlay.R — S1.3b (BL-069): varianti di precisione per le operazioni di sovrapposizione sui territori.
# Pre-registrazione: results/S1.3b/S1.3b_preregistration.md (commit d96fea1).
# Varianti: V0 grezzo (sf, GEOS di sistema); V0g grezzo con il pacchetto geos (GEOS 3.11.1); V1 arrotondamento sf a griglia g
# (st_set_precision -> WKB, materializzato); V2 snap-rounding GEOS a griglia g (geos_set_precision + *_prec);
# V3 annodamento esplicito senza griglia (fusione dei vertici entro tol, giunzioni a T inserite).
suppressPackageStartupMessages({ library(sf); library(geos); library(wk) })
`%||%` <- function(a, b) if (is.null(a)) b else a
G_CHECK <- c(V0 = 0, V0g = 0, V1 = 1e-6, V2 = 1e-6, V3 = 1e-7)   # g usato nelle tolleranze tau (pre-registrato)

#' Perimetro (anelli esterni + buchi) di ogni geometria; via GEOSLength, che per i poligoni restituisce il perimetro.
s3b_perim <- function(x) geos_length(if (inherits(x, "geos_geometry")) x else as_geos_geometry(x))

#' V1: coordinate arrotondate come le vede GEOS quando sf ha una precisione (round(x * 1/g) / (1/g)), materializzate.
s3b_round <- function(x, g) { y <- st_as_sfc(st_as_binary(st_set_precision(x, 1 / g))); st_precision(y) <- 0; y }

s3b_sfg_list <- function(x) lapply(seq_along(x), function(i) x[[i]])

#' Costruisce una variante: territori e regione nella rappresentazione della variante.
s3b_variant <- function(terr, reg, v, g = NULL, tol = NULL, round_terr = TRUE) {
  switch(v,
    V0  = list(v = v, kind = "sf", terr = terr, reg = reg, g = 0),
    V0g = list(v = v, kind = "geos", prec = FALSE, terr = as_geos_geometry(terr), reg = as_geos_geometry(reg), g = 0),
    V1  = { g <- g %||% 1e-6
            list(v = v, kind = "sf", terr = if (round_terr) s3b_round(terr, g) else terr, reg = s3b_round(reg, g), g = g) },
    V2  = { g <- g %||% 1e-6
            list(v = v, kind = "geos", prec = TRUE, terr = geos_set_precision(as_geos_geometry(terr), g),
                 reg = geos_set_precision(as_geos_geometry(reg), g), g = g) },
    V3  = { tol <- tol %||% 1e-7; nd <- s3b_node(c(list(reg[[1]]), s3b_sfg_list(terr)), tol)
            list(v = v, kind = "sf", terr = st_sfc(nd$geoms[-1]), reg = st_sfc(nd$geoms[1]), g = tol, info = nd$info, d = nd$d) },
    stop("variante sconosciuta ", v))
}

s3b_area  <- function(V, x) if (V$kind == "sf") as.numeric(st_area(x)) else geos_area(x)
s3b_union <- function(V) if (V$kind == "sf") st_union(V$terr) else if (V$prec)
  geos_unary_union_prec(geos_make_collection(V$terr), V$g) else geos_unary_union(geos_make_collection(V$terr))
s3b_symd  <- function(V, a, b) if (V$kind == "sf") st_sym_difference(a, b) else if (V$prec) geos_sym_difference_prec(a, b, V$g) else geos_sym_difference(a, b)
#' intersezione elemento per elemento (a[k] con b[k]); restituisce le aree
s3b_inter_area <- function(V, a, b) {
  if (!length(a)) return(numeric(0))
  if (V$kind == "sf") vapply(seq_along(a), function(k) sum(as.numeric(st_area(st_intersection(a[k], b[k])))), 0)
  else geos_area(if (V$prec) geos_intersection_prec(a, b, V$g) else geos_intersection(a, b))
}
s3b_sf <- function(V, x) if (V$kind == "sf") x else st_as_sfc(as_wkb(x))

#' Coppie i < j con interni sovrapposti (predicato, affidabile).
s3b_overlap_pairs <- function(ts) {
  pr <- st_relate(ts, ts, pattern = "2********")
  p <- do.call(rbind, lapply(seq_along(pr), function(q) { j <- pr[[q]]; j <- j[j > q]; if (length(j)) cbind(q, j) }))
  if (is.null(p)) matrix(integer(0), 0, 2) else p
}

#' Misure C-1..C-4 (e aree per C-5) di una variante. P = perimetri dei territori V0, Preg = perimetro regione V0,
#' pts = punti per la copertura (NULL = salta C-3), gchk = g delle tolleranze.
s3b_measure <- function(V, P, Preg, pts = NULL, gchk = G_CHECK[[V$v]], full = TRUE) {
  tm <- c(); t0 <- proc.time()[[3]]
  a <- s3b_area(V, V$terr); tau <- gchk * (Preg + sum(P))
  u <- s3b_union(V); ua <- sum(s3b_area(V, u)); sdA <- sum(s3b_area(V, s3b_symd(V, u, V$reg)))
  tm["union"] <- proc.time()[[3]] - t0; t0 <- proc.time()[[3]]
  ts <- s3b_sf(V, V$terr); rs <- s3b_sf(V, V$reg)
  pp <- s3b_overlap_pairs(ts)
  pa <- s3b_inter_area(V, V$terr[pp[, 1]], V$terr[pp[, 2]]); tij <- gchk * (P[pp[, 1]] + P[pp[, 2]])
  tm["pairs"] <- proc.time()[[3]] - t0
  out <- list(variant = V$v, g = V$g, gchk = gchk, n = length(a), sum_a = sum(a), union_a = ua, ov = sum(a) - ua, symd = sdA, tau = tau,
              n_pairs = nrow(pp), pairs = data.frame(i = pp[, 1], j = pp[, 2], area = pa, tau_ij = tij, psum = P[pp[, 1]] + P[pp[, 2]]),
              n_pairs_over = sum(pa > tij), area = a)
  if (!full) { out$time <- tm; return(out) }
  t0 <- proc.time()[[3]]
  ty <- as.character(st_geometry_type(ts)); out$types <- table(ty); out$n_empty <- sum(st_is_empty(ts))
  out$n_invalid <- sum(!st_is_valid(ts)); out$n_collapsed <- sum(a <= 0)
  out$reg_valid <- st_is_valid(rs)
  if (!is.null(pts)) {
    h <- lengths(st_intersects(pts, ts)); an <- which(h != 1)
    dmin <- vapply(an, function(k) {
      cand <- st_is_within_distance(pts[k], ts, 1e-5)[[1]]
      b <- as.numeric(st_distance(pts[k], st_boundary(rs)))
      if (length(cand)) b <- c(b, as.numeric(st_distance(pts[k], st_boundary(ts[cand]))))
      min(b) }, 0)
    out$cov <- c(n_pts = length(pts), cov0 = sum(h == 0), cov1 = sum(h == 1), cov2 = sum(h >= 2),
                 n_anom = length(an), n_anom_far = sum(dmin > 1e-5), dmax = if (length(dmin)) max(dmin) else 0)
  }
  tm["integrity_cov"] <- proc.time()[[3]] - t0; out$time <- tm
  out
}

# ---- V3: annodamento esplicito ------------------------------------------------------------------------------------

#' Vertici (senza il vertice di chiusura) di una lista di POLYGON/MULTIPOLYGON: g, p (parte), r (anello), i, x, y.
s3b_geoms_df <- function(geoms) {
  out <- vector("list", length(geoms))
  for (k in seq_along(geoms)) {
    gk <- geoms[[k]]; parts <- if (inherits(gk, "MULTIPOLYGON")) unclass(gk) else list(unclass(gk))
    rr <- list()
    for (p in seq_along(parts)) for (r in seq_along(parts[[p]])) {
      m <- parts[[p]][[r]]; nv <- nrow(m) - 1L
      rr[[length(rr) + 1L]] <- cbind(k, p, r, seq_len(nv), m[seq_len(nv), 1], m[seq_len(nv), 2])
    }
    out[[k]] <- do.call(rbind, rr)
  }
  v <- do.call(rbind, out)
  data.frame(g = as.integer(v[, 1]), p = as.integer(v[, 2]), r = as.integer(v[, 3]), i = v[, 4], x = v[, 5], y = v[, 6])
}

#' Fusione dei vertici entro tol (hash su celle di lato tol + union-find); rappresentante = indice minore.
s3b_fuse <- function(x, y, tol) {
  n <- length(x); kx <- floor(x / tol); ky <- floor(y / tol)
  key <- paste(kx, ky); uk <- unique(key); cid <- match(key, uk); mem <- split(seq_len(n), cid); len_mem <- lengths(mem)
  I <- list(); J <- list()
  for (o in list(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(1, -1))) {
    same <- all(o == 0)
    c2 <- if (same) cid else match(paste(kx + o[1], ky + o[2]), uk)
    ii <- which(!is.na(c2)); if (!length(ii)) next
    i_ <- rep(ii, len_mem[c2[ii]]); j_ <- unlist(mem[c2[ii]], use.names = FALSE)
    keep <- if (same) i_ < j_ else rep(TRUE, length(i_))
    I[[length(I) + 1L]] <- i_[keep]; J[[length(J) + 1L]] <- j_[keep]
  }
  I <- unlist(I); J <- unlist(J)
  if (is.null(I)) I <- J <- integer(0)
  d <- sqrt((x[I] - x[J])^2 + (y[I] - y[J])^2); k <- d <= tol; I <- I[k]; J <- J[k]; d <- d[k]
  lab <- seq_len(n)
  if (length(I)) repeat {
    old <- lab; m <- pmin(lab[I], lab[J]); o <- order(m, decreasing = TRUE)
    lab[I[o]] <- pmin(old[I[o]], m[o]); l2 <- lab; lab[J[o]] <- pmin(l2[J[o]], m[o])
    lab <- lab[lab]
    if (identical(lab, old)) break
  }
  list(rep = lab, I = I, J = J, d = d)
}

#' Giunzioni a T: vertici di una geometria entro tol dall'interno di un lato di un'altra geometria (hash dei lati su celle).
s3b_tjunctions <- function(v, tol, cell = 1) {
  n <- nrow(v); ring <- match(paste(v$g, v$p, v$r), unique(paste(v$g, v$p, v$r)))
  last <- c(ring[-1] != ring[-n], TRUE); first_of <- match(ring, ring)
  nxt <- seq_len(n) + 1L; nxt[last] <- first_of[last]
  ax <- v$x; ay <- v$y; bx <- v$x[nxt]; by <- v$y[nxt]
  cx0 <- floor((pmin(ax, bx) - tol) / cell); cx1 <- floor((pmax(ax, bx) + tol) / cell)
  cy0 <- floor((pmin(ay, by) - tol) / cell); cy1 <- floor((pmax(ay, by) + tol) / cell)
  ncx <- cx1 - cx0 + 1; ncy <- cy1 - cy0 + 1; ncell <- ncx * ncy
  e <- rep(seq_len(n), ncell); off <- sequence(ncell) - 1
  ecx <- cx0[e] + off %% ncx[e]; ecy <- cy0[e] + off %/% ncx[e]
  ek <- (ecx + 1e5) * 1e6 + (ecy + 1e5)
  vk <- (floor(v$x / cell) + 1e5) * 1e6 + (floor(v$y / cell) + 1e5)
  uk <- unique(ek); memE <- split(e, match(ek, uk)); vc <- match(vk, uk)
  vi <- which(!is.na(vc)); vi <- rep(vi, lengths(memE)[vc[vi]]); ej <- unlist(memE[vc[!is.na(vc)]], use.names = FALSE)
  k <- v$g[vi] != v$g[ej]; vi <- vi[k]; ej <- ej[k]
  dx <- bx[ej] - ax[ej]; dy <- by[ej] - ay[ej]; L2 <- dx^2 + dy^2
  t <- ((v$x[vi] - ax[ej]) * dx + (v$y[vi] - ay[ej]) * dy) / L2
  qx <- ax[ej] + t * dx; qy <- ay[ej] + t * dy; dd <- sqrt((v$x[vi] - qx)^2 + (v$y[vi] - qy)^2)
  is_end <- (v$x[vi] == ax[ej] & v$y[vi] == ay[ej]) | (v$x[vi] == bx[ej] & v$y[vi] == by[ej])
  hit <- L2 > 0 & t > 0 & t < 1 & dd <= tol & !is_end
  vi <- vi[hit]; ej <- ej[hit]; t <- t[hit]; dd <- dd[hit]
  if (length(vi)) {   # un vertice con le stesse coordinate gia' presente sul lato (inserimento doppio): una sola volta per lato
    kk <- !duplicated(data.frame(ej, v$x[vi], v$y[vi])); vi <- vi[kk]; ej <- ej[kk]; t <- t[kk]; dd <- dd[kk]
  }
  new <- data.frame(g = v$g[ej], p = v$p[ej], r = v$r[ej], i = v$i[ej] + t, x = v$x[vi], y = v$y[vi])
  list(new = new, n = length(vi), d = dd)
}

#' Ricostruzione: anelli ordinati, duplicati consecutivi rimossi, anelli con < 3 vertici distinti scartati (contati).
s3b_rebuild <- function(v, types) {
  v <- v[order(v$g, v$p, v$r, v$i), ]
  ig <- split(seq_len(nrow(v)), v$g); n_drop <- 0L
  geoms <- lapply(seq_along(types), function(k) {
    ii <- ig[[as.character(k)]]
    parts <- split(ii, v$p[ii])
    pl <- lapply(parts, function(pi) {
      rings <- split(pi, v$r[pi])
      rl <- lapply(rings, function(ri) {
        m <- cbind(v$x[ri], v$y[ri]); nm <- nrow(m)
        dup <- c(FALSE, m[-1, 1] == m[-nm, 1] & m[-1, 2] == m[-nm, 2]); m <- m[!dup, , drop = FALSE]
        while (nrow(m) > 1 && m[nrow(m), 1] == m[1, 1] && m[nrow(m), 2] == m[1, 2]) m <- m[-nrow(m), , drop = FALSE]
        if (nrow(m) < 3) { n_drop <<- n_drop + 1L; return(NULL) }
        rbind(m, m[1, ]) })
      if (is.null(rl[[1]])) return(NULL)
      unname(Filter(Negate(is.null), rl)) })
    pl <- unname(Filter(Negate(is.null), pl))
    if (!length(pl)) return(st_polygon())
    if (types[k] == "POLYGON" && length(pl) == 1) st_polygon(pl[[1]]) else st_multipolygon(pl) })
  list(geoms = geoms, n_drop = n_drop)
}

#' Annodamento: fusione entro tol, poi giunzioni a T entro tol. geoms[[1]] ha la precedenza come rappresentante.
s3b_node <- function(geoms, tol, tjunction = TRUE) {
  types <- vapply(geoms, function(gk) class(gk)[2], "")
  v <- s3b_geoms_df(geoms)
  if (tol <= 0) {                               # mutante M2: nessuna fusione, nessuna giunzione (ricostruzione identica)
    rb <- s3b_rebuild(v, types)
    return(list(geoms = rb$geoms, info = c(n_vertices = nrow(v), n_pairs_cross = 0, n_pairs_cross_identical = 0, n_pairs_cross_diff = 0,
                n_moved = 0, n_moved_geoms = 0, n_tjunction = 0, n_rings_dropped = rb$n_drop, d_max = 0), d = numeric(0), d_tj = numeric(0),
                moved_g = integer(0), tj_g = integer(0)))
  }
  fu <- s3b_fuse(v$x, v$y, tol); rp <- fu$rep
  cross <- v$g[fu$I] != v$g[fu$J]
  moved <- v$x != v$x[rp] | v$y != v$y[rp]
  v$x <- v$x[rp]; v$y <- v$y[rp]
  tj <- if (tjunction) s3b_tjunctions(v, tol) else list(new = v[0, ], n = 0L, d = numeric(0))
  rb <- s3b_rebuild(rbind(v, tj$new), types)
  info <- c(n_vertices = nrow(v), n_pairs_cross = sum(cross), n_pairs_cross_identical = sum(cross & fu$d == 0),
            n_pairs_cross_diff = sum(cross & fu$d > 0), n_moved = sum(moved), n_moved_geoms = length(unique(v$g[moved])),
            n_tjunction = tj$n, n_rings_dropped = rb$n_drop,
            d_max = if (any(cross & fu$d > 0)) max(fu$d[cross & fu$d > 0]) else 0)
  list(geoms = rb$geoms, info = info, d = fu$d[cross & fu$d > 0], d_tj = tj$d, moved_g = v$g[moved], tj_g = tj$new$g)
}

# ---- CP-2: perturbazione indipendente dei vertici -----------------------------------------------------------------
s3b_perturb <- function(terr, delta, seed) {
  set.seed(seed)
  st_sfc(lapply(seq_along(terr), function(i) {
    gk <- terr[[i]]; mp <- inherits(gk, "MULTIPOLYGON"); parts <- if (mp) unclass(gk) else list(unclass(gk))
    parts <- lapply(parts, function(pr) lapply(pr, function(m) { k <- nrow(m) - 1L
      m[seq_len(k), ] <- m[seq_len(k), ] + runif(2 * k, -delta, delta); m[k + 1L, ] <- m[1L, ]; m }))
    if (mp) st_multipolygon(parts) else st_polygon(parts[[1]]) }))
}
