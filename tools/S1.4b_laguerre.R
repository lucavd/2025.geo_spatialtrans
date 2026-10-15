# tools/S1.4b_laguerre.R — S1.4b (BL-098): prototipo del diagramma di potenza (Laguerre) e dei nuclei disgiunti.
# Pre-registrazione: results/S1.4b/S1.4b_preregistration.md (commit ccf1c0e). Prototipo: NON e' codice del pacchetto
# (l'integrazione in R/04b3, R/04b4 e' di S1.4c). Dipendenze: geometry (qhull), geos, sf.
# Convenzioni: µm; potenza pow_i(p) = |p - c_i|^2 - w_i; cella_i = {p : pow_i(p) <= pow_j(p) per ogni j}.
# Teorema usato da V-a: se i dischi (c_i, R_i) sono disgiunti e w_i = R_i^2, ogni disco sta nella propria cella.

#' Ambiente del prototipo (per i mutanti di C-mut: le funzioni si leggono dall'ambiente, come tools/S1.4_mutants.R)
lag_env <- function(file = "tools/S1.4b_laguerre.R") { e <- new.env(); sys.source(file, envir = e); e }

#' Coppie di vicini del diagramma di potenza = spigoli della triangolazione regolare, cioe' dell'inviluppo convesso
#' inferiore dei punti sollevati (x, y, x^2 + y^2 - w) (coordinate centrate). Restituisce list(pairs = matrice k x 2
#' con i < j, hidden = indici dei punti assenti dall'inviluppo inferiore, cioe' con cella vuota).
lag_neighbors <- function(x, y, w) {
  n <- length(x); stopifnot(n == length(y), n == length(w), n >= 1)
  mx <- mean(x); my <- mean(y); xc <- x - mx; yc <- y - my
  if (n <= 4) {                                   # qhull 3D richiede > 4 punti in posizione generale: tutte le coppie
    e <- if (n >= 2) t(utils::combn(n, 2)) else matrix(integer(0), 0, 2)    # (semipiani in piu' non cambiano la cella)
    return(list(pairs = e, hidden = integer(0), center = c(mx, my)))
  }
  z <- xc^2 + yc^2 - w
  h <- geometry::convhulln(cbind(xc, yc, z), options = "Qt n")
  low <- h$normals[, 3] < 0
  tri <- h$hull[low, , drop = FALSE]
  e <- rbind(tri[, 1:2], tri[, 2:3], tri[, c(1, 3)])
  e <- cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2]))
  e <- unique(e)
  list(pairs = e, hidden = setdiff(seq_len(n), unique(as.vector(tri))), center = c(mx, my))
}

# Ritaglio di Sutherland–Hodgman di un poligono convesso (matrice n x 2, senza ripetizione del primo vertice)
# con il semipiano a . p <= b.
.lag_clip <- function(P, a, b) {
  if (is.null(P) || nrow(P) == 0) return(P)
  s <- P[, 1] * a[1] + P[, 2] * a[2] - b
  inside <- s <= 0
  if (all(inside)) return(P)
  if (!any(inside)) return(P[0, , drop = FALSE])
  m <- nrow(P); nx <- c(2:m, 1L)
  out <- vector("list", 2 * m); k <- 0L
  for (i in seq_len(m)) {
    j <- nx[i]
    if (inside[i]) { k <- k + 1L; out[[k]] <- P[i, ] }
    if (inside[i] != inside[j]) {
      t <- s[i] / (s[i] - s[j]); k <- k + 1L; out[[k]] <- P[i, ] + t * (P[j, ] - P[i, ])
    }
  }
  do.call(rbind, out[seq_len(k)])
}

#' Celle convesse (non ritagliate alla regione) nel riquadro `box` = c(xmin, ymin, xmax, ymax).
#' Coordinate centrate su `nb$center` per la precisione; restituite nelle coordinate originali.
lag_cells <- function(x, y, w, box, nb = lag_neighbors(x, y, w)) {
  n <- length(x); mx <- nb$center[1]; my <- nb$center[2]
  xc <- x - mx; yc <- y - my; q <- xc^2 + yc^2 - w
  B <- rbind(c(box[1], box[2]), c(box[3], box[2]), c(box[3], box[4]), c(box[1], box[4])) - matrix(c(mx, my), 4, 2, byrow = TRUE)
  pr <- nb$pairs; adj <- split(c(pr[, 2], pr[, 1]), c(pr[, 1], pr[, 2]))
  cells <- vector("list", n)
  for (i in seq_len(n)) {
    P <- B; js <- adj[[as.character(i)]]
    if (i %in% nb$hidden || is.null(js)) { cells[[i]] <- if (i %in% nb$hidden) P[0, , drop = FALSE] else P; next }
    for (j in js) {
      # pow_i <= pow_j  <=>  2 p.(c_j - c_i) <= q_j - q_i
      P <- .lag_clip(P, 2 * c(xc[j] - xc[i], yc[j] - yc[i]), q[j] - q[i])
      if (nrow(P) == 0) break
    }
    cells[[i]] <- if (nrow(P)) P + matrix(c(mx, my), nrow(P), 2, byrow = TRUE) else P
  }
  cells
}

#' Diagramma di potenza ritagliato alla regione (geos).
#' snap = TRUE (default, rappresentazione F1 di S1.3d): regione e celle portate sulla griglia di `grid` µm
#' (`geos_set_precision`), le celle non contenute intersecate con la regione con snap-rounding (`*_prec`), cosi' che i
#' lati condivisi coincidano. snap = FALSE: celle esatte dell'algoritmo (solo per i test analitici), ritaglio senza griglia.
#' Frammenti (prototipo): restano al generatore (MULTIPOLYGON); la regola D-S1.3.2 e' di S1.4c.
#' @return list(cells = geos_geometry, area, n_hidden, hidden, n_clipped, n_multi, pairs, region = regione usata)
power_diagram <- function(x, y, w, region, grid = 1e-6, snap = TRUE) {
  reg <- geos::as_geos_geometry(region)
  if (length(reg) > 1) reg <- geos::geos_unary_union(geos::geos_make_collection(reg))
  if (snap) reg <- geos::geos_set_precision(reg, grid)
  bb <- c(geos::geos_xmin(reg), geos::geos_ymin(reg), geos::geos_xmax(reg), geos::geos_ymax(reg))
  pad <- max(1, 0.1 * max(bb[3] - bb[1], bb[4] - bb[2]))
  box <- c(min(bb[1], min(x)) - pad, min(bb[2], min(y)) - pad, max(bb[3], max(x)) + pad, max(bb[4], max(y)) + pad)
  nb <- lag_neighbors(x, y, w)
  cl <- lag_cells(x, y, w, box, nb)
  nv <- vapply(cl, nrow, 1L)
  ok <- nv >= 3
  g <- rep(geos::geos_empty("polygon"), length(x))
  if (any(ok)) {
    fid <- rep(which(ok), nv[ok]); X <- unlist(lapply(cl[ok], function(P) P[, 1])); Y <- unlist(lapply(cl[ok], function(P) P[, 2]))
    g[ok] <- geos::geos_make_polygon(X, Y, feature_id = fid)
  }
  if (snap) g[ok] <- geos::geos_set_precision(g[ok], grid)
  inside <- rep(FALSE, length(x))
  inside[ok] <- geos::geos_prepared_contains(reg, g[ok])
  todo <- which(ok & !inside)
  if (length(todo)) g[todo] <- if (snap) geos::geos_intersection_prec(g[todo], reg, grid) else geos::geos_intersection(g[todo], reg)
  ty <- geos::geos_type(g)
  hid <- sort(unique(c(nb$hidden, which(!ok | geos::geos_is_empty(g)))))
  list(cells = g, area = geos::geos_area(g), n_hidden = length(hid), hidden = hid,
       n_clipped = length(todo), n_multi = sum(ty == "multipolygon"), pairs = nb$pairs, region = reg)
}

#' Coppie (i < j) con R_i + R_j > gamma * d_ij (candidate da STRtree sui quadrati di lato 2 R / gamma).
near_pairs <- function(x, y, R, gamma = 0.99) {
  h <- R / gamma
  rect <- geos::geos_create_rectangle(x - h, y - h, x + h, y + h)
  tr <- geos::geos_basic_strtree(rect)
  q <- geos::geos_basic_strtree_query(tr, rect)
  i <- q$x; j <- q$tree; k <- i < j; i <- i[k]; j <- j[k]
  d <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
  bad <- R[i] + R[j] > gamma * d
  cbind(i = i[bad], j = j[bad], d = d[bad])
}

#' Fattore di scala dei nuclei (V-a/V-b): s_i = min(1, min_j gamma d_ij / (R_i + R_j), gamma d_reg_i / R_i),
#' con R = raggio circoscritto del poligono (k32 * r) e d_reg = distanza del centro dal bordo della regione.
#' Una passata; ogni vincolo e' soddisfatto perche' entrambi i membri della coppia sono scalati almeno di f_ij.
#' @return list(s, s_pair, s_reg, n_pairs)
nonoverlap_scale <- function(x, y, R, d_reg, gamma = 0.99) {
  pr <- near_pairs(x, y, R, gamma)
  s_pair <- rep(1, length(x))
  if (nrow(pr)) {
    f <- gamma * pr[, "d"] / (R[pr[, "i"]] + R[pr[, "j"]])
    m <- tapply(c(f, f), c(pr[, "i"], pr[, "j"]), min)
    id <- as.integer(names(m)); s_pair[id] <- pmin(s_pair[id], m)
  }
  s_reg <- pmin(1, gamma * d_reg / R)
  list(s = pmin(s_pair, s_reg), s_pair = s_pair, s_reg = s_reg, n_pairs = nrow(pr))
}

#' Poligoni dei nuclei circolari (n_vertices vertici, raggio circoscritto Rc) come geos_geometry.
nucleus_polys <- function(x, y, Rc, nv = 32L) {
  t <- 2 * pi * (seq_len(nv) - 1) / nv
  X <- x + outer(Rc, cos(t)); Y <- y + outer(Rc, sin(t))
  geos::geos_make_polygon(as.vector(t(X)), as.vector(t(Y)), feature_id = rep(seq_along(x), each = nv))
}

k32 <- sqrt(2 * pi / (32 * sin(2 * pi / 32)))

#' Nuclei di V-a/V-b dalle aree estratte (nucleus_area_target di derive_cell_geometry(), stesso flusso casuale).
#' @return data.frame(a_tgt, s, s_pair, s_reg, a_n, Rc) e n_pairs.
vab_nuclei <- function(x, y, a_tgt, region_boundary, gamma = 0.99) {
  r <- sqrt(a_tgt / pi); R <- k32 * r
  d_reg <- geos::geos_distance(geos::geos_make_point(x, y), region_boundary)
  sc <- nonoverlap_scale(x, y, R, d_reg, gamma)
  list(df = data.frame(a_tgt = a_tgt, s = sc$s, s_pair = sc$s_pair, s_reg = sc$s_reg, a_n = sc$s^2 * a_tgt, Rc = sc$s * R),
       n_pairs = sc$n_pairs)
}
