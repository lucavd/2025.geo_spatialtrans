# tools/S1.3_arbiter.R — arbitro di C-10a (addendum): cella di Voronoi esatta per ritaglio di semipiani
# (Sutherland–Hodgman) in R puro, coordinate locali al generatore, indipendente da GEOS e da deldir.
# K vicini raddoppiato fino al certificato d_(K+1) > 2 max||v|| (corretto dopo la revisione, RA-code-08); lati contati se > 1e-9 µm (BL-075).
hp_cell <- function(i, x, y, K0 = 32L) {
  # Certificato (dopo la revisione, RA-code-08): un generatore a distanza d dal generatore i non puo' tagliare la cella
  # se d > 2 * max||v|| (v = vertici della cella, coordinate locali). Si raddoppia K finche' d_(K+1) > 2 * max||v||
  # o finche' si usano tutti i punti (allora bounded = FALSE se la cella tocca il riquadro iniziale).
  n <- length(x); d2 <- (x - x[i])^2 + (y - y[i])^2; o <- order(d2); K <- min(K0, n - 1L)
  repeat {
    nb <- o[2:(K + 1)]; R <- 4 * sqrt(d2[nb[K]]) + 1
    P <- rbind(c(-R, -R), c(R, -R), c(R, R), c(-R, R))
    for (j in nb) {
      q <- c(x[j] - x[i], y[j] - y[i]); cst <- sum(q^2) / 2; f <- P %*% q - cst; m <- nrow(P); out <- NULL
      for (k in seq_len(m)) { k2 <- if (k == m) 1L else k + 1L; a <- P[k, ]; b <- P[k2, ]; fa <- f[k]; fb <- f[k2]
        if (fa <= 0) out <- rbind(out, a)
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0)) out <- rbind(out, a + (b - a) * fa / (fa - fb)) }
      P <- out; if (is.null(P)) break
    }
    rmax <- sqrt(max(rowSums(P^2)))
    certified <- K < n - 1L && sqrt(d2[o[K + 2]]) > 2 * rmax
    if (certified || K >= n - 1L) break
    K <- min(2L * K, n - 1L)
  }
  xs <- P[, 1]; ys <- P[, 2]
  area <- abs(sum(xs * c(ys[-1], ys[1]) - c(xs[-1], xs[1]) * ys)) / 2
  e <- sqrt(rowSums((P - P[c(2:nrow(P), 1), , drop = FALSE])^2))
  c(area = area, nsides = sum(e > 1e-9), min_edge = min(e), K = K, bounded = certified || rmax < R * 0.999)
}
#' Tile discordanti fra due motori (aree rel. > 1e-9 o lati diversi) fra quelli selezionati, con l'arbitro.
arbitrate <- function(tg, td, x, y, sel = rep(TRUE, length(tg))) {
  ag <- as.numeric(sf::st_area(tg)); ad <- as.numeric(sf::st_area(td))
  nsd <- function(t) vapply(t, function(q) { m <- unclass(q)[[1]]; sum(sqrt(rowSums(diff(m)^2)) > 1e-9) }, 1L); ng <- nsd(tg); nd <- nsd(td)   # lati > 1e-9 µm (BL-075)
  k <- which(sel & (abs(ag - ad) / ad > 1e-9 | ng != nd))
  if (!length(k)) return(data.frame(cell = integer(0), a_geos = numeric(0), a_deldir = numeric(0), a_hp = numeric(0), ns_geos = integer(0),
                                    ns_deldir = integer(0), ns_hp = integer(0), min_edge = numeric(0), bounded = logical(0), geos_ok = logical(0), deldir_ok = logical(0)))
  h <- t(vapply(k, hp_cell, numeric(5), x = x, y = y))
  data.frame(cell = k, a_geos = ag[k], a_deldir = ad[k], a_hp = h[, "area"], ns_geos = ng[k], ns_deldir = nd[k], ns_hp = h[, "nsides"],
             min_edge = h[, "min_edge"], bounded = h[, "bounded"] == 1,
             geos_ok = abs(ag[k] - h[, "area"]) / h[, "area"] <= 1e-9 & ng[k] == h[, "nsides"],
             deldir_ok = abs(ad[k] - h[, "area"]) / h[, "area"] <= 1e-9 & nd[k] == h[, "nsides"])
}
