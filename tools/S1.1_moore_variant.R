# S1.1 — controprova CP-1: contorno di Moore "letterale" (design v1.1 §5.1, passo 3)
# Traccia il bordo esterno di una componente passando per i CENTRI dei pixel
# (Moore-neighbour tracing; arresto quando dal pixel di partenza si ripete la
# prima mossa: il criterio di Jacob da solo non termina sulle linee larghe 1 px) e ne calcola l'area
# con la formula di Gauss (shoelace). Non gestisce i buchi, come nel design.
# Serve solo a dimostrare che il check di area C1 ha potere discriminante.

#' @param B matrice logica [riga, colonna] della componente (TRUE = pixel)
#' @return area in pixel^2 del poligono per i centri (0 per un pixel o una linea)
moore_centre_area <- function(B) {
  B <- rbind(FALSE, cbind(FALSE, B, FALSE), FALSE)          # bordo di fondo
  idx <- which(B, arr.ind = TRUE)
  if (nrow(idx) <= 1L) return(0)
  s <- idx[order(idx[, 1], idx[, 2])[1L], ]                 # primo pixel in ordine di riga
  # 8 vicini in senso orario (y verso il basso), a partire da ovest
  off <- matrix(c(0, -1, -1, -1, -1, 0, -1, 1, 0, 1, 1, 1, 1, 0, 1, -1), ncol = 2, byrow = TRUE)
  p <- s; b <- s + off[1L, ]                                # backtrack iniziale: ovest (fondo)
  q1 <- NULL
  path <- matrix(NA_real_, nrow = 4L * nrow(idx) + 8L, ncol = 2L); n <- 1L; path[1L, ] <- p
  for (iter in seq_len(4L * nrow(idx) + 8L)) {
    d0 <- which(off[, 1] == b[1] - p[1] & off[, 2] == b[2] - p[2])
    found <- FALSE
    for (k in 1:8) {
      d <- ((d0 - 1L + k) %% 8L) + 1L
      q <- p + off[d, ]
      if (B[q[1], q[2]]) { found <- TRUE; break }
    }
    if (!found) break
    if (is.null(q1)) q1 <- q
    else if (all(p == s) && all(q == q1)) break              # stessa prima mossa: anello chiuso
    b <- p + off[((d - 2L) %% 8L) + 1L, ]                    # vicino esaminato prima di q
    p <- q
    n <- n + 1L; path[n, ] <- p
  }
  P <- path[seq_len(n), , drop = FALSE]
  if (n > 1L && all(P[n, ] == s)) P <- P[-n, , drop = FALSE]  # ultimo = partenza
  n <- nrow(P)
  if (n < 3L) return(0)
  xr <- P[, 2]; yr <- P[, 1]
  abs(sum(xr * c(yr[-1L], yr[1L]) - c(xr[-1L], xr[1L]) * yr)) / 2
}
