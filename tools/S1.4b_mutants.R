# tools/S1.4b_mutants.R — S1.4b, C-mut: mutanti del prototipo Laguerre (pre-registrazione, commit ccf1c0e).
#  MU1 pesi ignorati (w = 0, cioe' Voronoi); MU2 segno del peso (pow = |p-c|^2 + w); MU3 non sovrapposizione rimossa;
#  MU4 vicini dai 4 piu' vicini euclidei invece della triangolazione regolare. MU5 (sollevamento dell'arbitro) e' in
#  tools/S1.4b_arbiter.py --mutant MU5.
source("tools/S1.4b_laguerre.R")
lag_mutate <- function(e, m = "none") {
  if (m == "none") return(e)
  if (m == "MU1") { f <- e$power_diagram; e$power_diagram <- function(x, y, w, region, ...) f(x, y, 0 * w, region, ...) }
  else if (m == "MU2") { f <- e$power_diagram; e$power_diagram <- function(x, y, w, region, ...) f(x, y, -w, region, ...) }
  else if (m == "MU3") { e$nonoverlap_scale <- function(x, y, R, d_reg, gamma = 0.99) {
      s_reg <- pmin(1, gamma * d_reg / R); list(s = s_reg, s_pair = rep(1, length(x)), s_reg = s_reg, n_pairs = 0L) } }
  else if (m == "MU4") { e$lag_neighbors <- function(x, y, w) {
      k <- 4; n <- length(x); pr <- NULL
      for (i in seq_len(n)) { d <- (x - x[i])^2 + (y - y[i])^2; d[i] <- Inf; j <- order(d)[1:k]; pr <- rbind(pr, cbind(i, j)) }
      pr <- unique(cbind(pmin(pr[, 1], pr[, 2]), pmax(pr[, 1], pr[, 2])))
      list(pairs = pr, hidden = integer(0), center = c(mean(x), mean(y))) } }
  else stop("mutante sconosciuto: ", m)
  # le funzioni originali hanno ambiente e (sys.source): power_diagram trova il lag_neighbors mutato in e
  e
}
