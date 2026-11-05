# pilot_metrics.R — performance inside one replicate (ADEMP §5.2, §7.2).
# est, truth: spots x types (same dimnames); est rows with NA = failed spots.

# Failed spot: dropped/unsolved (NA), or proportions not summing to 1 within `tol`.
failed_spots <- function(est, tol = 1e-6) {
  apply(est, 1, function(p) anyNA(p) || abs(sum(p) - 1) > tol)
}

rmse <- function(est, truth) sqrt(mean((est - truth)^2))

# Null baseline (DECISIONS 2026-11-03, Luca): share of each cell type among all training cells
# (full training split, not the capped reference), over all reference types.
# n_train: named vector of training cells per reference type.
null_props <- function(n_train) n_train / sum(n_train)

# Main analysis (ADEMP §7.2): failed spots get the null-baseline proportions.
fill_failed <- function(est, p_null) {
  bad <- failed_spots(est)
  if (any(bad)) est[bad, ] <- matrix(p_null[colnames(est)], sum(bad), ncol(est), byrow = TRUE)
  est
}

# Main run (ADEMP §7.2, Pierangela 2026-11-03): a spot fails if it has a missing value, a proportion
# below -tol, or a sum that differs from 1 by more than tol; values in [-tol, 0) are numerical zeros,
# set to 0, then the row is rescaled to 1. Failed rows become NA, so fill_failed() treats them alike.
# Returns the cleaned matrix with attribute "min_raw" (smallest raw value, for the residue check).
clean_props <- function(est, tol = 1e-5) {
  min_raw <- suppressWarnings(min(est, na.rm = TRUE))
  bad <- apply(est, 1, function(p) anyNA(p) || any(p < -tol) || abs(sum(p) - 1) > tol)
  est[bad, ] <- NA
  ok <- !bad
  if (any(ok)) {
    x <- est[ok, , drop = FALSE]; x[x < 0] <- 0
    est[ok, ] <- x / rowSums(x)
  }
  attr(est, "min_raw") <- min_raw
  est
}

# Mean over spots of the Jensen–Shannon divergence (base 2, between 0 and 1).
jsd <- function(est, truth) {
  kl <- function(p, q) { i <- p > 0; sum(p[i] * log2(p[i] / q[i])) }
  mean(vapply(seq_len(nrow(est)), function(s) {
    p <- est[s, ]; q <- truth[s, ]; m <- (p + q) / 2
    0.5 * kl(p, m) + 0.5 * kl(q, m)
  }, numeric(1)))
}
