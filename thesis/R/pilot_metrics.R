# pilot_metrics.R — performance inside one replicate (ADEMP §5.2, §7.2).
# est, truth: spots x types (same dimnames); est rows with NA = failed spots.

# Failed spot: dropped/unsolved (NA), or proportions not summing to 1 within `tol`.
failed_spots <- function(est, tol = 1e-6) {
  apply(est, 1, function(p) anyNA(p) || abs(sum(p) - 1) > tol)
}

rmse <- function(est, truth) sqrt(mean((est - truth)^2))

# Mean over spots of the Jensen–Shannon divergence (base 2, between 0 and 1).
jsd <- function(est, truth) {
  kl <- function(p, q) { i <- p > 0; sum(p[i] * log2(p[i] / q[i])) }
  mean(vapply(seq_len(nrow(est)), function(s) {
    p <- est[s, ]; q <- truth[s, ]; m <- (p + q) / 2
    0.5 * kl(p, m) + 0.5 * kl(q, m)
  }, numeric(1)))
}
