# arm_p.R — baseline gene means of arm P, after Splatter (Zappia et al. 2017, doi:10.1186/s13059-017-1305-0).
# The estimation reproduces splatEstMean and splatEstOutlier, and the simulation splatSimGeneMeans
# (splatter 1.37.1, R/splat-estimate.R and R/splat-simulate.R, GitHub Oshlack/splatter, read 2026-11-04).
# Input: one mean per gene (ADEMP §2.1: pooled relative profile of the training cells, non-marker genes).

# Splatter utils.R winsorize(): values below the q quantile set to it, above 1 - q set to it.
winsorize_q <- function(x, q) {
  lohi <- stats::quantile(x, c(q, 1 - q), na.rm = TRUE)
  if (diff(lohi) < 0) lohi <- rev(lohi)
  x[!is.na(x) & x < lohi[1]] <- lohi[1]
  x[!is.na(x) & x > lohi[2]] <- lohi[2]
  x
}

estimate_splat_means <- function(means) {
  means <- means[means != 0]
  # splatEstMean: gamma on winsorized means, maximum goodness-of-fit (Cramer-von Mises), else moments
  w <- winsorize_q(means, q = 0.1)
  fit <- fitdistrplus::fitdist(w, "gamma", method = "mge", gof = "CvM")
  method <- "mge_CvM"
  if (fit$convergence > 0) {
    fit <- fitdistrplus::fitdist(w, "gamma", method = "mme"); method <- "mme"
  }
  # splatEstOutlier: log mean above median + 2 MAD; factors = mean / median mean, log-normal fit
  lm <- log(means)
  outs <- which(lm > stats::median(lm) + 2 * stats::mad(lm))
  out <- list(prob = length(outs) / length(means), facLoc = NA_real_, facScale = NA_real_)
  if (length(outs) > 1) {
    f <- fitdistrplus::fitdist(means[outs] / stats::median(means), "lnorm")
    out$facLoc <- unname(f$estimate["meanlog"]); out$facScale <- unname(f$estimate["sdlog"])
  }
  list(shape = unname(fit$estimate["shape"]), rate = unname(fit$estimate["rate"]), gamma_method = method,
       out_prob = out$prob, out_facLoc = out$facLoc, out_facScale = out$facScale)
}

# splatSimGeneMeans: base ~ gamma; outlier genes (probability out_prob) get median(base) x log-normal factor.
# Splatter's getLNormFactors(n, sel.prob, neg.prob = 0, loc, scale): selected genes by
# rbinom(n, 1, sel.prob), factors rlnorm(loc, scale), and a factor below 1 is inverted (1/f), so with
# neg.prob = 0 every outlier factor is >= 1.
lnorm_factors <- function(n, sel_prob, loc, scale) {
  sel <- as.logical(stats::rbinom(n, 1, sel_prob))
  dir <- (-1)^stats::rbinom(sum(sel), 1, 0)   # neg.prob = 0: always +1, kept so the RNG stream matches
  f <- stats::rlnorm(sum(sel), loc, scale)
  dir[f < 1] <- -dir[f < 1]
  out <- rep(1, n); out[sel] <- f^dir
  out
}
simulate_splat_means <- function(n_genes, est) {
  base <- stats::rgamma(n_genes, shape = est$shape, rate = est$rate)
  if (!is.finite(est$out_facLoc)) return(base)
  facs <- lnorm_factors(n_genes, est$out_prob, est$out_facLoc, est$out_facScale)
  is_out <- facs != 1
  base[is_out] <- stats::median(base) * facs[is_out]
  base
}
