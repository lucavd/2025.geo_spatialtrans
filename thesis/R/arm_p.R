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

# ---------------------------------------------------------------------------------------------------
# Arm P in one replicate (DECISIONS 2026-11-04, Luca: P as a coherent de novo world; Pierangela:
# markers "some", φ by rank). Pre-registered settings, written 2026-11-05 before the first run:
#   - μ^P: n_genes base means from simulate_splat_means() with the tissue calibration
#     (04p_calibrate_p.R); n_genes = genes of the reference the run uses. Genes are named P00001...,
#     because they are not the real genes. Markers: each type k gets as many markers as found for
#     it in training, drawn at random among all P genes (without replacement within a type; a gene
#     drawn for two types gets both fold-changes, one per column); each fold-change is drawn with
#     replacement from the empirical fold-changes of the training markers of the SAME type
#     (Pierangela, 2026-11-05: per type, so each type keeps the strength of its own markers, e.g.
#     Plasma). Columns normalised to 1 (ADEMP §2.1).
#   - φ of the P genes (phi_by_rank): the P gene with the r-th smallest base mean takes the φ of the
#     real gene with the r-th smallest pooled training mean; ties at random. The real pooled mean is
#     the pooled training profile (pooled_mean(): Σ_k μ_E,gk T_k / Σ_k T_k, T_k = total counts of
#     the training cells of type k, which is exact because μ_E,gk = counts_gk / T_k).
#   - Method reference of P (simulate_reference_p): one simulated cell for each cell of the real
#     capped reference (§4), same cell type, same donor label, library size = total counts of that
#     real cell over the kept genes (the saved reference holds only kept genes; in colon the kept
#     genes carry 99.7–99.95% of the counts, column sums of μ_E: declared); counts from the common
#     noise model, NB(mean μ^P_gk × L, φ^P_g). NNLS profiles of P = mean counts per simulated
#     reference cell of each type (the real arms use all training cells; P has only the simulated
#     reference: declared). A new μ^P, φ^P and reference in every replicate.
#   - Order of random numbers inside a replicate: μ^P, then φ^P ties, then the reference, then the
#     spots of P (so the pilot streams of E, C− and the anchor are untouched if P is drawn last).
# Checks: tests/test_arm_p_o.R (pre-registered in its header).

# Pooled training mean of every real gene, from μ_E (genes x types) and the training library
# sizes by type (L_by_type, all training cells).
pooled_mean <- function(mu_E, L_by_type) {
  tot <- vapply(L_by_type[colnames(mu_E)], sum, numeric(1))
  as.vector(mu_E %*% tot) / sum(tot)
}

# μ^P for one replicate. calib = <tissue>_p_calibration.rds. Returns list(mu = genes x types,
# columns summing to 1; base = base means before the markers; markers = list of P genes per type).
simulate_mu_p <- function(calib, n_genes, types = calib$types) {
  stopifnot(all(types %in% names(calib$markers)))
  gid <- sprintf("P%05d", seq_len(n_genes))
  base <- stats::setNames(simulate_splat_means(n_genes, calib$baseline), gid)
  mu <- matrix(base, n_genes, length(types), dimnames = list(gid, types))
  mk <- stats::setNames(vector("list", length(types)), types)
  for (k in types) {
    n_k <- length(calib$markers[[k]])
    if (!n_k) { mk[[k]] <- character(0); next }
    fc <- calib$fold_changes[[k]]
    g <- sample.int(n_genes, n_k)
    mu[g, k] <- mu[g, k] * fc[sample.int(length(fc), n_k, replace = TRUE)]
    mk[[k]] <- gid[g]
  }
  list(mu = sweep(mu, 2, colSums(mu), "/"), base = base, markers = mk)
}

# φ for the P genes by rank of the mean (ADEMP §2.1). phi_real and mean_real: one per real gene;
# mean_p: base means of the P genes (same length). Ties broken at random.
phi_by_rank <- function(phi_real, mean_real, mean_p) {
  stopifnot(length(phi_real) == length(mean_real), length(mean_p) == length(phi_real))
  phi_sorted <- unname(phi_real)[order(rank(mean_real, ties.method = "random"))]
  stats::setNames(phi_sorted[rank(mean_p, ties.method = "random")], names(mean_p))
}

# Method reference of P: ref_meta = the real reference's cells (cell_type, donor_id), L = their
# library sizes. Returns list(ref = genes x cells counts, ref_meta, profiles = genes x types mean
# counts per cell), the same fields the methods read from the real reference.
simulate_reference_p <- function(mu_p, phi_p, ref_meta, L, cell_names = NULL, block = 1000) {
  stopifnot(length(L) == nrow(ref_meta), all(ref_meta$cell_type %in% colnames(mu_p)))
  x <- simulate_cells(mu_p, phi_p, ref_meta$cell_type, L, block)
  if (!is.null(cell_names)) colnames(x) <- cell_names
  types <- colnames(mu_p)
  ind <- Matrix::sparseMatrix(i = seq_len(ncol(x)), j = match(ref_meta$cell_type, types), x = 1,
                              dims = c(ncol(x), length(types)))
  n <- tabulate(match(ref_meta$cell_type, types), length(types))
  profiles <- sweep(as.matrix(x %*% ind), 2, n, "/")
  dimnames(profiles) <- list(rownames(mu_p), types)
  list(ref = x, ref_meta = ref_meta, profiles = profiles)
}

# Method reference of P rebuilt from a saved world (scripts/04_simulate.R, world_p.rds: mu, phi and
# ref_seed = the L'Ecuyer-CMRG state of the reference sub-stream). Same inputs -> identical reference,
# so the reference is not stored. L = Matrix::colSums(ref) of the real capped reference (ADEMP §2.1).
reference_p_from_world <- function(world, ref_meta, L, cell_names = NULL) {
  assign(".Random.seed", world$ref_seed, envir = .GlobalEnv)
  simulate_reference_p(world$mu, world$phi, ref_meta, L, cell_names)
}
