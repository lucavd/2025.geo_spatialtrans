# simulate_spots.R — pseudo-spot generator (PROTOCOLLO §5.2; ADEMP §2.2–2.4).
# One replicate = compositions + per-cell library sizes drawn once and shared by every arm and
# by the anchor (common random numbers, ADEMP §2.4); only μ changes between arms.

# Spot compositions: π_s ~ Dirichlet(alpha) over the spot types, then n_cells cells drawn
# multinomially (as in the Open Problems spatial-decomposition simulator). alpha is one value
# (symmetric) or one value per type, e.g. alpha_k = K * p_k with p_k the training share of type k
# (ADEMP §2.2, Pierangela 2026-10-31). Returns a spots x types matrix of cell counts; the truth is
# counts / n_cells.
draw_compositions <- function(n_spots, n_cells, types, alpha = 1) {
  k <- length(types)
  alpha <- rep_len(alpha, k)
  g <- matrix(rgamma(n_spots * k, shape = rep(alpha, each = n_spots)), n_spots, k)
  p <- g / rowSums(g)
  cnt <- t(vapply(seq_len(n_spots), function(s) as.vector(rmultinom(1, n_cells, p[s, ])), numeric(k)))
  dimnames(cnt) <- list(sprintf("spot%04d", seq_len(n_spots)), types)
  cnt
}

# One row per simulated cell: its spot, its type, and a library size drawn from the real library
# sizes of training cells of the same type (DECISIONS 2026-10-29, noise model).
draw_cells <- function(comp, L_by_type) {
  spot <- rep(rep(rownames(comp), ncol(comp)), as.vector(comp))
  type <- rep(rep(colnames(comp), each = nrow(comp)), as.vector(comp))
  L <- numeric(length(type))
  for (k in unique(type)) {
    i <- which(type == k)
    L[i] <- L_by_type[[k]][sample.int(length(L_by_type[[k]]), length(i), replace = TRUE)]
  }
  o <- order(spot)
  data.frame(spot = spot[o], type = type[o], L = L[o], stringsAsFactors = FALSE)
}

# Sum of the columns of x (genes x cells) by spot -> genes x spots, spots in `spot_names` order.
sum_by_spot <- function(x, spot, spot_names) {
  ind <- Matrix::sparseMatrix(i = seq_along(spot), j = match(spot, spot_names), x = 1,
                              dims = c(length(spot), length(spot_names)))
  out <- x %*% ind
  colnames(out) <- spot_names
  out
}

# Simulated spots for one arm: counts of gene g in a cell of type k ~ NB(mean = mu[g, k] * L,
# overdispersion phi[g]) (glmGamPoi convention: var = m + phi * m^2), summed over the cells of the
# spot. mu: genes x types (relative profile), phi: one per gene (same genes as mu).
# Done in blocks of cells to bound memory.
simulate_arm <- function(mu, phi, cells, spot_names, block = 2000) {
  stopifnot(all(cells$type %in% colnames(mu)), length(phi) == nrow(mu))
  size <- pmin(1 / pmax(phi, 1e-12), 1e12)   # phi ~ 0: NB with huge size = Poisson
  out <- Matrix::Matrix(0, nrow(mu), length(spot_names), sparse = TRUE,
                        dimnames = list(rownames(mu), spot_names))
  for (b in split(seq_len(nrow(cells)), ceiling(seq_len(nrow(cells)) / block))) {
    m <- mu[, cells$type[b], drop = FALSE] * rep(cells$L[b], each = nrow(mu))
    y <- rnbinom(length(m), size = rep(size, length(b)), mu = m)
    y <- Matrix::Matrix(matrix(y, nrow(mu)), sparse = TRUE)
    out <- out + sum_by_spot(y, cells$spot[b], spot_names)
  }
  out
}

# Real anchor: each spot sums real held-out cells with the same composition (cells drawn with
# replacement within type). held: genes x cells counts; held_type: type of each column.
anchor_spots <- function(held, held_type, comp) {
  spot_names <- rownames(comp)
  idx <- integer(0); spot <- character(0)
  for (k in colnames(comp)) {
    pool <- which(held_type == k)
    n <- comp[, k]
    pick <- pool[sample.int(length(pool), sum(n), replace = TRUE)]
    idx <- c(idx, pick); spot <- c(spot, rep(spot_names, n))
  }
  list(counts = sum_by_spot(held[, idx, drop = FALSE], spot, spot_names), cells_used = idx)
}

# Depth: a spot above the target is thinned at random (binomial, each UMI kept with probability
# target / total), as synthspot and Open Problems downsample; a spot below the target is kept.
# The target is cells per spot x UMI per cell (ADEMP §2.2), so depth and cells per spot stay crossed.
thin_to_depth <- function(x, target) {
  x <- methods::as(x, "CsparseMatrix")
  tot <- Matrix::colSums(x)
  p <- pmin(1, target / tot)
  j <- rep(seq_len(ncol(x)), diff(x@p))
  x@x <- as.numeric(rbinom(length(x@x), size = x@x, prob = p[j]))
  Matrix::drop0(x)
}

# Arm C−: the gene labels of μ are permuted within each cell type (PROTOCOLLO §5.2); column sums
# are unchanged, so library sizes are unchanged.
permute_mu <- function(mu) {
  out <- mu
  for (k in seq_len(ncol(mu))) out[, k] <- mu[sample.int(nrow(mu)), k]
  out
}
