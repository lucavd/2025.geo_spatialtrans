# aim1.R — measures of Aim 1, fidelity of the cell-type profiles μ (ADEMP §5.1; DECISIONS 2026-11-11).
# Profiles are genes x types, columns = relative means (sum of counts / sum of library sizes).
# Log-profile = log(1 + 10^4 μ) (ADEMP §5.1). No random numbers are drawn here except where stated
# (noise_ceiling_draw, simulate_type_cells), which use the current RNG state.

log_profile <- function(mu) log1p(1e4 * mu)

# Pooled relative profile per type (as μ_E in 03_reference.R): x genes x cells counts, type per
# column, L = library size of each cell on ALL genes (so columns sum to the share in the kept genes).
pooled_mu <- function(x, type, L, types) {
  stopifnot(ncol(x) == length(type), length(L) == length(type), all(types %in% type))
  keep <- type %in% types
  ind <- Matrix::sparseMatrix(i = which(keep), j = match(type[keep], types), x = 1,
                              dims = c(length(type), length(types)))
  s <- as.matrix(x %*% ind)
  out <- sweep(s, 2, as.numeric(tapply(L[keep], factor(type[keep], types), sum)), "/")
  dimnames(out) <- list(rownames(x), types)
  out
}

# Pearson and Spearman correlation per type between two profiles on the same genes and types.
fidelity_cor <- function(src, tgt) {
  stopifnot(identical(rownames(src), rownames(tgt)), identical(colnames(src), colnames(tgt)))
  a <- log_profile(src); b <- log_profile(tgt)
  data.frame(cell_type = colnames(src),
             pearson = vapply(seq_len(ncol(a)), function(k) cor(a[, k], b[, k]), 0),
             spearman = vapply(seq_len(ncol(a)), function(k) cor(a[, k], b[, k], method = "spearman"), 0))
}

# log2 fold-change of each gene in each type vs the mean of the other types, on the log-profile scale
# (ADEMP §5.1: "highest log2 fold-change vs the mean of the other types, on the log(1 + 10^4 μ) scale").
lfc_matrix <- function(mu) {
  lp <- log_profile(mu)
  K <- ncol(lp)
  out <- (lp - (rowSums(lp) - lp) / (K - 1)) / log(2)
  out
}

# Top-k genes per type (ties broken by gene order, so the result is deterministic).
top_k <- function(lfc, k) lapply(seq_len(ncol(lfc)), function(j) rownames(lfc)[order(-lfc[, j])[seq_len(k)]])

# Jaccard index of the top-k markers per type, for each k in ks. src and tgt: genes x types, same
# genes and types (the type set is the same, so "the other types" are the same in both).
jaccard_topk <- function(src, tgt, ks = c(10, 20, 50, 100, 200)) {
  stopifnot(identical(rownames(src), rownames(tgt)), identical(colnames(src), colnames(tgt)))
  ls <- lfc_matrix(src); lt <- lfc_matrix(tgt)
  out <- sapply(ks, function(k) {
    a <- top_k(ls, k); b <- top_k(lt, k)
    mapply(function(x, y) length(intersect(x, y)) / length(union(x, y)), a, b)
  })
  out <- matrix(out, ncol(src), length(ks), dimnames = list(colnames(src), paste0("jaccard_k", ks)))
  data.frame(cell_type = colnames(src), out, row.names = NULL)
}

# One random split-half of the noise ceiling (ADEMP §5.1): for each type, cells split at random into
# two halves of n_k/2 (odd n_k: the extra cell to a random half), pooled half-profiles, correlation r
# (Pearson and Spearman on the log-profiles), Spearman-Brown reliability 2r/(1+r), ceiling =
# sqrt(reliability) (0 if the reliability is negative). Uses the current RNG state.
noise_ceiling_draw <- function(x, type, L, types) {
  res <- lapply(types, function(k) {
    idx <- which(type == k)
    n <- length(idx)
    h <- sample(rep(1:2, length.out = n))           # equal halves; with odd n the extra cell is random
    m1 <- Matrix::rowSums(x[, idx[h == 1], drop = FALSE]) / sum(L[idx[h == 1]])
    m2 <- Matrix::rowSums(x[, idx[h == 2], drop = FALSE]) / sum(L[idx[h == 2]])
    a <- log_profile(m1); b <- log_profile(m2)
    r <- c(pearson = cor(a, b), spearman = cor(a, b, method = "spearman"))
    sqrt(pmax(2 * r / (1 + r), 0))
  })
  data.frame(cell_type = types, ceiling_pearson = vapply(res, `[`, 0, 1),
             ceiling_spearman = vapply(res, `[`, 0, 2))
}

# ---- Highly variable genes, Seurat VST (FindVariableFeatures, selection.method = "vst", defaults:
# loess span 0.3, clip.max = sqrt(number of cells), nfeatures 2,000; Seurat source, VST.default and
# SparseRowVarStd). Computed from sums so that the training counts can be streamed (brain):
#   pass 1: per gene sum and sum of squares -> mean, variance (n - 1);
#   loess of log10(variance) on log10(mean) over genes with variance > 0 -> expected variance;
#   pass 2: per gene sum over the non-zero counts of min((x - mean)/sd_exp, clip)^2, plus the zeros
#   ((n - nnz) * (mean/sd_exp)^2), divided by n - 1 = standardized variance; genes with sd_exp = 0 get 0;
#   top nfeatures by standardized variance.
vst_expected <- function(mean, variance, span = 0.3) {
  ok <- variance > 0
  fit <- stats::loess(log10(variance[ok]) ~ log10(mean[ok]), span = span)
  ve <- numeric(length(mean))
  ve[ok] <- 10^fit$fitted
  ve
}
# m: cells x genes dgCMatrix (one window); returns per-gene sum of clipped squared z over non-zeros
# and number of non-zeros.
vst_pass2 <- function(m, mean, sd, vmax) {
  m <- methods::as(m, "CsparseMatrix")
  g <- rep(seq_len(ncol(m)), diff(m@p))
  z <- (m@x - mean[g]) / sd[g]
  z[!is.finite(z)] <- 0
  m@x <- pmin(z, vmax)^2
  list(sumz2 = Matrix::colSums(m), nnz = diff(m@p))
}
vst_finish <- function(mean, sd, n, sumz2, nnz, nfeatures = 2000, genes) {
  sv <- (sumz2 + (n - nnz) * (mean / sd)^2) / (n - 1)
  sv[sd == 0 | !is.finite(sv)] <- 0
  list(hvg = genes[order(-sv)[seq_len(nfeatures)]], std_var = setNames(sv, genes))
}
# In-memory version (genes x cells), for the simulated reference of arm P and for the test.
vst_dense <- function(x, nfeatures = 2000, span = 0.3) {
  m <- Matrix::t(methods::as(x, "CsparseMatrix"))
  n <- nrow(m)
  mean <- Matrix::colSums(m) / n
  variance <- (Matrix::colSums(m^2) - n * mean^2) / (n - 1)
  sd <- sqrt(vst_expected(mean, variance, span))
  p2 <- vst_pass2(m, mean, sd, sqrt(n))
  vst_finish(mean, sd, n, p2$sumz2, p2$nnz, nfeatures, rownames(x))
}

# ---- Separability: mean silhouette width per type (Rousseeuw 1987), Euclidean distances in the space
# of the first 50 PCs (scater::calculatePCA defaults otherwise) of log-normalised counts
# (scuttle::normalizeCounts = logNormCounts: library-size factors, log2(x/sf + 1)) on the HVGs
# (ADEMP §5.1, as in Crowell et al. 2023). x: genes x cells counts.
silhouette_by_type <- function(x, type, hvg, ncomp = 50) {
  lx <- scuttle::normalizeCounts(x)
  pc <- scater::calculatePCA(lx, ncomponents = ncomp, subset_row = intersect(hvg, rownames(x)))
  f <- factor(type)
  s <- cluster::silhouette(as.integer(f), stats::dist(pc))
  out <- tapply(s[, "sil_width"], f, mean)
  stats::setNames(as.numeric(out), names(out))
}

# Cells from a profile with the common noise model (ADEMP §2.3): n_by_type cells of each type,
# library sizes drawn with replacement from the training library sizes of the type (L_by_type).
# Uses the current RNG state. Needs simulate_cells() from R/simulate_spots.R.
simulate_type_cells <- function(mu, phi, n_by_type, L_by_type) {
  type <- rep(names(n_by_type), n_by_type)
  L <- unlist(lapply(names(n_by_type), function(k) {
    v <- L_by_type[[k]]; v[sample.int(length(v), n_by_type[[k]], replace = TRUE)]
  }), use.names = FALSE)
  list(x = simulate_cells(mu, phi, type, L), type = type)
}

# Distribution measures of arm P (ADEMP §5.1): quantiles and two-sample Kolmogorov-Smirnov distance.
dist_quantiles <- function(v, probs = c(.1, .5, .9, .99, .999))
  stats::setNames(as.numeric(stats::quantile(v, probs, names = FALSE)), paste0("q", probs * 100))
ks_d <- function(a, b) unname(suppressWarnings(stats::ks.test(a, b)$statistic))

# Number of Monte Carlo draws (DECISIONS 2026-11-11): upper 95% limit of the SD (chi-square, n - 1 df,
# as 06b_pilot_nsim.R) and n = ceiling((SD_upper / target)^2), at least n_min.
sd_upper <- function(v) { v <- v[is.finite(v)]; n <- length(v); stats::sd(v) * sqrt((n - 1) / stats::qchisq(0.025, n - 1)) }
n_draws_rule <- function(sd_up, target = 0.005, n_min = 10) max(n_min, ceiling(max(sd_up, na.rm = TRUE)^2 / target^2))
