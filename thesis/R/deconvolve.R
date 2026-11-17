# deconvolve.R — wrappers for the CPU deconvolution methods (ADEMP §4). Each takes the spots
# (genes x spots counts) and the training reference (genes x cells counts, cell types, donors)
# and returns a spots x types matrix of proportions over all reference types (rows of dropped or
# unsolved spots are NA: failures, ADEMP §7.2), plus the run time in seconds.

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- expr
  list(prop = val, seconds = proc.time()[["elapsed"]] - t0)
}

# Complete a proportion matrix to all spots and all reference types (missing -> NA rows / 0 cols).
complete_prop <- function(p, spots, types) {
  out <- matrix(NA_real_, length(spots), length(types), dimnames = list(spots, types))
  p <- as.matrix(p)
  ok_r <- intersect(rownames(p), spots); ok_c <- intersect(colnames(p), types)
  out[ok_r, ] <- 0
  out[ok_r, ok_c] <- p[ok_r, ok_c]
  out
}

# NNLS baseline: spot counts ~ profiles %*% w, profiles = mean counts per cell of each type in the
# training split; proportions = w / sum(w).
run_nnls <- function(spots, profiles) {
  timed({
    A <- as.matrix(profiles[rownames(spots), , drop = FALSE])
    p <- t(apply(as.matrix(spots), 2, function(b) nnls::nnls(A, b)$x))
    colnames(p) <- colnames(A)
    p / rowSums(p)
  })
}

# RCTD (spacexr), full mode (ADEMP §4 proposal), fake coordinates, package defaults otherwise, except
# UMI_min_sigma = 100 (amendment 2026-11-17, before any result): choose_sigma_c() estimates the
# platform-effect sigma only on spots with more than UMI_min_sigma UMI (default 300), and at 2 cells
# the HD spots are thinned to exactly 124 (lymphoid) and 242 (brain) UMI, so RCTD stopped on every
# 2-cell dataset of the HD pilots ("choose_sigma_c determined a N_fit of 0"). 100 = UMI_min, so sigma is
# estimated on the spots RCTD deconvolves; same value as Bilous et al. 2026 (Nat Methods,
# doi:10.1038/s41592-026-03089-8) for low-UMI Xenium data. Where every spot has more than 300 UMI
# (colon, 10 cells) the eligible spots are the same as with the default. umi_min_sigma = 300 reproduces
# the pilots run before the amendment.
run_rctd <- function(spots, ref, ref_type, max_cores = 1, umi_min_sigma = 100) {
  timed({
    ct <- factor(ref_type); names(ct) <- colnames(ref)
    reference <- spacexr::Reference(ref, ct, nUMI = setNames(Matrix::colSums(ref), colnames(ref)))
    puck <- spacexr::SpatialRNA(NULL, spots, nUMI = setNames(Matrix::colSums(spots), colnames(spots)),
                                use_fake_coords = TRUE)
    obj <- suppressMessages(spacexr::create.RCTD(puck, reference, max_cores = max_cores,
                                                 UMI_min_sigma = umi_min_sigma))
    obj <- suppressMessages(spacexr::run.RCTD(obj, doublet_mode = "full"))
    as.matrix(spacexr::normalize_weights(obj@results$weights))
  })
}

# MuSiC (bulk method, second baseline): spots as bulk samples, training donors as samples.
# Basis computed once (DECISIONS 2026-11-12, Luca): music_prop() recomputes music_basis() on the
# reference at every call, but the reference is the same in every replicate and arm except P.
# In MuSiC 1.0.0 music_basis(non.zero = TRUE) computes every quantity on all genes expressed in the
# reference and only at the end keeps the rows of `markers` (= genes expressed in the spots, in
# music_prop()); M.S does not depend on `markers`. So the basis on all genes, subset to those genes,
# is the same basis. music_prop_basis() is the rest of music_prop() (branch ct.cov = FALSE, the
# default we use) unchanged. Identity with music_prop() is checked before use (tests/test_music_basis.R).
music_ref_sce <- function(ref, ref_type, ref_donor) {
  # MuSiC calls counts() unqualified: attach SingleCellExperiment (not loaded on PSOCK workers)
  suppressPackageStartupMessages(library(SingleCellExperiment))
  SingleCellExperiment::SingleCellExperiment(
    list(counts = ref), colData = S4Vectors::DataFrame(type = ref_type, donor = ref_donor))
}

music_basis_once <- function(ref, ref_type, ref_donor) {
  sce <- music_ref_sce(ref, ref_type, ref_donor)
  b <- suppressMessages(MuSiC::music_basis(sce, non.zero = TRUE, markers = NULL, clusters = "type",
                                           samples = "donor", verbose = FALSE))
  b$n_ref_genes <- nrow(ref)
  b
}

# Copy of the ct.cov = FALSE path of MuSiC::music_prop() 1.0.0 (markers = NULL, cell_size = NULL,
# default iter.max, nu, eps, centered, normalize), starting from a precomputed basis.
music_prop_basis <- function(bulk.mtx, basis, iter.max = 1000, nu = 0.0001, eps = 0.01,
                             centered = FALSE, normalize = FALSE) {
  bulk.gene <- rownames(bulk.mtx)[rowMeans(bulk.mtx) != 0]
  bulk.mtx <- bulk.mtx[bulk.gene, ]
  # music_prop: rows of the basis = intersect(bulk.gene, reference genes), in bulk.gene order
  cm.gene <- intersect(bulk.gene, rownames(basis$Disgn.mtx))
  if (length(cm.gene) < 0.2 * min(length(bulk.gene), basis$n_ref_genes)) stop("Too few common genes!")
  m.sc <- match(cm.gene, rownames(basis$Disgn.mtx)); m.bulk <- match(cm.gene, bulk.gene)
  D1 <- basis$Disgn.mtx[m.sc, ]
  M.S <- colMeans(basis$S, na.rm = TRUE)
  Yjg <- MuSiC:::relative.ab(bulk.mtx[m.bulk, ]); N.bulk <- ncol(bulk.mtx)
  Sigma <- basis$Sigma[m.sc, ]
  valid.ct <- (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) == 0) & (!is.na(M.S))
  if (sum(valid.ct) <= 1) stop("Not enough valid cell type!")
  D1 <- D1[, valid.ct]; M.S <- M.S[valid.ct]; Sigma <- Sigma[, valid.ct]
  Est.prop.weighted <- NULL
  for (i in 1:N.bulk) {
    nz <- Yjg[, i] != 0
    if (sum(!nz) > 0) {
      D1.temp <- D1[nz, ]; Yjg.temp <- Yjg[nz, i]; Sigma.temp <- Sigma[nz, ]
    } else {
      D1.temp <- D1; Yjg.temp <- Yjg[, i]; Sigma.temp <- Sigma
    }
    lm.D1.weighted <- MuSiC:::music.iter(Yjg.temp, D1.temp, M.S, Sigma.temp, iter.max = iter.max,
                                         nu = nu, eps = eps, centered = centered, normalize = normalize)
    Est.prop.weighted <- rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
  }
  colnames(Est.prop.weighted) <- colnames(D1)
  rownames(Est.prop.weighted) <- colnames(Yjg)
  Est.prop.weighted
}

# basis = NULL: original music_prop() call (pilots); otherwise the precomputed basis.
run_music <- function(spots, ref, ref_type, ref_donor, basis = NULL) {
  timed({
    if (is.null(basis)) {
      sce <- music_ref_sce(ref, ref_type, ref_donor)
      est <- suppressMessages(MuSiC::music_prop(bulk.mtx = as.matrix(spots), sc.sce = sce,
                                                clusters = "type", samples = "donor", verbose = FALSE))
      est$Est.prop.weighted
    } else {
      music_prop_basis(as.matrix(spots), basis)
    }
  })
}

# spatialDWLS (standalone package). Settings (ADEMP §4):
#  - normalisation: counts / total * 1e4; spots log2(x + 1), signature on the unlogged scale;
#  - markers: top n_markers genes per type by log2 fold-change of the mean normalised expression
#    (type vs mean of the other types), from the training reference;
#  - spot clusters (changed 2026-11-04, ADEMP §4): Leiden clustering with the settings of Spotless
#    (saeyslab/spotless-benchmark, subworkflows/deconvolution/spatialdwls/script_nf.R: nn.dims = 10,
#    nn.k = 4, cluster.res = 0.4, cluster.n_iter = 1000) on a shared-nearest-neighbour graph
#    (bluster::makeSNNGraph, default "rank" weights) of the first 10 principal components of the
#    log-normalised spots over the marker genes; igraph instead of Giotto. The colon pilot used
#    k-means (k = 10, 20 PCs), which had no source. Leiden uses the R RNG, so the replicate seed fixes it;
#  - n_cell = the true cells-per-spot level.
run_dwls <- function(spots, ref, ref_type, n_cell, n_markers = 100, nn_dims = 10, nn_k = 4,
                     cluster_res = 0.4, cluster_n_iter = 1000) {
  timed({
    norm <- function(x) { x <- methods::as(x, "CsparseMatrix"); Matrix::t(Matrix::t(x) / Matrix::colSums(x)) * 1e4 }
    ref_n <- norm(ref)
    types <- sort(unique(ref_type))
    ind <- Matrix::sparseMatrix(i = seq_along(ref_type), j = match(ref_type, types), x = 1)
    means <- as.matrix(ref_n %*% ind) / rep(tabulate(match(ref_type, types)), each = nrow(ref))
    colnames(means) <- types
    lfc <- vapply(seq_along(types), function(k)
      log2(means[, k] + 0.1) - log2(rowMeans(means[, -k, drop = FALSE]) + 0.1), numeric(nrow(means)))
    markers <- unique(unlist(lapply(seq_along(types), function(k)
      rownames(means)[order(lfc[, k], decreasing = TRUE)[seq_len(n_markers)]])))
    sig <- spatialDWLS::dwls_signature(ref_n, markers, ref_type)
    expr <- log2(norm(spots) + 1)
    pcs <- stats::prcomp(t(as.matrix(expr[markers, ])), rank. = nn_dims)$x
    g <- bluster::makeSNNGraph(pcs, k = nn_k)
    cl <- igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity",
                                                    resolution = cluster_res, n_iterations = cluster_n_iter))
    p <- spatialDWLS::dwls(expr, sig, clusters = cl, n_cell = n_cell)
    t(p)
  })
}
