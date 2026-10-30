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

# RCTD (spacexr), full mode (ADEMP §4 proposal), fake coordinates, package defaults otherwise.
run_rctd <- function(spots, ref, ref_type, max_cores = 1) {
  timed({
    ct <- factor(ref_type); names(ct) <- colnames(ref)
    reference <- spacexr::Reference(ref, ct, nUMI = setNames(Matrix::colSums(ref), colnames(ref)))
    puck <- spacexr::SpatialRNA(NULL, spots, nUMI = setNames(Matrix::colSums(spots), colnames(spots)),
                                use_fake_coords = TRUE)
    obj <- suppressMessages(spacexr::create.RCTD(puck, reference, max_cores = max_cores))
    obj <- suppressMessages(spacexr::run.RCTD(obj, doublet_mode = "full"))
    as.matrix(spacexr::normalize_weights(obj@results$weights))
  })
}

# MuSiC (bulk method, second baseline): spots as bulk samples, training donors as samples.
run_music <- function(spots, ref, ref_type, ref_donor) {
  timed({
    sce <- SingleCellExperiment::SingleCellExperiment(
      list(counts = ref), colData = S4Vectors::DataFrame(type = ref_type, donor = ref_donor))
    est <- suppressMessages(MuSiC::music_prop(bulk.mtx = as.matrix(spots), sc.sce = sce,
                                              clusters = "type", samples = "donor", verbose = FALSE))
    est$Est.prop.weighted
  })
}

# spatialDWLS (standalone package). Pilot settings (ADEMP §4, [PILOT]):
#  - normalisation: counts / total * 1e4; spots log2(x + 1), signature on the unlogged scale;
#  - markers: top n_markers genes per type by log2 fold-change of the mean normalised expression
#    (type vs mean of the other types), from the training reference;
#  - spot clusters: k-means (k_clusters, nstart 10) on the first 20 principal components of the
#    log-normalised spots over the marker genes;
#  - n_cell = the true cells-per-spot level.
run_dwls <- function(spots, ref, ref_type, n_cell, n_markers = 100, k_clusters = 10) {
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
    pcs <- stats::prcomp(t(as.matrix(expr[markers, ])), rank. = 20)$x
    cl <- stats::kmeans(pcs, centers = min(k_clusters, ncol(spots) - 1), nstart = 10)$cluster
    p <- spatialDWLS::dwls(expr, sig, clusters = cl, n_cell = n_cell)
    t(p)
  })
}
