# stream_counts.R — read selected cells (rows) of a CSR count matrix in an h5ad file window by
# window and hand each window to a callback, so that a large dataset (PsychAD brain, ~1 million
# nuclei) is never held in memory as a whole. Needs read_csr_rows() from R/read_h5ad.R.

# rows: 1-based h5ad rows to keep. FUN(m, r) receives m = cells x genes dgCMatrix of the kept
# rows of one window and r = their h5ad row numbers (same order as the rows of m).
# window: number of consecutive h5ad rows read at once (memory ~ window x genes per cell).
stream_counts <- function(file, group, rows, n_genes, FUN, window = 20000, verbose = TRUE) {
  indptr <- as.vector(rhdf5::h5read(file, paste0(group, "/indptr"), bit64conversion = "double"))
  rhdf5::H5close()
  rows <- sort(unique(rows))
  starts <- seq(min(rows), max(rows), by = window)
  for (w in seq_along(starts)) {
    first <- starts[w]
    last <- min(first + window - 1, max(rows))
    r <- rows[rows >= first & rows <= last]
    if (length(r) == 0) next
    m <- read_csr_rows(file, group, first, last, n_genes, indptr)
    FUN(m[r - first + 1, , drop = FALSE], r)
    rhdf5::H5close()
    if (verbose && w %% 10 == 0) cat(sprintf("    window %d / %d\n", w, length(starts)))
  }
  invisible(NULL)
}
