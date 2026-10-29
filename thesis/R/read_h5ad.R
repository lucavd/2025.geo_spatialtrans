# read_h5ad.R — read cell metadata (obs) and detected genes per cell from an h5ad file,
# without loading the count matrix (PsychAD has 1.3 million nuclei).
# Needs rhdf5 (Bioconductor).

# Columns of obs as plain R vectors. Categorical columns (codes + categories) become
# character; nullable columns (values + mask) get NA where masked.
read_obs_columns <- function(file, cols) {
  info <- rhdf5::h5ls(file, recursive = 3)
  read_one <- function(col) {
    p <- paste0("/obs/", col)
    kids <- info$name[info$group == p]
    if ("codes" %in% kids) {
      codes <- as.vector(rhdf5::h5read(file, paste0(p, "/codes")))
      cats <- as.vector(rhdf5::h5read(file, paste0(p, "/categories")))
      out <- cats[codes + 1L]
      out[codes < 0] <- NA
      out
    } else if ("values" %in% kids) {
      v <- as.vector(rhdf5::h5read(file, paste0(p, "/values")))
      if ("mask" %in% kids) v[as.logical(rhdf5::h5read(file, paste0(p, "/mask")))] <- NA
      v
    } else if (any(info$group == "/obs" & info$name == col)) {
      as.vector(rhdf5::h5read(file, p))
    } else {
      stop("Column '", col, "' not found in obs of ", basename(file))
    }
  }
  index_name <- rhdf5::h5readAttributes(file, "/obs")[["_index"]]
  out <- data.frame(cell_id = as.vector(rhdf5::h5read(file, paste0("/obs/", index_name))),
                    stringsAsFactors = FALSE)
  for (col in cols) out[[col]] <- read_one(col)
  rhdf5::H5close()
  out
}

# Gene ids (var index, Ensembl in CELLxGENE) and symbols (var/feature_name) for the count
# matrix in `group`: "/X" uses /var, "/raw/X" uses /raw/var. Symbols come from /var, which
# CELLxGENE keeps in the same order as /raw/var; the lengths are checked.
read_var_genes <- function(file, group) {
  var <- if (group == "/raw/X") "/raw/var" else "/var"
  index_name <- rhdf5::h5readAttributes(file, var)[["_index"]]
  ids <- as.vector(rhdf5::h5read(file, paste0(var, "/", index_name)))
  main_ids <- as.vector(rhdf5::h5read(file, paste0("/var/", rhdf5::h5readAttributes(file, "/var")[["_index"]])))
  codes <- as.vector(rhdf5::h5read(file, "/var/feature_name/codes"))
  cats <- as.vector(rhdf5::h5read(file, "/var/feature_name/categories"))
  rhdf5::H5close()
  symbols <- setNames(cats[codes + 1L], main_ids)[ids]
  data.frame(gene_id = ids, symbol = unname(symbols), stringsAsFactors = FALSE)
}

# Rows first..last (1-based, contiguous) of a CSR count matrix as a cells x genes dgCMatrix.
# Only the slice of data/indices for those rows is read. Positions are doubles because brain
# has more than 2^31 non-zeros.
read_csr_rows <- function(file, group, first, last, n_genes, indptr) {
  from <- indptr[first]                 # 0-based offset of the first entry of row `first`
  n <- indptr[last + 1] - from
  x <- as.vector(rhdf5::h5read(file, paste0(group, "/data"), start = from + 1, count = n))
  j <- as.vector(rhdf5::h5read(file, paste0(group, "/indices"), start = from + 1, count = n))
  per_row <- diff(indptr[first:(last + 1)])
  Matrix::sparseMatrix(i = rep.int(seq_along(per_row), per_row), j = j + 1L, x = as.numeric(x),
                       dims = c(length(per_row), n_genes))
}

# Number of detected genes (non-zero entries) per cell, from the row pointer of a CSR
# matrix ("/X" or "/raw/X"): cell i has indptr[i + 1] - indptr[i] non-zero genes.
# indptr is read as double because the total number of non-zeros can exceed 2^31.
n_genes_per_cell <- function(file, group) {
  enc <- rhdf5::h5readAttributes(file, group)[["encoding-type"]]
  if (!identical(enc, "csr_matrix")) stop(group, " is '", enc, "', expected csr_matrix (cells in rows)")
  indptr <- as.vector(rhdf5::h5read(file, paste0(group, "/indptr"), bit64conversion = "double"))
  rhdf5::H5close()
  diff(indptr)
}
