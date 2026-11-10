# tools/S2.6_check.R - S2.6, G-2: invarianti di ogni esecuzione di tools/S2.6_pieces.R. Uso: Rscript --vanilla tools/S2.6_check.R <dir1> [<dir2> ...]
# Scrive una riga CSV per esecuzione su stdout; esce con 1 se un invariante fallisce.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(Matrix))
dirs <- commandArgs(TRUE); bad <- 0L
cat("tag,mode,n_bins,n_na,u2_eq_c8,cover_once,cons_pieces,n_empty,sum_u2,ok\n")
for (d in dirs) {
  j <- jsonlite::fromJSON(file.path(d, "run.json"), simplifyVector = FALSE); g <- j$grid; nb <- g$n_row * g$n_col
  u2 <- readBin(file.path(d, "u2.bin"), "integer", nb + 1L, size = 4L); c8 <- readRDS(file.path(d, "c8.rds"))
  n8 <- (g$n_row %/% 4L) * (g$n_col %/% 4L); id <- seq_len(nb)
  b8 <- as.integer((((id - 1L) %/% g$n_col) %/% 4L) * (g$n_col %/% 4L) + ((id - 1L) %% g$n_col) %/% 4L + 1L)
  nna <- sum(is.na(u2)) + sum(is.na(c8@x))
  eq <- length(u2) == nb && ncol(c8) == n8 && all(abs(Matrix::colSums(c8) - tabulate(rep(b8, u2), n8)) == 0)
  cov1 <- TRUE; cons <- TRUE; nemp <- 0L
  if (identical(j$mode, "one")) {
    cons <- isTRUE(abs(sum(as.numeric(u2)) - (j$n_captured + j$n_background)) == 0)
  } else {
    cov1 <- isTRUE(j$cover_max == 1L)
    for (p in j$pieces) { if (isTRUE(p$empty)) { nemp <- nemp + 1L; next }; if (!isTRUE(p$cons_ok)) cons <- FALSE }
  }
  ok <- nna == 0 && eq && cov1 && cons; if (!ok) bad <- bad + 1L
  cat(sprintf("%s,%s,%d,%d,%s,%s,%s,%d,%.0f,%s\n", j$tag, j$mode, nb, nna, eq, cov1, cons, nemp, sum(as.numeric(u2)), ok))
}
quit(status = if (bad == 0L) 0L else 1L)
