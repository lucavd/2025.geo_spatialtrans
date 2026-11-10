.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
a <- commandArgs(TRUE); d1 <- a[1]; d2 <- a[2]
u1 <- readBin(file.path(d1, "u2.bin"), "integer", 1e8, size = 4L); u2 <- readBin(file.path(d2, "u2.bin"), "integer", 1e8, size = 4L)
j2 <- jsonlite::fromJSON(file.path(d2, "run.json"))
c1 <- readRDS(file.path(d1, "c8.rds")); c2 <- readRDS(file.path(d2, "c8.rds"))
cov <- u2 > 0 | TRUE
cat("len", length(u1), length(u2), " cover_max", j2$cover_max, " bins covered", j2$n_bins_covered, "\n")
g <- j2$grid; ii <- (seq_along(u1) - 1L) %/% g$n_col; jj <- (seq_along(u1) - 1L) %% g$n_col
x <- g$x0 + (jj + 0.5) * 2; y <- g$y0 + (ii + 0.5) * 2
U <- j2$pieces$U; if (is.list(U)) U <- do.call(rbind, U)
inU <- rep(FALSE, length(u1)); for (k in seq_len(nrow(U))) inU <- inU | (x >= U[k,1] & x < U[k,3] & y >= U[k,2] & y < U[k,4])
cat("bins in U", sum(inU), " identical u2 in U:", identical(u1[inU], u2[inU]), " sum one", sum(u1[inU]), " sum pieces", sum(u2[inU]), " outside U pieces", sum(u2[!inU]), "\n")
