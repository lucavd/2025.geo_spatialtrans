# test_simulate.R — tests of the pseudo-spot generator (R/simulate_spots.R) and of the metrics
# (R/pilot_metrics.R) on a toy example. Run from thesis/:  Rscript tests/test_simulate.R
source("R/simulate_spots.R"); source("R/pilot_metrics.R")
ok <- function(cond, msg) cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n")
set.seed(1)
G <- 200; types <- c("A", "B", "C")
mu <- matrix(rgamma(G * 3, 1), G, 3, dimnames = list(paste0("g", 1:G), types)); mu <- sweep(mu, 2, colSums(mu), "/")
phi <- rep(0.1, G); L <- list(A = rep(1000, 50), B = rep(2000, 50), C = rep(3000, 50))

comp <- draw_compositions(100, 5, types)
ok(all(rowSums(comp) == 5) && all(comp >= 0), "compositions: 5 cells per spot")
cells <- draw_cells(comp, L)
ok(nrow(cells) == 500 && all(table(factor(cells$type, types)) == colSums(comp)), "cells match the compositions")
x <- simulate_arm(mu, phi, cells, rownames(comp))
ok(all(dim(x) == c(G, 100)), "simulated spots: genes x spots")
expected <- as.vector(comp %*% c(1000, 2000, 3000))
ok(abs(mean(Matrix::colSums(x) / expected) - 1) < 0.05, "mean spot total = sum of library sizes (within 5%)")
y <- thin_to_depth(x, 4000)
ok(all(Matrix::colSums(y) <= Matrix::colSums(x)) && all(Matrix::colSums(y)[Matrix::colSums(x) <= 4000] == Matrix::colSums(x)[Matrix::colSums(x) <= 4000]),
   "thinning only lowers totals and keeps spots below the target")
ok(all(Matrix::colSums(y)[Matrix::colSums(x) > 4000] == 4000), "thinned spots exactly at the target")
ok(all(y <= x), "thinning never adds counts to a gene")
# Proportions kept on average (DECISIONS 2026-11-02): thinning one spot R times, the mean count of
# each gene must match x_g * t / T within 5 hypergeometric SEs of the mean
j <- which.max(Matrix::colSums(x)); xs <- x[, j, drop = FALSE]; Tt <- sum(xs); tg <- 4000; R <- 500
yb <- rowMeans(vapply(seq_len(R), function(r) as.vector(thin_to_depth(xs, tg)), numeric(G)))
p <- as.vector(xs) / Tt; se <- sqrt(tg * p * (1 - p) * (Tt - tg) / (Tt - 1) / R)
ok(all(abs(yb - tg * p)[p > 0] <= 5 * se[p > 0]) && all(yb[p == 0] == 0), "thinning keeps gene proportions on average")
mc <- permute_mu(mu)
ok(all(abs(colSums(mc) - colSums(mu)) < 1e-12) && !identical(mc, mu), "C−: permuted, column sums unchanged")
held <- Matrix::Matrix(matrix(rpois(G * 30, 2), G, 30), sparse = TRUE); ht <- rep(types, each = 10)
a <- anchor_spots(held, ht, comp)
ok(all(Matrix::colSums(a$counts) > 0) && length(a$cells_used) == 500 && all(table(factor(ht[a$cells_used], types)) == colSums(comp)),
   "anchor: real cells drawn with the spot compositions")
tr <- comp / 5
ok(rmse(tr, tr) == 0 && jsd(tr, tr) == 0, "metrics are 0 on the truth")
ok(abs(jsd(rbind(c(1, 0)), rbind(c(0, 1))) - 1) < 1e-12, "JSD of disjoint compositions = 1 (base 2)")
e <- tr; e[1, ] <- NA; e[2, ] <- e[2, ] * 2
ok(identical(unname(which(failed_spots(e))), 1:2), "failures: NA rows and sums != 1")
set.seed(2)
p <- c(A = 0.8, B = 0.15, C = 0.05)
cv <- draw_compositions(4000, 10, types, alpha = 3 * p)   # alpha_k = K * p_k
ok(all(abs(colMeans(cv / 10) - p) < 0.02), "vector alpha: mean composition = training shares")
# UMI-share truth (2026-11-03): rows sum to 1; equal library sizes give the cell-share truth;
# the anchor's spot vector matches its cells; a type with larger cells gets a larger share
tu <- umi_truth(cells$spot, cells$type, cells$L, rownames(comp), types)
eqL <- umi_truth(cells$spot, cells$type, rep(1, nrow(cells)), rownames(comp), types)
ok(all(abs(rowSums(tu) - 1) < 1e-12) && max(abs(eqL - comp / rowSums(comp))) < 1e-12,
   "UMI truth: rows sum to 1, equal L gives the cell truth")
ok(length(a$spot) == length(a$cells_used) && all(table(factor(a$spot, rownames(comp))) == rowSums(comp)),
   "anchor: one spot label per cell used")
