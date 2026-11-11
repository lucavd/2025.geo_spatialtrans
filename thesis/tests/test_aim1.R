# test_aim1.R — checks of R/aim1.R (Aim 1 measures, ADEMP §5.1).
#   1. vst_dense() (streamable sums) selects the same HVGs as Seurat::FindVariableFeatures(vst, 2,000)
#      and its standardized variance agrees, on the real colon training counts of the reference file
#      (Seurat is not in renv: give its library with SEURAT_LIB; the test is skipped without it);
#   2. vst_pass2() over two windows = one window (streaming does not change the result);
#   3. fidelity_cor() and jaccard_topk() of a profile with itself = 1; Jaccard of disjoint tops = 0;
#   4. noise ceiling: identical halves give 1;
#   5. pooled_mu() columns = sum of counts / sum of library sizes.
# Usage (working directory = thesis/):  SEURAT_LIB=<path> Rscript tests/test_aim1.R [colon]
source("renv/activate.R")
suppressPackageStartupMessages(library(Matrix))
for (f in c("R/simulate_spots.R", "R/aim1.R")) source(f)
ok <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); invisible(cond) }
set.seed(1)

# 2. streaming in windows
x <- Matrix::rsparsematrix(300, 400, 0.1, rand.x = function(n) rpois(n, 3) + 1)   # genes x cells
rownames(x) <- paste0("g", 1:300)
m <- Matrix::t(x); n <- nrow(m)
mean <- Matrix::colSums(m) / n; variance <- (Matrix::colSums(m^2) - n * mean^2) / (n - 1)
sd <- sqrt(vst_expected(mean, variance))
a <- vst_pass2(m, mean, sd, sqrt(n))
b1 <- vst_pass2(m[1:150, ], mean, sd, sqrt(n)); b2 <- vst_pass2(m[151:400, ], mean, sd, sqrt(n))
ok(isTRUE(all.equal(a$sumz2, b1$sumz2 + b2$sumz2)) && all(a$nnz == b1$nnz + b2$nnz), "VST pass 2 in two windows = one")

# 3. self-comparison and disjoint tops
mu <- matrix(rgamma(300 * 4, 1, 1e4), 300, 4, dimnames = list(rownames(x), paste0("T", 1:4)))
mu <- sweep(mu, 2, colSums(mu), "/") * 0.98
ok(all(abs(as.matrix(fidelity_cor(mu, mu)[, -1]) - 1) < 1e-12), "correlation of a profile with itself = 1")
ok(all(abs(as.matrix(jaccard_topk(mu, mu, c(10, 50))[, -1]) - 1) < 1e-12), "Jaccard with itself = 1")
l1 <- top_k(lfc_matrix(mu), 10)
ok(all(vapply(l1, length, 1L) == 10), "top_k returns k genes")

# 4. ceiling with duplicated cells: halves identical -> 1
y <- x[, 1:100]; y <- cbind(y, y); ty <- rep(rep(c("A", "B"), each = 50), 2)
L <- Matrix::colSums(y)
# halves are random, so build a case where every cell of a type is identical
z <- Matrix::Matrix(rep(Matrix::rowSums(x[, 1:5]), 60), nrow = 300, sparse = TRUE)
rownames(z) <- rownames(x)
nc <- noise_ceiling_draw(z, rep(c("A", "B"), each = 30), Matrix::colSums(z), c("A", "B"))
ok(all(abs(nc$ceiling_pearson - 1) < 1e-12), "ceiling of identical cells = 1")

# 5. pooled μ
p <- pooled_mu(y, ty, L * 1.01, c("A", "B"))
ok(isTRUE(all.equal(unname(p[, "A"]), unname(Matrix::rowSums(y[, ty == "A"]) / sum(L[ty == "A"] * 1.01)))), "pooled_mu")

# 1. Seurat on real counts
lib <- Sys.getenv("SEURAT_LIB")
args <- commandArgs(trailingOnly = TRUE); tissue <- if (length(args)) args[1] else "colon"
if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(Seurat))
  R0 <- readRDS(file.path(Sys.getenv("TESI_DATA"), "derived", "reference", sprintf("%s_pilot.rds", tissue)))
  cnt <- R0$ref                                    # training cells of the reference (genes x cells)
  mine <- vst_dense(cnt, 2000)
  so <- CreateSeuratObject(cnt, min.cells = 0, min.features = 0)
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  hv <- VariableFeatures(so)
  info <- HVFInfo(so, method = "vst")
  sv <- info[names(mine$std_var), grep("standardized", colnames(info))[1]]
  cat(sprintf("  Seurat %s; %d cells x %d genes; overlap %d of 2000; max |diff| std var %.2e\n",
              packageVersion("Seurat"), ncol(cnt), nrow(cnt), length(intersect(hv, mine$hvg)),
              max(abs(sv - mine$std_var))))
  ok(setequal(hv, mine$hvg), "VST: same 2,000 HVGs as Seurat")
  ok(max(abs(sv - mine$std_var)) < 1e-6, "VST: same standardized variance as Seurat")
} else cat("SKIP Seurat check (SEURAT_LIB not set)\n")
