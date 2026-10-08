# tools/S1.2_fixcheck.R — confronto prima/dopo le correzioni della revisione avversariale (S1.2)
# prima: /mnt/micron/geo_spatialtrans/S1.2/prefix/ (sim/ prodotte da 6316790; CSV estratti con git show 6316790:...); dopo: checkB/sim/ e results/S1.2/
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
P <- "/mnt/micron/geo_spatialtrans/S1.2/prefix"; N <- "/mnt/micron/geo_spatialtrans/S1.2/checkB/sim"
rows <- list()
for (f in list.files(file.path(P, "sim"))) {
  a <- readRDS(file.path(P, "sim", f)); b <- readRDS(file.path(N, f))
  ga <- as.matrix(a[, grep("^g_", names(a))]); gb <- as.matrix(b[, grep("^g_", names(b))])
  rows[[f]] <- data.frame(unit = f, n_rows = nrow(a), same_rows = nrow(a) == nrow(b),
    same_reps = identical(a$reps, b$reps), same_feasible = identical(a$feasible, b$feasible),
    same_n = isTRUE(all.equal(a$n_mean, b$n_mean)), max_abs_dg = suppressWarnings(max(abs(ga - gb), na.rm = TRUE)),
    identical_g = identical(ga, gb))
}
X <- do.call(rbind, rows)
ta <- read.csv(file.path(P, "S1.2_test_results.csv")); tb <- read.csv("results/S1.2/S1.2_test_results.csv")
k <- merge(ta, tb, by = c("check", "input", "metric"), suffixes = c("_prima", "_dopo"))
b1a <- read.csv(file.path(P, "S1.2_B1b.csv")); b1b <- read.csv("results/S1.2/S1.2_B1b.csv")
out <- data.frame(
  verifica = c("unita' di simulazione confrontate", "righe identiche (g bit per bit)", "max |Δg| su tutte le righe",
               "reps e fattibilita' identiche", "asserzioni check C comuni", "con valore identico", "B1b densita' seed 42 identica (30 ROI)"),
  valore = c(nrow(X), sprintf("%d/%d", sum(X$identical_g), nrow(X)), sprintf("%.3g", max(X$max_abs_dg)),
             sprintf("%d/%d", sum(X$same_reps & X$same_feasible), nrow(X)), nrow(k), sprintf("%d/%d", sum(k$value_prima == k$value_dopo), nrow(k)),
             sprintf("%d/%d", sum(b1a$dens_seed42 == b1b$dens_seed42), nrow(b1a))))
write.csv(out, "results/S1.2/S1.2_fixcheck.csv", row.names = FALSE)
write.csv(k[k$value_prima != k$value_dopo, ], "results/S1.2/S1.2_fixcheck_testdiff.csv", row.names = FALSE)
print(out); print(k[k$value_prima != k$value_dopo, c("check", "input", "metric", "value_prima", "value_dopo")])
