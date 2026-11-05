# S2.2, analisi post hoc PH-4 (dichiarata dopo G-2): la tolleranza «|z| < 5 per ogni gene» di G-2 (b) non tiene conto del numero di confronti
# (da 7 sottotipi x migliaia di geni a 84 sottotipi x migliaia di geni). Qui si guarda la distribuzione degli z dei profili ricostruiti
# (modello M) dai file del pannello: media, deviazione standard, numero di |z| > 5 contro l'atteso sotto N(0, 1), soglia di Bonferroni
# per un errore di famiglia di 0.001. Uso: Rscript --vanilla tools/S2.2_posthoc_ph4.R -> results/S2.2/S2.2_PH4.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
D <- "/mnt/micron/geo_spatialtrans/S2.2/panel"
out <- list()
for (A in paste0("A", 1:6)) for (sd in c(20261108L, 20261109L)) {
  prof <- load_cell_profiles(A); fs <- list.files(D, pattern = sprintf("^%s_.*_s%d_M\\.rds$", A, sd), full.names = TRUE)
  sums <- 0; nsub <- 0
  for (p in fs) { s <- readRDS(p)
    agg <- Matrix::sparseMatrix(i = seq_len(nrow(s$cells)), j = match(s$cells$subtype, prof$subtypes$subtype), x = 1, dims = c(nrow(s$cells), nrow(prof$subtypes)))
    sums <- sums + as.matrix(s$counts %*% agg); nsub <- nsub + Matrix::colSums(agg) }
  z <- c()
  for (i in which(nsub >= 200)) { Tt <- sum(sums[, i]); ex <- Tt * prof$P[i, ]; k <- ex >= 20; z <- c(z, (sums[k, i] - ex[k]) / sqrt(ex[k] * (1 - prof$P[i, k]))) }
  n <- length(z); k5 <- sum(abs(z) > 5); e5 <- n * 2 * pnorm(-5); zb <- qnorm(1 - 0.001 / (2 * n))
  out[[length(out) + 1]] <- data.frame(archetype = A, seed = sd, n_z = n, mean_z = mean(z), sd_z = sd(z), max_abs_z = max(abs(z)),
    n_abs_z_gt5 = k5, expected_gt5 = e5, p_poisson_ge = ppois(k5 - 1, e5, lower.tail = FALSE), bonferroni_z_0.001 = zb, pass_bonferroni = max(abs(z)) < zb,
    frac_abs_z_gt2 = mean(abs(z) > 2), expected_frac_gt2 = 2 * pnorm(-2))
}
tab <- do.call(rbind, out); write.csv(tab, "results/S2.2/S2.2_PH4.csv", row.names = FALSE); print(tab, digits = 3)
