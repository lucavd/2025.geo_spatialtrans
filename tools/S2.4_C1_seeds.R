# tools/S2.4_C1_seeds.R - S2.4, correzione C-1: verifica dei flussi L'Ecuyer-CMRG con semi consecutivi (docs/steps/S2.4.md).
# 8 000 estrazioni di Poisson con medie lam = 0.5 exp(N(0, 0.4^2)) (seme 99), 400 semi; per posizione z = (somma - 400 lam) / sqrt(400 lam).
# Uso: Rscript --vanilla tools/S2.4_C1_seeds.R

nb <- 8000; set.seed(99); lam <- 0.5 * exp(rnorm(nb, 0, 0.4))
run <- function(seeds, which = 2L, kind = "L'Ecuyer-CMRG") {
  acc <- numeric(nb); acc2 <- numeric(nb)
  for (s in seeds) { RNGkind(kind, "Inversion", "Rejection"); set.seed(s); x <- .Random.seed
    if (kind == "L'Ecuyer-CMRG") for (k in seq_len(which)) x <- parallel::nextRNGStream(x)
    assign(".Random.seed", x, envir = globalenv()); K <- rpois(nb, lam); acc <- acc + K; acc2 <- acc2 + K^2 }
  n <- length(seeds); m <- acc / n; v <- acc2 / n - m^2
  z <- (acc - n * lam) / sqrt(n * lam)
  c(maxz = max(abs(z)), var_over_mean = mean(v / lam))
}
print(run(501:900, 2L)); print(run(501:900, 1L)); print(run(sample(1e6, 400), 2L)); print(run(501:900, 0L, "Mersenne-Twister"))
