
# S1.4d post hoc PH-1 (dichiarata dopo i risultati): da dove viene l'eccesso di cellule quasi senza citoplasma (P-3, tumore
# del colon)? Bersagli del modello di taglia prima del tetto contro il reale; effetto di phi piu' bassi sul troncamento (B-4).
L <- readLines("tools/S1.4d_bench.R"); k <- grep('^if \\(STAGE == "gen"\\)', L); args <- c("none", "/tmp/s14d_dg"); eval(parse(text = L[1:(k - 1)]))
rows <- list()
for (A in c("A1", "A2", "A3", "A4")) { method <- ARCH$primary[ARCH$archetype == A]
  for (kk in 1:5) { o14 <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, kk))); te <- load_real(A, o14$roi_test, method)
    cc <- te$cells; ev <- evset(cc); f2 <- o14$fits$M2; aref <- f2[["a_ref"]]
    for (s in 1:10) { set.seed(o14$seed + 100000 * s); z <- rnorm(te$n)
      tgt <- exp(log(f2[["rho"]]^2 * aref) + f2[["beta"]] * (log(cc$a_t) - log(aref)) + f2[["sigma"]] * z); r <- tgt / cc$a_t
      rows[[length(rows) + 1]] <- data.frame(archetype = A, roi = o14$roi_test, seed = s, sigma = f2[["sigma"]], beta = f2[["beta"]],
        tgt_nc95 = mean(r[ev] >= 0.95), tgt_nc100 = mean(r[ev] >= 1), real_nc95 = mean((cc$a_n / cc$a_t)[ev] >= 0.95), real_nc100 = mean((cc$a_n / cc$a_t)[ev] >= 1),
        tgt_nc_q99 = unname(quantile(r[ev], 0.99)), real_nc_q99 = unname(quantile((cc$a_n / cc$a_t)[ev], 0.99)),
        b4_phi100 = mean(tgt > 1.00 * cc$a_t), b4_phi095 = mean(tgt > 0.95 * cc$a_t), b4_phi090 = mean(tgt > 0.90 * cc$a_t), b4_phi085 = mean(tgt > 0.85 * cc$a_t)) } } }
x <- do.call(rbind, rows); dir.create("results/S1.4d", showWarnings = FALSE); write.csv(x, "results/S1.4d/S1.4d_posthoc_nc.csv", row.names = FALSE)
roi <- aggregate(. ~ archetype + roi, data = x[, setdiff(names(x), "seed")], FUN = mean)
ag <- aggregate(. ~ archetype, data = roi[, setdiff(names(roi), "roi")], FUN = mean)
se <- aggregate(. ~ archetype, data = roi[, setdiff(names(roi), "roi")], FUN = function(v) sd(v) / sqrt(length(v)))
print(ag, digits = 3); cat("SE fra ROI:\n"); print(se[, c("archetype", "tgt_nc95", "real_nc95", "b4_phi095", "b4_phi090")], digits = 3)
