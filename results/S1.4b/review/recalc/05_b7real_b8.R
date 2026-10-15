.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

suppressPackageStartupMessages(library(spatstat.geom))
OUTD <- "results/S1.4b/review/recalc"
ph1 <- read.csv("results/S1.4/S1.4_ph1.csv", stringsAsFactors = FALSE)
fs <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.4b/loo", full.names = TRUE))
rows <- list()
for (f in fs) {
  L <- readRDS(f); A <- L$archetype; roi <- L$roi_test; me <- L$method
  te <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.4/real/%s_%s_%s.rds", A, roi, me)); cc <- te$cells
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  b7 <- mean(sqrt(cc$a_n[ev] / pi) > cc$d_gen[ev])
  # B-8: raggi da a_n se appaiato, altrimenti a_own
  a <- ifelse(cc$matched & is.finite(cc$a_n) & cc$a_n > 0, cc$a_n, cc$a_own); r <- sqrt(a / pi)
  kp <- sqrt(2 * pi / (32 * sin(2 * pi / 32)))
  P <- ppp(cc$x, cc$y, window = owin(range(cc$x) + c(-1, 1), range(cc$y) + c(-1, 1)), check = FALSE)
  cp <- closepairs(P, 2 * kp * max(r, na.rm = TRUE) / 0.99, what = "indices", twice = TRUE)
  d <- sqrt((cc$x[cp$i] - cc$x[cp$j])^2 + (cc$y[cp$i] - cc$y[cp$j])^2)
  hit_c <- tabulate(cp$i[r[cp$i] + r[cp$j] > d], te$n) > 0
  hit_r <- tabulate(cp$i[kp * (r[cp$i] + r[cp$j]) > 0.99 * d], te$n) > 0
  evj <- ev[cp$j]
  hit_c_ev <- tabulate(cp$i[r[cp$i] + r[cp$j] > d & evj], te$n) > 0
  p <- ph1[ph1$archetype == A & ph1$roi == roi, ]
  rows[[f]] <- data.frame(archetype = A, roi = roi, method = me, n_ev = sum(ev), ph1_n = if (nrow(p)) p$n else NA,
    b7_mine = b7, b7_loo = L$b7_real, b7_ph1 = if (nrow(p)) p$circ_gen_1 else NA,
    b8c_mine = mean(hit_c[ev]), b8c_loo = L$b8[["circle"]], b8r_mine = mean(hit_r[ev]), b8r_loo = L$b8[["rule"]],
    b8c_mine_evonly = mean(hit_c_ev[ev]), na_r = sum(!is.finite(r)))
}
O <- do.call(rbind, rows); write.csv(O, file.path(OUTD, "recalc_b7real_b8.csv"), row.names = FALSE)
options(width = 220); print(O, row.names = FALSE, digits = 5)
cat("max|b7 mine-loo|", max(abs(O$b7_mine - O$b7_loo)), " max|b7 mine-ph1|", max(abs(O$b7_mine - O$b7_ph1), na.rm = TRUE),
    " n mismatch", sum(O$n_ev != O$ph1_n, na.rm = TRUE), "\n")
cat("max|b8c mine-loo|", max(abs(O$b8c_mine - O$b8c_loo)), " max|b8r mine-loo|", max(abs(O$b8r_mine - O$b8r_loo)), "\n")
