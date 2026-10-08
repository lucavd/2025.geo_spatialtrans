# tools/S1.3_predictions.R — previsioni del check B di S1.3 ricavate dai nulli di R3 (prima di qualsiasi codice S1.3).
# I nulli CSR e RSA(d*) di R3 (tools/R3_run.R, stadio null) usano seed_centroids() e la stessa maschera/regola
# d'interna che S1.3 usera': le loro distanze dal reale sono la previsione per G1 e G2.
# Uso (root del repo): Rscript --vanilla tools/S1.3_predictions.R  -> results/S1.3/S1.3_predictions_from_R3.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
dir.create("results/S1.3", showWarnings = FALSE, recursive = TRUE)
ns <- read.csv("results/R3/R3_null_summary.csv")
rs <- read.csv("results/R3/R3_roi_summary.csv")
rs <- rs[rs$method == ifelse(rs$archetype %in% c("A4", "A6"), "spaceranger", "cellpose_rgb"), ]
stopifnot(nrow(rs) == 30)
key <- c("median_eq_r", "cv", "cv_loc", "median_ecc_T")
agg <- aggregate(ns[, key], ns[, c("archetype", "roi_id", "model")], median)
m <- merge(agg, rs[, c("archetype", "roi_id", key)], by = c("archetype", "roi_id"), suffixes = c("_sim", "_real"))
stopifnot(nrow(m) == 60)
m$rel_eq_r <- m$median_eq_r_sim / m$median_eq_r_real - 1
m$rel_cv <- m$cv_sim / m$cv_real - 1
m$rel_cv_loc <- m$cv_loc_sim / m$cv_loc_real - 1
m$d_ecc <- m$median_ecc_T_sim - m$median_ecc_T_real
m$pass_B1 <- abs(m$rel_eq_r) <= 0.05; m$pass_B2 <- abs(m$rel_cv) <= 0.10
m$pass_B3 <- abs(m$rel_cv_loc) <= 0.10; m$pass_B4 <- abs(m$d_ecc) < 0.05
write.csv(m, "results/S1.3/S1.3_predictions_from_R3_roi.csv", row.names = FALSE)
f <- function(v) sprintf("%+.3f..%+.3f", min(v), max(v))
out <- do.call(rbind, lapply(split(m, list(m$archetype, m$model), drop = TRUE), function(d) data.frame(
  archetype = d$archetype[1], model = d$model[1],
  B1_rel_eq_r = f(d$rel_eq_r), B1_n_pass = sum(d$pass_B1), B2_rel_cv = f(d$rel_cv), B2_n_pass = sum(d$pass_B2),
  B3_rel_cv_loc = f(d$rel_cv_loc), B3_n_pass = sum(d$pass_B3), B4_d_ecc = f(d$d_ecc), B4_n_pass = sum(d$pass_B4))))
out <- out[order(out$model, out$archetype), ]
write.csv(out, "results/S1.3/S1.3_predictions_from_R3.csv", row.names = FALSE)
print(out, row.names = FALSE)
