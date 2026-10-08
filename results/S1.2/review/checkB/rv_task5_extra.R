# rv_task5_extra.R — revisione avversariale S1.2 / checkB: B2 per variante di banda e relazione d* ~ nucleo / regola.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task5_extra.R
SO <- "results/S1.2/review/checkB/sens"; RV <- "results/S1.2/review/checkB"; R_GRID <- seq(0, 30, by = 0.5); I <- R_GRID > 0
S <- do.call(rbind, lapply(list.files(SO, pattern = "\\.rds$", full.names = TRUE), readRDS))
out <- list()
for (v in c("x1", "x0.5", "x2", "own")) for (A in paste0("A", 1:6)) {
  s <- S[S$archetype == A & S$variant == v, ]; rois <- sort(unique(s$roi_id))
  Rm <- sapply(rois, function(r) s$g[s$roi_id == r & s$source == "real_primary"]); lo <- apply(Rm, 1, min); hi <- apply(Rm, 1, max)
  for (m in c("PD", "CSR", "NUC")) { Gm <- sapply(rois, function(r) s$g[s$roi_id == r & s$source == m]); gs <- rowMeans(Gm)
    per <- apply(Gm, 2, function(g) mean((g >= lo & g <= hi)[I]))
    out[[length(out) + 1]] <- data.frame(variant = v, archetype = A, model = m, coverage = mean((gs >= lo & gs <= hi)[I]), cov_per_roi_max = max(per)) } }
b2 <- do.call(rbind, out); write.csv(b2, file.path(RV, "rv_task5_b2_variants.csv"), row.names = FALSE)
print(reshape(b2[b2$model == "PD", c("variant", "archetype", "coverage")], idvar = "archetype", timevar = "variant", direction = "wide"), digits = 2)
print(reshape(b2[b2$model == "CSR", c("variant", "archetype", "coverage")], idvar = "archetype", timevar = "variant", direction = "wide"), digits = 2)
ss <- read.csv(file.path(RV, "rv_task4_sens_summary.csv")); ss <- ss[ss$archetype != "tutti", ]
nuc <- c(A1 = 4.1, A2 = 6.3, A3 = 6.3, A4 = 4.2, A5 = 11.3, A6 = 4.8)
cr <- do.call(rbind, lapply(split(ss, ss$variant), function(x) data.frame(variant = x$variant[1],
  spearman_dstar_dnuc = cor(x$d_star, nuc[x$archetype], method = "spearman"), spearman_dstar_drule = cor(x$d_star, x$d_rule, method = "spearman"),
  pearson_dstar_dnuc = cor(x$d_star, nuc[x$archetype]), pearson_dstar_drule = cor(x$d_star, x$d_rule),
  ratio_dstar_dnuc = paste(sprintf("%.2f", x$d_star / nuc[x$archetype]), collapse = "/"), ratio_dstar_drule = paste(sprintf("%.2f", x$d_star / x$d_rule), collapse = "/"))))
write.csv(cr, file.path(RV, "rv_task5_dstar_corr.csv"), row.names = FALSE); print(cr, digits = 3)
