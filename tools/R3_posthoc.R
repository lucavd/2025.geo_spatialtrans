# tools/R3_posthoc.R — R3, analisi POST HOC (dopo aver visto i risultati; nessuna soglia, solo descrittivo).
# P1 CP-2 con StarDist in A4/A6 (sensibilita' D-4 pre-registrata, qui estesa a CP-2 e CP-3).
# P2 ipotesi del volume escluso: l'allineamento territorio–nucleo cresce con q = asse maggiore del nucleo /
#    distanza del generatore dal vicino piu' prossimo? (spiegherebbe allineamento in A1, A2, A4 e assenza in A5, A6).
# P3 effetto della selezione delle celle interne: eq_r mediano di tutti i territori (anche ritagliati) vs interne.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(ggplot2) })
source("tools/R3_voronoi_metrics.R")
OUT <- "/mnt/micron/geo_spatialtrans/R3"; RES <- "results/R3"; BASE_SEED <- 20261008
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]; rois$roi_idx <- seq_len(nrow(rois))
cells <- as.data.frame(read_parquet(file.path(OUT, "cells_all.parquet")))
nuc <- as.data.frame(read_parquet("results/R2/R2_nuclei_all.parquet", col_select = c("archetype", "roi_id", "method", "scale", "keep", "label", "major_um")))
nuc <- nuc[nuc$scale == 1 & nuc$keep, ]
gl <- do.call(rbind, lapply(split(unique(cells[, c("archetype", "roi_id", "method")]), seq_len(nrow(unique(cells[, c("archetype", "roi_id", "method")])))), function(k) {
  g <- as.data.frame(read_parquet(file.path(OUT, "real", sprintf("%s_%s_%s_gen.parquet", k$archetype, k$roi_id, k$method))))
  nn <- spatstat.geom::nndist(g$x, g$y)
  data.frame(archetype = k$archetype, roi_id = k$roi_id, method = k$method, idx = g$idx, label = g$label, d_nn = nn)
}))

cells <- merge(cells, gl, by = c("archetype", "roi_id", "method", "idx"))
cells <- merge(cells, nuc[, c("archetype", "roi_id", "method", "label", "major_um")], by = c("archetype", "roi_id", "method", "label"), all.x = TRUE)
ci <- cells[cells$interior, ]
# P1
p1 <- do.call(rbind, lapply(c("A4", "A6"), function(A) do.call(rbind, lapply(rois$roi_id[rois$archetype == A], function(roi) {
  do.call(rbind, lapply(c("spaceranger", "stardist_he"), function(me) {
    d <- ci[ci$archetype == A & ci$roi_id == roi & ci$method == me, ]
    e <- d[d$ecc_N >= 0.8 & d$ecc_T >= 0.5 & is.finite(d$theta_N), ]
    dth <- r3_axis_diff(e$theta_T, e$theta_N) * 180 / pi
    set.seed(BASE_SEED + 1000 * rois$roi_idx[rois$archetype == A & rois$roi_id == roi] + 400)
    perm <- replicate(2000, median(r3_axis_diff(e$theta_T, sample(e$theta_N)) * 180 / pi))
    data.frame(archetype = A, roi_id = roi, method = me, n_eligible = nrow(e), median_dtheta = median(dth), perm_p = mean(perm <= median(dth)),
               frac_cut = mean(d$frac_out > 0.05, na.rm = TRUE), nc_median = median(d$area_nuc / d$area), ecc_T_median = median(d$ecc_T))
  }))
}))))
write.csv(p1, file.path(RES, "R3_posthoc_P1_stardist.csv"), row.names = FALSE)
# P2
e <- ci[ci$ecc_N >= 0.8 & ci$ecc_T >= 0.5 & is.finite(ci$theta_N) & is.finite(ci$major_um), ]
e <- e[e$method == c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")[e$archetype] |
       (e$method == "stardist_he"), ]
e$dth <- r3_axis_diff(e$theta_T, e$theta_N) * 180 / pi; e$q <- e$major_um / e$d_nn
e$qbin <- cut(e$q, c(0, 0.25, 0.5, 0.75, 1, 1.5, Inf), right = FALSE)
p2 <- aggregate(dth ~ archetype + method + qbin, e, function(v) c(median = median(v), n = length(v)))
p2 <- do.call(data.frame, p2); names(p2) <- c("archetype", "method", "qbin", "median_dtheta", "n")
write.csv(p2, file.path(RES, "R3_posthoc_P2_alignment_vs_q.csv"), row.names = FALSE)
qa <- aggregate(q ~ archetype + method, e, median); names(qa)[3] <- "q_median"
write.csv(qa, file.path(RES, "R3_posthoc_P2_q_median.csv"), row.names = FALSE)
p2$lab <- paste(p2$archetype, p2$method)
g <- ggplot(p2[p2$n >= 30, ], aes(qbin, median_dtheta, group = lab, colour = archetype, linetype = method == "stardist_he")) +
  geom_line() + geom_point(aes(size = n)) + geom_hline(yintercept = 45, linetype = 2, colour = "grey50") +
  scale_size_area(max_size = 3, name = "n celle") + scale_linetype_discrete(name = "StarDist (sensibilità)") +
  labs(x = "q = asse maggiore del nucleo / distanza dal generatore più vicino", y = "mediana Δθ (°)",
       title = "Post hoc P2: allineamento territorio–nucleo in funzione dell'impaccamento (bin con n ≥ 30)") + theme_bw(base_size = 10)
ggsave(file.path(RES, "figures", "fig10_posthoc_alignment_q.png"), g, width = 9, height = 5.5, dpi = 150, bg = "white")
# P3
p3 <- do.call(rbind, lapply(split(cells, list(cells$archetype, cells$method), drop = TRUE), function(d)
  data.frame(archetype = d$archetype[1], method = d$method[1], eq_r_median_all = median(sqrt(d$area / pi)),
             eq_r_median_interior = median(sqrt(d$area[d$interior] / pi)), eq_r_of_mean_area = sqrt(mean(d$area) / pi),
             frac_interior = mean(d$interior), mean_area_interior_over_all = mean(d$area[d$interior]) / mean(d$area))))
write.csv(p3, file.path(RES, "R3_posthoc_P3_interior_bias.csv"), row.names = FALSE)
print(p1); print(qa); print(p3)
