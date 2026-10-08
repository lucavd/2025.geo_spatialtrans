# tools/S1.2_analysis.R — tabelle PASS/WARN/FAIL e figure di check B e controprove (S1.2)
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.2_analysis.R
# Input: /mnt/micron/geo_spatialtrans/S1.2/checkB/{real_pcf.rds, sim/*.rds}, results/S1.2/S1.2_B1b.csv, S1.2_manual_pcf.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(ggplot2); library(sf); library(spatstat.geom); library(png); library(arrow) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R")
OUT <- "/mnt/micron/geo_spatialtrans/S1.2/checkB"; RES <- "results/S1.2"; FIG <- file.path(RES, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ARCH <- read.csv(file.path(RES, "S1.2_archetype_params.csv"))
R_GRID <- seq(0, 30, by = 0.5); IR <- which(R_GRID > 0)            # r in (0, 30] per B2 (pre-registrato)
# D come pre-registrato: integrale dei trapezi su r = 0..30 (RA-checkB-07: la versione pre-revisione partiva da 0.5,
# deviazione non necessaria che cambiava l'esito di CP-1); D05 (r = 0.5..30) resta come analisi di sensibilita'.
I0 <- seq_along(R_GRID)
trap <- function(f, ii) sum((head(f, -1) + tail(f, -1)) / 2 * diff(R_GRID[ii]))
Dist   <- function(a, b) { stopifnot(all(is.finite(a)), all(is.finite(b))); trap((a[I0] - b[I0])^2, I0) }
Dist05 <- function(a, b) trap((a[IR] - b[IR])^2, IR)

real <- readRDS(file.path(OUT, "real_pcf.rds"))
rkey <- vapply(real, function(z) paste(z$archetype, z$roi_id), "")
names(real) <- rkey
sim <- do.call(rbind, lapply(list.files(file.path(OUT, "sim"), full.names = TRUE), readRDS))
gcols <- grep("^g_", names(sim), value = TRUE); stopifnot(length(gcols) == length(R_GRID))
sim$key <- paste(sim$archetype, sim$roi_id)
G <- as.matrix(sim[, gcols])
sim$D <- vapply(seq_len(nrow(sim)), function(i) if (sim$reps[i] > 0) Dist(G[i, ], real[[sim$key[i]]]$primary$g) else NA_real_, 0)
sim$D05 <- vapply(seq_len(nrow(sim)), function(i) if (sim$reps[i] > 0) Dist05(G[i, ], real[[sim$key[i]]]$primary$g) else NA_real_, 0)
write.csv(sim[, setdiff(names(sim), gcols)], file.path(RES, "S1.2_sim_models.csv"), row.names = FALSE)
cat(sprintf("[analysis] %d righe modello, %d ROI\n", nrow(sim), length(unique(sim$key))))

# curve lunghe per figure e revisione
curves <- rbind(
  do.call(rbind, lapply(real, function(z) rbind(
    data.frame(archetype = z$archetype, roi_id = z$roi_id, source = "reale primario", method = z$primary$method, r = R_GRID, g = z$primary$g),
    data.frame(archetype = z$archetype, roi_id = z$roi_id, source = "reale secondario", method = z$secondary$method, r = R_GRID, g = z$secondary$g)))),
  do.call(rbind, lapply(which(sim$model %in% c("CSR", "PD", "NUC")), function(i)
    data.frame(archetype = sim$archetype[i], roi_id = sim$roi_id[i], source = sim$model[i], method = sprintf("d=%.2f", sim$d[i]), r = R_GRID, g = G[i, ]))))
write.csv(curves, file.path(RES, "S1.2_pcf_curves.csv"), row.names = FALSE)

st3 <- function(x, pass, warn = -Inf) ifelse(x >= pass, "PASS", ifelse(x >= warn, "WARN", "FAIL"))
summ <- list()
add <- function(A, check, metric, value, threshold, status, predicted) summ[[length(summ) + 1]] <<-
  data.frame(archetype = A, check = check, metric = metric, value = value, threshold = threshold, status = status, predicted = predicted)

# ---- B1a -------------------------------------------------------------------------
pd <- sim[sim$model == "PD", ]
pd$area_mm2 <- vapply(pd$key, function(k) real[[k]]$area_mm2, 0)
pd <- merge(pd, ARCH[, c("archetype", "tot", "evid")], by = "archetype")
pd$in_int <- pd$dens_min >= pd$evid - 1 / pd$area_mm2 & pd$dens_max <= pd$tot + 1 / pd$area_mm2
write.csv(pd[, c("archetype", "roi_id", "area_mm2", "n_mean", "dens_mean", "dens_min", "dens_max", "evid", "tot", "in_int", "lost_mean", "feasible")],
          file.path(RES, "S1.2_B1a.csv"), row.names = FALSE)
for (A in ARCH$archetype) { s <- pd[pd$archetype == A, ]
  add(A, "B1a", "ROI con densita' PD nell'intervallo (20 repliche)", sprintf("%d/5", sum(s$in_int)), "5/5 (coerenza)", ifelse(all(s$in_int), "PASS", "FAIL"), "PASS") }

# ---- B1b -------------------------------------------------------------------------
b1 <- read.csv(file.path(RES, "S1.2_B1b.csv"))
b1$in_int <- b1$dens_seed42 >= b1$evid & b1$dens_seed42 <= b1$tot
write.csv(b1, file.path(RES, "S1.2_B1b.csv"), row.names = FALSE)
pred_b1b <- c(A1 = "PASS", A2 = "PASS", A3 = "PASS", A4 = "PASS (al limite)", A5 = "PASS", A6 = "FAIL")
for (A in ARCH$archetype) { s <- b1[b1$archetype == A, ]
  add(A, "B1b", "ROI graphclust con densita' nell'intervallo", sprintf("%d/5 (esclusa mediana %.3f)", sum(s$in_int), median(s$frac_excluded)), ">= 4/5",
      ifelse(sum(s$in_int) >= 4, "PASS", "FAIL"), pred_b1b[[A]]) }

# ---- B2: copertura dell'inviluppo reale ------------------------------------------
pred_b2 <- c(A1 = "FAIL", A2 = "FAIL", A3 = "FAIL", A4 = "WARN", A5 = "FAIL", A6 = "FAIL")
b2rows <- list()
for (A in ARCH$archetype) {
  ks <- rkey[startsWith(rkey, paste0(A, " "))]
  Rm <- sapply(ks, function(k) real[[k]]$primary$g)
  lo <- apply(Rm, 1, min); hi <- apply(Rm, 1, max)
  for (m in c("PD", "CSR", "NUC")) {
    gs <- colMeans(G[sim$archetype == A & sim$model == m, , drop = FALSE])
    cov <- mean(gs[IR] >= lo[IR] & gs[IR] <= hi[IR])
    b2rows[[length(b2rows) + 1]] <- data.frame(archetype = A, model = m, coverage = cov,
      cov_r_le10 = mean((gs >= lo & gs <= hi)[IR][R_GRID[IR] <= 10]), cov_r_gt10 = mean((gs >= lo & gs <= hi)[IR][R_GRID[IR] > 10]))
    if (m == "PD") add(A, "B2", "frazione di r in (0,30] dentro l'inviluppo (PD)", sprintf("%.2f", cov), ">= 0.80 PASS; 0.50-0.80 WARN", st3(cov, 0.8, 0.5), pred_b2[[A]])
  }
}
b2 <- do.call(rbind, b2rows); write.csv(b2, file.path(RES, "S1.2_B2.csv"), row.names = FALSE)

# ---- CP-1 --------------------------------------------------------------------------
cp1 <- merge(sim[sim$model == "PD", c("archetype", "roi_id", "D")], sim[sim$model == "CSR", c("archetype", "roi_id", "D")],
             by = c("archetype", "roi_id"), suffixes = c("_PD", "_CSR"))
cp1$reduction <- 1 - cp1$D_PD / cp1$D_CSR
write.csv(cp1, file.path(RES, "S1.2_CP1.csv"), row.names = FALSE)
pred_cp1 <- c(A1 = "sì", A2 = "sì", A3 = "sì", A4 = "sì", A5 = "sì", A6 = "no")
just <- c()
for (A in ARCH$archetype) { s <- cp1[cp1$archetype == A, ]
  j <- sum(s$D_PD < s$D_CSR) >= 4 && median(s$reduction) >= 0.25; just[A] <- j
  add(A, "CP-1", "ROI con D_PD < D_CSR; mediana riduzione", sprintf("%d/5; %.2f", sum(s$D_PD < s$D_CSR), median(s$reduction)),
      ">= 4/5 e >= 0.25", ifelse(j, "PASS", "FAIL"), pred_cp1[[A]]) }
add("tutti", "CP-1", "archetipi in cui PD e' giustificato", sprintf("%d/6", sum(just)), ">= 4/6", ifelse(sum(just) >= 4, "PASS", "FAIL"), "5/6 (no A6)")

# ---- CP-2 --------------------------------------------------------------------------
pred_cp2 <- c(A1 = "adeguata", A2 = "adeguata", A3 = "non adeguata", A4 = "non adeguata", A5 = "adeguata", A6 = "non adeguata")
gr <- sim[grepl("^G", sim$model), ]
cp2c <- do.call(rbind, lapply(split(gr, paste(gr$archetype, gr$d)), function(s)
  data.frame(archetype = s$archetype[1], d = s$d[1], n_roi_feasible = sum(s$feasible & s$reps == 20),
             D_med = if (all(s$feasible & s$reps == 20)) median(s$D) else NA, D_q25 = if (all(s$feasible)) quantile(s$D, .25) else NA,
             D_q75 = if (all(s$feasible)) quantile(s$D, .75) else NA)))
cp2c <- cp2c[order(cp2c$archetype, cp2c$d), ]
write.csv(cp2c, file.path(RES, "S1.2_CP2_curve.csv"), row.names = FALSE)
cp2 <- do.call(rbind, lapply(ARCH$archetype, function(A) {
  s <- cp2c[cp2c$archetype == A & !is.na(cp2c$D_med), ]; a <- ARCH[ARCH$archetype == A, ]
  dstar <- s$d[which.min(s$D_med)]; Dstar <- min(s$D_med)
  Drule <- median(sim$D[sim$archetype == A & sim$model == "PD"]); Dnuc <- median(sim$D[sim$archetype == A & sim$model == "NUC"])
  Dcsr <- median(sim$D[sim$archetype == A & sim$model == "CSR"])
  data.frame(archetype = A, d_star = dstar, D_star = Dstar, d_rule = a$d_rule, D_rule = Drule, ratio_rule = Drule / Dstar,
             d_nuc = a$d_nuc, D_nuc = Dnuc, ratio_nuc = Dnuc / Dstar, D_csr = Dcsr, d_max_feasible = max(s$d),
             frac_D_csr_left_at_dstar = Dstar / Dcsr)
}))
write.csv(cp2, file.path(RES, "S1.2_CP2.csv"), row.names = FALSE)
for (i in seq_len(nrow(cp2))) add(cp2$archetype[i], "CP-2", "D(d_regola)/D(d*); d* (µm)", sprintf("%.2f; d* %.1f (regola %.2f, nucleo %.1f: %.2f)", cp2$ratio_rule[i], cp2$d_star[i], cp2$d_rule[i], cp2$d_nuc[i], cp2$ratio_nuc[i]),
                                 "<= 1.25 adeguata", ifelse(cp2$ratio_rule[i] <= 1.25, "PASS", "FAIL"), pred_cp2[[cp2$archetype[i]]])

# ---- CP-3 --------------------------------------------------------------------------
cp3 <- do.call(rbind, lapply(real, function(z) data.frame(archetype = z$archetype, roi_id = z$roi_id,
  D_seg = Dist(z$primary$g, z$secondary$g), primary = z$primary$method, secondary = z$secondary$method,
  n_primary = z$primary$n, n_secondary = z$secondary$n)))
cp3 <- merge(cp3, sim[sim$model == "CSR", c("archetype", "roi_id", "D")], by = c("archetype", "roi_id"))
names(cp3)[names(cp3) == "D"] <- "D_CSR"
cp3 <- merge(cp3, sim[sim$model == "PD", c("archetype", "roi_id", "D")], by = c("archetype", "roi_id")); names(cp3)[names(cp3) == "D"] <- "D_PD"
write.csv(cp3, file.path(RES, "S1.2_CP3.csv"), row.names = FALSE)
pred_cp3 <- c(A1 = "incerto", A2 = "PASS", A3 = "PASS", A4 = "PASS", A5 = "PASS", A6 = "incerto")
p3 <- c()
for (A in ARCH$archetype) { s <- cp3[cp3$archetype == A, ]
  ok <- median(s$D_seg) < 0.5 * median(s$D_CSR); p3[A] <- ok
  add(A, "CP-3", "mediana D(primario, secondario) / mediana D_CSR; D_seg / D_PD", sprintf("%.2f; %.2f", median(s$D_seg) / median(s$D_CSR), median(s$D_seg) / median(s$D_PD)),
      "< 0.5", ifelse(ok, "PASS", "FAIL"), pred_cp3[[A]]) }
add("tutti", "CP-3", "archetipi con riferimento discriminante", sprintf("%d/6", sum(p3)), ">= 4/6 (previsione)", ifelse(sum(p3) >= 4, "PASS", "FAIL"), ">= 4/6")

S <- do.call(rbind, summ); write.csv(S, file.path(RES, "S1.2_checkB_summary.csv"), row.names = FALSE)
# sensibilita' (non pre-registrata): D integrato da r = 0.5
sens <- do.call(rbind, lapply(ARCH$archetype, function(A) {
  a <- sim[sim$archetype == A, ]; pdD <- a$D05[a$model == "PD"]; csD <- a$D05[a$model == "CSR"]
  z <- real[startsWith(rkey, paste0(A, " "))]; dseg <- vapply(z, function(q) Dist05(q$primary$g, q$secondary$g), 0)
  data.frame(archetype = A, variante = "D da r = 0.5", cp1_wins = sum(pdD < csD), cp1_red = median(1 - pdD / csD),
             cp1 = sum(pdD < csD) >= 4 && median(1 - pdD / csD) >= 0.25, cp3_ratio = median(dseg) / median(csD), cp3 = median(dseg) < 0.5 * median(csD))
}))
write.csv(sens, file.path(RES, "S1.2_CP_sensitivity_r05.csv"), row.names = FALSE)
print(S[, c("archetype", "check", "value", "status", "predicted")], right = FALSE)

# ---- figure ----------------------------------------------------------------------
thm <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
lab <- setNames(sprintf("%s (%s)", ARCH$archetype, ifelse(ARCH$primary == "spaceranger", "SR4", "Cellpose")), ARCH$archetype)
cr <- curves[curves$r > 0, ]
env <- do.call(rbind, lapply(split(cr[cr$source == "reale primario", ], paste(cr$archetype[cr$source == "reale primario"], cr$r[cr$source == "reale primario"])),
  function(s) data.frame(archetype = s$archetype[1], r = s$r[1], lo = min(s$g), hi = max(s$g), med = median(s$g))))
mm <- function(src) { s <- cr[cr$source == src, ]; aggregate(g ~ archetype + r, s, mean) }
sec <- aggregate(g ~ archetype + r, cr[cr$source == "reale secondario", ], median)
dstar_curves <- do.call(rbind, lapply(seq_len(nrow(cp2)), function(i) {
  ii <- which(sim$archetype == cp2$archetype[i] & sim$model == sprintf("G%04.1f", cp2$d_star[i]))
  data.frame(archetype = cp2$archetype[i], r = R_GRID[IR], g = colMeans(G[ii, IR, drop = FALSE]))
}))
lines_df <- rbind(cbind(mm("PD"), modello = "Poisson-disk (d = 2/3 eq_r)"), cbind(mm("CSR"), modello = "CSR"),
                  cbind(dstar_curves, modello = "hard-core d* (CP-2)"), cbind(sec, modello = "reale secondario (StarDist), mediana"))
p1 <- ggplot() + geom_ribbon(data = env, aes(r, ymin = lo, ymax = hi), fill = "grey75") +
  geom_line(data = env, aes(r, med), colour = "black", linewidth = 0.6) +
  geom_line(data = lines_df, aes(r, g, colour = modello, linetype = modello), linewidth = 0.6) +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
  facet_wrap(~archetype, ncol = 3, labeller = labeller(archetype = lab)) +
  scale_colour_manual(values = c("Poisson-disk (d = 2/3 eq_r)" = "#D55E00", "CSR" = "#0072B2", "hard-core d* (CP-2)" = "#009E73", "reale secondario (StarDist), mediana" = "grey30")) +
  scale_linetype_manual(values = c("Poisson-disk (d = 2/3 eq_r)" = "solid", "CSR" = "solid", "hard-core d* (CP-2)" = "dashed", "reale secondario (StarDist), mediana" = "dotted")) +
  coord_cartesian(ylim = c(0, 2.6)) + labs(x = "r (µm)", y = "g(r)", colour = NULL, linetype = NULL,
  title = "g(r): reale (grigio = min–max dei 5 ROI, nero = mediana) vs simulato (media di 5 ROI × 20 repliche)") + thm + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.2_fig1_pcf.png"), p1, width = 11, height = 7.5, dpi = 150, bg = "white")

cc <- cp2c[!is.na(cp2c$D_med), ]
vl <- rbind(data.frame(archetype = cp2$archetype, d = cp2$d_rule, tipo = "regola 2/3 eq_r"), data.frame(archetype = cp2$archetype, d = cp2$d_nuc, tipo = "diametro nucleare"),
            data.frame(archetype = cp2$archetype, d = cp2$d_star, tipo = "d*"))
p2 <- ggplot(cc, aes(d, D_med)) + geom_ribbon(aes(ymin = D_q25, ymax = D_q75), fill = "grey80") + geom_line() + geom_point(size = 0.8) +
  geom_vline(data = vl, aes(xintercept = d, colour = tipo, linetype = tipo)) + facet_wrap(~archetype, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("regola 2/3 eq_r" = "#D55E00", "diametro nucleare" = "#CC79A7", "d*" = "#009E73")) +
  labs(x = "distanza minima d (µm)", y = "D = ∫ (g_sim − g_reale)² dr  (mediana dei 5 ROI, banda IQR)", colour = NULL, linetype = NULL,
       title = "CP-2: discrepanza in funzione di d (solo d con tutte le repliche piazzate)") + thm + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.2_fig2_cp2.png"), p2, width = 11, height = 7, dpi = 150, bg = "white")

p3f <- ggplot(cp1, aes(D_CSR, D_PD)) + geom_abline(slope = 1, intercept = 0, colour = "grey50") + geom_abline(slope = 0.75, intercept = 0, colour = "grey50", linetype = "dashed") +
  geom_point(aes(colour = archetype), size = 2.5) + scale_x_log10() + scale_y_log10() +
  labs(x = "D_CSR", y = "D_PD", colour = NULL, title = "CP-1: Poisson-disk vs CSR per ROI (continua: uguaglianza; tratteggio: riduzione 25 %)") + thm
ggsave(file.path(FIG, "S1.2_fig3_cp1.png"), p3f, width = 7, height = 6, dpi = 150, bg = "white")

p4 <- ggplot(b1, aes(roi_id, dens_seed42)) + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = evid, ymax = tot), fill = "grey85") +
  geom_errorbar(aes(ymin = dens_min10, ymax = dens_max10), width = 0.2) + geom_point(aes(colour = in_int), size = 2.5) +
  geom_point(aes(y = expected), shape = 4) + facet_wrap(~archetype, scales = "free", ncol = 3) +
  scale_colour_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#D55E00"), labels = c(`TRUE` = "nell'intervallo", `FALSE` = "fuori")) +
  labs(x = "ROI", y = "cellule / mm² di tessuto", colour = NULL, title = "B1b: densita' simulata su mappe graphclust 8 µm (grigio = [evidenti, totale]; × = totale × (1 − esclusa))") + thm
ggsave(file.path(FIG, "S1.2_fig4_b1b.png"), p4, width = 11, height = 6.5, dpi = 150, bg = "white")

man <- read.csv(file.path(RES, "S1.2_manual_K.csv")); man <- man[!grepl("@", man$source) & man$r > 0, ]
srcl <- c(manual = "manuale (Luca)", cellpose_rgb = "Cellpose-SAM RGB", stardist_he = "StarDist HE", spaceranger = "Space Ranger 4")
man$fonte <- srcl[man$source]
p5 <- ggplot(man, aes(r, K_over_pir2, colour = fonte, linetype = fonte)) + geom_line(linewidth = 0.7) + geom_hline(yintercept = 1, colour = "grey50") + facet_wrap(~archetype, ncol = 3) +
  scale_colour_manual(values = c("manuale (Luca)" = "black", "Cellpose-SAM RGB" = "#E69F00", "StarDist HE" = "#56B4E9", "Space Ranger 4" = "#CC79A7")) +
  scale_linetype_manual(values = c("manuale (Luca)" = "solid", "Cellpose-SAM RGB" = "solid", "StarDist HE" = "dashed", "Space Ranger 4" = "solid")) +
  labs(x = "r (µm)", y = "K(r) / (π r²)  (pool delle 4 finestre R2b, senza lisciamento)", colour = NULL, linetype = NULL,
       title = "CP-3 descrittivo: punti manuali vs segmentatori nelle stesse finestre (1 = CSR; < 1 = deficit cumulato di coppie)") + thm + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.2_fig5_manual.png"), p5, width = 11, height = 6.5, dpi = 150, bg = "white")

# esempio visivo: ritaglio 150 × 150 µm del ROI 1 di ogni archetipo, reale primario vs PD (seed 1)
rois <- read.csv("results/R2/R2_rois_checked.csv"); nuc <- read_parquet("results/R2/R2_nuclei_all.parquet", col_select = c("archetype", "roi_id", "method", "scale", "keep", "x_um", "y_um"))
ex <- do.call(rbind, lapply(ARCH$archetype, function(A) {
  a <- ARCH[ARCH$archetype == A, ]; roi <- rois$roi_id[rois$archetype == A][1]; side_um <- rois$side_px[rois$archetype == A][1] * rois$um_per_px[rois$archetype == A][1]
  m <- readPNG(sprintf("results/R2/tissue_masks/%s_%s_valid_ds4.png", A, roi)); if (length(dim(m)) == 3) m <- m[, , 1]; m <- m > 0.5
  px <- side_um / ncol(m); idx <- which(m, arr.ind = TRUE)
  x0 <- side_um / 2 - 75; y0 <- side_um / 2 - 75
  sel <- idx[, 2] * px > x0 - 5 & (idx[, 2] - 1) * px < x0 + 155 & idx[, 1] * px > y0 - 5 & (idx[, 1] - 1) * px < y0 + 155
  reg <- extract_regions(data.frame(x = idx[sel, 2], y = idx[sel, 1], cl = 1L), pixel_size_um = px, min_region_area_um2 = 0, cluster_col = "cl", simplify_tol_um = 0, stride = 1L, verbose = FALSE)
  o <- seed_centroids(reg, data.frame(cell_type = "c", density = a$tot), data.frame(cluster_id = 1, cell_type = "c", fraction = 1), random_seed = 1, verbose = FALSE)
  d <- nuc[nuc$archetype == A & nuc$roi_id == roi & nuc$method == a$primary & nuc$scale == 1 & nuc$keep, ]
  inbox <- function(x, y) x >= x0 & x <= x0 + 150 & y >= y0 & y <= y0 + 150
  rbind(data.frame(archetype = A, fonte = sprintf("reale (%s)", ifelse(a$primary == "spaceranger", "SR4", "Cellpose")), x = d$x_um, y = d$y_um)[inbox(d$x_um, d$y_um), ],
        data.frame(archetype = A, fonte = "simulato PD", x = o$centroids$x, y = o$centroids$y)[inbox(o$centroids$x, o$centroids$y), ])
}))
ex$fonte2 <- ifelse(grepl("^reale", ex$fonte), "reale", "simulato PD (densita' totale manuale)")
p6 <- ggplot(ex, aes(x, y)) + geom_point(size = 0.5) + facet_grid(fonte2 ~ archetype) + coord_equal() + scale_y_reverse() +
  labs(x = "x (µm)", y = "y (µm)", title = "Ritaglio 150 × 150 µm al centro del ROI 1: centroidi reali vs simulati") + thm + theme(axis.text = element_text(size = 7))
ggsave(file.path(FIG, "S1.2_fig6_example.png"), p6, width = 13, height = 5.2, dpi = 150, bg = "white")
cat("[analysis] fatto\n")
