# tools/R3_figures.R — figure del report R3 (sfondo bianco, un pannello per archetipo). Legge le uscite di
# tools/R3_run.R e tools/R3_analysis.R. Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/R3_figures.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(ggplot2); library(sf) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R")
OUT <- "/mnt/micron/geo_spatialtrans/R3"; RES <- "results/R3"; FIG <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE)
PRIM <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
LAB <- c(A1 = "A1 epitelio (intestino)", A2 = "A2 tumore (CRC)", A3 = "A3 stroma (CRC)", A4 = "A4 linfonodo",
         A5 = "A5 corteccia", A6 = "A6 miocardio")
th <- theme_bw(base_size = 10) + theme(strip.background = element_rect(fill = "white"), panel.grid.minor = element_blank())
cells <- as.data.frame(read_parquet(file.path(OUT, "cells_all.parquet")))
cells <- cells[cells$method == PRIM[cells$archetype], ]
ci <- cells[cells$interior, ]; ci$lab <- factor(LAB[ci$archetype], LAB)
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px

# fig1: territori su un ritaglio 120x120 µm (ROI r1/f1 di ogni archetipo), interne vs ritagliate
crop_side <- c(A1 = 80, A2 = 100, A3 = 150, A4 = 60, A5 = 250, A6 = 250)
tiles_df <- do.call(rbind, lapply(names(PRIM), function(A) {
  roi <- if (A == "A4") "f1" else "r1"
  d <- cells[cells$archetype == A & cells$roi_id == roi, ]
  s <- crop_side[[A]]; cx <- median(d$x); cy <- median(d$y)
  sel <- abs(d$x - cx) < s / 2 + 60 & abs(d$y - cy) < s / 2 + 60
  dd <- deldir::deldir(d$x[sel], d$y[sel], rw = c(cx - s / 2 - 60, cx + s / 2 + 60, cy - s / 2 - 60, cy + s / 2 + 60), round = FALSE)
  tl <- deldir::tile.list(dd); ids <- which(sel)
  do.call(rbind, lapply(seq_along(tl), function(k) {
    t <- tl[[k]]; j <- ids[t$ptNum]
    if (abs(d$x[j] - cx) > s / 2 || abs(d$y[j] - cy) > s / 2) return(NULL)
    data.frame(archetype = A, id = paste(A, j), x = (t$x - cx), y = (t$y - cy), interior = d$interior[j])
  }))
}))
tiles_df$lab <- factor(LAB[tiles_df$archetype], LAB)
gen_df <- do.call(rbind, lapply(names(PRIM), function(A) { roi <- if (A == "A4") "f1" else "r1"
  d <- cells[cells$archetype == A & cells$roi_id == roi, ]; s <- crop_side[[A]]; cx <- median(d$x); cy <- median(d$y)
  k <- abs(d$x - cx) <= s / 2 & abs(d$y - cy) <= s / 2; data.frame(archetype = A, x = d$x[k] - cx, y = d$y[k] - cy) }))
gen_df$lab <- factor(LAB[gen_df$archetype], LAB)
p <- ggplot(tiles_df, aes(x, y, group = id)) + geom_polygon(aes(fill = interior), colour = "grey30", linewidth = 0.2) +
  geom_point(data = gen_df, aes(x, y), inherit.aes = FALSE, size = 0.4) +
  scale_fill_manual(values = c(`TRUE` = "#cfe3f3", `FALSE` = "#f3d9cf"), labels = c(`TRUE` = "interna", `FALSE` = "ritagliata"), name = NULL) +
  scale_y_reverse() + facet_wrap(~lab, scales = "free", nrow = 2) + theme(aspect.ratio = 1) +
  labs(x = "µm (dal centro del ritaglio)", y = "µm", title = "Voronoi sui centroidi nucleari reali (ROI r1 / f1; lato del ritaglio adattato alla densità)") + th
ggsave(file.path(FIG, "fig1_tiles.png"), p, width = 11, height = 7.5, dpi = 150, bg = "white")

# fig2: area normalizzata (A / media del ROI): reale vs CSR vs RSA (replica 1)
nf <- list.files(file.path(OUT, "null"), "_01\\.rds$", full.names = TRUE)
nd <- do.call(rbind, lapply(nf, function(f) { o <- readRDS(f); z <- o$cells[o$cells$interior, ]
  data.frame(archetype = o$archetype, roi_id = o$roi_id, source = o$model, a = z$area / mean(z$area), ecc_T = z$ecc_T) }))
rd <- do.call(rbind, lapply(split(ci, list(ci$archetype, ci$roi_id), drop = TRUE), function(z)
  data.frame(archetype = z$archetype[1], roi_id = z$roi_id[1], source = "reale", a = z$area / mean(z$area), ecc_T = z$ecc_T)))
ad <- rbind(rd, nd); ad$lab <- factor(LAB[ad$archetype], LAB); ad$source <- factor(ad$source, c("reale", "CSR", "RSA"))
cols <- c(reale = "black", CSR = "#d95f02", RSA = "#1b9e77")
p <- ggplot(ad, aes(a, colour = source)) + geom_density(linewidth = 0.6, adjust = 0.8) +
  stat_function(fun = function(y) dgamma(y, 3.5, 3.5), colour = "grey50", linetype = 2, inherit.aes = FALSE) +
  scale_colour_manual(values = cols, name = NULL) + coord_cartesian(xlim = c(0, 3)) + facet_wrap(~lab, nrow = 2) +
  labs(x = "area del territorio / media del ROI (celle interne, pool 5 ROI)", y = "densità",
       title = "Distribuzione delle aree normalizzate; tratteggio: gamma(7/2, 7/2) di Poisson-Voronoi (Ferenc & Néda 2007)") + th
ggsave(file.path(FIG, "fig2_area_norm.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")

# fig3: CP-1, CV locale per ROI: reale vs intervallo dei nulli
cp1 <- read.csv(file.path(RES, "R3_cp1_roi.csv")); ns <- read.csv(file.path(RES, "R3_null_summary.csv"))
ns$lab <- factor(LAB[ns$archetype], LAB); cp1$lab <- factor(LAB[cp1$archetype], LAB)
p <- ggplot(ns, aes(roi_id, cv_loc, colour = model)) + geom_point(position = position_dodge(0.5), size = 0.8, alpha = 0.6) +
  geom_point(data = cp1, aes(roi_id, cv_loc_real), inherit.aes = FALSE, shape = 4, size = 2.5, stroke = 1.1) +
  scale_colour_manual(values = cols[c("CSR", "RSA")], name = "nullo (20 repl.)") + facet_wrap(~lab, scales = "free", nrow = 2) +
  labs(x = "ROI", y = "CV locale dell'area (A · λ̂ locale)", title = "CP-1: CV locale reale (×) contro CSR e RSA(d*) nello stesso ROI, stesso n") + th
ggsave(file.path(FIG, "fig3_cp1_cvloc.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")

# fig4: eq_radius con intervallo C8 provvisorio e intervallo dalle densità manuali
man <- data.frame(archetype = names(PRIM), lo = sqrt(1e6 / (pi * c(12419, 8924, 3100, 28096, 1185, 988))),
                  hi = sqrt(1e6 / (pi * c(10291, 7376, 2311, 26844, 953, 906)))); man$lab <- factor(LAB[man$archetype], LAB)
p <- ggplot(ci, aes(eq_r)) + annotate("rect", xmin = 5, xmax = 25, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_histogram(bins = 60, fill = "grey40") + geom_rect(data = man, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
  fill = "#e7298a", alpha = 0.35) + facet_wrap(~lab, scales = "free", nrow = 2) +
  labs(x = "eq_radius del territorio (µm)", y = "celle interne",
       title = "eq_radius: grigio = intervallo C8 provvisorio [5, 25] µm; rosa = eq_r della media dalle densità manuali [totale, evidenti]") + th
ggsave(file.path(FIG, "fig4_eq_radius.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")

# fig5: rapporto dei raggi con il catalogo del design
des <- data.frame(archetype = c("A1", "A2", "A3", "A4"), r = c(0.45, 0.45, 0.45, 0.55)); des$lab <- factor(LAB[des$archetype], LAB)
p <- ggplot(ci[is.finite(ci$ratio), ], aes(ratio)) + geom_histogram(bins = 60, fill = "grey40") +
  geom_vline(data = des, aes(xintercept = r), colour = "#e7298a") + geom_vline(data = des, aes(xintercept = r * 0.8), colour = "#e7298a", linetype = 3) +
  geom_vline(data = des, aes(xintercept = r * 1.2), colour = "#e7298a", linetype = 3) + coord_cartesian(xlim = c(0, 1.2)) +
  facet_wrap(~lab, scales = "free_y", nrow = 2) +
  labs(x = "rapporto dei raggi √(area nucleo / area territorio)", y = "celle interne",
       title = "B-R3.2: rapporto dei raggi reale; rosa = nucleus_to_eq_ratio del catalogo (design §4.4) ± 20 %") + th
ggsave(file.path(FIG, "fig5_ratio.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")

# fig6: CP-2, Δθ fra asse del territorio e asse del nucleo
d <- ci[ci$ecc_N >= 0.8 & ci$ecc_T >= 0.5 & is.finite(ci$theta_N), ]
d$dth <- r3_axis_diff(d$theta_T, d$theta_N) * 180 / pi
p <- ggplot(d, aes(dth)) + geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 90, 5), fill = "grey40") +
  geom_hline(yintercept = 1 / 90, colour = "#d95f02", linetype = 2) + facet_wrap(~lab, nrow = 2) +
  labs(x = "Δθ asse territorio − asse nucleo (°)", y = "densità",
       title = "CP-2: allineamento territorio–nucleo (nuclei ecc ≥ 0.8, territori ecc ≥ 0.5); tratteggio = nullo uniforme") + th
ggsave(file.path(FIG, "fig6_cp2_dtheta.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")

# fig7: CP-2a/c, eccentricita' del territorio reale vs nulli
p <- ggplot(ad, aes(ecc_T, colour = source)) + geom_density(linewidth = 0.6) + scale_colour_manual(values = cols, name = NULL) +
  facet_wrap(~lab, nrow = 2) + labs(x = "eccentricità del territorio e_T", y = "densità",
  title = "Eccentricità dei territori: reale vs CSR vs RSA(d*) (celle interne, replica 1 per i nulli)") + th
ggsave(file.path(FIG, "fig7_eccT.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")

# fig8: CP-3, frazione del nucleo fuori dal proprio territorio
p <- ggplot(ci[is.finite(ci$frac_out), ], aes(pmax(frac_out, 1e-4))) + geom_histogram(bins = 50, fill = "grey40") +
  geom_vline(xintercept = 0.05, colour = "#e7298a") + scale_x_log10() + facet_wrap(~lab, scales = "free_y", nrow = 2) +
  labs(x = "frazione dei pixel del nucleo fuori dal proprio territorio (0 riportato a 1e-4)", y = "nuclei (interne)",
       title = "CP-3: il Voronoi taglia i nuclei reali? rosa = soglia del 5 %") + th
ggsave(file.path(FIG, "fig8_cp3_containment.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")

# fig9: CP-4, finestre R2b — manuale vs segmentatore
cc <- do.call(rbind, lapply(list.files(file.path(OUT, "calib"), "\\.parquet$", full.names = TRUE), function(f) as.data.frame(read_parquet(f))))
cc <- cc[cc$interior, ]; cc$archetype <- substr(cc$win_id, 1, 2); cc$lab <- factor(LAB[cc$archetype], LAB)
cc$fonte <- ifelse(cc$source == "manual_Luca", "manuale (Luca)", "segmentatore primario")
p <- ggplot(cc, aes(sqrt(area / pi), colour = fonte)) + stat_ecdf(linewidth = 0.6) + facet_wrap(~lab, scales = "free_x", nrow = 2) +
  scale_colour_manual(values = c(`manuale (Luca)` = "black", `segmentatore primario` = "#7570b3"), name = NULL) +
  labs(x = "eq_radius del territorio (µm)", y = "ECDF", title = "CP-4: Voronoi sui punti manuali vs sul segmentatore, stesse 4 finestre R2b per archetipo") + th
ggsave(file.path(FIG, "fig9_cp4_calibration.png"), p, width = 11, height = 6.5, dpi = 150, bg = "white")
cat("figure scritte:", length(list.files(FIG)), "\n")
