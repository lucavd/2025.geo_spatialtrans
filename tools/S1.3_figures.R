# tools/S1.3_figures.R — figure del report S1.3 (sfondo bianco, un pannello per archetipo dove pertinente).
# Legge le uscite di tools/S1.3_roi.R, tools/S1.3_cp.R, tools/S1.3_analysis.R. Le figure non definiscono verdetti.
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.3_figures.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(ggplot2); library(arrow) })
source("R/04b1_extract_regions.R"); source("R/04b3_tessellate_voronoi.R"); source("tools/R3_voronoi_metrics.R")
OUT <- commandArgs(TRUE)[1]; if (is.na(OUT)) OUT <- "/mnt/micron/geo_spatialtrans/S1.3"
RES <- "results/S1.3"; FIG <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
R3 <- "/mnt/micron/geo_spatialtrans/R3"
PRIM <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
LAB <- c(A1 = "A1 epitelio (intestino)", A2 = "A2 tumore (CRC)", A3 = "A3 stroma (CRC)", A4 = "A4 linfonodo", A5 = "A5 corteccia", A6 = "A6 miocardio")
GEN <- c(real = "reale", CSR = "G1 CSR", RSA = "G2 RSA(d*)", RSArule = "G3 RSA regola")
COL <- c("reale" = "black", "G1 CSR" = "#999999", "G2 RSA(d*)" = "#0072B2", "G3 RSA regola" = "#D55E00")
th <- theme_bw(base_size = 10) + theme(strip.background = element_rect(fill = "white"), panel.grid.minor = element_blank(), legend.position = "bottom")
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
roi1 <- function(A) if (A == "A4") "f1" else "r1"
save <- function(p, f, w, h) ggsave(file.path(FIG, f), plot = p, width = w, height = h, dpi = 150, bg = "white")

# fig1: territori su un ritaglio per archetipo — reale, G2, G3 (replica 1), stessa maschera
crop_side <- c(A1 = 80, A2 = 100, A3 = 150, A4 = 60, A5 = 300, A6 = 300)
tiles <- do.call(rbind, lapply(names(PRIM), function(A) {
  roi <- roi1(A); w <- r3_roi_window(A, roi, rois); s <- crop_side[[A]]; cx <- w$side_um / 2; cy <- w$side_um / 2
  sq <- st_sfc(st_polygon(list(rbind(c(cx - s/2, cy - s/2), c(cx + s/2, cy - s/2), c(cx + s/2, cy + s/2), c(cx - s/2, cy + s/2), c(cx - s/2, cy - s/2)))))
  reg_g <- st_intersection(st_sfc(w$poly[[1]]), sq); if (!length(reg_g)) return(NULL)
  reg <- list(region_df = data.frame(region_id = 1L), region_polygons = st_sfc(st_union(reg_g)[[1]]))
  pts <- list(real = as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", A, roi, PRIM[[A]]))))[, c("x", "y")],
              RSA = readRDS(file.path(OUT, "null", sprintf("%s_%s_RSA_01.rds", A, roi)))$cells[, c("x", "y")],
              RSArule = readRDS(file.path(OUT, "null", sprintf("%s_%s_RSArule_01.rds", A, roi)))$cells[, c("x", "y")])
  do.call(rbind, lapply(names(pts), function(g) {
    p <- pts[[g]]; inr <- lengths(st_intersects(st_as_sf(p, coords = c("x", "y")), reg$region_polygons)) > 0; p <- p[inr, ]
    if (nrow(p) < 3) return(NULL)
    tv <- tessellate_voronoi(data.frame(cell_id = seq_len(nrow(p)), region_id = 1L, x = p$x, y = p$y), reg, verbose = FALSE)
    co <- st_coordinates(st_cast(tv$cell_territories, "MULTIPOLYGON"))       # territori misti POLYGON/MULTIPOLYGON (regione MULTIPOLYGON)
    data.frame(archetype = LAB[[A]], gen = GEN[[g]], id = paste(co[, "L3"], co[, "L2"], co[, "L1"]), x = co[, "X"] - cx, y = co[, "Y"] - cy)
  }))
}))
tiles$gen <- factor(tiles$gen, GEN[c("real", "RSA", "RSArule")]); tiles$archetype <- factor(tiles$archetype, LAB)
p1 <- ggplot(tiles, aes(x, y, group = id)) + geom_polygon(fill = "grey95", colour = "grey20", linewidth = 0.15) +
  facet_wrap(~ gen + archetype, nrow = 3, scales = "free", labeller = function(l) label_value(l, multi_line = FALSE)) +
  scale_y_reverse() + labs(x = "µm", y = "µm") + th + theme(aspect.ratio = 1)   # ritagli quadrati, scale proprie per pannello
save(p1, "fig1_tessellations.png", 13, 7)

# fig2: QQ delle aree normalizzate (interne, replica 1, pool dei 5 ROI) — reale (pacchetto) vs G1–G3
qq <- do.call(rbind, lapply(names(PRIM), function(A) {
  real <- do.call(c, lapply(list.files(file.path(OUT, "real"), sprintf("^%s_.*_%s\\.rds$", A, PRIM[[A]]), full.names = TRUE), function(f) {
    d <- readRDS(f)$geos$cells; d <- d[d$interior, ]; d$area / mean(d$area) }))
  pr <- seq(0.005, 0.995, by = 0.005); qr <- quantile(real, pr)
  do.call(rbind, lapply(c("CSR", "RSA", "RSArule"), function(g) {
    s <- do.call(c, lapply(list.files(file.path(OUT, "null"), sprintf("^%s_.*_%s_01\\.rds$", A, g), full.names = TRUE), function(f) {
      d <- readRDS(f)$cells; d <- d[d$interior, ]; d$area / mean(d$area) }))
    data.frame(archetype = LAB[[A]], gen = GEN[[g]], q_real = qr, q_sim = quantile(s, pr))
  }))
}))
qq$archetype <- factor(qq$archetype, LAB); qq$gen <- factor(qq$gen, GEN[-1])
p2 <- ggplot(qq, aes(q_real, q_sim, colour = gen)) + geom_abline(linetype = 2, colour = "grey50") + geom_line(linewidth = 0.7) +
  facet_wrap(~archetype, nrow = 2, scales = "free") + scale_colour_manual(values = COL, name = NULL) +
  labs(x = "quantile reale, area / media (interne)", y = "quantile simulato") + th
save(p2, "fig2_qq_area.png", 10, 6.5)

# fig3: check B per ROI (riferimento R3, pre-registrato) con bande di tolleranza
mb <- read.csv(file.path(RES, "S1.3_checkB_roi.csv"))
long <- do.call(rbind, lapply(c("B1", "B2", "B3", "B4"), function(b) data.frame(archetype = mb$archetype, roi = mb$roi_id, gen = GEN[mb$model], metric = b, value = mb[[b]])))
mlab <- c(B1 = "B-1 eq_r mediano (rel.)", B2 = "B-2 CV globale (rel.)", B3 = "B-3 CV locale (rel.)", B4 = "B-4 eccentricità (Δ ass.)")
tol <- data.frame(metric = c("B1", "B2", "B3", "B4"), t = c(0.05, 0.10, 0.10, 0.05)); tol$metric <- factor(mlab[tol$metric], mlab)
long$metric <- factor(mlab[long$metric], mlab); long$gen <- factor(long$gen, GEN[-1])
p3 <- ggplot(long, aes(archetype, value, colour = gen)) +
  geom_rect(data = tol, aes(xmin = -Inf, xmax = Inf, ymin = -t, ymax = t), inherit.aes = FALSE, fill = "#E8F3E8") +
  geom_hline(yintercept = 0, colour = "grey50") + geom_point(position = position_dodge(width = 0.6), size = 1.6) +
  facet_wrap(~metric, scales = "free_y", nrow = 2) + scale_colour_manual(values = COL, name = NULL) + labs(x = NULL, y = "simulato vs reale (per ROI)") + th
save(p3, "fig3_checkB.png", 10, 6.5)

# fig4: CP-1 — area scoperta per regione: testo del design vs D-S1.3.1
cp1 <- read.csv(file.path(RES, "S1.3_cp1_regions.csv")); cp1 <- cp1[cp1$n_cells > 0, ]
d4 <- rbind(data.frame(n = cp1$n_cells, f = pmax(cp1$uncovered_design, 1e-17), v = "testo del design §5.3 (globale)"),
            data.frame(n = cp1$n_cells, f = pmax(abs(cp1$uncovered_pkg), 1e-17), v = "D-S1.3.1 (per regione)"))
p4 <- ggplot(d4, aes(n, f, colour = v)) + geom_point(alpha = 0.5, size = 1) + geom_hline(yintercept = c(0.01, 1e-9), linetype = 2, colour = "grey40") +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = c("#D55E00", "#0072B2"), name = NULL) +
  labs(x = "cellule nella regione", y = "frazione di area scoperta") + th
save(p4, "fig4_cp1.png", 7, 4.5)

# fig5: CP-2 — eccentricità mediana vs anisotropia imposta k
cp2 <- read.csv(file.path(RES, "S1.3_cp2_anisotropy.csv")); b0 <- median(cp2$median_ecc_T[cp2$k == 1])
p5 <- ggplot(cp2, aes(factor(k), median_ecc_T)) + geom_jitter(width = 0.1, size = 1, colour = "grey40") +
  stat_summary(fun = median, geom = "point", colour = "#0072B2", size = 3) +
  geom_hline(yintercept = b0 + c(0, 0.05), linetype = c(1, 2), colour = c("grey50", "#D55E00")) +
  labs(x = "rapporto d'anisotropia imposto k", y = "eccentricità mediana dei territori (interne)") + th
save(p5, "fig5_cp2.png", 6, 4)

# fig6: prestazioni dei due motori (C-10c)
pf <- read.csv(file.path(RES, "S1.3_perf.csv")); pf <- pf[is.finite(pf$elapsed_s), ]
p6 <- ggplot(pf, aes(n, elapsed_s, colour = backend)) + geom_line() + geom_point() + geom_hline(yintercept = 150, linetype = 2) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = c(geos = "#0072B2", deldir = "#D55E00"), name = "motore") +
  labs(x = "punti (quadrato, densità A4)", y = "tempo di tessellate_voronoi() (s)") + th
save(p6, "fig6_engines_time.png", 6.5, 4.5)

# fig7: smussatura su un ritaglio A1 (cs = 0, 1/3, 1) + lacune per archetipo (CP-3)
w <- r3_roi_window("A1", "r1", rois); s <- 60; cx <- w$side_um / 2; cy <- w$side_um / 2
sq <- st_sfc(st_polygon(list(rbind(c(cx - s/2, cy - s/2), c(cx + s/2, cy - s/2), c(cx + s/2, cy + s/2), c(cx - s/2, cy + s/2), c(cx - s/2, cy - s/2)))))
reg <- list(region_df = data.frame(region_id = 1L), region_polygons = st_sfc(st_union(st_intersection(st_sfc(w$poly[[1]]), sq))[[1]]))
g <- as.data.frame(read_parquet(file.path(R3, "real", "A1_r1_cellpose_rgb_gen.parquet")))[, c("x", "y")]
g <- g[lengths(st_intersects(st_as_sf(g, coords = c("x", "y")), reg$region_polygons)) > 0, ]
sm <- do.call(rbind, lapply(c(0, 1/3, 1), function(cs) {
  tv <- tessellate_voronoi(data.frame(cell_id = seq_len(nrow(g)), region_id = 1L, x = g$x, y = g$y), reg, corner_smoothing = cs, verbose = FALSE)
  co <- st_coordinates(st_cast(tv$cell_territories, "MULTIPOLYGON")); data.frame(cs = sprintf("corner_smoothing = %.2f (n_iter %d)", cs, tv$info$n_iter), id = paste(co[, "L3"], co[, "L2"], co[, "L1"]), x = co[, "X"], y = co[, "Y"])
}))
p7 <- ggplot(sm, aes(x, y, group = id)) + geom_polygon(fill = "#F3E9D2", colour = "grey20", linewidth = 0.2) + facet_wrap(~cs) + coord_equal() + scale_y_reverse() + labs(x = "µm", y = "µm") + th
save(p7, "fig7_smoothing.png", 11, 4)

cat("figure scritte in", FIG, "\n")
