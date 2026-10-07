# S1.1 — figure del report (dopo R/testing/test_S1.1.R)
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.1_figures.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(ggplot2) })
source("R/04b1_extract_regions.R")
IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"; OUTD <- "/mnt/micron/geo_spatialtrans/S1.1/outputs"
RES <- "results/S1.1"; FIG <- file.path(RES, "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
load_in  <- function(id) readRDS(file.path(IN, paste0(id, ".rds")))
load_out <- function(id) readRDS(file.path(OUTD, paste0(id, ".rds")))
pal <- function(n) grDevices::hcl.colors(max(n, 2), "Dark 3")

# fondo: tutte le celle di tessuto in grigio (cio' che resta grigio = escluso dal filtro)
draw_map <- function(id, main, lwd = 0.3, border = "grey20") {
  o <- load_in(id); z <- load_out(id); px <- o$pixel_size_um
  x <- o$clust$x; y <- o$clust$y
  s <- z$rf$info$stride
  plot(NA, xlim = c(min(x) - 1, max(x) - 1 + s) * px, ylim = rev(c(min(y) - 1, max(y) - 1 + s) * px),
       asp = 1, axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.85)
  rect((x - 1) * px, (y - 1) * px, (x - 1 + s) * px, (y - 1 + s) * px, col = "grey75", border = NA)
  if (z$rf$info$n_regions > 0) {
    cl <- z$rf$region_df$cluster_id; lv <- sort(unique(o$clust[[o$cluster_col]]))
    plot(z$rf$region_polygons, col = pal(length(lv))[match(cl, lv)], border = border, lwd = lwd, add = TRUE)
  }
  mtext(sprintf("%d regioni · %.1f%% area esclusa", z$rf$info$n_regions, 100 * z$rf$info$frac_area_excluded),
        side = 1, line = 0.2, cex = 0.6)
}
png_open <- function(f, w, h) { png(file.path(FIG, f), width = w, height = h, res = 150, bg = "white") }

# Fig 1 — avversari (poligoni esatti sopra la griglia dei pixel)
adv <- c("adv_checker", "adv_donut", "adv_donut_island", "adv_corner", "adv_pinch_hole", "adv_line",
         "adv_border", "adv_blobs_9_10", "adv_single_cluster", "adv_pixel35", "adv_L")
png_open("S1.1_fig1_adversarial.png", 2400, 1800); par(mfrow = c(3, 4), mar = c(1.5, 0.5, 2, 0.5))
for (id in adv) {
  o <- load_in(id); r <- extract_regions(o$clust, 1, 0, simplify_tol_um = 0, verbose = FALSE)
  x <- o$clust$x; y <- o$clust$y
  plot(NA, xlim = c(min(x) - 2, max(x) + 1), ylim = rev(c(min(y) - 2, max(y) + 1)), asp = 1, axes = FALSE,
       xlab = "", ylab = "", main = sprintf("%s\n%d regioni, area %s", id, nrow(r$region_df),
                                            paste(head(r$region_df$area_um2, 3), collapse = "/")), cex.main = 0.8)
  if (max(x) <= 60 && max(y) <= 60) abline(v = seq(min(x) - 2, max(x) + 1), h = seq(min(y) - 2, max(y) + 1), col = "grey92", lwd = 0.4)
  plot(r$region_polygons, col = adjustcolor(pal(4)[match(r$region_df$cluster_id, sort(unique(r$region_df$cluster_id)))], 0.7),
       border = "black", lwd = 0.8, add = TRUE)
}
plot.new(); text(0.5, 0.5, "grigio chiaro = griglia dei pixel\nbordo nero = poligono esatto\n(lati dei pixel, buchi inclusi)", cex = 0.9)
dev.off()

# Fig 2 — sintetici I1-I6
png_open("S1.1_fig2_synthetic.png", 2400, 1650); par(mfrow = c(2, 3), mar = c(1.5, 0.5, 2, 0.5))
for (id in c("I1_syn600_c1", "I2_syn600_c2", "I3_syn600_c3", "I4_syn600_c4", "I5_syn600_c2_labels", "I6_syn6800x6500_c2"))
  draw_map(id, id, lwd = 0.2)
dev.off()

# Fig 3 — mappe reali intere (graphclust 8 µm)
png_open("S1.1_fig3_real_full.png", 3000, 1800); par(mfrow = c(2, 3), mar = c(1.5, 0.5, 2, 0.5))
for (a in c("A1", "A2A3", "A4", "A5", "A6")) draw_map(paste0("real_full_", a), paste0("real_full_", a, " — ", load_in(paste0("real_full_", a))$meta$dataset), lwd = 0.05, border = NA)
plot.new(); text(0.5, 0.5, "colore = cluster graphclust\n(regione >= 100 µm²)\ngrigio = componenti escluse\n(1 bin = 64 µm²)", cex = 1)
dev.off()

# Fig 4 — ROI per archetipo (6 righe x 5 ROI)
rois <- read.csv("results/R2/R2_rois_checked.csv", stringsAsFactors = FALSE)
png_open("S1.1_fig4_real_roi_by_archetype.png", 2500, 3000); par(mfrow = c(6, 5), mar = c(1.2, 0.3, 1.6, 0.3))
for (k in seq_len(nrow(rois))) draw_map(sprintf("roi_%s_%s", rois$archetype[k], rois$roi_id[k]),
                                        sprintf("%s %s", rois$archetype[k], rois$roi_id[k]), lwd = 0.2)
dev.off()

# Fig 5 — descrittivo B
desc <- read.csv(file.path(RES, "S1.1_checkB_descriptive.csv"), stringsAsFactors = FALSE)
all <- readRDS(file.path(OUTD, "S1.1_regions_all.rds"))
reg <- all$regions
reg$panel <- ifelse(reg$group == "real_roi", sub("^roi_(A[1-6])_.*", "\\1", reg$id),
              ifelse(reg$group == "syn", sub("_.*", "", reg$id), NA))
reg <- reg[!is.na(reg$panel), ]
reg$kind <- ifelse(grepl("^A", reg$panel), "reale (ROI R2, 8 µm)", "sintetico (k-means, 1 µm)")
g1 <- ggplot(reg, aes(area_px_um2, colour = panel)) + stat_ecdf(aes(weight = area_px_um2), geom = "step") +
  scale_x_log10() + facet_wrap(~kind) + theme_bw(base_size = 10) +
  labs(x = "area della regione (µm², scala log)", y = "frazione cumulativa dell'area di tessuto",
       colour = NULL, title = "Distribuzione delle aree delle regioni pesata per area (regioni >= 100 µm²)")
ggsave(file.path(FIG, "S1.1_fig5a_checkB_ecdf.png"), g1, width = 9, height = 4, dpi = 150, bg = "white")
dd <- desc[desc$group %in% c("real_roi", "syn"), ]
dd$panel <- ifelse(dd$group == "real_roi", sub("^roi_(A[1-6])_.*", "\\1", dd$id), sub("_.*", "", dd$id))
dl <- rbind(data.frame(panel = dd$panel, metric = "frazione di area esclusa (< 100 µm²)", value = dd$frac_area_excluded),
            data.frame(panel = dd$panel, metric = "frazione di componenti di 1 pixel/bin", value = dd$frac_components_single_px),
            data.frame(panel = dd$panel, metric = "log10 mediana pesata dell'area (µm²)", value = log10(dd$area_weighted_median_um2)))
g2 <- ggplot(dl, aes(panel, value)) + geom_point(alpha = 0.7) + stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "firebrick") +
  facet_wrap(~metric, scales = "free_y") + theme_bw(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Frammentazione per archetipo (5 ROI, punti) e sintetici I1-I6; barra = mediana")
ggsave(file.path(FIG, "S1.1_fig5b_checkB_fragmentation.png"), g2, width = 10, height = 3.8, dpi = 150, bg = "white")

# Fig 6 — CP-1 Moore
mo <- read.csv(file.path(RES, "S1.1_cp1_moore.csv"))
mo$relerr <- (mo$area_moore_px - mo$n_px) / mo$n_px
sq <- data.frame(n_px = 2^(seq(2, 18, 0.25)))
sq$relerr <- ((sqrt(sq$n_px) - 1)^2 - sq$n_px) / sq$n_px
g3 <- ggplot(mo, aes(n_px, relerr, colour = id)) + geom_point(alpha = 0.6, size = 1) +
  geom_line(data = sq, aes(n_px, relerr), inherit.aes = FALSE, linetype = 2) +
  geom_hline(yintercept = c(-0.02, 0.02), colour = "grey40") + scale_x_log10() + theme_bw(base_size = 10) +
  labs(x = "pixel nella componente (log)", y = "errore relativo d'area (Moore sui centri)", colour = NULL,
       title = "CP-1: contorno di Moore sui centri vs area esatta",
       subtitle = "tratteggio = quadrato, (sqrt(n)-1)^2 / n - 1; linee grigie = ±2 % (soglia ROADMAP)")
ggsave(file.path(FIG, "S1.1_fig6_cp1_moore.png"), g3, width = 8, height = 4.5, dpi = 150, bg = "white")

# Fig 7 — CP-3 nullo
cp3 <- read.csv(file.path(RES, "S1.1_cp3_null.csv"))
c3l <- rbind(data.frame(dataset = cp3$dataset, mappa = "reale", metric = "regioni >= 100 µm²", value = cp3$n_regions_real),
             data.frame(dataset = cp3$dataset, mappa = "nullo (etichette permutate)", metric = "regioni >= 100 µm²", value = cp3$n_regions_null),
             data.frame(dataset = cp3$dataset, mappa = "reale", metric = "frazione di area esclusa", value = cp3$frac_excl_real),
             data.frame(dataset = cp3$dataset, mappa = "nullo (etichette permutate)", metric = "frazione di area esclusa", value = cp3$frac_excl_null))
g4 <- ggplot(c3l, aes(dataset, value, fill = mappa)) + geom_col(position = "dodge") + facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c("grey60", "steelblue")) + theme_bw(base_size = 10) + labs(x = NULL, y = NULL, fill = NULL,
  title = "CP-3: mappe graphclust reali vs nullo spaziale (stessa composizione)")
ggsave(file.path(FIG, "S1.1_fig7a_cp3_null.png"), g4, width = 9, height = 3.8, dpi = 150, bg = "white")
png_open("S1.1_fig7b_cp3_maps_A1.png", 2400, 1200); par(mfrow = c(1, 2), mar = c(1.5, 0.5, 2, 0.5))
draw_map("real_full_A1", "A1 reale", lwd = 0.05, border = NA); draw_map("null_full_A1", "A1 nullo (etichette permutate)", lwd = 0.05, border = NA)
dev.off()

# Fig 8 — costo della semplificazione
per <- read.csv(file.path(RES, "S1.1_per_input.csv"), stringsAsFactors = FALSE)
pp <- per[!is.na(per$c6_overlap_simpl), ]
pl <- rbind(data.frame(id = pp$id, group = pp$group, metric = "C1s: 95° percentile errore d'area per regione", value = pp$c1s_p95_region),
            data.frame(id = pp$id, group = pp$group, metric = "C6: sovrapposizione + scoperto / area", value = pp$c6_overlap_simpl + pp$c6_uncovered_simpl),
            data.frame(id = pp$id, group = pp$group, metric = "C5: 1 - concordanza raster", value = 1 - pp$c5_match_simpl))
thr <- data.frame(metric = unique(pl$metric), t = c(0.02, 0.01, 0.02))
g5 <- ggplot(pl, aes(value, reorder(id, value), colour = group)) + geom_point() + geom_vline(data = thr, aes(xintercept = t), linetype = 2) +
  facet_wrap(~metric, scales = "free_x") + theme_bw(base_size = 8) + labs(x = NULL, y = NULL, colour = NULL,
  title = "Costo della semplificazione (st_simplify, 0.5 µm); tratteggio = soglia pre-registrata")
ggsave(file.path(FIG, "S1.1_fig8_simplification_cost.png"), g5, width = 11, height = 9, dpi = 150, bg = "white")
cat("[figure] fatto\n")
