# tools/S1.4c_figures.R — S1.4c: figure del report (sfondo bianco, un pannello per archetipo). Uso: Rscript --vanilla tools/S1.4c_figures.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(ggplot2); library(arrow) })
RES <- "results/S1.4c"; FIG <- file.path(RES, "fig"); dir.create(FIG, showWarnings = FALSE); OUT <- "/mnt/micron/geo_spatialtrans/S1.4c"
th <- theme_bw(base_size = 11) + theme(plot.background = element_rect(fill = "white", colour = NA), strip.background = element_rect(fill = "grey92"))
col_v <- c("K-c" = "#B2182B", "K-1" = "#F4A582", "K-2" = "#2166AC")
fv <- read.csv(file.path(RES, "S1.4c_fold_verdicts.csv")); fv$variant <- factor(fv$variant, names(col_v))
lab <- c(B1 = "B-1 N/C", B2 = "B-2 A_N", B7 = "B-7 cerchio fuori", B8 = "B-8 cerchi sovrapposti", S1 = "S-1 a contatto*", S2 = "S-2 perimetro a contatto*", S3 = "S-3 solidita'*")
z <- fv[fv$metric %in% names(lab), ]; z$metric <- factor(lab[z$metric], lab); z$r <- z$dist_median / pmax(z$yard_median, 1e-4)
p1 <- ggplot(z, aes(metric, r, colour = variant)) + geom_hline(yintercept = 1, linetype = 2) +
  geom_point(position = position_dodge(width = 0.6), size = 1.6, alpha = 0.85) + scale_y_log10() + facet_wrap(~archetype, nrow = 2) +
  scale_colour_manual(values = col_v, name = NULL) +
  labs(x = NULL, y = "distanza simulato-reale / mediana del metro inter-ROI (log)",
       title = "Check B per fold (mediana di 10 semi): sotto la linea = PASS", subtitle = "* declassate a descrittive da CP-2 (dipendono dal segmentatore)") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(FIG, "S1.4c_fig2_checkB.png"), p1, width = 11, height = 7, dpi = 150, bg = "white")
o <- read.csv(file.path(RES, "S1.4c_obs_raw.csv")); rf <- read.csv(file.path(RES, "S1.4c_real_refs.csv"))
of <- aggregate(cbind(b7, b8, s1, m2) ~ archetype + variant + fold, o, median)
lv <- c(b7 = "B-7 cerchio di pari area fuori dal Voronoi", b8 = "B-8 cerchi di pari area sovrapposti", s1 = "S-1 nuclei a contatto (maschera, 2 px)", m2 = "D-1 nuclei con > 5 % di pixel fuori dal Voronoi")
sl <- do.call(rbind, lapply(names(lv), function(m) data.frame(archetype = of$archetype, variant = of$variant, metric = lv[[m]], value = of[[m]])))
rl <- do.call(rbind, lapply(names(lv), function(m) data.frame(archetype = rf$archetype, metric = lv[[m]], value = rf[[m]])))
sl$variant <- factor(sl$variant, names(col_v)); sl$metric <- factor(sl$metric, lv); rl$metric <- factor(rl$metric, lv)
p2 <- ggplot() + geom_point(data = sl, aes(archetype, value, colour = variant), position = position_dodge(width = 0.6), size = 1.6, alpha = 0.85) +
  geom_point(data = rl, aes(archetype, value), shape = 4, size = 2.2, colour = "black") + facet_wrap(~metric, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = col_v, name = "simulati (fold, mediana dei semi)") +
  labs(x = NULL, y = "frazione", title = "Contatto fra nuclei: simulati contro reali", subtitle = "x nere = ROI reali (segmentatore di forma)") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4c_fig3_contact.png"), p2, width = 10, height = 7, dpi = 150, bg = "white")
b4 <- do.call(rbind, lapply(list.files(file.path(OUT, "gen"), full.names = TRUE), function(f) { g <- readRDS(f)
  a <- aggregate(trunc_all ~ variant, g$b4, median); a$archetype <- g$archetype; a$fold <- g$fold; a }))
b4$variant <- factor(b4$variant, names(col_v))
p3 <- ggplot(b4, aes(archetype, trunc_all, colour = variant)) + geom_hline(yintercept = 0.05, linetype = 2) +
  geom_point(position = position_dodge(width = 0.6), size = 2) + scale_colour_manual(values = col_v, name = NULL) +
  labs(x = NULL, y = "frazione di nuclei troncati (per fold, mediana dei semi)", title = "B-4: troncamento dell'area estratta (soglia 5 %)") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4c_fig4_trunc.png"), p3, width = 7.5, height = 4.5, dpi = 150, bg = "white")
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
cp <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) do.call(rbind, lapply(c("cellpose_rgb", "stardist_he"), function(sg) {
  z <- as.data.frame(read_parquet(file.path(OUT, "cp2", sprintf("%s_%s_%s_cp2.parquet", rois$archetype[i], rois$roi_id[i], sg)))); z <- z[z$ev, ]
  data.frame(archetype = rois$archetype[i], roi = rois$roi_id[i], seg = sg, metric = c("S-1 nuclei a contatto", "S-3 solidita' (mediana)"),
             value = c(mean(z$n_c2 > 0), median(z$sol))) }))))
p4 <- ggplot(cp, aes(archetype, value, colour = seg, group = interaction(roi, seg))) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_wrap(~metric, scales = "free_y") + scale_colour_manual(values = c(cellpose_rgb = "#1B9E77", stardist_he = "#7570B3"), name = NULL) +
  labs(x = NULL, y = NULL, title = "CP-2: riferimenti di forma con due segmentatori sugli stessi ROI", subtitle = "un punto per ROI; stesso insieme di nuclei (centroide a >= 10 µm dal tessuto non valido)") +
  th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4c_fig5_cp2.png"), p4, width = 9, height = 4.5, dpi = 150, bg = "white")
cat("figure scritte in", FIG, "\n")
