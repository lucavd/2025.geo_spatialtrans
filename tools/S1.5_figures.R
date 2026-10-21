# tools/S1.5_figures.R - S1.5: figure del report (sfondo bianco, un pannello per tessuto dove pertinente).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(ggplot2); library(sf); library(png); library(grid) })
RES <- "results/S1.5"; FIG <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE)
TIS <- c(A1 = "intestino", A2 = "tumore del colon", A3 = "stroma", A4 = "linfonodo", A5 = "cervello", A6 = "cuore")
COL <- c(confermata = "#2b6cb0", smentita = "#dd6b20", indeterminata = "#718096")
th <- theme_bw(base_size = 9) + theme(plot.background = element_rect(fill = "white", colour = NA), panel.grid.minor = element_blank(),
                                      plot.title = element_text(size = 9, hjust = 0), strip.background = element_rect(fill = "grey95"))
V <- read.csv(file.path(RES, "S1.5_verdicts.csv")); V <- V[V$winner_arm & V$id %in% paste0("P-", 1:7), ]
V$tissue <- factor(V$tissue, levels = TIS); V$x <- as.numeric(V$tissue)
lab <- c(`P-1` = "P-1 area mediana del nucleo (µm²)", `P-2` = "P-2 rapporto nucleo/territorio", `P-3` = "P-3 nuclei troncati",
         `P-4` = "P-4 cellule quasi senza citoplasma", `P-5` = "P-5 nuclei a contatto (finestre)", `P-6` = "P-6 CV dei territori, sim/reale - 1",
         `P-7` = "P-7 eccentricità, sim - reale")
V$panel <- factor(lab[V$id], levels = lab)
p1 <- ggplot(V) + geom_rect(aes(xmin = x - 0.35, xmax = x + 0.35, ymin = int_lo, ymax = int_hi), fill = "grey85") +
  geom_errorbar(aes(x = x, ymin = lo2, ymax = hi2, colour = outcome), width = 0.25) + geom_point(aes(x = x, y = estimate, colour = outcome), size = 1.6) +
  facet_wrap(~panel, scales = "free_y", ncol = 4) + scale_x_continuous(breaks = 1:6, labels = TIS) + scale_colour_manual(values = COL, name = "esito") +
  labs(x = NULL, y = NULL, title = "Simulatore completo contro i riferimenti: stima ± 2 errori standard (punto e barra), intervallo di rilevanza (grigio)") + th +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = c(0.88, 0.2))
ggsave(file.path(FIG, "S1.5_fig1_predictions.png"), p1, width = 10, height = 5.2, dpi = 150, bg = "white")
# PH-1: contatti per braccio
w <- read.csv(file.path(RES, "S1.5_posthoc_windows.csv")); ps <- aggregate(s1 ~ arm + archetype + seed, w, mean)
ag <- do.call(rbind, lapply(split(ps, list(ps$arm, ps$archetype)), function(d) data.frame(arm = d$arm[1], archetype = d$archetype[1], m = mean(d$s1), se = sd(d$s1) / sqrt(nrow(d)))))
ag$tissue <- factor(TIS[ag$archetype], levels = TIS)
ag$arm <- factor(c(realgen = "generatori veri", sim_segdens = "semina, densità dei segmentatori", preserve_area = "semina, densità manuale")[ag$arm],
                 levels = c("generatori veri", "semina, densità dei segmentatori", "semina, densità manuale"))
B <- data.frame(tissue = factor(TIS, levels = TIS), lo = c(0.366, 0.217, 0.083, 0.350, 0, 0), hi = c(0.650, 0.567, 0.283, 0.750, 0.217, 0.070))
p2 <- ggplot(ag) + geom_rect(data = B, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi), fill = "grey85", inherit.aes = FALSE) +
  geom_pointrange(aes(x = arm, y = m, ymin = m - 2 * se, ymax = m + 2 * se, colour = arm), size = 0.3) + facet_wrap(~tissue, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("#2b6cb0", "#9f7aea", "#dd6b20"), name = NULL) +
  labs(x = NULL, y = "frazione di nuclei a contatto", title = "I contatti in eccesso vengono dalla semina omogenea, non dalla densità: con i generatori veri tornano nell'intervallo manuale (grigio)") +
  th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
ggsave(file.path(FIG, "S1.5_fig2_contacts.png"), p2, width = 10, height = 3.4, dpi = 150, bg = "white")
# esempio: H&E della finestra e livello simulato della stessa finestra (seme 1)
cr <- read.csv("/mnt/micron/geo_spatialtrans/S1.5/fig/crops.csv"); gl <- list()
for (k in seq_len(nrow(cr))) {
  A <- cr$archetype[k]; L <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.5/layers/%s_%s_preserve_area.rds", A, cr$roi_id[k]))
  bb <- st_as_sfc(st_bbox(c(xmin = cr$x0_um[k], ymin = cr$y0_um[k], xmax = cr$x0_um[k] + cr$side_um[k], ymax = cr$y0_um[k] + cr$side_um[k])))
  keep <- lengths(st_intersects(L$cell_territories, bb)) > 0
  Tt <- st_intersection(st_sf(geometry = L$cell_territories[keep]), bb); Nn <- st_intersection(st_sf(pale = L$cell_df$nucleus_pale[keep], geometry = L$cell_nuclei[keep]), bb)
  img <- readPNG(cr$png[k])
  gh <- ggplot() + annotation_raster(img, 0, 1, 0, 1) + coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void() +
    labs(title = sprintf("%s, H&E (%.0f µm)", TIS[[A]], cr$side_um[k])) + theme(plot.title = element_text(size = 8))
  gs <- ggplot() + geom_sf(data = Tt, fill = "#f7e8ef", colour = "grey40", linewidth = 0.15) +
    geom_sf(data = Nn, aes(fill = pale), colour = NA) + scale_fill_manual(values = c(`FALSE` = "#3c366b", `TRUE` = "#b794f4"), guide = "none") +
    scale_y_reverse() + coord_sf(expand = FALSE, datum = NA) + theme_void() + labs(title = "simulato, seme 1") + theme(plot.title = element_text(size = 8))
  gl[[length(gl) + 1]] <- gh; gl[[length(gl) + 1]] <- gs
}
png(file.path(FIG, "S1.5_fig3_example.png"), width = 12, height = 4.4, units = "in", res = 150, bg = "white")
grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 6)))
for (k in 1:6) { print(gl[[2 * k - 1]], vp = viewport(layout.pos.row = 1, layout.pos.col = k)); print(gl[[2 * k]], vp = viewport(layout.pos.row = 2, layout.pos.col = k)) }
dev.off()
# CP-1: rapporto stimata / vera
K <- rbind(transform(read.csv(file.path(RES, "S1.5_CP1_scale.csv")), check = "come pre-registrato"),
           transform(read.csv(file.path(RES, "S1.5_CP1_scale_est_b.csv")), check = "con controllo (post hoc)"))
K <- K[K$variant == "pipeline", ]; K$tissue <- factor(TIS[K$archetype], levels = TIS); K$esito <- ifelse(K$refused %in% TRUE, "rifiuto", "numero restituito")
K$y <- ifelse(K$refused %in% TRUE, 0.03, K$ratio)
p4 <- ggplot(K, aes(x = factor(scale_true), y = y)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.2, fill = "grey85") +
  geom_jitter(aes(shape = esito, colour = esito), width = 0.15, height = 0, size = 1.4) + facet_grid(check ~ tissue) +
  scale_y_log10(breaks = c(0.03, 0.1, 0.3, 1, 3), labels = c("rifiuto", "0.1", "0.3", "1", "3")) + scale_shape_manual(values = c(`numero restituito` = 16, rifiuto = 4)) +
  scale_colour_manual(values = c(`numero restituito` = "#dd6b20", rifiuto = "#718096")) +
  labs(x = "scala vera (µm/px)", y = "scala stimata / vera", title = "La stima della scala dall'immagine ritrova la scala solo nel linfonodo e a volte restituisce numeri sbagliati (5 ROI per tessuto)") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.5_fig4_scale.png"), p4, width = 10, height = 4.6, dpi = 150, bg = "white")
cat("figure scritte\n")
