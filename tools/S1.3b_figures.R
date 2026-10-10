# tools/S1.3b_figures.R — figure del report S1.3b (sfondo bianco; un pannello per archetipo dove pertinente).
REAL_ARGS <- commandArgs(TRUE)
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
suppressMessages(source("tools/S1.3_roi.R")); rm(commandArgs); source("tools/S1.3b_overlay.R")
suppressPackageStartupMessages(library(ggplot2))
src <- readLines("tools/S1.3b_run.R"); eval(parse(text = src[grep("^regen <- function", src):(grep("^reg_features <- function", src) - 1)]))
RES <- "results/S1.3b"; FIG <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11) + theme(plot.background = element_rect(fill = "white", colour = NA)))
tab <- read.csv(file.path(RES, "S1.3b_null_table.csv"))
VL <- c("V0", "V0g", "V1", "V2", "V3")
# Fig 1 — errori V0 per archetipo e variante
f1 <- do.call(rbind, lapply(VL, function(v) { d <- tab[tab$variant == v, ]
  u <- if (v %in% c("V0", "V0g")) d$ov / d$A > 1e-9 | d$symd / d$A > 1e-9 else abs(d$ov) > d$tau | d$symd > d$tau
  rbind(data.frame(variant = v, archetype = d$archetype, tipo = "unione sbagliata (nulli)", n = as.integer(u)),
        data.frame(variant = v, archetype = d$archetype, tipo = "intersezione di coppia sbagliata (coppie)", n = d$pairs_over)) }))
f1 <- aggregate(n ~ variant + archetype + tipo, f1, sum); f1$variant <- factor(f1$variant, VL)
p1 <- ggplot(f1, aes(variant, n, fill = tipo)) + geom_col(position = position_dodge(0.8), width = 0.75) + facet_wrap(~archetype, nrow = 1) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) + labs(x = NULL, y = "conteggio su 300 nulli per archetipo", fill = NULL,
  title = "Errori silenziosi di sovrapposizione per variante (1 800 nulli di S1.3)", subtitle = "V0/V0g: soglie 1e-9·A (unione) e 1e-6·(Pi+Pj) (coppie); V1–V3: tau pre-registrata") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG, "S1.3b_fig1_errori_per_variante.png"), p1, width = 11, height = 4.2, dpi = 150, bg = "white")
# Fig 2 — la coppia 766/1987 (A3 r3 CSR 3): celle e intersezione dichiarata da V0
G <- regen(which(rois$archetype == "A3" & rois$roi_id == "r3"), "CSR", 3L); t <- G$terr
g0 <- st_intersection(t[766], t[1987]); V2 <- s3b_variant(t, G$reg, "V2"); g2 <- st_as_sfc(as_wkb(geos_intersection_prec(V2$terr[766], V2$terr[1987], 1e-6)))
bb <- st_bbox(c(t[766], t[1987])); nb <- t[lengths(st_intersects(t, st_as_sfc(bb))) > 0]
d1 <- rbind(st_sf(ruolo = "territori vicini", geometry = nb), st_sf(ruolo = "cella 766", geometry = t[766]), st_sf(ruolo = "cella 1987", geometry = t[1987]))
d1$ruolo <- factor(d1$ruolo, c("territori vicini", "cella 766", "cella 1987"))
v12 <- s3b_geoms_df(list(t[[766]], t[[1987]])); fu <- s3b_fuse(v12$x, v12$y, 1e-6); cr <- v12$g[fu$I] != v12$g[fu$J]
sh <- data.frame(x = v12$x[fu$I[cr]], y = v12$y[fu$I[cr]], tipo = ifelse(fu$d[cr] == 0, "vertice identico", sprintf("vertice diverso (%.1e µm)", fu$d[cr])))
p2a <- ggplot() + geom_sf(data = d1, aes(fill = ruolo), colour = "grey30", linewidth = 0.2) +
  geom_sf(data = st_sf(geometry = st_collection_extract(g0, "POLYGON")), fill = NA, colour = "red", linewidth = 0.9, linetype = "22") +
  geom_point(data = sh, aes(x, y, shape = tipo), size = 3) + scale_shape_manual(values = c(1, 4)) +
  scale_fill_manual(values = c("grey92", "#56B4E9", "#E69F00")) + coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], datum = NA) +
  labs(title = sprintf("V0: st_intersection(766, 1987) = %.2f µm² (tratteggio rosso)", sum(as.numeric(st_area(g0)))),
       subtitle = sprintf("V2 (snap-rounding 1e-6 µm): %.2g µm²", sum(as.numeric(st_area(g2)))), fill = NULL, shape = NULL, x = NULL, y = NULL)
ggsave(file.path(FIG, "S1.3b_fig2_coppia_766_1987.png"), p2a, width = 7.5, height = 6, dpi = 150, bg = "white")
# Fig 3 — C-5: max |Δa| / (g·P_i) per nullo
f3 <- tab[tab$variant %in% c("V1", "V2", "V3"), c("archetype", "variant", "c5_ratio")]
n0 <- tapply(f3$c5_ratio == 0, f3$variant, sum)
p3 <- ggplot(f3, aes(variant, c5_ratio)) + geom_boxplot(outlier.size = 0.6) + geom_hline(yintercept = 1, colour = "red", linetype = 2) +
  facet_wrap(~archetype, nrow = 1) + scale_y_log10() + labs(x = NULL, y = "max |Δa| / (g·P_i) per nullo (log)",
  title = "C-5: spostamento d'area per cella rispetto alla tolleranza (linea rossa = 1)", subtitle = sprintf("V1, V2: g = 1e-6 µm; V3: g = 1e-7 µm. Nulli senza alcuna cella spostata (rapporto 0, fuori scala): V1 %d, V2 %d, V3 %d su 1 800", n0[["V1"]], n0[["V2"]], n0[["V3"]]))
ggsave(file.path(FIG, "S1.3b_fig3_c5.png"), p3, width = 11, height = 3.8, dpi = 150, bg = "white")
# Fig 4 — CP-2: FAIL per δ e variante, silenziosi vs eccezioni
ct <- read.csv(file.path(RES, "S1.3b_cp2.csv"))
f4 <- rbind(data.frame(variant = ct$variant, delta = ct$delta, tipo = "eccezione GEOS (input invalido)", n = as.integer(!is.na(ct$error))),
            data.frame(variant = ct$variant, delta = ct$delta, tipo = "errore silenzioso (C-1)", n = as.integer(ct$fail & is.na(ct$error))))
f4 <- aggregate(n ~ variant + delta + tipo, f4, sum); f4$delta <- factor(sprintf("δ = %g µm", f4$delta), sprintf("δ = %g µm", c(1e-12, 1e-10, 1e-8)))
p4 <- ggplot(f4, aes(variant, n, fill = tipo)) + geom_col() + facet_wrap(~delta, nrow = 1) + scale_fill_manual(values = c("grey50", "#D55E00")) +
  geom_hline(yintercept = 60 * 15 / 1800, linetype = 3) + labs(x = NULL, y = "nulli con FAIL (su 60)", fill = NULL,
  title = "CP-2: stress con vertici perturbati indipendentemente", subtitle = "linea punteggiata = basale atteso dagli errori V0 (15/1 800 × 60 = 0.5)") + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.3b_fig4_cp2.png"), p4, width = 8, height = 3.8, dpi = 150, bg = "white")
cat("figure fatte\n")
