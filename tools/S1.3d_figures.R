# tools/S1.3d_figures.R — S1.3d: figure del report (sfondo bianco, un pannello per archetipo dove pertinente).
# Uso: Rscript --vanilla tools/S1.3d_figures.R [outdir] ; legge results/S1.3d (o S13D_OUTR), scrive results/S1.3d/figures/.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(ggplot2); library(sf); library(arrow); library(jsonlite); library(polyclip) })
a <- commandArgs(TRUE); OUTD <- if (length(a)) a[1] else "/mnt/micron/geo_spatialtrans/S1.3d"
IN <- Sys.getenv("S13D_OUTR", "results/S1.3d"); FIG <- file.path(IN, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
G_GRID <- 1e-6; source("tools/S1.3d_candidates.R")
theme_set(theme_bw(base_size = 10) + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white", colour = NA)))
CC <- c(V0 = "grey45", F1 = "#1b9e77", F2 = "#7570b3", F3 = "#e7298a", F4 = "#66a61e", F5 = "#d95f02", F6 = "#e6ab02")
O1 <- read.csv(file.path(IN, "S1.3d_o1_long.csv")); O2 <- read.csv(file.path(IN, "S1.3d_o2_long.csv"))
fl <- function(x) pmax(x, 1e-12)
# ---- Fig. 1: O1, massimo errore/tau per nullo (N), per archetipo ---------------------------------------------------
d <- O1[O1$set == "N" & O1$applicable %in% TRUE, ]
d$ratio <- pmax(d$pair_max_ratio, d$cell_max_ratio, d$union_err / d$tauU, d$symd_err / d$tauU, na.rm = TRUE)
d$cand <- factor(d$cand, levels = names(CC))
p <- ggplot(d, aes(cand, fl(ratio), colour = cand)) + geom_hline(yintercept = 1, linetype = 2) + geom_jitter(width = 0.25, height = 0, size = 0.5, alpha = 0.6) +
  scale_y_log10() + scale_colour_manual(values = CC, guide = "none") + facet_wrap(~archetype, nrow = 2) +
  labs(x = NULL, y = "max errore / τ per nullo (unione, diff. simm., coppie, celle)", title = "O1 sui 1 800 nulli: errore massimo rispetto alla tolleranza")
ggsave(file.path(FIG, "fig1_o1_ratio.png"), p, width = 9, height = 5.5, dpi = 150)
# ---- Fig. 2: O2 (2 µm) e O3, massimo errore/tau per caso (insieme O2 di N), per archetipo ------------------------------------
d2 <- O2[O2$set == "N" & O2$op %in% c("O2", "O3") & !(O2$op == "O2" & O2$s == 8), ]; d2$what <- ifelse(d2$op == "O2", "O2 (2 µm)", "O3 (2→8 µm)")
d2$cand <- factor(d2$cand, levels = names(CC))
p <- ggplot(d2, aes(cand, fl(max_ratio), colour = cand, shape = what)) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = 1 / sqrt(2), linetype = 3) +
  geom_jitter(width = 0.25, height = 0, size = 0.9, alpha = 0.7) + scale_y_log10() + scale_colour_manual(values = CC, guide = "none") + facet_wrap(~archetype, nrow = 2) +
  labs(x = NULL, y = "max errore / τ per caso", shape = NULL, title = "O2 e O3 sui 120 nulli dell'insieme O2 (tratteggio: τ; punteggiato: 1/√2, limite dell'arrotondamento)")
ggsave(file.path(FIG, "fig2_o2_ratio.png"), p, width = 9, height = 5.5, dpi = 150)
# ---- Fig. 3: copertura V0 e fallimenti V0 per archetipo -------------------------------------------------------------------
v <- O1[O1$set == "N" & O1$cand == "V0", ]
tb <- do.call(rbind, lapply(split(v, v$archetype), function(x) data.frame(archetype = x$archetype[1],
      what = c("copertura V0 invalida", "errore silenzioso V0 (O1)", "errore con copertura valida"),
      frac = c(mean(x$v0_valid %in% FALSE), mean(x$silent %in% TRUE), mean(x$silent %in% TRUE & x$v0_valid %in% TRUE)))))
p <- ggplot(tb, aes(archetype, frac, fill = what)) + geom_col(position = position_dodge(0.8), width = 0.75) + scale_fill_manual(values = c("grey60", "#d95f02", "black")) +
  labs(x = NULL, y = "frazione dei 300 nulli", fill = NULL, title = "V0: copertura invalida (GEOS CoverageValidator) e errori silenziosi di O1") + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig3_v0_coverage.png"), p, width = 8, height = 4.2, dpi = 150)
# ---- Fig. 4: prestazioni su S ------------------------------------------------------------------------------------------------
pf <- file.path(IN, "S1.3d_perf.csv")
if (file.exists(pf)) { P <- read.csv(pf); P$cand <- factor(P$cand, levels = names(CC))
  L <- rbind(data.frame(P[, c("cand", "op")], what = "tempo (s, mediana di 3)", val = P$t_total), data.frame(P[, c("cand", "op")], what = "memoria di picco (GB)", val = P$rss_gb))
  p <- ggplot(L, aes(cand, val, fill = cand)) + geom_col() + scale_fill_manual(values = CC, guide = "none") + facet_grid(what ~ op, scales = "free_y") +
    labs(x = NULL, y = NULL, title = "S (449 536 celle, 16 mm²): preparazione + operazione, processo separato")
  ggsave(file.path(FIG, "fig4_perf.png"), p, width = 9, height = 5, dpi = 150) }
# ---- Fig. 5: un errore di F5 (polyclip) e una cella F4 invalida -------------------------------------------------------------
case_plot <- function(key, i, ix, iy, s, title) {
  m <- fromJSON(file.path(OUTD, "geom", paste0(key, "_meta.json"))); t <- df_to_sfc(as.data.frame(read_parquet(file.path(OUTD, "geom", paste0(key, "_cells.parquet")))), m$n)
  x0 <- m$gs$X0 + s * ix; y0 <- m$gs$Y0 + s * iy; b <- st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
  ex <- st_intersection(t[i], st_sfc(b)); R <- PC(rings_of(t[[i]]), list(list(x = c(x0, x0 + s, x0 + s, x0), y = c(y0, y0, y0 + s, y0 + s))), "intersection")
  f5 <- st_sfc(lapply(R, function(r) st_polygon(list(rbind(cbind(r$x, r$y), c(r$x[1], r$y[1]))))))
  ae <- sum(as.numeric(st_area(ex))); a5 <- pc_area(R)
  ggplot() + geom_sf(data = st_sf(geometry = t[i]), fill = "grey90", colour = "grey40") + geom_sf(data = st_sf(geometry = st_sfc(b)), fill = NA, colour = "black", linetype = 2) +
    geom_sf(data = st_sf(geometry = ex), fill = "#1b9e77", alpha = 0.5, colour = NA) + geom_sf(data = st_sf(geometry = f5), fill = NA, colour = "#d95f02", linewidth = 0.8) +
    coord_sf(xlim = c(x0 - 0.6 * s, x0 + 1.6 * s), ylim = c(y0 - 0.6 * s, y0 + 1.6 * s), datum = NA) +
    labs(title = title, subtitle = sprintf("pezzo esatto %.6f µm² (verde) · polyclip %.6f µm² (arancio)", ae, a5), x = NULL, y = NULL)
}
v5 <- O2[O2$cand == "F5" & O2$op == "O2" & O2$s == 2 & O2$fail %in% TRUE, ]
if (nrow(v5)) { best <- NULL
  for (k in setdiff(unique(v5$key), "S_cperf")) { r <- readRDS(list.files(file.path(OUTD, "res"), pattern = paste0("^", k, "_o12?\\.rds$"), full.names = TRUE)[1]); vv <- r$viol[["o2_F5_2"]]
    if (is.null(vv)) next; vv$at <- ifelse(is.na(vv$area_t), 0, vv$area_t); q <- vv[which.max(abs(vv$area - vv$at)), ]; q$key <- k; q$d <- abs(q$area - q$at)
    if (is.null(best) || q$d > best$d) best <- q }
  kk <- as.integer(strsplit(best$k, " ")[[1]])
  p <- case_plot(best$key, kk[1], kk[2], kk[3], 2, sprintf("F5 su %s (input valido): cella %d, bin (%d, %d)", best$key, kk[1], kk[2], kk[3]))
  ggsave(file.path(FIG, "fig5_f5_case.png"), p, width = 5.5, height = 5.5, dpi = 150) }
c5f <- list.files(file.path(OUTD, "res"), pattern = "_o12?\\.rds$", full.names = TRUE)
for (f in c5f) { r <- readRDS(f); if (!is.null(r$c5$info$n_invalid) && r$c5$info$n_invalid > 0) {
  m <- fromJSON(file.path(OUTD, "geom", paste0(r$key, "_meta.json"))); t4 <- df_to_sfc(as.data.frame(read_parquet(file.path(OUTD, "geom", paste0(r$key, "_f4cells.parquet")))), m$n)
  reg <- df_to_sfc(transform(as.data.frame(read_parquet(file.path(OUTD, "geom", paste0(r$key, "_region.parquet")))), cell = 1L), 1L)
  b <- which(!st_is_valid(t4))[1]; rs <- st_is_valid(t4[b], reason = TRUE); xy <- as.numeric(regmatches(rs, gregexpr("[-0-9.]+", rs))[[1]]); xy <- tail(xy, 2)
  w <- 6; p <- ggplot() + geom_sf(data = st_sf(geometry = reg), fill = "grey95", colour = "grey30") + geom_sf(data = st_sf(geometry = t4[b]), fill = "#66a61e", alpha = 0.4, colour = "#66a61e") +
    annotate("point", x = xy[1], y = xy[2], colour = "red", size = 2) + coord_sf(xlim = xy[1] + c(-w, w), ylim = xy[2] + c(-w, w), datum = NA) +
    labs(title = sprintf("F4 su %s: cella %d invalida", r$key, b), subtitle = rs, x = NULL, y = NULL)
  ggsave(file.path(FIG, "fig6_f4_invalid.png"), p, width = 5.5, height = 5.5, dpi = 150); break } }
cat("figure scritte in", FIG, "\n"); print(list.files(FIG))
