# 08_aim2_summary.R — descriptive first look at Aim 2 (CPU methods), written 2026-11-26 after the main
# runs of 05_deconvolve.R for colon and lymph node. Not the T5 analysis: no ranking, no H3 test.
# Inputs: <TESI_DATA>/sim/main/<tissue>_deconv_metrics.csv (one row per replicate x condition x
# dataset x method; failed spots already filled with the null baseline, DECISIONS 2026-11-03).
# Outputs in <TESI_DATA>/derived/aim2/:
#   aim2_summary.csv — mean RMSE, MCSE, failure rate, n by tissue x cells x dataset x method;
#   aim2_differences.csv — paired RMSE differences arm − anchor (same replicate), mean and MCSE.
# Figure in thesis/plots/ (small, tracked in git so it reaches the PC with a pull):
#   aim2_first_look.pdf — A4: (A) mean RMSE ± 2 MCSE by arm, null baseline dotted;
#                         (B) O − anchor and E′ − anchor, mean ± 2 MCSE.
# MCSE = SD / sqrt(n replicates) (Morris et al. 2019).
#
# Usage (working directory = thesis/):
#   Rscript scripts/08_aim2_summary.R [tissues, e.g. colon,lymphoid]

args <- commandArgs(trailingOnly = TRUE)
tissues <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else c("colon", "lymphoid")
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
out_dir <- file.path(data_dir, "derived", "aim2")
plot_dir <- "plots"
for (x in c(out_dir, plot_dir)) dir.create(x, recursive = TRUE, showWarnings = FALSE)

m <- do.call(rbind, lapply(tissues, function(t)
  read.csv(file.path(data_dir, "sim", "main", paste0(t, "_deconv_metrics.csv")))))
stopifnot(!anyNA(m$rmse))
tissue_lab <- c(colon = "Colon", lymphoid = "Lymph node", brain = "Brain")
arm_lev <- c("anchor", "O", "Eprime", "E", "P", "Cminus")
arm_lab <- c(anchor = "Anchor", O = "O", Eprime = "E′", E = "E", P = "P", Cminus = "C−")
meth_lev <- c("RCTD", "spatialDWLS", "MuSiC", "NNLS")
# Okabe–Ito colours (colour-blind safe); shapes give a second encoding of the method
meth_col <- c(RCTD = "#009E73", spatialDWLS = "#E69F00", MuSiC = "#0072B2", NNLS = "#CC79A7")
meth_shp <- c(RCTD = 16, spatialDWLS = 17, MuSiC = 15, NNLS = 18)

# Summary table
key <- c("tissue", "n_cells", "dataset", "method")
s <- aggregate(cbind(rmse, fail_rate) ~ tissue + n_cells + dataset + method, m, mean)
s$rmse_sd <- aggregate(rmse ~ tissue + n_cells + dataset + method, m, sd)$rmse
s$n <- aggregate(rmse ~ tissue + n_cells + dataset + method, m, length)$rmse
s$mcse <- s$rmse_sd / sqrt(s$n)
s <- s[order(s$tissue, s$n_cells, match(s$dataset, arm_lev), s$method), ]
write.csv(s, file.path(out_dir, "aim2_summary.csv"), row.names = FALSE)

# Paired differences arm − anchor (same replicate, condition and method)
a <- m[m$dataset == "anchor", c("tissue", "rep", "n_cells", "method", "rmse")]
names(a)[5] <- "rmse_anchor"
d <- merge(m[m$dataset != "anchor", c("tissue", "rep", "n_cells", "dataset", "method", "rmse")], a)
d$diff <- d$rmse - d$rmse_anchor
dd <- aggregate(diff ~ tissue + n_cells + dataset + method, d, mean)
dd$sd <- aggregate(diff ~ tissue + n_cells + dataset + method, d, sd)$diff
dd$n <- aggregate(diff ~ tissue + n_cells + dataset + method, d, length)$diff
dd$mcse <- dd$sd / sqrt(dd$n)
write.csv(dd, file.path(out_dir, "aim2_differences.csv"), row.names = FALSE)

# Plot style (after fig5_A4): white, left/bottom axes only, no grid, small sans text, bold titles
th <- theme_classic(base_size = 8, base_family = "sans") +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
        plot.tag = element_text(face = "bold", size = 12),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(colour = "black"),
        legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(3, "mm"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
lab_cells <- function(x) paste(x, "cells per spot")
prep <- function(z) {
  z$tissue <- factor(tissue_lab[z$tissue], levels = tissue_lab[tissues])
  z$cells <- factor(lab_cells(z$n_cells), levels = lab_cells(c(2, 10)))
  z
}
pd <- position_dodge(width = 0.6)

sa <- prep(s[s$method %in% meth_lev, ])
sa$arm <- factor(arm_lab[sa$dataset], levels = arm_lab[arm_lev])
sa$method <- factor(sa$method, levels = meth_lev)
nb <- prep(aggregate(rmse ~ tissue + n_cells, s[s$method == "null baseline", ], mean))   # one line per panel
pA <- ggplot(sa, aes(arm, rmse, colour = method, shape = method)) +
  geom_hline(data = nb, aes(yintercept = rmse), linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(ymin = rmse - 2 * mcse, ymax = rmse + 2 * mcse), width = 0, linewidth = 0.4, position = pd) +
  geom_point(size = 1.6, position = pd) +
  facet_grid(cells ~ tissue) +
  scale_colour_manual(values = meth_col) + scale_shape_manual(values = meth_shp) +
  labs(x = NULL, y = "RMSE of cell-type proportions",
       title = "Error by source of μ (dotted: null baseline)") + th

db <- prep(dd[dd$dataset %in% c("O", "Eprime") & dd$method %in% meth_lev, ])
db$arm <- factor(ifelse(db$dataset == "O", "O − anchor", "E′ − anchor"),
                 levels = c("O − anchor", "E′ − anchor"))
db$method <- factor(db$method, levels = meth_lev)
pB <- ggplot(db, aes(method, diff, colour = method, shape = arm)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(ymin = diff - 2 * mcse, ymax = diff + 2 * mcse), width = 0, linewidth = 0.4, position = pd) +
  geom_point(size = 1.8, position = pd, fill = "white", stroke = 0.6) +
  facet_grid(cells ~ tissue) +
  scale_colour_manual(values = meth_col, guide = "none") +
  scale_shape_manual(values = c(21, 16)) +
  labs(x = NULL, y = "Paired RMSE difference",
       title = "Reference arms against the anchor (same replicate)") + th

p <- (pA / pB) + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))
ggsave(file.path(plot_dir, "aim2_first_look.pdf"), p, width = 210, height = 260, units = "mm",
       device = if (capabilities("cairo")) cairo_pdf else pdf, bg = "white")
cat("Written:", file.path(out_dir, c("aim2_summary.csv", "aim2_differences.csv")),
    file.path(plot_dir, "aim2_first_look.pdf"), sep = "\n")
cat("Replicates per cell:", range(s$n), "\n")
