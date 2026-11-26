# 05b_diag_rctd_spots.R — diagnostic (read-only), written 2026-11-26 after the lymphoid main run of
# 05_deconvolve.R: RCTD failed on 75% of the spots of arm P at 2 cells (other arms <= 4%).
# Question: does RCTD drop those spots before fitting? create.RCTD() (spacexr) keeps a spot only if
# nUMI >= UMI_min (100) and it has >= counts_MIN (10) counts on gene_list_bulk, the genes it selects
# from the reference as differing between types (restrict_counts()). A dropped spot has no estimate,
# so 05 counts it as failed. Expectation: in P the selected genes carry a smaller share of the spot UMI
# than in the real-gene arms, so more spots fall below 10 counts.
# Runs only create.RCTD() (no fit) on saved datasets with the references of 05; writes nothing.
# Also prints the RCTD failure rate by arm from <tissue>_deconv_metrics.csv.
#
# Usage (working directory = thesis/):
#   Rscript scripts/05b_diag_rctd_spots.R <tissue> [reps, e.g. 1,2,3]

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "lymphoid"
reps <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else 1:3
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
suppressPackageStartupMessages(library(Matrix))
for (f in c("R/simulate_spots.R", "R/arm_p.R", "R/deconvolve.R")) source(f)
R0 <- readRDS(file.path(data_dir, "derived", "reference", sprintf("%s_main.rds", tissue)))
sim_dir <- file.path(data_dir, "sim", "main", tissue)
L_ref <- Matrix::colSums(R0$ref)
hd_depth <- c(colon = 174, lymphoid = 62, brain = 121)
keys <- sprintf("c%02d_d%04d", c(2, 10), hd_depth[[tissue]])
options(width = 150)

# 1. RCTD failure rate by arm and cells (all replicates), and spread over replicates
m <- read.csv(file.path(dirname(sim_dir), paste0(tissue, "_deconv_metrics.csv")))
r <- m[m$method == "RCTD", ]
cat("RCTD failure rate, mean over replicates\n")
print(round(tapply(r$fail_rate, list(r$dataset, r$n_cells), mean), 3))
cat("RCTD failure rate in P, quantiles over replicates\n")
for (nc in c(2, 10)) { cat(nc, "cells: "); print(round(quantile(r$fail_rate[r$dataset == "P" & r$n_cells == nc]), 3)) }

# 2. create.RCTD() on saved datasets: selected genes, their UMI share, spots dropped
one <- function(spots, ref, ref_type) {
  ct <- factor(ref_type); names(ct) <- colnames(ref)
  reference <- spacexr::Reference(ref, ct, nUMI = setNames(Matrix::colSums(ref), colnames(ref)))
  puck <- spacexr::SpatialRNA(NULL, spots, nUMI = setNames(Matrix::colSums(spots), colnames(spots)),
                              use_fake_coords = TRUE)
  obj <- suppressMessages(spacexr::create.RCTD(puck, reference, max_cores = 1, UMI_min_sigma = 100))
  g <- obj@internal_vars$gene_list_bulk
  on_g <- Matrix::colSums(spots[g, , drop = FALSE]); tot <- Matrix::colSums(spots)
  data.frame(n_genes_bulk = length(g), n_genes_reg = length(obj@internal_vars$gene_list_reg),
             spot_umi = median(tot), umi_on_genes = median(on_g), share_on_genes = median(on_g / tot),
             below10 = mean(on_g < 10), kept = ncol(obj@spatialRNA@counts) / ncol(spots))
}
out <- NULL
for (i in reps) for (key in keys) {
  rd <- file.path(sim_dir, sprintf("r%03d", i))
  x <- readRDS(file.path(rd, paste0(key, ".rds")))
  RP <- reference_p_from_world(readRDS(file.path(rd, "world_p.rds")), R0$ref_meta, L_ref, colnames(R0$ref))
  est <- readRDS(file.path(rd, paste0("est_", key, ".rds")))$metrics
  for (d in names(x$counts)) {
    ref <- if (d == "P") RP$ref else R0$ref
    z <- tryCatch(one(x$counts[[d]], ref, R0$ref_meta$cell_type),
                  error = function(e) { cat("error", i, key, d, conditionMessage(e), "\n"); NULL })
    if (is.null(z)) next
    fr <- est$fail_rate[est$dataset == d & est$method == "RCTD"]
    out <- rbind(out, cbind(rep = i, cond = key, dataset = d, z, rctd_fail_rate = fr))
  }
}
cat("\ncreate.RCTD() on saved datasets (medians over spots)\n")
print(out, digits = 3, row.names = FALSE)
cat("\nMeans over replicates\n")
print(aggregate(cbind(n_genes_bulk, share_on_genes, below10, kept, rctd_fail_rate) ~ cond + dataset, out, mean), digits = 3)
