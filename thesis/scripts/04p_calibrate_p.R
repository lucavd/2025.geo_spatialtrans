# 04p_calibrate_p.R — calibration of arm P on the main training split, once per tissue
# (ADEMP §2.1 "P parameters"; DECISIONS 2026-11-03 and 2026-11-04, Pierangela).
#
# Pre-registered settings (written 2026-11-04 before the first run):
#   - data: main training split only (split == "training"), cell types with >= 50 training cells
#     (type_kept); genes = ADEMP §2.6.2 filter (count > 1 in >= 10 training cells); the held-out group
#     is never read;
#   - markers: scuttle::logNormCounts, then scran::findMarkers(groups = cell type, test.type = "t" (Welch),
#     direction = "up", lfc = 1, pval.type = "some", block = donor); marker of type k = FDR < 0.05;
#     (changed 2026-11-04 after the first colon run: pval.type = "all" gave 0 markers to 7 of 18 types,
#     3 of them spot types, whose P profiles would be identical; OSCA 11.2.2.3 suggests "some" as the
#     compromise; DECISIONS 2026-11-04, Pierangela);
#     check: number of markers per type also at FDR < 0.01;
#   - baseline means (input kept from ADEMP v1.0): pooled relative profile of all training cells (sum of
#     counts of gene g / sum of all counts, as mu_E in ADEMP §2.3), in counts per 10,000, on the genes that
#     are a marker of no type;
#   - estimation (refined 2026-11-04, Pierangela): the procedure of Splatter (Zappia et al. 2017; source
#     R/splat-estimate.R, splatter 1.37.1, functions splatEstMean and splatEstOutlier), applied to these
#     means: (a) gamma: means winsorized at the 10th/90th percentiles, then fitdistrplus::fitdist(method =
#     "mge", gof = "CvM"), method of moments if it does not converge; (b) expression outliers: genes with
#     log mean > median + 2 MAD of the log means; out.prob = their share; factors = mean / median of the
#     means, log-normal fit (fitdistrplus, MLE) -> out.facLoc, out.facScale. Difference from Splatter,
#     declared: Splatter averages library-size-normalised cells (each cell scaled to the median library
#     size), we use the pooled profile, because the generator is defined on it (ADEMP §2.3);
#     in a replicate (Splatter splatSimGeneMeans): base mean ~ gamma; with probability out.prob a gene's mean
#     becomes median(base means) x a log-normal factor; then markers x fold-change; columns normalised to 1;
#   - fold-change of marker g of type k: mu_E[g, k] / mean over the other types of mu_E[g, j], with mu_E the
#     pooled profile of each type (ADEMP §2.3), the same definition as the top-k markers of ADEMP §5.1;
#   - checks printed: cells, donors, types, genes; markers per type at 0.05 and 0.01; gamma and outlier
#     parameters with empirical vs simulated quantiles of the gene means (10, 50, 90, 99, 99.9%, from one
#     draw of the fitted model); fold-change summary and share of fold-changes below 1 (expected ~0, since
#     markers are tested as up against every other type).
# Output: <TESI_DATA>/derived/p_calibration/<tissue>_p_calibration.rds (markers, n per type, gamma
# parameters, fold-changes per type, gene list, settings, session info).
#
# Usage (working directory = thesis/):  Rscript scripts/04p_calibrate_p.R colon

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
stopifnot(tissue %in% c("brain", "colon", "lymphoid"))
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
out_dir <- file.path(data_dir, "derived", "p_calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source("R/read_h5ad.R"); source("R/stream_counts.R"); source("R/arm_p.R")
suppressMessages({ library(SingleCellExperiment); library(scran); library(scuttle) })

seed <- 42               # study seed (Pierangela, 2026-11-03); findMarkers is deterministic
gene_min_cells <- 10     # ADEMP §2.6.2
fdr_main <- 0.05; fdr_check <- 0.01
counts_group <- c(brain = "/X", colon = "/raw/X", lymphoid = "/raw/X")[[tissue]]
t0 <- proc.time()[["elapsed"]]
el <- function() (proc.time()[["elapsed"]] - t0) / 60
set.seed(seed)

# --- training cells and counts (same reading as 04m_fit_scdesign3.R) -------------------------------
cells <- read.csv(file.path(data_dir, "derived", "split", paste0(tissue, "_cells.csv.gz")),
                  stringsAsFactors = FALSE)
cells <- cells[cells$split == "training" & cells$type_kept, ]
cells <- cells[order(cells$row), ]
cat(sprintf("[%s] training: %d cells, %d donors, %d types\n", tissue, nrow(cells),
            length(unique(cells$donor_id)), length(unique(cells$cell_type))))
file <- list.files(file.path(data_dir, "raw"), pattern = paste0("^", tissue, "_.*\\.h5ad$"), full.names = TRUE)
genes <- read_var_genes(file, counts_group)
mats <- list()
stream_counts(file, counts_group, cells$row, nrow(genes), function(m, r) mats[[length(mats) + 1]] <<- m)
X <- Matrix::t(do.call(rbind, mats)); rm(mats); invisible(gc())   # genes x cells
rownames(X) <- genes$gene_id; colnames(X) <- cells$cell_id
keep <- Matrix::rowSums(X > 1) >= gene_min_cells
X <- X[keep, ]
cat(sprintf("  read in %.1f min; genes kept %d of %d\n", el(), nrow(X), length(keep)))
stopifnot(all(X@x >= 0), all(X@x == round(X@x)))

# --- markers --------------------------------------------------------------------------------------
sce <- SingleCellExperiment(list(counts = X),
                            colData = data.frame(cell_type = cells$cell_type, donor = cells$donor_id))
sce <- logNormCounts(sce)
mk <- findMarkers(sce, groups = sce$cell_type, test.type = "t", direction = "up", lfc = 1,
                  pval.type = "some", block = sce$donor)
types <- names(mk)
markers <- lapply(mk, function(d) rownames(d)[!is.na(d$FDR) & d$FDR < fdr_main])
n_markers <- data.frame(cell_type = types,
                        n_fdr05 = vapply(mk, function(d) sum(d$FDR < fdr_main, na.rm = TRUE), 1L),
                        n_fdr01 = vapply(mk, function(d) sum(d$FDR < fdr_check, na.rm = TRUE), 1L))
cat(sprintf("  markers in %.1f min\n", el())); print(n_markers, row.names = FALSE)
all_markers <- unique(unlist(markers))
cat(sprintf("  genes marker of at least one type: %d; of more than one type: %d\n",
            length(all_markers), sum(table(unlist(markers)) > 1)))

# --- pooled profiles mu_E (ADEMP §2.3) and baseline gamma -----------------------------------------
ind <- Matrix::sparseMatrix(i = seq_len(ncol(X)), j = match(cells$cell_type, types), x = 1,
                            dims = c(ncol(X), length(types)))
sums <- as.matrix(X %*% ind); colnames(sums) <- types
mu_E <- sweep(sums, 2, colSums(sums), "/")
base <- Matrix::rowSums(X) / sum(X) * 1e4                       # pooled profile, CP10k
nonmk <- setdiff(rownames(X), all_markers)
b <- base[nonmk]
b <- b[b != 0]
est <- estimate_splat_means(b)
cat(sprintf("  Splatter-style fit on %d non-marker genes (CP10k): gamma shape %.4f, rate %.4f (%s); outliers prob %.4f, facLoc %.4f, facScale %.4f\n",
            length(b), est$shape, est$rate, est$gamma_method, est$out_prob, est$out_facLoc, est$out_facScale))
q <- c(0.1, 0.5, 0.9, 0.99, 0.999)
sim_b <- simulate_splat_means(length(b), est)
print(rbind(empirical = quantile(b, q), simulated = quantile(sim_b, q)), digits = 3)

# --- empirical fold-changes of the markers --------------------------------------------------------
fc <- lapply(types, function(k) {
  g <- markers[[k]]
  if (!length(g)) return(numeric(0))
  other <- rowMeans(mu_E[g, setdiff(types, k), drop = FALSE])
  setNames(mu_E[g, k] / other, g)
})
names(fc) <- types
fc_all <- unlist(fc)
cat(sprintf("  fold-changes: n %d, median %.2f, 10%% %.2f, 90%% %.2f, share below 1: %.4f, infinite: %d\n",
            length(fc_all), median(fc_all[is.finite(fc_all)]), quantile(fc_all[is.finite(fc_all)], 0.1),
            quantile(fc_all[is.finite(fc_all)], 0.9), mean(fc_all < 1), sum(!is.finite(fc_all))))

res <- list(tissue = tissue, types = types, genes = rownames(X), markers = markers, n_markers = n_markers,
            baseline = est, fold_changes = fc,
            settings = list(seed = seed, gene_min_cells = gene_min_cells, fdr = fdr_main, fdr_check = fdr_check,
                            lfc = 1, pval.type = "some", block = "donor", scale = "CP10k"),
            n_cells = ncol(X), minutes = el(), session = sessionInfo())
out <- file.path(out_dir, paste0(tissue, "_p_calibration.rds"))
saveRDS(res, out)
cat(sprintf("  saved %s in %.1f min\n", out, el()))
