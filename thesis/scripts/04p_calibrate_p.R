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
# Amendment 2026-11-11 (DECISIONS, Pierangela; before the brain calibration and the lymphoid check):
#   findMarkers on at most 10,000 training cells per type (random within type, seed 42; smaller types whole),
#   in every tissue, because the brain training (897,347 nuclei) does not fit in memory (same reasoning as the
#   20,000-cell sample of φ); gene filter, pooled profile and mu_E stay on all training cells (streamed sums).
#   Pre-registered check in lymphoid (4 types above 10,000 cells): markers with the cap vs all cells, Jaccard
#   of the marker sets (FDR < 0.05) >= 0.9 for every type, otherwise the cap rises to 20,000 in every tissue.
#   Outcome 2026-11-11: the check failed (cd4_cm 0.600; 5 more types 0.750-0.857), so the cap is 20,000 in
#   every tissue; the check at 20,000 is repeated as a reported value only (DECISIONS 2026-11-11).
# Usage (working directory = thesis/):  Rscript scripts/04p_calibrate_p.R <tissue> [all]

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
# cells per type for findMarkers (DECISIONS 2026-11-11); 2nd argument "all" = no cap, written to
# <tissue>_p_calibration_all.rds, used only by the pre-registered lymphoid check (tests/test_p_cap.R)
marker_cap <- if (length(args) >= 2 && args[2] == "all") Inf else 20000   # 10,000 failed the lymphoid check (DECISIONS 2026-11-11)
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
# Cells for findMarkers (DECISIONS 2026-11-11): at most marker_cap cells per type, drawn at random within
# the type (seed 42); every type is kept, a type with <= marker_cap cells is taken whole. Gene filter,
# pooled profile and mu_E are sums over ALL training cells, accumulated window by window (the full matrix
# is never held in memory; PsychAD brain did not fit in 96 GB).
tps <- sort(unique(cells$cell_type))
sub_rows <- sort(unlist(lapply(split(cells$row, cells$cell_type), function(r)
  if (length(r) <= marker_cap) r else r[sample.int(length(r), marker_cap)]), use.names = FALSE))
G0 <- nrow(genes)
detect <- numeric(G0)
sums_all <- matrix(0, G0, length(tps), dimnames = list(genes$gene_id, tps))
mats <- list()
row_type <- setNames(cells$cell_type, cells$row)
stream_counts(file, counts_group, cells$row, G0, function(m, r) {
  stopifnot(all(m@x >= 0), all(m@x == round(m@x)))
  detect <<- detect + Matrix::colSums(m > 1)
  j <- match(row_type[as.character(r)], tps)
  ind <- Matrix::sparseMatrix(i = seq_along(j), j = j, x = 1, dims = c(length(j), length(tps)))
  sums_all <<- sums_all + as.matrix(Matrix::crossprod(m, ind))
  s <- r %in% sub_rows
  if (any(s)) mats[[length(mats) + 1]] <<- m[s, , drop = FALSE]
})
keep <- detect >= gene_min_cells
X <- Matrix::t(do.call(rbind, mats)); rm(mats); invisible(gc())   # genes x cells (marker cells)
rownames(X) <- genes$gene_id
X <- X[keep, ]
sub <- cells[match(sub_rows, cells$row), ]
colnames(X) <- sub$cell_id
cat(sprintf("  read in %.1f min; genes kept %d of %d; marker cells %d of %d (cap %s per type)\n", el(),
            nrow(X), length(keep), ncol(X), nrow(cells), format(marker_cap)))
stopifnot("every type among the marker cells" = setequal(unique(sub$cell_type), tps))

# --- markers --------------------------------------------------------------------------------------
sce <- SingleCellExperiment(list(counts = X),
                            colData = data.frame(cell_type = sub$cell_type, donor = sub$donor_id))
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
stopifnot(identical(sort(types), tps))
sums <- sums_all[keep, types, drop = FALSE]                      # all training cells
mu_E <- sweep(sums, 2, colSums(sums), "/")
base <- rowSums(sums) / sum(sums) * 1e4                          # pooled profile, CP10k
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
                            lfc = 1, pval.type = "some", block = "donor", scale = "CP10k",
                            marker_cap = marker_cap),
            n_cells = nrow(cells), n_marker_cells = ncol(X), minutes = el(), session = sessionInfo())
out <- file.path(out_dir, paste0(tissue, "_p_calibration", if (is.finite(marker_cap)) "" else "_all", ".rds"))
saveRDS(res, out)
cat(sprintf("  saved %s in %.1f min\n", out, el()))
