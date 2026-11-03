# 04m_fit_scdesign3.R — fit scDesign3 for arm M, once per tissue, on the main training split
# (ADEMP §2.1 "Arm M"; DECISIONS 2026-11-02, 2026-11-03: M is a secondary analysis).
#
# Settings written 2026-11-03 before the first run (Pierangela: scDesign3 defaults):
#   - data: main training split only (split == "training"), cell types with >= 50 training cells
#     (type_kept); the held-out group is never read;
#   - genes: the ADEMP §2.6.2 filter, count > 1 in >= 10 training cells (same as 03_reference.R);
#   - marginal: negative binomial (family_use = "nb", scDesign3 default), mu_formula = "cell_type + donor"
#     (cell type and donor as covariates, as for φ), sigma_formula = "1" (default: one dispersion per gene),
#     usebam = FALSE (default);
#   - library size: no explicit offset; scDesign3's default model has none, and the donor term absorbs
#     donor-level depth (declared; spot depth is then set by the common thinning, ADEMP §2.2);
#   - copula: Gaussian (default), corr_formula = "1" (one gene-gene correlation for all cells),
#     important_feature = "all" (default);
#   - parallel: BiocParallel SnowParam (Windows has no forking), workers from TESI_WORKERS.
# Output: <TESI_DATA>/derived/scdesign3/<tissue>_<tag>.rds with data, marginal fits, copula, times.
# Optional 2nd argument = number of genes drawn at random (fixed seed) for a timing test.
#
# Usage (working directory = thesis/):  Rscript scripts/04m_fit_scdesign3.R colon [n_genes]

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
n_test <- if (length(args) >= 2) as.integer(args[2]) else NA
stopifnot(tissue %in% c("brain", "colon", "lymphoid"))
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
out_dir <- file.path(data_dir, "derived", "scdesign3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source("R/read_h5ad.R"); source("R/stream_counts.R")
suppressMessages({ library(scDesign3); library(SingleCellExperiment) })

seed <- 42               # study seed (Pierangela, 2026-11-03)
gene_min_cells <- 10     # ADEMP §2.6.2
workers <- as.integer(Sys.getenv("TESI_WORKERS", min(30, parallel::detectCores() - 2)))
counts_group <- c(brain = "/X", colon = "/raw/X", lymphoid = "/raw/X")[[tissue]]
t0 <- proc.time()[["elapsed"]]
el <- function() (proc.time()[["elapsed"]] - t0) / 60

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
if (!is.na(n_test)) { set.seed(seed); X <- X[sort(sample(nrow(X), n_test)), ] }
cat(sprintf("  read in %.1f min; genes %d (filter kept %d)\n", el(), nrow(X), sum(keep)))

sce <- SingleCellExperiment(list(counts = X),
                            colData = data.frame(cell_type = factor(cells$cell_type),
                                                 donor = factor(cells$donor_id), row.names = cells$cell_id))
# Parallel back-end: forking (mcmapply) on Linux; on Windows SnowParam for the marginals and a serial
# copula, because SnowParam fails to send the copula data to the workers (colon test, 2026-11-03).
on_unix <- .Platform$OS.type == "unix"
par_m <- if (on_unix) "mcmapply" else "bpmapply"
bp <- if (on_unix) NULL else BiocParallel::SnowParam(workers = workers, RNGseed = seed)
tag <- if (is.na(n_test)) "main" else paste0("test", n_test)
set.seed(seed)
dat <- construct_data(sce, assay_use = "counts", celltype = "cell_type", pseudotime = NULL, spatial = NULL,
                      other_covariates = "donor", corr_by = "1")
t_m <- el()
marginal <- fit_marginal(data = dat, predictor = "gene", mu_formula = "cell_type + donor", sigma_formula = "1",
                         family_use = "nb", n_cores = workers, usebam = FALSE,
                         parallelization = par_m, BPPARAM = bp)
t_m <- el() - t_m
cat(sprintf("  marginals: %d genes in %.1f min (%.2f s per gene with %d workers)\n",
            nrow(X), t_m, 60 * t_m / nrow(X), workers))
saveRDS(list(data = dat, marginal = marginal, workers = workers, minutes = t_m),
        file.path(out_dir, paste0(tissue, "_", tag, "_marginal.rds")))   # kept if the copula fails
t_c <- el()
set.seed(seed)
copula <- fit_copula(sce = sce, assay_use = "counts", input_data = dat$dat, marginal_list = marginal,
                     family_use = "nb", copula = "gaussian", important_feature = "all",
                     n_cores = if (on_unix) workers else 1, parallelization = "mcmapply")
t_c <- el() - t_c
cat(sprintf("  copula in %.1f min\n", t_c))

saveRDS(list(data = dat, marginal = marginal, copula = copula, genes = rownames(X), seed = seed, workers = workers,
             minutes = c(marginal = t_m, copula = t_c, total = el()), date = Sys.time(),
             version = as.character(packageVersion("scDesign3"))),
        file.path(out_dir, paste0(tissue, "_", tag, ".rds")))
cat(sprintf("  saved %s_%s.rds; total %.1f min\n", tissue, tag, el()))
