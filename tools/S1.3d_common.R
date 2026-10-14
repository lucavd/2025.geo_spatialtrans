# tools/S1.3d_common.R — S1.3d (BL-084): caricamento comune. Pre-registrazione results/S1.3d/S1.3d_preregistration.md (ee6208b).
# Rigenera i nulli esattamente come tools/S1.3b_run.R (stessa regen()), carica gli helper di S1.3b e i pacchetti dei candidati.
.s13d_args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) .s13d_args
source("tools/S1.3_roi.R")          # rois, ARCH, gen_points, s13_tessellate, c123, one_region, R3, TOL, BASE_SEED, NREP
rm(commandArgs)
source("tools/S1.3b_overlay.R")     # s3b_* (varianti V0..V3 di S1.3b)
suppressPackageStartupMessages({ library(polyclip); library(exactextractr); library(terra) })
OUT13 <- "/mnt/micron/geo_spatialtrans/S1.3"; SEL_SEED <- 20261010; MODELS <- c("CSR", "RSA", "RSArule")
G_GRID <- 1e-6                                   # g (pre-registrato) per tutte le tolleranze
regen <- function(i, model, rp) {                 # identica a tools/S1.3b_run.R
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]
  w <- r3_roi_window(A, roi, rois)
  n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  p <- gen_points(i, model, rp, w, n_obs, a)
  s <- s13_tessellate(p$x, p$y, w, backend = "geos")
  list(A = A, roi = roi, w = w, p = p, s = s, terr = s$tv$cell_territories, reg = st_sfc(w$poly[[1]]))
}
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)
fname <- function(A, roi, model, rp) sprintf("%s_%s_%s_%02d.rds", A, roi, model, rp)
