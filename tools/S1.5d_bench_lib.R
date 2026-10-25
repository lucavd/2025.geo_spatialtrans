# tools/S1.5d_bench_lib.R - S1.5d: funzioni condivise da calibrazione dell'eterogeneita', banco, test e controllo della calibrazione.
# compart_for_d(A, roi, fold): compartimenti ricalibrati (results/S1.5d/S1.5d_compartments_<A>.csv, tools/S1.5d_calibration_types.py) della piega;
# het_for(A, roi, family, fold): eterogeneita' calibrata (results/S1.5d/S1.5d_het_<family>.csv, tools/S1.5d_calibration_pos.py) della piega.
# Predefinito: lascia-fuori della ROI simulata (docs/steps/S1.5d.md, «Banco»).
compart_for_d <- function(A, roi, fold = c("loo", "pooled")) {
  fold <- match.arg(fold)
  d <- utils::read.csv(sprintf("results/S1.5d/S1.5d_compartments_%s.csv", A), stringsAsFactors = FALSE)
  key <- if (fold == "loo") paste0("loo_", roi) else "pooled"
  e <- d[d$fold == key, c("cell_type", "ell_um", "s")]
  if (!nrow(e)) stop("compart_for_d(): piega ", key, " assente per ", A, call. = FALSE)
  rownames(e) <- NULL; e
}
het_for <- function(A, roi, family, fold = c("loo", "pooled")) {
  fold <- match.arg(fold)
  d <- utils::read.csv(sprintf("results/S1.5d/S1.5d_het_%s.csv", family), stringsAsFactors = FALSE)
  key <- if (fold == "loo") paste0("loo_", roi) else "pooled"
  e <- d[d$archetype == A & d$fold == key, c("void_frac", "void_ell_um", "dens_sigma", "dens_ell_um")]
  if (nrow(e) != 1L) stop("het_for(): piega ", key, " assente o duplicata per ", A, " (", family, ")", call. = FALSE)
  rownames(e) <- NULL; e
}
