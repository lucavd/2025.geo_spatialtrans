# tools/S1.5e_bench_lib.R - S1.5e: parametri dei vuoti calibrati per famiglia e piega (results/S1.5e/S1.5e_het_<famiglia>.csv,
# tools/S1.5e_calibration_pos.py), con le colonne void_mode, void_aspect, dens_cap. Predefinito: lascia-fuori della ROI simulata.
het_for_e <- function(A, roi, family, fold = c("loo", "pooled")) {
  fold <- match.arg(fold)
  d <- utils::read.csv(sprintf("results/S1.5e/S1.5e_het_%s.csv", family), stringsAsFactors = FALSE)
  key <- if (fold == "loo") paste0("loo_", roi) else "pooled"
  e <- d[d$archetype == A & d$fold == key, c("void_frac", "void_ell_um", "dens_sigma", "dens_ell_um", "void_mode", "void_aspect", "dens_cap")]
  if (nrow(e) != 1L) stop("het_for_e(): piega ", key, " assente o duplicata per ", A, " (", family, ")", call. = FALSE)
  e$dens_cap <- as.logical(e$dens_cap); rownames(e) <- NULL; e
}
