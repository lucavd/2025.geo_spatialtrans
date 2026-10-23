# tools/S1.5c_bench_lib.R - S1.5c: funzioni condivise da banco, test e controllo della calibrazione.
# compart_for(A, roi, fold): tabella dei compartimenti (cell_type, ell_um, s) della piega; predefinito lascia-fuori del ROI simulato
# (definizione: docs/steps/S1.5c.md, «Banco»). Legge results/S1.5c/S1.5c_compartments_<A>.csv (tools/S1.5c_calibration.py).
compart_for <- function(A, roi, fold = c("loo", "pooled")) {
  fold <- match.arg(fold)
  d <- utils::read.csv(sprintf("results/S1.5c/S1.5c_compartments_%s.csv", A), stringsAsFactors = FALSE)
  key <- if (fold == "loo") paste0("loo_", roi) else "pooled"
  e <- d[d$fold == key, c("cell_type", "ell_um", "s")]
  if (!nrow(e)) stop("compart_for(): piega ", key, " assente per ", A, call. = FALSE)
  rownames(e) <- NULL; e
}
