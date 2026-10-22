# tools/S1.5b_calcheck.R - S1.5b, cancello G-4: ricalcolo in R dei modelli di taglia dagli stessi file della calibrazione (Python).
# (1) modello di un solo gruppo = preset di S1.5 in intestino, tumore, stroma, cervello (<= 5e-4); (2) modelli dei gruppi = calibrazione (<= 1e-9).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
M <- "/mnt/micron/geo_spatialtrans/R6"; cal <- read.csv("results/S1.5b/S1.5b_type_calibration.csv", stringsAsFactors = FALSE)
EPI <- c("epiteliale normale" = "epiteliale", "epiteliale tumorale" = "epiteliale"); out <- list()
for (A in paste0("A", 1:6)) {
  lab <- as.data.frame(read_parquet(sprintf("%s/labels/%s_exp_umi20_seed20261022_nuclei_ph1.parquet", M, A)))
  lab$type <- ifelse(lab$cell_type %in% names(EPI), EPI[lab$cell_type], lab$cell_type)
  d <- do.call(rbind, lapply(sort(unique(lab$roi_id)), function(r) {
    t <- as.data.frame(read_parquet(sprintf("%s/terr/%s_%s.parquet", M, A, r))); l <- lab[lab$roi_id == r, ]
    k <- match(sprintf("%.6f %.6f", t$x, t$y), sprintf("%.6f %.6f", l$x_um, l$y_um))
    if (anyNA(k)) stop(A, " ", r, ": territori senza etichetta")
    data.frame(roi = r, a_n = t$a_n, a_t = t$a_t, ev = t$interior & t$matched & t$a_n > 0, type = l$type[k]) }))
  d <- d[d$ev, ]; tp <- .TP_TABLE[.TP_TABLE$archetype == A, ]; aref <- 1e6 / tp$density_model
  f1 <- fit_nucleus_size_model(d$a_n, d$a_t, aref, "full")
  out[[length(out) + 1]] <- data.frame(archetype = A, what = "single_vs_preset", group = NA, rho = f1[["rho"]], beta = f1[["beta"]], sigma = f1[["sigma"]],
    ref_rho = tp$rho, ref_beta = tp$beta, ref_sigma = tp$sigma, tol = 5e-4, applies = !A %in% c("A4", "A6"))
  cA <- cal[cal$archetype == A, ]
  if (length(unique(cA$size_group)) > 1) for (g in sort(unique(cA$size_group))) {
    ty <- cA$cell_type[cA$size_group == g]; e <- d[d$type %in% ty, ]; fg <- fit_nucleus_size_model(e$a_n, e$a_t, aref, "full"); r0 <- cA[cA$size_group == g, ][1, ]
    out[[length(out) + 1]] <- data.frame(archetype = A, what = "group_vs_python", group = g, rho = fg[["rho"]], beta = fg[["beta"]], sigma = fg[["sigma"]],
      ref_rho = r0$rho_group, ref_beta = r0$beta_group, ref_sigma = r0$sigma_group, tol = 1e-9, applies = TRUE)
  }
}
O <- do.call(rbind, out)
O$maxdiff <- pmax(abs(O$rho - O$ref_rho), abs(O$beta - O$ref_beta), abs(O$sigma - O$ref_sigma))
O$pass <- !O$applies | O$maxdiff <= O$tol
write.csv(O, "results/S1.5b/S1.5b_G4_calcheck.csv", row.names = FALSE); print(O[, c("archetype", "what", "group", "maxdiff", "tol", "applies", "pass")])
cat(sprintf("%s G-4: %d/%d confronti applicabili entro tolleranza\n", if (all(O$pass)) "PASS" else "FAIL", sum(O$pass & O$applies), sum(O$applies)))
