# test_music_basis.R — check of the MuSiC basis computed once per tissue (DECISIONS 2026-11-12, Luca),
# written 2026-11-17 before any use in the main run. On datasets written by 04_simulate.R, the proportions
# of run_music(basis = music_basis_once(...)) must be identical to those of the original
# MuSiC::music_prop() call (run_music(basis = NULL)); also reported: the largest absolute difference and
# the time saved (basis once + all datasets, against music_prop() on every dataset).
# Arm P has its own reference in every replicate, so its basis is rebuilt per replicate (checked too).
# Usage (working directory = thesis/):  Rscript tests/test_music_basis.R [tissue] [pilot|main] [n_rep]
args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
ref_mode <- if (length(args) >= 2) args[2] else "pilot"
n_rep <- if (length(args) >= 3) as.integer(args[3]) else 2L
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), file.exists("thesis.Rproj"))
suppressPackageStartupMessages(library(Matrix))
for (f in c("R/simulate_spots.R", "R/arm_p.R", "R/deconvolve.R")) source(f)
ok <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); invisible(cond) }

R0 <- readRDS(file.path(data_dir, "derived", "reference", sprintf("%s_%s.rds", tissue, ref_mode)))
sim_dir <- file.path(data_dir, "sim", ref_mode, tissue)
reps <- head(sort(list.files(sim_dir, pattern = "^r[0-9]+$")), n_rep)
stopifnot("no replicates on disk" = length(reps) > 0)
L_ref <- Matrix::colSums(R0$ref)

t_basis <- system.time(B <- music_basis_once(R0$ref, R0$ref_meta$cell_type, R0$ref_meta$donor_id))[["elapsed"]]
rows <- list()
for (rp in reps) {
  rd <- file.path(sim_dir, rp)
  keys <- sub("\\.rds$", "", list.files(rd, pattern = "^c[0-9]+_d[0-9]+\\.rds$"))
  BP <- NULL; tBP <- 0
  for (key in keys) {
    x <- readRDS(file.path(rd, paste0(key, ".rds")))
    for (d in names(x$counts)) {
      if (d == "P") {
        RP <- reference_p_from_world(readRDS(file.path(rd, "world_p.rds")), R0$ref_meta, L_ref, colnames(R0$ref))
        ref <- RP$ref
        if (is.null(BP)) tBP <- system.time(BP <- music_basis_once(ref, RP$ref_meta$cell_type, RP$ref_meta$donor_id))[["elapsed"]]
        b <- BP; tb <- tBP; tBP <- 0   # P basis time counted once per replicate
      } else { ref <- R0$ref; b <- B; tb <- 0 }
      a0 <- run_music(x$counts[[d]], ref, R0$ref_meta$cell_type, R0$ref_meta$donor_id)
      a1 <- run_music(x$counts[[d]], ref, R0$ref_meta$cell_type, R0$ref_meta$donor_id, basis = b)
      cat(rp, key, d, identical(a0$prop, a1$prop), round(a0$seconds, 1), round(a1$seconds, 1), "\n")
      rows[[length(rows) + 1]] <- data.frame(rep = rp, cond = key, dataset = d,
        identical = identical(a0$prop, a1$prop),
        max_abs_diff = if (identical(dim(a0$prop), dim(a1$prop))) max(abs(a0$prop - a1$prop)) else NA,
        sec_music_prop = a0$seconds, sec_basis_once = a1$seconds + tb)
    }
  }
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE, digits = 3)
cat(sprintf("\nBasis of the real reference: %.1f s (once per tissue)\n", t_basis))
cat(sprintf("Total over %d datasets: music_prop() %.0f s; basis once %.0f s (incl. real basis) -> %.0f%% of the time\n",
            nrow(tab), sum(tab$sec_music_prop), sum(tab$sec_basis_once) + t_basis,
            100 * (sum(tab$sec_basis_once) + t_basis) / sum(tab$sec_music_prop)))
ok(all(tab$identical), sprintf("proportions identical to music_prop() in %d of %d datasets (max |diff| %.3g)",
                              sum(tab$identical), nrow(tab), max(tab$max_abs_diff, na.rm = TRUE)))
