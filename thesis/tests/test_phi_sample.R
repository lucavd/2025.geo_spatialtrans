# test_phi_sample.R — pre-registered check of the φ sample (DECISIONS 2026-10-29, refined 2026-10-30;
# ADEMP §2.3). Written before the comparison was run.
# In colon all training cells fit in memory, so φ_g fitted on the 20,000-cell stratified sample is
# compared with φ_g fitted on all training cells (same genes, same design ~ cell_type + donor).
# The sample is adequate if
#   1. the Pearson correlation of log φ across genes is >= 0.99, and
#   2. the median of φ_sample / φ_all is between 0.95 and 1.05 (no systematic bias).
# If either fails, the sample rises to 50,000 cells in every tissue (written as an amendment).
# Needs: Rscript scripts/03_reference.R colon pilot      -> colon_pilot.rds
#        Rscript scripts/03_reference.R colon pilot all  -> colon_pilot_phiall.rds
# Run from thesis/:  Rscript tests/test_phi_sample.R [tissue] [mode]

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
mode <- if (length(args) >= 2) args[2] else "pilot"
dir <- file.path(Sys.getenv("TESI_DATA"), "derived", "reference")
s <- readRDS(file.path(dir, sprintf("%s_%s.rds", tissue, mode)))
a <- readRDS(file.path(dir, sprintf("%s_%s_phiall.rds", tissue, mode)))
stopifnot("Different gene lists" = identical(names(s$phi), names(a$phi)))

ok <- s$phi > 1e-6 & a$phi > 1e-6   # log scale needs φ > 0; genes at the lower bound reported apart
r <- cor(log(s$phi[ok]), log(a$phi[ok]))
ratio <- median(s$phi[ok] / a$phi[ok])
cat(sprintf("[%s %s] φ on %d sampled cells vs %d cells; %d genes (%d with φ ~ 0 in either fit)\n",
            tissue, mode, s$settings$phi_cells, a$settings$phi_cells, length(ok), sum(!ok)))
cat(sprintf("  correlation of log φ: %.4f  -> %s (>= 0.99)\n", r, if (r >= 0.99) "PASS" else "FAIL"))
cat(sprintf("  median ratio sample/all: %.3f -> %s (0.95-1.05)\n", ratio,
            if (ratio >= 0.95 && ratio <= 1.05) "PASS" else "FAIL"))
cat(sprintf("  genes with ratio outside 0.5-2: %d\n", sum(s$phi[ok] / a$phi[ok] < 0.5 | s$phi[ok] / a$phi[ok] > 2)))
