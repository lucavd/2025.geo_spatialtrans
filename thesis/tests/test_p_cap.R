# test_p_cap.R — pre-registered check of the marker cap of arm P (DECISIONS 2026-11-11, written before the run).
# Lymphoid has 4 training types above 10,000 cells and fits in memory: markers (FDR < 0.05) with the cap
# (<tissue>_p_calibration.rds) vs all cells (<tissue>_p_calibration_all.rds). The cap is adequate if the
# Jaccard index of the two marker sets is >= 0.9 for every type; otherwise the cap rises to 20,000 in every
# tissue. Also reported: markers per type in both runs; gene list, baseline and fold-changes of the types
# below the cap must be computed on the same cells (gene list identical).
# Usage (working directory = thesis/):  Rscript scripts/04p_calibrate_p.R lymphoid
#                                       Rscript scripts/04p_calibrate_p.R lymphoid all
#                                       Rscript tests/test_p_cap.R [lymphoid]
args <- commandArgs(trailingOnly = TRUE); tissue <- if (length(args)) args[1] else "lymphoid"
d <- file.path(Sys.getenv("TESI_DATA"), "derived", "p_calibration")
cap <- readRDS(file.path(d, paste0(tissue, "_p_calibration.rds")))
all <- readRDS(file.path(d, paste0(tissue, "_p_calibration_all.rds")))
stopifnot(is.finite(cap$settings$marker_cap), is.infinite(all$settings$marker_cap))
cat(sprintf("Cap %s cells per type%s\n", format(cap$settings$marker_cap), if (cap$settings$marker_cap == 20000) " (2nd check, reported only: DECISIONS 2026-11-11)" else ""))
ok <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); invisible(cond) }
ok(identical(cap$genes, all$genes), "same kept genes (filter on all training cells)")
ok(identical(cap$types, all$types), "same types")
jac <- function(a, b) if (!length(union(a, b))) 1 else length(intersect(a, b)) / length(union(a, b))
n_cells <- table(read.csv(file.path(Sys.getenv("TESI_DATA"), "derived", "split", paste0(tissue, "_cells.csv.gz")),
                          stringsAsFactors = FALSE)[, c("split", "type_kept", "cell_type")] |>
                   subset(split == "training" & type_kept, select = cell_type))
tab <- data.frame(cell_type = cap$types, n_train = as.integer(n_cells[cap$types]),
                  markers_cap = lengths(cap$markers)[cap$types], markers_all = lengths(all$markers)[cap$types],
                  jaccard = mapply(jac, cap$markers[cap$types], all$markers[cap$types]))
print(tab[order(-tab$n_train), ], row.names = FALSE, digits = 3)
ok(all(tab$jaccard >= 0.9), sprintf("Jaccard >= 0.9 for every type (min %.3f, %s)", min(tab$jaccard),
                                    tab$cell_type[which.min(tab$jaccard)]))
