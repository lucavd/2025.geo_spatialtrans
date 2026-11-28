# 09a_export_c2l.R — exports the inputs of cell2location (Python, GPU) to plain files that Python reads
# without R (ADEMP §4, §6; DECISIONS 2026-11-04). No new random numbers, no change to the datasets:
# the same spots and references read by 05_deconvolve.R are written as MatrixMarket counts.
#
# Written to <TESI_DATA>/c2l/<mode>/<tissue>/:
#   ref_real/            counts.mtx (genes x cells), genes.txt, cells.csv (cell, cell_type, donor_id)
#   r<rep>/ref_P/        same files for the reference of arm P, rebuilt from world_p.rds as in 05
#   r<rep>/<cond>/<dataset>/  counts.mtx (genes x spots), genes.txt, spots.txt, n_cells.txt
#   types.txt            reference cell types, in the order of 05_deconvolve.R
# Files already present are not rewritten, so the script can be rerun.
#
# Usage (working directory = thesis/):
#   Rscript scripts/09a_export_c2l.R <tissue> [pilot|main] [rep_from] [rep_to]

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
ref_mode <- if (length(args) >= 2) args[2] else "main"
rep_from <- if (length(args) >= 3) as.integer(args[3]) else 1L
rep_to <- if (length(args) >= 4) as.integer(args[4]) else rep_from
stopifnot(tissue %in% c("colon", "lymphoid", "brain"), ref_mode %in% c("pilot", "main"))

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
suppressPackageStartupMessages(library(Matrix))
for (f in c("R/simulate_spots.R", "R/arm_p.R")) source(f)

R0 <- readRDS(file.path(data_dir, "derived", "reference", sprintf("%s_%s.rds", tissue, ref_mode)))
sim_dir <- file.path(data_dir, "sim", ref_mode, tissue)
out_dir <- file.path(data_dir, "c2l", ref_mode, tissue)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
L_ref <- Matrix::colSums(R0$ref)
hd_depth <- c(colon = 174, lymphoid = 62, brain = 121)
cond_keys <- sprintf("c%02d_d%04d", c(2, 10), hd_depth[[tissue]])   # as in 05_deconvolve.R

write_counts <- function(m, d, col_file, col_df = NULL) {
  if (file.exists(file.path(d, "counts.mtx"))) return(invisible(FALSE))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  m <- as(as(m, "CsparseMatrix"), "generalMatrix")
  if (any(m@x != round(m@x))) stop("non-integer counts in ", d)   # cell2location needs raw counts
  writeLines(rownames(m), file.path(d, "genes.txt"))
  if (is.null(col_df)) writeLines(colnames(m), file.path(d, col_file))
  else write.csv(col_df, file.path(d, col_file), row.names = FALSE)
  tmp <- file.path(d, "counts.mtx.tmp"); Matrix::writeMM(m, tmp); file.rename(tmp, file.path(d, "counts.mtx"))
  invisible(TRUE)
}
ref_df <- function(ref, meta) data.frame(cell = colnames(ref), cell_type = meta$cell_type, donor_id = meta$donor_id)

writeLines(R0$types_ref, file.path(out_dir, "types.txt"))
write_counts(R0$ref, file.path(out_dir, "ref_real"), "cells.csv", ref_df(R0$ref, R0$ref_meta))

for (i in rep_from:rep_to) {
  rd <- file.path(sim_dir, sprintf("r%03d", i))
  od <- file.path(out_dir, sprintf("r%03d", i))
  RP <- reference_p_from_world(readRDS(file.path(rd, "world_p.rds")), R0$ref_meta, L_ref, colnames(R0$ref))
  if (ncol(RP$ref) != ncol(R0$ref) || !identical(RP$ref_meta, R0$ref_meta)) stop("D2")   # check of 05
  write_counts(RP$ref, file.path(od, "ref_P"), "cells.csv", ref_df(RP$ref, RP$ref_meta))
  for (key in cond_keys) {
    x <- readRDS(file.path(rd, paste0(key, ".rds")))
    for (d in names(x$counts)) {
      ref_genes <- if (d == "P") rownames(RP$ref) else rownames(R0$ref)
      if (!identical(rownames(x$counts[[d]]), ref_genes)) stop("D1: ", d)   # check of 05
      dd <- file.path(od, key, d)
      write_counts(x$counts[[d]], dd, "spots.txt")
      writeLines(as.character(x$meta$n_cells), file.path(dd, "n_cells.txt"))
    }
  }
  cat(sprintf("[%s] rep %d exported to %s\n", tissue, i, od))
}
