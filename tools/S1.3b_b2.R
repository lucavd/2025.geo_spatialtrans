# tools/S1.3b_b2.R — S1.3b, B-2: i verdetti del check B di S1.3 ricalcolati con le aree di ogni variante.
# Per ogni variante: copia degli RDS dei nulli di S1.3 con $summary e $cells (aree, interne) della variante; S1.3_analysis.R
# invariato, eseguito in un worktree git temporaneo (scrive in results/S1.3 del worktree, non del repo). V0 = controllo.
# Uso (root del repo): Rscript --vanilla tools/S1.3b_b2.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
OUTB <- commandArgs(TRUE)[1]; if (is.na(OUTB)) OUTB <- "/mnt/micron/geo_spatialtrans/S1.3b"
OUT13 <- "/mnt/micron/geo_spatialtrans/S1.3"; B2 <- file.path(OUTB, "b2"); dir.create(B2, showWarnings = FALSE)
REPO <- normalizePath("."); WT <- file.path(B2, "wt")
if (!dir.exists(WT)) { stopifnot(system(sprintf("git worktree add --detach %s HEAD", WT)) == 0)
  file.symlink(file.path(REPO, "renv/library"), file.path(WT, "renv/library")) }
for (v in c("V0", "V1", "V2", "V3")) {
  dv <- file.path(B2, v); dir.create(file.path(dv, "null"), recursive = TRUE, showWarnings = FALSE)
  for (s in c("real", "cp3")) if (!file.exists(file.path(dv, s))) file.symlink(file.path(OUT13, s), file.path(dv, s))
  for (f in list.files(file.path(OUT13, "null"), pattern = "\\.rds$")) {
    o <- readRDS(file.path(OUT13, "null", f)); b <- readRDS(file.path(OUTB, "null", f)); m <- b$res[[v]]
    o$summary <- m$summary
    if (!is.null(o$cells)) { stopifnot(nrow(o$cells) == nrow(m$cells)); o$cells$area <- m$cells$area; o$cells$interior <- m$cells$interior }
    saveRDS(o, file.path(dv, "null", f))
  }
  st <- system(sprintf("cd %s && Rscript --vanilla tools/S1.3_analysis.R %s > %s/analysis_%s.log 2>&1", WT, dv, B2, v)); stopifnot(st == 0)
  file.copy(file.path(WT, "results/S1.3/S1.3_verdicts.csv"), file.path(B2, sprintf("verdicts_%s.csv", v)), overwrite = TRUE)
  cat(v, "fatto\n")
}
