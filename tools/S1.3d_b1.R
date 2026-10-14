# tools/S1.3d_b1.R — S1.3d, B-1: i verdetti del check B di S1.3 ricalcolati con le aree di ogni candidato (come tools/S1.3b_b2.R).
# Per ogni candidato: copie degli RDS dei nulli di S1.3 con $summary e $cells (aree, interne) del candidato (calcolati in cands_one da
# r3_summary sulla tabella per cella di S1.3); tools/S1.3_analysis.R invariato in un worktree git temporaneo. V0 = controllo (deve
# riprodurre S1.3). Uso (root del repo): Rscript --vanilla tools/S1.3d_b1.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
OUTD <- commandArgs(TRUE)[1]; if (is.na(OUTD)) OUTD <- "/mnt/micron/geo_spatialtrans/S1.3d"
OUT13 <- "/mnt/micron/geo_spatialtrans/S1.3"; B1 <- file.path(OUTD, "b1"); dir.create(B1, showWarnings = FALSE); RES <- file.path(OUTD, "res")
REPO <- normalizePath("."); WT <- file.path(B1, "wt")
if (!dir.exists(WT)) { stopifnot(system(sprintf("git worktree add --detach %s HEAD", WT)) == 0)
  file.symlink(file.path(REPO, "renv/library"), file.path(WT, "renv/library")) }
CANDS <- c("V0", "F1", "F2", "F3", "F4", "F5")
nf <- list.files(file.path(OUT13, "null"), pattern = "\\.rds$")
for (v in CANDS) {
  dv <- file.path(B1, v); dir.create(file.path(dv, "null"), recursive = TRUE, showWarnings = FALSE)
  for (s in c("real", "cp3")) if (!file.exists(file.path(dv, s))) file.symlink(file.path(OUT13, s), file.path(dv, s))
  for (f in nf) {
    key <- sub("\\.rds$", "", f); rf <- file.path(RES, paste0(key, c("_o12.rds", "_o1.rds"))); rf <- rf[file.exists(rf)][1]
    stopifnot(!is.na(rf)); r <- readRDS(rf); m <- r$b1[[v]]; stopifnot(!is.null(m))
    o <- readRDS(file.path(OUT13, "null", f)); o$summary <- m$summary
    if (!is.null(o$cells)) { stopifnot(nrow(o$cells) == nrow(m$cells)); o$cells$area <- m$cells$area; o$cells$interior <- m$cells$interior }
    saveRDS(o, file.path(dv, "null", f))
  }
  st <- system(sprintf("cd %s && Rscript --vanilla tools/S1.3_analysis.R %s > %s/analysis_%s.log 2>&1", WT, dv, B1, v)); stopifnot(st == 0)
  file.copy(file.path(WT, "results/S1.3/S1.3_verdicts.csv"), file.path(B1, sprintf("verdicts_%s.csv", v)), overwrite = TRUE)
  cat(v, "fatto\n")
}
ref <- read.csv("results/S1.3/S1.3_verdicts.csv"); kB <- ref$section == "B" & !is.na(ref$confirmed); kb <- ref$section == "B"
out <- do.call(rbind, lapply(CANDS, function(v) { x <- read.csv(file.path(B1, sprintf("verdicts_%s.csv", v)))
  data.frame(cand = v, all_identical = identical(x[, c("section", "id", "value", "outcome", "confirmed")], ref[, c("section", "id", "value", "outcome", "confirmed")]),
             pred_same = sum(x$outcome[kB] == ref$outcome[kB] & x$confirmed[kB] == ref$confirmed[kB]), n_pred = sum(kB),
             value_same = sum(x$value[kb] == ref$value[kb]), n_rows_B = sum(kb)) }))
print(out); write.csv(out, "results/S1.3d/S1.3d_b1.csv", row.names = FALSE)
