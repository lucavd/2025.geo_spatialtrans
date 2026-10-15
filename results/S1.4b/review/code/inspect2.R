
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
F <- read.csv("results/S1.4b/review/code/recomputed_folds.csv")
for (k in 1:5) { o <- readRDS(sprintf("/tmp/code_vc/A6_fold%d.rds", k)); cat("A6 fold", k, "yard d_sp med/max", round(median(o$yard$d_sp),4), round(max(o$yard$d_sp),4),
  " ks_nc med/max", round(median(o$yard$ks_nc),4), round(max(o$yard$ks_nc),4), "\n") }
print(F[F$A == "A6" & F$variant != "V-ab", c("variant","k","ks_nc","d_sp","v1","v3","at_p_vs_at","poi_shift_med")], digits = 4)
# margine minimo |d - soglia| / soglia su tutte le celle PASS/WARN della tabella, per variante (fragilita' del conteggio)
